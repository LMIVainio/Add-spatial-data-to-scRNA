# This is a simple way of estimating the spatial position in you single cell dataset using the AspWood dataset (not provided here)
# The method sums up each cells general gene expression location and takes the maximum value and saves it as a location data estimator in
# Seurat@meta.data

# Set up libraries
library(readxl)
library(dplyr)
library(Seurat)

## Load your seurat data
your_seurat_object <- LoadSeuratRds("path to your Seurat object")

# Read AspWood data
AspWood <- read_excel("path to your AspWood data set excel file", sheet = "Annotated genes")
AspWood <- column_to_rownames(AspWood, var = "Genes")

# Get an individual cell expression data (in this case SCT normalized expression data)
assay_dataset <- GetAssayData(your_seurat_object, assay = "SCT", slot = "data")

# Select tree as separate entity
tree_name <- "T1"
AspWood %>% select(contains(tree_name)) -> tree_asp

# Loop trough each cell
for (cell_n in 1:NCOL(assay_dataset)){
  # Get the expression data of individual cell
  cell_expressions <- assay_dataset[, cell_n]
  # Generate an empty matrix containing tree specific column names (section numbers)
  cell_location <- matrix(0, nrow = 1, ncol = ncol(tree_asp),
                                dimnames = list("Cell total expession",
                                                colnames(tree_asp)))
  # Loop through each gene in cell
  for (gene_n in 1:NROW(cell_expressions)){
    # get the gene name and format to aspwood
    gene_name <- names(cell_expressions)[gene_n]
    # get the expression level of that gene
    gene_expression <- cell_expressions[gene_n]
    # Skip to next gene if expression is 0 at this cell
    if (gene_expression == 0){next}
    # take a subset from AspWood containing only the location data for
    # gene n and multiply it by expression value 
    expression_location <- tree_asp %>% filter(row.names(tree_asp) %in% gene_name) * gene_expression
    # if not found in asp_wood, accumulate the distribution by adding it to previous gene
    if (nrow(expression_location) != 0) {
      cell_location <- cell_location + expression_location
    }
  }
  # Divide cell_location by it's max value
  cell_location <- cell_location/max(cell_location)
  # Add location estimation to your metadata by cell_name
  # get cell name
  cell_name <- colnames(assay_dataset)[cell_n]
  #AddMetaData(object = your_seurat_object, metadata = cell_location, key = cell_name)
  your_seurat_object@meta.data[cell_name, "cell_location"] <- as.vector(which.max(cell_location))
  # Print progress if necessary
  #print(paste("Cell", cell_n, "out of", NCOL(assay_dataset)))
} 
# Save seurat object
SaveSeuratRds(your_seurat_object, file = "path to your Seurat object")
# Finally plot your data
DimPlot(your_seurat_object, group.by = 'cell_location')
