#install.packages('Seurat')
library(Seurat)

# Load the SAN dataset
san.data <- Read10X(data.dir = "C:/Users/siddm/Desktop/Loyola Academics/Loyola Senior Year Semester 2/COMP 383 Comp Bio - Wheeler/scRNA-Seq-Variation-Pipeline/mouse_heart_GEO_data/SAN_GEO")
dim(san.data)

# Initialize the Seurat object with the raw data
# remove genes expressed less than 3 cells
# remove cells with less than 200 gene expressed
san <- CreateSeuratObject(counts = san.data, min.cells = 3, min.features = 200)
san

# create a subset with cells that have the gene value Mybphl > 0
san.ourgene = subset(x = san, subset = Mybphl > 0)
san.ourgene

cellsWithMYBPHL <- GetAssayData(object = san.ourgene, assay = "RNA", slot = "data")
write.csv(cellsWithMYBPHL, "cellsWithMYBPHL.csv")