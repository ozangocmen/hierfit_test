#installing 
install.packages("devtools")
devtools::install_github("yasinkaymaz/HieRFIT")
library(HieRFIT)

##recommended packages
install.packages("DiagrammeR")
install.packages("e1071")
install.packages("alluvial")
install.packages("ggalluvial")

library(DiagrammeR)
library(e1071)
library(alluvial)
library(ggalluvial)

#Seurat package
library(Seurat)

#loading and updating object
setwd("C:/Users/OZAN/Downloads")
pbmc <- readRDS("pbmc3k_final.Rds")
pbmc <- UpdateSeuratObject(pbmc)

#frequency of cell types in the data
table(pbmc@meta.data$ClusterNames_0.6)
classLabels <- pbmc@meta.data$ClusterNames_0.6
head(classLabels)

#reading the tree file 
setwd("C:/Users/OZAN/Desktop")
treeTable <- read.delim("pbmc3k_taxa.txt", header = F)
PlotTopoTree(obj = treeTable)

#generating the model
refmod <- CreateHieR(RefData = pbmc[["RNA"]]@data,
                     ClassLabels = pbmc@meta.data$ClusterNames_0.6,
                     Tree = treeTable,
                     species = "hsapiens")
SaveHieRMod(refMod = refmod, filePrefix = "PBMC3K_HierMod")

#classification accuracy of the mode
PlotTopoNodeAcc(refMod = refmod)

###################  ##########################################################
### 2nd SESSION ###  #Projecting the reference cell type on to a query dataset#
###################  ##########################################################

setwd("C:/Users/OZAN/Desktop/datasets/filtered_feature_bc_matrix")
list.files()
new.pbmc.data <- Read10X("C:/Users/OZAN/Desktop/datasets/filtered_feature_bc_matrix")
newPBMC <- CreateSeuratObject(counts = new.pbmc.data, project = "pbmc10k", min.cells = 3, min.features = 200)


#analysis steps
newPBMC <- NormalizeData(newPBMC)

newPBMC <- FindVariableFeatures(newPBMC, selection.method = "vst", nfeatures = 2000)
newPBMC <- ScaleData(newPBMC)
newPBMC <- RunPCA(newPBMC)
newPBMC <- FindNeighbors(newPBMC, dims = 1:10)
newPBMC <- FindClusters(newPBMC, resolution = 1)

setwd("C:/Users/OZAN/Desktop")
refmod <- readRDS( "PBMC3K_HierMod.RDS")
hierObj <- HieRFIT(Query = newPBMC[["RNA"]]@data, refMod = refmod)
head(hierObj@Evaluation)

#Exploring the projection results:
PlotBarStats(HieRobj = hierObj)
PlotTopoStats(HieRobj = hierObj)

newPBMC@meta.data$res.1
CrossCheck(HieRobj = hierObj, Prior = newPBMC@meta.data$res.1)

sessionInfo()
