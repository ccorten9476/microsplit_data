#Read in data - sample 1
library(readr)
barcodes <- read_csv("GSM4594094_M11_barcodes.csv", col_names = FALSE)
matrix <- read_csv("GSM4594094_M11_dcm.csv", col_names = FALSE)
features<- read_csv("GSM4594094_M11_genes.csv", col_names = FALSE)

#Merge data
library(data.table)
matrix <- as.data.frame(matrix)
matrix <- transpose(matrix)
features_vector <- features$X1
barcodes_vector <- barcodes$X1
row.names(matrix) <- features_vector
colnames(matrix) <- barcodes_vector
#Calculate percentage of zeros in data matrix
total.zeros <- sum(matrix==0)
total.counts <- nrow(matrix) * ncol(matrix)
(total.zeros/total.counts) * 100
rm(total.counts)
rm(total.zeros)
#Create Seurat Object and inspect the data
library(Seurat)
seuratObj <- CreateSeuratObject(counts=matrix, min.cells=3) #Has already been filtered for genes >200
head(Cells(seuratObj))
Cells(seuratObj)
Features(seuratObj)
dim(seuratObj)
seuratObj$sample <- "Sample1"
metaData <- seuratObj@meta.data
seuratObj@meta.data
head(metaData)
head(seuratObj@meta.data)

#
library(ggplot2)
p1 <- ggplot(metaData,aes("",nCount_RNA)) +
  geom_jitter(height=0,width=0.3,color="lightblue") +
  geom_violin(fill="gray80",alpha=0.5) + 
  scale_y_continuous(trans="log2",breaks=c(500,1000,2000,
                                           50000,75000,
                                           100000,150000)) +
  labs(title="Total UMI counts per cell") +
  theme_classic()
p1

p2 <- ggplot(metaData,aes("",nFeature_RNA)) +
  geom_jitter(height=0,width=0.3,color="lightblue") +
  geom_violin(fill="gray80",alpha=0.5) + 
  scale_y_continuous(trans="log2",breaks=c(200,300,400,
                                           5000,7500,10000)) +
  labs(title="Number of genes per cell") +
  theme_classic()
p2
library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

#Normalize the data
seuratObj <- NormalizeData(object=seuratObj)
seuratObj@assays$RNA@layers$counts[1:5,1:10]

#Find highly variable genes
seuratObj <- FindVariableFeatures(seuratObj)
length(VariableFeatures(seuratObj)) 
head(VariableFeatures(seuratObj)) 

VariableFeaturePlot(seuratObj)
p <- VariableFeaturePlot(seuratObj)
top20 <- head(VariableFeatures(seuratObj),20)
LabelPoints(plot=p,points=top20,repel=TRUE)
p2 <- LabelPoints(plot=p,points=top20,repel=TRUE)
#All genes seem highly variable so we take all of them for further analysis

#Scaling the data
seuratObj <- ScaleData(seuratObj, features = rownames(seuratObj))
seuratObj@assays$RNA@layers$scale.data[1:5,1:4]
dim(seuratObj@assays$RNA@layers$scale.data)
#PCA
seuratObj <- RunPCA(seuratObj,nfeatures.print=10)

seuratObj@reductions$pca@cell.embeddings[1:5,1:5]
seuratObj@reductions$pca@feature.loadings[1:5,1:5]
DimPlot(seuratObj,reduction="pca") + theme_few() + NoLegend() + labs(x = "PCA1",y = "PCA2",)
graph2ppt(file="PCA_MICROSPLITplot_SP1.pptx", width=6, height=5)   

#Check the PC's
DimHeatmap(seuratObj,dims=1:12,cells=500,balanced=TRUE)
ElbowPlot(seuratObj,ndims=50)

#Find Neighbours
seuratObj <- FindNeighbors(seuratObj,dims=1:4)
names(seuratObj@graphs)
head(seuratObj@meta.data)
#Clustering based on the found neighbors
seuratObj <- FindClusters(seuratObj,resolution=0.1,graph.name="RNA_snn") 
library(clustree)
clustree(seuratObj@meta.data,prefix="RNA_snn_res.1")
#UMAP clustering
seuratObj <- RunUMAP(seuratObj,dims=1:4,n.neighbors=20)
library(ggthemes)
DimPlot(seuratObj, label = "TRUE", label.size = 5)+ theme_few() + NoLegend() + labs(x = "UMAP1",y = "UMAP2",)
library(export)
graph2ppt(file="UMAP_MICROSPLITplot.pptx", width=6, height=5)   
