#packages I used in the following code.
library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)

#folder "I-A10" is a folder that contains an output from the cell ranger count function called "filtered_feature_bc_matrix.
#I used filtered over the raw matrix because in 10X documentation it said that raw matrix includes "all valid barcodes from GEMs (Gel Bead-In EMulsions) captured in the data". 
#"because most GEMs do not actually contain cells, it follows that most barcodes in the data do not correspond to cells, but rather background noise (e.g. GEMs with free-floating mRNA from lysed or dead cells)."
#https://kb.10xgenomics.com/hc/en-us/articles/360001892491-What-is-the-difference-between-the-filtered-and-raw-gene-barcode-matrix-
IA10.data <- Read10X(data.dir = "I-A10/") #matrix of data. genes are the rows and cells are the cols

#creates a seurat object. hyperparameter are following the Seurat tutorial
IA10 <- CreateSeuratObject(counts = IA10.data, project = "I-A10", min.cells = 3, min.features = 200)
#none of the cells in the dataset had mitochondrial DNA. Cells with large amounts of mitochondrial DNA are usually low quality/ dead cells
IA10[["percent.mt"]] <- PercentageFeatureSet(IA10, pattern = "^MT-")

#visulaization of quality control feautures
VlnPlot(IA10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
IA10plot1 <- FeatureScatter(IA10, feature1 = "nCount_RNA", feature2 = "percent.mt")
IA10plot2 <- FeatureScatter(IA10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
IA10plot1 + IA10plot2

#subsetting the Seurat object to get rid of low quality cells. We filter cells that have unique feature counts over 2,500 or less than 200. 
IA10 <- subset(IA10, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Using the default parameters from the Seurat tutorial "y default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result."
IA10 <- NormalizeData(IA10, normalization.method = "LogNormalize", scale.factor = 10000)

#finds 2000 feautures that exihibit the highest cell-to-cell variation
IA10 <- FindVariableFeatures(IA10, selection.method = "vst", nfeatures = 2000)
#visulizations of the variable features
IA10top10 <- head(VariableFeatures(IA10), 10)
IA10plot3 <- VariableFeaturePlot(IA10)
IA10plot4 <- LabelPoints(plot = IA10plot3, points = IA10top10, repel = TRUE) 

#scaling the data - Shifts the expression of each gene, so that the mean expression across cells is 0 + Scales the expression of each gene, so that the variance across cells is 1
IA10 <- ScaleData(IA10, features = rownames(IA10))

#runs a pca on the Seurat object. 
IA10 <- RunPCA(IA10, features = VariableFeatures(object = IA10))
#vizualize the PCA result
VizDimLoadings(IA10, dims = 1:2, reduction = "pca")
DimPlot(IA10, reduction = "pca") 
DimHeatmap(IA10, dims = 1, cells = 500,balanced = TRUE)
DimHeatmap(IA10, dims = 1:15, cells = 500, balanced = TRUE)

#determining the dimesionailty of the data
#jack straw plot identifies PCs with a strong enrichment of feautures (genes) with low p-values
IA10 <- JackStraw(IA10, num.replicate = 100, dims = 20)
IA10 <- ScoreJackStraw(IA10, dims = 1:20)
JackStrawPlot(IA10, dims = 1:20)
#other way is to use an elbow plot which ranks the PCs based on the percentage of variance explained by each one. Where you see an "elbow" is the heuristic to set number of clusters
ElbowPlot(IA10)
#I ended up using 10 because that is where there was approximatley an "elbow" in the elbow plot and was the default paramenter used in many of the later downstream functions

#Clustering the Cells the method seurat uses " embesd cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.
#contructing the KNN graph
IA10 <- FindNeighbors(IA10, dims = 1:10)
#implements the "louvian algorithm" to group cells together. the resolution parameter "sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters." 
#the default resolution value was 0.8 so that is what I used for this dataset. Using this resolution outputs 11 "communities" or clusters  
IA10 <- FindClusters(IA10) 

IA10 <- RunUMAP(IA10, dims = 1:10)
DimPlot(IA10, reduction = "umap")

#finding differentially expressed genes within this sample
#parameters : min.pct only tests genes that in 25% of either the clusters. used value from seurat tutorial. logfc.threshold : diplays genes that have that minimum threshold. default value is 0.25
IA10.markers <- FindAllMarkers(IA10, min.pct = 0.25, logfc.threshold = 0.25)

#pathway analysis using enrichR
library(enrichR)
#dataset I decided to use was the Azimuth Cell Type Dataset. Was accurate in predicting the cell types in the PBMC Seurat tutorial.And from the Variable Feature Plot it looked like a lot of the most differential expressed genes were canoncial cell type markers for blood cells so I though the sample was from blood. 
#Was not sure if it was peripheral blood or whole blood, but Azimuth dataset in enrichR was the only one giving my statistically sig. results
dbs <- listEnrichrDbs()
dbs <- c("Azimuth_Cell_Types_2021")

#creating the list that will hold the values of the GSEA (Gene Set Enrichment Analysis)
IA10.clusters.enriched <- vector(mode = "list", length = 12)

#loop that would go cluster by cluster to do the "enrichr" function and do the GSEA. Would imput genes positively expressed >0.25 and have p_val_adj < 0.05. 
for(i in 0:10) {
  IA10.clusters.enriched[i+1] <- list(IA10.markers %>%
                                       filter(cluster == i & avg_log2FC >= 0.25 & p_val_adj < 0.05))
  IA10.clusters.enriched[i+1] <- enrichr(as.data.frame(IA10.clusters.enriched[i+1])$gene, dbs)
}

#I then manually investigated each cluster to see which was the most significant according to the results stored in the "IA10.clusters.enriched list. I inputted those values into the list below
IA10clusternames <- c("CD4+ Naive T", "CD4+ T Cell", "CD8 T", "TYROBP+ CD74+ Layer 1-6 Microglia", "CD4+ Central Memory T 3", "CD16+ Monocyte", "CD16+ Monocyte", "Natural Killer T", "B Cell", "CD16+ Monocyte", "Proliferating NK/T")

#setting the cluster numbers to their celltypes according to the Azimuth dataset.
names(IA10clusternames) <- levels(IA10) #alligns the cluster name in the "IA10clusternames list with the cluster numbers.
IA10 <- RenameIdents(IA10, IA10clusternames) #identities are now based on the "IA10clusternames list"
DimPlot(IA10, reduction = "umap", label = TRUE, pt.size = 0.5) #new UMAP clustering with celltypes


#repeat the same process for GO enrichment.
dbs <- c("GO_Biological_Process_2021")
IA10.clusters.GOenriched <- vector(mode = "list", length = 11)
for(i in 0:10) {
  IA10.clusters.GOenriched[i+1] <- list(IA10.markers %>%
                                         filter(cluster == i & avg_log2FC >= 0.25 & p_val_adj < 0.05))
  IA10.clusters.GOenriched[i+1] <- enrichr(as.data.frame(IA10.clusters.GOenriched[i+1])$gene, dbs)
}
IA10GOclusternames <- c("SRP-dependent cotranslational protein targeting to membrane", "SRP-dependent cotranslational protein targeting to membrane", "cytoplasmic translation", "regulation of mRNA stability", "regulation of mRNA splicing, via spliceosome", "positive regulation of transcription, DNA-templated", "neutrophil degranulation", "positive regulation of mRNA catabolic process", "cotranslational protein targeting to membrane", "Fc receptor mediated stimulatory signaling pathway", "mRNA processing")
names(IA10GOclusternames) <- levels(IA10) 
IA10 <- RenameIdents(IA10, IA10GOclusternames) 
DimPlot(IA10, reduction = "umap", label = TRUE, pt.size = 0.5)





