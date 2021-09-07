#packages I used in the following code.
library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)

#folder "N-A9 is a folder that contains an output from the cell ranger count function called "filtered_feature_bc_matrix.
#I used filtered over the raw matrix because in 10X documentation it said that raw matrix includes "all valid barcodes from GEMs (Gel Bead-In EMulsions) captured in the data". 
#"because most GEMs do not actually contain cells, it follows that most barcodes in the data do not correspond to cells, but rather background noise (e.g. GEMs with free-floating mRNA from lysed or dead cells)."
#https://kb.10xgenomics.com/hc/en-us/articles/360001892491-What-is-the-difference-between-the-filtered-and-raw-gene-barcode-matrix-
NA9.data <- Read10X(data.dir = "N-A9/") #matrix of data. genes are the rows and cells are the cols

#creates a seurat object. hyperparameter are following the Seurat tutorial
NA9 <- CreateSeuratObject(counts = NA9.data,project = "NA-9", min.cells = 3, min.features = 200)
#none of the cells in the dataset had mitochondrial DNA. Cells with large amounts of mitochondrial DNA are usually low quality/ dead cells
NA9[["percent.mt"]] <- PercentageFeatureSet(NA9, pattern = "^MT-")

#visulaization of quality control feautures
NA9plot1 <- VlnPlot(NA9, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
NA9plot2 <- FeatureScatter(NA9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#subsetting the Seurat object to get rid of low quality cells. We filter cells that have unique feature counts over 2,500 or less than 200. 
NA9 <- subset(NA9, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Using the default parameters from the Seurat tutorial "y default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result."
NA9 <- NormalizeData(NA9, normalization.method = "LogNormalize", scale.factor = 10000)

#finds 2000 feautures that exihibit the highest cell-to-cell variation
NA9 <- FindVariableFeatures(NA9, selection.method = "vst", nfeatures = 2000)
#visulizations of the variable features
NA9top10 <- head(VariableFeatures(NA9), 10) #list of the top10 most highly variable genes
NA9plot3 <- VariableFeaturePlot(NA9)
NA9plot4 <- LabelPoints(plot = NA9plot3, points = NA9top10, repel = TRUE) 

#scaling the data - Shifts the expression of each gene, so that the mean expression across cells is 0 + Scales the expression of each gene, so that the variance across cells is 1
NA9 <- ScaleData(NA9, features = rownames(NA9))

#runs a pca on the Seurat object. 
NA9 <- RunPCA(NA9, features = VariableFeatures(object = NA9))
#vizualize the PCA result
VizDimLoadings(NA9, dims = 1:2, reduction = "pca")
DimPlot(NA9, reduction = "pca") 
DimHeatmap(NA9, dims = 1, cells = 500,balanced = TRUE)
DimHeatmap(NA9, dims = 1:15, cells = 500, balanced = TRUE)

#determining the dimesionailty of the data
#jack straw plot identifies PCs with a strong enrichment of feautures (genes) with low p-values
NA9 <- JackStraw(NA9, num.replicate = 100, dims = 20)
NA9 <- ScoreJackStraw(NA9, dims = 1:20)
JackStrawPlot(NA9, dims = 1:20)
#other way is to use an elbow plot which ranks the PCs based on the percentage of variance explained by each one. Where you see an "elbow" is the heuristic to set number of clusters
ElbowPlot(NA9)
#I ended up using 10 because that is where there was approximatley an "elbow" in the elbow plot and was the default paramenter used in many of the later downstream functions

#Clustering the Cells the method seurat uses " embesd cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.
#contructing the KNN graph
NA9 <- FindNeighbors(NA9, dims = 1:10)
#implements the "louvian algorithm" to group cells together. the resolution parameter "sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters." 
#the default resolution value was 0.8 so that is what I used for this dataset. Using this resolution outputs 12 "communities" or clusters  
NA9 <- FindClusters(NA9) 

NA9 <- RunUMAP(NA9, dims = 1:10)
NA9plot5 <- DimPlot(NA9, reduction = "umap")

#finding differentially expressed genes within this sample
#parameters : min.pct only tests genes that in 25% of either the clusters. used value from seurat tutorial. logfc.threshold : diplays genes that have that minimum threshold. default value is 0.25
NA9.markers <- FindAllMarkers(NA9, min.pct = 0.25, logfc.threshold = 0.25)

#pathway analysis using enrichR
library(enrichR)
#dataset I decided to use was the Azimuth Cell Type Dataset. Was accurate in predicting the cell types in the PBMC Seurat tutorial.And from the Variable Feature Plot it looked like a lot of the most differential expressed genes were canoncial cell type markers for blood cells so I though the sample was from blood. 
#Was not sure if it was peripheral blood or whole blood, but Azimuth dataset in enrichR was the only one giving my statistically sig. results
dbs <- listEnrichrDbs()
dbs <- c("Azimuth_Cell_Types_2021")

#creating the list that will hold the values of the GSEA (Gene Set Enrichment Analysis)
NA9.clusters.enriched <- vector(mode = "list", length = 12)

#loop that would go cluster by cluster to do the "enrichr" function and do the GSEA. Would imput genes positively expressed >0.25 and have p_val_adj < 0.05. 
for(i in 0:11) {
  NA9.clusters.enriched[i+1] <- list(NA9.markers %>%
                                   filter(cluster == i & avg_log2FC >= 0.25 & p_val_adj < 0.05))
  NA9.clusters.enriched[i+1] <- enrichr(as.data.frame(NA9.clusters.enriched[i+1])$gene, dbs)
}

#I then manually investigated each cluster to see which was the most significant according to the results stored in the "NA9.clusters.enriched list. I inputted those values into the list below
NA9clusternames <- c("CD4+ Central Memory T", "CD4+ T", "TYROBP+ CD74+ Layer 1-6 Microglia","CD4+ T", "CD56-dim Natural Killer 3", "Monocyte", "B Cell", "Monocyte", "CD16+ Monocyte", "CD16+ Monocyte","Plasmacytoid Dendritic Cell", "Proliferating NK/T" )

#setting the cluster numbers to their celltypes according to the Azimuth dataset.
names(NA9clusternames) <- levels(NA9) #alligns the cluster name in the "NA9clusternames list with the cluster numbers.
NA9 <- RenameIdents(NA9, NA9clusternames) #identities are now based on the "NA9clusternames list"
NA9plot6 <- DimPlot(NA9, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) #new UMAP clustering with celltypes


#repeat the same process for GO enrichment.
dbs <- c("GO_Biological_Process_2021")
NA9.clusters.GOenriched <- vector(mode = "list", length = 12)
for(i in 0:11) {
  NA9.clusters.GOenriched[i+1] <- list(NA9.markers %>%
                                       filter(cluster == i & avg_log2FC >= 0.25 & p_val_adj < 0.05))
  NA9.clusters.GOenriched[i+1] <- enrichr(as.data.frame(NA9.clusters.GOenriched[i+1])$gene, dbs)
}
NA9GOclusternames <- c("SRP-dependent cotranslational protein targeting to membrane","SRP-dependent cotranslational protein targeting to membrane", "mRNA splicing, via spliceosome", "cotranslational protein targeting to membrane", "cytokine-mediated signaling pathway", "regulation of mRNA stability", "cytoplasmic translation", "positive regulation of mRNA catabolic process", "neutrophil degranulation", "neutrophil degranulation", "endoplasmic reticulum to cytosol transport", "DNA metabolic process")
NA9 <- FindClusters(NA9) #resetting things so it is 
names(NA9GOclusternames) <- levels(NA9) 
NA9 <- RenameIdents(NA9, NA9GOclusternames) 
NA9plot7 <- DimPlot(NA9, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)





