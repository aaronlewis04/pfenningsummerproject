library(Seurat)
library(SeuratData)
library(patchwork)

#the psuedocode for these functions and their respective parameters are in the other scripts. 
#loading the NA9 dataset and creating its Seurat object
NA9.data <- Read10X(data.dir = "N-A9/") #matrix of data. genes are the rows and cells are the cols
NA9 <- CreateSeuratObject(counts = NA9.data,project = "NA-9", min.cells = 3, min.features = 200)
NA9[["percent.mt"]] <- PercentageFeatureSet(NA9, pattern = "^MT-")
NA9 <- subset(NA9, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#loading the IA10 dataset and creating its Seurat object 
IA10.data <- Read10X(data.dir = "I-A10/") #matrix of data. genes are the rows and cells are the cols
IA10 <- CreateSeuratObject(counts = IA10.data,project = "I-A10", min.cells = 3, min.features = 200)
IA10[["percent.mt"]] <- PercentageFeatureSet(IA10, pattern = "^MT-")
IA10 <- subset(IA10, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#merging the two datasets. Adding metadata that specifies the sample each comes from
NA9 <- AddMetaData(NA9, metadata="NA9", col.name="Sample")
IA10 <-AddMetaData(IA10, metadata="IA10", col.name="Sample")
rada <- merge(NA9, y = IA10, add.cell.ids = c("NA9", "IA10"), project = "rada") 

#following the seurat integration tutorial : https://satijalab.org/seurat/articles/integration_introduction.html 
#now splitting the dataset into a list of two seurat objects
rada.list <- SplitObject(rada, split.by = "Sample")
#normalize and identify variable features for each dataset independently
rada.list <- lapply(X = rada.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
rada.features <- SelectIntegrationFeatures(object.list = rada.list)
#Identify anchors - pairs of cells from each dataset that are contained within each other's neighborhoods (also known as mutual nearest neighbors). Used later to integrate the two datasets 
rada.anchors <- FindIntegrationAnchors(object.list = rada.list, anchor.features = rada.features)
#integrates the two datasets
rada <- IntegrateData(anchorset = rada.anchors)
#there are now two assays in the Seurat object, one with the integrated dataset and on with the original dataset. Here we specify that we will perform downstream analysis on the integrated dataset. The original dataset is in the "RNA" assay. 
DefaultAssay(rada) <- "integrated"
# Run the standard workflow for visualization and clustering. Following the parameter values of the Seurat tutorial. 
rada <- ScaleData(rada, verbose = FALSE)
rada <- RunPCA(rada, npcs = 30, verbose = FALSE)
rada <- RunUMAP(rada, reduction = "pca", dims = 1:30)
rada <- FindNeighbors(rada, reduction = "pca", dims = 1:30)
rada <- FindClusters(rada, resolution = 0.5)
#UMAP visualizations
p1 <- DimPlot(rada, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(rada, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

#enrichr GSEA analysis
library(enrichR)
#finding the markers for each cluster within the integrated data
rada.integrated.markers <- FindAllMarkers(rada, min.pct = 0.25, logfc.threshold = 0.25)
dbs <- c("Azimuth_Cell_Types_2021")
rada.clusters.enriched <- vector(mode = "list", length = 11)
library(dplyr)
for(i in 0:10) {
  rada.clusters.enriched[i+1] <- list(rada.integrated.markers %>%
                                         filter(cluster == i & avg_log2FC >= 0.25 & p_val_adj < 0.05))
  rada.clusters.enriched[i+1] <- enrichr(as.data.frame(rada.clusters.enriched[i+1])$gene, dbs)
}
rada$number.cluster <- Idents(rada) #saving the original number clusters for later analysis
#manually inputted the results into this next function
rada <- RenameIdents(rada, `0` = "CD4+ Central Memory T 3", `1` = "CD4+ Central Memory T 3", `2` = "CD4+ T Cell",
                                   `3` = "Monocyte", `4` = "TYROBP+ CD74+ Layer 1-6 Microglia", `5` = "B Cell", `6` = "CD16+ Monocyte", `7` = "CD16+ Monocyte", `8` = "Proliferating NK/T", `9` = "Plasmacytoid Dendritic Cell",
                                   `10` = "Hematopoietic Stem And Progenitor Cell")
p3 <- DimPlot(rada, label = TRUE) #dim plot with the new identities 
rada$celltype <- Idents(rada) #saving the Azimuth celltypes for later analysis

#running GO on it
dbs <- c("GO_Biological_Process_2021")
rada.clusters.GOenriched <- vector(mode = "list", length = 11)
for(i in 0:10) {
  rada.clusters.GOenriched[i+1] <- list(rada.integrated.markers %>%
                                        filter(cluster == i & avg_log2FC >= 0.25 & p_val_adj < 0.05))
  rada.clusters.GOenriched[i+1] <- enrichr(as.data.frame(rada.clusters.GOenriched[i+1])$gene, dbs)
}
rada <- FindClusters(rada, resolution = 0.5) #resets the identities to numbers. Easiest way I found to then rename the identities to GO pathways in the next step.
#manually inputted the results into this next function
rada <- RenameIdents(rada, `0` = "SRP-dependent cotranslational protein targeting to membrane", `1` = "cytokine-mediated signaling pathway", `2` = "SRP-dependent cotranslational protein targeting to membrane",
                     `3` = "regulation of mRNA stability", `4` = "mRNA splicing, via spliceosome", `5` = "cotranslational protein targeting to membrane", `6` = "neutrophil degranulation", `7` = "Fc receptor mediated stimulatory signaling pathway", `8` = "mRNA processing", `9` = "protein phosphorylation",
                     `10` = "positive regulation of transcription, DNA-templated")

rada$GOprocess <- Idents(rada)
p4 <- DimPlot(rada, label = TRUE)
#differential expression of NA9 sample compared to IA10 sample
#Seurat Github and tutorial instructs to set the assay back to the original unintegegrated results to perform differential expression analysis. 
DefaultAssay(rada) <- "RNA"
Idents(rada) <- rada$orig.ident
diffexpgenes <- FindMarkers(rada, ident.1 = "NA-9", ident.2 = "I-A10") #differenital expression on the two datasets. 106 diff. exp. genes

diffexpgenes$gene <- rownames(diffexpgenes)
library(ggrepel)
#creating a volcano plot from the differential expression results 
p5 <- ggplot(data=diffexpgenes, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(data = diffexpgenes[1:10,],aes(x = avg_log2FC, y = -log10(p_val_adj),label=gene, max.overlaps = 20)) +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggtitle("NA-9 Sample vs I-A10 Sample") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(plot.title = element_text(face="bold"))

#violin plots for the expression of Ube3a the gene of interest in the integrated dataset
p6 <- VlnPlot(rada, features = c("Ube3a"), split.by = "Sample", group.by = "number.cluster",
                 pt.size = 0, combine = FALSE)
p7 <- VlnPlot(rada, features = c("Ube3a"), split.by = "Sample", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
p8 <- VlnPlot(rada, features = c("Ube3a"), split.by = "Sample", group.by = "GOprocess",
              pt.size = 0, combine = FALSE)

