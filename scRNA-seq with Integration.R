library(Seurat)
library(patchwork)
library(tidyverse)
library(dplyr)
library(clustree)
library(gridExtra)
library(SingleCellExperiment)
library(SingleR)
library(celldex)


Data <- Read10X("C:/Users/Joy/Desktop/SingleCell/New")

sdf <- CreateSeuratObject(Data, project = "PBMC")
head(sdf@meta.data)

head(colnames(sdf))

# map suffix to sample name (from the GEO series description)
suffix_map <- c(
  "1"="nCoV1","2"="nCoV2","3"="Flu1","4"="Flu2","5"="Normal1",
  "6"="Flu3","7"="Flu4","8"="Flu5","9"="nCoV3","10"="nCoV4",
  "11"="nCoV5","12"="nCoV6","13"="Normal2","14"="Normal3","15"="nCoV7",
  "16"="nCoV8","17"="nCoV9","18"="nCoV10","19"="Normal4","20"="nCoV11"
)


barcodes <- colnames(sdf)
suffixes <- sub(".*-(\\d+)$", "\\1", barcodes)

# create mapped sample names
sample_labels <- suffix_map[suffixes]

# IMPORTANT: assign barcode names so Seurat can match them
names(sample_labels) <- barcodes

# add metadata safely
sdf <- AddMetaData(sdf, metadata = sample_labels, col.name = "sample")
View(sdf@meta.data)


# define condition from sample strings
sdf$condition <- ifelse(grepl("^nCoV", sdf$sample), "COVID",
                        ifelse(grepl("^Normal", sdf$sample), "Healthy", "Flu"))

table(sdf$sample, sdf$condition)
View(sdf@meta.data)


samples_to_keep <- c("nCoV1", "nCoV2", "Flu1", "Flu2", "Normal1", "Normal2")

sdf_subset <- subset(sdf, subset = sample %in% samples_to_keep)

table(sdf_subset$sample)
table(sdf_subset$condition)
head(sdf_subset@meta.data)

# QC and filtering per sample
sdf_subset[["percent.mt"]] <- PercentageFeatureSet(sdf_subset, pattern = "^MT[-\\.]")

VlnPlot(sdf, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)

FeatureScatter(sdf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

sdf_subset


# Filtering
sdf_subset <- subset(sdf_subset, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
sdf_subset

#Normalizing
sdf_subset <- SCTransform(sdf_subset,
                   vars.to.regress = "percent.mt",
                   variable.features.n = 3000,
                   vst.flavor = "v2",
                   verbose = T)

head(sdf_subset@meta.data)

sdf_subset <- RunPCA(sdf_subset, assay = "SCT", verbose = T)
ElbowPlot(sdf_subset)

sdf_subset <- RunTSNE(sdf_subset, reduction = "pca", dims = 1:15)
TSNEPlot(sdf_subset)

sdf_subset <- FindNeighbors(sdf_subset, dims = 1:15)
sdf_subset <- FindClusters(sdf_subset, resolution = 0.4)

plot1 <- DimPlot(sdf_subset, reduction = "tsne", label = T, pt.size = 0.5, 
        group.by = "sample")+
  ggtitle("Before Integration for sample")

plot1

plot2 <- DimPlot(sdf_subset, reduction = "tsne", label = T, pt.size = 0.5, 
                 group.by = "condition")+
  ggtitle("Before Integration for condition")


plot1 + plot2


#############################################################
#-----perform integration to correct for batch effects------#
#############################################################

# split by sample
seurat_list <- SplitObject(sdf_subset, split.by = "sample")

# run SCTransform on each sample (regress percent.mt as before)
for (i in seq_along(seurat_list)) {
  seurat_list[[i]] <- SCTransform(
    seurat_list[[i]],
    vars.to.regress = "percent.mt",
    variable.features.n = 3000,
    vst.flavor = "v2",
    verbose = FALSE
  )
}

#Select features to use for integration
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

#Prep for SCT integration
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features, 
                                  verbose = T)


#Find anchors using CCA
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "cca",   # use CCA for anchor finding as you requested
  dims = 1:15
)

#ntegrate data (SCT)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                            dims = 1:15)


#Run PCA on the integrated assay and check elbow
# make sure we use the integrated assay for downstream steps
DefaultAssay(integrated) <- "integrated"

integrated <- RunPCA(integrated, verbose = T)
ElbowPlot(integrated)  # pick how many PCs to use (often 1:15 or 1:30)

#Run t-SNE (using PCA dims) and plot
# choose dims based on ElbowPlot, e.g. dims = 1:15
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:15, 
                      check_duplicates = FALSE)

DimPlot(integrated, reduction = "tsne", group.by = "sample", pt.size = 0.5) + 
  ggtitle("By sample (integrated)")


#Neighbors and clustering
integrated <- FindNeighbors(integrated, dims = 1:15)
integrated <- FindClusters(integrated, resolution = 0.4)

head(integrated@meta.data)

plot3 <- DimPlot(integrated, reduction = "tsne", label = T, pt.size = 0.5,
                 group.by = "sample")+
  ggtitle("After Integration for Sample")


plot4 <- DimPlot(integrated, reduction = "tsne", label = T, pt.size = 0.5,
                 group.by = "condition")+
  ggtitle("After Integration for condition")


plot1 + plot3
plot2 + plot4


#Prepare: set assay and identity
# use SCT for marker finding (integrated assay was used for PCA/clustering)
DefaultAssay(integrated) <- "SCT"
Idents(integrated) <- "seurat_clusters"   # make sure identities are cluster IDs

#Find marker genes for every cluster
# Prepare the SCT assay for marker detection
integrated <- PrepSCTFindMarkers(integrated)


# Now run FindAllMarkers
markers <- FindAllMarkers(integrated,
                          only.pos = FALSE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25,
                          test.use = "wilcox")


head(markers)

# save markers list for review
write.csv(markers, file = "all_clusters_markers.csv", row.names = FALSE)


# get top 5 markers per cluster
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# save top5 markers list for review
write.csv(top5, file = "top5_markers_per_cluster.csv", row.names = FALSE)


# view a few gene plots on the tSNE
FeaturePlot(integrated, features = c("CD3D","MS4A1","CD14"), reduction = "tsne",
            ncol = 3)
VlnPlot(integrated, features = c("CD3D","MS4A1","CD14"), group.by = "condition", 
        pt.size = 0.5)


#Cell Annotation
# prepare reference and test
ref <- HumanPrimaryCellAtlasData()   # good general reference for human blood
sce <- as.SingleCellExperiment(integrated)

# per-cell prediction
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
integrated$Annotation <- pred$labels

head(integrated)


##Rename clustering
## show cluster levels (current identities)
levels(Idents(integrated))

# contingency table: cluster x SingleR label
table(integrated$seurat_clusters, integrated$Annotation)

# Find dominant cell type per cluster
cluster_annotation <- as.data.frame(table(integrated$seurat_clusters, integrated$Annotation)) %>%
  group_by(Var1) %>%
  slice_max(Freq) %>%
  dplyr::select(Cluster = Var1, Annotation = Var2)


# Create a vector with new labels
new_cluster_ids <- cluster_annotation$Annotation
names(new_cluster_ids) <- cluster_annotation$Cluster

# Rename clusters
integrated <- RenameIdents(integrated, new_cluster_ids)

# Check
levels(integrated)

# Visualize
DimPlot(integrated, reduction = "tsne", label = TRUE, pt.size = 0.5)

# saving the result
saveRDS(integrated, file = "C:/Users/Joy/Desktop/SingleCell/New/integrated_data.rds")