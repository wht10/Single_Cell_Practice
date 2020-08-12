#####
# Standard Workflow

devtools::install_github('satijalab/seurat-data')
library(SeuratData)
AvailableData()

InstallData("panc8")

data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}


reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")]

for (i in 1:length(reference.list)){
  all.genes <- rownames(reference.list[[i]])
  reference.list[[i]] <- ScaleData(reference.list[[i]], features = all.genes)
}

for (i in 1:length(reference.list)){
  reference.list[[i]] <- RunPCA(reference.list[[i]], features = VariableFeatures(object = reference.list[[i]]))
}


ElbowPlot(reference.list[[1]])
ElbowPlot(reference.list[[2]])

pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)


library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
p1 + p2

pancreas.query <- pancreas.list[["fluidigmc1"]]
pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype, 
                            dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)

pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)

table(pancreas.query$predicted.id)

VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")

#####
# SCTransform 

options(future.globals.maxSize = 4000 * 1024^2)

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = FALSE)
}

pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = FALSE)

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
plots <- DimPlot(pancreas.integrated, group.by = c("tech", "celltype"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                  override.aes = list(size = 3)))

#####
#Reference Based


InstallData("pbmcsca")
data("pbmcsca")
pbmc.list <- SplitObject(pbmcsca, split.by = "Method")
for (i in names(pbmc.list)) {
  pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
}
pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset <- which(names(pbmc.list) == "10x Chromium (v3)")

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", 
                                       anchor.features = pbmc.features, reference = reference_dataset)
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT")

pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(object = pbmc.integrated, dims = 1:30)
plots <- DimPlot(pbmc.integrated, group.by = c("Method", "CellType"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 4, byrow = TRUE, 
                                                                     override.aes = list(size = 2.5)))