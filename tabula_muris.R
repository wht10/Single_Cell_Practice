Brain_Myeloid_Count <- read.csv("Tabula_Muris/FACS/Brain_Myeloid-counts.csv")
Brain_Non_Myeloid_Count <- read.csv("Tabula_Muris/FACS/Brain_Non-Myeloid-counts.csv")
meta <- read.csv("Tabula_Muris/annotations_facs.csv")

identical(Brain_Myeloid_Count[,1],Brain_Non_Myeloid_Count[,1])

rownames(Brain_Myeloid_Count) <- Brain_Myeloid_Count[,1]
rownames(Brain_Non_Myeloid_Count) <- Brain_Non_Myeloid_Count[,1]
Brain_Myeloid_Count <- Brain_Myeloid_Count[,-1]
Brain_Non_Myeloid_Count <- Brain_Non_Myeloid_Count[,-1]

library(tidyverse)
brain_meta <- meta %>% filter(tissue == "Brain_Myeloid" | tissue == "Brain_Non-Myeloid")
rownames(brain_meta) <- brain_meta$cell
i <- which(brain_meta$tissue == "Brain_Myeloid")
y <- which(brain_meta$tissue == "Brain_Non-Myeloid")
mye.ann <- Brain_Myeloid_Count[,rownames(brain_meta)[i]]
nonmye.ann <- Brain_Non_Myeloid_Count[,rownames(brain_meta)[y]]

downsamp_Myeloid <- mye.ann[,sample(colnames(mye.ann[,1:ncol(mye.ann)]), size = 600, replace = F)]

Brain_Count <- downsamp_Myeloid %>% add_column(nonmye.ann)

brain_meta_down <- brain_meta[colnames(Brain_Count),]
identical(rownames(brain_meta_down),colnames(Brain_Count))

brain_meta_down_sel <- brain_meta_down %>% select(cell_ontology_class,subtissue,tissue,mouse.id,mouse.sex)

library(Seurat)
tm.brain.down <- CreateSeuratObject(Brain_Count,meta.data = brain_meta_down_sel, project = "tm.brain.down")

tm.brain.down[["percent.mt"]] <- PercentageFeatureSet(tm.brain.down,pattern = "^Mt")
VlnPlot(tm.brain.down, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


tm.brain.down <- SCTransform(tm.brain.down, vars.to.regress = "percent.mt")

tm.brain.down <- NormalizeData(tm.brain.down, normalization.method = "LogNormalize", scale.factor = 10000)

tm.brain.down <- RunPCA(tm.brain.down)

ElbowPlot(tm.brain.down, ndims = 50)

tm.brain.down <- FindNeighbors(tm.brain.down, dims = 1:6)

tm.brain.down <- tm.brain.down %>% RunUMAP(dims = 1:6)

tm.brain.down <- FindClusters(tm.brain.down,resolution = .24)

Idents(tm.brain.down) <- tm.brain.down$cell_ontology_class

DimPlot(tm.brain.down, label = T, repel = T)

VlnPlot(tm.brain.down, features = "Camkk2")

neuron.de <- FindMarkers(tm.brain.down, ident.1 = "neuron", features = "Camkk2")
neuron.de
neuron.micro.de <- FindMarkers(tm.brain.down, ident.1 = "neuron", ident.2 = "microglial cell", features = "Camkk2")
neuron.micro.de

microglia.de <- FindMarkers(tm.brain.down, ident.1 = "microglial cell", features = "Camkk2", logfc.threshold = 0)
microglia.de

FeaturePlot(tm.brain.down, features = "Camkk2", min.cutoff =  "q10", max.cutoff = "q90")

counts <- GetAssayData(tm.brain.down, slot = "data", assay = "SCT")

camkk2.counts <- counts["Camkk2",]

library(kableExtra)

camkk2 <- data.frame(counts = camkk2.counts, celltype = tm.brain.down$cell_ontology_class)
camkk2 %>% group_by(celltype) %>% summarise(avg = mean(counts),
                                            median = median(counts), 
                                            med_over_zero = median(counts != 0),
                                            percent_zero = mean(counts == 0),
                                            percent_over_zero = mean(counts != 0)) %>% kable(
                                              col.names = c("Celltype", "Average Expression",
                                                            "Median Expression", "Median Expression of Non-Zero Values",
                                                            "% of Cells with Zero Expression", 
                                                            "% of Cells with non-Zero Expression")) %>% kable_styling()

saveRDS(tm.brain.down, file = "TabulaMuris_Brain_FACS_myeloiddown.Robj")

