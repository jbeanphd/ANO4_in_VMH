# Jonathan Bean, PhD
# 2022-07-26

#libraries
libs <- c('RcppThread','dplyr','tidyr','ggplot2','DataCombine', 'extrafont','gdata','ggrepel','Seurat')

lapply(libs, require, character.only = TRUE)


# data downloaded from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172204
# published in eLife
# Affinati et al. 2021 
# https://pubmed.ncbi.nlm.nih.gov/34018926/
# https://elifesciences.org/articles/69065



#Import Data
expression_matrix_ms <- ReadMtx(
  mtx = "../../Mouse_GSE172204/GSE172204_matrix.mtx.gz", features = "../../Mouse_GSE172204/GSE172204_features.tsv.gz",
  cells = "../../Mouse_GSE172204/GSE172204_barcodes.tsv.gz", feature.column = 1
)
metadata <- read.csv('../../Mouse_GSE172204/GSE172204_metadata.csv.gz', row.names = 1)
seurat_object_ms <- CreateSeuratObject(counts = expression_matrix_ms, meta.data = metadata)


# filter out cells with missing data
filtered_seurat <- subset(x = seurat_object_ms, 
                          subset= (Run != 'NA'))

# normalize data, regress batch out
# significant batch effects without regressing out Run. Data not shown
SCTseurat <- SCTransform(filtered_seurat, vars.to.regress = 'Run')

# dimension reductions PCA and UMAP
SCTseurat <- RunPCA(SCTseurat, assay = 'SCT')
SCTseurat <- RunUMAP(SCTseurat, assay = 'SCT', dims = (1:30))

# corresponding to figure 1 - supplement 1 H
# different shape but labels from metadata do cluster together
UMAPPlot(SCTseurat,group.by = "Run", label = TRUE, shuffle = TRUE)
UMAPPlot(SCTseurat,group.by = "Cell_type", label = TRUE)

print(x = SCTseurat[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# set active ident to cell type
filtered_seurat@active.ident <- factor(filtered_seurat$Cell_type)

# filter for neurons
neurons <- subset(filtered_seurat, idents = "Neuron")


# save objects for future use
saveRDS(seurat_object_ms, 'data/seurat_object_ms.rds')
saveRDS(filtered_seurat, 'data/filtered_seurat.rds')
saveRDS(neurons, 'data/neurons.rds')

# load object of interest
neurons <- readRDS('data/neurons.rds')
neurons@active.ident <- factor(neurons$Neuron_cluster)

neurons <- FindVariableFeatures(neurons, assay = 'RNA')
neurons <- SCTransform(neurons, conserve.memory = FALSE, assay = 'RNA', vars.to.regress = 'Run')
var_genes <- intersect(
  neurons@assays$RNA@var.features, neurons@assays$SCT@var.features
)


neurons <- RunPCA(neurons, assay = 'SCT')
ElbowPlot(neurons, ndims = 50)
# I don't see an asymptote, so I'll use all available dimensions 
neurons <- RunUMAP(neurons, assay = 'SCT', dims = 1:50)

UMAPPlot(neurons, group.by = "Run", label = TRUE)
UMAPPlot(neurons, group.by = "Neuron_cluster", label = TRUE)



# compare RNA and SCT assays for batch effects
DefaultAssay(neurons) <- "RNA"
neurons <- ScaleData(neurons, features = var_genes, vars.to.regress = "Run")
neurons <- RunPCA(neurons, npcs = 100)

# custome function from
# https://github.com/alanrupp/affinati-elife-2021
knee_test <- function(object) {
  n_pc = ncol(object@reductions$pca@cell.embeddings)
  total_var <- sum(object@reductions$pca@stdev^2)
  percent_var <- cumsum(object@reductions$pca@stdev^2)/total_var * n_pc
  diminishing <- which(percent_var - lag(percent_var) < 1)
  return(min(diminishing) - 1)
}

neurons@reductions$pca@misc$sig_pcs <- knee_test(neurons)
neurons <- RunUMAP(neurons, dims = 1:neurons@reductions$pca@misc$sig_pcs)
umap_object <- neurons[["umap"]]
pca_object <- neurons[["pca"]]
UMAPPlot(neurons, group.by = "Run", label = TRUE, shuffle = TRUE)
# compare RNA and SCT assays for batch effects
DefaultAssay(neurons) <- "SCT"
neurons <- RunPCA(neurons, npcs = 100)
neurons@reductions$pca@misc$sig_pcs <- knee_test(neurons)
neurons <- RunUMAP(neurons, dims = 1:neurons@reductions$pca@misc$sig_pcs)
UMAPPlot(neurons, group.by = "Run", label = TRUE, shuffle = TRUE)

# corresponding to figure 1B
UMAPPlot(neurons, group.by = "Neuron_cluster", label = TRUE, shuffle = TRUE)


#filtered_seurat <- readRDS('data/filtered_seurat.rds')

#subset vmh neurons
vmh <- subset(filtered_seurat, VMH_cluster != "")

DefaultAssay(vmh) <- "RNA"
vmh <- FindVariableFeatures(vmh)
vmh <- SCTransform(vmh, conserve.memory = FALSE, vars.to.regress = 'Run')
vmh <- FindVariableFeatures(vmh)

var_genes <- intersect(
  vmh@assays$RNA@var.features, vmh@assays$SCT@var.features
)
# compare RNA and SCT assays for batch effects, none provided Run is regressed out. 
DefaultAssay(vmh) <- "RNA"
vmh <- ScaleData(vmh, features = var_genes, vars.to.regress = "Run")
vmh <- RunPCA(vmh, npcs = 100)
vmh@reductions$pca@misc$sig_pcs <- knee_test(vmh)
ElbowPlot(vmh, ndims = 100)
vmh <- RunUMAP(vmh, dims = 1:vmh@reductions$pca@misc$sig_pcs)


DimPlot(vmh, group = 'Run', label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE)

DimPlot(vmh, group = 'VMH_cluster', label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) +
  theme_void() +
  theme(legend.position = "none")
#ggsave(filename = 'figures/VMH_clusters.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)

umap_object <- vmh[["umap"]]
pca_object <- vmh[["pca"]]
# compare RNA and SCT assays for batch effects
DefaultAssay(vmh) <- "SCT"
vmh <- RunPCA(vmh, npcs = 100)
vmh@reductions$pca@misc$sig_pcs <- knee_test(vmh)
vmh <- RunUMAP(vmh, dims = 1:vmh@reductions$pca@misc$sig_pcs)


DimPlot(vmh, group = 'Run', label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE)

# corresponding to figure 2A
DimPlot(vmh, group = 'VMH_cluster', label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) +
  theme_void() +
  theme(legend.position = "none")




###
DefaultAssay(vmh) <- "RNA"
#vmh <- ScaleData(vmh, features = var_genes, vars.to.regress = "Run")
vmh <- RunPCA(vmh, npcs = 100)
vmh@reductions$pca@misc$sig_pcs <- knee_test(vmh)
ElbowPlot(vmh, ndims = 100)
vmh <- RunUMAP(vmh, dims = 1:vmh@reductions$pca@misc$sig_pcs)

DimPlot(vmh, group = 'Run', label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE)

DimPlot(vmh, group = 'VMH_cluster', label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) 
#  theme_void() +
#  theme(legend.position = "none")


DefaultAssay(vmh) <- "SCT"

FeaturePlot(vmh, "Ano4", pt.size = 0.5) +
  scale_color_gradient(name = expression(italic("Ano4")),
                       low = "grey", high = "red") +
  theme_void() +
  theme(legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12))
#ggsave(filename = 'figures/Ano4_in_VMH_clusters_SCT2.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)



FeaturePlot(vmh, "Gck", pt.size = 0.5) +
  scale_color_gradient(name = expression(italic("Gck")),
                       low = "grey", high = "red") +
  theme_void() +
  theme(legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12))
#ggsave(filename = 'figures/Gck_in_VMH_clusters_SCT2.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)


FeaturePlot(vmh, "Adcyap1", pt.size = 0.5) +
  scale_color_gradient(name = expression(italic("Adcyap1")),
                       low = "grey", high = "red") +
  theme_void() +
  theme(legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12))
#ggsave(filename = 'figures/Adcyap1_in_VMH_clusters_SCT2.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)



FeaturePlot(vmh, "Nos1", pt.size = 0.5) +
  scale_color_gradient(name = expression(italic("Nos1")),
                       low = "grey", high = "red") +
  theme_void() +
  theme(legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12))
#ggsave(filename = 'figures/Nos1_in_VMH_clusters_SCT2.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)


FeaturePlot(vmh, "Esr1", pt.size = 0.5) +
  scale_color_gradient(name = expression(italic("Esr1")),
                       low = "grey", high = "red") +
  theme_void() +
  theme(legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12))
#ggsave(filename = 'figures/Esr1_in_VMH_clusters_SCT2.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)


FeaturePlot(vmh, "Cckbr", pt.size = 0.5) +
  scale_color_gradient(name = expression(italic("Cckbr")),
                       low = "grey", high = "red") +
  theme_void() +
  theme(legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12))
#ggsave(filename = 'figures/Cckbr_in_VMH_clusters_SCT2.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)



Ano4 <- subset(filtered_seurat, VMH_cluster != "")
Ano4 <-  subset(Ano4, subset = Ano4 > 1)

Ano4 <- FindVariableFeatures(Ano4)
Ano4 <- ScaleData(Ano4, features = Ano4@assays$RNA@var.features, vars.to.regress = "Run")
Ano4 <- RunPCA(Ano4, npcs = 100)
Ano4@reductions$pca@misc$sig_pcs <- knee_test(vmh)
ElbowPlot(Ano4, ndims = 100)
Ano4 <- RunUMAP(Ano4, dims = 1:Ano4@reductions$pca@misc$sig_pcs)





UMAPPlot(Ano4, shuffle = TRUE, group.by = 'VMH_cluster', label = TRUE)


DimPlot(Ano4, group = 'VMH_cluster', label = TRUE, label.size = 4, pt.size = 0.5, repel = TRUE) +
  theme_void() +
  theme(legend.position = "none")
#ggsave(filename = 'figures/Ano4_clusters.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)


FeaturePlot(Ano4, "Ano4", pt.size = 0.5, max.cutoff = 10) +
  scale_color_gradient(name = expression(italic("Ano4")),
                       low = "grey", high = "red") +
  theme_void() +
  theme(legend.key.width = unit(0.1, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12))


RidgePlot(vmh, features = 'Ano4', group.by = 'VMH_cluster', assay = 'SCT') +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 12))

#ggsave(filename = 'figures/Ano4_in_VMH_clusters_RP_SCT.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)




saveRDS(Ano4, file = 'data/Ano4.rds')

Ano4 <- readRDS('data/Ano4.rds')

Ano4_meta <- Ano4@meta.data
Ano4_meta <- Ano4_meta %>% 
  mutate(Ano4_cluster = ifelse(
    VMH_cluster == 'Dlk1_1', 'Dlk1_1,2',
    ifelse(
      VMH_cluster == 'Dlk1_2', 'Dlk1_1,2',
      ifelse(
        VMH_cluster == 'Esr1_1', 'Esr1_1,2',
        ifelse(
          VMH_cluster == 'Esr1_2', 'Esr1_1,2',
          ifelse(
            VMH_cluster == 'Esr1_3', 'Esr1_3,4',
            ifelse(
              VMH_cluster == 'Esr1_4', 'Esr1_3,4',
              ifelse(
                VMH_cluster == 'Fezf1_1', 'Fezf1_1:3',
                ifelse(
                  VMH_cluster == 'Fezf1_2', 'Fezf1_1:3',
                  ifelse(
                    VMH_cluster == 'Fezf1_3', 'Fezf1_1:3',
                    ifelse(
                      VMH_cluster == 'Foxp2_1', 'Foxp2_1,2',
                      ifelse(
                        VMH_cluster == 'Foxp2_2', 'Foxp2_1,2',
                        ifelse(
                          VMH_cluster == 'Foxp2_3', 'Foxp2_3',
                          ifelse(
                            VMH_cluster == 'Lepr_1', 'Lepr_1,3',
                            ifelse(
                              VMH_cluster == 'Lepr_3', 'Lepr_1,3',
                              ifelse(
                                VMH_cluster == 'Lepr_2', 'Lepr_2',
                                ifelse(
                                  VMH_cluster == 'Lepr_4', 'Lepr_4,5',
                                  ifelse(
                                    VMH_cluster == 'Lepr_5', 'Lepr_4,5',
                                    ifelse(
                                      VMH_cluster == 'Lepr_6', 'Lepr_6,9',
                                      ifelse(
                                        VMH_cluster == 'Lepr_9', 'Lepr_6,9',
                                        ifelse(
                                          VMH_cluster == 'Lepr_7', 'Lepr_7,8',
                                          ifelse(
                                            VMH_cluster == 'Lepr_8', 'Lepr_7,8',
                                            ifelse(
                                              VMH_cluster == 'Nfib_1', 'Nfib_1:3',
                                              ifelse(
                                                VMH_cluster == 'Nfib_2', 'Nfib_1:3',
                                                ifelse(
                                                  VMH_cluster == 'Nfib_3', 'Nfib_1:3', 'Missed')))))))))))))))))))))))))

Ano4@meta.data <- Ano4_meta                                                                                 

Ano4 <- SCTransform(Ano4, vars.to.regress = 'Run')

UMAPPlot(Ano4, group.by = 'Ano4_cluster')
DotPlot(Ano4, features = c('Ano4','Gck','Adcyap1','Nos1','Esr1','Cckbr'), group.by = 'Ano4_cluster', assay = 'SCT') +
  theme_classic() +
  theme(axis.title = element_text(size = 10, colour = 'black'),
        legend.text = element_text(size = 8, colour = 'black'),
        legend.title = element_text(size = 8, color = 'black', angle = 90, vjust = 1.5),
        axis.text = element_text(size = 8, colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, colour = 'black'))
#ggsave(filename = 'figures/Ano4_dotplot_SCT.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)



DotPlot(vmh, features = c('Ano4','Gck','Adcyap1','Nos1','Cckbr'), group.by = 'VMH_cluster', assay = 'SCT') +
  theme_classic() +
  theme(axis.title = element_text(size = 10, colour = 'black'),
        legend.text = element_text(size = 8, colour = 'black'),
        legend.title = element_text(size = 8, color = 'black', angle = 90, vjust = 1.5),
        axis.text = element_text(size = 8, colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, colour = 'black'))
#ggsave(filename = 'figures/VMH_neuron_dotplot_SCT2.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)


DotPlot(vmh, features = c('Ano4','Nr5a1'), group.by = 'VMH_cluster', assay = 'SCT') +
  theme_classic() +
  theme(axis.title = element_text(size = 10, colour = 'black'),
        legend.text = element_text(size = 8, colour = 'black'),
        legend.title = element_text(size = 8, color = 'black', angle = 90, vjust = 1.5),
        axis.text = element_text(size = 8, colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, colour = 'black'))
#ggsave(filename = 'figures/VMH_neuron_dotplot_SCT.tiff', device = 'tiff', units = 'in', width = 4, height = 4, dpi = 300)

gc()

