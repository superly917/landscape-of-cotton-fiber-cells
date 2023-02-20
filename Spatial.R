library("Seurat")
library("Matrix")
library("dplyr")
library("tibble")
library("ggplot2")
library("patchwork")
library("cowplot")
library("dplyr")
library("gridExtra")
library("RColorBrewer")
library("MAST")

# Create seurat object
dir = c('spaceranger/global_embryo/outs','spaceranger/torpedo_embryo/outs','spaceranger/cytoledon_embryo/outs')
names(dir) = c('global_embryo', 'torpedo_embryo', 'cytoledon_embryo')
seu_list <- list()
for(i in 1:length(dir)){
seu_list[[i]] <-Seurat::Load10X_Spatial(data.dir = dir[i])
seu_list[[i]]@meta.data$orig.ident <-names(dir)[i]
}

seurat_ob <- merge(seu_list[[1]], y=c(seu_list[[2]], seu_list[[3]]))

# QC
plot1 <- VlnPlot(seurat_ob, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seurat_ob, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

plot1 <- VlnPlot(seurat_ob, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seurat_ob, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

mito.genes <- grep(pattern = "^(MT|mt)(-|_)", x = rownames(seurat_ob), value = T,perl=T)
raw.counts = GetAssayData(seurat_ob, slot = "counts")
percent.mito <- Matrix::colSums(raw.counts[mito.genes, ])/Matrix::colSums(raw.counts)
seurat_ob <- AddMetaData(object = seurat_ob, metadata = percent.mito, col.name = "percent.mito")
plot1 <- VlnPlot(seurat_ob, features = "percent.mito", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seurat_ob, features = "percent.mito") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

DefaultAssay(brain.merge) <- "SCT"
seurat_ob <- SCTransform(seurat_ob, assay = "Spatial", vars.to.regress = "nCount_Spatial", verbose = FALSE, return.only.var.genes = FALSE)

# Clustering
seurat_ob@meta.data$batchid = seurat_ob@meta.data$orig.ident
seurat_ob <- RunMnn.Seurat(seurat_ob, features = VariableFeatures(seurat_ob) ,
							   assay = DefaultAssay(seurat_ob),
							   batch = batchid)
seurat_ob <- FindNeighbors(seurat_ob, reduction = "mnn", dims = 1:30,
							 features = VariableFeatures(seurat_ob),
							 nn.eps = 0, force.recalc = T, verbose = F)
seurat_ob <- FindClusters(object = seurat_ob, resolution = 0.4, algorithm = 1, verbose = F)
seurat_ob <- RunTSNE(seurat_ob, reduction = "mnn",
                                   dim.embed = 2, tsne.method = "Rtsne",
                                   features = VariableFeatures(singlecell_ob), perplexity = 30,
                                   dims = 1:30, num_threads = 10,
                                   max_iter = 2000, check_duplicates = F)
p1 <- DimPlot(seurat_ob, reduction = "tsne", label = TRUE)
p2 <- SpatialDimPlot(seurat_ob, label = TRUE, label.size = 3)
p1 + p2

# Find marker
global_DEGs = FindAllMarkers(object = seurat_ob,only.pos = T,test.use = 'bimod', logfc.threshold = 0, min.pct = 0.25)
global_DEGs = global_DEGs %>% mutate( gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>% select( gene, everything())
write.table(global_DEGs,"all_DEGs_for_all_clusters.xls", col.names =T, row.names = F,sep = "\t")

topn_markers  = global_DEGs %>% group_by(cluster) %>% 
              #filter(p_val < opt$pvalue ) %>%
              arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
              # filter(gene_diff > pct_fold_cutoff)  %>% 
              top_n(10,gene_diff)

write.table(topn_markers, "top10_for_each_clusters.xls", col.names =T,row.names = F,sep = "\t")


topn_markers <- data.frame()
markers2vis <- read.table("top10_for_each_clusters.xls", sep = "\t", header = T)
topn_markers <- markers2vis %>%
	group_by(cluster) %>%
	arrange(p_val, desc(avg_logFC), desc(gene_diff)) %>%
	top_n(10, .data[['avg_logFC']]) %>%
	arrange(cluster) %>%
	mutate(folder_suffix = paste0("cluster", cluster)) %>%
	select(cluster, gene, folder_suffix)
markers2vis4heatmap <- unique(as.vector(topn_markers$gene))

ggheat <- DoHeatmap(
	object = seurat_ob,
	features = markers2vis4heatmap,
	group.by = "clusters",
	group.bar = T,
	label = F
) + theme(axis.text.y = element_text(size = 4))
ggheat + guides(fill = guide_colorbar(title.position = "top", order = 1), color = guide_legend(order = 2, override.aes = list(alpha = 1)))
ggsave("topmarker_gene_heatmap.pdf")


# Diffexp
res = Seurat::FindMarkers(seurat_ob, ident.1 = "cytoledon_embryo",
				  ident.2 = "torpedo_embryo",
				  group.by = "sampleid",
				  test.use = "MAST",
				  only.pos = F )
res = res %>% tibble::rownames_to_column(var = "gene") %>%
	  dplyr::rename( pvalue = p_val, padj = p_val_adj)
res = res %>% dplyr::select(gene,everything())  
res$FoldChange = exp(1)^res$avg_logFC
res = res %>% select(-avg_logFC)
res$log2FoldChange = log2(res$FoldChange)
write.table(res, "diffexp_genes.xls", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na="")

