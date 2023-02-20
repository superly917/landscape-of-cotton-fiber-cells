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
library("DoubletFinder")

# Create seurat object
dir = c('cellranger/global_embryo/outs/filtered_feature_bc_matrix','cellranger/torpedo_embryo/outs/filtered_feature_bc_matrix','cellranger/cytoledon_embryo/outs/filtered_feature_bc_matrix')
names(dir) = c('global_embryo', 'torpedo_embryo', 'cytoledon_embryo')
seu_list <- list()
for(i in 1:length(dir)){
countMatrixSparse <-Seurat::Read10X(data.dir = dir[i])
seu_list[[i]] <- CreateSeuratObject( countMatrixSparse[["Gene Expression"]], names.field = 2, assay = "RNA", names.delim = "-" )

seu_list[[i]]@meta.data$orig.ident <-names(dir)[i]
}

seurat_ob <- merge(seu_list[[1]], y=c(seu_list[[2]], seu_list[[3]]))

# QC
plot1 <- VlnPlot(seurat_ob, features = "nCount_RNA", pt.size = 0.1,x.lab.rot = T,group.by = "sampleid") + ggtitle("nUMI")+labs(xlab = "") + theme( plot.title = element_text(hjust = 0.5))
ggsave("nUMI_beforeQC_vlnplot.pdf",plot=plot1)
ggsave("nUMI_beforeQC_vlnplot.png",plot=plot1)

plot2 <- VlnPlot(seurat_ob, features = "nFeature_RNA", pt.size = 0.1,x.lab.rot = T,group.by = "sampleid") + ggtitle("nGene")+labs(xlab = "") + theme( plot.title = element_text(hjust = 0.5))
ggsave("nGene_beforeQC_vlnplot.pdf",plot=plot2)
ggsave("nGene_beforeQC_vlnplot.png",plot=plot2)

DefaultAssay(seurat_ob)="RNA"
UMIs_per_cell = Matrix::colSums(raw.counts)
genes_per_cell = Matrix::colSums(raw.counts>0)
df = data.frame(UMIs_per_cell=UMIs_per_cell, genes_per_cell=genes_per_cell)
suppressPackageStartupMessages(library(MASS))
df = df[order(df$UMIs_per_cell),]
pdf(file.path(output_dir,"outliers.pdf"))
plot(df, log='xy')
m <- rlm(genes_per_cell~UMIs_per_cell,data=df)
p.level = 1e-3
suppressWarnings(pb <- data.frame(predict(m, interval='prediction',
level = 1-p.level,
type="response")))
polygon(c(df$UMIs_per_cell, rev(df$UMIs_per_cell)), c(pb$lwr, rev(pb$upr)), col=adjustcolor(2,alpha=0.1), border = NA)
outliers <- rownames(df)[df$genes_per_cell > pb$upr | df$genes_per_cell < pb$lwr]
points(df[outliers,],col=2, bg="red",pch = 16,cex=0.6)
title("Outlier of Cells")
dev.off()

seurat_ob = subset(seurat_ob, cells = colnames(raw.counts)[!colnames(raw.counts) %in% outliers])

CellRemover = function( object, parameters,lower_limit, high_limit, sdfold ){
    lower_threshold = vector()
    upper_threshold = vector()
    param_list = list()
    for ( paramx in parameters ){
        param_value = object@meta.data[,paramx]
        param_vector = param_value[param_value>0]
        mean4paramx = mean(log10(param_vector))
        sd4paramx = sd(log10(param_vector))
        upper_bound <- 10^(mean4paramx + sdfold*sd4paramx)
        lower_bound <- 10^(mean4paramx - sdfold*sd4paramx)
        if ( is.null( lower_limit[paramx]) | is.na( lower_limit[paramx]) |lower_limit[paramx] %in% c("NA","NULL") ){
            lower_threshold = lower_bound
        }else{
            lower_threshold = as.numeric(lower_limit[paramx])
        }
        if ( is.null( high_limit[paramx]) | is.na( high_limit[paramx]) | high_limit[paramx] %in% c("NA","NULL") ){
            upper_threshold = upper_bound
        }else{
            upper_threshold = as.numeric(high_limit[paramx])
        }
        suppressWarnings({
			threshold =  paste0( paramx, " > ", lower_threshold , " & " ,  paramx, " < " , upper_threshold ) 
			df = object@meta.data
			desired_cells= subset(df, eval( parse(text=threshold)))
			object = object[, rownames(desired_cells)]
            # object = SubsetData(object, subset.name = paramx,
            # low.threshold = lower_threshold, high.threshold = upper_threshold)
        })
    }
    return( object )
}

seurat_by_sample = future_lapply(SplitObject(seurat_ob, split.by = "sampleid" ),
function(x)CellRemover(object = x, parameters = c( "nFeature_RNA", "nCount_RNA",),lower_limit = c( nFeature_RNA = "NULL", nCount_RNA = "NULL"),
high_limit =  c( nFeature_RNA = "NULL", nCount_RNA = "NULL"), sdfold = 2 ))
merged_seurat = seurat_by_sample[[1]]
for( idx in 2:length(seurat_by_sample) ){
	merged_seurat = merge(x = merged_seurat, y = seurat_by_sample[[idx]],
	do.scale = F, do.center = F, do.normalize = F)
}
seurat_ob = merged_seurat


RemoveDoublets <-function(
  object,
  identity,
  doublet.rate,
  pN=0.25,
  PCs=1:30,
  use.SCT=FALSE,
  num.cores=1,
  quietly=TRUE
){

  tic("get sweep parameters")
  # calculate parameters
  if (quietly==TRUE){
    invisible(capture.output(sweep.res.list <- paramSweep_v3(object, PCs = PCs, sct=use.SCT, num.cores=num.cores)))
    invisible(capture.output(sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)))
    ff <- tempfile()
    png(filename=ff)
    invisible(capture.output(bcmvn <- 
find.pK(sweep.stats)))
    dev.off()
    unlink(ff)
  }else{
    sweep.res.list <- paramSweep_v3(object, PCs = PCs, sct=use.SCT, num.cores=num.cores)
    sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
    ff <- tempfile()
    png(filename=ff)
    bcmvn <- 
find.pK(sweep.stats)
    dev.off()
    unlink(ff)
  }
  toc()

  # choose parameters
  maxBCmetric    <- max(bcmvn$BCmetric, na.rm = TRUE)
  pK <- as.numeric(as.character(bcmvn[bcmvn$BCmetric==maxBCmetric, ]$pK))

  # compute doublet scores
  tic("Removing doublets")
  annotations    <- object@meta.data$identity  ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi       <- round(doublet.rate*length(colnames(x = object)))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj   <- round(nExp_poi*(1-homotypic.prop))
  seu.scored     <- doubletFinder_v3(object, PCs =PCs, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = use.SCT)
  toc()

  # pick out doublets
  cname <-colnames(seu.scored[[]])
  DF<-cname[grep('^DF',cname)]
  seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")

  # remove doublets
  seu.removed <- subset(seu.scored, subset = doublet != 1)
  return(list(removed=seu.removed, original=seu.scored))
}
rate=function(num){
	if (num >10000) rate= 0.1 else rate= num*0.0000077
	return(rate)
}

seurat_ob = SetIdent( seurat_ob, value = "sampleid")
obj = SplitObject(seurat_ob,split.by = "sampleid")
obj_rm=list()
for( i in names(obj)){
	obj[[i]] <- NormalizeData(obj[[i]])
	obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
	obj[[i]] <- ScaleData(obj[[i]])
	obj[[i]] <- RunPCA(obj[[i]])
	obj[[i]] <- RunUMAP(obj[[i]], dims = 1:10)
	obj_rm[[i]] = RemoveDoublets(obj[[i]], doublet.rate=rate(dim(obj[[i]])[2]),  num.cores=4)
}
removed= lapply(obj_rm,FUN=function(x) x = x$removed)
if ( length(removed) > 1 ) {
	seurat_ob = merge(removed[[1]],do.call(c,removed[-1]))
} else {
	seurat_ob = removed[[1]]
}


plot1 <- VlnPlot(seurat_ob, features = "nCount_RNA", pt.size = 0.1,x.lab.rot = T,group.by = "sampleid") + ggtitle("nUMI")+labs(xlab = "") + theme( plot.title = element_text(hjust = 0.5))
ggsave("nUMI_afterQC_vlnplot.pdf",plot=plot1)
ggsave("nUMI_afterQC_vlnplot.png",plot=plot1)

plot2 <- VlnPlot(seurat_ob, features = "nFeature_RNA", pt.size = 0.1,x.lab.rot = T,group.by = "sampleid") + ggtitle("nGene")+labs(xlab = "") + theme( plot.title = element_text(hjust = 0.5))
ggsave("nGene_afterQC_vlnplot.pdf",plot=plot2)
ggsave("nGene_afterQC_vlnplot.png",plot=plot2)


seurat_ob <- NormalizeData(object = seurat_ob,
						normalization.method = 'LogNormalize',scale.factor = 10000)
vars2regress =c("nCount_RNA", "percent.mito")
seurat_ob = FindVariableFeatures(object= seurat_ob, loess.span = 0.3,
					clip.max = "auto", mean.function = "FastExpMean",
					dispersion.function = "FastLogVMR", num.bin = 20,
					nfeature = 4000, binning.method = "equal_width" )
seurat_ob <- ScaleData(object = seurat_ob, features = rownames(seurat_ob),
					vars.to.regress = vars2regress, verbose = T )


# Clustering
seurat_ob <- RunPCA(seurat_ob, features = VariableFeatures(seurat_ob) ,
							   assay = DefaultAssay(seurat_ob),
							   npcs = 30)
seurat_ob <- FindNeighbors(seurat_ob, reduction = "pca", dims = 1:30,
							 features = VariableFeatures(seurat_ob),
							 nn.eps = 0, force.recalc = T, verbose = F)
seurat_ob <- FindClusters(object = seurat_ob, resolution = 0.4, algorithm = 1, verbose = F)
seurat_ob <- RunTSNE(seurat_ob, reduction = "pca",
                                   dim.embed = 2, tsne.method = "Rtsne",
                                   features = VariableFeatures(singlecell_ob), perplexity = 30,
                                   dims = 1:30, num_threads = 10,
                                   max_iter = 2000, check_duplicates = F)
p1 <- DimPlot(seurat_ob, reduction = "tsne", label = TRUE)
ggsave('tsne_plot.pdf', plot=p)
ggsave('tsne_plot.png', plot=p)

# Find marker
global_DEGs = FindAllMarkers(object = seurat_ob,only.pos = T,test.use = 'bimod', logfc.threshold = 0, min.pct = 0.25)
global_DEGs = global_DEGs %>% mutate( gene_diff = round(global_DEGs$pct.1 / global_DEGs$pct.2, 3)) %>% select( gene, everything())
write.table(global_DEGs,"all_DEGs_for_all_clusters.xls", col.names =T, row.names = F,sep = "\t")

topn_markers  = global_DEGs %>% group_by(cluster) %>% 
              arrange(p_val,desc(avg_logFC),desc(gene_diff)) %>%
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

