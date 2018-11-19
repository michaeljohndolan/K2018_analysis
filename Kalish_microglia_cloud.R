#Code to load up the full, aggregated dataset and extract microglia and compare gene expression
#accross timepoints. 
#Run on the cloud 
library(Seurat)
library(Matrix)
library(ggplot2)

#Set path and load up the aggregated and completely processed object 
path<-"/data/mike/Kalish_Analysis/data/data/"
setwd(path)
all.kalish<-readRDS(file = "all.kalish.rds")

#Plot tSNE with labels 
TSNEPlot(object = all.kalish, pt.size = 0.05, do.label = TRUE)

#Examine the expression of different microglial markers and compare with macrophage markers
FeaturePlot(object = all.kalish, features.plot = c("Cx3cr1", "Mrc1", "Snap25", "Flt1", "Mbp" ,"Pdgfra"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=0.5)

#So cluster 12 appears to be microglia, while cluster 16 is macrophages. How many microglia are there
#accross the different timepoints. ~200-510 microglia per timepoint. 
table(all.kalish@meta.data[,c(3,5)])

#Extract these microglia for further analysis, only ~1800 single cells. 
all.mgls<-SubsetData(object = all.kalish, ident.use = 12, subset.raw = TRUE, do.clean = TRUE )

#Repeat analyses on this microglial dataset 
all.mgls <- NormalizeData(object = all.mgls, normalization.method = "LogNormalize", scale.factor = 10000)
all.mgls<-FindVariableGenes(object = all.mgls, mean.function = ExpMean, dispersion.function = LogVMR
                            ,do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
all.mgls <- ScaleData(object = all.mgls, genes.use = all.mgls@var.genes, display.progress = TRUE
                      ,vars.to.regress = c("percent.mito",'nUMI'), do.par = TRUE, num.cores = 4)
all.mgls<- RunPCA(object = all.mgls, pc.genes = all.mgls@var.genes,pcs.compute = 30 , do.print = TRUE, pcs.print = 1:5
                  ,genes.print = 5)
PCElbowPlot(object = all.mgls)

#Find clusters
all.mgls <- FindClusters(object = all.mgls, reduction.type = "pca", dims.use = 1:16, 
                         resolution = 0.8, print.output = 0, save.SNN = TRUE)

## Run and plot tSNE
all.mgls <- RunTSNE(object = all.mgls, dims.use = 1:16, do.fast = TRUE)
TSNEPlot(object = all.mgls, pt.size = 1, do.label = TRUE)

#Save the full datset 
saveRDS(all.mgls, "all.mgls.rds")

#Plot some different features of this dataset
FeaturePlot(object = all.mgls, features.plot = c("C1qa", "Itgam"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=1)
FeaturePlot(object = all.mgls, features.plot = c("percent.mito"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=1)

#Switch identities to the different timepoints, plot and calculate genes that define P5
all.mgls <- SetAllIdent(object = all.mgls, id = "orig.ident")
TSNEPlot(object = all.mgls, pt.size = 1, do.label = TRUE)
P5_markers<-FindMarkers(object = all.mgls,ident.1 = "P5")

#Calcuate a difference metric and then sort by difference and fold change 
P5_markers$geneName<-row.names(P5_markers)
P5_markers<-mutate(P5_markers, diff=pct.1-pct.2)
P5_markers<-filter(P5_markers, p_val_adj<0.05)
P5_markers.ordered<-arrange(P5_markers, desc(diff), avg_logFC)
View(P5_markers.ordered)

#Examine the expression of several genes in the main and microglial datasets
FeaturePlot(object = all.mgls, features.plot = c("Apoe"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=1)
FeaturePlot(object = all.kalish, features.plot = c("Vgf"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=1)

#Specifically subset and examine the P5 dataset alone.
P5<-SubsetData(object = all.mgls, ident.use = "P5", subset.raw = TRUE, do.clean = TRUE )

#Repeat analyses on this microglial dataset 
P5 <- NormalizeData(object = P5, normalization.method = "LogNormalize", scale.factor = 10000)
P5<-FindVariableGenes(object = P5, mean.function = ExpMean, dispersion.function = LogVMR
                      ,do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
P5 <- ScaleData(object = P5, genes.use = P5@var.genes, display.progress = TRUE
                ,vars.to.regress = c("percent.mito",'nUMI'), do.par = TRUE, num.cores = 4)
P5<- RunPCA(object = P5, pc.genes = P5@var.genes,pcs.compute = 30 , do.print = TRUE, pcs.print = 1:5
            ,genes.print = 5)
PCElbowPlot(object = P5)

#Find clusters
P5 <- FindClusters(object = P5, reduction.type = "pca", dims.use = 1:15, 
                   resolution = 0.8, print.output = 0, save.SNN = TRUE)

## Run and plot tSNE
P5 <- RunTSNE(object = P5, dims.use = 1:15, do.fast = TRUE)
TSNEPlot(object = P5, pt.size = 1, do.label = TRUE)

#Identify the subsets within P5, see C1qa, C1qb and C1qc are all enriched in cluster 1
P5_subtype_markers<-FindAllMarkers(object = P5)
P5_subtype_markers<-mutate(P5_subtype_markers, diff=pct.1-pct.2)
P5_subtype_markers<-filter(P5_subtype_markers, p_val_adj<0.05)
P5_subtype_markers<-arrange(P5_subtype_markers, desc(diff), avg_logFC)
View(P5_subtype_markers)

#Known microglial genes with roles in pruning/lysosomes and phagocytosis
DotPlot(object = P5, genes.plot =  c("C1qa", "C1qb", "C1qc", "Grn", "Sirpa", "Itgam", "Laptm5", "Hmgb2") , plot.legend = TRUE)

#AD genes in P5 microglia
DotPlot(object = P5, genes.plot =  c("Tyrobp", "Apoe", "Trem2") , plot.legend = TRUE)

#Previously unknown interesting genes: CCL3 (chemokine, implicated in retinal damage), Vsir (inhibitory immune checkpoint)
#Fcer1g (viral), pathogen recognition (Fcgr3), Trem2, Apoe, fragile X interactor Cyfip1, 
DotPlot(object=P5, genes.plot=c("Ccl3", "Vsir", "Cyfip1", "Mpeg1", "Hpgds"), plot.legend = TRUE)

#What genes are highly correlated with C1qa in this P5 dataset 
matrix<-P5@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["C1qa",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations<-na.omit(correlations)


#Create a multipage PDF with markers for neurons, astros, oligos etc.
pdf(file = "/data/mike/Outputs/test.pdf")
FeaturePlot(object = all.mgls, features.plot = c("Apoe"), cols.use = c("grey", "blue") 
            ,reduction.use = "tsne", pt.size=1)
dev.off()