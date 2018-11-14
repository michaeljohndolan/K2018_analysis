#Code to examine all.mgls from Kalish2018 in more detail. 
library(Seurat)
library(Matrix)
library(dplyr)

#Read in the all.mgl dataset
setwd("~/Google Drive (mdolan@broadinstitute.org)/Misc_Projects/Thalamus_scRNAseq")
all.mgls<-readRDS(file = "all.mgls.rds")

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

#Identify the subsets within P5, see C1qa, C1qb and C1qc are all enriched in cluster 1 (216)
P5_subtype_markers<-FindAllMarkers(object = P5)
P5_subtype_markers<-mutate(P5_subtype_markers, diff=pct.1-pct.2)
P5_subtype_markers<-filter(P5_subtype_markers, p_val_adj<0.05)
P5_subtype_markers<-arrange(P5_subtype_markers, desc(diff), avg_logFC)
View(P5_subtype_markers)
View(filter(P5_subtype_markers, pct.2<0.2 ,p_val_adj<0.05))


#Known microglial genes with roles in pruning/lysosomes and phagocytosis
DotPlot(object = P5, genes.plot =  c("C1qa", "C1qb", "C1qc", "Grn", "Sirpa", "Itgam", "Laptm5", "Hmgb2") , plot.legend = TRUE)

#AD genes in P5 microglia
DotPlot(object = P5, genes.plot =  c("Tyrobp", "Apoe", "Trem2") , plot.legend = TRUE)

#Previously unknown interesting genes: CCL3 (chemokine, implicated in retinal damage), Vsir (inhibitory immune checkpoint)
#Fcer1g (viral), pathogen recognition (Fcgr3), Trem2, Apoe, fragile X interactor Cyfip1, 
DotPlot(object=P5, genes.plot=c("Ccl3", "Vsir", "Cyfip1", "Mpeg1", "Hpgds"), plot.legend = TRUE)

DotPlot(object=P5, genes.plot=c("Cxcl2", "Hspa1b", "Hspa1a", "Tnfaip3", "Cd86"), plot.legend = TRUE)


#What genes are highly correlated with C1qa in this P5 dataset 
matrix<-P5@data
matrix_mod<-as.matrix(matrix)
gene<-as.numeric(matrix_mod["C1qa",])
correlations<-apply(matrix_mod,1,function(x){cor(gene,x)})
correlations<-na.omit(correlations)
View(correlations)

#Next step, compare these pruning microglia to all microglia in the Kalish dataset 
#Do some gene ontology, look at P5 neurons and find a "pruned neuron" signature, look at Tim's dataset. 
#Examine the other microglia at P5, what are they doing. Read Tim's paper 

