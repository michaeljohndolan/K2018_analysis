#Code to merge the samples for each timepoint and process, analyze and save all the data. 
#This was run on google cloud 

#Load up required libraries and custom functions
library(Seurat)
library(Matrix)
library(dplyr)
QC.plotter<-function(object) {
  data<-object@meta.data
  median.mito<-median(data$percent.mito); print(paste0("Median mito is: ", median.mito))
  median.nUMI<-median(data$nUMI); print(paste0("Median nUMI is: ", median.nUMI))
  
  g<-ggplot(data, aes(x=nUMI, y=percent.mito))
  g<-g+geom_point(size=0.5, alpha=0.05)
  g<-g+coord_cartesian(expand = F)
  g<-g+geom_hline(yintercept = median.mito, color="red")
  g<-g+geom_vline(xintercept = median.nUMI, color="green")
  g
}
name.extract<-function(name) {
  name<-sapply(strsplit(name, "\\Q.\\E"), "[", 1)
  name<-paste0(sapply(strsplit(name, "_"), "[", 2), "_", sapply(strsplit(name, "_"), "[", 3))
  name
}

#Set count matrix data path (tsv files, not included in repo as very large)
path<-"/data/mike/Kalish_Analysis/data/data/"
setwd(path)
samples<-list.files()

#List out the different timepoints
time.points<-c("P5_", "P10_", "P16_", "P21_")

#Loops through each timepoint, creates a Seurat object for first sample and iteratively merges the samples
for(i in 1:4){ 
  P.X<-time.points[i]
  P.X.tsv<-grep(pattern = P.X, x = samples, value = TRUE, fixed = TRUE)
  for(k in 1:length(P.X.tsv)) {
    if(k==1) {
      object<-read.table(file =P.X.tsv[k] , header = TRUE, row.names = 1, sep = '\t')
      object<-t(object)
      object<-CreateSeuratObject(raw.data = object , min.cells = 3, min.genes = 200, 
                                 project =P.X.tsv[k] )
    } #Will initialize a Seurat object 
    if(k>1) {
      temp<-read.table(file =P.X.tsv[k] , header = TRUE, row.names = 1, sep = '\t')
      temp<-t(temp)
      temp<-CreateSeuratObject(raw.data = temp , min.cells = 3, min.genes = 200, 
                               project =P.X.tsv[k])
      object<-MergeSeurat(object1 = object, object2 = temp, project = P.X, do.normalize = FALSE
                          , add.cell.id1 = k ,add.cell.id2 = P.X.tsv[k])
      rm(temp)
    } #Will merge subsequent Seurat samples
  }
  saveRDS(object = object, file = paste0(P.X, "aggregated.rds"))
  rm(object)
}

##Fix up the directory organization. Then load in each object into a single large seurat object for QC  (maintaining timepoint info)
#the run analysis and clustering. Then pull out all the microglia in a separate cluster and analyse/differential expression.
P5<-readRDS("P5_aggregated.rds"); P5@meta.data$orig.ident<-"P5"
P10<-readRDS("P10_aggregated.rds"); P10@meta.data$orig.ident<-"P10"
P16<-readRDS("P16_aggregated.rds"); P16@meta.data$orig.ident<-"P16"
P21<-readRDS("P21_aggregated.rds"); P21@meta.data$orig.ident<-"P21"

#Aggregate the Seurat objects to a single object. 
all.kalish<-MergeSeurat(object1 = P5, object2 = P10, do.normalize = FALSE, add.cell.id1 = "P5", add.cell.id2 = "P10")
all.kalish<-MergeSeurat(object1 = all.kalish, object2 = P16, do.normalize = FALSE, add.cell.id2 = "P16")
all.kalish<-MergeSeurat(object1 = all.kalish, object2 = P21, do.normalize = FALSE, add.cell.id2 = "P21")

#First part of cell QC, calculate the percent mito genes expressed for each cell and add this to the metadata slot. 
mito.genes <- grep(pattern = "^mt.", x = rownames(x = all.kalish@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(all.kalish@raw.data[mito.genes, ])/Matrix::colSums(all.kalish@raw.data)
all.kalish<-AddMetaData(object = all.kalish, metadata = percent.mito, col.name = "percent.mito")

#Plot cellQC parameters from the metadata slot with violin plot and geneplot. 
VlnPlot(object = all.kalish, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, point.size.use=0.005)
par(mfrow = c(1, 2))
GenePlot(object = all.kalish, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.05)
GenePlot(object = all.kalish, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.05)
median(all.kalish@meta.data$percent.mito)
median(all.kalish@meta.data$nUMI)
median(all.kalish@meta.data$nGene)

#Filter the dataset by cells
nrow(all.kalish@meta.data)
all.kalish<- FilterCells(object = all.kalish, subset.names = c("nGene", "nUMI","percent.mito"), 
                         low.thresholds = c(200, 500, -Inf), high.thresholds = c(2500, 15000, 0.1))
nrow(all.kalish@meta.data)

#Normalize the data and find the variable genes
all.kalish <- NormalizeData(object = all.kalish, normalization.method = "LogNormalize", scale.factor = 10000)
all.kalish<-FindVariableGenes(object = all.kalish, mean.function = ExpMean, dispersion.function = LogVMR, 
                              do.plot = TRUE, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

#Scale the data accross the different genes and regress out variations due to nUMI (high expression) and per.cent mito
# (genes correlated with mt expression). Might play around with this b/c presumable microglia during pruning use up 
# a lot of energy. 
all.kalish <- ScaleData(object = all.kalish, genes.use = all.kalish@var.genes, display.progress = TRUE, 
                        vars.to.regress = c("percent.mito",'nUMI'), do.par = TRUE, num.cores = 4)

#Run the PCA. Keep PCs 1-30. 
all.kalish<- RunPCA(object = all.kalish, pc.genes = all.kalish@var.genes,pcs.compute = 30 , do.print = TRUE, pcs.print = 1:5, 
                    genes.print = 5)
PCElbowPlot(object = all.kalish)

#Cluster the data with 15 PCs
all.kalish <- FindClusters(object = all.kalish, reduction.type = "pca", dims.use = 1:18, 
                           resolution = 0.6, print.output = 0, save.SNN = TRUE)
#the object in case GC crashes
saveRDS(all.kalish, "all.kalish.rds")

## Run and plot tSNE
all.kalish <- RunTSNE(object = all.kalish, dims.use = 1:18, do.fast = TRUE)
TSNEPlot(object = all.kalish, pt.size = 0.05, do.label = TRUE)

#Save the final object with tSNE calculated. 
saveRDS(all.kalish, "all.kalish.rds")

