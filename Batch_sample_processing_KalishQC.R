#Code to run through all the samples and do the same QC as the Kalish paper. 

#Cells with fewer than 500 or more than 15,000 UMI counts were excluded. Note that this is not a great way to do QC 
#accross different timepoints imo. The data were log-normalized and scaled to 10,000 transcripts per cell.
#Variable genes were identified using the following parameters: x.low.cutoff: 0.0125; x.high.cutoff: 3; y.cutoff: 0.5. We limited the
#analysis to the top 30 principal components. Clustering resolution
#was set to 0.6. 

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
  name<-paste0(sapply(strsplit(name, "_"), "[", 1), "_", sapply(strsplit(name, "_"), "[", 2))
  name
  }

#Set count matrix data path (tsv files, not included in repo as very large)
path<-"/Users/mdolan/Desktop/Kalish2018_Data/"
setwd(path)
samples<-list.files()

#Loop through the samples and run the Kalish QC and save the RDS file and print 
for(i in 1:length(samples)) {

  #Load the data and initialize the object 
  object_name<-name.extract(samples[i])
  object<-read.table(samples[i], header=TRUE, row.names=1, sep = '\t')
  object<-t(object)
  object<-CreateSeuratObject(raw.data = object , min.cells = 3, min.genes = 200, 
                             project =  object_name)
  
  #Calculate the percent mito for QC 
  mito.genes <- grep(pattern = "^mt.", x = rownames(x = object@raw.data), value = TRUE)
  percent.mito <- Matrix::colSums(object@raw.data[mito.genes, ])/Matrix::colSums(object@raw.data)
  object<-AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito")
  
  #Filter the cells as per Kalish. No nGene filtering! 
  object<- FilterCells(object =object, subset.names = c("nUMI"), 
                         low.thresholds = c(500), high.thresholds = 15000)
  
  saveRDS(object, file = paste0(object_name, "_Kalishfiltered.rds"))

}

  #Actually just want to do the below on the final dataset 
  #Normalize the filtered dataset 
  object<- NormalizeData(object = object, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
  
  #Detect the most variable genes, play around with these to see what the higher average expression genes are? 
  object<- FindVariableGenes(object =  object  , mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  
  #Scale the data accross the different genes and regress out variations due to nUMI (high expression) and per.cent mito
  object <- ScaleData(object = object, vars.to.regress = c("nUMI", "percent.mito"))
  
  
  #Run the PCA. Keep PCs 1-30, as per Kalish 
  object<- RunPCA(object = object, pc.genes = object@var.genes, do.print = TRUE, pcs.print = 1:5, 
                    genes.print = 5)

  
  #what are significant PCs? 
  
  #Cluster the data 
  object<- FindClusters(object = object, reduction.type = "pca", dims.use = 1:10, 
                          resolution = 0.6, print.output = 0, save.SNN = FALSE)
  
  #Use the clustering to inform a tSNE on the significant PCs
  object<- RunTSNE(object = object, dims.use = 1:30, do.fast = TRUE)
  TSNEPlot(object = object, do.label = TRUE, pt.size = 0.25)
 