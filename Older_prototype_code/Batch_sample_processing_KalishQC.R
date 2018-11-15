#Code to run through all the samples and do the same cell QC as the Kalish paper. 

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
  name<-paste0(sapply(strsplit(name, "_"), "[", 2), "_", sapply(strsplit(name, "_"), "[", 3))
  name
}


#Set count matrix data path (tsv files, not included in repo as very large)
path<-"/Users/mdolan/Desktop/Kalish2018_Data/"
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

#Create aggregated dataset 

#Begin preprocessing on aggregated dataset.
#Normalize the data as per Kalish protocol. 
P5<- NormalizeData(object = P5, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
  
#Detect the most variable genes, play around with these to see what the higher average expression genes are? 
P5<- FindVariableGenes(object =  P5  , mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  
#Scale the data accross the different genes and regress out variations due to nUMI (high expression) and per.cent mito
P5 <- ScaleData(object = P5, block.size = 100)
  
  
#Run the PCA. Keep PCs 1-30, as per Kalish 
P5<- RunPCA(object = P5, pc.genes = P5@var.genes, do.print = TRUE, pcs.print = 1:5, 
                    genes.print = 5)

  
#what are significant PCs? 
PCElbowPlot(P5, num.pc = 1:30)  

#Cluster the data 
object<- FindClusters(object = object, reduction.type = "pca", dims.use = 1:10, 
                          resolution = 0.6, print.output = 0, save.SNN = FALSE)
  
#Use the clustering to inform a tSNE on the significant PCs
object<- RunTSNE(object = object, dims.use = 1:30, do.fast = TRUE)
TSNEPlot(object = object, do.label = TRUE, pt.size = 0.25)
 