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

#Set count matrix data path (tsv files, not included in repo as very large)
path<-"/Users/mdolan/Desktop/Kalish2018_Data/"

