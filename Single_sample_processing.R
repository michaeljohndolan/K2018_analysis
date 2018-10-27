#Initial code to dissect the interactions between cell-types in the developing thalamus. Performing analysis 
#with just one sample. 

#Load up required libraries 
library(Seurat)
library(Matrix)
library(dplyr)

#Set count matrix data path (tsv files, not included in repo as very large)
path<-"/Users/michaeljohndolan/Desktop/Kalish_data/"
setwd(path)

#Read in a single tsv file for analysis, careful b/c first col is row names. File is 1.5GB! Transpose to get 
#genes x columns for Seurat 
P5_1<-read.table(file = "GSM2971228_P5_1.counts.tsv", header = T, row.names = 1, sep = '\t')
P5_1<-t(P5_1)

#Initialize the Seurat Object, contains 14373 cells
P5_1.Seu<-CreateSeuratObject(raw.data = P5_1 , min.cells = 3, min.genes = 200, 
                             project = "P5_1")

#Perform the cell QC, calculate the percent mito genes expressed for each cell and add this to the metadata slot. 
mito.genes <- grep(pattern = "^mt.", x = rownames(x = P5_1.Seu@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(P5_1.Seu@raw.data[mito.genes, ])/Matrix::colSums(P5_1.Seu@raw.data)
P5_1.Seu<-AddMetaData(object = P5_1.Seu, metadata = percent.mito, col.name = "percent.mito")

#Lets start to plot some cellQC parameters from the metadata slot with violin plot and geneplot. May need to write my own functions
#for plotting this much data at once! 
nrow(P5_1.Seu@meta.data)
VlnPlot(object = P5_1.Seu, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = P5_1.Seu, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.05)
GenePlot(object = P5_1.Seu, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.05)

#Using this data we filter out the cells that have high percent mito and abnormally high (by eye)
#nUMIs and nGenes. May need to write some code to do this in an automatic manner or maybe not.
nrow(P5_1.Seu@meta.data)
P5_1.Seu<- FilterCells(object = P5_1.Seu, subset.names = c("nUMI", "nGene", "percent.mito"), 
                    low.thresholds = c(0, 200, -Inf), high.thresholds = c(10000, 3000, 0.4))
nrow(P5_1.Seu@meta.data)

#Normalize the filtered dataset 
P5_1.Seu <- NormalizeData(object = P5_1.Seu, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

#Detect the most variable genes, play around with these to see what the higher average expression genes are? 
P5_1.Seu <- FindVariableGenes(object = P5_1.Seu , mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3.5, y.cutoff = 0.5)

#Scale the data accross the different genes and regress out variations due to nUMI (high expression) and per.cent mito
# (genes correlated with mt expression). Might play around with this b/c presumable microglia during pruning use up 
# a lot of energy. 
P5_1.Seu <- ScaleData(object = P5_1.Seu, vars.to.regress = c("nUMI", "percent.mito"))











