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

#Initialize the Seurat Object 
P5_1.Seu<-CreateSeuratObject(raw.data = P5_1 , min.cells = 3, min.genes = 200, 
                             project = "P5_1")

#Perform the cell QC, calculate the percent mito 
mito.genes <- grep(pattern = "^MT-", x = rownames(x = P5_1.Seu@raw.data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)





