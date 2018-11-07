#Initial code to dissect the interactions between cell-types in the developing thalamus. Performing analysis 
#with just one sample. 

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
path<-"/Users/michaeljohndolan/Desktop/Kalish_data/"
setwd(path)

#Read in a single tsv file for analysis, careful b/c first col is row names. File is 1.5GB! Transpose to get 
#genes x columns for Seurat 
P5_1<-read.table(file = "GSM2971232_P10_1.counts.tsv", header = T, row.names = 1, sep = '\t')
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
VlnPlot(object = P5_1.Seu, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,  point.size.use=0.005)
par(mfrow = c(1, 2))
GenePlot(object = P5_1.Seu, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.05)
GenePlot(object = P5_1.Seu, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.05)
median(P5_1.Seu@meta.data$percent.mito)
median(P5_1.Seu@meta.data$nUMI)
median(P5_1.Seu@meta.data$nGene)

#Using this data we filter out the cells that have high percent mito and abnormally high (by eye)
#nUMIs and nGenes. May need to write some code to do this in an automatic manner or maybe not.
#There are initially 14373 cells, giving us 2.3GB object. 
nrow(P5_1.Seu@meta.data)
P5_1.Seu<- FilterCells(object = P5_1.Seu, subset.names = "nUMI", 
                    low.thresholds = 500, high.thresholds = 15000)
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

#Run the PCA. Keep PCs 1-17. 
P5_1.Seu<- RunPCA(object = P5_1.Seu, pc.genes = P5_1.Seu@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
VizPCA(object = P5_1.Seu, pcs.use = 1:2)
PCElbowPlot(object = P5_1.Seu)

#Cluster the data 
P5_1.Seu<- FindClusters(object = P5_1.Seu, reduction.type = "pca", dims.use = 1:19, 
                     resolution = 0.8, print.output = 0, save.SNN = TRUE)

#Use the clustering to inform a tSNE on the significant PCs
P5_1.Seu<- RunTSNE(object = P5_1.Seu, dims.use = 1:19, do.fast = TRUE)
TSNEPlot(object = P5_1.Seu, do.label = TRUE, pt.size = 0.25)

#Find differentially expressed genes 
All.markers<-FindAllMarkers(object = P5_1.Seu, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#Identify the microglia
FeaturePlot(object = P5_1.Seu, features.plot = c("Cx3cr1", "C1qa"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne", pt.size=0.8)

#Extract microglia cluster only 
filter(P5_1.Seu@meta.data, res.0.8==10)
P10_1_mgl <- SubsetData(object = P5_1.Seu,ident.use = 9 )
saveRDS(P10_1_mgl, "P10_1_mgl.rds")
