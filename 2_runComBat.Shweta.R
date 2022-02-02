############
##
## This script runs ComBat on the log2TPM patient samples 
##
############
suppressMessages(library(dplyr))
suppressMessages(library(plyr))

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(sva)) ## for ComBat
suppressMessages(library(doParallel))
#registerDoParallel(cores=20)
#rm(list=ls()); 
#source("~/.Rprofile");

setwd("/Users/schavan/Box/schavan/a2idea-related/")
rawd <- fread("testing/abundance_genelevel_tpm.csv",data.table = F); dim(rawd);head(rawd)

##Before batch correction
##Use Kallisto files directly?
##Use only protein coding genes
##Avoid Zeros by adding 1 to TPM?

Ns <- ncol(rawd) - 1 ## get out the number of samples in this file
Ns
pheno <- fread("testing/SampleSheet.txt")


# rawd <- fread("testing/log2TPM.N35.tsv",data.table = F); dim(rawd);head(rawd)
# 
# ## This code constructs the phenotype data.frame for the patients
# 
# pheno=read.csv("testing/SampleCancerType_relabel_3jan21.csv", as.is=T)
# pheno$tag <- gsub("^([^_]+)_.*", "\\1", pheno$providerID)
# pheno$tag[ -grep("^F", pheno$tag) ] <- paste0("X", pheno$tag[ -grep("^F", pheno$tag) ])
# pheno <- pheno[ pheno$providerID != "failed", ]

## Here we are making sure pheno and rawd have the same column ordering
pheno$index = 0
for(i in 1:nrow(pheno)) {
	pheno$index[i] <- match(pheno$V1[i], colnames(rawd))
}
rm(i)

rawd <- rawd[, c(1,pheno$index) ]; dim(rawd)

colnames(pheno)[1] <- "sampleID"
pheno$tag = as.factor(pheno$sampleID)
pheno$Diagnosis <- as.factor(pheno$V5)
pheno$Batch  <- as.factor(pheno$V3)
pheno$Provider  <- as.factor(pheno$V6)
pheno$DosePoint <- as.factor(pheno$V2)
pheno <- select(pheno, -c(index,V2,V3,V5,V6))

mat <- as.matrix(rawd[,-c(1)])
rownames(mat) <- rawd[,1]
#=======================================================================


## Trim the expression matrix
## Remove genes that have less than 0.5 counts across 25% of the samples
scoreRow <- function(x, cutOff=0.25, pct=0.25) {
	j <- sum(x > cutOff)/length(x)
	ret = FALSE
	if(j >= pct) ret = TRUE
	return(ret)
}
j <- apply(mat, 1, scoreRow)
fmat <- mat[j,] ## filtered matrix
rm(j, scoreRow)



#==========================================================================================
## Function to make PCA plots
getPlot <- function(expr_mat, feat, title=NULL) {
	## expr_mat: matrix (genes as rows)
	## feat: the feature from 'pheno' you want to plot

	df_pca <- prcomp( t(expr_mat) )
	sm_df <- summary(df_pca)$importance
	
	tmp <- data.frame(Sample=rownames(df_pca$x), PC1=df_pca$x[,1], PC2=df_pca$x[,2])
	tmp <- inner_join(x=pheno, y=tmp, by=c("tag"="Sample"))

	if(is.null(title)) title <- feat;

	if(feat == "labels") {
		ret <- ggplot(tmp) + 
			geom_point(aes(x=PC1, y=PC2), color="red", size=1) +
			geom_text_repel(aes(x=PC1, y=PC2, label=sampleID)) + 
			xlab(paste0("PC1, ",round(sm_df[2,1]*100,1),"%")) + 
			ylab(paste0("PC2, ",round(sm_df[2,2]*100,1),"%")) + 
			ggtitle(title) + 
			theme_bw()
	} else {
		tmp <- select(tmp, sampleID, PC1, PC2, feat)
		colnames(tmp)[4] <- "feature"
		ret <- ggplot(tmp) + 
			geom_point(aes(x=PC1, y=PC2, color=feature), alpha=0.5, size=3) + 
			xlab(paste0("PC1, ",round(sm_df[2,1]*100,1),"%")) + 
			ylab(paste0("PC2, ",round(sm_df[2,2]*100,1),"%")) + 
			ggtitle(title) + 
			theme_bw() +
			theme(legend.title=element_blank())
	
	}
	return(ret)
}
#==========================================================================================


## Get a PCA plot with and without removing genes
p0_noFilter <- getPlot(mat, feat="labels", title="Before ComBat, no filter")
p0_withFilter <- getPlot(fmat, feat="labels", title="Before ComBat, with filter")


## These are constants you'll need to perform ComBat correction
M <- select(pheno, -sampleID, -V4, -tag, -Provider)
rownames(M) <- pheno$tag
modcombat <- model.matrix(~1, data=M)
batch <- M$Batch

## Perform batch correction on unfiltered gene set
combat_edata <- ComBat(dat=mat, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F)

## Perform batch correction on filtered gene set
fcombat_edata <- ComBat(dat=fmat, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F)


## PCA plots
p01_noFilterPostCombat <- getPlot(combat_edata, feat="labels", title="Post ComBat, no filter")
p01_withFilterPostCombat <- getPlot(fcombat_edata, feat="labels", title="Post ComBat, with filter")

p1_noFilter <- getPlot(combat_edata, feat="Diagnosis", title="Diagnosis, Post ComBat, no filter")
p1_withFilter <- getPlot(fcombat_edata, feat="Diagnosis", title="Diagnosis, Post ComBat, with filter")

p2_noFilter <- getPlot(combat_edata, feat="DosePoint", title="DosePoint, Post ComBat, no filter")
p2_withFilter <- getPlot(fcombat_edata, feat="DosePoint", title="DosePoint, Post ComBat, with filter")

p3_noFilter <- getPlot(combat_edata, feat="Batch", title="Batch, Post ComBat, no filter")
p3_withFilter <- getPlot(fcombat_edata, feat="Batch", title="Batch, Post ComBat, with filter")

#p4_noFilter <- getPlot(combat_edata, feat="Provider", title="Provider, Post ComBat, no filter")
#p4_withFilter <- getPlot(fcombat_edata, feat="Provider", title="Provider, Post ComBat, with filter")

outF=paste0("testing/pca_plots.N",Ns,".pdf")
pdf(file=outF, width=10, height=10)
plot(p0_noFilter)
plot(p0_withFilter)
plot(p01_noFilterPostCombat)
plot(p01_withFilterPostCombat)
plot(p1_noFilter)
plot(p1_withFilter)
plot(p2_noFilter)
plot(p2_withFilter)
plot(p3_noFilter)
plot(p3_withFilter)
plot(p4_noFilter)
plot(p4_withFilter)
dev.off()



## write out both matrices (with and without filtering)
out <- data.frame(round(mat, 2), stringsAsFactors=F)
outF=paste0("testing/comBat.noGeneFilter.log2TPM.N",Ns,".tsv")
write.table(out, file=outF, sep="\t", row.names=F)
rm(out)

out <- data.frame(round(fmat, 2), stringsAsFactors=F)
outF=paste0("testing/comBat.withGeneFilter.log2TPM.N",Ns,".tsv")
write.table(out, file=outF, sep="\t", row.names=F)
#rm(out)

out$Gene = row.names(out)
dt = as.data.table(out)
dt %>% filter(Gene %in% c('BIRC5', 'CTAG1B', 'PRAME', 'SSX2', 'WT1'))
