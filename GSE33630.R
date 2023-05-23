# Title	Normal thyrocytes vs papillary vs anaplastic thyroid carcinomas
# Organism	Homo sapiens
# Experiment type	Expression profiling by array
# Summary	We profiled the gene expression of 11 anaplastic thyroid carcinomas (ATC), 49 papillary thyroid carcinomas (PTC) and 45 normal thyroids (N)
# 
# Overall design	We hibridized a series of anaplastic thyroid carcinomas (ATC) and papillary thyroid carcinomas (PTC) onto Affymetrix U133 Plus 2.0 arrays. ATCs were obtained from different hospitals in France and Belgium. Paired RNA samples of PTCs and non-tumoral thyroid tissues were obtained from Ukraine via the Chernobyl Tissue Bank (www.chernobyltissuebank.com). Diagnoses were confirmed by the members of the International Pathology Panel of the Chernobyl Tissue Bank.

#Platforms (1)	
#GPL570	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
#RMA Normalized

# source("http://www.bioconductor.org/biocLite.R")
#BiocManager::install("bioconductor")
#
library("BiocManager")
 
if(!require(GEOquery)) library(GEOquery)
# ###   source("GEOparche") EN LA UNI
if(!require(Biobase)) library(Biobase)
# 
 
## EN CASA ########################################################################
gse1 = getGEO("GSE33630") #, AnnotGPL = T, getGPL = T) #ACA CAMBIAR GSE
es1=gse1$GSE33630_series_matrix.txt.gz
e1 = exprs(es1) #ACA CAMBIAR GSE  
dim(e1) 
# [1] 54675   105 - es lo que dice en GEO al ver los valores de una muestra:
feno1=pData(es1) ## Muestras
geno1=fData(es1) ## Genes
# geno$ID[1:10]
# geno$`Gene Symbol`[1:10]

## DEJAR SOLO MUESTRAS CONTROL y PTC
## SACAMOS ATC!
tejido1 = as.factor(feno1$characteristics_ch1)
summary(tejido1)
iATC1 = which(grepl("ATC",tejido1))
e1 = e1[,-iATC1]
feno1 = feno1[-iATC1,]

## FACTOR // EXPRESION DIFERENCIAL CON LIMMA
 
tejido1 = as.factor(feno1$characteristics_ch1) # ahora sin ATC
summary(tejido1)
iPTC1 = which(grepl("PTC",tejido1))

fc1 = rep(0, ncol(e1)) 
fc1[iPTC1] = 1
fc1 = as.factor(fc1)
summary(fc1)

library(limma)
head(e1)
dim(e1)
length(fc1)

design1 = model.matrix(~0 + fc1)
colnames(design1) <- c("CONTROL", "PTC")
fit1 = lmFit(e1, design1) # CREA EL MODELO LINEAL
contrast.matrix1 = makeContrasts("CONTROL-PTC", levels = design1)
PTC_fits1 <- contrasts.fit(fit1, contrast.matrix1)

# Given a microarray linear model fit, compute moderated t-statistics, 
# moderated F-statistic, and log-odds of differential expression 
# by empirical Bayes moderation of the standard errors towards a common value

ebFit1 <- eBayes(PTC_fits1)
tT1 = topTable(ebFit1, number=nrow(e1), coef=1, adjust.method = "BH", p.value=0.01) 
dim(tT1) # 18527
## EJEMPLOS uno FC positivo, otro FC negativo
head(tT1)
boxplot(e1[rownames(tT1)[1],]~fc1) # logFC=4.2; >0, sub expr en ptc
boxplot(e1[rownames(tT1)[3],]~fc1) # logFC=-2.07; <0, sobre expr en ptc

probes1=rownames(tT1) # 18527
conGeneSymbol1 = filter(geno1, ID %in% probes1)
conGeneSymbol1 = subset(conGeneSymbol1, select=c(ID,`Gene Symbol`))
dim(conGeneSymbol1) # 18527 2
length(unique(conGeneSymbol1$`Gene Symbol`)) # 10645
SYMBS1 = unique(conGeneSymbol1$`Gene Symbol`)
GEO1file = filter(GENES, GENSYMBOL %in% SYMBS1) # 
dim(GEO1file) # 9628 2
write.csv(GEO1file,paste(DIR, "GEO1.csv",sep=""))

## Cuando FC<0, sub expresados en PTC
## Cuando FC>0, sobre expresados en PTC
## Separarlos

iSOBREgeo1 = which(tT1$logFC>0) # sobre-expr en PTC: 11105 probes
probesGEO1sobre = rownames(e1[iSOBREgeo1,])

genSOBRE1 = e1[probesGEO1sobre[1],]
boxplot(genSOBRE1~fc1,  ylab = "Expression Level", xlab="Control (0) Vs. PTC (1)", main = paste('GEO1: ',probesGEO1sobre[1] ), cex=0.5)

iSUBgeo1 = which(tT1$logFC<0) # under-expr en PTC: 7422 genes
probesGEO1sub = rownames(e1[iSUBgeo1,])

genSUB1 = e1[probesGEO1sub[1],]
boxplot(genSUB1~fc1,  ylab = "Expression Level", xlab="Control (0) Vs. PTC (1)", main = paste('GEO1: ',probesGEO1sub[1] ), cex=0.5)


A1sobre = filter(geno1, ID %in% probesGEO1sobre)
A1sobre = subset(A1sobre, select=c(ID,`Gene Symbol`))
dim(A1sobre) # 11105 2
length(unique(A1sobre$`Gene Symbol`)) # 7749
B1sobre = unique(A1sobre$`Gene Symbol`)
B1sobre = filter(GENES, GENSYMBOL %in% B1sobre) # 
dim(B1sobre) # 6565 2
B1sobre=cbind(B1sobre, UpDown="UpPTC")

A1sub = filter(geno1, ID %in% probesGEO1sub)
A1sub = subset(A1sub, select=c(ID,`Gene Symbol`))
dim(A1sub) # 7422 2
length(unique(A1sub$`Gene Symbol`)) # 5730
B1sub = unique(A1sub$`Gene Symbol`)
B1sub = filter(GENES, GENSYMBOL %in% B1sub) # 
dim(B1sub) # 4957 2
B1sub = cbind(B1sub, UpDown="DownPTC")

GEO1completo=rbind(B1sobre,B1sub)
write.csv(GEO1completo,paste(DIR, "GEO1b.csv",sep=""))

write.csv(e1, paste(DIR, "e1.csv"))
save(e1, file = "e1.RData")