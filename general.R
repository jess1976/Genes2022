##### Author: JESSICA CARBALLIDO
##### Bahía Blanca, Buenos Aires, Argentina
##### 2021
##### Differential Expression: An overrated analysis?
##### 
##### Expression Data (XENA BROWSER) - counts (raw), phenotype (samples) and genes
##### https://xenabrowser.net/datapages/
##### 
##### Cancer types: Thyroid cancer, Prostate cancer, Non-Small cell Lung cancer and Breast cancer
##### 
##### Genesets annotated to each type of cancer:
#####   GENESET KEGG HUMAN 2021 (enrichr) 
##### 
##### Imported using IMPORT DATASET from the RStudio environment

#### Genes' information: shared in all the cases of study
#### GENES #################################################################################################
gencode.v22.annotation.gene <- read.delim("~/Desktop/DATOS/TCGA.BRCA/gencode.v22.annotation.gene.probeMap")
GENES=data.frame(gencode.v22.annotation.gene$id, gencode.v22.annotation.gene$gene)
rm(gencode.v22.annotation.gene)

colnames(GENES)=c("ENSID","GENSYMBOL")
rownames(GENES)=GENES[,1]
save(GENES, file = "GENES.RData")

#### MUESTRAS ###############################################################################################

# ---- muestrasTH --------------------------------------------------------------
# Thyroid
TCGA.THCA.GDC_phenotype <- read.delim("~/Desktop/DATOS/TCGA.THCA/TCGA-THCA.GDC_phenotype.tsv")
fenoTH=TCGA.THCA.GDC_phenotype
rm(TCGA.THCA.GDC_phenotype)

length(fenoTH$sample_type.samples) #615
summary(as.factor(fenoTH$sample_type.samples))
sum(fenoTH$sample_type.samples=="Primary Tumor")
sum(fenoTH$sample_type.samples=="Solid Tissue Normal")

muestrasTH = data.frame(fenoTH$submitter_id.samples,as.factor(fenoTH$sample_type.samples))
head(muestrasTH)
dim(muestrasTH)
colnames(muestrasTH)=c("SAMPLE","TYPE")
rownames(muestrasTH)=muestrasTH$SAMPLE
head(muestrasTH)

# Dejo solamente TUMOR Y NORMAL 
muestrasTH = subset(muestrasTH, TYPE=="Primary Tumor" | TYPE=="Solid Tissue Normal")
dim(muestrasTH) # quedan 607 muestras - 100 normal, 507 tumor

muestrasTH$TYPE=as.factor(muestrasTH$TYPE) 
muestrasTH$TYPE=droplevels(muestrasTH$TYPE)
muestrasTH$SAMPLE=as.character(muestrasTH$SAMPLE)
head(muestrasTH)
class(muestrasTH$TYPE)
summary(muestrasTH$TYPE)
sum(summary(muestrasTH$TYPE))

## FALTA VER CUALES DE ESTAS ESTAN REALMENTE EN LA MATRIZ DE EXPRESION

# ---- muestrasBR --------------------------------------------------------------
# Breast
TCGA.BRCA.GDC_phenotype.tsv <- read.delim("~/Desktop/DATOS/TCGA.BRCA/TCGA-BRCA.GDC_phenotype.tsv.gz")
fenoBR=TCGA.BRCA.GDC_phenotype.tsv
rm(TCGA.BRCA.GDC_phenotype.tsv)

length(fenoBR$sample_type.samples) # 1284
summary(as.factor(fenoBR$sample_type.samples))
sum(fenoBR$sample_type.samples=="Primary Tumor")
sum(fenoBR$sample_type.samples=="Solid Tissue Normal")

muestrasBR = data.frame(fenoBR$submitter_id.samples,as.factor(fenoBR$sample_type.samples))
head(muestrasBR)
dim(muestrasBR)
colnames(muestrasBR)=c("SAMPLE","TYPE")
rownames(muestrasBR)=muestrasBR$SAMPLE
head(muestrasBR)

# Dejo solamente TUMOR Y NORMAL 
muestrasBR = subset(muestrasBR, TYPE=="Primary Tumor" | TYPE=="Solid Tissue Normal")
dim(muestrasBR) # quedan 1276 muestras - 162 normal, 1114 tumor

muestrasBR$TYPE=as.factor(muestrasBR$TYPE) 
muestrasBR$TYPE=droplevels(muestrasBR$TYPE)
muestrasBR$SAMPLE=as.character(muestrasBR$SAMPLE)
head(muestrasBR)
class(muestrasBR$TYPE)
summary(muestrasBR$TYPE)
sum(summary(muestrasBR$TYPE))

## FALTA VER CUALES DE ESTAS ESTAN REALMENTE EN LA MATRIZ DE EXPRESION

# ---- muestrasLU --------------------------------------------------------------
# Lung
TCGA.LUAD.GDC_phenotype.tsv <- read.delim("~/Desktop/DATOS/TCGA.LUAD/TCGA-LUAD.GDC_phenotype.tsv.gz")
fenoLU=TCGA.LUAD.GDC_phenotype.tsv
rm(TCGA.LUAD.GDC_phenotype.tsv)

length(fenoLU$sample_type.samples) # 877
summary(as.factor(fenoLU$sample_type.samples))
sum(fenoLU$sample_type.samples=="Primary Tumor")
sum(fenoLU$sample_type.samples=="Solid Tissue Normal")

muestrasLU = data.frame(fenoLU$submitter_id.samples,as.factor(fenoLU$sample_type.samples))
head(muestrasLU)
dim(muestrasLU)
colnames(muestrasLU)=c("SAMPLE","TYPE")
rownames(muestrasLU)=muestrasLU$SAMPLE
head(muestrasLU)
# Dejo solamente TUMOR Y NORMAL 
muestrasLU = subset(muestrasLU, TYPE=="Primary Tumor" | TYPE=="Solid Tissue Normal")
dim(muestrasLU) # quedan 871 muestras - 274 normal, 597 tumor

muestrasLU$TYPE=as.factor(muestrasLU$TYPE) 
muestrasLU$TYPE=droplevels(muestrasLU$TYPE)
muestrasLU$SAMPLE=as.character(muestrasLU$SAMPLE)
head(muestrasLU)
class(muestrasLU$TYPE)
summary(muestrasLU$TYPE)
sum(summary(muestrasLU$TYPE))

## FALTA VER CUALES DE ESTAS ESTAN REALMENTE EN LA MATRIZ DE EXPRESION


# ---- muestrasPR --------------------------------------------------------------
# Prostate
TCGA.PRAD.GDC_phenotype.tsv <- read.delim("~/Desktop/DATOS/TCGA.PRAD/TCGA-PRAD.GDC_phenotype.tsv.gz")
fenoPR=TCGA.PRAD.GDC_phenotype.tsv
rm(TCGA.PRAD.GDC_phenotype.tsv)

length(fenoPR$sample_type.samples) # 623
summary(as.factor(fenoPR$sample_type.samples))
sum(fenoPR$sample_type.samples=="Primary Tumor")
sum(fenoPR$sample_type.samples=="Solid Tissue Normal")

muestrasPR = data.frame(fenoPR$submitter_id.samples,as.factor(fenoPR$sample_type.samples))
head(muestrasPR)
dim(muestrasPR)
colnames(muestrasPR)=c("SAMPLE","TYPE")
rownames(muestrasPR)=muestrasPR$SAMPLE
head(muestrasPR)
# Dejo solamente TUMOR Y NORMAL 
muestrasPR = subset(muestrasPR, TYPE=="Primary Tumor" | TYPE=="Solid Tissue Normal")
dim(muestrasPR) # quedan 871 muestras - 274 normal, 597 tumor

muestrasPR$TYPE=as.factor(muestrasPR$TYPE) 
muestrasPR$TYPE=droplevels(muestrasPR$TYPE)
muestrasPR$SAMPLE=as.character(muestrasPR$SAMPLE)
head(muestrasPR)
class(muestrasPR$TYPE)
summary(muestrasPR$TYPE)
sum(summary(muestrasPR$TYPE))

## FALTA VER CUALES DE ESTAS ESTAN REALMENTE EN LA MATRIZ DE EXPRESION

#### COUNTS ###############################################################################################
# Descarga matriz de conteos  --- BAJAR RAW DATA - FPKM NO SIRVE PARA EXPR. DIF.
# 
#  OJO!!!!!!!!!!!! NO USAR FPKM PARA EXPRESION DIFERENCIAL
#  USAR RAW COUNTS
#  IMPORTAR USANDO: Import Dataset de pesta?a ENVIRNMENT - RSTUDIO
#  Data from the same sample but from different vials/portions/analytes/aliquotes is averaged; data from different samples is combined into genomicMatrix; all data is then log2(x+1) transformed.
# 

# ---- conteosTH ----

TCGA.THCA.htseq_counts.tsv <- read.delim("~/Desktop/DATOS/TCGA.THCA/TCGA-THCA.htseq_counts.tsv.gz")
CONTEOS1=TCGA.THCA.htseq_counts.tsv
rm(TCGA.THCA.htseq_counts.tsv)
dim(CONTEOS1) # hay 569 muestras, son menos q las de muestrasTH (607 tumor+normal), hay q tomar solo estas
conteosTH=CONTEOS1

# Los nombres de las muestras estan en formato: TCGA.J8.A3YE.01A
# Hay q ponerlos en formato: TCGA-BJ-A28W-01A --- reemplazamos . por -
# Ojo que la columna 1 viene con los Ensembl_ID
head(conteosTH)[1:5]
aux = colnames(conteosTH)[2:ncol(conteosTH)]
aux = gsub(".","-", aux, fixed=T)
colnames(conteosTH)[2:ncol(conteosTH)]=aux
head(conteosTH)[1:5]

# # MODIFICO MUESTRAS dejando solo las q aparecen en la matriz de expresion
muestrasTH=muestrasTH[colnames(conteosTH[,2:ncol(conteosTH)]),]
dim(muestrasTH)
summary(muestrasTH)
## son 58 normal y 502 tumor.
## aprovechamos a balancear los niveles del factor, y dejamos 58 en cada clase
normalesTH=which(muestrasTH$TYPE=="Solid Tissue Normal")
tumoresTH=which(muestrasTH$TYPE=="Primary Tumor")
## selecciono 50 subconjuntos de muestras de tumores, con la long de muestras normales
tumoresTH.mat=matrix(data=NA, nrow=50, ncol=length(normalesTH))
for (i in 1:50)
{
  tumoresTH.mat[i,]=sample(tumoresTH, length(normalesTH), replace = FALSE)
}

length(normalesTH)
dim(tumoresTH.mat)

# Ahora, uso la primera columna de la matriz CONTEOS original para nombrar a las filas,
# y la elimino
rownames(conteosTH)=CONTEOS1[,1]
head(conteosTH)[1:5]
conteosTH=conteosTH[,-1]
dim(conteosTH)

### Ahora filtro CONTEOS para quedarme solo con los genes del geneset KEGG

kegg2021TH <- read.table("~/Desktop/DATOS/TCGA.THCA/GENESETS/New folder/kegg2021.txt", quote="\"", comment.char="")
dim(kegg2021TH)
L1=as.vector(kegg2021TH[,1])
class(L1)
length(L1) # 37

cuales=which(GENES$GENSYMBOL %in% L1)
ens=rownames(GENES[cuales,])

conteosTH = conteosTH[ens,]
dim(conteosTH)

## genero lista de 50 matrices de conteos: cada nueva matriz de conteos tiene 37 genes y 116 muestras (58 normal, 58 tumor q van cambiando aprovechando las muestras extra de tumor q habia)
conteosTH.list= list()
for (j in 1:50)
{ 
  seleccionMuestrasTH=c(rownames(muestrasTH[normalesTH,]),rownames(muestrasTH[tumoresTH.mat[j,],]))
  conteosTH.list[[j]]=as.matrix(conteosTH[,seleccionMuestrasTH])
} 
length(conteosTH.list)

matExprTH.list = conteosTH.list
length(matExprTH.list)
rm(conteosTH.list)

# genero el factor de clases q indica q las primeras 58 muestras son sanas y las ult son enfermas
FCth = muestrasTH[colnames(matExprTH.list[[1]]),"TYPE"]
summary(FCth)

# ---- conteosBR ----
TCGA.BRCA.htseq_counts.tsv <- read.delim("~/Desktop/DATOS/TCGA.BRCA/TCGA-BRCA.htseq_counts.tsv.gz")
CONTEOS2=TCGA.BRCA.htseq_counts.tsv
rm(TCGA.BRCA.htseq_counts.tsv)
dim(CONTEOS2) # hay 1217 muestras + 1 columna ENSEMBLID, son menos q las de muestrasBR (1284 tumor+normal), hay q tomar solo estas
conteosBR=CONTEOS2

# Los nombres de las muestras estan en formato: TCGA.J8.A3YE.01A
# Hay q ponerlos en formato: TCGA-BJ-A28W-01A --- reemplazamos . por -
# Ojo que la columna 1 viene con los Ensembl_ID
head(conteosBR)[1:5]
aux = colnames(conteosBR)[2:ncol(conteosBR)]
aux = gsub(".","-", aux, fixed=T)
colnames(conteosBR)[2:ncol(conteosBR)]=aux
head(conteosBR)[1:5]

# # Filtro COLUMNAS: MODIFICO MUESTRAS dejando solo las q aparecen en la matriz de expresion
muestrasBR=muestrasBR[colnames(conteosBR[,2:ncol(conteosBR)]),]
dim(muestrasBR) # 1217    2
summary(muestrasBR)
## son 113 normal y 1097 tumor.
 
## aprovechamos a balancear los niveles del factor, y dejamos 113 en cada clase
normalesBR=which(muestrasBR$TYPE=="Solid Tissue Normal")
tumoresBR=which(muestrasBR$TYPE=="Primary Tumor")

## selecciono 50 subconjuntos de muestras de tumores, con la long de muestras normales
tumoresBR.mat=matrix(data=NA, nrow=50, ncol=length(normalesBR))
for (i in 1:50)
{
  tumoresBR.mat[i,]=sample(tumoresBR, length(normalesBR), replace = FALSE)
}

length(normalesBR)
dim(tumoresBR.mat)

# Ahora, uso la primera columna de la matriz CONTEOS original para nombrar a las filas
# y la elimino
rownames(conteosBR)=CONTEOS2[,1]
head(conteosBR)[1:5]
conteosBR=conteosBR[,-1]
dim(conteosBR)

### Ahora filtro filas de CONTEOS para quedarme solo con los genes del geneset KEGG
kegg2021BR = read.table("~/Desktop/DATOS/TCGA.BRCA/GENESETS/New folder/kegg2021.txt", quote="\"", comment.char="")
dim(kegg2021BR)
L1=as.vector(kegg2021BR[,1])
class(L1)
length(L1) # 147

cuales=which(GENES$GENSYMBOL %in% L1)
ens=rownames(GENES[cuales,])

conteosBR = conteosBR[ens,]
dim(conteosBR)

## genero lista de 50 matrices de conteos: cada nueva matriz de conteos tiene 148 genes y 226 muestras (113 normal, 113 tumor q van cambiando aprovechando las muestras extra de tumor q habia)
conteosBR.list= list()
for (j in 1:50)
{ 
  seleccionMuestrasBR=c(rownames(muestrasBR[normalesBR,]),rownames(muestrasBR[tumoresBR.mat[j,],]))
  conteosBR.list[[j]]=as.matrix(conteosBR[,seleccionMuestrasBR])
} 
length(conteosBR.list)

matExprBR.list = conteosBR.list
length(matExprBR.list)
rm(conteosBR.list)

FCbr = muestrasBR[colnames(matExprBR.list[[1]]),"TYPE"]
summary(FCbr)

# ---- conteosLU ----

TCGA.LUAD.htseq_counts.tsv <- read.delim("~/Desktop/DATOS/TCGA.LUAD/TCGA-LUAD.htseq_counts.tsv.gz")
CONTEOS3 = TCGA.LUAD.htseq_counts.tsv
rm(TCGA.LUAD.htseq_counts.tsv)
dim(CONTEOS3) # hay 585 muestras + 1 columna ENSEMBLID, son menos q las de muestrasBR (1284 tumor+normal), hay q tomar solo estas
conteosLU=CONTEOS3

# Los nombres de las muestras estan en formato: TCGA.J8.A3YE.01A
# Hay q ponerlos en formato: TCGA-BJ-A28W-01A --- reemplazamos . por -
# Ojo que la columna 1 viene con los Ensembl_ID
head(conteosLU)[1:5]
aux = colnames(conteosLU)[2:ncol(conteosLU)]
aux = gsub(".","-", aux, fixed=T)
colnames(conteosLU)[2:ncol(conteosLU)]=aux
head(conteosLU)[1:5]

# # Filtro COLUMNAS: MODIFICO MUESTRAS dejando solo las q aparecen en la matriz de expresion
muestrasLU=muestrasLU[colnames(conteosLU[,2:ncol(conteosLU)]),]
dim(muestrasLU) # 585    2
summary(muestrasLU)
## son 59 normal y 524 tumor (2 NA).
## aprovechamos a balancear los niveles del factor, y dejamos 59 en cada clase
normalesLU=which(muestrasLU$TYPE=="Solid Tissue Normal")
tumoresLU=which(muestrasLU$TYPE=="Primary Tumor")

## selecciono 50 subconjuntos de muestras de tumores, con la long de muestras normales
tumoresLU.mat=matrix(data=NA, nrow=50, ncol=length(normalesLU))
for (i in 1:50)
{
  tumoresLU.mat[i,]=sample(tumoresLU, length(normalesLU), replace = FALSE)
}

length(normalesLU)
dim(tumoresLU.mat)
any(tumoresLU.mat=="NA")


# Ahora, uso la primera columna de la matriz CONTEOS original para nombrar a las filas
# y la elimino
rownames(conteosLU)=CONTEOS3[,1]
head(conteosLU)[1:5]
conteosLU=conteosLU[,-1]
dim(conteosLU)


### Ahora filtro filas de CONTEOS para quedarme solo con los genes del geneset KEGG

kegg2021LU= read.table("~/Desktop/DATOS/TCGA.LUAD/GENESETS/New folder/KEGG2021.txt", quote="\"", comment.char="")  
dim(kegg2021LU)
L1=as.vector(kegg2021LU[,1])
class(L1)
length(L1) # 72

cuales=which(GENES$GENSYMBOL %in% L1)
ens=rownames(GENES[cuales,])

conteosLU = conteosLU[ens,]
dim(conteosLU)

## genero lista de 50 matrices de conteos: cada nueva matriz de conteos tiene  73 genes (un gen tiene 2 sondas) y 118 muestras (59 normal, 59 tumor q van cambiando aprovechando las muestras extra de tumor q habia)
conteosLU.list= list()
for (j in 1:50)
{ 
  seleccionMuestrasLU=c(rownames(muestrasLU[normalesLU,]),rownames(muestrasLU[tumoresLU.mat[j,],]))
  conteosLU.list[[j]]=as.matrix(conteosLU[,seleccionMuestrasLU])
} 
length(conteosLU.list)

matExprLU.list = conteosLU.list
length(matExprLU.list)
rm(conteosLU.list)


FClu = muestrasLU[colnames(matExprLU.list[[1]]),"TYPE"]
summary(FClu)

# ---- conteosPR ----
TCGA.PRAD.htseq_counts.tsv <- read.delim("~/Desktop/DATOS/TCGA.PRAD/TCGA-PRAD.htseq_counts.tsv.gz")
CONTEOS4 = TCGA.PRAD.htseq_counts.tsv
rm(TCGA.PRAD.htseq_counts.tsv)
dim(CONTEOS4) # hay 551 muestras + 1 columna ENSEMBLID, son menos q las de muestrasPR 
conteosPR=CONTEOS4

# Los nombres de las muestras estan en formato: TCGA.J8.A3YE.01A
# Hay q ponerlos en formato: TCGA-BJ-A28W-01A --- reemplazamos . por -
# Ojo que la columna 1 viene con los Ensembl_ID
head(conteosPR)[1:5]
aux = colnames(conteosPR)[2:ncol(conteosPR)]
aux = gsub(".","-", aux, fixed=T)
colnames(conteosPR)[2:ncol(conteosPR)]=aux
head(conteosPR)[1:5]

# # Filtro COLUMNAS: MODIFICO MUESTRAS dejando solo las q aparecen en la matriz de expresion
muestrasPR=muestrasPR[colnames(conteosPR[,2:ncol(conteosPR)]),]
dim(muestrasPR) # 551    2
summary(muestrasPR)
## son 52 normal y 498 tumor.

## aprovechamos a balancear los niveles del factor, y dejamos 52 en cada clase
normalesPR = which(muestrasPR$TYPE=="Solid Tissue Normal")
tumoresPR = which(muestrasPR$TYPE=="Primary Tumor")

## selecciono 50 subconjuntos de muestras de tumores, con la long de muestras normales
tumoresPR.mat=matrix(data=NA, nrow=50, ncol=length(normalesPR))
for (i in 1:50)
{
  tumoresPR.mat[i,]=sample(tumoresPR, length(normalesPR), replace = FALSE)
}

length(normalesPR)
dim(tumoresPR.mat)
any(tumoresPR.mat=="NA")


# Ahora, uso la primera columna de la matriz CONTEOS original para nombrar a las filas
# y la elimino
rownames(conteosPR)=CONTEOS4[,1]
head(conteosPR)[1:5]
conteosPR=conteosPR[,-1]
dim(conteosPR)

### Ahora filtro filas de CONTEOS para quedarme solo con los genes del geneset KEGG

kegg2021PR <- read.table("~/Desktop/DATOS/TCGA.PRAD/GENESETS/New folder/KEGG2021.txt", quote="\"", comment.char="")
dim(kegg2021PR)
L1=as.vector(kegg2021PR[,1])
class(L1)
length(L1) # 97

cuales=which(GENES$GENSYMBOL %in% L1)
ens=rownames(GENES[cuales,])

conteosPR = conteosPR[ens,]
dim(conteosPR)

## genero lista de 50 matrices de conteos: cada nueva matriz de conteos tiene 98 genes (dos sondas para un gen) y 104 muestras (52 normal, 52 tumor q van cambiando aprovechando las muestras extra de tumor q habia)
conteosPR.list= list()
for (j in 1:50)
{ 
  seleccionMuestrasPR=c(rownames(muestrasPR[normalesPR,]),rownames(muestrasPR[tumoresPR.mat[j,],]))
  conteosPR.list[[j]]=as.matrix(conteosPR[,seleccionMuestrasPR])
} 
length(conteosPR.list)

matExprPR.list = conteosPR.list
length(matExprPR.list)
rm(conteosPR.list)

FCpr = muestrasPR[colnames(matExprPR.list[[1]]),"TYPE"]
summary(FCpr)



### Differential Expression Analysis ----
## 
# BiocManager::install("edgeR")
# library(edgeR)

# ---- DE TH ----
 
disenio=model.matrix(~0+FCth)
colnames(disenio)= levels(FCth)
head(disenio)
apply(disenio, 2, sum)

## ARMO LISTA DE GENES DIF EXPR en cada matriz de expresion de la lista en MatExprTH.list
noTH.list=list()
siTH.list=list()
for (i in 1:50)
{ 
  matExprTH=matExprTH.list[[i]]
  dge = DGEList(counts=matExprTH, group=FCth)
  dge <- calcNormFactors(dge)
  dge <- estimateGLMCommonDisp(dge, disenio)
  dge <- estimateGLMTagwiseDisp(dge, disenio)
  
  etest=exactTest(dge)
#  View(etest[["table"]])
  ttTH <-topTags(etest, n=nrow(matExprTH), adjust.method = "none") 
  # si no ajusto pval, toma los q son pval menor a 0.5
  # es una restriccion sobre el FDR (p value ajustado si pongo restriccion en pval)
  
  dim(ttTH)
  
  tablaEDif = etest[["table"]]
  tablaORDth <- tablaEDif[order(tablaEDif$PValue), ]

  SIth=rownames(tablaORDth[which(tablaORDth$PValue<=0.05),])
  NOth=rownames(tablaORDth[which(tablaORDth$PValue>0.05),])
  
  siTH.list[[i]]=SIth
  noTH.list[[i]]=NOth
}
# vector SIEMPREde y NUNCAde, con la interseccion de los q fueron DE y la interseccion de los no DE respect.
SIEMPREde=siTH.list[[1]]
NUNCAde=noTH.list[[1]]
for (i in 2:50)
{ 
  SIEMPREde=intersect(SIEMPREde, siTH.list[[i]])
  NUNCAde=intersect(NUNCAde, noTH.list[[i]])
  }
length(SIEMPREde) # 3
length(NUNCAde) # 27

# BOXPLOTS: uno de SIEMPREde y uno de NUNCAde
unDE=SIEMPREde[1]
unNODE=NUNCAde[length(NUNCAde)]
nombreDE=GENES[unDE, 2]
nombreNODE=GENES[unNODE, 2]
par(mfrow=c(1,2))
boxplot(matExprTH[unDE,]~FCth, ylim=c(0,15), ylab = "Expression Level", xlab="", main='KEGG Human 2021', sub= nombreDE, cex=0.5)
# peor
boxplot(matExprTH[unNODE,]~FCth, ylim=c(0,15), ylab = "Expression Level", xlab="", main='Thyroid Cancer', sub= nombreNODE, cex=0.5)

write.csv(SIEMPREde, "deTH.csv")
write.csv(NUNCAde, "NOdeTH.csv")
SIEMPREdeTH=SIEMPREde
NUNCAdeTH=NUNCAde

# ---- DE BR ----

disenio=model.matrix(~0+FCbr)
colnames(disenio)= levels(FCbr)
head(disenio)
apply(disenio, 2, sum)

## ARMO LISTA DE GENES DIF EXPR en cada matriz de expresion de la lista en MatExprBR.list
noBR.list=list()
siBR.list=list()
for (i in 1:50)
{ 
  matExprBR=matExprBR.list[[i]]
  dge = DGEList(counts=matExprBR, group=FCbr)
  dge <- calcNormFactors(dge)
  dge <- estimateGLMCommonDisp(dge, disenio)
  dge <- estimateGLMTagwiseDisp(dge, disenio)
  
  etest=exactTest(dge)
  #  View(etest[["table"]])
  ttBR <-topTags(etest, n=nrow(matExprBR), adjust.method = "none") 
  # si no ajusto pval, toma los q son pval menor a 0.5
  # es una restriccion sobre el FDR (p value ajustado si pongo restriccion en pval)
  
  dim(ttBR)
  
  tablaEDif = etest[["table"]]
  tablaORDbr <- tablaEDif[order(tablaEDif$PValue), ]

  SIbr=rownames(tablaORDbr[which(tablaORDbr$PValue<=0.05),])
  NObr=rownames(tablaORDbr[which(tablaORDbr$PValue>0.05),])
  
  siBR.list[[i]]=SIbr
  noBR.list[[i]]=NObr
}
# vector SIEMPREde y NUNCAde, con la interseccion de los q fueron DE y la interseccion de los no DE respect.
SIEMPREde=siBR.list[[1]]
NUNCAde=noBR.list[[1]]
for (i in 2:50)
{ 
  SIEMPREde=intersect(SIEMPREde, siBR.list[[i]])
  NUNCAde=intersect(NUNCAde, noBR.list[[i]])
}
length(SIEMPREde) #
length(NUNCAde) # 

# BOXPLOTS: uno de SIEMPREde y uno de NUNCAde
unDE=SIEMPREde[1]
unNODE=NUNCAde[length(NUNCAde)]
nombreDE=GENES[unDE, 2]
nombreNODE=GENES[unNODE, 2]
par(mfrow=c(1,2))
boxplot(matExprBR[unDE,]~FCbr, ylim=c(0,15), ylab = "Expression Level", xlab="", main='KEGG Human 2021', sub= nombreDE, cex=0.5)
# peor
boxplot(matExprBR[unNODE,]~FCbr, ylim=c(0,15), ylab = "Expression Level", xlab="", main='Breast Cancer', sub= nombreNODE, cex=0.5)

write.csv(SIEMPREde, "deBR.csv")
write.csv(NUNCAde, "NOdeBR.csv")
SIEMPREdeBR=SIEMPREde
NUNCAdeBR=NUNCAde

# ---- DE LU ----

disenio=model.matrix(~0+FClu)
colnames(disenio)= levels(FClu)
head(disenio)
apply(disenio, 2, sum)

## ARMO LISTA DE GENES DIF EXPR en cada matriz de expresion de la lista en MatExprLU.list
noLU.list=list()
siLU.list=list()
for (i in 1:50)
{ 
  matExprLU=matExprLU.list[[i]]
  dge = DGEList(counts=matExprLU, group=FClu)
  dge <- calcNormFactors(dge)
  dge <- estimateGLMCommonDisp(dge, disenio)
  dge <- estimateGLMTagwiseDisp(dge, disenio)
  
  etest=exactTest(dge)
  #  View(etest[["table"]])
  ttLU <-topTags(etest, n=nrow(matExprLU), adjust.method = "none") 
  # si no ajusto pval, toma los q son pval menor a 0.5
  # es una restriccion sobre el FDR (p value ajustado si pongo restriccion en pval)
  
  dim(ttLU)
  
  tablaEDif = etest[["table"]]
  tablaORDlu <- tablaEDif[order(tablaEDif$PValue), ]
  
  SIlu=rownames(tablaORDlu[which(tablaORDlu$PValue<=0.05),])
  NOlu=rownames(tablaORDlu[which(tablaORDlu$PValue>0.05),])
  
  siLU.list[[i]]=SIlu
  noLU.list[[i]]=NOlu
}
# vector SIEMPREde y NUNCAde, con la interseccion de los q fueron DE y la interseccion de los no DE respect.
SIEMPREde=siLU.list[[1]]
NUNCAde=noLU.list[[1]]
for (i in 2:50)
{ 
  SIEMPREde=intersect(SIEMPREde, siLU.list[[i]])
  NUNCAde=intersect(NUNCAde, noLU.list[[i]])
}
length(SIEMPREde) #
length(NUNCAde) # 

# BOXPLOTS: uno de SIEMPREde y uno de NUNCAde
unDE=SIEMPREde[1]
unNODE=NUNCAde[length(NUNCAde)]
nombreDE=GENES[unDE, 2]
nombreNODE=GENES[unNODE, 2]
par(mfrow=c(1,2))
boxplot(matExprLU[unDE,]~FClu, ylim=c(0,15), ylab = "Expression Level", xlab="", main='KEGG Human 2021', sub= nombreDE, cex=0.5)
# peor
boxplot(matExprLU[unNODE,]~FClu, ylim=c(0,15), ylab = "Expression Level", xlab="", main='Lung Cancer', sub= nombreNODE, cex=0.5)

write.csv(SIEMPREde, "deLU.csv")
write.csv(NUNCAde, "NOdeLU.csv")
SIEMPREdeLU=SIEMPREde
NUNCAdeLU=NUNCAde

# ---- DE PR ----

disenio=model.matrix(~0+FCpr)
colnames(disenio)= levels(FCpr)
head(disenio)
apply(disenio, 2, sum)

## ARMO LISTA DE GENES DIF EXPR en cada matriz de expresion de la lista en MatExprLU.list
noPR.list=list()
siPR.list=list()
for (i in 1:50)
{ 
  matExprPR=matExprPR.list[[i]]
  dge = DGEList(counts=matExprPR, group=FCpr)
  dge <- calcNormFactors(dge)
  dge <- estimateGLMCommonDisp(dge, disenio)
  dge <- estimateGLMTagwiseDisp(dge, disenio)
  
  etest=exactTest(dge)
  #  View(etest[["table"]])
  ttPR <-topTags(etest, n=nrow(matExprPR), adjust.method = "none") 
  # si no ajusto pval, toma los q son pval menor a 0.5
  # es una restriccion sobre el FDR (p value ajustado si pongo restriccion en pval)
  
  dim(ttPR)
  
  tablaEDif = etest[["table"]]
  tablaORDpr <- tablaEDif[order(tablaEDif$PValue), ]
  
  SIpr=rownames(tablaORDpr[which(tablaORDpr$PValue<=0.05),])
  NOpr=rownames(tablaORDpr[which(tablaORDpr$PValue>0.05),])
  
  siPR.list[[i]]=SIpr
  noPR.list[[i]]=NOpr
}
# vector SIEMPREde y NUNCAde, con la interseccion de los q fueron DE y la interseccion de los no DE respect.
SIEMPREde=siPR.list[[1]]
NUNCAde=noPR.list[[1]]
for (i in 2:50)
{ 
  SIEMPREde=intersect(SIEMPREde, siPR.list[[i]])
  NUNCAde=intersect(NUNCAde, noPR.list[[i]])
}
length(SIEMPREde) #
length(NUNCAde) # 

# BOXPLOTS: uno de SIEMPREde y uno de NUNCAde
unDE=SIEMPREde[1]
unNODE=NUNCAde[length(NUNCAde)]
nombreDE=GENES[unDE, 2]
nombreNODE=GENES[unNODE, 2]
par(mfrow=c(1,2))
boxplot(matExprPR[unDE,]~FCpr, ylim=c(0,15), ylab = "Expression Level", xlab="", main='KEGG Human 2021', sub= nombreDE, cex=0.5)
# peor
boxplot(matExprPR[unNODE,]~FCpr, ylim=c(0,15), ylab = "Expression Level", xlab="", main='Prostate Cancer', sub= nombreNODE, cex=0.5)

write.csv(SIEMPREde, "dePR.csv")
write.csv(NUNCAde, "NOdePR.csv")
SIEMPREdePR=SIEMPREde
NUNCAdePR=NUNCAde

#### ---- Resultados Regresión Logística y Clasificación Random Forest
library(caret)

### Logistic Regression and Classification ----
# BiocManager::install("ggpubr")
library("ggpubr")
library("effsize") # para fc. cohen.d (tamaño del efecto)
library(caret)


# ---- GLM TH ---- 
## Datos
dim(matExprTH.list[[1]])
length(FCth)
summary(FCth)
maxPRED = min(floor(length(FCth)/10), length(SIEMPREdeTH)) 
predSI = SIEMPREdeTH[1:maxPRED] 
 
## Regresion Logistica, un modelo para cada matriz en matExprTH.list variando los genes NOde (ya 
## que tenemos mas noDE que DE)
## y voy registrando AIC y log.lik
## 
regTH=data.frame()
for (i in 1:50)
{
  cambianNOde = sample(1:length(NUNCAdeTH),maxPRED, replace = FALSE)
  predNO = NUNCAdeTH[cambianNOde]
  datos = data.frame(t(matExprTH.list[[i]]), FCth)
  subSI=datos[,predSI]
  subNO= datos[,predNO]
  
  modelo.DE <- glm(FCth ~ ., data = subSI, family = binomial(logit), maxit=100)
  modelo.noDE <- glm(FCth ~ ., data = subNO, family = binomial(logit), maxit=100)
  de=glance(modelo.DE)
  no.de=glance(modelo.noDE)
  
  regTH = rbind(regTH, c(AIC.DE = de$AIC , AIC.notDE = no.de$AIC , log.lik.DE = de$logLik , log.lik.notDE = no.de$logLik))
  
}
colnames(regTH) = c("AIC.DE", "AIC.notDE", "log.lik.DE", "log.lik.notDE")
promediosTH = apply(regTH, 2, mean)
write.csv(regTH, "detalleRegresionTH.csv")
# write.csv(promediosTH, "mediasRegresionTH.csv")

# test para ver si la diferencia entre AICs es significativa
# veo si la distr es normal
# SHAPIRO-WILK: The data is normal if the p-value is above 0.05. 
shapiro.test(regTH$AIC.DE) # p-value = 0.7955 -- NORMAL
densityplot(regTH$AIC.DE)
shapiro.test(regTH$AIC.notDE) # p-value = 3.163e-05  -- NO NORMAL
densityplot(regTH$AIC.notDE)
# uso test no  parametrico 
wilcox.test(regTH$AIC.DE,regTH$AIC.notDE,alternative = "two.sided",paired=FALSE, conf.int = 0.95)
#  p-value < 2.2e-16
# existe una probabilidad menor a 0.05 de que  la dif. entre las medianas se deba al azar
# wilcox usa MEDIANAS
medianasTH=apply(regTH, 2, median)
medianasTH
cohen.d(regTH$AIC.DE,regTH$AIC.notDE)
# Cohen's d
# 
# d estimate: -5.25084 (large)
# 95 percent confidence interval:
#     lower     upper 
# -6.087750 -4.413929 



## ---- RF TH ----

  Ind <- createDataPartition(y = FCth, p = 0.7,list = FALSE)
  FCthTrain=FCth[Ind]
  FCthTest=FCth[-Ind]
  tablaACCKAP.TH = data.frame()

for (i in 1:50)
{ 
matExprTH = matExprTH.list[[i]]
  
datosTH.si = data.frame(t(matExprTH[SIEMPREdeTH,]), FCth)
datosTH.no = data.frame(t(matExprTH[NUNCAdeTH[1:length(SIEMPREdeTH)],]), FCth)
  
EntrenamientoTH.si <- datosTH.si[Ind,]
dim(EntrenamientoTH.si)
TestTH.si <- datosTH.si[-Ind,]
dim(TestTH.si)
nrow(EntrenamientoTH.si) + nrow(TestTH.si)
  
control_train <- trainControl(method = "repeatedcv", number = 5,
                                repeats = 10)
  
RF.TH.si = train(FCth ~ ., data = EntrenamientoTH.si, method = "rf", trControl = control_train, tuneLength = 2, tuneGrid = expand.grid(mtry=c(1,2)))
# print(RF.TH.si)
# RF.TH.si$results
pred.TH.si = predict(RF.TH.si, TestTH.si)
CM.TH.si = confusionMatrix(table(TestTH.si[,"FCth"], pred.TH.si))
CM.TH.si$overall[c(1,2)]

# print(varImp(RF.TH.si))
# plot(varImp(RF.TH.si))
# se sabe q FCth pertenece al data.frame EntrenamientoTH.si y es la predictora 
# mtry es la cantidad de variables q se usaran para dividir en cada rama
  
EntrenamientoTH.no <- datosTH.no[Ind,]
dim(EntrenamientoTH.no)
TestTH.no <- datosTH.no[-Ind,]
dim(TestTH.no)
nrow(EntrenamientoTH.no) + nrow(TestTH.no)
  
control_train <- trainControl(method = "repeatedcv", number = 5,
                                repeats = 10)
  
RF.TH.no = train(FCth ~ ., data = EntrenamientoTH.no, method = "rf", trControl = control_train, tuneLength = 2, tuneGrid = expand.grid(mtry=c(1,2)))
#  print(RF.TH.no)
#  RF.TH.no$results
pred.TH.no = predict(RF.TH.no, TestTH.no)
CM.TH.no = confusionMatrix(table(TestTH.no[,"FCth"], pred.TH.no))
CM.TH.no$overall[c(1,2)]

tablaACCKAP.TH = rbind(tablaACCKAP.TH, c(CM.TH.si$overall[c(1,2)],CM.TH.no$overall[c(1,2)])) 
}
  colnames(tablaACCKAP.TH)=c("ACC.TH.si", "KAPPA.TH.si", "ACC.TH.no", "KAPPA.TH.no")
  
  apply(tablaACCKAP.TH, 2, mean)
  apply(tablaACCKAP.TH, 2, median)
  test1=wilcox.test(tablaACCKAP.TH$ACC.TH.si, tablaACCKAP.TH$ACC.TH.no, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
  test1$p.value
  test2=wilcox.test(tablaACCKAP.TH$KAPPA.TH.si, tablaACCKAP.TH$KAPPA.TH.no, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
  test2$p.value
  cohen.d(tablaACCKAP.TH$ACC.TH.si, tablaACCKAP.TH$ACC.TH.no)
  cohen.d(tablaACCKAP.TH$KAPPA.TH.si, tablaACCKAP.TH$KAPPA.TH.no)
  
#  print(varImp(RF.TH.no))
#  plot(varImp(RF.TH.no))
# se sabe q FCth pertenece al data.frame EntrenamientoTH.si y es la predictora 
# mtry es la cantidad de variables q se usaran para dividir en cada rama
  
# ---- GLM BR ---- 
##Datos
dim(matExprBR)
length(FCbr)
summary(FCbr)
maxPRED = min(floor(length(FCbr)/10), length(SIEMPREdeBR)) 
predSI = SIEMPREdeBR[1:maxPRED] 

## Regresion Logistica, un modelo para cada matriz en matExprBR.list variando los genes NOde (ya 
## que tenemos mas noDE que DE)
## y voy registrando AIC y log.lik
## 
regBR=data.frame()
for (i in 1:50)
{
  cambianNOde = sample(1:length(NUNCAdeBR),maxPRED, replace = FALSE)
  predNO = NUNCAdeBR[cambianNOde]
  datos = data.frame(t(matExprBR.list[[i]]), FCbr)
  subSI=datos[,predSI]
  subNO= datos[,predNO]
  
  modelo.DE <- glm(FCbr ~ ., data = subSI, family = binomial(logit), maxit=100)
  modelo.noDE <- glm(FCbr ~ ., data = subNO, family = binomial(logit), maxit=100)
  de=glance(modelo.DE)
  no.de=glance(modelo.noDE)
  
  regBR = rbind(regBR, c(AIC.DE = de$AIC , AIC.notDE = no.de$AIC , log.lik.DE = de$logLik , log.lik.notDE = no.de$logLik))
  
}
colnames(regBR) = c("AIC.DE", "AIC.notDE", "log.lik.DE", "log.lik.notDE")
promediosBR = apply(regBR, 2, mean)
write.csv(regBR, "detalleRegresionBR.csv")
# write.csv(promediosBR, "mediasRegresionBR.csv")

# test para ver si la diferencia entre AICs es significativa
# veo si la distr es normal
# SHAPIRO-WILK: The data is normal if the p-value is above 0.05. 
shapiro.test(regBR$AIC.DE) # p-value = 0.07 -- NORMAL alfa menor a 0.05 no es normal
densityplot(regBR$AIC.DE)
shapiro.test(regBR$AIC.notDE) # p-value = 0.1339  -- NORMAL
densityplot(regBR$AIC.notDE)

# igual uso test no  parametrico 
wilcox.test(regBR$AIC.DE,regBR$AIC.notDE,alternative = "two.sided",paired=FALSE, conf.int = 0.95)
#  p-value =  1.552e-14

# existe una probabilidad menor a 0.05 de que  la dif. entre las medianas se deba al azar
# wilcox usa MEDIANAS
medianasBR=apply(regBR, 2, median)
medianasBR
cohen.d(regBR$AIC.DE,regBR$AIC.notDE)
# Cohen's d
# 
# d estimate: -2.342016 (large)
# 95 percent confidence interval:
#     lower     upper 
# -2.857310 -1.826722 

## ---- RF BR ----

  Ind <- createDataPartition(y = FCbr, p = 0.7,list = FALSE)
  FCbrTrain=FCbr[Ind]
  FCbrTest=FCbr[-Ind]
  tablaACCKAP.BR = data.frame()
  
for (i in 1:50)
{
matExprBR = matExprBR.list[[i]]

datosBR.si = data.frame(t(matExprBR[SIEMPREdeBR,]), FCbr)
datosBR.no = data.frame(t(matExprBR[NUNCAdeBR[1:length(SIEMPREdeBR)],]), FCbr)

EntrenamientoBR.si <- datosBR.si[Ind,]
dim(EntrenamientoBR.si)
TestBR.si <- datosBR.si[-Ind,]
dim(TestBR.si)
nrow(EntrenamientoBR.si) + nrow(TestBR.si)

control_train <- trainControl(method = "repeatedcv", number = 5,
                              repeats = 10)


RF.BR.si = train(FCbr ~ ., data = EntrenamientoBR.si, method = "rf", trControl = control_train, tuneLength = 2, tuneGrid = expand.grid(mtry=c(2,4, 10, 20)))
# print(RF.BR.si)
# RF.BR.si$results
pred.BR.si = predict(RF.BR.si, TestBR.si)
CM.BR.si = confusionMatrix(table(TestBR.si[,"FCbr"],pred.BR.si))
CM.BR.si$overall[c(1,2)]

# print(varImp(RF.BR.si))
# plot(varImp(RF.BR.si))
# se sabe q FCbr pertenece al data.frame EntrenamientoBR.si y es la predictora 
# mtry es la cantidad de variables q se usaran para dividir en cada rama

EntrenamientoBR.no <- datosBR.no[Ind,]
dim(EntrenamientoBR.no)
TestBR.no <- datosBR.no[-Ind,]
dim(TestBR.no)
nrow(EntrenamientoBR.no) + nrow(TestBR.no)

control_train <- trainControl(method = "repeatedcv", number = 5,
                              repeats = 10)

RF.BR.no = train(FCbr ~ ., data = EntrenamientoBR.no, method = "rf", trControl = control_train, tuneLength = 2, tuneGrid = expand.grid(mtry=c(2, 4, 10, 20)))
# print(RF.BR.no)
# RF.BR.no$results
pred.BR.no = predict(RF.BR.no, TestBR.no)
CM.BR.no = confusionMatrix(table(TestBR.no[,"FCbr"], pred.BR.no))

tablaACCKAP.BR = rbind(tablaACCKAP.BR, c(CM.BR.si$overall[c(1,2)], CM.BR.no$overall[c(1,2)]))
}
  colnames(tablaACCKAP.BR)=c("ACC.BR.si", "KAPPA.BR.si", "ACC.BR.no", "KAPPA.BR.no")
  
  apply(tablaACCKAP.BR, 2, mean)
  apply(tablaACCKAP.BR, 2, median)
  test1=wilcox.test(tablaACCKAP.BR$ACC.BR.si, tablaACCKAP.BR$ACC.BR.no, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
  test1$p.value
  test2=wilcox.test(tablaACCKAP.BR$KAPPA.BR.si, tablaACCKAP.BR$KAPPA.BR.no, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
  test2$p.value
  cohen.d(tablaACCKAP.BR$ACC.BR.si, tablaACCKAP.BR$ACC.BR.no)
  cohen.d(tablaACCKAP.BR$KAPPA.BR.si, tablaACCKAP.BR$KAPPA.BR.no)
  
    
# print(varImp(RF.BR.no))
# plot(varImp(RF.BR.no))
# se sabe q FCbr pertenece al data.frame EntrenamientoBR.si y es la predictora 
# mtry es la cantidad de variables q se usaran para dividir en cada rama


# ---- GLM LU ---- 
##Datos
dim(matExprLU)
length(FClu)
summary(FClu)
maxPRED = min(floor(length(FClu)/10), length(SIEMPREdeLU)) 
predSI = SIEMPREdeLU[1:maxPRED] 

## Regresion Logistica, un modelo para cada matriz en matExprLU.list variando los genes NOde (ya 
## que tenemos mas noDE que DE)
## y voy registrando AIC y log.lik
## 
regLU=data.frame()
for (i in 1:50)
{
  cambianNOde = sample(1:length(NUNCAdeLU),maxPRED, replace = FALSE)
  predNO = NUNCAdeLU[cambianNOde]
  datos = data.frame(t(matExprLU.list[[i]]), FClu)
  subSI=datos[,predSI]
  subNO= datos[,predNO]
  
  modelo.DE <- glm(FClu ~ ., data = subSI, family = binomial(logit), maxit=100)
  modelo.noDE <- glm(FClu ~ ., data = subNO, family = binomial(logit), maxit=100)
  de=glance(modelo.DE)
  no.de=glance(modelo.noDE)
  
  regLU = rbind(regLU, c(AIC.DE = de$AIC , AIC.notDE = no.de$AIC , log.lik.DE = de$logLik , log.lik.notDE = no.de$logLik))
  
}
colnames(regLU) = c("AIC.DE", "AIC.notDE", "log.lik.DE", "log.lik.notDE")
promediosLU = apply(regLU, 2, mean)
write.csv(regLU, "detalleRegresionLU.csv")
write.csv(promediosLU, "mediasRegresionLU.csv")


# test para ver si la diferencia entre AICs es significativa
# veo si la distr es normal
# SHAPIRO-WILK: The data is normal if the p-value is above 0.05. 
shapiro.test(regLU$AIC.DE) # p-value = 0.0006825 -- no es NORMAL, alfa menor a 0.05 no es normal
densityplot(regLU$AIC.DE)
shapiro.test(regLU$AIC.notDE) # p-value = p-value = 0.2591  -- NORMAL
densityplot(regLU$AIC.notDE)

# igual uso test no  parametrico 
wilcox.test(regLU$AIC.DE,regLU$AIC.notDE, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
#  p-value =  1.824e-14
# existe una probabilidad menor a 0.05 de que  la dif. entre las medianas se deba al azar
# wilcox usa MEDIANAS
medianasLU=apply(regLU, 2, median)
medianasLU
cohen.d(regBR$AIC.DE,regLU$AIC.notDE)
# Cohen's d
# 
# d estimate: -1.942799 (large)
# 95 percent confidence interval:
#     lower     upper 
# -2.424302 -1.461295 

## ---- RF LU ----

  Ind <- createDataPartition(y = FClu, p = 0.7,list = FALSE)
  FCluTrain=FClu[Ind]
  FCluTest=FClu[-Ind]
  tablaACCKAP.LU = data.frame()
  
for (i in 1:50)
{
matExprLU = matExprLU.list[[i]]

datosLU.si = data.frame(t(matExprLU[SIEMPREdeLU,]), FClu)
datosLU.no = data.frame(t(matExprLU[NUNCAdeLU[1:length(SIEMPREdeLU)],]), FClu)

EntrenamientoLU.si <- datosLU.si[Ind,]
dim(EntrenamientoLU.si)
TestLU.si <- datosLU.si[-Ind,]
dim(TestLU.si)
nrow(EntrenamientoLU.si) + nrow(TestLU.si)

control_train <- trainControl(method = "repeatedcv", number = 5,
                              repeats = 10)


RF.LU.si = train(FClu ~ ., data = EntrenamientoLU.si, method = "rf", trControl = control_train, tuneLength = 2, tuneGrid = expand.grid(mtry=c(1, 2, 4, 8)))
# print(RF.LU.si)
# RF.LU.si$results
pred.LU.si = predict(RF.LU.si, TestLU.si)
CM.LU.si = confusionMatrix(table(TestLU.si[,"FClu"],pred.LU.si))
CM.LU.si$overall[c(1,2)]

# print(varImp(RF.LU.si))
# plot(varImp(RF.LU.si))
# se sabe q FClu pertenece al data.frame EntrenamientoLU.si y es la predictora 
# mtry es la cantidad de variables q se usaran para dividir en cada rama

EntrenamientoLU.no <- datosLU.no[Ind,]
dim(EntrenamientoLU.no)
TestLU.no <- datosLU.no[-Ind,]
dim(TestLU.no)
nrow(EntrenamientoLU.no) + nrow(TestLU.no)

control_train <- trainControl(method = "repeatedcv", number = 5,
                              repeats = 10)

RF.LU.no = train(FClu ~ ., data = EntrenamientoLU.no, method = "rf", trControl = control_train, tuneLength = 2, tuneGrid = expand.grid(mtry=c(1, 2, 4, 8)))
# print(RF.LU.no)
# RF.LU.no$results
pred.LU.no = predict(RF.LU.no, TestLU.no)
CM.LU.no = confusionMatrix(table(TestLU.no[,"FClu"], pred.LU.no))

tablaACCKAP.LU =  rbind(tablaACCKAP.LU, c(CM.LU.si$overall[c(1,2)], CM.LU.no$overall[c(1,2)]))
}
  colnames(tablaACCKAP.LU)=c("ACC.LU.si", "KAPPA.LU.si", "ACC.LU.no", "KAPPA.LU.no")
  
  apply(tablaACCKAP.LU, 2, mean)
  apply(tablaACCKAP.LU, 2, median)
  test1=wilcox.test(tablaACCKAP.LU$ACC.LU.si, tablaACCKAP.LU$ACC.LU.no, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
  test1$p.value
  test2=wilcox.test(tablaACCKAP.LU$KAPPA.LU.si, tablaACCKAP.LU$KAPPA.LU.no, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
  test2$p.value
  cohen.d(tablaACCKAP.LU$ACC.LU.si, tablaACCKAP.LU$ACC.LU.no)
  cohen.d(tablaACCKAP.LU$KAPPA.LU.si, tablaACCKAP.LU$KAPPA.LU.no)
  
# print(varImp(RF.LU.no))
# plot(varImp(RF.LU.no))
# se sabe q FClu pertenece al data.frame EntrenamientoLU.no y es la predictora 
# mtry es la cantidad de variables q se usaran para dividir en cada rama

# ---- GLM PR ---- 
##Datos
dim(matExprPR)
length(FCpr)
summary(FCpr)
maxPRED = min(floor(length(FCpr)/10), length(SIEMPREdePR)) 
predSI = SIEMPREdePR[1:maxPRED] 

## Regresion Logistica, un modelo para cada matriz en matExprPR.list variando los genes NOde (ya 
## que tenemos mas noDE que DE)
## y voy registrando AIC y log.lik
## 
regPR=data.frame()
for (i in 1:50)
{
  cambianNOde = sample(1:length(NUNCAdePR),maxPRED, replace = FALSE)
  predNO = NUNCAdePR[cambianNOde]
  datos = data.frame(t(matExprPR.list[[i]]), FCpr)
  subSI=datos[,predSI]
  subNO= datos[,predNO]
  
  modelo.DE <- glm(FCpr ~ ., data = subSI, family = binomial(logit), maxit=100)
  modelo.noDE <- glm(FCpr ~ ., data = subNO, family = binomial(logit), maxit=100)
  de=glance(modelo.DE)
  no.de=glance(modelo.noDE)
  
  regPR = rbind(regPR, c(AIC.DE = de$AIC , AIC.notDE = no.de$AIC , log.lik.DE = de$logLik , log.lik.notDE = no.de$logLik))
  
}
colnames(regPR) = c("AIC.DE", "AIC.notDE", "log.lik.DE", "log.lik.notDE")
promediosPR = apply(regPR, 2, mean)
write.csv(regPR, "detalleRegresionPR.csv")
write.csv(promediosPR, "mediasRegresionPR.csv")


# test para ver si la diferencia entre AICs es significativa
# veo si la distr es normal
# SHAPIRO-WILK: The data is normal if the p-value is above 0.05. 
shapiro.test(regPR$AIC.DE) # p-value = 0.718 -- NORMAL, alfa menor a 0.05 no es normal
densityplot(regPR$AIC.DE)
shapiro.test(regPR$AIC.notDE) # p-value = p-value = p-value = 0.718  -- NORMAL
densityplot(regPR$AIC.notDE)

# igual uso test no  parametrico porq algunos no tienen distr normal
wilcox.test(regPR$AIC.DE,regPR$AIC.notDE, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
# p-value < 2.2e-16
# existe una probabilidad menor a 0.05 de que  la dif. entre las medianas se deba al azar
# wilcox usa MEDIANAS
medianasPR=apply(regPR, 2, median)
medianasPR
cohen.d(regPR$AIC.DE,regPR$AIC.notDE)
# Cohen's d
# 
# d estimate: -3.168979 (large)
# 95 percent confidence interval:
#     lower     upper 
# -3.765020 -2.572937 

## ---- RF PR ----

  Ind <- createDataPartition(y = FCpr, p = 0.7,list = FALSE)
  FCprTrain=FCpr[Ind]
  FCprTest=FCpr[-Ind]
  tablaACCKAP.PR = data.frame()
  
for (i in 1:50)
{
matExprPR = matExprPR.list[[i]]

datosPR.si = data.frame(t(matExprPR[SIEMPREdePR,]), FCpr)
datosPR.no = data.frame(t(matExprPR[NUNCAdePR[1:length(SIEMPREdePR)],]), FCpr)

EntrenamientoPR.si <- datosPR.si[Ind,]
# dim(EntrenamientoPR.si)
TestPR.si <- datosPR.si[-Ind,]
# dim(TestPR.si)
# nrow(EntrenamientoPR.si) + nrow(TestPR.si)

control_train <- trainControl(method = "repeatedcv", number = 5,
                              repeats = 10)

RF.PR.si = train(FCpr ~ ., data = EntrenamientoPR.si, method = "rf", trControl = control_train, tuneLength = 2, tuneGrid = expand.grid(mtry=c(1, 2, 4, 7)))
# print(RF.PR.si)
# RF.PR.si$results
pred.PR.si = predict(RF.PR.si, TestPR.si)
CM.PR.si = confusionMatrix(table(TestPR.si[,"FCpr"],pred.PR.si))

# print(varImp(RF.PR.si))
# plot(varImp(RF.PR.si))
# se sabe q FCpr pertenece al data.frame EntrenamientoPR.si y es la predictora 
# mtry es la cantidad de variables q se usaran para dividir en cada rama

EntrenamientoPR.no <- datosPR.no[Ind,]
dim(EntrenamientoPR.no)
TestPR.no <- datosPR.no[-Ind,]
dim(TestPR.no)
nrow(EntrenamientoPR.no) + nrow(TestPR.no)

control_train <- trainControl(method = "repeatedcv", number = 5,
                              repeats = 10)

RF.PR.no = train(FCpr ~ ., data = EntrenamientoPR.no, method = "rf", trControl = control_train, tuneLength = 2, tuneGrid = expand.grid(mtry=c(1, 2, 4, 7)))
# print(RF.PR.no)
# RF.PR.no$results
pred.PR.no = predict(RF.PR.no, TestPR.no)
CM.PR.no = confusionMatrix(table(TestPR.no[,"FCpr"], pred.PR.no))

tablaACCKAP.PR =  rbind(tablaACCKAP.PR, c(CM.PR.si$overall[c(1,2)], CM.PR.no$overall[c(1,2)]))

# print(varImp(RF.PR.no))
# plot(varImp(RF.PR.no))
# se sabe q FClu pertenece al data.frame EntrenamientoPR.no y es la predictora 
# mtry es la cantidad de variables q se usaran para dividir en cada rama
}
colnames(tablaACCKAP.PR)=c("ACC.PR.si", "KAPPA.PR.si", "ACC.PR.no", "KAPPA.PR.no")

apply(tablaACCKAP.PR, 2, mean)
apply(tablaACCKAP.PR, 2, median)
test1=wilcox.test(tablaACCKAP.PR$ACC.PR.si, tablaACCKAP.PR$ACC.PR.no, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
test1$p.value
test2=wilcox.test(tablaACCKAP.PR$KAPPA.PR.si, tablaACCKAP.PR$KAPPA.PR.no, alternative = "two.sided",paired=FALSE, conf.int = 0.95)
test2$p.value
cohen.d(tablaACCKAP.PR$ACC.PR.si, tablaACCKAP.PR$ACC.PR.no)
cohen.d(tablaACCKAP.PR$KAPPA.PR.si, tablaACCKAP.PR$KAPPA.PR.no)


##### FIN
 