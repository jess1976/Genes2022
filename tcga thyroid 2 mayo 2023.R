## bajando a mano los archivos descargados con XENA BROWSER
## Sacar de general.R

## SAMPLES INFO
TCGA.THCA.GDC_phenotype.tsv.gz <- read.delim("~/Dropbox/Research/-- CODIGOS/GENES journal/TCGA DATA/TCGA-THCA.GDC_phenotype.tsv.gz.tsv")
# TCGA.THCA.GDC_phenotype.tsv.gz <- read.delim("~/Library/CloudStorage/Dropbox/Research/CODIGOS/GENES journal/TCGA DATA/TCGA-THCA.GDC_phenotype.tsv.gz.tsv")
# View(TCGA.THCA.GDC_phenotype.tsv.gz)
DF.samples=TCGA.THCA.GDC_phenotype.tsv.gz
rm(TCGA.THCA.GDC_phenotype.tsv.gz)
summary(as.factor(DF.samples$sample_type.samples))
# Metastatic       Primary Tumor Solid Tissue Normal 
# 8                 507                 100 
## SACAR METASTATIC
DF.samples=subset(DF.samples, sample_type.samples!="Metastatic")
dim(DF.samples)
summary(as.factor(DF.samples$sample_type.samples))

summary(as.factor(DF.samples$primary_diagnosis.diagnoses))
# Carcinoma, NOS           Follicular adenocarcinoma, NOS 
# 1                                        1 
# Follicular carcinoma, minimally invasive     Nonencapsulated sclerosing carcinoma 
# 1                                        4 
# Oxyphilic adenocarcinoma            Papillary adenocarcinoma, NOS 
# 1                                      444 
# Papillary carcinoma, columnar cell  Papillary carcinoma, follicular variant 
# 43                                      120 
### DEJAR PAP ADENO,NOS y PAP CARC,FOLLICULAR V
DF.samples=subset(DF.samples, primary_diagnosis.diagnoses == "Papillary adenocarcinoma, NOS" | primary_diagnosis.diagnoses == "Papillary carcinoma, follicular variant")
dim(DF.samples)
summary(as.factor(DF.samples$sample_type.samples))
muestrasFINAL = DF.samples$submitter_id.samples

## GENES INFO
gencode.v22.annotation.gene <- read.delim("~/Desktop/DATOS/TCGA.BRCA/gencode.v22.annotation.gene.probeMap")
GENES=data.frame(gencode.v22.annotation.gene$id, gencode.v22.annotation.gene$gene)
rm(gencode.v22.annotation.gene)
colnames(GENES)=c("ENSID","GENSYMBOL")
rownames(GENES)=GENES[,1]

## MATRIZ DE EXPRESION

TCGA.THCA.htseq_counts.tsv.gz <- read.delim("~/Dropbox/Research/-- CODIGOS/GENES journal/TCGA DATA/TCGA-THCA.htseq_counts.tsv.gz.tsv")
# TCGA.THCA.htseq_counts.tsv.gz <- read.delim("~/Library/CloudStorage/Dropbox/Research/CODIGOS/GENES journal/TCGA DATA/TCGA-THCA.htseq_counts.tsv.gz.tsv")

# View(TCGA.THCA.htseq_counts.tsv.gz)
CONTEOS=TCGA.THCA.htseq_counts.tsv.gz
# rm(TCGA.THCA.htseq_counts.tsv.gz)
dim(CONTEOS)
# 1ra columna es el ENSEMBL ID
# Sacarla y usarla para nombre de fila
rownames(CONTEOS)=CONTEOS[,1]
CONTEOS = CONTEOS[,-1]
CONTEOS = as.matrix(CONTEOS)
# Nombres de muestras en CONTEOS separarn codigos con "." y en fenotype se separan codigos con "-". Pasamos todo al formato mas usado que es con "-" ("TCGA-DJ-A2PN-01A" ).
# 
aux = colnames(CONTEOS)
aux = gsub(".","-", aux, fixed=T)
colnames(CONTEOS)=aux

### VERIFIAR CUALES DE LAS MUESTRAS EN CONTEOS LAS TENGO EN DF.SAMPLES

muestrasEnCOMUN = intersect(colnames(CONTEOS), DF.samples$submitter_id.samples)

# Me quedo en el DF de muestras, solo con aquellas que aparecen en CONTEOS
DF.samples=subset(DF.samples, DF.samples$submitter_id.samples %in% muestrasEnCOMUN)
dim(DF.samples)
summary(as.factor(DF.samples$sample_type.samples))
## Reordeno y filtro las columnas de CONTEOS para que queden en igual orden que las filas de DF.samples
CONTEOS = CONTEOS[,DF.samples$submitter_id.samples]

iControles=which(DF.samples$sample_type.samples == "Solid Tissue Normal")
iPTCs=which(DF.samples$sample_type.samples == "Primary Tumor")

# ## ARMAR EL RESAMPLEO DE DATOS PTC PARA BALANCEAR LAS CLASES
# ## selecciono 50 subconjuntos de muestras de tumores, con la long de muestras normales
# iPTC.mat=matrix(data=NA, nrow=30, ncol=length(iControles))
# for (i in 1:30)
# {
#   iPTC.mat[i,]=sample(iPTCs, length(iControles), replace = FALSE)
# }
# length(iControles)
# dim(iPTC.mat)
# 
# ## ARMAR EL FACTOR DE CLASES PARA INDICAR SI MUESTRA ES SANA (0) O PTC (1)
# #### EXPRESION DIFERENCIAL - edgeR
# 
# FC = rep(0, ncol(CONTEOS))  
# FC[iPTCs] = 1
# FC = as.factor(FC)
# summary(FC)
# 
# library(edgeR)
# 
# disenio=model.matrix(~0+FC)
# colnames(disenio)= levels(FC)
# head(disenio)
# apply(disenio, 2, sum)
# 
# dge = DGEList(counts=CONTEOS, group=FC)
# dge <- calcNormFactors(dge)
# dge <- estimateGLMCommonDisp(dge, disenio)
# dge <- estimateGLMTagwiseDisp(dge, disenio)
# etest=exactTest(dge)
# #  View(etest[["table"]])
# ttTH <-topTags(etest, n=nrow(CONTEOS), adjust.method = "BH", sort.by="PValue", p.value = 0.01) # el 0.01 de tope de p.value en realidad aplica al FDR
# 
# # si no ajusto pval, toma los q son pval menor a 0.5
# # es una restriccion sobre el FDR (p value ajustado si pongo restriccion en pval)
# #dim(ttTH)
# ## ME QUEDO CON PVAL<=0.01 y FC>2
# tablaDEGs = ttTH$table ## todos FDR<0.01
# dim(tablaDEGs) # Son 8267 genes
# # write.csv(rownames(tablaDEGs),"TCGA sinFC.csv") # puse solo ensembls
# ## si filtro por logFC -- no hay en comun con GEO si agrego este filtro en ambos
# # tablaDEGs = subset(tablaDEGs, abs(logFC)>=2)
# # dim(tablaDEGs) # quedan 188 genes
# 
# ensTCGA=rownames(tablaDEGs)
# library(dplyr)
# ensSymbTCGA = filter(GENES, ENSID %in% ensTCGA)
# 
# # BOXPLOTS ejemplos
# unENSid = ensTCGA[1]
# i = which(rownames(CONTEOS) == unENSid)
# gen = CONTEOS[i,]
# boxplot(gen~FC,  ylab = "Expression Level", xlab="Control (0) Vs. PTC (1)", main = paste('TCGA-THCA: ', GENES[ENSid,]$GENSYMBOL), cex=0.5)
# 
# unENSid = ensTCGA[9]
# i = which(rownames(CONTEOS) == unENSid)
# gen = CONTEOS[i,]
# boxplot(gen~FC,  ylab = "Expression Level", xlab="Control (0) Vs. PTC (1)", main = paste('TCGA-THCA: ', GENES[ENSid,]$GENSYMBOL), cex=0.5)
# ## ylim=c(0,13),
# # FIN BOXPLOTS ejemplos
# 
# ## Cuando FC<0, sub expresados en PTC
# ## Cuando FC>0, sobre expresados en PTC
# ## Separarlos
# 
# iSOBRE = which(tablaDEGs$logFC>0) # sobre-expr en PTC: 2330 genes
# unUpPTC = ensTCGA[iSOBRE[1]]
# i = which(rownames(CONTEOS) == unUpPTC)
# gen = CONTEOS[i,]
# boxplot(gen~FC,  ylab = "Expression Level", xlab="Control (0) Vs. PTC (1)", main = paste('TCGA-THCA: ', GENES[unUpPTC,]$GENSYMBOL), cex=0.5)
# 
# iSUB = which(tablaDEGs$logFC<0) # under-expr en PTC: 5937 genes
# unDownPTC = ensTCGA[iSUB[1]]
# i = which(rownames(CONTEOS) == unDownPTC)
# gen = CONTEOS[i,]
# boxplot(gen~FC,  ylab = "Expression Level", xlab="Control (0) Vs. PTC (1)", main = paste('TCGA-THCA: ', GENES[unDownPTC,]$GENSYMBOL), cex=0.5)
# 
# TCGAfile = cbind(ensSymbTCGA, UpDown = "DOWNptc")
# TCGAfile[iSOBRE,]$UpDown="UPptc"
# DIR = "~/Dropbox/Research/-- CODIGOS/GENES journal/"
# setwd(DIR)
# write.csv(TCGAfile, paste(DIR, "DE TCGA.csv", sep=""))
# 
# 
# # GTEx
# # https://gtexportal.org/home/datasets
# # Gene read counts by tissue
# # gene_reads_2017-06-05_v8_thyroid.gct.gz
# 
# write.csv(CONTEOS, paste(DIR, "CONTEOS.csv"))
# save(CONTEOS, file = "CONTEOS.RData")
# save(FC, file = "FC.RData")
# save(DF.samples, file="DF samples.RData")
# 

## NUEVO - revision MAYO 2023
# 1 ENSG00000122420.8  PTGFR
# 2 ENSG00000130226.15 DPP6
# 3 ENSG00000145864.11   GABRB2
# 4 ENSG00000172667.9  ZMAT3

los4 = c("ENSG00000172667.9","ENSG00000145864.11","ENSG00000130226.15","ENSG00000122420.8")

tumoresLos4 = CONTEOS[los4,iPTCs]
dim(tumoresLos4)
View(tumoresLos4)
# 4 457

Disc.Los4.1 = cut(tumoresLos4[1,], summary(tumoresLos4[1,]), labels=LETTERS[1:5])
Disc.Los4.1[1]="A"
summary(Disc.Los4.1)
Disc.Los4.2 = cut(tumoresLos4[2,], summary(tumoresLos4[2,]), labels=LETTERS[1:5])
Disc.Los4.2[1]="A"
summary(Disc.Los4.2)
Disc.Los4.3 = cut(tumoresLos4[3,], summary(tumoresLos4[3,]), labels=LETTERS[1:5])
Disc.Los4.3[which(is.na(Disc.Los4.3))]="A"
summary(Disc.Los4.3)
Disc.Los4.4 = cut(tumoresLos4[4,], summary(tumoresLos4[4,]), labels=LETTERS[1:5])
Disc.Los4.4[which(is.na(Disc.Los4.4))]="A"
summary(Disc.Los4.4)



# CLINICOS de interes
DF.tumores=DF.samples[iPTCs,]
dim(DF.tumores)
# 457 114
View(DF.tumores)

## SEXO
summary(factor(DF.tumores$gender.demographic))
# female   male 
# 375    137 
Sex=factor(DF.tumores$gender.demographic)
which(is.na(Sex)) # no hay NA

summary(aov(tumoresLos4[1,] ~  Sex))
summary(aov(tumoresLos4[2,] ~  Sex))
summary(aov(tumoresLos4[3,] ~  Sex))
summary(aov(tumoresLos4[4,] ~  Sex))


## EDAD al momento del diagnóstico
summary(DF.tumores$age_at_initial_pathologic_diagnosis)
ageDIAG=DF.tumores$age_at_initial_pathologic_diagnosis
cor(tumoresLos4[1,], ageDIAG, method="kendall")
cor.test(tumoresLos4[1,], ageDIAG, method = "kendall")
cor(tumoresLos4[2,], ageDIAG, method="kendall")
cor.test(tumoresLos4[2,], ageDIAG, method = "kendall")
cor(tumoresLos4[3,], ageDIAG, method="kendall")
cor.test(tumoresLos4[3,], ageDIAG, method = "kendall")
cor(tumoresLos4[4,], ageDIAG, method="kendall")
cor.test(tumoresLos4[4,], ageDIAG, method = "kendall")

cut(ageDIAG, breaks=summary(ageDIAG), labels=c("A","B","C","D","E"))
  # intervalos (breaks) - los cuartiles
  # labels - A (15,34] B (34,46] C (46,46.7] D (46.7,57] E (57,88]
discret.ageDIAG = cut(ageDIAG, breaks=summary(ageDIAG), labels=c("A","B","C","D","E"))
which(is.na(discret.ageDIAG)) # 285 344


## Localizacion primaria
primLoc = factor(DF.tumores$primary_thyroid_gland_neoplasm_location_anatomic_site)
summary(factor(primLoc))
#       Bilateral    Isthmus  Left lobe Right lobe 
# 6         74         19        159        199 

summary(aov(tumoresLos4[1,] ~ primLoc))
boxplot((tumoresLos4[1,] ~ primLoc))
summary(aov(tumoresLos4[2,] ~ primLoc))
boxplot((tumoresLos4[2,] ~ primLoc))
summary(aov(tumoresLos4[3,] ~ primLoc)) ## SI!
boxplot((tumoresLos4[3,] ~ primLoc))
summary(aov(tumoresLos4[4,] ~ primLoc))
boxplot((tumoresLos4[4,] ~ primLoc))

library(PMCMRplus)
PMCMRplus::kwAllPairsNemenyiTest(tumoresLos4[3,] ~ primLoc, method="Chisq")

## RAZA
raza=factor(DF.tumores$race.demographic)
PMCMRplus::kwAllPairsNemenyiTest(tumoresLos4[3,] ~ raza, method="Chisq")
chisq.test(table(Disc.Los4.1, raza))
chisq.test(table(Disc.Los4.2, raza)) # SI!!
chisq.test(table(Disc.Los4.3, raza))
chisq.test(table(Disc.Los4.4, raza)) # SI
summary(table(Disc.Los4.2, raza)) ## equivalente a chisq.test

summary(aov(tumoresLos4[1,] ~  raza))
summary(aov(tumoresLos4[2,] ~  raza)) # SI!
summary(aov(tumoresLos4[3,] ~  raza))
summary(aov(tumoresLos4[4,] ~  raza))

## ESTADIO 
## DF.tumores$tumor_stage.diagnoses

estadio = factor(DF.tumores$tumor_stage.diagnoses)

# junto stage iv, iva, ivc
a = which(estadio=="stage iva")
b = which(estadio=="stage ivc")
c = which(estadio=="stage iv")
estadio[c(a,b,c)]="stage iv"
estadio = factor(estadio)
summary(estadio)

summary(aov(tumoresLos4[1,] ~  estadio))
summary(aov(tumoresLos4[2,] ~  estadio)) ## SI!
summary(aov(tumoresLos4[3,] ~  estadio)) ## SI!
summary(aov(tumoresLos4[4,] ~  estadio)) ## SI!



## GRade low/high
summary(factor(DF.tumores$tumor_grade.diagnoses))


## modelo general
for (i in 1:4)
{
data.lm = lm(tumoresLos4[i,] ~  estadio+Sex+primLoc+raza+ageDIAG)
data.av = aov(data.lm)
print(summary(data.av))
summary(data.lm)
print(summary(data.lm)$r.squared)
}

