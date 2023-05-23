# https://gtexportal.org/home/datasets
# Gene read counts by tissue
# gene_reads_2017-06-05_v8_thyroid.gct.gz

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)

# filtroNOMBRE filtra los nombres de las filas que tienen formato:
# "A1BG|1"
# Está separado con | el SYMBOL del ENTREZ ID
# filtroNOMBRE devuelve el gene SYMBOL
 
filtroNOMBRE = function(nombreLargo){
  Splitted = strsplit(nombreLargo,"|")[[1]]
  iBARRA = which(Splitted=="|")
  Splitted = Splitted[1:(iBARRA-1)]
  nombreGEN=c()
  for (i in 1:length(Splitted))
    nombreGEN=paste(nombreGEN, Splitted[i], sep="")
  return(nombreGEN)
}



# datQuery = GDCquery(project="TCGA-THCA", 
#                     data.category = "Gene expression",
#                     data.type = "Gene expression quantification",
#                     platform = "Illumina HiSeq", 
#                     file.type  = "results",
#                     experimental.strategy = "RNA-Seq",
#                     legacy = TRUE)

datQuery = GDCquery(project="TCGA-THCA", 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type  = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)

GDCdownload(datQuery, method = "api", files.per.chunk = 30) 
data <- GDCprepare(datQuery, summarizedExperiment = TRUE)
class(data)

##### INFO MUESTRAS #### 
# View(colData(data))
View(data.frame(colData(data)))
DF_samples=data.frame(colData(data))
# FULL NAME
DF_samples$barcode[1]
# SAMPLE NAME
DF_samples$sample[1]
# PATIENT
DF_samples$sample[1]
# TYPE OF SAMPLE
DF_samples$shortLetterCode[1]
DF_samples$definition[1]
DF_samples$sample_type[1]
summary(factor(DF_samples$sample_type))
# Metastatic       Primary Tumor Solid Tissue Normal 
# 8                 505                  59 
summary(factor(DF_samples$primary_diagnosis))
#### DEJAR: primary diagnosis

# Carcinoma, NOS           Follicular adenocarcinoma, NOS 
# 1                                        1 
# Follicular carcinoma, minimally invasive     Nonencapsulated sclerosing carcinoma 
# 1                                        4 
# Oxyphilic adenocarcinoma            Papillary adenocarcinoma, NOS 
# 1                                      412 
# Papillary carcinoma, columnar cell  Papillary carcinoma, follicular variant 
# 40                                      112 


iPA = which(DF_samples$primary_diagnosis=="Papillary adenocarcinoma, NOS")
iPC = which(DF_samples$primary_diagnosis=="Papillary carcinoma, follicular variant")
namePA = DF_samples$barcode[iPA]
namePC = DF_samples$barcode[iPC]
### OJO! Estos no son todo MUESTRAS DE TUMOR
### Hay muestras normal tissue q, como son sacadas de paciente cancer, en primary diagnosis tienen alguno de estos dos.... por eso era mejor usar GTEx...

PTC = c(namePA, namePC)
matrizF = matriz[, PTC]


##### INFO GENES #### 
DF_genes=data.frame(rowData(data))
DF_genes$gene_id[1]
DF_genes$entrezgene[1]
DF_genes$ensembl_gene_id[1]

##### CONTEOS ##### 
matriz = assay(data)
dim(matriz)   #  [1] 19947 genes   572 muestras
rownames(matriz)
colnames(matriz)
matrizF = matriz[,PTC]
dim(matrizF)

# which(grepl("HMOX1",rownames(matriz)))
# which(rownames(matriz)=="HMOX1") # 2152
# dim(matriz) # 19947   572
# strsplit(rownames(matriz)[1], "|")

#########################  Filtering & DE analysis con edgeR
# BiocManager::install("edgeR")
library(edgeR)

dataFilt <- TCGAanalyze_Filtering(tabDF = matrizF,
                                  method = "quantile", 
                                  qnt.cut =  0.25)
dim(dataFilt)
# [1] 14960  524

# which(grepl("HINT1",rownames(dataFilt))) # 1630


#### get subtype SAMPLE information

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))
length(samplesNT)
# 56

# selection of tumor samples "TP"

samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                 typesample = c("TP"))
length(samplesTP)
# 460

##################### Diff.expr.analysis (DEA)  ANALISIS ESTADISTICO UNIVARIADO
##################### 

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                mat2 = dataFilt[,samplesTP], 
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 2,
                            method = "exactTest")
# IMPORTANTE
# por el tema DESBALANCE DE CLASES (entre #NT y #PT)
# DECIDIR SI RE-SAMPLEAMOS LOS TUMORES Y HACEMOS VARIAS CORRIDAS, QUEDANDONOS 
# CON LOS Q APARECIERON DE EN TODAS LAS CORRIDAS

# method is 'glmLRT' (1) or 'exactTest' (2) used for edgeR (1) Fit a negative binomial generalized log-linear model to the read counts for each gene (2) Compute genewise exact tests for differences in the means between two groups of negative-binomially distributed counts.

head(dataDEGs)
orden=order(dataDEGs$PValue)
dataDEGsORDENADO=dataDEGs[orden,]
head(dataDEGsORDENADO)
rownames(dataDEGsORDENADO[1:20,])
dim(dataDEGs) # 2162 genes DE con FDR<0.01, logFC.cut = 1, method = glmRT

dim(dataDEGs) # 775 genes DE con FDR<0.01, logFC.cut = 2, method = glmRT

dim(dataDEGs) # 778 genes DE con FDR<0.01, logFC.cut = 2, method = exactTest

###no hace falta bajando results_normalized
###genesDE2 = unlist(lapply(rownames(dataDEGsORDENADO), filtroNOMBRE))
# todos los genesDE estan en genesDE2

factorNORMALES.TRUE = factor(colnames(dataFilt) %in% samplesNT)
summary(factorNORMALES.TRUE) # calse1: tumor clase2: normal
elMAS = rownames(dataDEGsORDENADO[3,])# logFC 8.992538, sobreexpr en 


boxplot(log(dataFilt[elMAS,]), factorNORMALES.TRUE)
out = (boxplot(dataFilt[elMAS,], factorNORMALES.TRUE)$out)

out <- boxplot.stats(dat$hwy)$out
out_ind <- which(dat$hwy %in% c(out))
out_ind


################################# RANDOM FOREST PARA feature selection MULTI-VARIADO

# BiocManager::install("randomForest")
library("randomForest")

#reshape de la matriz de conteos para pasar al clasificador --- filas genes a columnas.

nuevaM=t(dataFilt)
dim(nuevaM)  # 1215 samples x 14960 features

# armar factor de clases 0 NORMAL 1 TUMOR   -- 113 muestras de cada una

normal=samplesNT
tumor=samplesTP[1:113]

nuevaM=nuevaM[c(normal,tumor),]

fY=as.factor(c(rep(0,113),rep(1,113)))
#habia genes repetidos en mas de una columna
#renombre todas las col agregando un numero de col al comienzo
nombresCOL=colnames(nuevaM)

nombresCOL=paste(nombresCOL,1:14955)

colnames(nuevaM)=nombresCOL # si quedan con un espacio no le gusta al classif
#My guess (strongly supported by the example below) is that randomForest() can't handle non-syntactic variable/column names, i.e. ones with spaces or punctuation other than dots in them. You could try names(flight) <- make.names(names(flight)) to fix this.
#
# con make.names los pone como el clasif quiere
colnames(nuevaM)=make.names(colnames(nuevaM))
nuevaM=cbind(nuevaM,fY)

nuevaMrf = randomForest(x=nuevaM, y=fY, importance=TRUE, proximity=FALSE)
#The loss function is mse for regression and gini-impurity for classification. 
importancia=data.frame(importance(nuevaMrf))

 max(importancia$MeanDecreaseGini)
#[1] 1.882229
 which.max(importancia$MeanDecreaseGini)
#[1] 749

sort(importancia$MeanDecreaseGini, decreasing = TRUE)

rfORD=importancia[order(importancia$MeanDecreaseGini, decreasing = TRUE),]
rfORD[1:20,]

print(nuevaMrf)






################ clinical information
# creo q no es necesario hacerlo aparte porque se descargó dentro de colData()
# 
# dataClin <- GDCquery_clinic(project="TCGA-THCA", type = "clinical")
# dim(dataClin)
# # [1] 507 69 # 507 registros de pacientes
# 
# dataClin[1,]$primary_diagnosis
# # dataClin[1,]$bcr_patient_barcode ## idem anterior
# dataClin[1:20,]$primary_diagnosis
# 
# dataClin[1,]$submitter_id
# dataClin[1,]$classification_of_tumor
# dataClin[1,]$progression_or_recurrence
# dataClin[1,]$age_at_ini # son dias, dividir por 365 para sacar en años
# 
# which(dataClin$submitter_id=="TCGA-BJ-A28W") # 345
# # Asi en la pagina: "TCGA-BJ-A28W-01A" 
# # Pagina: https://xenabrowser.net/datapages/?dataset=TCGA-THCA.GDC_phenotype.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# # Así en dataClin: "TCGA-BJ-A28W"
# View(dataClin[345,])	
# which(dataClin$submitter_id==(colnames(matriz)[1])) # ninguno, en la matriz los nombres estan mas largos
# which(grepl("TCGA-BJ-A28W",colnames(matriz)))
# [1] 219 307 # Dos... 
# colnames(matriz)[c(219,307)]
# [1] "TCGA-BJ-A28W-11A-11R-A32Y-07" "TCGA-BJ-A28W-01A-11R-A32Y-07"
#### Creo q 11A es tejido normal, y 01A es tumor...
#### en el cuarto grupo de letras

