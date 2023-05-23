# 
# .ComunesSub <- read.csv("~/Dropbox/Research/-- CODIGOS/GENES journal/ ComunesSub.csv")
# View(.ComunesSub)
# .ComunesSobre <- read.csv("~/Dropbox/Research/-- CODIGOS/GENES journal/ ComunesSobre.csv")
# View(.ComunesSobre)

### Logistic Regression and Classification ----
# BiocManager::install("ggpubr")
# BiocManager::install("MuMin")
library("ggpubr")
library("effsize") # para fc. cohen.d (tamaño del efecto)
library(caret)
library(broom)
# library(MuMin)
library(randomForest)
library(org.Hs.eg.db)
library(STRINGdb)

insDIR = "~/Dropbox/Research/-- CODIGOS/GENES journal/"
setwd(insDIR)
# 
# If you want to load such an .Rdata file into your environment, 
# simply do load(file = "data.Rdata")
# Then, the object is available in your workspace with its old name. 
load(file="CONTEOS.RData") # CONTEOS
load(file="FC.RData") # FC
load("CSUB.RData") # ComunesSub
load("CSOBRE.RData") # ComunesSobre
load("DF samples.RData") # DF.samples

summary(FC)
dim(CONTEOS)
length(ComunesSub) # 40
length(ComunesSobre) # 55
los95 = c(ComunesSub, ComunesSobre)
## Filtrar CONTEOS para dejar solo GENES en ComunesSobre y ComunesSub
CONTEOS = CONTEOS[los95,]
dim(CONTEOS)

# RECURSIVE FEATURE ELIMINATION
# SUPUESTAMENTE 30 runs y elijo los q se repiten
# define the control using a random forest selection function
# control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# # run the RFE algorithm
# results <- rfe(t(CONTEOS), FC, sizes=c(1:10), rfeControl=control)
# Preds = filter(GENES, ENSID %in% predictors(results))
# save(Preds, file = "Preds.RData")
load("Preds.RData")
colS = keytypes(org.Hs.eg.db)[c(3,6,24,26)]
Preds2 = select(org.Hs.eg.db, keys=Preds$GENSYMBOL, columns=colS, keytype="SYMBOL")
# save(Preds2, file = "Preds2.RData")

Preds[1,]$ENSID %in% ComunesSobre # F
Preds[2,]$ENSID %in% ComunesSobre # T
Preds[3,]$ENSID %in% ComunesSobre # F RAROOOOOOOOOOOO
Preds[4,]$ENSID %in% ComunesSobre # F

boxplot(CONTEOS[Preds[1,]$ENSID,]~FC, main=Preds[1,]$GENSYMBOL, ylab="Gene Expression", xlab = "Control (0) vs. PTC (1)", col = c("green", "red"))
boxplot(CONTEOS[Preds[2,]$ENSID,]~FC, main=Preds[2,]$GENSYMBOL,ylab="Gene Expression", xlab = "Control (0) vs. PTC (1)", col = c("green", "red"))
boxplot(CONTEOS[Preds[3,]$ENSID,]~FC, main=Preds[3,]$GENSYMBOL,ylab="Gene Expression", xlab = "Control (0) vs. PTC (1)", col = c("green", "red"))
boxplot(CONTEOS[Preds[4,]$ENSID,]~FC,  main=Preds[4,]$GENSYMBOL,ylab="Gene Expression", xlab = "Control (0) vs. PTC (1)", col = c("green", "red"))


preds3 = merge(Preds, Preds2, by.x="GENSYMBOL", by.y="SYMBOL")
save(preds3, file = "preds3.RData")
los4 = Preds$GENSYMBOL[1:4]


string_db = STRINGdb$new(version="11.5", species=9606, score_threshold = 200)# , network_type="functional")
stIDS = string_db$mp(Preds$GENSYMBOL)
string_db$plot_network(stIDS)
tabla = string_db$get_enrichment(stIDS)

# chisq.test(sub1hmoxDF.t$`HIGH/LOW`, sub1hmoxDF.t$patTUMOR)
# # ---- Comparación expresion VS CLINICAL STAGE
# hmoxDF.t$FclSTAGE = factor(hmoxDF.t$FclSTAGE)
# kruskal.test(hmoxDF.t$expresion, hmoxDF.t$FclSTAGE)
# # pval: 0.43
# # No hay diferencia significativa entre las medias de los grupos
# boxplot(hmoxDF.t$expresion ~ factor(hmoxDF.t$FclSTAGE))




listaRFE=predictors(results)
for (i in 1:10)
{
  results <- rfe(t(CONTEOS), FC, sizes=c(1:10), rfeControl=control)
  listaRFE = intersec(listaRFE, predictors(results))
  }

# summarize the results
print(results)
# list the chosen features
SeleccionRFE = predictors(results)

save(SeleccionRFE, file="SeleccionRFE.RData")
# plot the results
plot(results, type=c("g", "o"))

load("GENES.Rdata")
library(dplyr)
Preds = filter(GENES, ENSID %in% SeleccionRFE)

## ENRIQUECIMIENTO DE ESOS GENES
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(clusterProfiler)
library(DOSE)

keytypes(org.Hs.eg.db)
kk = enrichGO(gene = Preds$GENSYMBOL, OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL",ont = "BP", pvalueCutoff=0.05)

View(head(kk))

DF.GO = as.data.frame(kk)
nrow(DF.GO) ## Cant de terminos enriquecidos con cutoff 0.05 q se usa para limitar el P.ADJ

## espera ENTREZ ID
kk = enrichKEGG(gene=Preds$GENSYMBOL, pvalueCutoff=0.05)



# modelo.GLM <- glm(FC ~ ., data = data.frame(t(CONTEOS[1:1000,])), family = binomial(logit), maxit=100)
# glance(modelo.GLM)
# tidy(modelo.GLM)
# summary(modelo.GLM)

# 
# Max.Vars = floor(length(FC)/10)
# dredge(modelo.GLM, m.lim=c(0,Max.Vars)) # MuMin elige el mejor modelo probando todas las combinaciones de predictoras

Comunes = c(ComunesSub, ComunesSobre)
los95 = filter(GENES, ENSID %in% Comunes)
kk = enrichGO(gene=.Comunes$GENSYMBOL, OrgDb = 'org.Hs.eg.db', keyType = "SYMBOL",ont = "BP", pvalueCutoff=0.01)



ensids <- c("ENSG00000130720", "ENSG00000103257", "ENSG00000156414", 
            "ENSG00000144644", "ENSG00000159307", "ENSG00000144485")
keytypes(org.Hs.eg.db)
colS = keytypes(org.Hs.eg.db)[c(3,6,24,26)]
select(org.Hs.eg.db, keys=ensids, columns=colS, keytype="ENSEMBL")

todos = select(org.Hs.eg.db, keys=.Comunes$GENSYMBOL, columns=colS, keytype="SYMBOL")
entrezTODOS = unique(todos$ENTREZID)
kk = enrichGO(gene=entrezTODOS, OrgDb = 'org.Hs.eg.db', keyType = "ENTREZID",ont = "BP", pvalueCutoff=1)

write.csv(entrezTODOS, "entrezTODOS.csv")
Comunes
kk = enrichKEGG(gene=entrezTODOS, pvalueCutoff=0.05)

#### REVISION MAJOR CHANGES ####
#### 
#### 
library(edgeR)
# USAR CONTEOS RESUMIDA A los95
# calcular los FC
disenio=model.matrix(~0+FC)
colnames(disenio)= levels(FC)
head(disenio)
apply(disenio, 2, sum)

dge = DGEList(counts=CONTEOS, group=FC)
dge <- calcNormFactors(dge)
dge <- estimateGLMCommonDisp(dge, disenio)
dge <- estimateGLMTagwiseDisp(dge, disenio)
etest=exactTest(dge)
#  View(etest[["table"]])
ttTH <-topTags(etest, n=nrow(CONTEOS), adjust.method = "BH", sort.by="PValue", p.value = 0.01) # el 0.01 de tope de p.value en realidad aplica al FDR

# si no ajusto pval, toma los q son pval menor a 0.5
# es una restriccion sobre el FDR (p value ajustado si pongo restriccion en pval)
#dim(ttTH)
## ME QUEDO CON PVAL<=0.01 y FC>2
tablaDEGs = ttTH$table ## todos FDR<0.01
View(tablaDEGs)
tablaDEGs = cbind(tablaDEGs, ENS = rownames(tablaDEGs))

ensTCGA=rownames(tablaDEGs)
library(dplyr)
ensSymbTCGA = filter(GENES, ENSID %in% ensTCGA)

# GENES de "tcga thyroid 2.R"
View(ensSymbTCGA)
genes.of.interest = c("DPP6", "GABRB2", "PTGFR", "ZMAT3")
mapeo = subset(ensSymbTCGA, GENSYMBOL %in% genes.of.interest)
tabla.paper = merge(mapeo, tablaDEGs,by.x="ENSID", by.y="ENS")
write.csv(tabla.paper, "table_FoldCHANGE.csv")

