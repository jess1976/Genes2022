
DE.TCGA <- read.csv("~/Dropbox/Research/-- CODIGOS/GENES journal/DE TCGA.csv")
View(DE.TCGA)
GEO1b <- read.csv("~/Dropbox/Research/-- CODIGOS/GENES journal/GEO1b.csv")
# View(GEO1b)
GEO2b <- read.csv("~/Dropbox/Research/-- CODIGOS/GENES journal/GEO2b.csv")
# View(GEO2b)

## OJO: Hay probes q generaron filas con igual GENSYMBOL o ENSID
## POR ESO USO COLUMNA X que es distinta en cada fila
length(which(duplicated(GEO1b$GENSYMBOL)))
length(which(duplicated(GEO1b$ENSID)))

A = DE.TCGA$X # 8267
B = GEO1b$X # 11522
C = GEO2b$X # 2823
length(intersect(A,B)) # 1488
length(intersect(A,C)) # 421
length(intersect(B,C)) # 2691

Comunes = intersect(A,B) 
Comunes = intersect(Comunes, C) 
length(Comunes) # 421
nombres421 = filter(GENES, ENSID %in% Comunes)
write.csv(nombres421, paste(DIR, "Comunes.csv"))

library(Vennerable)
conjs = list(TCGA_THCA=A, GEO33630=B, GEO29265=C)
w <- compute.Venn(Venn(Sets=conjs), doWeights=F, type = "circles") 
library(grid)
grid.newpage()
gp <- VennThemes(w) 
# TITULOS CONJS
gp[["SetText"]][["Set1"]]$fontsize = 12
gp[["SetText"]][["Set2"]]$fontsize = 12
gp[["SetText"]][["Set3"]]$fontsize = 12
# tamaño numeros dentro de los sets
gp[["FaceText"]]$`001`$fontsize = 12
gp[["FaceText"]]$`010`$fontsize = 12
gp[["FaceText"]]$`100`$fontsize = 12
gp[["FaceText"]]$`011`$fontsize = 12
gp[["FaceText"]]$`110`$fontsize = 12
gp[["FaceText"]]$`101`$fontsize = 12
gp[["FaceText"]]$`111`$fontsize = 12

plot(w, gp = gp)

## SOBRE EXPRESADOS EN PTC
Asobre = subset(DE.TCGA, UpDown=="UPptc", select="X") # 2330
Bsobre = subset(GEO1b, UpDown=="UpPTC", select="X") # 6565
Csobre = subset(GEO2b, UpDown=="UpPTC", select="X") # 1501
ComunesSobre = intersect(Asobre$X, Bsobre$X)
ComunesSobre = intersect(ComunesSobre, Csobre$X)
length(ComunesSobre) # 55
nombres55 = filter(GENES, ENSID %in% ComunesSobre)
write.csv(nombres55, paste(DIR, "ComunesSobre.csv"))

conjs = list(TCGA_THCAup=Asobre[,1], GEO33630up=Bsobre[,1], GEO29265up=Csobre[,1])
w <- compute.Venn(Venn(Sets=conjs), doWeights=F, type = "circles") 
grid.newpage()
gp <- VennThemes(w) 
# TITULOS CONJS
gp[["SetText"]][["Set1"]]$fontsize = 12
gp[["SetText"]][["Set2"]]$fontsize = 12
gp[["SetText"]][["Set3"]]$fontsize = 12
# tamaño numeros dentro de los sets
gp[["FaceText"]]$`001`$fontsize = 12
gp[["FaceText"]]$`010`$fontsize = 12
gp[["FaceText"]]$`100`$fontsize = 12
gp[["FaceText"]]$`011`$fontsize = 12
gp[["FaceText"]]$`110`$fontsize = 12
gp[["FaceText"]]$`101`$fontsize = 12
gp[["FaceText"]]$`111`$fontsize = 12

plot(w, gp = gp)


## SUB EXPRESADOS EN PTC
Asub = subset(DE.TCGA, UpDown=="DOWNptc", select="X") # 5937
Bsub = subset(GEO1b, UpDown=="DownPTC", select="X") # 4957
Csub = subset(GEO2b, UpDown=="DownPTC", select="X") # 1322
ComunesSub = intersect(Asub$X, Bsub$X)
ComunesSub = intersect(ComunesSub, Csub$X)
length(ComunesSub) # 40
nombres40 = filter(GENES, ENSID %in% ComunesSub)
write.csv(nombres40, paste(DIR, "ComunesSub.csv"))

conjs = list(TCGA_THCAdown=Asub[,1], GEO33630down=Bsub[,1], GEO29265down=Csub[,1])
w <- compute.Venn(Venn(Sets=conjs), doWeights=F, type = "circles") 
grid.newpage()
gp <- VennThemes(w) 
# TITULOS CONJS
gp[["SetText"]][["Set1"]]$fontsize = 12
gp[["SetText"]][["Set2"]]$fontsize = 12
gp[["SetText"]][["Set3"]]$fontsize = 12
# tamaño numeros dentro de los sets
gp[["FaceText"]]$`001`$fontsize = 12
gp[["FaceText"]]$`010`$fontsize = 12
gp[["FaceText"]]$`100`$fontsize = 12
gp[["FaceText"]]$`011`$fontsize = 12
gp[["FaceText"]]$`110`$fontsize = 12
gp[["FaceText"]]$`101`$fontsize = 12
gp[["FaceText"]]$`111`$fontsize = 12

plot(w, gp = gp)

## DIFERENCIA ENTRE LOS COMUNES GENERALS Y CUANDO SEPARO EN UP Y DOWN
length(which(Comunes %in% ComunesSub)) # 40
length(which(Comunes %in% ComunesSobre)) # 55

iComunes1 = which(Comunes %in% ComunesSub)
iComunes2 = which(Comunes %in% ComunesSobre)
iComunes =c(iComunes1,iComunes2)
noCoinciden = Comunes[-iComunes] # 326

## Verificar con graficos por las dudas
noCoinciden[1]
# lo busque en los 3 datasets, en los dos de GEO aparece down, en TCGA aparece UP
noCoinciden[2]
# En TCGA, DOWN; en GEO esta UP y DOWN (por haber sido distintas sondas)
noCoinciden[3]
# EN TCGA y GEO2, DOWN; en GEO1 UP y DOWN
noCoinciden[4]
# EN TCGA y GEO2, DOWN; en GEO1 UP 

length(A) # 8267
sacarA = A %in% noCoinciden
A = DE.TCGA[!sacarA,]
dim(A) # 7941 4
write.csv(A, paste(DIR, "DE TCGA clean.csv"))

length(B) # 11522
sacarB = B %in% noCoinciden
B = GEO1b[!sacarB,]
length(B) # 11196 4
write.csv(B, paste(DIR, "GEO1b clean.csv"))

length(C) # 2823
sacarC = C %in% noCoinciden
C = GEO2b[!sacarC,]
dim(C) # 2497 4
write.csv(C, paste(DIR, "GEO2b clean.csv"))
