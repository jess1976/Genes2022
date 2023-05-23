.DE.TCGA.clean <- read.csv("~/Dropbox/Research/-- CODIGOS/GENES journal/ DE TCGA clean.csv")
View(.DE.TCGA.clean)

.GEO1b.clean <- read.csv("~/Dropbox/Research/-- CODIGOS/GENES journal/ GEO1b clean.csv")
View(.GEO1b.clean)

.GEO2b.clean <- read.csv("~/Dropbox/Research/-- CODIGOS/GENES journal/ GEO2b clean.csv")
View(.GEO2b.clean)

A = .DE.TCGA.clean$X # 7941
B = .GEO1b.clean$X # 11196
C = .GEO2b.clean$X # 2497

length(intersect(A,B)) # 7941
length(intersect(A,C)) # 2497
length(intersect(B,C)) # 2497

Comunes = intersect(A,B) 
Comunes = intersect(Comunes, C) 
length(Comunes) # 2497
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

plot(w, gp = gp) # Fig 1

## SOBRE EXPRESADOS EN PTC
Asobre = subset(.DE.TCGA.clean, UpDown=="UPptc", select="X") # 2330
Bsobre = subset(.GEO1b.clean, UpDown=="UpPTC", select="X") # 6565
Csobre = subset(.GEO2b.clean, UpDown=="UpPTC", select="X") # 1501
ComunesSobre = intersect(Asobre$X, Bsobre$X)
ComunesSobre = intersect(ComunesSobre, Csobre$X)
length(ComunesSobre) # 55
nombres55 = filter(GENES, ENSID %in% ComunesSobre)
write.csv(nombres55, paste(DIR, "ComunesSobre.csv"))
save(ComunesSobre, file= "CSOBRE.RData")

conjs = list(TCGA_THCA=Asobre[,1], GEO33630=Bsobre[,1], GEO29265=Csobre[,1])
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

plot(w, gp = gp) # fig 2a

## SUB EXPRESADOS EN PTC
Asub = subset(.DE.TCGA.clean, UpDown=="DOWNptc", select="X") # 5937
Bsub = subset(.GEO1b.clean, UpDown=="DownPTC", select="X") # 4957
Csub = subset(.GEO2b.clean, UpDown=="DownPTC", select="X") # 1322
ComunesSub = intersect(Asub$X, Bsub$X)
ComunesSub = intersect(ComunesSub, Csub$X)
length(ComunesSub) # 40
nombres40 = filter(GENES, ENSID %in% ComunesSub)
write.csv(nombres40, paste(DIR, "ComunesSub.csv"))
save(ComunesSub, file = "CSUB.RData")

conjs = list(TCGA_THCA=Asub[,1], GEO33630=Bsub[,1], GEO29265=Csub[,1])
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

plot(w, gp = gp) # fig 2b

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
A = A[!sacarA]
length(A) # 7941
write.csv(A, paste(DIR, "DE TCGA clean.csv"))

length(B) # 11522
sacarB = B %in% noCoinciden
B = B[!sacarB]
length(B) # 11196
write.csv(B, paste(DIR, "GEO1b clean.csv"))

length(C) # 2823
sacarC = C %in% noCoinciden
C = C[!sacarC]
length(C) # 2497
write.csv(C, paste(DIR, "GEO2b clean.csv"))

