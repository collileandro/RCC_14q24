## libraries ##

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(EDASeq)
library(survival)
library(survminer)

RemOut <- function(x, na.rm = TRUE,...) {
x[!x %in% boxplot.stats(x)$out]}

############## Download TCGA data ####################

library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

GDCdownload(query)
data <- GDCprepare(query)

dataPrep <- TCGAanalyze_Preprocessing(object = data, 
                                      cor.cut = 0.6,
                                      datatype = "HTSeq - Counts")                      

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent") 
                                      
normal <- TCGAquery_SampleTypes(colnames(dataPrep), typesample = c("NT"))
tumor <- TCGAquery_SampleTypes(colnames(dataPrep), typesample = c("TP"))

############################## Sub-selecting Genes Expression ##############################

## DPF3 = ENSG00000205683 ##
DPF3 <- as.data.frame(dataNorm[grep("ENSG00000205683", rownames(dataNorm)),] )
DPF3$name = rownames(DPF3)
DPF3$full.RNA = rownames(DPF3)
DPF3$code = NA
for (i in 1:nrow(DPF3)) { DPF3$code[i] <- substring(DPF3[i,2], 14, 15) }
for (i in 1:nrow(DPF3)) { DPF3[i,2] <- substring(DPF3[i,2], 1, 12) }
#DPF3 <- DPF3[which(DPF3$code=='01'), 1:3]
colnames(DPF3) <- c("DPF3", "name", "full.RNA")
#hist(log(DPF3))    

## IRX2 = ENSG00000170561 ##
IRX2 <- as.data.frame(dataNorm[grep("ENSG00000170561", rownames(dataNorm)),] )
IRX2$name = rownames(IRX2)
IRX2$full.RNA = rownames(IRX2)
IRX2$code = NA
for (i in 1:nrow(IRX2)) { IRX2$code[i] <- substring(IRX2[i,2], 14, 15) }
for (i in 1:nrow(IRX2)) { IRX2[i,2] <- substring(IRX2[i,2], 1, 12) }
#IRX2 <- IRX2[which(IRX2$code=='01'), 1:3]
colnames(IRX2) <- c("IRX2", "name", "full.RNA")
#hist(log(IRX2))    

## IRX5 = ENSG00000176842 ##
IRX5 <- as.data.frame(dataNorm[grep("ENSG00000176842", rownames(dataNorm)),] )
IRX5$name = rownames(IRX5)
IRX5$full.RNA = rownames(IRX5)
IRX5$code = NA
for (i in 1:nrow(IRX5)) { IRX5$code[i] <- substring(IRX5[i,2], 14, 15) }
for (i in 1:nrow(IRX5)) { IRX5[i,2] <- substring(IRX5[i,2], 1, 12) }
#IRX5 <- IRX5[which(IRX5$code=='01'), 1:3]
colnames(IRX5) <- c("IRX5", "name", "full.RNA")
#hist(log(IRX5)) 


data.gene <- merge(DPF3, IRX2, by.x="full.RNA", by.y="full.RNA")
data.gene <- merge(data.gene, IRX5, by.x="full.RNA", by.y="full.RNA")
data.gene <- data.gene[, c("name", "full.RNA", "DPF3", "IRX2", "IRX5")]

data.gene$Sample <- substr(data.gene$full.RNA, 1, 15)

############################## Genotyping Data ##############################

geno <- read.table("SNP_plus_RNA.data", as.is=T, header=T)

geno.rs4903064 <- geno[, c("Individual", "rs4903064.73279420.T.C")]


############################## HYPOXIA SCORE Data ##############################

hypo.KIRC <- read.csv("Hypoxia_score_cell_rep_2016.csv", as.is=T, header=T)

## hypoxia data is from https://www.sciencedirect.com/science/article/pii/S2211124716301279?via%3Dihub#fig2

############################# Comb RNA, SNP and hypo ###############################

data <- merge(data.gene, geno.rs4903064, by.x="name", by.y="Individual")
data <- merge(data, hypo.KIRC, by.x="Sample", by.y="Patient.tumor", all=T)

data <- mutate(data, rs4903064.col = case_when(rs4903064.73279420.T.C == 0 ~ "yellow", 
										 rs4903064.73279420.T.C == 1 ~ "orange", 
										 rs4903064.73279420.T.C == 2 ~ "orange"))

############################# DPF3 vs IRX2/5 #################################

############ normal ############

data.normal <- data[which(data$full.RNA %in% normal),]

data.normal$DPF3 <- log10(data.normal$DPF3)
data.normal$IRX2 <- log10(data.normal$IRX2)
data.normal$IRX5 <- log10(data.normal$IRX5)

par(mfrow=c(2,2))

signif(cor.test(data.normal$DPF3, data.normal$IRX2)$p.value, digits=2) -> p
round(cor.test(data.normal$DPF3, data.normal$IRX2)$estimate,2) -> est
title = paste0("p = ", p, "; Estimate = ", est)
plot(data.normal$IRX2, data.normal$DPF3, pch = 19, col="blue", ylab="log10(DPF3 Normalized Counts)", xlab="log10(IRX2 Normalized Counts)", main="DPF3 vs IRX2 - Normal Renal (TCGA, n=71)")
title(title, line=0.5, cex.main=0.9)
abline(lm(DPF3 ~ IRX2, data.normal), col="blue")

signif(cor.test(data.normal$DPF3, data.normal$IRX5)$p.value, digits=2) -> p
round(cor.test(data.normal$DPF3, data.normal$IRX5)$estimate,2) -> est
title = paste0("p = ", p, "; Estimate = ", est)
plot(data.normal$IRX5, data.normal$DPF3, pch = 19, col="blue", ylab="log10(DPF3 Normalized Counts)", xlab="log10(IRX5 Normalized Counts)", main="DPF3 vs IRX5 - Normal Renal (TCGA, n=71)")
abline(lm(DPF3 ~ IRX5, data.normal), col="blue")
title(title, line=0.5, cex.main=0.9)


summary(lm(DPF3 ~ IRX2 + IRX5 + rs4903064.73279420.T.C, data.normal))

#par(mfrow=c(1,2))

plot(data.normal$IRX2, data.normal$DPF3, pch = 19, col=data.normal$rs4903064.col, , ylab="log10(DPF3 Normalized Counts)", xlab="log10(IRX2 Normalized Counts)", main="DPF3 vs IRX2 - Normal Renal (TCGA, n=71)")
abline(lm(DPF3 ~ IRX2, data.normal[which(data.normal$rs4903064.col %in% c("red", "orange")),]), col="orange")
#abline(lm(DPF3 ~ IRX2, data.normal[which(data.normal$rs4903064.col=="orange"),]), col="orange")
abline(lm(DPF3 ~ IRX2, data.normal[which(data.normal$rs4903064.col=="yellow"),]), col="yellow")
IRX2.coef.CCTC <- lm(DPF3 ~ IRX2, data.normal[which(data.normal$rs4903064.col %in% c("red", "orange")),])$coefficients[2]
IRX2.coef.TT <- lm(DPF3 ~ IRX2, data.normal[which(data.normal$rs4903064.col=="yellow"),])$coefficients[2]
title = paste0("beta.TT=", round(IRX2.coef.TT, 2), "; beta.TC.CC=", round(IRX2.coef.CCTC, 2))
title(title, line=0.5, cex.main=0.9)

plot(data.normal$IRX5, data.normal$DPF3, pch = 19, col=data.normal$rs4903064.col, , ylab="log10(DPF3 Normalized Counts)", xlab="log10(IRX5 Normalized Counts)", main="DPF3 vs IRX5 - Normal Renal (TCGA, n=71)")
abline(lm(DPF3 ~ IRX5, data.normal[which(data.normal$rs4903064.col %in% c("red", "orange")),]), col="orange")
#abline(lm(DPF3 ~ IRX5, data.normal[which(data.normal$rs4903064.col=="orange"),]), col="orange")
abline(lm(DPF3 ~ IRX5, data.normal[which(data.normal$rs4903064.col=="yellow"),]), col="yellow")
IRX5.coef.CCTC <- lm(DPF3 ~ IRX5, data.normal[which(data.normal$rs4903064.col %in% c("red", "orange")),])$coefficients[2]
IRX5.coef.TT <- lm(DPF3 ~ IRX5, data.normal[which(data.normal$rs4903064.col=="yellow"),])$coefficients[2]
title = paste0("beta.TT=", round(IRX5.coef.TT, 2), "; beta.TC.CC=", round(IRX5.coef.CCTC, 2))
title(title, line=0.5, cex.main=0.9)

dev.print(pdf, "Plot_Cor_DPF3_IRX_com_alelos_rs4903064.NORMAL.pdf")

summary(lm(DPF3 ~ IRX2 + IRX5 + rs4903064.73279420.T.C, data.normal))

summary(lm(DPF3 ~ IRX2 * IRX5 * rs4903064.73279420.T.C, data.normal))



#### 3D PLOT #####

library("plot3D")

par(mfrow=c(1,3))

## ALL ##

x <- data.normal$IRX5
y <- data.normal$IRX2
z <- data.normal$DPF3
# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 2, 
    theta = 20, phi = 20, ticktype = "detailed",
    xlab = "IRX5 Expression", ylab = "IRX2 Expression", zlab = "DPF3 Expression",  
    surf = list(x = x.pred, y = y.pred, z = z.pred,  
    facets = NA, fit = fitpoints), main = "Normal KIRC (n=71)")
    
## rs4903064 CT.CC ##

x <- data.normal[which(data.normal$rs4903064.col %in% c("red", "orange")),]$IRX5
y <- data.normal[which(data.normal$rs4903064.col %in% c("red", "orange")),]$IRX2
z <- data.normal[which(data.normal$rs4903064.col %in% c("red", "orange")),]$DPF3
# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 2, 
    theta = 20, phi = 20, ticktype = "detailed",
    xlab = "IRX5 Expression", ylab = "IRX2 Expression", zlab = "DPF3 Expression",  
    surf = list(x = x.pred, y = y.pred, z = z.pred,  
    facets = NA, fit = fitpoints), main = "rs4903064 TC/CC Normal KIRC") 

## rs4903064 TT ##
    
x <- data.normal[which(data.normal$rs4903064.col %in% c("yellow")),]$IRX5
y <- data.normal[which(data.normal$rs4903064.col %in% c("yellow")),]$IRX2
z <- data.normal[which(data.normal$rs4903064.col %in% c("yellow")),]$DPF3
# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 2, 
    theta = 20, phi = 20, ticktype = "detailed",
    xlab = "IRX5 Expression", ylab = "IRX2 Expression", zlab = "DPF3 Expression",  
    surf = list(x = x.pred, y = y.pred, z = z.pred,  
    facets = NA, fit = fitpoints), main = "rs4903064 TT Normal KIRC")           

dev.print(pdf, "Plot_3D_Todos_e_estratificado_por_alelos_rs4903064.pdf")

## FINAL PLOT NORMAL ##

#par(mfrow=c(2,2))

layout(matrix(c(1,3,2,3), 2, 2, byrow = TRUE))

signif(cor.test(data.normal$DPF3, data.normal$IRX2)$p.value, digits=2) -> p
round(cor.test(data.normal$DPF3, data.normal$IRX2)$estimate,2) -> est
title = paste0("p = ", p, "; Estimate = ", est)
plot(data.normal$IRX2, data.normal$DPF3, pch = 19, col="blue", ylab="log10(DPF3)", xlab="log10(IRX2)", main="DPF3 vs IRX2")
title(title, line=0.5, cex.main=0.9)
abline(lm(DPF3 ~ IRX2, data.normal), col="blue")

signif(cor.test(data.normal$DPF3, data.normal$IRX5)$p.value, digits=2) -> p
round(cor.test(data.normal$DPF3, data.normal$IRX5)$estimate,2) -> est
title = paste0("p = ", p, "; Estimate = ", est)
plot(data.normal$IRX5, data.normal$DPF3, pch = 19, col="blue", ylab="log10(DPF3)", xlab="log10(IRX5)", main="DPF3 vs IRX5")
abline(lm(DPF3 ~ IRX5, data.normal), col="blue")
title(title, line=0.5, cex.main=0.9)

x <- data.normal$IRX5
y <- data.normal$IRX2
z <- data.normal$DPF3
# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 2, clab = c("log10(DPF3)"),
    theta = 20, phi = 20, colkey=list(side = 1), # ticktype = "detailed",
    xlab = "log10(IRX5)", ylab = "log10(IRX2", zlab = "log10(DPF3)",  
    surf = list(x = x.pred, y = y.pred, z = z.pred,
    facets = NA, fit = fitpoints), main = "DPF3 vs IRX2 vs IRX5")

    
dev.print(pdf, "Plot_DPF3_IRX2-5.pdf")

############ tumor ############

data.tumor <- data[which(data$full.RNA %in% tumor),]
data.tumor <- data.tumor[complete.cases(data.tumor),]

data.tumor$DPF3 <- log10(data.tumor$DPF3+1)
data.tumor$IRX2 <- log10(data.tumor$IRX2+1)
data.tumor$IRX5 <- log10(data.tumor$IRX5+1)

par(mfrow=c(1,2))
signif(cor.test(data.tumor$DPF3, data.tumor$IRX2)$p.value, digits=2) -> p
round(cor.test(data.tumor$DPF3, data.tumor$IRX2)$estimate,2) -> est
title = paste0("p = ", p, "; Estimate = ", est)
plot(data.tumor$IRX2, data.tumor$DPF3, pch = 19, col="blue", ylab="DPF3 Normalized Counts", xlab="IRX2 Normalized Counts", main="Tumor Renal (TCGA, n=538)")
abline(lm(DPF3 ~ IRX2, data.tumor), col="blue")
title(title, line=-1)

signif(cor.test(data.tumor$DPF3, data.tumor$IRX5)$p.value, digits=2) -> p
round(cor.test(data.tumor$DPF3, data.tumor$IRX5)$estimate,2) -> est
title = paste0("p = ", p, "; Estimate = ", est)
plot(data.tumor$IRX5, data.tumor$DPF3, pch = 19, col="blue", ylab="DPF3 Normalized Counts", xlab="IRX5 Normalized Counts", main="Tumor Renal (TCGA, n=538)")
abline(lm(DPF3 ~ IRX5, data.tumor), col="blue")
title(title, line=-1)

par(mfrow=c(1,2))

plot(data.tumor$IRX2, data.tumor$DPF3, pch = 19, col=data.tumor$rs4903064.col, , ylab="DPF3 Normalized Counts", xlab="IRX2 Normalized Counts", main="Tumor Renal (TCGA, n=538)")
abline(lm(DPF3 ~ IRX2, data.tumor[which(data.tumor$rs4903064.col=="red"),]), col="red")
abline(lm(DPF3 ~ IRX2, data.tumor[which(data.tumor$rs4903064.col=="orange"),]), col="orange")
abline(lm(DPF3 ~ IRX2, data.tumor[which(data.tumor$rs4903064.col=="yellow"),]), col="yellow")
plot(data.tumor$IRX5, data.tumor$DPF3, pch = 19, col=data.tumor$rs4903064.col, , ylab="DPF3 Normalized Counts", xlab="IRX5 Normalized Counts", main="Tumor Renal (TCGA, n=538)")
abline(lm(DPF3 ~ IRX5, data.tumor[which(data.tumor$rs4903064.col=="red"),]), col="red")
abline(lm(DPF3 ~ IRX5, data.tumor[which(data.tumor$rs4903064.col=="orange"),]), col="orange")
abline(lm(DPF3 ~ IRX5, data.tumor[which(data.tumor$rs4903064.col=="yellow"),]), col="yellow")

summary(lm(DPF3 ~ IRX2 * IRX5 * rs4903064.73279420.T.C, data.normal))

############################# DPF3 vs HYPOXIA #################################


data <- merge(data.gene, geno.rs4903064, by.x="name", by.y="Individual")
data <- merge(data, hypo.KIRC, by.x="Sample", by.y="Patient.tumor")

data <- mutate(data, rs4903064.col = case_when(rs4903064.73279420.T.C == 0 ~ "yellow", 
										 rs4903064.73279420.T.C == 1 ~ "orange", 
										 rs4903064.73279420.T.C == 2 ~ "red"))

############ tumor ############

data.tumor <- data[which(data$full.RNA %in% tumor),]
data.tumor <- data.tumor[complete.cases(data.tumor),]

par(mfrow=c(1,2))

plot(data.tumor$hypoxia, data.tumor$DPF3, pch = 19, col="blue", ylab="DPF3 Normalized Counts", xlab="Hypoxia Score", main="Tumor Renal (TCGA)")
abline(lm(DPF3 ~ hypoxia, data.tumor), col="blue")

plot(data.tumor$hypoxia, data.tumor$DPF3, pch = 19, col=data.tumor$rs4903064.col, , ylab="DPF3 Normalized Counts", xlab="Hypoxia Score", main="DPF3 vs Hypoxia Score (TCGA KIRC)")
abline(lm(DPF3 ~ hypoxia, data.tumor[which(data.tumor$rs4903064.col=="red"),]), col="red")
abline(lm(DPF3 ~ hypoxia, data.tumor[which(data.tumor$rs4903064.col=="orange"),]), col="orange")
abline(lm(DPF3 ~ hypoxia, data.tumor[which(data.tumor$rs4903064.col=="yellow"),]), col="yellow")
hypo.coef.TT <- lm(DPF3 ~ hypoxia, data.tumor[which(data.tumor$rs4903064.col=="yellow"),])$coefficients[2]
hypo.coef.CT <- lm(DPF3 ~ hypoxia, data.tumor[which(data.tumor$rs4903064.col=="orange"),])$coefficients[2]
hypo.coef.CC <- lm(DPF3 ~ hypoxia, data.tumor[which(data.tumor$rs4903064.col=="red"),])$coefficients[2]
title = paste0("beta.TT=", round(hypo.coef.TT, 2), "; beta.TC=", round(hypo.coef.CT, 2), "; beta.CC=", round(hypo.coef.CC, 2))
title(title, line=0.5, cex.main=0.9)
legend("topleft", legend=c("TT", "TC", "CC"), title="rs4903064", 
       col=c("yellow", "orange", "red"), pch=19, cex=0.8,
       box.lty=0)

dev.print(pdf, "Plot_DPF3_hypoxia.pdf")




summary(lm(DPF3 ~ rs4903064.73279420.T.C + hypoxia + rs4903064.73279420.T.C * hypoxia + IRX2 * IRX5 , data.tumor))
