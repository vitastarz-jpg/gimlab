setwd("C:/Users/dygks/OneDrive/πŸ≈¡ »≠∏È/Namsu/MSN/Research/US PROPEL_pilot RNA seq/03_Rrunning")

rm(list=ls())
getwd
#Loading libraries #1
library(anomalize)
library(forecast)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readxl)

#Loading libraries #2
library(tsibbledata)
library(feasts)
library(tibble)

#Loading libraries #3
library(data.table); library(dplyr); library(matrixStats); library(pheatmap); library(corrplot); library(tidyr); library(progress); library(ggplot2)

#Loading libraries #4
#install.packages("BiocManager")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("enrichplot")
library(AnnotationDbi); library(org.Hs.eg.db)
library(ggplot2); library(ggpubr); library(gridExtra)
library(DOSE); library(enrichplot); library(enrichR)

#Loading libraries #5
library(rpart); library(rpart.plot); library(randomForest); set.seed(312)

get_mode <- function(v) {  uniqv <- unique(v);   uniqv[which.max(tabulate(match(v, uniqv)))]}
replace_na <- function(x) {  x[x %in% c("NaN", "missing")] <- NA;   return(x)}
custom_col <- colorRampPalette(c("#03AED2", "#F7D940", "#FF0080"))(200)
####################################################################################
PROPEL_ML1 <- as.data.frame(read_excel("PROPEL/PROPEL US data_0728.xlsx", sheet=1)); names(PROPEL_ML1)[1] <- "ID"; PROPEL_ML1$ID <- substr(PROPEL_ML1$ID,1,9); PROPEL_ML1$ID <- gsub("\\-","_",PROPEL_ML1$ID)
PROPEL_ML2 <- as.data.frame(read_excel("PROPEL/PROPEL US data_0728.xlsx", sheet=2)); names(PROPEL_ML2)[1] <- "ID"; PROPEL_ML2$ID <- substr(PROPEL_ML2$ID,1,9); PROPEL_ML2$ID <- gsub("\\-","_",PROPEL_ML2$ID)
PROPEL_ML3 <- as.data.frame(read_excel("PROPEL/PROPEL US data_0728.xlsx", sheet=3)); names(PROPEL_ML3)[1] <- "ID"; PROPEL_ML3$ID <- substr(PROPEL_ML3$ID,1,9); PROPEL_ML3$ID <- gsub("\\-","_",PROPEL_ML3$ID)
PROPEL_ML4 <- as.data.frame(read_excel("PROPEL/PROPEL US data_0728.xlsx", sheet=4)); names(PROPEL_ML4)[1] <- "ID"; PROPEL_ML4$ID <- substr(PROPEL_ML4$ID,1,9); PROPEL_ML4$ID <- gsub("\\-","_",PROPEL_ML4$ID)
#Effect, non-effect group Íµ¨Î∂Ñ
PROPEL03c <- as.data.frame(read_excel("PROPEL/PROPEL_US_data_20240519.xlsx", sheet=4))
PROPEL03c <- dplyr::filter(PROPEL03c, grepl("PROPEL",ID))
rownames(PROPEL03c) <- PROPEL03c$ID
PROPEL03c[] <- lapply(PROPEL03c[], as.factor)
## 3Î≤?, 8Î≤àÏó¥??Ä C1B, C2B?ù∏?ç∞ A??Ä Í∞ôÏ?Ä Í∑∏Î£π?úºÎ°? ?ï†?ãπ?êòÍ∏? ?ïåÎ¨∏Ïóê ?Ç≠?†ú?ï®.
PROPEL03c <- PROPEL03c[,-c(3,8)]
head(PROPEL03c)

#Í≤∞Ï∏°Í∞? Ï≤òÎ¶¨
PROPEL_ML1 <- dplyr::mutate(PROPEL_ML1, across(everything(), replace_na))
PROPEL_ML2 <- dplyr::mutate(PROPEL_ML2, across(everything(), replace_na))
PROPEL_ML3 <- dplyr::mutate(PROPEL_ML3, across(everything(), replace_na))
PROPEL_ML4 <- dplyr::mutate(PROPEL_ML4, across(everything(), replace_na))
PROPEL_ML1[,1:2]; PROPEL_ML2[,1:2]; PROPEL_ML3[,1:2]; PROPEL_ML4[,1:2]
head (PROPEL_ML1)
head (PROPEL_ML2)

NACOUNT1 <- vector(); for(i in 1:ncol(PROPEL_ML1)){  NACOUNT1[i] <- sum(is.na(PROPEL_ML1[,i]))}
NACOUNT2 <- vector(); for(i in 1:ncol(PROPEL_ML2)){  NACOUNT2[i] <- sum(is.na(PROPEL_ML2[,i]))}
NACOUNT3 <- vector(); for(i in 1:ncol(PROPEL_ML3)){  NACOUNT3[i] <- sum(is.na(PROPEL_ML3[,i]))}
NACOUNT4 <- vector(); for(i in 1:ncol(PROPEL_ML4)){  NACOUNT4[i] <- sum(is.na(PROPEL_ML4[,i]))}
names(NACOUNT1) <- colnames(PROPEL_ML1)
names(NACOUNT2) <- colnames(PROPEL_ML2)
names(NACOUNT3) <- colnames(PROPEL_ML3)
names(NACOUNT4) <- colnames(PROPEL_ML4)
write.table(as.data.frame(NACOUNT1), file="clipboard",sep="\t",quote=FALSE,col.names=NA)
write.table(as.data.frame(NACOUNT2), file="clipboard",sep="\t",quote=FALSE,col.names=NA)
write.table(as.data.frame(NACOUNT3), file="clipboard",sep="\t",quote=FALSE,col.names=NA)
write.table(as.data.frame(NACOUNT4), file="clipboard",sep="\t",quote=FALSE,col.names=NA)

#ID ?ó¥ ?†ú?ô∏?ïòÍ≥? ?ÇòÎ®∏Ï?Ä ?ó¥ ?à´?ûêÎ°? Î≥Ä?ôò
##NAÍ∞íÎèÑ ?à´?ûêÎ°? Ï≤òÎ¶¨?ê®. ?ù¥ÎØ? ?ïû?óê?Ñú Î™®Îëê ?à´?ûêÎ°? Î≥Ä?ôò?ñàÍ∏? ?ïåÎ¨∏Ïóê Í∞Ä?ä•
PROPEL_ML1[,-1] <- lapply(PROPEL_ML1[,-1], as.numeric)
PROPEL_ML2[,-1] <- lapply(PROPEL_ML2[,-1], as.numeric)
PROPEL_ML3[,-1] <- lapply(PROPEL_ML3[,-1], as.numeric)
PROPEL_ML4[,-1] <- lapply(PROPEL_ML4[,-1], as.numeric)
str(PROPEL_ML1)
str(PROPEL_ML2)
str(PROPEL_ML4)
#Med op current, med ?ó¥ ?Ç≠?†ú (d/t Í≤∞Ï∏°Í∞? ?Üë)
PROPEL_ML1 <- dplyr::select(PROPEL_ML1, -meds_op_current, -curr_meds)
colnames (PROPEL_ML1)

#Í≤∞Ï∏°Í∞íÏùÑ ?èâÍ∑? Í∞íÏúºÎ°? ??ÄÏ≤?
PROPEL_ML1$age <- ifelse(is.na(PROPEL_ML1$age), mean(PROPEL_ML1$age, na.rm = TRUE), PROPEL_ML1$age)
PROPEL_ML1$qst_vdt_c_avg1 <- ifelse(is.na(PROPEL_ML1$qst_vdt_c_avg1), mean(PROPEL_ML1$qst_vdt_c_avg1, na.rm = TRUE), PROPEL_ML1$qst_vdt_c_avg1)
PROPEL_ML1$qst_vdt_p_avg1 <- ifelse(is.na(PROPEL_ML1$qst_vdt_p_avg1), mean(PROPEL_ML1$qst_vdt_p_avg1, na.rm = TRUE), PROPEL_ML1$qst_vdt_p_avg1)
PROPEL_ML1$lbp_unemployed <- ifelse(is.na(PROPEL_ML1$lbp_unemployed), get_mode(PROPEL_ML1$lbp_unemployed), PROPEL_ML1$lbp_unemployed)
PROPEL_ML1$lbp_workcomp   <- ifelse(is.na(PROPEL_ML1$lbp_workcomp  ), get_mode(PROPEL_ML1$lbp_workcomp  ), PROPEL_ML1$lbp_workcomp  )
#Í≤∞Ï∏°Í∞? ??ÄÏ≤? ?êò?óà?äîÏßÄ ?ôï?ù∏
head(PROPEL_ML1$age)
sum(is.na(PROPEL_ML1$lbp_unemployed))

PROPEL_ML2$age <- ifelse(is.na(PROPEL_ML2$age), mean(PROPEL_ML2$age, na.rm = TRUE), PROPEL_ML2$age)
PROPEL_ML2$qst_vdt_c_avg1 <- ifelse(is.na(PROPEL_ML2$qst_vdt_c_avg1), mean(PROPEL_ML2$qst_vdt_c_avg1, na.rm = TRUE), PROPEL_ML2$qst_vdt_c_avg1)
PROPEL_ML2$qst_vdt_p_avg1 <- ifelse(is.na(PROPEL_ML2$qst_vdt_p_avg1), mean(PROPEL_ML2$qst_vdt_p_avg1, na.rm = TRUE), PROPEL_ML2$qst_vdt_p_avg1)

PROPEL_ML4$qst_vdt_c_avg1 <- ifelse(is.na(PROPEL_ML4$qst_vdt_c_avg1), mean(PROPEL_ML4$qst_vdt_c_avg1, na.rm = TRUE), PROPEL_ML4$qst_vdt_c_avg1)
PROPEL_ML4$qst_vdt_p_avg1 <- ifelse(is.na(PROPEL_ML4$qst_vdt_p_avg1), mean(PROPEL_ML4$qst_vdt_p_avg1, na.rm = TRUE), PROPEL_ML4$qst_vdt_p_avg1)
sum(is.na(PROPEL_ML4$qst_vdt_p_avg1))

#E or NE group ?ó¨Î∂ÄÎ•? ML file?óê Î≥ëÌï©
PROPEL_ML1_C1A <- merge(x=PROPEL_ML1,y=PROPEL03c[,c(1,2)], by="ID"); rownames(PROPEL_ML1_C1A) <- PROPEL_ML1_C1A$ID; PROPEL_ML1_C1A <- PROPEL_ML1_C1A[,-1]
PROPEL_ML2_C1A <- merge(x=PROPEL_ML2,y=PROPEL03c[,c(1,2)], by="ID"); rownames(PROPEL_ML2_C1A) <- PROPEL_ML2_C1A$ID; PROPEL_ML2_C1A <- PROPEL_ML2_C1A[,-1]
PROPEL_ML3_C1A <- merge(x=PROPEL_ML3,y=PROPEL03c[,c(1,2)], by="ID"); rownames(PROPEL_ML3_C1A) <- PROPEL_ML3_C1A$ID; PROPEL_ML3_C1A <- PROPEL_ML3_C1A[,-1]
PROPEL_ML4_C1A <- merge(x=PROPEL_ML4,y=PROPEL03c[,c(1,2)], by="ID"); rownames(PROPEL_ML4_C1A) <- PROPEL_ML4_C1A$ID; PROPEL_ML4_C1A <- PROPEL_ML4_C1A[,-1]
PROPEL_ML1_C1C <- merge(x=PROPEL_ML1,y=PROPEL03c[,c(1,3)], by="ID"); rownames(PROPEL_ML1_C1C) <- PROPEL_ML1_C1C$ID; PROPEL_ML1_C1C <- PROPEL_ML1_C1C[,-1]
PROPEL_ML2_C1C <- merge(x=PROPEL_ML2,y=PROPEL03c[,c(1,3)], by="ID"); rownames(PROPEL_ML2_C1C) <- PROPEL_ML2_C1C$ID; PROPEL_ML2_C1C <- PROPEL_ML2_C1C[,-1]
PROPEL_ML3_C1C <- merge(x=PROPEL_ML3,y=PROPEL03c[,c(1,3)], by="ID"); rownames(PROPEL_ML3_C1C) <- PROPEL_ML3_C1C$ID; PROPEL_ML3_C1C <- PROPEL_ML3_C1C[,-1]
PROPEL_ML4_C1C <- merge(x=PROPEL_ML4,y=PROPEL03c[,c(1,3)], by="ID"); rownames(PROPEL_ML4_C1C) <- PROPEL_ML4_C1C$ID; PROPEL_ML4_C1C <- PROPEL_ML4_C1C[,-1]
PROPEL_ML1_C1D <- merge(x=PROPEL_ML1,y=PROPEL03c[,c(1,4)], by="ID"); rownames(PROPEL_ML1_C1D) <- PROPEL_ML1_C1D$ID; PROPEL_ML1_C1D <- PROPEL_ML1_C1D[,-1]
PROPEL_ML2_C1D <- merge(x=PROPEL_ML2,y=PROPEL03c[,c(1,4)], by="ID"); rownames(PROPEL_ML2_C1D) <- PROPEL_ML2_C1D$ID; PROPEL_ML2_C1D <- PROPEL_ML2_C1D[,-1]
PROPEL_ML3_C1D <- merge(x=PROPEL_ML3,y=PROPEL03c[,c(1,4)], by="ID"); rownames(PROPEL_ML3_C1D) <- PROPEL_ML3_C1D$ID; PROPEL_ML3_C1D <- PROPEL_ML3_C1D[,-1]
PROPEL_ML4_C1D <- merge(x=PROPEL_ML4,y=PROPEL03c[,c(1,4)], by="ID"); rownames(PROPEL_ML4_C1D) <- PROPEL_ML4_C1D$ID; PROPEL_ML4_C1D <- PROPEL_ML4_C1D[,-1]
PROPEL_ML1_C1E <- merge(x=PROPEL_ML1,y=PROPEL03c[,c(1,5)], by="ID"); rownames(PROPEL_ML1_C1E) <- PROPEL_ML1_C1E$ID; PROPEL_ML1_C1E <- PROPEL_ML1_C1E[,-1]
PROPEL_ML2_C1E <- merge(x=PROPEL_ML2,y=PROPEL03c[,c(1,5)], by="ID"); rownames(PROPEL_ML2_C1E) <- PROPEL_ML2_C1E$ID; PROPEL_ML2_C1E <- PROPEL_ML2_C1E[,-1]
PROPEL_ML3_C1E <- merge(x=PROPEL_ML3,y=PROPEL03c[,c(1,5)], by="ID"); rownames(PROPEL_ML3_C1E) <- PROPEL_ML3_C1E$ID; PROPEL_ML3_C1E <- PROPEL_ML3_C1E[,-1]
PROPEL_ML4_C1E <- merge(x=PROPEL_ML4,y=PROPEL03c[,c(1,5)], by="ID"); rownames(PROPEL_ML4_C1E) <- PROPEL_ML4_C1E$ID; PROPEL_ML4_C1E <- PROPEL_ML4_C1E[,-1]
PROPEL_ML1_C2A <- merge(x=PROPEL_ML1,y=PROPEL03c[,c(1,6)], by="ID"); rownames(PROPEL_ML1_C2A) <- PROPEL_ML1_C2A$ID; PROPEL_ML1_C2A <- PROPEL_ML1_C2A[,-1]
PROPEL_ML2_C2A <- merge(x=PROPEL_ML2,y=PROPEL03c[,c(1,6)], by="ID"); rownames(PROPEL_ML2_C2A) <- PROPEL_ML2_C2A$ID; PROPEL_ML2_C2A <- PROPEL_ML2_C2A[,-1]
PROPEL_ML3_C2A <- merge(x=PROPEL_ML3,y=PROPEL03c[,c(1,6)], by="ID"); rownames(PROPEL_ML3_C2A) <- PROPEL_ML3_C2A$ID; PROPEL_ML3_C2A <- PROPEL_ML3_C2A[,-1]
PROPEL_ML4_C2A <- merge(x=PROPEL_ML4,y=PROPEL03c[,c(1,6)], by="ID"); rownames(PROPEL_ML4_C2A) <- PROPEL_ML4_C2A$ID; PROPEL_ML4_C2A <- PROPEL_ML4_C2A[,-1]
PROPEL_ML1_C2C <- merge(x=PROPEL_ML1,y=PROPEL03c[,c(1,7)], by="ID"); rownames(PROPEL_ML1_C2C) <- PROPEL_ML1_C2C$ID; PROPEL_ML1_C2C <- PROPEL_ML1_C2C[,-1]
PROPEL_ML2_C2C <- merge(x=PROPEL_ML2,y=PROPEL03c[,c(1,7)], by="ID"); rownames(PROPEL_ML2_C2C) <- PROPEL_ML2_C2C$ID; PROPEL_ML2_C2C <- PROPEL_ML2_C2C[,-1]
PROPEL_ML3_C2C <- merge(x=PROPEL_ML3,y=PROPEL03c[,c(1,7)], by="ID"); rownames(PROPEL_ML3_C2C) <- PROPEL_ML3_C2C$ID; PROPEL_ML3_C2C <- PROPEL_ML3_C2C[,-1]
PROPEL_ML4_C2C <- merge(x=PROPEL_ML4,y=PROPEL03c[,c(1,7)], by="ID"); rownames(PROPEL_ML4_C2C) <- PROPEL_ML4_C2C$ID; PROPEL_ML4_C2C <- PROPEL_ML4_C2C[,-1]
PROPEL_ML1_C2D <- merge(x=PROPEL_ML1,y=PROPEL03c[,c(1,8)], by="ID"); rownames(PROPEL_ML1_C2D) <- PROPEL_ML1_C2D$ID; PROPEL_ML1_C2D <- PROPEL_ML1_C2D[,-1]
PROPEL_ML2_C2D <- merge(x=PROPEL_ML2,y=PROPEL03c[,c(1,8)], by="ID"); rownames(PROPEL_ML2_C2D) <- PROPEL_ML2_C2D$ID; PROPEL_ML2_C2D <- PROPEL_ML2_C2D[,-1]
PROPEL_ML3_C2D <- merge(x=PROPEL_ML3,y=PROPEL03c[,c(1,8)], by="ID"); rownames(PROPEL_ML3_C2D) <- PROPEL_ML3_C2D$ID; PROPEL_ML3_C2D <- PROPEL_ML3_C2D[,-1]
PROPEL_ML4_C2D <- merge(x=PROPEL_ML4,y=PROPEL03c[,c(1,8)], by="ID"); rownames(PROPEL_ML4_C2D) <- PROPEL_ML4_C2D$ID; PROPEL_ML4_C2D <- PROPEL_ML4_C2D[,-1]
PROPEL_ML1_C2E <- merge(x=PROPEL_ML1,y=PROPEL03c[,c(1,9)], by="ID"); rownames(PROPEL_ML1_C2E) <- PROPEL_ML1_C2E$ID; PROPEL_ML1_C2E <- PROPEL_ML1_C2E[,-1]
PROPEL_ML2_C2E <- merge(x=PROPEL_ML2,y=PROPEL03c[,c(1,9)], by="ID"); rownames(PROPEL_ML2_C2E) <- PROPEL_ML2_C2E$ID; PROPEL_ML2_C2E <- PROPEL_ML2_C2E[,-1]
PROPEL_ML3_C2E <- merge(x=PROPEL_ML3,y=PROPEL03c[,c(1,9)], by="ID"); rownames(PROPEL_ML3_C2E) <- PROPEL_ML3_C2E$ID; PROPEL_ML3_C2E <- PROPEL_ML3_C2E[,-1]
PROPEL_ML4_C2E <- merge(x=PROPEL_ML4,y=PROPEL03c[,c(1,9)], by="ID"); rownames(PROPEL_ML4_C2E) <- PROPEL_ML4_C2E$ID; PROPEL_ML4_C2E <- PROPEL_ML4_C2E[,-1]
head(PROPEL_ML2_C1A)


#Radnom forest??Ä Decision Tree ML Î∂ÑÏÑù
getwd()
dir.exists("Output_20240729")
#Íµ¨Ï°∞ ?ôï?ù∏?ïò?ó¨ yÍ∞? (C1A?ù¥ N?ù∏Í∞Ä EN)?ù¥ Î≤îÏ£º?òï?ù¥ ?ïÑ?ãà?ùºÎ©? ?ïÑ?ûò ?ãù ?Ç¨?ö©?ïò?ó¨ Î≤îÏ£º?òï?úºÎ°? Î≥Ä?ôò ?ïÑ?öî?öî
str(PROPEL_ML1_C1A$C1A)
PROPEL_ML1_C1A$C1A <- as.factor(PROPEL_ML1_C1A$C1A)
PROPEL_ML2_C1A$C1A <- as.factor(PROPEL_ML2_C1A$C1A)
PROPEL_ML1_C1C$C1C <- as.factor(PROPEL_ML1_C1C$C1C)
PROPEL_ML2_C1C$C1C <- as.factor(PROPEL_ML2_C1C$C1C)

pdf("Output_20240729/03_ML1_C1A.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML1_1AAd1 <- rpart(C1A ~ ., data = PROPEL_ML1_C1A); rpart.plot(PROPEL_ML1_1AAd1, digits=3, type=2, extra=1); PROPEL_ML1_C1Ar1 <- randomForest(C1A ~ ., PROPEL_ML1_C1A); varImpPlot(PROPEL_ML1_C1Ar1); dev.off()
pdf("Output_20240729/03_ML2_C1A.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML2_1AAd1 <- rpart(C1A ~ ., data = PROPEL_ML2_C1A); rpart.plot(PROPEL_ML2_1AAd1, digits=3, type=2, extra=1); PROPEL_ML2_C1Ar1 <- randomForest(C1A ~ ., PROPEL_ML2_C1A); varImpPlot(PROPEL_ML2_C1Ar1); dev.off()
pdf("Output_20240729/03_ML3_C1A.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML3_1AAd1 <- rpart(C1A ~ ., data = PROPEL_ML3_C1A); rpart.plot(PROPEL_ML3_1AAd1, digits=3, type=2, extra=1); PROPEL_ML3_C1Ar1 <- randomForest(C1A ~ ., PROPEL_ML3_C1A); varImpPlot(PROPEL_ML3_C1Ar1); dev.off()
pdf("Output_20240729/03_ML4_C1A.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML4_1AAd1 <- rpart(C1A ~ ., data = PROPEL_ML4_C1A); rpart.plot(PROPEL_ML4_1AAd1, digits=3, type=2, extra=1); PROPEL_ML4_C1Ar1 <- randomForest(C1A ~ ., PROPEL_ML4_C1A); varImpPlot(PROPEL_ML4_C1Ar1); dev.off()
pdf("Output_20240729/03_ML1_C1C.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML1_1CAd1 <- rpart(C1C ~ ., data = PROPEL_ML1_C1C); rpart.plot(PROPEL_ML1_1CAd1, digits=3, type=2, extra=1); PROPEL_ML1_C1Cr1 <- randomForest(C1C ~ ., PROPEL_ML1_C1C); varImpPlot(PROPEL_ML1_C1Cr1); dev.off()
pdf("Output_20240729/03_ML2_C1C.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML2_1CAd1 <- rpart(C1C ~ ., data = PROPEL_ML2_C1C); rpart.plot(PROPEL_ML2_1CAd1, digits=3, type=2, extra=1); PROPEL_ML2_C1Cr1 <- randomForest(C1C ~ ., PROPEL_ML2_C1C); varImpPlot(PROPEL_ML2_C1Cr1); dev.off()
pdf("Output_20240729/03_ML3_C1C.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML3_1CAd1 <- rpart(C1C ~ ., data = PROPEL_ML3_C1C); rpart.plot(PROPEL_ML3_1CAd1, digits=3, type=2, extra=1); PROPEL_ML3_C1Cr1 <- randomForest(C1C ~ ., PROPEL_ML3_C1C); varImpPlot(PROPEL_ML3_C1Cr1); dev.off()
pdf("Output_20240729/03_ML4_C1C.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML4_1CAd1 <- rpart(C1C ~ ., data = PROPEL_ML4_C1C); rpart.plot(PROPEL_ML4_1CAd1, digits=3, type=2, extra=1); PROPEL_ML4_C1Cr1 <- randomForest(C1C ~ ., PROPEL_ML4_C1C); varImpPlot(PROPEL_ML4_C1Cr1); dev.off()
pdf("Output_20240729/03_ML1_C1D.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML1_1DAd1 <- rpart(C1D ~ ., data = PROPEL_ML1_C1D); rpart.plot(PROPEL_ML1_1DAd1, digits=3, type=2, extra=1); PROPEL_ML1_C1Dr1 <- randomForest(C1D ~ ., PROPEL_ML1_C1D); varImpPlot(PROPEL_ML1_C1Dr1); dev.off()
pdf("Output_20240729/03_ML2_C1D.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML2_1DAd1 <- rpart(C1D ~ ., data = PROPEL_ML2_C1D); rpart.plot(PROPEL_ML2_1DAd1, digits=3, type=2, extra=1); PROPEL_ML2_C1Dr1 <- randomForest(C1D ~ ., PROPEL_ML2_C1D); varImpPlot(PROPEL_ML2_C1Dr1); dev.off()
pdf("Output_20240729/03_ML3_C1D.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML3_1DAd1 <- rpart(C1D ~ ., data = PROPEL_ML3_C1D); rpart.plot(PROPEL_ML3_1DAd1, digits=3, type=2, extra=1); PROPEL_ML3_C1Dr1 <- randomForest(C1D ~ ., PROPEL_ML3_C1D); varImpPlot(PROPEL_ML3_C1Dr1); dev.off()
pdf("Output_20240729/03_ML4_C1D.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML4_1DAd1 <- rpart(C1D ~ ., data = PROPEL_ML4_C1D); rpart.plot(PROPEL_ML4_1DAd1, digits=3, type=2, extra=1); PROPEL_ML4_C1Dr1 <- randomForest(C1D ~ ., PROPEL_ML4_C1D); varImpPlot(PROPEL_ML4_C1Dr1); dev.off()
pdf("Output_20240729/03_ML1_C1E.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML1_1EAd1 <- rpart(C1E ~ ., data = PROPEL_ML1_C1E); rpart.plot(PROPEL_ML1_1EAd1, digits=3, type=2, extra=1); PROPEL_ML1_C1Er1 <- randomForest(C1E ~ ., PROPEL_ML1_C1E); varImpPlot(PROPEL_ML1_C1Er1); dev.off()
pdf("Output_20240729/03_ML2_C1E.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML2_1EAd1 <- rpart(C1E ~ ., data = PROPEL_ML2_C1E); rpart.plot(PROPEL_ML2_1EAd1, digits=3, type=2, extra=1); PROPEL_ML2_C1Er1 <- randomForest(C1E ~ ., PROPEL_ML2_C1E); varImpPlot(PROPEL_ML2_C1Er1); dev.off()
pdf("Output_20240729/03_ML3_C1E.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML3_1EAd1 <- rpart(C1E ~ ., data = PROPEL_ML3_C1E); rpart.plot(PROPEL_ML3_1EAd1, digits=3, type=2, extra=1); PROPEL_ML3_C1Er1 <- randomForest(C1E ~ ., PROPEL_ML3_C1E); varImpPlot(PROPEL_ML3_C1Er1); dev.off()
pdf("Output_20240729/03_ML4_C1E.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML4_1EAd1 <- rpart(C1E ~ ., data = PROPEL_ML4_C1E); rpart.plot(PROPEL_ML4_1EAd1, digits=3, type=2, extra=1); PROPEL_ML4_C1Er1 <- randomForest(C1E ~ ., PROPEL_ML4_C1E); varImpPlot(PROPEL_ML4_C1Er1); dev.off()
pdf("Output_20240729/03_ML1_C2A.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML1_2AAd1 <- rpart(C2A ~ ., data = PROPEL_ML1_C2A); rpart.plot(PROPEL_ML1_2AAd1, digits=3, type=2, extra=1); PROPEL_ML1_C2Ar1 <- randomForest(C2A ~ ., PROPEL_ML1_C2A); varImpPlot(PROPEL_ML1_C2Ar1); dev.off()
pdf("Output_20240729/03_ML2_C2A.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML2_2AAd1 <- rpart(C2A ~ ., data = PROPEL_ML2_C2A); rpart.plot(PROPEL_ML2_2AAd1, digits=3, type=2, extra=1); PROPEL_ML2_C2Ar1 <- randomForest(C2A ~ ., PROPEL_ML2_C2A); varImpPlot(PROPEL_ML2_C2Ar1); dev.off()
pdf("Output_20240729/03_ML3_C2A.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML3_2AAd1 <- rpart(C2A ~ ., data = PROPEL_ML3_C2A); rpart.plot(PROPEL_ML3_2AAd1, digits=3, type=2, extra=1); PROPEL_ML3_C2Ar1 <- randomForest(C2A ~ ., PROPEL_ML3_C2A); varImpPlot(PROPEL_ML3_C2Ar1); dev.off()
pdf("Output_20240729/03_ML4_C2A.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML4_2AAd1 <- rpart(C2A ~ ., data = PROPEL_ML4_C2A); rpart.plot(PROPEL_ML4_2AAd1, digits=3, type=2, extra=1); PROPEL_ML4_C2Ar1 <- randomForest(C2A ~ ., PROPEL_ML4_C2A); varImpPlot(PROPEL_ML4_C2Ar1); dev.off()
pdf("Output_20240729/03_ML1_C2C.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML1_2CAd1 <- rpart(C2C ~ ., data = PROPEL_ML1_C2C); rpart.plot(PROPEL_ML1_2CAd1, digits=3, type=2, extra=1); PROPEL_ML1_C2Cr1 <- randomForest(C2C ~ ., PROPEL_ML1_C2C); varImpPlot(PROPEL_ML1_C2Cr1); dev.off()
pdf("Output_20240729/03_ML2_C2C.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML2_2CAd1 <- rpart(C2C ~ ., data = PROPEL_ML2_C2C); rpart.plot(PROPEL_ML2_2CAd1, digits=3, type=2, extra=1); PROPEL_ML2_C2Cr1 <- randomForest(C2C ~ ., PROPEL_ML2_C2C); varImpPlot(PROPEL_ML2_C2Cr1); dev.off()
pdf("Output_20240729/03_ML3_C2C.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML3_2CAd1 <- rpart(C2C ~ ., data = PROPEL_ML3_C2C); rpart.plot(PROPEL_ML3_2CAd1, digits=3, type=2, extra=1); PROPEL_ML3_C2Cr1 <- randomForest(C2C ~ ., PROPEL_ML3_C2C); varImpPlot(PROPEL_ML3_C2Cr1); dev.off()
pdf("Output_20240729/03_ML4_C2C.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML4_2CAd1 <- rpart(C2C ~ ., data = PROPEL_ML4_C2C); rpart.plot(PROPEL_ML4_2CAd1, digits=3, type=2, extra=1); PROPEL_ML4_C2Cr1 <- randomForest(C2C ~ ., PROPEL_ML4_C2C); varImpPlot(PROPEL_ML4_C2Cr1); dev.off()
pdf("Output_20240729/03_ML1_C2D.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML1_2DAd1 <- rpart(C2D ~ ., data = PROPEL_ML1_C2D); rpart.plot(PROPEL_ML1_2DAd1, digits=3, type=2, extra=1); PROPEL_ML1_C2Dr1 <- randomForest(C2D ~ ., PROPEL_ML1_C2D); varImpPlot(PROPEL_ML1_C2Dr1); dev.off()
pdf("Output_20240729/03_ML2_C2D.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML2_2DAd1 <- rpart(C2D ~ ., data = PROPEL_ML2_C2D); rpart.plot(PROPEL_ML2_2DAd1, digits=3, type=2, extra=1); PROPEL_ML2_C2Dr1 <- randomForest(C2D ~ ., PROPEL_ML2_C2D); varImpPlot(PROPEL_ML2_C2Dr1); dev.off()
pdf("Output_20240729/03_ML3_C2D.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML3_2DAd1 <- rpart(C2D ~ ., data = PROPEL_ML3_C2D); rpart.plot(PROPEL_ML3_2DAd1, digits=3, type=2, extra=1); PROPEL_ML3_C2Dr1 <- randomForest(C2D ~ ., PROPEL_ML3_C2D); varImpPlot(PROPEL_ML3_C2Dr1); dev.off()
pdf("Output_20240729/03_ML4_C2D.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML4_2DAd1 <- rpart(C2D ~ ., data = PROPEL_ML4_C2D); rpart.plot(PROPEL_ML4_2DAd1, digits=3, type=2, extra=1); PROPEL_ML4_C2Dr1 <- randomForest(C2D ~ ., PROPEL_ML4_C2D); varImpPlot(PROPEL_ML4_C2Dr1); dev.off()
pdf("Output_20240729/03_ML1_C2E.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML1_2EAd1 <- rpart(C2E ~ ., data = PROPEL_ML1_C2E); rpart.plot(PROPEL_ML1_2EAd1, digits=3, type=2, extra=1); PROPEL_ML1_C2Er1 <- randomForest(C2E ~ ., PROPEL_ML1_C2E); varImpPlot(PROPEL_ML1_C2Er1); dev.off()
pdf("Output_20240729/03_ML2_C2E.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML2_2EAd1 <- rpart(C2E ~ ., data = PROPEL_ML2_C2E); rpart.plot(PROPEL_ML2_2EAd1, digits=3, type=2, extra=1); PROPEL_ML2_C2Er1 <- randomForest(C2E ~ ., PROPEL_ML2_C2E); varImpPlot(PROPEL_ML2_C2Er1); dev.off()
pdf("Output_20240729/03_ML3_C2E.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML3_2EAd1 <- rpart(C2E ~ ., data = PROPEL_ML3_C2E); rpart.plot(PROPEL_ML3_2EAd1, digits=3, type=2, extra=1); PROPEL_ML3_C2Er1 <- randomForest(C2E ~ ., PROPEL_ML3_C2E); varImpPlot(PROPEL_ML3_C2Er1); dev.off()
pdf("Output_20240729/03_ML4_C2E.pdf", width=15,height=8); par(mfrow=c(1,2));PROPEL_ML4_2EAd1 <- rpart(C2E ~ ., data = PROPEL_ML4_C2E); rpart.plot(PROPEL_ML4_2EAd1, digits=3, type=2, extra=1); PROPEL_ML4_C2Er1 <- randomForest(C2E ~ ., PROPEL_ML4_C2E); varImpPlot(PROPEL_ML4_C2Er1); dev.off()
file.exists("Output_20240729/03_ML1_C1A.pdf")

####################################################################################
#Load RNA-seq data, then make dataframe as expression levels
##data frame?ùÑ ÎßåÎì§?ñ¥ Î™®Îì† RNA-seq count ?åå?ùº?ùÑ Î≥ëÌï©
RNAseq_ALL01 <- data.frame()
head(RNAseq_ALL01)
IBS_HC_01_01_S51 <- as.data.frame(fread("RNAseq/IBS-HC-01-01_S51.counts", header = FALSE, sep = "\t")); RNAseq_ALL01 <- IBS_HC_01_01_S51         ; rm(IBS_HC_01_01_S51)
IBS_HC_02_01_S29 <- as.data.frame(fread("RNAseq/IBS-HC-02-01_S29.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[, 4] <- IBS_HC_02_01_S29[,3]; rm(IBS_HC_02_01_S29)
IBS_HC_03_01_S13 <- as.data.frame(fread("RNAseq/IBS-HC-03-01_S13.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[, 5] <- IBS_HC_03_01_S13[,3]; rm(IBS_HC_03_01_S13)
IBS_HC_04_01_S05 <- as.data.frame(fread("RNAseq/IBS-HC-04-01_S5.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[, 6] <- IBS_HC_04_01_S05[,3]; rm(IBS_HC_04_01_S05)
IBS_HC_05_01_S43 <- as.data.frame(fread("RNAseq/IBS-HC-05-01_S43.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[, 7] <- IBS_HC_05_01_S43[,3]; rm(IBS_HC_05_01_S43)
IBS_HC_06_01_S36 <- as.data.frame(fread("RNAseq/IBS-HC-06-01_S36.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[, 8] <- IBS_HC_06_01_S36[,3]; rm(IBS_HC_06_01_S36)
IBS_HC_07_01_S04 <- as.data.frame(fread("RNAseq/IBS-HC-07-01_S4.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[, 9] <- IBS_HC_07_01_S04[,3]; rm(IBS_HC_07_01_S04)
IBS_HC_08_01_S12 <- as.data.frame(fread("RNAseq/IBS-HC-08-01_S12.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,10] <- IBS_HC_08_01_S12[,3]; rm(IBS_HC_08_01_S12)
IBS_HC_09_01_S56 <- as.data.frame(fread("RNAseq/IBS-HC-09-01_S56.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,11] <- IBS_HC_09_01_S56[,3]; rm(IBS_HC_09_01_S56)
IBS_HC_10_01_S49 <- as.data.frame(fread("RNAseq/IBS-HC-10-01_S49.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,12] <- IBS_HC_10_01_S49[,3]; rm(IBS_HC_10_01_S49)
IBS_HC_11_01_S44 <- as.data.frame(fread("RNAseq/IBS-HC-11-01_S44.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,13] <- IBS_HC_11_01_S44[,3]; rm(IBS_HC_11_01_S44)
IBS_HC_12_01_S37 <- as.data.frame(fread("RNAseq/IBS-HC-12-01_S37.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,14] <- IBS_HC_12_01_S37[,3]; rm(IBS_HC_12_01_S37)
IBS_HC_13_01_S21 <- as.data.frame(fread("RNAseq/IBS-HC-13-01_S21.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,15] <- IBS_HC_13_01_S21[,3]; rm(IBS_HC_13_01_S21)
IBS_HC_14_01_S57 <- as.data.frame(fread("RNAseq/IBS-HC-14-01_S57.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,16] <- IBS_HC_14_01_S57[,3]; rm(IBS_HC_14_01_S57)
IBS_HC_15_01_S50 <- as.data.frame(fread("RNAseq/IBS-HC-15-01_S50.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,17] <- IBS_HC_15_01_S50[,3]; rm(IBS_HC_15_01_S50)
IBS_HC_16_01_S28 <- as.data.frame(fread("RNAseq/IBS-HC-16-01_S28.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,18] <- IBS_HC_16_01_S28[,3]; rm(IBS_HC_16_01_S28)
IBS_HC_17_01_S20 <- as.data.frame(fread("RNAseq/IBS-HC-17-01_S20.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,19] <- IBS_HC_17_01_S20[,3]; rm(IBS_HC_17_01_S20)
IBS_HC_18_01_S35 <- as.data.frame(fread("RNAseq/IBS-HC-18-01_S35.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,20] <- IBS_HC_18_01_S35[,3]; rm(IBS_HC_18_01_S35)
IBS_HC_19_01_S27 <- as.data.frame(fread("RNAseq/IBS-HC-19-01_S27.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,21] <- IBS_HC_19_01_S27[,3]; rm(IBS_HC_19_01_S27)
IBS_HC_20_01_S42 <- as.data.frame(fread("RNAseq/IBS-HC-20-01_S42.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,22] <- IBS_HC_20_01_S42[,3]; rm(IBS_HC_20_01_S42)
IBS_HC_21_01_S58 <- as.data.frame(fread("RNAseq/IBS-HC-21-01_S58.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,23] <- IBS_HC_21_01_S58[,3]; rm(IBS_HC_21_01_S58)
PROPEL_01_01_S26 <- as.data.frame(fread("RNAseq/PROPEL-01-01_S26.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,24] <- PROPEL_01_01_S26[,3]; rm(PROPEL_01_01_S26)
PROPEL_01_07_S01 <- as.data.frame(fread("RNAseq/PROPEL-01-07_S1.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,25] <- PROPEL_01_07_S01[,3]; rm(PROPEL_01_07_S01)
PROPEL_02_01_S32 <- as.data.frame(fread("RNAseq/PROPEL-02-01_S32.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,26] <- PROPEL_02_01_S32[,3]; rm(PROPEL_02_01_S32)
PROPEL_02_07_S19 <- as.data.frame(fread("RNAseq/PROPEL-02-07_S19.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,27] <- PROPEL_02_07_S19[,3]; rm(PROPEL_02_07_S19)
PROPEL_03_01_S24 <- as.data.frame(fread("RNAseq/PROPEL-03-01_S24.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,28] <- PROPEL_03_01_S24[,3]; rm(PROPEL_03_01_S24)
PROPEL_03_07_S09 <- as.data.frame(fread("RNAseq/PROPEL-03-07_S9.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,29] <- PROPEL_03_07_S09[,3]; rm(PROPEL_03_07_S09)
PROPEL_04_01_S16 <- as.data.frame(fread("RNAseq/PROPEL-04-01_S16.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,30] <- PROPEL_04_01_S16[,3]; rm(PROPEL_04_01_S16)
PROPEL_05_01_S08 <- as.data.frame(fread("RNAseq/PROPEL-05-01_S8.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,31] <- PROPEL_05_01_S08[,3]; rm(PROPEL_05_01_S08)
PROPEL_06_07_S11 <- as.data.frame(fread("RNAseq/PROPEL-06-07_S11.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,32] <- PROPEL_06_07_S11[,3]; rm(PROPEL_06_07_S11)
PROPEL_07_01_S60 <- as.data.frame(fread("RNAseq/PROPEL-07-01_S60.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,33] <- PROPEL_07_01_S60[,3]; rm(PROPEL_07_01_S60)
PROPEL_07_07_S17 <- as.data.frame(fread("RNAseq/PROPEL-07-07_S17.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,34] <- PROPEL_07_07_S17[,3]; rm(PROPEL_07_07_S17)
PROPEL_09_01_S53 <- as.data.frame(fread("RNAseq/PROPEL-09-01_S53.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,35] <- PROPEL_09_01_S53[,3]; rm(PROPEL_09_01_S53)
PROPEL_10_01_S59 <- as.data.frame(fread("RNAseq/PROPEL-10-01_S59.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,36] <- PROPEL_10_01_S59[,3]; rm(PROPEL_10_01_S59)
PROPEL_10_07_S25 <- as.data.frame(fread("RNAseq/PROPEL-10-07_S25.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,37] <- PROPEL_10_07_S25[,3]; rm(PROPEL_10_07_S25)
PROPEL_11_01_S38 <- as.data.frame(fread("RNAseq/PROPEL-11-01_S38.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,38] <- PROPEL_11_01_S38[,3]; rm(PROPEL_11_01_S38)
PROPEL_12_01_S46 <- as.data.frame(fread("RNAseq/PROPEL-12-01_S46.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,39] <- PROPEL_12_01_S46[,3]; rm(PROPEL_12_01_S46)
PROPEL_12_07_S33 <- as.data.frame(fread("RNAseq/PROPEL-12-07_S33.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,40] <- PROPEL_12_07_S33[,3]; rm(PROPEL_12_07_S33)
PROPEL_13_01_S39 <- as.data.frame(fread("RNAseq/PROPEL-13-01_S39.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,41] <- PROPEL_13_01_S39[,3]; rm(PROPEL_13_01_S39)
PROPEL_13_07_S40 <- as.data.frame(fread("RNAseq/PROPEL-13-07_S40.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,42] <- PROPEL_13_07_S40[,3]; rm(PROPEL_13_07_S40)
PROPEL_14_01_S31 <- as.data.frame(fread("RNAseq/PROPEL-14-01_S31.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,43] <- PROPEL_14_01_S31[,3]; rm(PROPEL_14_01_S31)
PROPEL_15_01_S23 <- as.data.frame(fread("RNAseq/PROPEL-15-01_S23.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,44] <- PROPEL_15_01_S23[,3]; rm(PROPEL_15_01_S23)
PROPEL_15_07_S03 <- as.data.frame(fread("RNAseq/PROPEL-15-07_S3.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,45] <- PROPEL_15_07_S03[,3]; rm(PROPEL_15_07_S03)
PROPEL_16_01_S07 <- as.data.frame(fread("RNAseq/PROPEL-16-01_S7.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,46] <- PROPEL_16_01_S07[,3]; rm(PROPEL_16_01_S07)
PROPEL_16_07_S55 <- as.data.frame(fread("RNAseq/PROPEL-16-07_S55.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,47] <- PROPEL_16_07_S55[,3]; rm(PROPEL_16_07_S55)
PROPEL_17_01_S15 <- as.data.frame(fread("RNAseq/PROPEL-17-01_S15.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,48] <- PROPEL_17_01_S15[,3]; rm(PROPEL_17_01_S15)
PROPEL_17_07_S47 <- as.data.frame(fread("RNAseq/PROPEL-17-07_S47.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,49] <- PROPEL_17_07_S47[,3]; rm(PROPEL_17_07_S47)
PROPEL_18_07_S48 <- as.data.frame(fread("RNAseq/PROPEL-18-07_S48.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,50] <- PROPEL_18_07_S48[,3]; rm(PROPEL_18_07_S48)
PROPEL_19_01_S45 <- as.data.frame(fread("RNAseq/PROPEL-19-01_S45.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,51] <- PROPEL_19_01_S45[,3]; rm(PROPEL_19_01_S45)
PROPEL_19_07_S54 <- as.data.frame(fread("RNAseq/PROPEL-19-07_S54.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,52] <- PROPEL_19_07_S54[,3]; rm(PROPEL_19_07_S54)
PROPEL_20_01_S52 <- as.data.frame(fread("RNAseq/PROPEL-20-01_S52.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,53] <- PROPEL_20_01_S52[,3]; rm(PROPEL_20_01_S52)
PROPEL_20_07_S41 <- as.data.frame(fread("RNAseq/PROPEL-20-07_S41.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,54] <- PROPEL_20_07_S41[,3]; rm(PROPEL_20_07_S41)
PROPEL_21_01_S30 <- as.data.frame(fread("RNAseq/PROPEL-21-01_S30.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,55] <- PROPEL_21_01_S30[,3]; rm(PROPEL_21_01_S30)
PROPEL_21_07_S02 <- as.data.frame(fread("RNAseq/PROPEL-21-07_S2.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,56] <- PROPEL_21_07_S02[,3]; rm(PROPEL_21_07_S02)
PROPEL_23_01_S22 <- as.data.frame(fread("RNAseq/PROPEL-23-01_S22.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,57] <- PROPEL_23_01_S22[,3]; rm(PROPEL_23_01_S22)
PROPEL_23_07_S10 <- as.data.frame(fread("RNAseq/PROPEL-23-07_S10.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,58] <- PROPEL_23_07_S10[,3]; rm(PROPEL_23_07_S10)
PROPEL_24_01_S14 <- as.data.frame(fread("RNAseq/PROPEL-24-01_S14.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,59] <- PROPEL_24_01_S14[,3]; rm(PROPEL_24_01_S14)
PROPEL_24_07_S34 <- as.data.frame(fread("RNAseq/PROPEL-24-07_S34.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,60] <- PROPEL_24_07_S34[,3]; rm(PROPEL_24_07_S34)
PROPEL_25_01_S06 <- as.data.frame(fread("RNAseq/PROPEL-25-01_S6.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,61] <- PROPEL_25_01_S06[,3]; rm(PROPEL_25_01_S06)
PROPEL_25_07_S18 <- as.data.frame(fread("RNAseq/PROPEL-25-07_S18.counts", header = FALSE, sep = "\t")); RNAseq_ALL01[,62] <- PROPEL_25_07_S18[,3]; rm(PROPEL_25_07_S18)
PROPEL_27_01_S14 <- as.data.frame(fread("RNAseq/PROPEL02701_S14.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,63] <- PROPEL_27_01_S14[,3]; rm(PROPEL_27_01_S14)
PROPEL_28_01_S15 <- as.data.frame(fread("RNAseq/PROPEL02801_S15.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,64] <- PROPEL_28_01_S15[,3]; rm(PROPEL_28_01_S15)
PROPEL_29_01_S01 <- as.data.frame(fread("RNAseq/PROPEL02901_S1.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,65] <- PROPEL_29_01_S01[,3]; rm(PROPEL_29_01_S01)
PROPEL_29_07_S02 <- as.data.frame(fread("RNAseq/PROPEL02907_S2.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,66] <- PROPEL_29_07_S02[,3]; rm(PROPEL_29_07_S02)
PROPEL_30_01_S16 <- as.data.frame(fread("RNAseq/PROPEL03001_S16.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,67] <- PROPEL_30_01_S16[,3]; rm(PROPEL_30_01_S16)
PROPEL_31_01_S03 <- as.data.frame(fread("RNAseq/PROPEL03101_S3.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,68] <- PROPEL_31_01_S03[,3]; rm(PROPEL_31_01_S03)
PROPEL_31_07_S04 <- as.data.frame(fread("RNAseq/PROPEL03107_S4.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,69] <- PROPEL_31_07_S04[,3]; rm(PROPEL_31_07_S04)
PROPEL_32_01_S05 <- as.data.frame(fread("RNAseq/PROPEL03201_S5.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,70] <- PROPEL_32_01_S05[,3]; rm(PROPEL_32_01_S05)
PROPEL_32_07_S06 <- as.data.frame(fread("RNAseq/PROPEL03207_S6.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,71] <- PROPEL_32_07_S06[,3]; rm(PROPEL_32_07_S06)
PROPEL_33_01_S07 <- as.data.frame(fread("RNAseq/PROPEL03301_S7.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,72] <- PROPEL_33_01_S07[,3]; rm(PROPEL_33_01_S07)
PROPEL_34_01_S08 <- as.data.frame(fread("RNAseq/PROPEL03401_S8.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,73] <- PROPEL_34_01_S08[,3]; rm(PROPEL_34_01_S08)
PROPEL_34_07_S09 <- as.data.frame(fread("RNAseq/PROPEL03407_S9.counts"  , header = FALSE, sep = "\t")); RNAseq_ALL01[,74] <- PROPEL_34_07_S09[,3]; rm(PROPEL_34_07_S09)
PROPEL_38_01_S10 <- as.data.frame(fread("RNAseq/PROPEL03801_S10.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,75] <- PROPEL_38_01_S10[,3]; rm(PROPEL_38_01_S10)
PROPEL_38_07_S11 <- as.data.frame(fread("RNAseq/PROPEL03807_S11.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,76] <- PROPEL_38_07_S11[,3]; rm(PROPEL_38_07_S11)
PROPEL_39_01_S17 <- as.data.frame(fread("RNAseq/PROPEL03901_S17.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,77] <- PROPEL_39_01_S17[,3]; rm(PROPEL_39_01_S17)
PROPEL_40_01_S12 <- as.data.frame(fread("RNAseq/PROPEL04001_S12.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,78] <- PROPEL_40_01_S12[,3]; rm(PROPEL_40_01_S12)
PROPEL_40_07_S13 <- as.data.frame(fread("RNAseq/PROPEL04007_S13.counts" , header = FALSE, sep = "\t")); RNAseq_ALL01[,79] <- PROPEL_40_07_S13[,3]; rm(PROPEL_40_07_S13)
head(RNAseq_ALL01)
names(RNAseq_ALL01) <- c("ENSEMBL","SYMBOL","IBS_HC_01_01_S51","IBS_HC_02_01_S29","IBS_HC_03_01_S13","IBS_HC_04_01_S05","IBS_HC_05_01_S43","IBS_HC_06_01_S36","IBS_HC_07_01_S04","IBS_HC_08_01_S12","IBS_HC_09_01_S56","IBS_HC_10_01_S49","IBS_HC_11_01_S44","IBS_HC_12_01_S37","IBS_HC_13_01_S21","IBS_HC_14_01_S57","IBS_HC_15_01_S50","IBS_HC_16_01_S28","IBS_HC_17_01_S20","IBS_HC_18_01_S35","IBS_HC_19_01_S27","IBS_HC_20_01_S42","IBS_HC_21_01_S58","PROPEL_01_01_S26","PROPEL_01_07_S01","PROPEL_02_01_S32","PROPEL_02_07_S19","PROPEL_03_01_S24","PROPEL_03_07_S09","PROPEL_04_01_S16","PROPEL_05_01_S08","PROPEL_06_07_S11","PROPEL_07_01_S60","PROPEL_07_07_S17","PROPEL_09_01_S53","PROPEL_10_01_S59","PROPEL_10_07_S25","PROPEL_11_01_S38","PROPEL_12_01_S46","PROPEL_12_07_S33","PROPEL_13_01_S39","PROPEL_13_07_S40","PROPEL_14_01_S31","PROPEL_15_01_S23","PROPEL_15_07_S03","PROPEL_16_01_S07","PROPEL_16_07_S55","PROPEL_17_01_S15","PROPEL_17_07_S47","PROPEL_18_07_S48","PROPEL_19_01_S45","PROPEL_19_07_S54","PROPEL_20_01_S52","PROPEL_20_07_S41","PROPEL_21_01_S30","PROPEL_21_07_S02","PROPEL_23_01_S22","PROPEL_23_07_S10","PROPEL_24_01_S14","PROPEL_24_07_S34","PROPEL_25_01_S06","PROPEL_25_07_S18","PROPEL_27_01_S14","PROPEL_28_01_S15","PROPEL_29_01_S01","PROPEL_29_07_S02","PROPEL_30_01_S16","PROPEL_31_01_S03","PROPEL_31_07_S04","PROPEL_32_01_S05","PROPEL_32_07_S06","PROPEL_33_01_S07","PROPEL_34_01_S08","PROPEL_34_07_S09","PROPEL_38_01_S10","PROPEL_38_07_S11","PROPEL_39_01_S17","PROPEL_40_01_S12","PROPEL_40_07_S13")
head(RNAseq_ALL01); dim(RNAseq_ALL01)

PROPEL01 <- as.data.frame(read_excel("PROPEL/PROPEL_US_data_20240519.xlsx", sheet=2))
#PROPELÏ∞∏Í?Ä?ûê Î≤àÌò∏ÎßåÏùÑ Ï∂îÏ∂ú!-> length(unique(substr(PROPEL01$subject_id,1,9)))
PROPEL01[1:5,1:5]; dim(PROPEL01); length(unique(substr(PROPEL01$subject_id,1,9)))
dplyr::filter(PROPEL01, substr(subject_id,12,12)==1)[,1:11]
PROPEL01[1:5,1:5]; dim(PROPEL01) #238 152

#"-"?ùÑ "_"Î°? Î≥Ä?ôò
PROPEL01$record_id  <- gsub("\\-","_", PROPEL01$record_id)
PROPEL01$subject_id <- gsub("\\-","_", PROPEL01$subject_id)
#RNAseq_ALL01?óê?Ñú 1,2Î≤àÏß∏ ?ó¥?ùÑ ?Ñ†?Éù?ïò?ó¨ ?ç∞?ù¥?Ñ∞ ?îÑ?†à?ûÑ?ùÑ ÎßåÎì¨
rownames(RNAseq_ALL01) <- RNAseq_ALL01$ENSEMBL
RNAseq_GENEs <- RNAseq_ALL01[,1:2]
head(RNAseq_GENEs)

#edgeR ?å®?Ç§ÏßÄ?óê?Ñú ?Ç¨?ö©?ïò?äî DEGlist ?ç∞?ù¥?Ñ∞ Íµ¨Ï°∞Î°? Î≥Ä?ôò
library(edgeR)
#Count ?ïÑ?ìú?óê ??Ä?û•?ïò?ó¨ count ?àòÎß? ?Çò?ò§Í≤? ?ï®.
RNAseq_ALL02 <- DGEList(counts=RNAseq_ALL01[,3:79])
head(RNAseq_ALL02)
#DEGlist count ?ç∞?ù¥?Ñ∞Î•? ?†ïÍ∑úÌôî?ïò?ó¨ ÎπÑÍµê Í∞Ä?ä•?ïòÍ≤? ?ï®.
RNAseq_ALL02 <- calcNormFactors(RNAseq_ALL02)
##?†ïÍ∑úÌôî Í≥ÑÏàò ?ôï?ù∏
RNAseq_ALL02$samples
##Counts per Million?ùÑ Í≥ÑÏÇ∞?ïò?ó¨ ?ùΩÍ∏? ?àòÎ•? Î∞±Îßå?ã®?úÑÎ°? ?†ïÍ∑úÌôî (Î≥¥Ï†ï)
RNAseq_ALL02 <- cpm(RNAseq_ALL02, normalized.lib.sizes=TRUE)
#log Í∞íÏúºÎ°? Î≥Ä?ôò > "0"?ù¥ ?ûà?ùÑ Í≤ΩÏö∞ Î°úÍ∑∏ Í≥ÑÏÇ∞ Î∂àÍ?Ä?ä•, Î∂ÑÌè¨ ?†ïÍ∑úÌôî Î∞? ?ù¥?ÉÅÏπ? Ï§ÑÏù¥Í∏? ?úÑ?ï®.
###!!! ?ïÑ?ûò ?ç∞?ù¥?Ñ∞Î•? ?ôú?ö©?ï¥?Ñú T1?ùò baseline Í∞íÏúºÎ°? ML Î∂ÑÏÑù?óê ?Ñ£?ñ¥Î≥¥Í∏∞!!!!!
RNAseq_ALL02 <- as.data.frame(log10(RNAseq_ALL02+1))
head(RNAseq_ALL02); dim(RNAseq_ALL02) #?ú†?†Ñ?ûê ?àò: 66028,?Éò?îå ?àò: 77
#?ú†?†Ñ?ûê?óê ?î∞Î•? Î∞úÌòÑ?üâ heatmap
pheatmap(RNAseq_ALL02[1:100,], show_rownames = FALSE)

########################################################################################
#ML Î∂ÑÏÑù
RNA_T1_samples <- grep("_01_", colnames(RNAseq_ALL02), value = TRUE)
length(RNA_T1_samples)  # Î™? Í∞úÏùò T1 ?Éò?îå?ù¥ ?ûà?äîÏßÄ ?ôï?ù∏

# ?òà: PROPEL_01_01_S26 -> PROPEL_01
convert_to_ID <- function(sample_name) {
  sub("(_01_.*)", "", sample_name)
}
RNA_T1_IDs <- sapply(RNA_T1_samples, convert_to_ID)
RNA_T1_expr <- RNAseq_ALL02[, RNA_T1_samples]
colnames(RNA_T1_expr) <- RNA_T1_IDs
RNA_T1_expr <- as.data.frame(t(RNA_T1_expr))  # ?ñâ: ?Éò?îå, ?ó¥: ?ú†?†Ñ?ûê
RNA_T1_expr$ID <- rownames(RNA_T1_expr)
ML2_with_RNA <- merge(PROPEL_ML2, RNA_T1_expr, by = "ID")

colnames(PROPEL03c)
PROPEL03c$ID <- rownames(PROPEL03c)
ML2_C2D <- merge(ML2_with_RNA, PROPEL03c[, c("ID", "C2D")], by = "ID")
ML2_C2D$C2D <- factor(ML2_C2D$C2D, levels = c("N", "E"))  # ÎπÑÎ∞ò?ùëÍµ? N, Î∞òÏùëÍµ? E

ML2_C2E <- merge(ML2_with_RNA, PROPEL03c[, c("ID", "C2E")], by = "ID")
ML2_C2E$C2E <- factor(ML2_C2E$C2E, levels = c("N", "E"))

library(randomForest)
set.seed(123)
#?ú†?†Ñ?ûê ?àòÍ∞Ä ?ÑàÎ¨? ÎßéÏïÑ ?óê?ü¨
rf_model <- randomForest(C2D ~ ., data = ML2_C2D[,-1], importance = TRUE, ntree = 500)

head(PROPEL03c$ID)
head(colnames(RNAseq_ALL02))

T1_samples <- colnames(RNAseq_ALL02)[grepl("_01_", colnames(RNAseq_ALL02))]
RNA_ID_map <- data.frame(
  SampleID = T1_samples,
  ID = gsub("(_01_.*$)", "", T1_samples)  # ?òà: "PROPEL_01_01_S26" ?Üí "PROPEL_01"
)

rna_df <- as.data.frame(t(RNAseq_ALL02[, T1_samples]))
rna_df$SampleID <- rownames(rna_df)

# ID Ï∂îÍ?Ä
rna_df <- merge(rna_df, RNA_ID_map, by = "SampleID")
ML2_with_RNA <- merge(PROPEL_ML2, rna_df[, -1], by = "ID")  # SampleID ?†úÍ±?

ML2_C2D <- merge(ML2_with_RNA, PROPEL03c[, c("ID", "C2D")], by = "ID")
table(ML2_C2D$C2D)


#?ÉÅ?úÑ 100Í∞úÎßå
# 1. ?èâÍ∑? Î∞úÌòÑ?üâ Í∏∞Ï?Ä ?ÉÅ?úÑ 100Í∞? ?ú†?†Ñ?ûê ?Ñ†?Éù
rna_only <- rna_df[, !(colnames(rna_df) %in% c("SampleID", "ID"))]
gene_means <- colMeans(rna_only)
top100_genes <- names(sort(gene_means, decreasing = TRUE))[1:100]

# 2. ?ïÑ?öî?ïú Ïª¨ÎüºÎß? Ï∂îÏ∂ú?ïò?ó¨ ?ã§?ãú Î≥ëÌï©
rna_top100 <- rna_df[, c("ID", top100_genes)]
ML2_with_RNA_top100 <- merge(PROPEL_ML2, rna_top100, by = "ID")
ML2_C2D_top100 <- merge(ML2_with_RNA_top100, PROPEL03c[, c("ID", "C2D")], by = "ID")

# 3. ?ûú?ç§?è¨?†à?ä§?ä∏ ?ã§?ñâ
library(randomForest)
rf_model <- randomForest(C2D ~ ., data = ML2_C2D_top100[,-1], importance = TRUE, ntree = 500)

# Ï§ëÏöî ?ú†?†Ñ?ûê ?ôï?ù∏
importance(rf_model)
# 1. Ï¢ÖÏÜçÎ≥Ä?àò factor Î≥Ä?ôò
ML2_C2D_top100$C2D <- as.factor(ML2_C2D_top100$C2D)

# 2. NAÍ∞Ä ?ûà?äîÏßÄ ?ôï?ù∏
sum(is.na(ML2_C2D_top100))  # Í≤∞Í≥ºÍ∞Ä 0?ù¥ ?ïÑ?ãàÎ©? NA Ï°¥Ïû¨

# NAÍ∞Ä ?ûà?ã§Î©? ?†úÍ±?
ML2_C2D_top100 <- na.omit(ML2_C2D_top100)

# 3. Î≥Ä?àòÎ™? Î¨∏Ï†ú ?ï¥Í≤?
colnames(ML2_C2D_top100) <- make.names(colnames(ML2_C2D_top100), unique = TRUE)

# 4. ?ã§?ãú Î™®Îç∏ ?ã§?ñâ
rf_model <- randomForest(C2D ~ ., data = ML2_C2D_top100[,-1], importance = TRUE, ntree = 500)

importance(rf_model)

varImpPlot(rf_model)

colnames(ML2_C2D_top100)[grepl("X_|Qual", colnames(ML2_C2D_top100))]
# Î∂àÌïÑ?öî?ïú Î≥Ä?àò ?†úÍ±?
ML2_C2D_top100_clean <- ML2_C2D_top100[, !grepl("X_|Qual", colnames(ML2_C2D_top100))]

# ?ã§?ãú Î™®Îç∏ ?Éù?Ñ±
ML2_C2D_top100_clean$C2D <- as.factor(ML2_C2D_top100_clean$C2D)
rf_model_clean <- randomForest(C2D ~ ., data = ML2_C2D_top100_clean[,-1], importance = TRUE, ntree = 500)
importance(rf_model)

varImpPlot(rf_model)

# importance Í∞ùÏ≤¥ Î∂àÎü¨?ò§Í∏?
imp_raw <- importance(rf_model_clean)

# ENSEMBL ?ù¥Î¶? Ï∂îÏ∂ú ?Üí .?à´?ûê ?†úÍ±?
ensembl_names <- rownames(imp_raw)
ensembl_clean <- gsub("\\.\\d+$", "", ensembl_names)

# SYMBOL ?†ïÎ≥? Î∂ôÏù¥Í∏?
RNAseq_GENEs$ENSEMBL_clean <- gsub("\\.\\d+$", "", RNAseq_GENEs$ENSEMBL)
gene_map <- RNAseq_GENEs[, c("ENSEMBL_clean", "SYMBOL")]

# SYMBOL Îß§Ïπ≠
symbol_names <- gene_map$SYMBOL[match(ensembl_clean, gene_map$ENSEMBL_clean)]

# SYMBOL ?óÜ?äî Í≤ΩÏö∞?äî ?õê?ûò ?ù¥Î¶? ?ú†ÏßÄ
final_labels <- ifelse(is.na(symbol_names), ensembl_names, symbol_names)

# rownames ÍµêÏ≤¥
rownames(imp_raw) <- final_labels

# ?ù¥?†ú varImpPlot ?ã§?ñâ
varImpPlot(rf_model_clean, n.var = 20, main = "Top 20 Important Variables")

########################################################################################
#ALL02 ?åå?ùº?óê?Ñú PROPEL_01_01_S26 -> PROPEL_01_01Î°? Î≥ÄÍ≤? 
RNAseq_ALL03 <- dplyr::select(RNAseq_ALL02, contains("PROPEL")); colnames(RNAseq_ALL03) <- substr(colnames(RNAseq_ALL03),1,12)
head(RNAseq_ALL03)
unique(substr(colnames(RNAseq_ALL03),8,9))
#PROPEL01 subject ID??Ä ALL03?óê ?ùºÏπòÌïò?äî IDÎ•? ?Ñ†?Éù?ïò?ó¨ COMMON IDÎ°? ??Ä?û• > ALL03?óê Î∂ÑÏÑù?ï† timepointÎß? ??Ä?û•, timepoint 02~06?ùÑ ?Ç≠?†ú?ïòÍ∏? ?úÑ?ï®.
head(PROPEL01$subject_id)
head(RNAseq_ALL03)
COMMON_ID <- intersect(PROPEL01$subject_id, colnames(RNAseq_ALL03))
PROPEL03 <- dplyr::filter(PROPEL01, subject_id %in% COMMON_ID)
str(PROPEL03)
head(PROPEL03)
RNAseq_ALL03 <- dplyr::select(RNAseq_ALL03, COMMON_ID)
head(RNAseq_ALL03)

dim(PROPEL03) #54,720
dim(RNAseq_ALL03) # 66028,54

#2Î≤? ?ù¥?ÉÅ Í≤ÄÏ≤? Ï±ÑÏ∑®Î•? ?ïú ??Ä?ÉÅ?ûêÎß? Ï∂îÏ∂ú.
COMMON_ID03 <- as.data.frame(table(substr(COMMON_ID,1,9)))
COMMON_ID03 <- dplyr::filter(COMMON_ID03, Freq == 2)
head(COMMON_ID03); dim(COMMON_ID03)
#?Åù?óê 01, 07Î∂ôÏó¨?Ñú  timepoint Íµ¨Î∂Ñ
RNAseq_ALL03 <- dplyr::select(RNAseq_ALL03, c(paste0(COMMON_ID03$Var1,"_01"),paste0(COMMON_ID03$Var1,"_07")))
head(RNAseq_ALL03); dim(RNAseq_ALL03) #66028, 40

#PROPEL03Í≥? RNAseq_ALL03?óê?Ñú subject ID Ï§ëÎ≥µ?êò?äî Í≤ÉÎßå Ï∂îÏ∂ú 
PROPEL03 <- dplyr::filter(PROPEL03, subject_id %in% colnames(RNAseq_ALL03))
PROPEL03[1:5, 1:19]; dim(PROPEL03) #40 152
PROPEL03

#PROPEL03?óê Effect or Non-effect group ?ó¨Î∂Ä Ï∂îÍ?Ä
PROPEL03c <- as.data.frame(read_excel("PROPEL/PROPEL_US_data_20240519.xlsx", sheet=4))
PROPEL03c <- dplyr::filter(PROPEL03c, grepl("PROPEL",ID))
head(PROPEL03c); tail(PROPEL03c)

#ID?ó¥?ùÑ ?ñâ ?ù¥Î¶ÑÏúºÎ°? ?Ñ§?†ï?ïòÍ≥? ID ?ó¥ ?Ç≠?†ú > Ï∞®Îì±Î∞úÌòÑ Î∞? PCAÍ∞Ä Í∞Ä?ä•?ï®.
rownames(PROPEL03c) <- PROPEL03c$ID; PROPEL03c <- PROPEL03c[,-1]
#Î™®Îì† ?ó¥?ùÑ Î≤îÏ£º?òï ?ç∞?ù¥?Ñ∞Î°? Î≥Ä?ôò> Ï¢ÖÏÜçÎ≥Ä?àò?ù∏ E or NE group?ù¥ Î≤îÏ£º?òï?ù¥Í∏? ?ïåÎ¨??
PROPEL03c[] <- lapply(PROPEL03c[], as.factor)
str(PROPEL03c)
head(PROPEL03c)

#?Åù?óê 01, 07?ùÑ ?Ç≠?†ú?ïòÎ©¥ÏÑú T1Í≥? T7 ?ç∞?ù¥?Ñ∞ Íµ¨Î∂Ñ?ïò?ó¨ ??Ä?û•
RNAseq_ALL03_T1 <- dplyr::select(RNAseq_ALL03, ends_with("_01")); colnames(RNAseq_ALL03_T1) <- substr(colnames(RNAseq_ALL03_T1),1,9); dim(RNAseq_ALL03_T1)
RNAseq_ALL03_T7 <- dplyr::select(RNAseq_ALL03, ends_with("_07")); colnames(RNAseq_ALL03_T7) <- substr(colnames(RNAseq_ALL03_T7),1,9); dim(RNAseq_ALL03_T7)
RNAseq_ALL03_T1 <- as.matrix(dplyr::select(RNAseq_ALL03_T1, rownames(PROPEL03c))); dim(RNAseq_ALL03_T1) #66028,20
RNAseq_ALL03_T7 <- as.matrix(dplyr::select(RNAseq_ALL03_T7, rownames(PROPEL03c))); dim(RNAseq_ALL03_T7) #66028,20
table(colnames(RNAseq_ALL03_T1)==rownames(PROPEL03c)) #20
table(colnames(RNAseq_ALL03_T7)==rownames(PROPEL03c)) #20

#Effect or Non-effect group?ùÑ BPI Î∂ÑÎ•ò Í∏∞Ï?Ä?óê ?î∞?ùº ?ïÑ?Ñ∞Îß? > ALL03?óê?Ñú Ï§ëÎ≥µ?êò?äî ?ó¥Îß? Ï∂îÏ∂ú
RNAseq_ALL03_E1 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C1A == "E"))))); dim(RNAseq_ALL03_E1) #66028,18
RNAseq_ALL03_N1 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C1A == "N"))))); dim(RNAseq_ALL03_N1) #66028,22
RNAseq_ALL03_E3 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C1C == "E"))))); dim(RNAseq_ALL03_E3) #66028, 8
RNAseq_ALL03_N3 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C1C == "N"))))); dim(RNAseq_ALL03_N3) #66028,32
RNAseq_ALL03_E4 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C1D == "E"))))); dim(RNAseq_ALL03_E4) #66028,24
RNAseq_ALL03_N4 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C1D == "N"))))); dim(RNAseq_ALL03_N4) #66028,16
RNAseq_ALL03_E5 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C1E == "E"))))); dim(RNAseq_ALL03_E5) #66028,26
RNAseq_ALL03_N5 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C1E == "N"))))); dim(RNAseq_ALL03_N5) #66028,14

RNAseq_ALL03_E6 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C2A == "E"))))); dim(RNAseq_ALL03_E6) #66028,18
RNAseq_ALL03_N6 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C2A == "N"))))); dim(RNAseq_ALL03_N6) #66028,22
RNAseq_ALL03_E8 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C2C == "E"))))); dim(RNAseq_ALL03_E8) #66028,24
RNAseq_ALL03_N8 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C2C == "N"))))); dim(RNAseq_ALL03_N8) #66028,16
RNAseq_ALL03_E9 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C2D == "E"))))); dim(RNAseq_ALL03_E9) #66028,26
RNAseq_ALL03_N9 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C2D == "N"))))); dim(RNAseq_ALL03_N9) #66028,14
RNAseq_ALL03_E0 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C2E == "E"))))); dim(RNAseq_ALL03_E0) #66028,22
RNAseq_ALL03_N0 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C2E == "N"))))); dim(RNAseq_ALL03_N0) #66028,18

#Î™®Îì† ??Ä?ÉÅ?ûê Ï§ëÏû¨ ?†Ñ?õÑ ?ôï?ù∏
RNAseq_ALL03_A1 <- as.matrix(dplyr::select(RNAseq_ALL03, contains(rownames(dplyr::filter(PROPEL03c, C3A == "E"))))); dim(RNAseq_ALL03_A1) #66028,40

head(RNAseq_ALL03_E1)
#01??Ä T1?úºÎ°? 07??Ä T7?ù¥?ùº?äî Î≤îÏ£º?òï Î≥Ä?àòÎ•? ?Éù?Ñ±
E1f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_E1),11,12)))
head(E1f)
N1f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_N1),11,12)))
E3f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_E3),11,12)))
N3f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_N3),11,12)))
E4f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_E4),11,12)))
N4f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_N4),11,12)))
E5f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_E5),11,12)))
N5f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_N5),11,12)))
E6f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_E6),11,12)))
N6f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_N6),11,12)))
E8f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_E8),11,12)))
N8f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_N8),11,12)))

E9f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_E9),11,12)))
N9f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_N9),11,12)))
E0f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_E0),11,12)))
N0f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_N0),11,12)))
A1f <- as.factor(paste0("T",substr(colnames(RNAseq_ALL03_A1),11,12)))

#?ó∞Íµ? Î™©Ìëú 3 Í∞ÅÍ∞Å Í∑∏Î£π Í∞? T1-T7?ùò Ï∞®Ïù¥ Î∂ÑÏÑù
#?úÑ?óê?Ñú ?ú†?†Ñ?ûê Î∞úÌòÑ?üâ ?àú?úºÎ°? ?†ïÎ¶¨Ìïú ?Ç¥?ö©?ù¥ ?óÜ?äî?ç∞ 1?ñâ?óê ?ûà?äî ?ú†?†Ñ?ûêÎß? ÎπÑÍµê?ïú ?ù¥?ú†?äî?
t.test(RNAseq_ALL03_E1[1, ] ~ E1f)

t.test(RNAseq_ALL03_N1[1, ] ~ N1f)

t.test(RNAseq_ALL03_E3[1, ] ~ E3f)
t.test(RNAseq_ALL03_N3[1, ] ~ N3f)
t.test(RNAseq_ALL03_E4[1, ] ~ E4f)
t.test(RNAseq_ALL03_N4[1, ] ~ N4f)
t.test(RNAseq_ALL03_E5[1, ] ~ E5f)
t.test(RNAseq_ALL03_N5[1, ] ~ N5f)
t.test(RNAseq_ALL03_E6[1, ] ~ E6f)
t.test(RNAseq_ALL03_N6[1, ] ~ N6f)
t.test(RNAseq_ALL03_E8[1, ] ~ E8f)
t.test(RNAseq_ALL03_N8[1, ] ~ N8f)
t.test(RNAseq_ALL03_E9[1, ] ~ E9f)
t.test(RNAseq_ALL03_N9[1, ] ~ N9f)
t.test(RNAseq_ALL03_E0[1, ] ~ E0f)
t.test(RNAseq_ALL03_N0[1, ] ~ N0f)
t.test(RNAseq_ALL03_A1[1, ] ~ A1f)

#####¿ÃπÃ ∆ƒ¿œ ¿˙¿Â«ÿ µ◊¿∏¥œ æ»«ÿµµ µ ~601±Ó¡ˆ
########################?ó∞Íµ¨Î™©?ëú 1(?†ÑÏ≤? ??Ä?ÉÅ?ûê)?ùò ??Ä?ïú FC?äî ?óÜ?ùå.##################
#?ù¥ÎØ? ??Ä?û•?êò?ñ¥ ?ûà?úº?ãà ?ï† ?ïÑ?öî ?óÜ?ùå!!
PV_E1 <- vector(); FC_E1 <- vector(); for (i in 1:nrow(RNAseq_ALL03_E1)){TTEST <- t.test(RNAseq_ALL03_E1[i, ] ~ E1f); PV_E1[i] <- TTEST$p.value; FC_E1[i] <- diff(TTEST$estimate)}
PV_N1 <- vector(); FC_N1 <- vector(); for (i in 1:nrow(RNAseq_ALL03_N1)){TTEST <- t.test(RNAseq_ALL03_N1[i, ] ~ N1f); PV_N1[i] <- TTEST$p.value; FC_N1[i] <- diff(TTEST$estimate)}
PV_E3 <- vector(); FC_E3 <- vector(); for (i in 1:nrow(RNAseq_ALL03_E3)){TTEST <- t.test(RNAseq_ALL03_E3[i, ] ~ E3f); PV_E3[i] <- TTEST$p.value; FC_E3[i] <- diff(TTEST$estimate)}
PV_N3 <- vector(); FC_N3 <- vector(); for (i in 1:nrow(RNAseq_ALL03_N3)){TTEST <- t.test(RNAseq_ALL03_N3[i, ] ~ N3f); PV_N3[i] <- TTEST$p.value; FC_N3[i] <- diff(TTEST$estimate)}
PV_E4 <- vector(); FC_E4 <- vector(); for (i in 1:nrow(RNAseq_ALL03_E4)){TTEST <- t.test(RNAseq_ALL03_E4[i, ] ~ E4f); PV_E4[i] <- TTEST$p.value; FC_E4[i] <- diff(TTEST$estimate)}
PV_N4 <- vector(); FC_N4 <- vector(); for (i in 1:nrow(RNAseq_ALL03_N4)){TTEST <- t.test(RNAseq_ALL03_N4[i, ] ~ N4f); PV_N4[i] <- TTEST$p.value; FC_N4[i] <- diff(TTEST$estimate)}
PV_E5 <- vector(); FC_E5 <- vector(); for (i in 1:nrow(RNAseq_ALL03_E5)){TTEST <- t.test(RNAseq_ALL03_E5[i, ] ~ E5f); PV_E5[i] <- TTEST$p.value; FC_E5[i] <- diff(TTEST$estimate)}
PV_N5 <- vector(); FC_N5 <- vector(); for (i in 1:nrow(RNAseq_ALL03_N5)){TTEST <- t.test(RNAseq_ALL03_N5[i, ] ~ N5f); PV_N5[i] <- TTEST$p.value; FC_N5[i] <- diff(TTEST$estimate)}
PV_E6 <- vector(); FC_E6 <- vector(); for (i in 1:nrow(RNAseq_ALL03_E6)){TTEST <- t.test(RNAseq_ALL03_E6[i, ] ~ E6f); PV_E6[i] <- TTEST$p.value; FC_E6[i] <- diff(TTEST$estimate)}
PV_N6 <- vector(); FC_N6 <- vector(); for (i in 1:nrow(RNAseq_ALL03_N6)){TTEST <- t.test(RNAseq_ALL03_N6[i, ] ~ N6f); PV_N6[i] <- TTEST$p.value; FC_N6[i] <- diff(TTEST$estimate)}
PV_E8 <- vector(); FC_E8 <- vector(); for (i in 1:nrow(RNAseq_ALL03_E8)){TTEST <- t.test(RNAseq_ALL03_E8[i, ] ~ E8f); PV_E8[i] <- TTEST$p.value; FC_E8[i] <- diff(TTEST$estimate)}
PV_N8 <- vector(); FC_N8 <- vector(); for (i in 1:nrow(RNAseq_ALL03_N8)){TTEST <- t.test(RNAseq_ALL03_N8[i, ] ~ N8f); PV_N8[i] <- TTEST$p.value; FC_N8[i] <- diff(TTEST$estimate)}
PV_E9 <- vector(); FC_E9 <- vector(); for (i in 1:nrow(RNAseq_ALL03_E9)){TTEST <- t.test(RNAseq_ALL03_E9[i, ] ~ E9f); PV_E9[i] <- TTEST$p.value; FC_E9[i] <- diff(TTEST$estimate)}
PV_N9 <- vector(); FC_N9 <- vector(); for (i in 1:nrow(RNAseq_ALL03_N9)){TTEST <- t.test(RNAseq_ALL03_N9[i, ] ~ N9f); PV_N9[i] <- TTEST$p.value; FC_N9[i] <- diff(TTEST$estimate)}
PV_E0 <- vector(); FC_E0 <- vector(); for (i in 1:nrow(RNAseq_ALL03_E0)){TTEST <- t.test(RNAseq_ALL03_E0[i, ] ~ E0f); PV_E0[i] <- TTEST$p.value; FC_E0[i] <- diff(TTEST$estimate)}
PV_N0 <- vector(); FC_N0 <- vector(); for (i in 1:nrow(RNAseq_ALL03_N0)){TTEST <- t.test(RNAseq_ALL03_N0[i, ] ~ N0f); PV_N0[i] <- TTEST$p.value; FC_N0[i] <- diff(TTEST$estimate)}
PV_A1 <- vector(); FC_A1 <- vector(); for (i in 1:nrow(RNAseq_ALL03_A1)){TTEST <- t.test(RNAseq_ALL03_A1[i, ] ~ A1f); PV_A1[i] <- TTEST$p.value; FC_A1[i] <- diff(TTEST$estimate)}
PVFC_DEG03 <- data.frame(PV_E1,FC_E1,PV_N1,FC_N1,PV_E3,FC_E3,PV_N3,FC_N3,PV_E4,FC_E4,PV_N4,FC_N4,PV_E5,FC_E5,PV_N5,FC_N5,PV_E6,FC_E6,PV_N6,FC_N6,PV_E8,FC_E8,PV_N8,FC_N8,PV_E9,FC_E9,PV_N9,FC_N9,PV_E0,FC_E0,PV_N0,FC_N0, PV_A1, FC_A1)
rownames(PVFC_DEG03) <- rownames(RNAseq_ALL03_E1)
head(PVFC_DEG03)
write.table(PVFC_DEG03, file="PVFC_DEG03.txt",sep="\t",quote=FALSE,col.names=NA)

#Í∞? ÏßëÎã®?ùò T1 Î∞? T7?óê?Ñú RNA Î∞úÌòÑ?üâ Ï∞®Ïù¥ Î∞? FC Íµ¨ÌïòÍ∏?  (?ó∞Íµ¨Î™©?ëú 2Ôº?)
##?ù¥ÎØ? ?åå?ùº??Ä Running ?è¥?çî?óê ??Ä?û•?êò?ñ¥ ?ûà?ùå.
### ?ûú?ç§?ïòÍ≤? Ï≤´Î≤àÏß? ?ñâ?óê ?ûà?äî ?ú†?†Ñ?ûê?ùò Î∞úÌòÑ?üâ Ï∞®Ïù¥?ùº ?Å∞ ?ùòÎØ? ?óÜ?ùå.
head (RNAseq_ALL03_T1)
t.test(RNAseq_ALL03_T1[1, ] ~ PROPEL03c$C1A)
t.test(RNAseq_ALL03_T1[1, ] ~ PROPEL03c$C1C)
t.test(RNAseq_ALL03_T1[1, ] ~ PROPEL03c$C1D)
t.test(RNAseq_ALL03_T1[1, ] ~ PROPEL03c$C1E)
t.test(RNAseq_ALL03_T1[1, ] ~ PROPEL03c$C2A)
t.test(RNAseq_ALL03_T1[1, ] ~ PROPEL03c$C2C)
t.test(RNAseq_ALL03_T1[1, ] ~ PROPEL03c$C2D)
t.test(RNAseq_ALL03_T1[1, ] ~ PROPEL03c$C2E)
t.test(RNAseq_ALL03_T7[1, ] ~ PROPEL03c$C1A)
t.test(RNAseq_ALL03_T7[1, ] ~ PROPEL03c$C1C)
t.test(RNAseq_ALL03_T7[1, ] ~ PROPEL03c$C1D)
t.test(RNAseq_ALL03_T7[1, ] ~ PROPEL03c$C1E)
t.test(RNAseq_ALL03_T7[1, ] ~ PROPEL03c$C2A)
t.test(RNAseq_ALL03_T7[1, ] ~ PROPEL03c$C2C)
t.test(RNAseq_ALL03_T7[1, ] ~ PROPEL03c$C2D)
t.test(RNAseq_ALL03_T7[1, ] ~ PROPEL03c$C2E)
PV_T1_1A <- vector(); FC_T1_1A <- vector(); for (i in 1:nrow(RNAseq_ALL03_T1)){TTEST <- t.test(RNAseq_ALL03_T1[i, ] ~ PROPEL03c$C1A); PV_T1_1A[i] <- TTEST$p.value; FC_T1_1A[i] <- diff(TTEST$estimate)}
PV_T1_1C <- vector(); FC_T1_1C <- vector(); for (i in 1:nrow(RNAseq_ALL03_T1)){TTEST <- t.test(RNAseq_ALL03_T1[i, ] ~ PROPEL03c$C1C); PV_T1_1C[i] <- TTEST$p.value; FC_T1_1C[i] <- diff(TTEST$estimate)}
PV_T1_1D <- vector(); FC_T1_1D <- vector(); for (i in 1:nrow(RNAseq_ALL03_T1)){TTEST <- t.test(RNAseq_ALL03_T1[i, ] ~ PROPEL03c$C1D); PV_T1_1D[i] <- TTEST$p.value; FC_T1_1D[i] <- diff(TTEST$estimate)}
PV_T1_1E <- vector(); FC_T1_1E <- vector(); for (i in 1:nrow(RNAseq_ALL03_T1)){TTEST <- t.test(RNAseq_ALL03_T1[i, ] ~ PROPEL03c$C1E); PV_T1_1E[i] <- TTEST$p.value; FC_T1_1E[i] <- diff(TTEST$estimate)}
PV_T1_2A <- vector(); FC_T1_2A <- vector(); for (i in 1:nrow(RNAseq_ALL03_T1)){TTEST <- t.test(RNAseq_ALL03_T1[i, ] ~ PROPEL03c$C2A); PV_T1_2A[i] <- TTEST$p.value; FC_T1_2A[i] <- diff(TTEST$estimate)}
PV_T1_2C <- vector(); FC_T1_2C <- vector(); for (i in 1:nrow(RNAseq_ALL03_T1)){TTEST <- t.test(RNAseq_ALL03_T1[i, ] ~ PROPEL03c$C2C); PV_T1_2C[i] <- TTEST$p.value; FC_T1_2C[i] <- diff(TTEST$estimate)}
PV_T1_2D <- vector(); FC_T1_2D <- vector(); for (i in 1:nrow(RNAseq_ALL03_T1)){TTEST <- t.test(RNAseq_ALL03_T1[i, ] ~ PROPEL03c$C2D); PV_T1_2D[i] <- TTEST$p.value; FC_T1_2D[i] <- diff(TTEST$estimate)}
PV_T1_2E <- vector(); FC_T1_2E <- vector(); for (i in 1:nrow(RNAseq_ALL03_T1)){TTEST <- t.test(RNAseq_ALL03_T1[i, ] ~ PROPEL03c$C2E); PV_T1_2E[i] <- TTEST$p.value; FC_T1_2E[i] <- diff(TTEST$estimate)}
PV_T7_1A <- vector(); FC_T7_1A <- vector(); for (i in 1:nrow(RNAseq_ALL03_T7)){TTEST <- t.test(RNAseq_ALL03_T7[i, ] ~ PROPEL03c$C1A); PV_T7_1A[i] <- TTEST$p.value; FC_T7_1A[i] <- diff(TTEST$estimate)}
PV_T7_1C <- vector(); FC_T7_1C <- vector(); for (i in 1:nrow(RNAseq_ALL03_T7)){TTEST <- t.test(RNAseq_ALL03_T7[i, ] ~ PROPEL03c$C1C); PV_T7_1C[i] <- TTEST$p.value; FC_T7_1C[i] <- diff(TTEST$estimate)}
PV_T7_1D <- vector(); FC_T7_1D <- vector(); for (i in 1:nrow(RNAseq_ALL03_T7)){TTEST <- t.test(RNAseq_ALL03_T7[i, ] ~ PROPEL03c$C1D); PV_T7_1D[i] <- TTEST$p.value; FC_T7_1D[i] <- diff(TTEST$estimate)}
PV_T7_1E <- vector(); FC_T7_1E <- vector(); for (i in 1:nrow(RNAseq_ALL03_T7)){TTEST <- t.test(RNAseq_ALL03_T7[i, ] ~ PROPEL03c$C1E); PV_T7_1E[i] <- TTEST$p.value; FC_T7_1E[i] <- diff(TTEST$estimate)}
PV_T7_2A <- vector(); FC_T7_2A <- vector(); for (i in 1:nrow(RNAseq_ALL03_T7)){TTEST <- t.test(RNAseq_ALL03_T7[i, ] ~ PROPEL03c$C2A); PV_T7_2A[i] <- TTEST$p.value; FC_T7_2A[i] <- diff(TTEST$estimate)}
PV_T7_2C <- vector(); FC_T7_2C <- vector(); for (i in 1:nrow(RNAseq_ALL03_T7)){TTEST <- t.test(RNAseq_ALL03_T7[i, ] ~ PROPEL03c$C2C); PV_T7_2C[i] <- TTEST$p.value; FC_T7_2C[i] <- diff(TTEST$estimate)}
PV_T7_2D <- vector(); FC_T7_2D <- vector(); for (i in 1:nrow(RNAseq_ALL03_T7)){TTEST <- t.test(RNAseq_ALL03_T7[i, ] ~ PROPEL03c$C2D); PV_T7_2D[i] <- TTEST$p.value; FC_T7_2D[i] <- diff(TTEST$estimate)}
PV_T7_2E <- vector(); FC_T7_2E <- vector(); for (i in 1:nrow(RNAseq_ALL03_T7)){TTEST <- t.test(RNAseq_ALL03_T7[i, ] ~ PROPEL03c$C2E); PV_T7_2E[i] <- TTEST$p.value; FC_T7_2E[i] <- diff(TTEST$estimate)}
PVFC_DEG02 <- data.frame(PV_T1_1A,FC_T1_1A,PV_T1_1C,FC_T1_1C,PV_T1_1D,FC_T1_1D,PV_T1_1E,FC_T1_1E,PV_T1_2A,FC_T1_2A,PV_T1_2C,FC_T1_2C,PV_T1_2D,FC_T1_2D,PV_T1_2E,FC_T1_2E,PV_T7_1A,FC_T7_1A,PV_T7_1C,FC_T7_1C,PV_T7_1D,FC_T7_1D,PV_T7_1E,FC_T7_1E,PV_T7_2A,FC_T7_2A,PV_T7_2C,FC_T7_2C,PV_T7_2D,FC_T7_2D,PV_T7_2E,FC_T7_2E)
rownames(PVFC_DEG02) <- rownames(RNAseq_ALL03_T1)
write.table(PVFC_DEG02, file="PVFC_DEG02.txt",sep="\t",quote=FALSE,col.names=NA)

#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#?ó¨Í∏∞ÏÑúÎ∂Ä?Ñ∞ ?ãú?ûë

PVFC_DEG02 <- as.data.frame(fread("PVFC_DEG02.txt.gz", sep = "\t")); rownames(PVFC_DEG02) <- PVFC_DEG02$V1; PVFC_DEG02 <- PVFC_DEG02[,-1]
PVFC_DEG02[is.na(PVFC_DEG02)] <- 1
head(PVFC_DEG02)

#### A1?ù¥ Í≥ÑÏÜç ?Ç¨?ùºÏßÄ?äî ?ò§Î•? ?ïÑ?ûò ÏΩîÎìú
PVFC_DEG03 <- as.data.frame(fread("PVFC_DEG03.txt.gz", sep = "\t")); rownames(PVFC_DEG03) <- PVFC_DEG03$V1; PVFC_DEG03 <- PVFC_DEG03[,-1]

PVFC_DEG03[is.na(PVFC_DEG03)] <- 1
head(PVFC_DEG03)


table(rownames(PVFC_DEG02)==rownames(PVFC_DEG03))
PVFC_DEG <- cbind(PVFC_DEG02,PVFC_DEG03)
head(PVFC_DEG)

summary(PVFC_DEG$FC_E1)



############################################################################################################
#Gene number range: 200-220
#?ó∞Íµ? Î™©Ìëú 3 Í∑∏Î£π ?ãπ T1 vs T7?ùò ?ú†?†Ñ?ûê Î∞úÌòÑ?üâ Ï∞®Ïù¥ Î∂ÑÏÑù Í≥ÑÏÜç
#?ú†?ùò?ïú Î≥Ä?ôî ?ûà?äî ?ú†?†Ñ?ûê Ï∂îÏ∂ú, P-value <0.05, FC
LIST_E1 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E1 < 0.05 & abs(FC_E1) > 0.145), PV_E1,FC_E1); LIST_E1$ENSEMBL <- rownames(LIST_E1); names(LIST_E1)[1:2] <- c("PV","FC"); dim(LIST_E1) #221,3
LIST_N1 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N1 < 0.05 & abs(FC_N1) > 0.040), PV_N1,FC_N1); LIST_N1$ENSEMBL <- rownames(LIST_N1); names(LIST_N1)[1:2] <- c("PV","FC"); dim(LIST_N1) #207,3
LIST_E3 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E3 < 0.05 & abs(FC_E3) > 0.070), PV_E3,FC_E3); LIST_E3$ENSEMBL <- rownames(LIST_E3); names(LIST_E3)[1:2] <- c("PV","FC"); dim(LIST_E3) #202,3
LIST_N3 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N3 < 0.05 & abs(FC_N3) > 0.110), PV_N3,FC_N3); LIST_N3$ENSEMBL <- rownames(LIST_N3); names(LIST_N3)[1:2] <- c("PV","FC"); dim(LIST_N3) #204,3
LIST_E4 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E4 < 0.05 & abs(FC_E4) > 0.153), PV_E4,FC_E4); LIST_E4$ENSEMBL <- rownames(LIST_E4); names(LIST_E4)[1:2] <- c("PV","FC"); dim(LIST_E4) #211,3
LIST_N4 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N4 < 0.05 & abs(FC_N4) > 0.035), PV_N4,FC_N4); LIST_N4$ENSEMBL <- rownames(LIST_N4); names(LIST_N4)[1:2] <- c("PV","FC"); dim(LIST_N4) #211,3
LIST_E5 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E5 < 0.05 & abs(FC_E5) > 0.113), PV_E5,FC_E5); LIST_E5$ENSEMBL <- rownames(LIST_E5); names(LIST_E5)[1:2] <- c("PV","FC"); dim(LIST_E5) #213,3
LIST_N5 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N5 < 0.05 & abs(FC_N5) > 0.045), PV_N5,FC_N5); LIST_N5$ENSEMBL <- rownames(LIST_N5); names(LIST_N5)[1:2] <- c("PV","FC"); dim(LIST_N5) #206,3
LIST_E6 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E6 < 0.05 & abs(FC_E6) > 0.131), PV_E6,FC_E6); LIST_E6$ENSEMBL <- rownames(LIST_E6); names(LIST_E6)[1:2] <- c("PV","FC"); dim(LIST_E6) #217,3
LIST_N6 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N6 < 0.05 & abs(FC_N6) > 0.035), PV_N6,FC_N6); LIST_N6$ENSEMBL <- rownames(LIST_N6); names(LIST_N6)[1:2] <- c("PV","FC"); dim(LIST_N6) #215,3
LIST_E8 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E8 < 0.05 & abs(FC_E8) > 0.126), PV_E8,FC_E8); LIST_E8$ENSEMBL <- rownames(LIST_E8); names(LIST_E8)[1:2] <- c("PV","FC"); dim(LIST_E8) #215,3
LIST_N8 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N8 < 0.05 & abs(FC_N8) > 0.052), PV_N8,FC_N8); LIST_N8$ENSEMBL <- rownames(LIST_N8); names(LIST_N8)[1:2] <- c("PV","FC"); dim(LIST_N8) #212,3
LIST_E9 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E9 < 0.05 & abs(FC_E9) > 0.140), PV_E9,FC_E9); LIST_E9$ENSEMBL <- rownames(LIST_E9); names(LIST_E9)[1:2] <- c("PV","FC"); dim(LIST_E9) #215,3
LIST_N9 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N9 < 0.05 & abs(FC_N9) > 0.035), PV_N9,FC_N9); LIST_N9$ENSEMBL <- rownames(LIST_N9); names(LIST_N9)[1:2] <- c("PV","FC"); dim(LIST_N9) #218,3
LIST_E0 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E0 < 0.05 & abs(FC_E0) > 0.122), PV_E0,FC_E0); LIST_E0$ENSEMBL <- rownames(LIST_E0); names(LIST_E0)[1:2] <- c("PV","FC"); dim(LIST_E0) #214,3
LIST_N0 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N0 < 0.05 & abs(FC_N0) > 0.045), PV_N0,FC_N0); LIST_N0$ENSEMBL <- rownames(LIST_N0); names(LIST_N0)[1:2] <- c("PV","FC"); dim(LIST_N0) #207,3
LIST_A1 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_A1 < 0.05 & abs(FC_A1) > 0.045), PV_A1,FC_A1); LIST_A1$ENSEMBL <- rownames(LIST_A1); names(LIST_A1)[1:2] <- c("PV","FC"); dim(LIST_A1) #1310,3

#?úÑ?óê ?ûà?äî ?ú†?†Ñ?ûê Î¶¨Ïä§?ä∏?ì§?ùÑ Î™®Îëê ?ï©ÏπòÎäî ?ûë?óÖ
LIST_EN <- Reduce(union, list(rownames(LIST_E1),rownames(LIST_N1),rownames(LIST_E3),rownames(LIST_N3),rownames(LIST_E4),rownames(LIST_N4),rownames(LIST_E5),rownames(LIST_N5),rownames(LIST_E6),rownames(LIST_N6),rownames(LIST_E8),rownames(LIST_N8),rownames(LIST_E9),rownames(LIST_N9),rownames(LIST_E0),rownames(LIST_N0), rownames(LIST_A1)))
RNAseq_ALL04 <- RNAseq_ALL03[which(rownames(RNAseq_ALL03) %in% LIST_EN),]
RNAseq_ALL04$ENSEMBL <- rownames(RNAseq_ALL04)
RNAseq_ALL04 <- merge(x=RNAseq_ALL04,y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL04) <- RNAseq_ALL04$ENSEMBL; RNAseq_ALL04 <- RNAseq_ALL04[,-1]
colnames(RNAseq_ALL04)
RNAseq_ALL03_E1x <- as.data.frame(RNAseq_ALL03_E1); RNAseq_ALL03_E1x$ENSEMBL <- rownames(RNAseq_ALL03_E1x); RNAseq_ALL03_E1x <- merge(x=RNAseq_ALL03_E1x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_N1x <- as.data.frame(RNAseq_ALL03_N1); RNAseq_ALL03_N1x$ENSEMBL <- rownames(RNAseq_ALL03_N1x); RNAseq_ALL03_N1x <- merge(x=RNAseq_ALL03_N1x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_E3x <- as.data.frame(RNAseq_ALL03_E3); RNAseq_ALL03_E3x$ENSEMBL <- rownames(RNAseq_ALL03_E3x); RNAseq_ALL03_E3x <- merge(x=RNAseq_ALL03_E3x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_N3x <- as.data.frame(RNAseq_ALL03_N3); RNAseq_ALL03_N3x$ENSEMBL <- rownames(RNAseq_ALL03_N3x); RNAseq_ALL03_N3x <- merge(x=RNAseq_ALL03_N3x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_E4x <- as.data.frame(RNAseq_ALL03_E4); RNAseq_ALL03_E4x$ENSEMBL <- rownames(RNAseq_ALL03_E4x); RNAseq_ALL03_E4x <- merge(x=RNAseq_ALL03_E4x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_N4x <- as.data.frame(RNAseq_ALL03_N4); RNAseq_ALL03_N4x$ENSEMBL <- rownames(RNAseq_ALL03_N4x); RNAseq_ALL03_N4x <- merge(x=RNAseq_ALL03_N4x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_E5x <- as.data.frame(RNAseq_ALL03_E5); RNAseq_ALL03_E5x$ENSEMBL <- rownames(RNAseq_ALL03_E5x); RNAseq_ALL03_E5x <- merge(x=RNAseq_ALL03_E5x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_N5x <- as.data.frame(RNAseq_ALL03_N5); RNAseq_ALL03_N5x$ENSEMBL <- rownames(RNAseq_ALL03_N5x); RNAseq_ALL03_N5x <- merge(x=RNAseq_ALL03_N5x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_E6x <- as.data.frame(RNAseq_ALL03_E6); RNAseq_ALL03_E6x$ENSEMBL <- rownames(RNAseq_ALL03_E6x); RNAseq_ALL03_E6x <- merge(x=RNAseq_ALL03_E6x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_N6x <- as.data.frame(RNAseq_ALL03_N6); RNAseq_ALL03_N6x$ENSEMBL <- rownames(RNAseq_ALL03_N6x); RNAseq_ALL03_N6x <- merge(x=RNAseq_ALL03_N6x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_E8x <- as.data.frame(RNAseq_ALL03_E8); RNAseq_ALL03_E8x$ENSEMBL <- rownames(RNAseq_ALL03_E8x); RNAseq_ALL03_E8x <- merge(x=RNAseq_ALL03_E8x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_N8x <- as.data.frame(RNAseq_ALL03_N8); RNAseq_ALL03_N8x$ENSEMBL <- rownames(RNAseq_ALL03_N8x); RNAseq_ALL03_N8x <- merge(x=RNAseq_ALL03_N8x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_E9x <- as.data.frame(RNAseq_ALL03_E9); RNAseq_ALL03_E9x$ENSEMBL <- rownames(RNAseq_ALL03_E9x); RNAseq_ALL03_E9x <- merge(x=RNAseq_ALL03_E9x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_N9x <- as.data.frame(RNAseq_ALL03_N9); RNAseq_ALL03_N9x$ENSEMBL <- rownames(RNAseq_ALL03_N9x); RNAseq_ALL03_N9x <- merge(x=RNAseq_ALL03_N9x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_E0x <- as.data.frame(RNAseq_ALL03_E0); RNAseq_ALL03_E0x$ENSEMBL <- rownames(RNAseq_ALL03_E0x); RNAseq_ALL03_E0x <- merge(x=RNAseq_ALL03_E0x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_N0x <- as.data.frame(RNAseq_ALL03_N0); RNAseq_ALL03_N0x$ENSEMBL <- rownames(RNAseq_ALL03_N0x); RNAseq_ALL03_N0x <- merge(x=RNAseq_ALL03_N0x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_A1x <- as.data.frame(RNAseq_ALL03_A1); RNAseq_ALL03_A1x$ENSEMBL <- rownames(RNAseq_ALL03_A1x); RNAseq_ALL03_A1x <- merge(x=RNAseq_ALL03_A1x, y=RNAseq_GENEs, by="ENSEMBL")
RNAseq_ALL03_E1x <- merge(x=RNAseq_ALL03_E1x, y=LIST_E1, by="ENSEMBL"); rownames(RNAseq_ALL03_E1x) <- RNAseq_ALL03_E1x$ENSEMBL; RNAseq_ALL03_E1x <- RNAseq_ALL03_E1x[,-1]
RNAseq_ALL03_N1x <- merge(x=RNAseq_ALL03_N1x, y=LIST_N1, by="ENSEMBL"); rownames(RNAseq_ALL03_N1x) <- RNAseq_ALL03_N1x$ENSEMBL; RNAseq_ALL03_N1x <- RNAseq_ALL03_N1x[,-1]
RNAseq_ALL03_E3x <- merge(x=RNAseq_ALL03_E3x, y=LIST_E3, by="ENSEMBL"); rownames(RNAseq_ALL03_E3x) <- RNAseq_ALL03_E3x$ENSEMBL; RNAseq_ALL03_E3x <- RNAseq_ALL03_E3x[,-1]
RNAseq_ALL03_N3x <- merge(x=RNAseq_ALL03_N3x, y=LIST_N3, by="ENSEMBL"); rownames(RNAseq_ALL03_N3x) <- RNAseq_ALL03_N3x$ENSEMBL; RNAseq_ALL03_N3x <- RNAseq_ALL03_N3x[,-1]
RNAseq_ALL03_E4x <- merge(x=RNAseq_ALL03_E4x, y=LIST_E4, by="ENSEMBL"); rownames(RNAseq_ALL03_E4x) <- RNAseq_ALL03_E4x$ENSEMBL; RNAseq_ALL03_E4x <- RNAseq_ALL03_E4x[,-1]
RNAseq_ALL03_N4x <- merge(x=RNAseq_ALL03_N4x, y=LIST_N4, by="ENSEMBL"); rownames(RNAseq_ALL03_N4x) <- RNAseq_ALL03_N4x$ENSEMBL; RNAseq_ALL03_N4x <- RNAseq_ALL03_N4x[,-1]
RNAseq_ALL03_E5x <- merge(x=RNAseq_ALL03_E5x, y=LIST_E5, by="ENSEMBL"); rownames(RNAseq_ALL03_E5x) <- RNAseq_ALL03_E5x$ENSEMBL; RNAseq_ALL03_E5x <- RNAseq_ALL03_E5x[,-1]
RNAseq_ALL03_N5x <- merge(x=RNAseq_ALL03_N5x, y=LIST_N5, by="ENSEMBL"); rownames(RNAseq_ALL03_N5x) <- RNAseq_ALL03_N5x$ENSEMBL; RNAseq_ALL03_N5x <- RNAseq_ALL03_N5x[,-1]
RNAseq_ALL03_E6x <- merge(x=RNAseq_ALL03_E6x, y=LIST_E6, by="ENSEMBL"); rownames(RNAseq_ALL03_E6x) <- RNAseq_ALL03_E6x$ENSEMBL; RNAseq_ALL03_E6x <- RNAseq_ALL03_E6x[,-1]
RNAseq_ALL03_N6x <- merge(x=RNAseq_ALL03_N6x, y=LIST_N6, by="ENSEMBL"); rownames(RNAseq_ALL03_N6x) <- RNAseq_ALL03_N6x$ENSEMBL; RNAseq_ALL03_N6x <- RNAseq_ALL03_N6x[,-1]
RNAseq_ALL03_E8x <- merge(x=RNAseq_ALL03_E8x, y=LIST_E8, by="ENSEMBL"); rownames(RNAseq_ALL03_E8x) <- RNAseq_ALL03_E8x$ENSEMBL; RNAseq_ALL03_E8x <- RNAseq_ALL03_E8x[,-1]
RNAseq_ALL03_N8x <- merge(x=RNAseq_ALL03_N8x, y=LIST_N8, by="ENSEMBL"); rownames(RNAseq_ALL03_N8x) <- RNAseq_ALL03_N8x$ENSEMBL; RNAseq_ALL03_N8x <- RNAseq_ALL03_N8x[,-1]
RNAseq_ALL03_E9x <- merge(x=RNAseq_ALL03_E9x, y=LIST_E9, by="ENSEMBL"); rownames(RNAseq_ALL03_E9x) <- RNAseq_ALL03_E9x$ENSEMBL; RNAseq_ALL03_E9x <- RNAseq_ALL03_E9x[,-1]
RNAseq_ALL03_N9x <- merge(x=RNAseq_ALL03_N9x, y=LIST_N9, by="ENSEMBL"); rownames(RNAseq_ALL03_N9x) <- RNAseq_ALL03_N9x$ENSEMBL; RNAseq_ALL03_N9x <- RNAseq_ALL03_N9x[,-1]
RNAseq_ALL03_E0x <- merge(x=RNAseq_ALL03_E0x, y=LIST_E0, by="ENSEMBL"); rownames(RNAseq_ALL03_E0x) <- RNAseq_ALL03_E0x$ENSEMBL; RNAseq_ALL03_E0x <- RNAseq_ALL03_E0x[,-1]
RNAseq_ALL03_N0x <- merge(x=RNAseq_ALL03_N0x, y=LIST_N0, by="ENSEMBL"); rownames(RNAseq_ALL03_N0x) <- RNAseq_ALL03_N0x$ENSEMBL; RNAseq_ALL03_N0x <- RNAseq_ALL03_N0x[,-1]
RNAseq_ALL03_A1x <- merge(x=RNAseq_ALL03_A1x, y=LIST_A1, by="ENSEMBL"); rownames(RNAseq_ALL03_A1x) <- RNAseq_ALL03_A1x$ENSEMBL; RNAseq_ALL03_A1x <- RNAseq_ALL03_A1x[,-1]


#æÍµµ ≤œ ø¿∑°∞…∏≤∏≤
dbs <- listEnrichrDbs()
ERR_E1 <- enrichr(RNAseq_ALL03_E1x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_N1 <- enrichr(RNAseq_ALL03_N1x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_E3 <- enrichr(RNAseq_ALL03_E3x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_N3 <- enrichr(RNAseq_ALL03_N3x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_E4 <- enrichr(RNAseq_ALL03_E4x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_N4 <- enrichr(RNAseq_ALL03_N4x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_E5 <- enrichr(RNAseq_ALL03_E5x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_N5 <- enrichr(RNAseq_ALL03_N5x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_E6 <- enrichr(RNAseq_ALL03_E6x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_N6 <- enrichr(RNAseq_ALL03_N6x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_E8 <- enrichr(RNAseq_ALL03_E8x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_N8 <- enrichr(RNAseq_ALL03_N8x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_E9 <- enrichr(RNAseq_ALL03_E9x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_N9 <- enrichr(RNAseq_ALL03_N9x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_E0 <- enrichr(RNAseq_ALL03_E0x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_N0 <- enrichr(RNAseq_ALL03_N0x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR_A1 <- enrichr(RNAseq_ALL03_A1x$SYMBOL, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
#E1[[1]]??Ä GOÎ•? E1[[2]]?äî KEGG Í≤∞Í≥ºÎ•? ?ùòÎØ∏Ìï®?ï®
ERR_E1_GO <- ERR_E1[[1]]; ERR_E1_KG <- ERR_E1[[2]]
ERR_N1_GO <- ERR_N1[[1]]; ERR_N1_KG <- ERR_N1[[2]]
ERR_E3_GO <- ERR_E3[[1]]; ERR_E3_KG <- ERR_E3[[2]]
ERR_N3_GO <- ERR_N3[[1]]; ERR_N3_KG <- ERR_N3[[2]]
ERR_E4_GO <- ERR_E4[[1]]; ERR_E4_KG <- ERR_E4[[2]]
ERR_N4_GO <- ERR_N4[[1]]; ERR_N4_KG <- ERR_N4[[2]]
ERR_E5_GO <- ERR_E5[[1]]; ERR_E5_KG <- ERR_E5[[2]]
ERR_N5_GO <- ERR_N5[[1]]; ERR_N5_KG <- ERR_N5[[2]]
ERR_E6_GO <- ERR_E6[[1]]; ERR_E6_KG <- ERR_E6[[2]]
ERR_N6_GO <- ERR_N6[[1]]; ERR_N6_KG <- ERR_N6[[2]]
ERR_E8_GO <- ERR_E8[[1]]; ERR_E8_KG <- ERR_E8[[2]]
ERR_N8_GO <- ERR_N8[[1]]; ERR_N8_KG <- ERR_N8[[2]]
ERR_E9_GO <- ERR_E9[[1]]; ERR_E9_KG <- ERR_E9[[2]]
ERR_N9_GO <- ERR_N9[[1]]; ERR_N9_KG <- ERR_N9[[2]]
ERR_E0_GO <- ERR_E0[[1]]; ERR_E0_KG <- ERR_E0[[2]]
ERR_N0_GO <- ERR_N0[[1]]; ERR_N0_KG <- ERR_N0[[2]]
ERR_A1_GO <- ERR_A1[[1]]; ERR_A1_KG <- ERR_A1[[2]]
#GO Ï§? p.value ÎØ∏Îßå?ù∏ GO term?ùÑ ?ôï?ù∏
#Í≤∞Í≥º?óê?Ñú [1] 20 (GO term), 4 (colum)
ERR_E1_GOx <- dplyr::filter(ERR_E1_GO, P.value < 0.001)[,c(1,3,9)]; ERR_E1_GOx$CLASS <- "GO"; dim(ERR_E1_GOx) #20,4
ERR_N1_GOx <- dplyr::filter(ERR_N1_GO, P.value < 0.050)[,c(1,3,9)]; ERR_N1_GOx$CLASS <- "GO"; dim(ERR_N1_GOx) #12,4
ERR_E3_GOx <- dplyr::filter(ERR_E3_GO, P.value < 0.050)[,c(1,3,9)]; ERR_E3_GOx$CLASS <- "GO"; dim(ERR_E3_GOx) #16,4
ERR_N3_GOx <- dplyr::filter(ERR_N3_GO, P.value < 0.001)[,c(1,3,9)]; ERR_N3_GOx$CLASS <- "GO"; dim(ERR_N3_GOx) #23,4

ERR_E4_GOx <- dplyr::filter(ERR_E4_GO, P.value < 0.001)[,c(1,3,9)]; ERR_E4_GOx$CLASS <- "GO"; dim(ERR_E4_GOx) #12,4
ERR_N4_GOx <- dplyr::filter(ERR_N4_GO, P.value < 0.050)[,c(1,3,9)]; ERR_N4_GOx$CLASS <- "GO"; dim(ERR_N4_GOx) # 6,4
ERR_E5_GOx <- dplyr::filter(ERR_E5_GO, P.value < 0.001)[,c(1,3,9)]; ERR_E5_GOx$CLASS <- "GO"; dim(ERR_E5_GOx) #13,4
ERR_N5_GOx <- dplyr::filter(ERR_N5_GO, P.value < 0.050)[,c(1,3,9)]; ERR_N5_GOx$CLASS <- "GO"; dim(ERR_N5_GOx) #10,4

ERR_E6_GOx <- dplyr::filter(ERR_E6_GO, P.value < 0.001)[,c(1,3,9)]; ERR_E6_GOx$CLASS <- "GO"; dim(ERR_E6_GOx) #
ERR_N6_GOx <- dplyr::filter(ERR_N6_GO, P.value < 0.050)[,c(1,3,9)]; ERR_N6_GOx$CLASS <- "GO"; dim(ERR_N6_GOx) #
ERR_E8_GOx <- dplyr::filter(ERR_E8_GO, P.value < 0.001)[,c(1,3,9)]; ERR_E8_GOx$CLASS <- "GO"; dim(ERR_E8_GOx) #
ERR_N8_GOx <- dplyr::filter(ERR_N8_GO, P.value < 0.005)[,c(1,3,9)]; ERR_N8_GOx$CLASS <- "GO"; dim(ERR_N8_GOx) #

ERR_E9_GOx <- dplyr::filter(ERR_E9_GO, P.value < 0.001)[,c(1,3,9)]; ERR_E9_GOx$CLASS <- "GO"; dim(ERR_E9_GOx)
ERR_N9_GOx <- dplyr::filter(ERR_N9_GO, P.value < 0.050)[,c(1,3,9)]; ERR_N9_GOx$CLASS <- "GO"; dim(ERR_N9_GOx)
ERR_E0_GOx <- dplyr::filter(ERR_E0_GO, P.value < 0.050)[,c(1,3,9)]; ERR_E0_GOx$CLASS <- "GO"; dim(ERR_E0_GOx)
ERR_N0_GOx <- dplyr::filter(ERR_N0_GO, P.value < 0.001)[,c(1,3,9)]; ERR_N0_GOx$CLASS <- "GO"; dim(ERR_N0_GOx)
ERR_A1_GOx <- dplyr::filter(ERR_A1_GO, P.value < 0.001)[,c(1,3,9)]; ERR_A1_GOx$CLASS <- "GO"; dim(ERR_A1_GOx)

ERR_E1_KGx <- dplyr::filter(ERR_E1_KG, P.value < 0.050)[,c(1,3,9)]; ERR_E1_KGx$CLASS <- "KG"; dim(ERR_E1_KGx)
ERR_N1_KGx <- dplyr::filter(ERR_N1_KG, P.value < 0.050)[,c(1,3,9)]; ERR_N1_KGx$CLASS <- "KG"; dim(ERR_N1_KGx)
ERR_E3_KGx <- dplyr::filter(ERR_E3_KG, P.value < 0.050)[,c(1,3,9)]; ERR_E3_KGx$CLASS <- "KG"; dim(ERR_E3_KGx)
ERR_N3_KGx <- dplyr::filter(ERR_N3_KG, P.value < 0.050)[,c(1,3,9)]; ERR_N3_KGx$CLASS <- "KG"; dim(ERR_N3_KGx)

ERR_E4_KGx <- dplyr::filter(ERR_E4_KG, P.value < 0.050)[,c(1,3,9)]; ERR_E4_KGx$CLASS <- "KG"; dim(ERR_E4_KGx)
#N4?äî ?ú†?ùò?ïú KGÍ∞Ä ?óÜ?ùå.
dim(ERR_N4_KGx)
ERR_N4_KGx <- dplyr::filter(ERR_N4_KG, P.value < 0.050)[,c(1,3,9)]; ERR_N4_KGx$CLASS <- "KG"; dim(ERR_N4_KGx)
ERR_E5_KGx <- dplyr::filter(ERR_E5_KG, P.value < 0.050)[,c(1,3,9)]; ERR_E5_KGx$CLASS <- "KG"; dim(ERR_E5_KGx)
ERR_N5_KGx <- dplyr::filter(ERR_N5_KG, P.value < 0.050)[,c(1,3,9)]; ERR_N5_KGx$CLASS <- "KG"; dim(ERR_N5_KGx)

ERR_E6_KGx <- dplyr::filter(ERR_E6_KG, P.value < 0.050)[,c(1,3,9)]; ERR_E6_KGx$CLASS <- "KG"; dim(ERR_E6_KGx)
ERR_N6_KGx <- dplyr::filter(ERR_N6_KG, P.value < 0.050)[,c(1,3,9)]; ERR_N6_KGx$CLASS <- "KG"; dim(ERR_N6_KGx)
ERR_E8_KGx <- dplyr::filter(ERR_E8_KG, P.value < 0.050)[,c(1,3,9)]; ERR_E8_KGx$CLASS <- "KG"; dim(ERR_E8_KGx)
ERR_N8_KGx <- dplyr::filter(ERR_N8_KG, P.value < 0.050)[,c(1,3,9)]; ERR_N8_KGx$CLASS <- "KG"; dim(ERR_N8_KGx)

ERR_E9_KGx <- dplyr::filter(ERR_E9_KG, P.value < 0.050)[,c(1,3,9)]; ERR_E9_KGx$CLASS <- "KG"; dim(ERR_E9_KGx)
#N9?äî ?ú†?ùò?ïú KGÍ∞Ä ?óÜ?ùå.
dim(ERR_N9_KGx)
ERR_N9_KGx <- dplyr::filter(ERR_N9_KG, P.value < 0.050)[,c(1,3,9)]; ERR_N9_KGx$CLASS <- "KG"; dim(ERR_N9_KGx)
ERR_E0_KGx <- dplyr::filter(ERR_E0_KG, P.value < 0.050)[,c(1,3,9)]; ERR_E0_KGx$CLASS <- "KG"; dim(ERR_E0_KGx)
ERR_N0_KGx <- dplyr::filter(ERR_N0_KG, P.value < 0.050)[,c(1,3,9)]; ERR_N0_KGx$CLASS <- "KG"; dim(ERR_N0_KGx)
ERR_A1_KGx <- dplyr::filter(ERR_A1_KG, P.value < 0.050)[,c(1,3,9)]; ERR_A1_KGx$CLASS <- "KG"; dim(ERR_A1_KGx)

write.table(rbind(ERR_E1_GOx,ERR_E1_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_N1_GOx,ERR_N1_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_E3_GOx,ERR_E3_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_N3_GOx,ERR_N3_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_E4_GOx,ERR_E4_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_N4_GOx,ERR_N4_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_E5_GOx,ERR_E5_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_N5_GOx,ERR_N5_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_E6_GOx,ERR_E6_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_N6_GOx,ERR_N6_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_E8_GOx,ERR_E8_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_N8_GOx,ERR_N8_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_E9_GOx,ERR_E9_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_N9_GOx,ERR_N9_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_E0_GOx,ERR_E0_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_N0_GOx,ERR_N0_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(rbind(ERR_A1_GOx,ERR_A1_KGx), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)

#ENSEMBL IDÎ•? ENTREZ IDÎ°? Î≥Ä?ôò
RNAseq_ALL03_E1x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_E1x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_N1x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_N1x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_E3x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_E3x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_N3x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_N3x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_E4x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_E4x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_N4x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_N4x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_E5x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_E5x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_N5x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_N5x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_E6x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_E6x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_N6x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_N6x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_E8x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_E8x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_N8x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_N8x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_E9x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_E9x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_N9x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_N9x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_E0x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_E0x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_N0x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_N0x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")
RNAseq_ALL03_A1x$ENTREZID <- mapIds(org.Hs.eg.db, keys = substr(rownames(RNAseq_ALL03_A1x),1,15),column = "ENTREZID",keytype = "ENSEMBL",multiVals = "first")

#NA Í∞íÏùÑ Î™®Îëê ÏßÄ?õå?Ñú ENTREZID?óê Îß§Ïπ≠?êòÏßÄ ?ïä??Ä ?ú†?†Ñ?ûêÎ•? ?Ç≠?†ú
RNAseq_ALL03_E1x <- dplyr::filter(RNAseq_ALL03_E1x, !is.na(ENTREZID))
RNAseq_ALL03_N1x <- dplyr::filter(RNAseq_ALL03_N1x, !is.na(ENTREZID))
RNAseq_ALL03_E3x <- dplyr::filter(RNAseq_ALL03_E3x, !is.na(ENTREZID))
RNAseq_ALL03_N3x <- dplyr::filter(RNAseq_ALL03_N3x, !is.na(ENTREZID))
RNAseq_ALL03_E4x <- dplyr::filter(RNAseq_ALL03_E4x, !is.na(ENTREZID))
RNAseq_ALL03_N4x <- dplyr::filter(RNAseq_ALL03_N4x, !is.na(ENTREZID))
RNAseq_ALL03_E5x <- dplyr::filter(RNAseq_ALL03_E5x, !is.na(ENTREZID))
RNAseq_ALL03_N5x <- dplyr::filter(RNAseq_ALL03_N5x, !is.na(ENTREZID))
RNAseq_ALL03_E6x <- dplyr::filter(RNAseq_ALL03_E6x, !is.na(ENTREZID))
RNAseq_ALL03_N6x <- dplyr::filter(RNAseq_ALL03_N6x, !is.na(ENTREZID))
RNAseq_ALL03_E8x <- dplyr::filter(RNAseq_ALL03_E8x, !is.na(ENTREZID))
RNAseq_ALL03_N8x <- dplyr::filter(RNAseq_ALL03_N8x, !is.na(ENTREZID))
RNAseq_ALL03_E9x <- dplyr::filter(RNAseq_ALL03_E9x, !is.na(ENTREZID))
RNAseq_ALL03_N9x <- dplyr::filter(RNAseq_ALL03_N9x, !is.na(ENTREZID))
RNAseq_ALL03_E0x <- dplyr::filter(RNAseq_ALL03_E0x, !is.na(ENTREZID))
RNAseq_ALL03_N0x <- dplyr::filter(RNAseq_ALL03_N0x, !is.na(ENTREZID))
RNAseq_ALL03_A1x <- dplyr::filter(RNAseq_ALL03_A1x, !is.na(ENTREZID))

#ENTREZ IDÎ•? keyÎ°? ?ïòÍ≥? FC Í∞íÏùÑ valueÎ°? Í∞ÄÏßÄ?äî ?òï?ÉúÎ°? Î≥ÄÍ≤?
RNAseq_ALL03_E1y <- as.vector(RNAseq_ALL03_E1x$FC); names(RNAseq_ALL03_E1y) <- RNAseq_ALL03_E1x$ENTREZID
RNAseq_ALL03_N1y <- as.vector(RNAseq_ALL03_N1x$FC); names(RNAseq_ALL03_N1y) <- RNAseq_ALL03_N1x$ENTREZID
RNAseq_ALL03_E3y <- as.vector(RNAseq_ALL03_E3x$FC); names(RNAseq_ALL03_E3y) <- RNAseq_ALL03_E3x$ENTREZID
RNAseq_ALL03_N3y <- as.vector(RNAseq_ALL03_N3x$FC); names(RNAseq_ALL03_N3y) <- RNAseq_ALL03_N3x$ENTREZID
RNAseq_ALL03_E4y <- as.vector(RNAseq_ALL03_E4x$FC); names(RNAseq_ALL03_E4y) <- RNAseq_ALL03_E4x$ENTREZID
RNAseq_ALL03_N4y <- as.vector(RNAseq_ALL03_N4x$FC); names(RNAseq_ALL03_N4y) <- RNAseq_ALL03_N4x$ENTREZID
RNAseq_ALL03_E5y <- as.vector(RNAseq_ALL03_E5x$FC); names(RNAseq_ALL03_E5y) <- RNAseq_ALL03_E5x$ENTREZID
RNAseq_ALL03_N5y <- as.vector(RNAseq_ALL03_N5x$FC); names(RNAseq_ALL03_N5y) <- RNAseq_ALL03_N5x$ENTREZID
RNAseq_ALL03_E6y <- as.vector(RNAseq_ALL03_E6x$FC); names(RNAseq_ALL03_E6y) <- RNAseq_ALL03_E6x$ENTREZID
RNAseq_ALL03_N6y <- as.vector(RNAseq_ALL03_N6x$FC); names(RNAseq_ALL03_N6y) <- RNAseq_ALL03_N6x$ENTREZID
RNAseq_ALL03_E8y <- as.vector(RNAseq_ALL03_E8x$FC); names(RNAseq_ALL03_E8y) <- RNAseq_ALL03_E8x$ENTREZID
RNAseq_ALL03_N8y <- as.vector(RNAseq_ALL03_N8x$FC); names(RNAseq_ALL03_N8y) <- RNAseq_ALL03_N8x$ENTREZID
RNAseq_ALL03_E9y <- as.vector(RNAseq_ALL03_E9x$FC); names(RNAseq_ALL03_E9y) <- RNAseq_ALL03_E9x$ENTREZID
RNAseq_ALL03_N9y <- as.vector(RNAseq_ALL03_N9x$FC); names(RNAseq_ALL03_N9y) <- RNAseq_ALL03_N9x$ENTREZID
RNAseq_ALL03_E0y <- as.vector(RNAseq_ALL03_E0x$FC); names(RNAseq_ALL03_E0y) <- RNAseq_ALL03_E0x$ENTREZID
RNAseq_ALL03_N0y <- as.vector(RNAseq_ALL03_N0x$FC); names(RNAseq_ALL03_N0y) <- RNAseq_ALL03_N0x$ENTREZID
RNAseq_ALL03_A1y <- as.vector(RNAseq_ALL03_A1x$FC); names(RNAseq_ALL03_A1y) <- RNAseq_ALL03_A1x$ENTREZID

#enrichDNG ?ï®?àò?äî ENTREZ ID Î¶¨Ïä§?ä∏Î•? ?ûÖ?†•?ïò?ó¨ ÏßàÎ≥ë Í¥Ä?†® ?ú†?†Ñ?ûê Î∂ÑÏÑù?ùÑ ?àò?ñâ?ïú?ã§.
RNAseq_ALL03_E1z <- enrichDGN(names(RNAseq_ALL03_E1y))
RNAseq_ALL03_N1z <- enrichDGN(names(RNAseq_ALL03_N1y))
RNAseq_ALL03_E3z <- enrichDGN(names(RNAseq_ALL03_E3y))
RNAseq_ALL03_N3z <- enrichDGN(names(RNAseq_ALL03_N3y))
RNAseq_ALL03_E4z <- enrichDGN(names(RNAseq_ALL03_E4y))
RNAseq_ALL03_N4z <- enrichDGN(names(RNAseq_ALL03_N4y))
RNAseq_ALL03_E5z <- enrichDGN(names(RNAseq_ALL03_E5y))
RNAseq_ALL03_N5z <- enrichDGN(names(RNAseq_ALL03_N5y))
RNAseq_ALL03_E6z <- enrichDGN(names(RNAseq_ALL03_E6y))
RNAseq_ALL03_N6z <- enrichDGN(names(RNAseq_ALL03_N6y))
RNAseq_ALL03_E8z <- enrichDGN(names(RNAseq_ALL03_E8y))
RNAseq_ALL03_N8z <- enrichDGN(names(RNAseq_ALL03_N8y))
RNAseq_ALL03_E9z <- enrichDGN(names(RNAseq_ALL03_E9y))
RNAseq_ALL03_N9z <- enrichDGN(names(RNAseq_ALL03_N9y))
RNAseq_ALL03_E0z <- enrichDGN(names(RNAseq_ALL03_E0y))
RNAseq_ALL03_N0z <- enrichDGN(names(RNAseq_ALL03_N0y))
RNAseq_ALL03_A1z <- enrichDGN(names(RNAseq_ALL03_A1y))

#setReadable?ùÑ ?ù¥?ö©?ïò?ó¨ ENTREZ IDÎ°? ?êú ?ú†?†Ñ?ûê ?†ïÎ≥¥Î?? ?ú†?†Ñ?ûê ?ã¨Î≥ºÎ°ú Î≥Ä?ôò, org.Hs.eg.db?äî ?ù∏Í∞? ?ú†?†Ñ?ûê ?†ïÎ≥? ?ç∞?ù¥?Ñ∞Î≤†Ïù¥?ä§?ûÑ.
RNAseq_ALL03_E1z <- setReadable(RNAseq_ALL03_E1z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_N1z <- setReadable(RNAseq_ALL03_N1z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_E3z <- setReadable(RNAseq_ALL03_E3z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_N3z <- setReadable(RNAseq_ALL03_N3z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_E4z <- setReadable(RNAseq_ALL03_E4z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_N4z <- setReadable(RNAseq_ALL03_N4z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_E5z <- setReadable(RNAseq_ALL03_E5z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_N5z <- setReadable(RNAseq_ALL03_N5z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_E6z <- setReadable(RNAseq_ALL03_E6z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_N6z <- setReadable(RNAseq_ALL03_N6z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_E8z <- setReadable(RNAseq_ALL03_E8z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_N8z <- setReadable(RNAseq_ALL03_N8z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_E9z <- setReadable(RNAseq_ALL03_E9z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_N9z <- setReadable(RNAseq_ALL03_N9z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_E0z <- setReadable(RNAseq_ALL03_E0z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_N0z <- setReadable(RNAseq_ALL03_N0z, 'org.Hs.eg.db', 'ENTREZID')
RNAseq_ALL03_A1z <- setReadable(RNAseq_ALL03_A1z, 'org.Hs.eg.db', 'ENTREZID')

#1800*1800
#cnetplot??Ä ?ú†?†Ñ?ûê Í∏∞Îä•Í∞? ?Ñ§?ä∏?õå?Å¨ ?ãúÍ∞ÅÌôî ?ïò?äî ?ó≠?ï†
#enrichDGN?óê?Ñú Ï∂úÎ†•?êú ?ú†?†Ñ?ûê??Ä Í¥Ä?†®?êú ?ÉùÎ¨ºÌïô?†Å Í≥ºÏ†ï?ùÑ ?ëú?ãú
#colorEdge=TURE?äî ?Ö∏?ìú Í∞? ?ó∞Í≤∞Ïùò ?Éâ?ùÑ ?ôú?Ñ±?ôî?ïò?ó¨ Í∞ôÏ?Ä Í∑∏Î£π?óê ?Üç?ïò?äî Í∏∞Îä•?ÅºÎ¶? ?Éâ?ÉÅ?ùÑ Í≥µÏú†?ïòÍ≤åÎÅî?ï®.
cp1 <- cnetplot(RNAseq_ALL03_E1z, colorEdge = TRUE, foldChange=RNAseq_ALL03_E1y, node_label = "all")
cp2 <- cnetplot(RNAseq_ALL03_N1z, colorEdge = TRUE, foldChange=RNAseq_ALL03_N1y, node_label = "all")
cp3 <- cnetplot(RNAseq_ALL03_E3z, colorEdge = TRUE, foldChange=RNAseq_ALL03_E3y, node_label = "all")
cp4 <- cnetplot(RNAseq_ALL03_N3z, colorEdge = TRUE, foldChange=RNAseq_ALL03_N3y, node_label = "all")
cp1 <- cnetplot(RNAseq_ALL03_E4z, colorEdge = TRUE, foldChange=RNAseq_ALL03_E4y, node_label = "all")
cp2 <- cnetplot(RNAseq_ALL03_N4z, colorEdge = TRUE, foldChange=RNAseq_ALL03_N4y, node_label = "all")
cp3 <- cnetplot(RNAseq_ALL03_E5z, colorEdge = TRUE, foldChange=RNAseq_ALL03_E5y, node_label = "all")
cp4 <- cnetplot(RNAseq_ALL03_N5z, colorEdge = TRUE, foldChange=RNAseq_ALL03_N5y, node_label = "all")
cp1 <- cnetplot(RNAseq_ALL03_E6z, colorEdge = TRUE, foldChange=RNAseq_ALL03_E6y, node_label = "all")
cp2 <- cnetplot(RNAseq_ALL03_N6z, colorEdge = TRUE, foldChange=RNAseq_ALL03_N6y, node_label = "all")
cp3 <- cnetplot(RNAseq_ALL03_E8z, colorEdge = TRUE, foldChange=RNAseq_ALL03_E8y, node_label = "all")
cp4 <- cnetplot(RNAseq_ALL03_N8z, colorEdge = TRUE, foldChange=RNAseq_ALL03_N8y, node_label = "all")
cp1 <- cnetplot(RNAseq_ALL03_E9z, colorEdge = TRUE, foldChange=RNAseq_ALL03_E9y, node_label = "all")
cp2 <- cnetplot(RNAseq_ALL03_N9z, colorEdge = TRUE, foldChange=RNAseq_ALL03_N9y, node_label = "all")
cp3 <- cnetplot(RNAseq_ALL03_E0z, colorEdge = TRUE, foldChange=RNAseq_ALL03_E0y, node_label = "all")
cp4 <- cnetplot(RNAseq_ALL03_N0z, colorEdge = TRUE, foldChange=RNAseq_ALL03_N0y, node_label = "all")
cp5 <- cnetplot(RNAseq_ALL03_A1z, colorEdge = TRUE, foldChange=RNAseq_ALL03_A1y, node_label = "all")
ggsave("Output_20240729/32_cnetplot_A1.pdf", plot=cp5, width=15,height=15)
ggsave("Output_20240729/32_cnetplot_E0.pdf", plot=cp3, width=15,height=15)

##?†ÑÏ≤? ??Ä?ÉÅ?ûê ??Ä?ÉÅ T1 vs T7 KEGG dot plot
ggplot(ERR_A1_KGx, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_point(aes(color = -log10(P.value)), size = 5) +
  coord_flip() +
  labs(
    x = "KEGG Pathway",
    y = "-log10(P-value)",
    title = "KEGG Enrichment Dot Plot (A1 Group)"
  ) +
  theme_minimal() +
  scale_color_gradient(low = "skyblue", high = "darkblue")

LIST_A1 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_A1 < 0.05 & abs(FC_A1) > 0.045), PV_A1,FC_A1)
LIST_A1$ENSEMBL <- rownames(LIST_A1)
names(LIST_A1)[1:2] <- c("PV","FC")

A1_volcano <- data.frame(
  logFC = PVFC_DEG03$FC_A1,
  pvalue = PVFC_DEG03$PV_A1,
  gene = rownames(PVFC_DEG03)
)

A1_volcano$logP <- -log10(A1_volcano$pvalue)
A1_volcano$significance <- "Not Sig"
A1_volcano$significance[A1_volcano$pvalue < 0.05 & A1_volcano$logFC > 0.1] <- "Up"
A1_volcano$significance[A1_volcano$pvalue < 0.05 & A1_volcano$logFC < -0.1] <- "Down"

# Volcano plot Í∑∏Î¶¨Í∏?
library(ggplot2)

ggplot(A1_volcano, aes(x = logFC, y = logP)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot (A1 group)",
    x = "Fold Change",
    y = "-Log10(p-value)"
  ) +
  theme_minimal()

#pvalueÎ•? 0.01 ?ù¥?ïòÎ°? ?ÇÆÏ∂∞ÏÑú ?ôï?ù∏?ù∏
############################################################################################################
LIST_E1 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E1 < 0.001 & abs(FC_E1) > 0.15), PV_E1,FC_E1); LIST_E1$ENSEMBL <- rownames(LIST_E1); dim(LIST_E1) #5,41
LIST_N1 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N1 < 0.01  & abs(FC_N1) > 0.09), PV_N1,FC_N1); LIST_N1$ENSEMBL <- rownames(LIST_N1); dim(LIST_N1) #4,41
LIST_E3 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E3 < 0.01  & abs(FC_E3) > 0.13), PV_E3,FC_E3); LIST_E3$ENSEMBL <- rownames(LIST_E3); dim(LIST_E3) #6,41
LIST_N3 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N3 < 0.001 & abs(FC_N3) > 0.11), PV_N3,FC_N3); LIST_N3$ENSEMBL <- rownames(LIST_N3); dim(LIST_N3) #6,41
LIST_E4 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E4 < 0.001 & abs(FC_E4) > 0.15), PV_E4,FC_E4); LIST_E4$ENSEMBL <- rownames(LIST_E4); dim(LIST_E4) #5,41
LIST_N4 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N4 < 0.01  & abs(FC_N4) > 0.09), PV_N4,FC_N4); LIST_N4$ENSEMBL <- rownames(LIST_N4); dim(LIST_N4) #4,41
LIST_E5 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E5 < 0.01  & abs(FC_E5) > 0.13), PV_E5,FC_E5); LIST_E5$ENSEMBL <- rownames(LIST_E5); dim(LIST_E5) #6,41
LIST_N5 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N5 < 0.001 & abs(FC_N5) > 0.11), PV_N5,FC_N5); LIST_N5$ENSEMBL <- rownames(LIST_N5); dim(LIST_N5) #6,41
LIST_E6 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E6 < 0.001 & abs(FC_E6) > 0.15), PV_E6,FC_E6); LIST_E6$ENSEMBL <- rownames(LIST_E6); dim(LIST_E6) #5,41
LIST_N6 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N6 < 0.01  & abs(FC_N6) > 0.09), PV_N6,FC_N6); LIST_N6$ENSEMBL <- rownames(LIST_N6); dim(LIST_N6) #4,41
LIST_E8 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E8 < 0.01  & abs(FC_E8) > 0.13), PV_E8,FC_E8); LIST_E8$ENSEMBL <- rownames(LIST_E8); dim(LIST_E8) #6,41
LIST_N8 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N8 < 0.001 & abs(FC_N8) > 0.11), PV_N8,FC_N8); LIST_N8$ENSEMBL <- rownames(LIST_N8); dim(LIST_N8) #6,41
LIST_E9 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E9 < 0.001 & abs(FC_E9) > 0.15), PV_E9,FC_E9); LIST_E9$ENSEMBL <- rownames(LIST_E9); dim(LIST_E9) #5,41
LIST_N9 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N9 < 0.01  & abs(FC_N9) > 0.09), PV_N9,FC_N9); LIST_N9$ENSEMBL <- rownames(LIST_N9); dim(LIST_N9) #4,41
LIST_E0 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_E0 < 0.01  & abs(FC_E0) > 0.13), PV_E0,FC_E0); LIST_E0$ENSEMBL <- rownames(LIST_E0); dim(LIST_E0) #6,41
LIST_N0 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_N0 < 0.001 & abs(FC_N0) > 0.11), PV_N0,FC_N0); LIST_N0$ENSEMBL <- rownames(LIST_N0); dim(LIST_N0) #6,41
LIST_A1 <- dplyr::select(dplyr::filter(PVFC_DEG, PV_A1 < 0.001 & abs(FC_A1) > 0.11), PV_A1,FC_A1); LIST_A1$ENSEMBL <- rownames(LIST_A1); dim(LIST_A1) #6,41
LIST_EN <- Reduce(union, list(rownames(LIST_E1),rownames(LIST_N1),rownames(LIST_E3),rownames(LIST_N3),rownames(LIST_E4),rownames(LIST_N4),rownames(LIST_E5),rownames(LIST_N5),rownames(LIST_E6),rownames(LIST_N6),rownames(LIST_E8),rownames(LIST_N8),rownames(LIST_E9),rownames(LIST_N9),rownames(LIST_E0),rownames(LIST_N0)))
#Î™®Îì† ?ú†?ùò?ïú ?ú†?†Ñ?ûê Î¶¨Ïä§?ä∏Îß? Ï∂îÏ∂ú
RNAseq_ALL04 <- RNAseq_ALL03[which(rownames(RNAseq_ALL03) %in% LIST_EN),]
RNAseq_ALL04$ENSEMBL <- rownames(RNAseq_ALL04)
RNAseq_ALL04 <- merge(x=RNAseq_ALL04,y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL04) <- RNAseq_ALL04$ENSEMBL; RNAseq_ALL04 <- RNAseq_ALL04[,-1]

#?ú†?†Ñ?ûê ?ã¨Î≥? ?ç∞?ù¥?Ñ∞ Ï∂îÍ?Ä
SYMBOL_E1 <- merge(x=LIST_E1,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_E1) <- SYMBOL_E1$ENSEMBL; SYMBOL_E1 <- SYMBOL_E1[,-1]; names(SYMBOL_E1)[1:2] <- c("PV","FC")
SYMBOL_N1 <- merge(x=LIST_N1,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_N1) <- SYMBOL_N1$ENSEMBL; SYMBOL_N1 <- SYMBOL_N1[,-1]; names(SYMBOL_N1)[1:2] <- c("PV","FC")
SYMBOL_E3 <- merge(x=LIST_E3,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_E3) <- SYMBOL_E3$ENSEMBL; SYMBOL_E3 <- SYMBOL_E3[,-1]; names(SYMBOL_E3)[1:2] <- c("PV","FC")
SYMBOL_N3 <- merge(x=LIST_N3,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_N3) <- SYMBOL_N3$ENSEMBL; SYMBOL_N3 <- SYMBOL_N3[,-1]; names(SYMBOL_N3)[1:2] <- c("PV","FC")
SYMBOL_E4 <- merge(x=LIST_E4,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_E4) <- SYMBOL_E4$ENSEMBL; SYMBOL_E4 <- SYMBOL_E4[,-1]; names(SYMBOL_E4)[1:2] <- c("PV","FC")
SYMBOL_N4 <- merge(x=LIST_N4,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_N4) <- SYMBOL_N4$ENSEMBL; SYMBOL_N4 <- SYMBOL_N4[,-1]; names(SYMBOL_N4)[1:2] <- c("PV","FC")
SYMBOL_E5 <- merge(x=LIST_E5,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_E5) <- SYMBOL_E5$ENSEMBL; SYMBOL_E5 <- SYMBOL_E5[,-1]; names(SYMBOL_E5)[1:2] <- c("PV","FC")
SYMBOL_N5 <- merge(x=LIST_N5,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_N5) <- SYMBOL_N5$ENSEMBL; SYMBOL_N5 <- SYMBOL_N5[,-1]; names(SYMBOL_N5)[1:2] <- c("PV","FC")
SYMBOL_E6 <- merge(x=LIST_E6,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_E6) <- SYMBOL_E6$ENSEMBL; SYMBOL_E6 <- SYMBOL_E6[,-1]; names(SYMBOL_E6)[1:2] <- c("PV","FC")
SYMBOL_N6 <- merge(x=LIST_N6,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_N6) <- SYMBOL_N6$ENSEMBL; SYMBOL_N6 <- SYMBOL_N6[,-1]; names(SYMBOL_N6)[1:2] <- c("PV","FC")
SYMBOL_E8 <- merge(x=LIST_E8,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_E8) <- SYMBOL_E8$ENSEMBL; SYMBOL_E8 <- SYMBOL_E8[,-1]; names(SYMBOL_E8)[1:2] <- c("PV","FC")
SYMBOL_N8 <- merge(x=LIST_N8,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_N8) <- SYMBOL_N8$ENSEMBL; SYMBOL_N8 <- SYMBOL_N8[,-1]; names(SYMBOL_N8)[1:2] <- c("PV","FC")
SYMBOL_E9 <- merge(x=LIST_E9,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_E9) <- SYMBOL_E9$ENSEMBL; SYMBOL_E9 <- SYMBOL_E9[,-1]; names(SYMBOL_E9)[1:2] <- c("PV","FC")
SYMBOL_N9 <- merge(x=LIST_N9,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_N9) <- SYMBOL_N9$ENSEMBL; SYMBOL_N9 <- SYMBOL_N9[,-1]; names(SYMBOL_N9)[1:2] <- c("PV","FC")
SYMBOL_E0 <- merge(x=LIST_E0,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_E0) <- SYMBOL_E0$ENSEMBL; SYMBOL_E0 <- SYMBOL_E0[,-1]; names(SYMBOL_E0)[1:2] <- c("PV","FC")
SYMBOL_N0 <- merge(x=LIST_N0,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_N0) <- SYMBOL_N0$ENSEMBL; SYMBOL_N0 <- SYMBOL_N0[,-1]; names(SYMBOL_N0)[1:2] <- c("PV","FC")
SYMBOL_A1 <- merge(x=LIST_A1,y=RNAseq_GENEs, by="ENSEMBL"); rownames(SYMBOL_A1) <- SYMBOL_A1$ENSEMBL; SYMBOL_A1 <- SYMBOL_A1[,-1]; names(SYMBOL_A1)[1:2] <- c("PV","FC")


SYMBOL_EN <- rbind(SYMBOL_E1,SYMBOL_N1,SYMBOL_E3,SYMBOL_N3)
SYMBOL_EN <- SYMBOL_EN[!duplicated(SYMBOL_EN),]
head(SYMBOL_EN)

pheatmap(RNAseq_ALL04[,1:40], labels_row = RNAseq_ALL04$SYMBOL, cutree_cols=2, color=custom_col, fontsize=4)

#Ï∞®Îì±Î∞úÌòÑ?êú ?ú†?†Ñ?ûê?ì§Îß? Ï∂îÏ∂ú?ïò?ó¨ ENSEMBL ?ù¥Î¶ÑÍ≥º Î≥ëÌï©
RNAseq_ALL03_E1x <- as.data.frame(RNAseq_ALL03_E1[which(rownames(RNAseq_ALL03_E1) %in% rownames(LIST_E1)),]); RNAseq_ALL03_E1x$ENSEMBL <- rownames(RNAseq_ALL03_E1x); RNAseq_ALL03_E1x <- merge(x=RNAseq_ALL03_E1x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_E1x) <- RNAseq_ALL03_E1x$ENSEMBL; RNAseq_ALL03_E1x <- RNAseq_ALL03_E1x[,-1]
RNAseq_ALL03_N1x <- as.data.frame(RNAseq_ALL03_N1[which(rownames(RNAseq_ALL03_N1) %in% rownames(LIST_N1)),]); RNAseq_ALL03_N1x$ENSEMBL <- rownames(RNAseq_ALL03_N1x); RNAseq_ALL03_N1x <- merge(x=RNAseq_ALL03_N1x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_N1x) <- RNAseq_ALL03_N1x$ENSEMBL; RNAseq_ALL03_N1x <- RNAseq_ALL03_N1x[,-1]
RNAseq_ALL03_E3x <- as.data.frame(RNAseq_ALL03_E3[which(rownames(RNAseq_ALL03_E3) %in% rownames(LIST_E3)),]); RNAseq_ALL03_E3x$ENSEMBL <- rownames(RNAseq_ALL03_E3x); RNAseq_ALL03_E3x <- merge(x=RNAseq_ALL03_E3x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_E3x) <- RNAseq_ALL03_E3x$ENSEMBL; RNAseq_ALL03_E3x <- RNAseq_ALL03_E3x[,-1]
RNAseq_ALL03_N3x <- as.data.frame(RNAseq_ALL03_N3[which(rownames(RNAseq_ALL03_N3) %in% rownames(LIST_N3)),]); RNAseq_ALL03_N3x$ENSEMBL <- rownames(RNAseq_ALL03_N3x); RNAseq_ALL03_N3x <- merge(x=RNAseq_ALL03_N3x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_N3x) <- RNAseq_ALL03_N3x$ENSEMBL; RNAseq_ALL03_N3x <- RNAseq_ALL03_N3x[,-1]
RNAseq_ALL03_E4x <- as.data.frame(RNAseq_ALL03_E4[which(rownames(RNAseq_ALL03_E4) %in% rownames(LIST_E4)),]); RNAseq_ALL03_E4x$ENSEMBL <- rownames(RNAseq_ALL03_E4x); RNAseq_ALL03_E4x <- merge(x=RNAseq_ALL03_E4x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_E4x) <- RNAseq_ALL03_E4x$ENSEMBL; RNAseq_ALL03_E4x <- RNAseq_ALL03_E4x[,-1]
RNAseq_ALL03_N4x <- as.data.frame(RNAseq_ALL03_N4[which(rownames(RNAseq_ALL03_N4) %in% rownames(LIST_N4)),]); RNAseq_ALL03_N4x$ENSEMBL <- rownames(RNAseq_ALL03_N4x); RNAseq_ALL03_N4x <- merge(x=RNAseq_ALL03_N4x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_N4x) <- RNAseq_ALL03_N4x$ENSEMBL; RNAseq_ALL03_N4x <- RNAseq_ALL03_N4x[,-1]
RNAseq_ALL03_E5x <- as.data.frame(RNAseq_ALL03_E5[which(rownames(RNAseq_ALL03_E5) %in% rownames(LIST_E5)),]); RNAseq_ALL03_E5x$ENSEMBL <- rownames(RNAseq_ALL03_E5x); RNAseq_ALL03_E5x <- merge(x=RNAseq_ALL03_E5x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_E5x) <- RNAseq_ALL03_E5x$ENSEMBL; RNAseq_ALL03_E5x <- RNAseq_ALL03_E5x[,-1]
RNAseq_ALL03_N5x <- as.data.frame(RNAseq_ALL03_N5[which(rownames(RNAseq_ALL03_N5) %in% rownames(LIST_N5)),]); RNAseq_ALL03_N5x$ENSEMBL <- rownames(RNAseq_ALL03_N5x); RNAseq_ALL03_N5x <- merge(x=RNAseq_ALL03_N5x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_N5x) <- RNAseq_ALL03_N5x$ENSEMBL; RNAseq_ALL03_N5x <- RNAseq_ALL03_N5x[,-1]
RNAseq_ALL03_E6x <- as.data.frame(RNAseq_ALL03_E6[which(rownames(RNAseq_ALL03_E6) %in% rownames(LIST_E6)),]); RNAseq_ALL03_E6x$ENSEMBL <- rownames(RNAseq_ALL03_E6x); RNAseq_ALL03_E6x <- merge(x=RNAseq_ALL03_E6x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_E6x) <- RNAseq_ALL03_E6x$ENSEMBL; RNAseq_ALL03_E6x <- RNAseq_ALL03_E6x[,-1]
RNAseq_ALL03_N6x <- as.data.frame(RNAseq_ALL03_N6[which(rownames(RNAseq_ALL03_N6) %in% rownames(LIST_N6)),]); RNAseq_ALL03_N6x$ENSEMBL <- rownames(RNAseq_ALL03_N6x); RNAseq_ALL03_N6x <- merge(x=RNAseq_ALL03_N6x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_N6x) <- RNAseq_ALL03_N6x$ENSEMBL; RNAseq_ALL03_N6x <- RNAseq_ALL03_N6x[,-1]
RNAseq_ALL03_E8x <- as.data.frame(RNAseq_ALL03_E8[which(rownames(RNAseq_ALL03_E8) %in% rownames(LIST_E8)),]); RNAseq_ALL03_E8x$ENSEMBL <- rownames(RNAseq_ALL03_E8x); RNAseq_ALL03_E8x <- merge(x=RNAseq_ALL03_E8x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_E8x) <- RNAseq_ALL03_E8x$ENSEMBL; RNAseq_ALL03_E8x <- RNAseq_ALL03_E8x[,-1]
RNAseq_ALL03_N8x <- as.data.frame(RNAseq_ALL03_N8[which(rownames(RNAseq_ALL03_N8) %in% rownames(LIST_N8)),]); RNAseq_ALL03_N8x$ENSEMBL <- rownames(RNAseq_ALL03_N8x); RNAseq_ALL03_N8x <- merge(x=RNAseq_ALL03_N8x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_N8x) <- RNAseq_ALL03_N8x$ENSEMBL; RNAseq_ALL03_N8x <- RNAseq_ALL03_N8x[,-1]
RNAseq_ALL03_E9x <- as.data.frame(RNAseq_ALL03_E9[which(rownames(RNAseq_ALL03_E9) %in% rownames(LIST_E9)),]); RNAseq_ALL03_E9x$ENSEMBL <- rownames(RNAseq_ALL03_E9x); RNAseq_ALL03_E9x <- merge(x=RNAseq_ALL03_E9x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_E9x) <- RNAseq_ALL03_E9x$ENSEMBL; RNAseq_ALL03_E9x <- RNAseq_ALL03_E9x[,-1]
RNAseq_ALL03_N9x <- as.data.frame(RNAseq_ALL03_N9[which(rownames(RNAseq_ALL03_N9) %in% rownames(LIST_N9)),]); RNAseq_ALL03_N9x$ENSEMBL <- rownames(RNAseq_ALL03_N9x); RNAseq_ALL03_N9x <- merge(x=RNAseq_ALL03_N9x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_N9x) <- RNAseq_ALL03_N9x$ENSEMBL; RNAseq_ALL03_N9x <- RNAseq_ALL03_N9x[,-1]
RNAseq_ALL03_E0x <- as.data.frame(RNAseq_ALL03_E0[which(rownames(RNAseq_ALL03_E0) %in% rownames(LIST_E0)),]); RNAseq_ALL03_E0x$ENSEMBL <- rownames(RNAseq_ALL03_E0x); RNAseq_ALL03_E0x <- merge(x=RNAseq_ALL03_E0x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_E0x) <- RNAseq_ALL03_E0x$ENSEMBL; RNAseq_ALL03_E0x <- RNAseq_ALL03_E0x[,-1]
RNAseq_ALL03_N0x <- as.data.frame(RNAseq_ALL03_N0[which(rownames(RNAseq_ALL03_N0) %in% rownames(LIST_N0)),]); RNAseq_ALL03_N0x$ENSEMBL <- rownames(RNAseq_ALL03_N0x); RNAseq_ALL03_N0x <- merge(x=RNAseq_ALL03_N0x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_N0x) <- RNAseq_ALL03_N0x$ENSEMBL; RNAseq_ALL03_N0x <- RNAseq_ALL03_N0x[,-1]
RNAseq_ALL03_A1x <- as.data.frame(RNAseq_ALL03_A1[which(rownames(RNAseq_ALL03_A1) %in% rownames(LIST_A1)),]); RNAseq_ALL03_A1x$ENSEMBL <- rownames(RNAseq_ALL03_A1x); RNAseq_ALL03_A1x <- merge(x=RNAseq_ALL03_A1x, y=RNAseq_GENEs, by="ENSEMBL"); rownames(RNAseq_ALL03_A1x) <- RNAseq_ALL03_A1x$ENSEMBL; RNAseq_ALL03_A1x <- RNAseq_ALL03_A1x[,-1]

#Í∑∏Î£π Î≥? ?ú†?†Ñ?ûê Î∞úÌòÑ ?ñë?ÉÅ 
E1pm <- pheatmap(dplyr::select(RNAseq_ALL03_E1x, -SYMBOL), labels_row = RNAseq_ALL03_E1x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
N1pm <- pheatmap(dplyr::select(RNAseq_ALL03_N1x, -SYMBOL), labels_row = RNAseq_ALL03_N1x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
E3pm <- pheatmap(dplyr::select(RNAseq_ALL03_E3x, -SYMBOL), labels_row = RNAseq_ALL03_E3x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
N3pm <- pheatmap(dplyr::select(RNAseq_ALL03_N3x, -SYMBOL), labels_row = RNAseq_ALL03_N3x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
E4pm <- pheatmap(dplyr::select(RNAseq_ALL03_E4x, -SYMBOL), labels_row = RNAseq_ALL03_E4x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=8)
N4pm <- pheatmap(dplyr::select(RNAseq_ALL03_N4x, -SYMBOL), labels_row = RNAseq_ALL03_N4x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
E5pm <- pheatmap(dplyr::select(RNAseq_ALL03_E5x, -SYMBOL), labels_row = RNAseq_ALL03_E5x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=8)
N5pm <- pheatmap(dplyr::select(RNAseq_ALL03_N5x, -SYMBOL), labels_row = RNAseq_ALL03_N5x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
E6pm <- pheatmap(dplyr::select(RNAseq_ALL03_E6x, -SYMBOL), labels_row = RNAseq_ALL03_E6x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
N6pm <- pheatmap(dplyr::select(RNAseq_ALL03_N6x, -SYMBOL), labels_row = RNAseq_ALL03_N6x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
E8pm <- pheatmap(dplyr::select(RNAseq_ALL03_E8x, -SYMBOL), labels_row = RNAseq_ALL03_E8x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=7)
N8pm <- pheatmap(dplyr::select(RNAseq_ALL03_N8x, -SYMBOL), labels_row = RNAseq_ALL03_N8x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
E9pm <- pheatmap(dplyr::select(RNAseq_ALL03_E9x, -SYMBOL), labels_row = RNAseq_ALL03_E9x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
N9pm <- pheatmap(dplyr::select(RNAseq_ALL03_N9x, -SYMBOL), labels_row = RNAseq_ALL03_N9x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
E0pm <- pheatmap(dplyr::select(RNAseq_ALL03_E0x, -SYMBOL), labels_row = RNAseq_ALL03_E0x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=8)
N0pm <- pheatmap(dplyr::select(RNAseq_ALL03_N0x, -SYMBOL), labels_row = RNAseq_ALL03_N0x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
A1pm <- pheatmap(dplyr::select(RNAseq_ALL03_A1x, -SYMBOL), labels_row = RNAseq_ALL03_A1x$SYMBOL, cutree_cols=2, color=custom_col, fontsize=15)
dim(RNAseq_ALL03_A1x)
dim(RNAseq_ALL03_E0x)

ggsave("Output_20240729/31_Heatmap_E1.pdf", plot = E1pm$gtable, width = 8, height = 4)
ggsave("Output_20240729/31_Heatmap_N1.pdf", plot = N1pm$gtable, width = 8, height = 4)
ggsave("Output_20240729/31_Heatmap_E3.pdf", plot = E3pm$gtable, width = 5, height = 4)
ggsave("Output_20240729/31_Heatmap_N3.pdf", plot = N3pm$gtable, width =10, height = 4)
ggsave("Output_20240729/31_Heatmap_E4.pdf", plot = E4pm$gtable, width = 8, height = 4)
ggsave("Output_20240729/31_Heatmap_N4.pdf", plot = N4pm$gtable, width = 8, height = 4)
ggsave("Output_20240729/31_Heatmap_E5.pdf", plot = E5pm$gtable, width = 5, height = 4)
ggsave("Output_20240729/31_Heatmap_N5.pdf", plot = N5pm$gtable, width =10, height = 4)
ggsave("Output_20240729/31_Heatmap_E6.pdf", plot = E6pm$gtable, width = 8, height = 4)
ggsave("Output_20240729/31_Heatmap_N6.pdf", plot = N6pm$gtable, width = 8, height = 4)
ggsave("Output_20240729/31_Heatmap_E8.pdf", plot = E8pm$gtable, width = 5, height = 4)
ggsave("Output_20240729/31_Heatmap_N8.pdf", plot = N8pm$gtable, width =10, height = 4)
ggsave("Output_20240729/31_Heatmap_E9.pdf", plot = E9pm$gtable, width = 8, height = 4)
ggsave("Output_20240729/31_Heatmap_N9.pdf", plot = N9pm$gtable, width = 8, height = 4)
ggsave("Output_20240729/31_Heatmap_E0.pdf", plot = E0pm$gtable, width = 5, height = 4)
ggsave("Output_20240729/31_Heatmap_N0.pdf", plot = N0pm$gtable, width =10, height = 4)

#?ú†?†Ñ?ûê ?ã¨Î≥? ?ãπ ID?ùò valueÎ•? ?ôï?ù∏?ï† ?àò ?ûàÍ≤åÎÅî Íµ¨Ï°∞Î•? Î≥Ä?ôò
colnames(RNAseq_ALL04)
RNAseq_ALL05 <- reshape2::melt(RNAseq_ALL04, id.vars = c("SYMBOL")); names(RNAseq_ALL05)[2] <- "ID"
head(RNAseq_ALL05)

#Í∑∏Î£πÎ≥? timepoint Ï∂îÍ?Ä
RNAseq_ALL05$Time <- substr(RNAseq_ALL05$ID,12,12)
RNAseq_ALL05$ID <- substr(RNAseq_ALL05$ID,1,9)
PROPEL03c$ID <- rownames(PROPEL03c)

#?ïÑ?ûò ?ò§Î•?, IDÍ∞Ä ?ùºÏπòÌïòÏßÄ ?ïäÍ±∞ÎÇò RNA_seq ?åå?ùº?óê?Ñú ?ÑàÎ¨? ÎßéÏïÑ Ï§ëÎ≥µ ?ò§Î•? ?ïà?êò?èÑ Í∑∏ÎÉ• ?ï®...
RNAseq_ALL05 <- merge(x=RNAseq_ALL05, y=PROPEL03c[,c(1,3,11)], by="ID",all.x=TRUE)

head(RNAseq_ALL05)
#?ã§?ãú ?ãú?ûë
RNAseq_ALL05$C1A <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C1A)
RNAseq_ALL05$C1C <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C1C)
RNAseq_ALL05$C1D <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C1D)
RNAseq_ALL05$C1E <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C1E)
RNAseq_ALL05$C2A <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C2A)
RNAseq_ALL05$C2C <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C2C)
RNAseq_ALL05$C2D <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C2D)
RNAseq_ALL05$C2E <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C2E)
RNAseq_ALL05$C3A <- paste0(RNAseq_ALL05$Time,RNAseq_ALL05$C3A)
                           
#Í∞? timepoint Î≥? ?Éò?îå ?àò 20/20 ?ôï?ù∏ Symbol [1,3]??Ä Ï≤´Î≤àÏß? ?ñâ?óê 3Î≤àÏß∏ ?ó¥?ù¥ ?ã¨Î≥ºÏúºÎ°? ?ì±Î°ùÎêò?ñ¥ ?ûàÍ∏? ?ïåÎ¨∏Ïóê Ï∂îÏ∂ú?ï®.
#Ï≤´Î≤àÏß? ?ñâ??Ä Í∞Ä?û• ?ú†?ùò?ïú ?ã¨Î≥ºÏù¥Í∏? ?ïåÎ¨∏Ïóê ?ùòÎØ∏Í?Ä ?ûà?ùå.
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[1,3])$C1A)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N1[1,3])$C1A)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E3[1,3])$C1C)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N3[1,3])$C1C)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E4[1,3])$C1C)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N4[1,3])$C1C)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E5[1,3])$C1C)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N5[1,3])$C1C)

table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E6[1,3])$C2A)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N6[1,3])$C2A)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E8[1,3])$C2C)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N8[1,3])$C2C) # ?óÜ?ùå
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E9[1,3])$C2D)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N9[1,3])$C2D)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E0[1,3])$C2E)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N0[1,3])$C2E)
table(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_A1[1,3])$C3A) # ?óÜ?ùå?ùå

#?ú†?ùò?ïú ?ú†?†Ñ?ûê?ù∏ÏßÄ ?ôï?ù∏?ùÑ ?úÑ?ïú Í≥ºÏ†ï
#[1,3]?óê ?úÑÏπòÌïú THOC2 ?ú†?†Ñ?ûê?ùò value Ï∂îÏ∂ú
dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[1,3])
head(SYMBOL_E1)
#SYMBOL_E1Í∞Ä Í∞Ä?û• ?ú†?ùò?ïú ?àú?Ñú??ÄÎ°? ?Çò?ó¥?êò?óà?úº?ãà Ï≤? ?ñâ?óê ?ûà?äî ?ú†?†Ñ?ûêÍ∞Ä Í∞Ä?û• ?ú†?ùò?ïú ?ú†?†Ñ?ûêÍ∞Ä ÎßûÏùå

#C1A Í∑∏Î£π?óê?Ñú THOC2 ?ú†?†Ñ?ûê?ùò Î∞úÌòÑ?üâ Ï∞®Ïù¥Í∞Ä ?ú†?ùò?ïúÏßÄ ?ôï?ù∏
t.test(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[1,3] & C1A == "1E")$value,dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[1,3] & C1A == "7E")$value)
table(dplyr::filter(RNAseq_ALL05, SYMBOL == SYMBOL_E1[1,3])$C1A)

#ÎπÑÎ™®?àò?ùº?Ñú ensg numberÎ°? THOC2Î•? Î∂àÎü¨?ò® ?ã§?ùå t-test ?ú†?ùò?ïú Ï∞®Ïù¥Í∞Ä ?ûà?äîÏßÄ ?ôï?ù∏?ö©?úºÎ°? Î∂ÑÏÑù
t.test(RNAseq_ALL03_E1[which(rownames(RNAseq_ALL03_E1) == "ENSG00000125676.19"), ] ~ E1f)
RNAseq_ALL03_E1[which(rownames(RNAseq_ALL03_E1) == "ENSG00000125676.19"),]
RNAseq_GENEs[which(rownames(RNAseq_GENEs) == "ENSG00000125676.19"),]
RNAseq_ALL04[which(rownames(RNAseq_ALL04) == "ENSG00000125676.19"),]
dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[1,3] & C1A == "1E")$value
dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[1,3] & C1A == "7E")$value

#C1A Í∑∏Î£π?óê?Ñú ?ú†?ùò?ïú ?ú†?†Ñ?ûê?ì§?ùò timepoint?óê ?î∞Î•? Ï∞®Ïù¥ Î∂ÑÏÑù
E1p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[1,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_E1[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E1p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[2,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_E1[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E1p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[3,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_E1[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E1p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[4,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_E1[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E1p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E1[5,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_E1[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_E1.pdf", plot=grid.arrange(E1p1,E1p2,E1p3,E1p4,E1p5, ncol=5), width=15,height=5)

N1p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N1[1,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_N1[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N1p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N1[2,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_N1[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N1p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N1[3,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_N1[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N1p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N1[4,3]), aes(x=C1A, y=value, fill=C1A)) + geom_boxplot() + ggtitle(SYMBOL_N1[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N1.pdf", plot=grid.arrange(N1p1,N1p2,N1p3,N1p4, ncol=4), width=12,height=5)

#C1C Í∑∏Î£π?óê?Ñú ?ú†?ùò?ïú ?ú†?†Ñ?ûê timepoint?óê ?î∞?ùº Î∂ÑÏÑù?Ñù
E3p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E3[ 1,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_E3[ 1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E3p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E3[ 2,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_E3[ 2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E3p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E3[ 3,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_E3[ 3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E3p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E3[ 4,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_E3[ 4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E3p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E3[ 5,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_E3[ 5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E3p6 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E3[ 6,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_E3[ 6,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_E3.pdf", plot=grid.arrange(E3p1,E3p2,E3p3,E3p4,E3p5,E3p6, ncol=6), width=18,height=5)

N3p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N3[ 1,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_N3[ 1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N3p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N3[ 2,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_N3[ 2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N3p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N3[ 3,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_N3[ 3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N3p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N3[ 4,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_N3[ 4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N3p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N3[ 5,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_N3[ 5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N3p6 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N3[ 6,3]), aes(x=C1C, y=value, fill=C1C)) + geom_boxplot() + ggtitle(SYMBOL_N3[ 6,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N3.pdf", plot=grid.arrange(N3p1,N3p2,N3p3,N3p4,N3p5,N3p6, ncol=6), width=18,height=5)

#C1D
E4p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E4[ 1,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_E4[ 1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E4p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E4[ 2,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_E4[ 2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E4p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E4[ 3,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_E4[ 3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E4p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E4[ 4,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_E4[ 4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E4p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E4[ 5,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_E4[ 5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E4p6 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E4[ 6,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_E4[ 6,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_E4.pdf", plot=grid.arrange(E4p1,E4p2,E4p3,E4p4,E4p5,E4p6, ncol=6), width=18,height=5)

N4p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N4[ 1,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_N4[ 1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N4p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N4[ 2,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_N4[ 2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N4p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N4[ 3,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_N4[ 3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N4p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N4[ 4,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_N4[ 4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N4p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N4[ 5,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_N4[ 5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N4p6 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N4[ 6,3]), aes(x=C1D, y=value, fill=C1D)) + geom_boxplot() + ggtitle(SYMBOL_N4[ 6,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N4.pdf", plot=grid.arrange(N4p1,N4p2,N4p3,N4p4,N4p5,N4p6, ncol=6), width=18,height=5)

#C1E
E5p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E5[ 1,3]), aes(x=C1E, y=value, fill=C1E)) + geom_boxplot() + ggtitle(SYMBOL_E5[ 1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E5p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E5[ 2,3]), aes(x=C1E, y=value, fill=C1E)) + geom_boxplot() + ggtitle(SYMBOL_E5[ 2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E5p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E5[ 3,3]), aes(x=C1E, y=value, fill=C1E)) + geom_boxplot() + ggtitle(SYMBOL_E5[ 3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E5p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E5[ 4,3]), aes(x=C1E, y=value, fill=C1E)) + geom_boxplot() + ggtitle(SYMBOL_E5[ 4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E5p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E5[ 5,3]), aes(x=C1E, y=value, fill=C1E)) + geom_boxplot() + ggtitle(SYMBOL_E5[ 5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E5p6 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E5[ 6,3]), aes(x=C1E, y=value, fill=C1E)) + geom_boxplot() + ggtitle(SYMBOL_E5[ 6,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_E5.pdf", plot=grid.arrange(E5p1,E5p2,E5p3,E5p4,E5p5,E5p6, ncol=6), width=18,height=5)

N5p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N5[ 1,3]), aes(x=C1E, y=value, fill=C1E)) + geom_boxplot() + ggtitle(SYMBOL_N5[ 1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N5p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N5[ 2,3]), aes(x=C1E, y=value, fill=C1E)) + geom_boxplot() + ggtitle(SYMBOL_N5[ 2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N5.pdf", plot=grid.arrange(N5p1,N5p2, ncol=2), width=18,height=5)

#C2A
E6p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E6[1,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_E6[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E6p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E6[2,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_E6[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E6p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E6[3,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_E6[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E6p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E6[4,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_E6[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E6p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E6[5,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_E6[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_E6.pdf", plot=grid.arrange(E6p1,E6p2,E6p3,E6p4,E6p5, ncol=5), width=15,height=5)

N6p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N6[1,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_N6[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N6p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N6[2,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_N6[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N6p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N6[3,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_N6[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N6p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N6[4,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_N6[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N6p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N6[5,3]), aes(x=C2A, y=value, fill=C2A)) + geom_boxplot() + ggtitle(SYMBOL_N6[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N6.pdf", plot=grid.arrange(N6p1,N6p2,N6p3,N6p4,N6p5, ncol=5), width=12,height=5)

#C2C
E8p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E8[1,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_E8[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E8p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E8[2,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_E8[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E8p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E8[3,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_E8[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E8p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E8[4,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_E8[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E8p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E8[5,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_E8[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_E8.pdf", plot=grid.arrange(E8p1,E8p2,E8p3,E8p4,E8p5, ncol=5), width=15,height=5)

#?óÜ?ùå
N8p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N8[1,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_N8[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N8p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N8[2,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_N8[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N8p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N8[3,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_N8[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N8p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N8[4,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_N8[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N8p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N8[5,3]), aes(x=C2C, y=value, fill=C2C)) + geom_boxplot() + ggtitle(SYMBOL_N8[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N8.pdf", plot=grid.arrange(N8p1,N8p2,N8p3,N8p4,N8p5, ncol=5), width=12,height=5)

#C2D
E9p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E9[1,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_E9[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E9p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E9[2,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_E9[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E9p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E9[3,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_E9[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E9p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E9[4,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_E9[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E9p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E9[5,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_E9[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_E9.pdf", plot=grid.arrange(E9p1,E9p2,E9p3,E9p4,E9p5, ncol=5), width=15,height=5)

N9p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N9[1,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_N9[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N9p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N9[2,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_N9[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N9p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N9[3,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_N9[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N9p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N9[4,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_N9[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N9p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N9[5,3]), aes(x=C2D, y=value, fill=C2D)) + geom_boxplot() + ggtitle(SYMBOL_N9[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N9.pdf", plot=grid.arrange(N9p1,N9p2,N9p3,N9p4,N9p5, ncol=5), width=12,height=5)

#C2E
E0p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E0[1,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_E0[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E0p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E0[2,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_E0[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E0p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E0[3,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_E0[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E0p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E0[4,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_E0[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
E0p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_E0[5,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_E0[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1E","7E")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_E0.pdf", plot=grid.arrange(E0p1,E0p2,E0p3,E0p4,E0p5, ncol=5), width=15,height=5)

N0p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N0[1,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_N0[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N0p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N0[2,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_N0[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N0p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N0[3,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_N0[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N0p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N0[4,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_N0[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
N0p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_N0[5,3]), aes(x=C2E, y=value, fill=C2E)) + geom_boxplot() + ggtitle(SYMBOL_N0[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1N","7N")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N0.pdf", plot=grid.arrange(N0p1, ncol=1), width=12,height=5)

#C3A ?ò§Î•? #C3AÎ•? Ï∞æÏùÑ ?àò ?óÜ?ùå.
A1p1 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_A1[1,3]), aes(x=C3A, y=value, fill=C3A)) + geom_boxplot() + ggtitle(SYMBOL_A1[1,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1A","7A")), method = "t.test", label = "p.format")
A1p2 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_A1[2,3]), aes(x=C3A, y=value, fill=C3A)) + geom_boxplot() + ggtitle(SYMBOL_A1[2,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1A","7A")), method = "t.test", label = "p.format")
A1p3 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_A1[3,3]), aes(x=C3A, y=value, fill=C3A)) + geom_boxplot() + ggtitle(SYMBOL_A1[3,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1A","7A")), method = "t.test", label = "p.format")
A1p4 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_A1[4,3]), aes(x=C3A, y=value, fill=C3A)) + geom_boxplot() + ggtitle(SYMBOL_A1[4,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1A","7A")), method = "t.test", label = "p.format")
A1p5 <- ggplot(dplyr::filter(RNAseq_ALL05,SYMBOL==SYMBOL_A1[5,3]), aes(x=C3A, y=value, fill=C3A)) + geom_boxplot() + ggtitle(SYMBOL_A1[5,3]) + theme(legend.position = "none") + stat_compare_means(comparisons = list(c("1A","7A")), method = "t.test", label = "p.format")
ggsave("Output_20240729/31_boxplot_N0.pdf", plot=grid.arrange(A1p1, A1p2, A1p3, A1p4, A1p5, ncol=1), width=12,height=5)


#####################RECUR_LIST, WHOGR_LIST ?ôï?ù∏ ?ïÑ?öî
#RESUR, WHOGRÍ∞Ä Í∑∏Î£π Í∞? timepoint 1Í≥? 7?ùÑ ÎπÑÍµê?ïú Í≤ÉÏù∏Í∞Ä?
#GO bar plot Í∑∏Î¶¥ ?àò ?ûà?ùå.
library(org.Hs.eg.db)
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
head(RECUR_LIST)
EGO1 <- enrichGO(gene = RECUR_LIST, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
EGO2 <- enrichGO(gene = WHOGR_LIST, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
EGO1
EGO2
barplot(EGO2)
dotplot(EGO2)

#GO Î∞? KEGG Î∂ÑÏÑù
library(enrichR)
dbs <- listEnrichrDbs()
ERR1 <- enrichr(RECUR_LIST, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR2 <- enrichr(WHOGR_LIST, databases = c("GO_Biological_Process_2018", "KEGG_2019_Human"))
ERR1_GO <- ERR1[[1]]; ERR1_KG <- ERR1[[2]]
ERR2_GO <- ERR2[[1]]; ERR2_KG <- ERR2[[2]]
ERR1_GOx <- dplyr::filter(ERR1_GO, P.value < 0.05 & Combined.Score > 350)[,c(1,3,8:9)]
ERR2_GOx <- dplyr::filter(ERR2_GO, P.value < 0.05 & Combined.Score > 350)[,c(1,3,8:9)]
dplyr::filter(ERR1_KG, P.value < 0.05)
dplyr::filter(ERR2_KG, P.value < 0.05)
write.table(ERR1_GOx, file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
write.table(ERR2_GOx, file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)

#?ÉÅ?úÑ 20Í∞úÏùò termÎß? ?ëú?ãú?ãú
plotEnrich(ERR1[[1]], showTerms = 20)
plotEnrich(ERR1[[2]], showTerms = 20)

#?ÉÅ?úÑ 10Í∞úÏùò reactome pathway ?ôï?ù∏
library(ReactomePA)
RPA1 <- enrichPathway(gene = RECUR_LIST, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
RPA2 <- enrichPathway(gene = WHOGR_LIST, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH")
barplot(RPA1, showCategory = 10)
dotplot(RPA1, showCategory = 10)

#KEGGÎ•? ?ù¥?ö©?ï¥?Ñú Í≤ΩÎ°ú Î∂ÑÏÑù
library(KEGGREST)
EntrezLIST1 <- bitr(RECUR_LIST, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
EntrezLIST2 <- bitr(RECUR_LIST, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
KEGG1 <- enrichKEGG(gene = EntrezLIST1, organism = 'hsa', pvalueCutoff = 0.05)
KEGG2 <- enrichKEGG(gene = EntrezLIST2, organism = 'hsa', pvalueCutoff = 0.05)
barplot(KEGG1)
dotplot(KEGG1)

#GO Î∂ÑÏÑù plot ?ôï?ù∏
library(gProfileR)
gPR1 <- gprofiler2::gost(query = RECUR_LIST, organism = "hsapiens")
gPR2 <- gprofiler2::gost(query = WHOGR_LIST, organism = "hsapiens")
gPR1x <- gPR1$result
gPR2x <- gPR2$result
ggplot(RECUR_enrichx, aes(x = reorder(term_name, -p_value), y = -log10(p_value))) +  geom_bar(stat = "identity") +  coord_flip() +  labs(title = "Top 10 significant GO and KEGG terms", x = "GO/KEGG Term", y = "-log10(p-value)")


#=====================================================================================================

#?ó∞Íµ? Î™©Ìëú 2. T1 Î∞? T7 Í∞ÅÍ∞Å Î∞òÏùëÍµ∞Í≥º ÎπÑÎ∞ò?ùëÍµ∞Í∞Ñ RNA-seq Î∂ÑÏÑù
DEG_NUMBERS <- data.frame()
#RNA-seq ?ç∞?ù¥?Ñ∞??Ä p-value, FCÎ•? ?úÑ?óê Îπ? ?ç∞?ù¥?Ñ∞ ?îÑ?†à?ûÑ?óê Ï∂îÍ?Ä?ïò?ó¨ up-regulation, down regulation ?ú†?†Ñ?ûê Í∞úÏàò ??Ä?û• 
Rs03_T1_1A <- as.data.frame(RNAseq_ALL03_T1); Rs03_T1_1A$PV <- PVFC_DEG02$PV_T1_1A; Rs03_T1_1A$logPV <- -log10(PVFC_DEG02$PV_T1_1A); Rs03_T1_1A$FC <- PVFC_DEG02$FC_T1_1A; Rs03_T1_1Ax <- dplyr::filter(Rs03_T1_1A, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 1,1] <- nrow(dplyr::filter(Rs03_T1_1Ax,FC < 0)); DEG_NUMBERS[ 1,2] <- nrow(dplyr::filter(Rs03_T1_1Ax,FC > 0))
Rs03_T1_1C <- as.data.frame(RNAseq_ALL03_T1); Rs03_T1_1C$PV <- PVFC_DEG02$PV_T1_1C; Rs03_T1_1C$logPV <- -log10(PVFC_DEG02$PV_T1_1C); Rs03_T1_1C$FC <- PVFC_DEG02$FC_T1_1C; Rs03_T1_1Cx <- dplyr::filter(Rs03_T1_1C, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 2,1] <- nrow(dplyr::filter(Rs03_T1_1Cx,FC < 0)); DEG_NUMBERS[ 2,2] <- nrow(dplyr::filter(Rs03_T1_1Cx,FC > 0))
Rs03_T1_1D <- as.data.frame(RNAseq_ALL03_T1); Rs03_T1_1D$PV <- PVFC_DEG02$PV_T1_1D; Rs03_T1_1D$logPV <- -log10(PVFC_DEG02$PV_T1_1D); Rs03_T1_1D$FC <- PVFC_DEG02$FC_T1_1D; Rs03_T1_1Dx <- dplyr::filter(Rs03_T1_1D, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 3,1] <- nrow(dplyr::filter(Rs03_T1_1Dx,FC < 0)); DEG_NUMBERS[ 3,2] <- nrow(dplyr::filter(Rs03_T1_1Dx,FC > 0))
Rs03_T1_1E <- as.data.frame(RNAseq_ALL03_T1); Rs03_T1_1E$PV <- PVFC_DEG02$PV_T1_1E; Rs03_T1_1E$logPV <- -log10(PVFC_DEG02$PV_T1_1E); Rs03_T1_1E$FC <- PVFC_DEG02$FC_T1_1E; Rs03_T1_1Ex <- dplyr::filter(Rs03_T1_1E, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 4,1] <- nrow(dplyr::filter(Rs03_T1_1Ex,FC < 0)); DEG_NUMBERS[ 4,2] <- nrow(dplyr::filter(Rs03_T1_1Ex,FC > 0))
Rs03_T1_2A <- as.data.frame(RNAseq_ALL03_T1); Rs03_T1_2A$PV <- PVFC_DEG02$PV_T1_2A; Rs03_T1_2A$logPV <- -log10(PVFC_DEG02$PV_T1_2A); Rs03_T1_2A$FC <- PVFC_DEG02$FC_T1_2A; Rs03_T1_2Ax <- dplyr::filter(Rs03_T1_2A, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 5,1] <- nrow(dplyr::filter(Rs03_T1_2Ax,FC < 0)); DEG_NUMBERS[ 5,2] <- nrow(dplyr::filter(Rs03_T1_2Ax,FC > 0))
Rs03_T1_2C <- as.data.frame(RNAseq_ALL03_T1); Rs03_T1_2C$PV <- PVFC_DEG02$PV_T1_2C; Rs03_T1_2C$logPV <- -log10(PVFC_DEG02$PV_T1_2C); Rs03_T1_2C$FC <- PVFC_DEG02$FC_T1_2C; Rs03_T1_2Cx <- dplyr::filter(Rs03_T1_2C, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 6,1] <- nrow(dplyr::filter(Rs03_T1_2Cx,FC < 0)); DEG_NUMBERS[ 6,2] <- nrow(dplyr::filter(Rs03_T1_2Cx,FC > 0))
Rs03_T1_2D <- as.data.frame(RNAseq_ALL03_T1); Rs03_T1_2D$PV <- PVFC_DEG02$PV_T1_2D; Rs03_T1_2D$logPV <- -log10(PVFC_DEG02$PV_T1_2D); Rs03_T1_2D$FC <- PVFC_DEG02$FC_T1_2D; Rs03_T1_2Dx <- dplyr::filter(Rs03_T1_2D, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 7,1] <- nrow(dplyr::filter(Rs03_T1_2Dx,FC < 0)); DEG_NUMBERS[ 7,2] <- nrow(dplyr::filter(Rs03_T1_2Dx,FC > 0))
Rs03_T1_2E <- as.data.frame(RNAseq_ALL03_T1); Rs03_T1_2E$PV <- PVFC_DEG02$PV_T1_2E; Rs03_T1_2E$logPV <- -log10(PVFC_DEG02$PV_T1_2E); Rs03_T1_2E$FC <- PVFC_DEG02$FC_T1_2E; Rs03_T1_2Ex <- dplyr::filter(Rs03_T1_2E, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 8,1] <- nrow(dplyr::filter(Rs03_T1_2Ex,FC < 0)); DEG_NUMBERS[ 8,2] <- nrow(dplyr::filter(Rs03_T1_2Ex,FC > 0))
Rs03_T7_1A <- as.data.frame(RNAseq_ALL03_T7); Rs03_T7_1A$PV <- PVFC_DEG02$PV_T7_1A; Rs03_T7_1A$logPV <- -log10(PVFC_DEG02$PV_T7_1A); Rs03_T7_1A$FC <- PVFC_DEG02$FC_T7_1A; Rs03_T7_1Ax <- dplyr::filter(Rs03_T7_1A, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[ 9,1] <- nrow(dplyr::filter(Rs03_T7_1Ax,FC < 0)); DEG_NUMBERS[ 9,2] <- nrow(dplyr::filter(Rs03_T7_1Ax,FC > 0))
Rs03_T7_1C <- as.data.frame(RNAseq_ALL03_T7); Rs03_T7_1C$PV <- PVFC_DEG02$PV_T7_1C; Rs03_T7_1C$logPV <- -log10(PVFC_DEG02$PV_T7_1C); Rs03_T7_1C$FC <- PVFC_DEG02$FC_T7_1C; Rs03_T7_1Cx <- dplyr::filter(Rs03_T7_1C, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[10,1] <- nrow(dplyr::filter(Rs03_T7_1Cx,FC < 0)); DEG_NUMBERS[10,2] <- nrow(dplyr::filter(Rs03_T7_1Cx,FC > 0))
Rs03_T7_1D <- as.data.frame(RNAseq_ALL03_T7); Rs03_T7_1D$PV <- PVFC_DEG02$PV_T7_1D; Rs03_T7_1D$logPV <- -log10(PVFC_DEG02$PV_T7_1D); Rs03_T7_1D$FC <- PVFC_DEG02$FC_T7_1D; Rs03_T7_1Dx <- dplyr::filter(Rs03_T7_1D, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[11,1] <- nrow(dplyr::filter(Rs03_T7_1Dx,FC < 0)); DEG_NUMBERS[11,2] <- nrow(dplyr::filter(Rs03_T7_1Dx,FC > 0))
Rs03_T7_1E <- as.data.frame(RNAseq_ALL03_T7); Rs03_T7_1E$PV <- PVFC_DEG02$PV_T7_1E; Rs03_T7_1E$logPV <- -log10(PVFC_DEG02$PV_T7_1E); Rs03_T7_1E$FC <- PVFC_DEG02$FC_T7_1E; Rs03_T7_1Ex <- dplyr::filter(Rs03_T7_1E, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[12,1] <- nrow(dplyr::filter(Rs03_T7_1Ex,FC < 0)); DEG_NUMBERS[12,2] <- nrow(dplyr::filter(Rs03_T7_1Ex,FC > 0))
Rs03_T7_2A <- as.data.frame(RNAseq_ALL03_T7); Rs03_T7_2A$PV <- PVFC_DEG02$PV_T7_2A; Rs03_T7_2A$logPV <- -log10(PVFC_DEG02$PV_T7_2A); Rs03_T7_2A$FC <- PVFC_DEG02$FC_T7_2A; Rs03_T7_2Ax <- dplyr::filter(Rs03_T7_2A, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[13,1] <- nrow(dplyr::filter(Rs03_T7_2Ax,FC < 0)); DEG_NUMBERS[13,2] <- nrow(dplyr::filter(Rs03_T7_2Ax,FC > 0))
Rs03_T7_2C <- as.data.frame(RNAseq_ALL03_T7); Rs03_T7_2C$PV <- PVFC_DEG02$PV_T7_2C; Rs03_T7_2C$logPV <- -log10(PVFC_DEG02$PV_T7_2C); Rs03_T7_2C$FC <- PVFC_DEG02$FC_T7_2C; Rs03_T7_2Cx <- dplyr::filter(Rs03_T7_2C, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[14,1] <- nrow(dplyr::filter(Rs03_T7_2Cx,FC < 0)); DEG_NUMBERS[14,2] <- nrow(dplyr::filter(Rs03_T7_2Cx,FC > 0))
Rs03_T7_2D <- as.data.frame(RNAseq_ALL03_T7); Rs03_T7_2D$PV <- PVFC_DEG02$PV_T7_2D; Rs03_T7_2D$logPV <- -log10(PVFC_DEG02$PV_T7_2D); Rs03_T7_2D$FC <- PVFC_DEG02$FC_T7_2D; Rs03_T7_2Dx <- dplyr::filter(Rs03_T7_2D, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[15,1] <- nrow(dplyr::filter(Rs03_T7_2Dx,FC < 0)); DEG_NUMBERS[15,2] <- nrow(dplyr::filter(Rs03_T7_2Dx,FC > 0))
Rs03_T7_2E <- as.data.frame(RNAseq_ALL03_T7); Rs03_T7_2E$PV <- PVFC_DEG02$PV_T7_2E; Rs03_T7_2E$logPV <- -log10(PVFC_DEG02$PV_T7_2E); Rs03_T7_2E$FC <- PVFC_DEG02$FC_T7_2E; Rs03_T7_2Ex <- dplyr::filter(Rs03_T7_2E, PV < 10^-2 & abs(FC) > 0.12); DEG_NUMBERS[16,1] <- nrow(dplyr::filter(Rs03_T7_2Ex,FC < 0)); DEG_NUMBERS[16,2] <- nrow(dplyr::filter(Rs03_T7_2Ex,FC > 0))
colnames(DEG_NUMBERS) <- c("NE_less","NE_more")
rownames(DEG_NUMBERS) <- c("T1_1A","T1_1C","T1_1D","T1_1E","T1_2A","T1_2C","T1_2D","T1_2E","T7_1A","T7_1C","T7_1D","T7_1E","T7_2A","T7_2C","T7_2D","T7_2E")
write.table(DEG_NUMBERS, file="clipboard",sep="\t",quote=FALSE,col.names=NA)

#?ú†?ùò?ïú Ï∞®Ïù¥Î•? Î≥¥Ïù¥?äî ?ú†?†Ñ?ûê ?ïÑ?Ñ∞ÎßÅÌïú ?í§ Î∞úÌòÑ Ï¶ùÍ?Ä, Í∞êÏÜå ?ú†?†Ñ?ûê ?ôï?ù∏
Rs03_T1_1Ay <- dplyr::filter(Rs03_T1_1A, PV < 10^-2 & abs(FC) > 0.15); nrow(dplyr::filter(Rs03_T1_1Ay,FC < 0)); nrow(dplyr::filter(Rs03_T1_1Ay,FC > 0))
Rs03_T1_1Cy <- dplyr::filter(Rs03_T1_1C, PV < 10^-3 & abs(FC) > 0.20); nrow(dplyr::filter(Rs03_T1_1Cy,FC < 0)); nrow(dplyr::filter(Rs03_T1_1Cy,FC > 0))
Rs03_T1_1Dy <- dplyr::filter(Rs03_T1_1D, PV < 10^-3 & abs(FC) > 0.19); nrow(dplyr::filter(Rs03_T1_1Dy,FC < 0)); nrow(dplyr::filter(Rs03_T1_1Dy,FC > 0))
Rs03_T1_1Ey <- dplyr::filter(Rs03_T1_1E, PV < 10^-2 & abs(FC) > 0.18); nrow(dplyr::filter(Rs03_T1_1Ey,FC < 0)); nrow(dplyr::filter(Rs03_T1_1Ey,FC > 0))
Rs03_T1_2Ay <- dplyr::filter(Rs03_T1_2A, PV < 10^-3 & abs(FC) > 0.15); nrow(dplyr::filter(Rs03_T1_2Ay,FC < 0)); nrow(dplyr::filter(Rs03_T1_2Ay,FC > 0))
Rs03_T1_2Cy <- dplyr::filter(Rs03_T1_2C, PV < 10^-2 & abs(FC) > 0.18); nrow(dplyr::filter(Rs03_T1_2Cy,FC < 0)); nrow(dplyr::filter(Rs03_T1_2Cy,FC > 0))
Rs03_T1_2Dy <- dplyr::filter(Rs03_T1_2D, PV < 10^-3 & abs(FC) > 0.15); nrow(dplyr::filter(Rs03_T1_2Dy,FC < 0)); nrow(dplyr::filter(Rs03_T1_2Dy,FC > 0))
Rs03_T1_2Ey <- dplyr::filter(Rs03_T1_2E, PV < 10^-2 & abs(FC) > 0.12); nrow(dplyr::filter(Rs03_T1_2Ey,FC < 0)); nrow(dplyr::filter(Rs03_T1_2Ey,FC > 0))
Rs03_T7_1Ay <- dplyr::filter(Rs03_T7_1A, PV < 10^-2 & abs(FC) > 0.12); nrow(dplyr::filter(Rs03_T7_1Ay,FC < 0)); nrow(dplyr::filter(Rs03_T7_1Ay,FC > 0))
Rs03_T7_1Cy <- dplyr::filter(Rs03_T7_1C, PV < 10^-3 & abs(FC) > 0.15); nrow(dplyr::filter(Rs03_T7_1Cy,FC < 0)); nrow(dplyr::filter(Rs03_T7_1Cy,FC > 0))
Rs03_T7_1Dy <- dplyr::filter(Rs03_T7_1D, PV < 10^-2 & abs(FC) > 0.12); nrow(dplyr::filter(Rs03_T7_1Dy,FC < 0)); nrow(dplyr::filter(Rs03_T7_1Dy,FC > 0))
Rs03_T7_1Ey <- dplyr::filter(Rs03_T7_1E, PV < 10^-2 & abs(FC) > 0.12); nrow(dplyr::filter(Rs03_T7_1Ey,FC < 0)); nrow(dplyr::filter(Rs03_T7_1Ey,FC > 0))
Rs03_T7_2Ay <- dplyr::filter(Rs03_T7_2A, PV < 10^-2 & abs(FC) > 0.12); nrow(dplyr::filter(Rs03_T7_2Ay,FC < 0)); nrow(dplyr::filter(Rs03_T7_2Ay,FC > 0))
Rs03_T7_2Cy <- dplyr::filter(Rs03_T7_2C, PV < 10^-2 & abs(FC) > 0.12); nrow(dplyr::filter(Rs03_T7_2Cy,FC < 0)); nrow(dplyr::filter(Rs03_T7_2Cy,FC > 0))
Rs03_T7_2Dy <- dplyr::filter(Rs03_T7_2D, PV < 10^-2 & abs(FC) > 0.12); nrow(dplyr::filter(Rs03_T7_2Dy,FC < 0)); nrow(dplyr::filter(Rs03_T7_2Dy,FC > 0))
Rs03_T7_2Ey <- dplyr::filter(Rs03_T7_2E, PV < 10^-2 & abs(FC) > 0.18); nrow(dplyr::filter(Rs03_T7_2Ey,FC < 0)); nrow(dplyr::filter(Rs03_T7_2Ey,FC > 0))

#?úÑ?óê?Ñú ?ÜµÍ≥ÑÏ†Å?úºÎ°? ?ú†?ùò?ïú ?ú†?†Ñ?ûê?ì§Îß? ?Ç®Í∏∞Í≥† ENSG ?ÑòÎ≤? Ï∂îÍ?Ä
Rs03_T1_1Az <- dplyr::filter(Rs03_T1_1A, PV < 0.05); dim(Rs03_T1_1Az); Rs03_T1_1Az$ENSEMBL <- rownames(Rs03_T1_1Az); Rs03_T1_1Az <- merge(x=Rs03_T1_1Az, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T1_1Cz <- dplyr::filter(Rs03_T1_1C, PV < 0.05); dim(Rs03_T1_1Cz); Rs03_T1_1Cz$ENSEMBL <- rownames(Rs03_T1_1Cz); Rs03_T1_1Cz <- merge(x=Rs03_T1_1Cz, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T1_1Dz <- dplyr::filter(Rs03_T1_1D, PV < 0.05); dim(Rs03_T1_1Dz); Rs03_T1_1Dz$ENSEMBL <- rownames(Rs03_T1_1Dz); Rs03_T1_1Dz <- merge(x=Rs03_T1_1Dz, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T1_1Ez <- dplyr::filter(Rs03_T1_1E, PV < 0.05); dim(Rs03_T1_1Ez); Rs03_T1_1Ez$ENSEMBL <- rownames(Rs03_T1_1Ez); Rs03_T1_1Ez <- merge(x=Rs03_T1_1Ez, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T1_2Az <- dplyr::filter(Rs03_T1_2A, PV < 0.05); dim(Rs03_T1_2Az); Rs03_T1_2Az$ENSEMBL <- rownames(Rs03_T1_2Az); Rs03_T1_2Az <- merge(x=Rs03_T1_2Az, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T1_2Cz <- dplyr::filter(Rs03_T1_2C, PV < 0.05); dim(Rs03_T1_2Cz); Rs03_T1_2Cz$ENSEMBL <- rownames(Rs03_T1_2Cz); Rs03_T1_2Cz <- merge(x=Rs03_T1_2Cz, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T1_2Dz <- dplyr::filter(Rs03_T1_2D, PV < 0.05); dim(Rs03_T1_2Dz); Rs03_T1_2Dz$ENSEMBL <- rownames(Rs03_T1_2Dz); Rs03_T1_2Dz <- merge(x=Rs03_T1_2Dz, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T1_2Ez <- dplyr::filter(Rs03_T1_2E, PV < 0.05); dim(Rs03_T1_2Ez); Rs03_T1_2Ez$ENSEMBL <- rownames(Rs03_T1_2Ez); Rs03_T1_2Ez <- merge(x=Rs03_T1_2Ez, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T7_1Az <- dplyr::filter(Rs03_T7_1A, PV < 0.05); dim(Rs03_T7_1Az); Rs03_T7_1Az$ENSEMBL <- rownames(Rs03_T7_1Az); Rs03_T7_1Az <- merge(x=Rs03_T7_1Az, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T7_1Cz <- dplyr::filter(Rs03_T7_1C, PV < 0.05); dim(Rs03_T7_1Cz); Rs03_T7_1Cz$ENSEMBL <- rownames(Rs03_T7_1Cz); Rs03_T7_1Cz <- merge(x=Rs03_T7_1Cz, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T7_1Dz <- dplyr::filter(Rs03_T7_1D, PV < 0.05); dim(Rs03_T7_1Dz); Rs03_T7_1Dz$ENSEMBL <- rownames(Rs03_T7_1Dz); Rs03_T7_1Dz <- merge(x=Rs03_T7_1Dz, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T7_1Ez <- dplyr::filter(Rs03_T7_1E, PV < 0.05); dim(Rs03_T7_1Ez); Rs03_T7_1Ez$ENSEMBL <- rownames(Rs03_T7_1Ez); Rs03_T7_1Ez <- merge(x=Rs03_T7_1Ez, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T7_2Az <- dplyr::filter(Rs03_T7_2A, PV < 0.05); dim(Rs03_T7_2Az); Rs03_T7_2Az$ENSEMBL <- rownames(Rs03_T7_2Az); Rs03_T7_2Az <- merge(x=Rs03_T7_2Az, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T7_2Cz <- dplyr::filter(Rs03_T7_2C, PV < 0.05); dim(Rs03_T7_2Cz); Rs03_T7_2Cz$ENSEMBL <- rownames(Rs03_T7_2Cz); Rs03_T7_2Cz <- merge(x=Rs03_T7_2Cz, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T7_2Dz <- dplyr::filter(Rs03_T7_2D, PV < 0.05); dim(Rs03_T7_2Dz); Rs03_T7_2Dz$ENSEMBL <- rownames(Rs03_T7_2Dz); Rs03_T7_2Dz <- merge(x=Rs03_T7_2Dz, y=RNAseq_GENEs, by="ENSEMBL")
Rs03_T7_2Ez <- dplyr::filter(Rs03_T7_2E, PV < 0.05); dim(Rs03_T7_2Ez); Rs03_T7_2Ez$ENSEMBL <- rownames(Rs03_T7_2Ez); Rs03_T7_2Ez <- merge(x=Rs03_T7_2Ez, y=RNAseq_GENEs, by="ENSEMBL")

#Volcano plot 
pdf(file = "Output_20240729/02_V_Rs03_T1_1A.pdf",width=6,height=6); with(Rs03_T1_1A, plot(FC, logPV, pch=20, cex=0.8, main = "T1_1A")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.15, col = "darkgreen", lty = 3); with(subset(Rs03_T1_1A, logPV > 2 & FC >  0.15), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T1_1A, logPV > 2 & FC < -0.15), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T1_1C.pdf",width=6,height=6); with(Rs03_T1_1C, plot(FC, logPV, pch=20, cex=0.8, main = "T1_1C")); abline(h=3, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.20, col = "darkgreen", lty = 3); with(subset(Rs03_T1_1C, logPV > 3 & FC >  0.20), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T1_1C, logPV > 3 & FC < -0.20), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T1_1D.pdf",width=6,height=6); with(Rs03_T1_1D, plot(FC, logPV, pch=20, cex=0.8, main = "T1_1D")); abline(h=3, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.19, col = "darkgreen", lty = 3); with(subset(Rs03_T1_1D, logPV > 3 & FC >  0.19), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T1_1D, logPV > 3 & FC < -0.19), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T1_1E.pdf",width=6,height=6); with(Rs03_T1_1E, plot(FC, logPV, pch=20, cex=0.8, main = "T1_1E")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.18, col = "darkgreen", lty = 3); with(subset(Rs03_T1_1E, logPV > 2 & FC >  0.18), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T1_1E, logPV > 2 & FC < -0.18), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T1_2A.pdf",width=6,height=6); with(Rs03_T1_2A, plot(FC, logPV, pch=20, cex=0.8, main = "T1_2A")); abline(h=3, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.15, col = "darkgreen", lty = 3); with(subset(Rs03_T1_2A, logPV > 3 & FC >  0.15), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T1_2A, logPV > 3 & FC < -0.15), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T1_2C.pdf",width=6,height=6); with(Rs03_T1_2C, plot(FC, logPV, pch=20, cex=0.8, main = "T1_2C")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.18, col = "darkgreen", lty = 3); with(subset(Rs03_T1_2C, logPV > 2 & FC >  0.18), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T1_2C, logPV > 2 & FC < -0.18), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T1_2D.pdf",width=6,height=6); with(Rs03_T1_2D, plot(FC, logPV, pch=20, cex=0.8, main = "T1_2D")); abline(h=3, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.15, col = "darkgreen", lty = 3); with(subset(Rs03_T1_2D, logPV > 3 & FC >  0.15), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T1_2D, logPV > 3 & FC < -0.15), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T1_2E.pdf",width=6,height=6); with(Rs03_T1_2E, plot(FC, logPV, pch=20, cex=0.8, main = "T1_2E")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.12, col = "darkgreen", lty = 3); with(subset(Rs03_T1_2E, logPV > 2 & FC >  0.12), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T1_2E, logPV > 2 & FC < -0.12), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T7_1A.pdf",width=6,height=6); with(Rs03_T7_1A, plot(FC, logPV, pch=20, cex=0.8, main = "T7_1A")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.12, col = "darkgreen", lty = 3); with(subset(Rs03_T7_1A, logPV > 2 & FC >  0.12), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T7_1A, logPV > 2 & FC < -0.12), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T7_1C.pdf",width=6,height=6); with(Rs03_T7_1C, plot(FC, logPV, pch=20, cex=0.8, main = "T7_1C")); abline(h=3, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.15, col = "darkgreen", lty = 3); with(subset(Rs03_T7_1C, logPV > 3 & FC >  0.15), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T7_1C, logPV > 3 & FC < -0.15), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T7_1D.pdf",width=6,height=6); with(Rs03_T7_1D, plot(FC, logPV, pch=20, cex=0.8, main = "T7_1D")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.12, col = "darkgreen", lty = 3); with(subset(Rs03_T7_1D, logPV > 2 & FC >  0.12), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T7_1D, logPV > 2 & FC < -0.12), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T7_1E.pdf",width=6,height=6); with(Rs03_T7_1E, plot(FC, logPV, pch=20, cex=0.8, main = "T7_1E")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.12, col = "darkgreen", lty = 3); with(subset(Rs03_T7_1E, logPV > 2 & FC >  0.12), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T7_1E, logPV > 2 & FC < -0.12), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T7_2A.pdf",width=6,height=6); with(Rs03_T7_2A, plot(FC, logPV, pch=20, cex=0.8, main = "T7_2A")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.12, col = "darkgreen", lty = 3); with(subset(Rs03_T7_2A, logPV > 2 & FC >  0.12), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T7_2A, logPV > 2 & FC < -0.12), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T7_2C.pdf",width=6,height=6); with(Rs03_T7_2C, plot(FC, logPV, pch=20, cex=0.8, main = "T7_2C")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.12, col = "darkgreen", lty = 3); with(subset(Rs03_T7_2C, logPV > 2 & FC >  0.12), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T7_2C, logPV > 2 & FC < -0.12), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T7_2D.pdf",width=6,height=6); with(Rs03_T7_2D, plot(FC, logPV, pch=20, cex=0.8, main = "T7_2D")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.12, col = "darkgreen", lty = 3); with(subset(Rs03_T7_2D, logPV > 2 & FC >  0.12), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T7_2D, logPV > 2 & FC < -0.12), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
pdf(file = "Output_20240729/02_V_Rs03_T7_2E.pdf",width=6,height=6); with(Rs03_T7_2E, plot(FC, logPV, pch=20, cex=0.8, main = "T7_2E")); abline(h=2, col = "blue", lty = 3); abline(v=0.12, col = "red", lty = 3); abline(v=-0.18, col = "darkgreen", lty = 3); with(subset(Rs03_T7_2E, logPV > 2 & FC >  0.18), points(FC, logPV, pch=20, cex=0.8, col="red")); with(subset(Rs03_T7_2E, logPV > 2 & FC < -0.18), points(FC, logPV, pch=20, cex=0.8, col="darkgreen")); dev.off()
dev.off()

# ƒ√∑≥∏Ì »Æ¿Œ
colnames(Rs03_T1_2Dy)
head(Rs03_T1_2Dy)

# T1 ∞·∞˙ √ﬂ√‚
Rs03_T1_2Dy$ENSEMBL <- rownames(Rs03_T1_2Dy)
Rs03_T1_2D_out <- merge(Rs03_T1_2Dy[, c("ENSEMBL","PV","logPV","FC")], 
                        RNAseq_GENEs, by="ENSEMBL")
write.table(Rs03_T1_2D_out[, c("ENSEMBL","SYMBOL","PV","logPV","FC")],
            file="clipboard", sep="\t", quote=FALSE, row.names=FALSE)

# T7 ∞·∞˙ √ﬂ√‚
Rs03_T7_2Dy$ENSEMBL <- rownames(Rs03_T7_2Dy)
Rs03_T7_2D_out <- merge(Rs03_T7_2Dy[, c("ENSEMBL","PV","logPV","FC")],
                        RNAseq_GENEs, by="ENSEMBL")
write.table(Rs03_T7_2D_out[, c("ENSEMBL","SYMBOL","PV","logPV","FC")],
            file="clipboard", sep="\t", quote=FALSE, row.names=FALSE)
# T1 ∞·∞˙ ∆ƒ¿œ ¿˙¿Â
write.table(Rs03_T1_2D_out[, c("ENSEMBL","SYMBOL","PV","logPV","FC")],
            file="Rs03_T1_C2D_results.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# T7 ∞·∞˙ ∆ƒ¿œ ¿˙¿Â
write.table(Rs03_T7_2D_out[, c("ENSEMBL","SYMBOL","PV","logPV","FC")],
            file="Rs03_T7_C2D_results.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
# T1 πﬂ«ˆ∞™ + ±◊∑Ï ¡§∫∏
t1_expr <- as.data.frame(t(RNAseq_ALL03_T1))
t1_expr$ID <- rownames(t1_expr)
t1_expr$Group <- ifelse(PROPEL03c$C2D == "E", "E", "N")

write.csv(t1_expr, "Rs03_T1_expression.csv", row.names=FALSE)

# T7 πﬂ«ˆ∞™ + ±◊∑Ï ¡§∫∏  
t7_expr <- as.data.frame(t(RNAseq_ALL03_T7))
t7_expr$ID <- rownames(t7_expr)
t7_expr$Group <- ifelse(PROPEL03c$C2D == "E", "E", "N")

write.csv(t7_expr, "Rs03_T7_expression.csv", row.names=FALSE)
#?ã§Î•? Î≥Ä?àò??Ä ?ÉÅÍ¥Ä?Ñ± ?ôï?ù∏?
##?ç∞?ù¥?Ñ∞ ?†ÑÏ≤òÎ¶¨
PROPEL03x <- PROPEL03
PROPEL03x <- dplyr::select(PROPEL03x, -other_employ, -curr_meds, -cig_smoking, -cig_use_years); dim(PROPEL03x)
colnames(PROPEL03x)[15:16] <- c("weight","height")
PROPEL03x[,4:148] <- lapply(PROPEL03x[,4:148], as.factor); str(PROPEL03x)
PROPEL03x[,c(5,9,14:17,62:86,110:114,136:148)] <- lapply(PROPEL03x[,c(5,9,14:17,62:86,110:114,136:148)], as.numeric); str(PROPEL03x)
str(PROPEL03x[,  1: 90])
str(PROPEL03x[, 91:148])

#?ãúÍ∞ÑÎ?Ä Î≥ÑÎ°ú Íµ¨Î∂ÑÎ∂?
PROPEL03T1 <- dplyr::filter(PROPEL03x, session == 1); rownames(PROPEL03T1) <- PROPEL03T1$record_id; PROPEL03T1 <- PROPEL03T1[,-1]
PROPEL03T7 <- dplyr::filter(PROPEL03x, session == 7); rownames(PROPEL03T7) <- PROPEL03T7$record_id; PROPEL03T7 <- PROPEL03T7[,-1]
table(rownames(PROPEL03T1)==rownames(PROPEL03c))
table(rownames(PROPEL03T7)==rownames(PROPEL03c))
PROPEL03T1 <- cbind(PROPEL03c,PROPEL03T1)
PROPEL03T7 <- cbind(PROPEL03c,PROPEL03T7)


PROPEL03d <- PROPEL03c; PROPEL03d$record_id <- rownames(PROPEL03d)
PROPEL03y <- merge(x=PROPEL03d, y=PROPEL03x, by="record_id", all.y=TRUE)
rownames(PROPEL03y) <- PROPEL03y$subject_id; PROPEL03y <- PROPEL03y[,-c(1,3,8,12:13)]

COLOR <- colorRampPalette(c("#0571B0", "#FFFCCC", "#D01C8B"))(500)
annotation_colors = list(
  session  = c("1"="#94FFD8", "7"="#FF76CE"),
  C1A      = c("E"="#006769", "N"="#FF0080"),
  C1C      = c("E"="#006769", "N"="#FF0080"),
  C1D      = c("E"="#006769", "N"="#FF0080"),
  C1E      = c("E"="#006769", "N"="#FF0080"),
  C2A      = c("E"="#006769", "N"="#FF0080"),
  C2C      = c("E"="#006769", "N"="#FF0080"),
  C2D      = c("E"="#006769", "N"="#FF0080"),
  C2E      = c("E"="#006769", "N"="#FF0080"),
  HLQ4  = c("28"="#E1CCEC", "22"="#850F8D"),
  HLQ7  = c("28"="#FED8B1", "22"="#6F4E37"),
  education_level  = c("28"="#FFE8C5", "22"="#FF0000"),
  income  = c("28"="#FFFF80", "22"="#FF0080"),
  mens_his  = c("1"="#0571B0", "2"="#CA0020", "3"="#00DFA2", "4"="#F6FA70", "5"="#344955"),
  cognition_TScore = c("28"="#FFF67E", "22"="#416D19"),
  fatigue_TScore   = c("28"="#DFA878", "22"="#6C3428"),
  FC   = c( "2"="#0571B0", "-2"="#CA0020"),
  logPV= c("10"="#BABABA", "1"="#111111"),
  age  = c("28"="#E0F4FF", "22"="#39A7FF"),
  gender  = c("2"="#CA0020", "1"="#0571B0"))



pheatmap(dplyr::select(Rs03_T1_1Ay, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_1Ay$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1A, session), annotation_row = dplyr::select(Rs03_T1_1Ay, logPV, FC), cutree_rows = 4, cutree_cols = 3)



table(PROPEL03T1$mens_his)
table(PROPEL03T7$mens_his)


c(fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1A, session)
str(PROPEL03T1)
table(PROPEL03$education_level)

#?ã§Î•? Î≥Ä?àò?ì§Í≥? Î≥ëÌï©?ï©
PROPEL03y_C1A <- cbind(PROPEL03y[,1],PROPEL03y[,-1:-8]); names(PROPEL03y_C1A)[1] <- "C1A"
PROPEL03y_C1C <- cbind(PROPEL03y[,2],PROPEL03y[,-1:-8]); names(PROPEL03y_C1C)[1] <- "C1C"
PROPEL03y_C1D <- cbind(PROPEL03y[,3],PROPEL03y[,-1:-8]); names(PROPEL03y_C1D)[1] <- "C1D"
PROPEL03y_C1E <- cbind(PROPEL03y[,4],PROPEL03y[,-1:-8]); names(PROPEL03y_C1E)[1] <- "C1E"
PROPEL03y_C2A <- cbind(PROPEL03y[,5],PROPEL03y[,-1:-8]); names(PROPEL03y_C2A)[1] <- "C2A"
PROPEL03y_C2C <- cbind(PROPEL03y[,6],PROPEL03y[,-1:-8]); names(PROPEL03y_C2C)[1] <- "C2C"
PROPEL03y_C2D <- cbind(PROPEL03y[,7],PROPEL03y[,-1:-8]); names(PROPEL03y_C2D)[1] <- "C2D"
PROPEL03y_C2E <- cbind(PROPEL03y[,8],PROPEL03y[,-1:-8]); names(PROPEL03y_C2E)[1] <- "C2E"

#Íµ¨Ï≤¥?†Å?úºÎ°? ?ñ¥?ñ§Í≤ÉÏùÑ ?ïò?äîÍ±¥Ï?Ä..?
pdf("Output_20240729/03_ML_C1A.pdf", width=15,height=8); par(mfrow=c(1,2)); PROPEL03y_C1Ad1 <- rpart(C1A ~ ., data = PROPEL03y_C1A); rpart.plot(PROPEL03y_C1Ad1, digits=3, type=2, extra=1); PROPEL03y_C1Ar1 <- randomForest(C1A ~ ., PROPEL03y_C1A); varImpPlot(PROPEL03y_C1Ar1); dev.off()
pdf("Output_20240729/03_ML_C1C.pdf", width=15,height=8); par(mfrow=c(1,2)); PROPEL03y_C1Cd1 <- rpart(C1C ~ ., data = PROPEL03y_C1C); rpart.plot(PROPEL03y_C1Cd1, digits=3, type=2, extra=1); PROPEL03y_C1Cr1 <- randomForest(C1C ~ ., PROPEL03y_C1C); varImpPlot(PROPEL03y_C1Cr1); dev.off()
pdf("Output_20240729/03_ML_C1D.pdf", width=15,height=8); par(mfrow=c(1,2)); PROPEL03y_C1Dd1 <- rpart(C1D ~ ., data = PROPEL03y_C1D); rpart.plot(PROPEL03y_C1Dd1, digits=3, type=2, extra=1); PROPEL03y_C1Dr1 <- randomForest(C1D ~ ., PROPEL03y_C1D); varImpPlot(PROPEL03y_C1Dr1); dev.off()
pdf("Output_20240729/03_ML_C1E.pdf", width=15,height=8); par(mfrow=c(1,2)); PROPEL03y_C1Ed1 <- rpart(C1E ~ ., data = PROPEL03y_C1E); rpart.plot(PROPEL03y_C1Ed1, digits=3, type=2, extra=1); PROPEL03y_C1Er1 <- randomForest(C1E ~ ., PROPEL03y_C1E); varImpPlot(PROPEL03y_C1Er1); dev.off()
pdf("Output_20240729/03_ML_C2A.pdf", width=15,height=8); par(mfrow=c(1,2)); PROPEL03y_C2Ad1 <- rpart(C2A ~ ., data = PROPEL03y_C2A); rpart.plot(PROPEL03y_C2Ad1, digits=3, type=2, extra=1); PROPEL03y_C2Ar1 <- randomForest(C2A ~ ., PROPEL03y_C2A); varImpPlot(PROPEL03y_C2Ar1); dev.off()
pdf("Output_20240729/03_ML_C2C.pdf", width=15,height=8); par(mfrow=c(1,2)); PROPEL03y_C2Cd1 <- rpart(C2C ~ ., data = PROPEL03y_C2C); rpart.plot(PROPEL03y_C2Cd1, digits=3, type=2, extra=1); PROPEL03y_C2Cr1 <- randomForest(C2C ~ ., PROPEL03y_C2C); varImpPlot(PROPEL03y_C2Cr1); dev.off()
pdf("Output_20240729/03_ML_C2D.pdf", width=15,height=8); par(mfrow=c(1,2)); PROPEL03y_C2Dd1 <- rpart(C2D ~ ., data = PROPEL03y_C2D); rpart.plot(PROPEL03y_C2Dd1, digits=3, type=2, extra=1); PROPEL03y_C2Dr1 <- randomForest(C2D ~ ., PROPEL03y_C2D); varImpPlot(PROPEL03y_C2Dr1); dev.off()
pdf("Output_20240729/03_ML_C2E.pdf", width=15,height=8); par(mfrow=c(1,2)); PROPEL03y_C2Ed1 <- rpart(C2E ~ ., data = PROPEL03y_C2E); rpart.plot(PROPEL03y_C2Ed1, digits=3, type=2, extra=1); PROPEL03y_C2Er1 <- randomForest(C2E ~ ., PROPEL03y_C2E); varImpPlot(PROPEL03y_C2Er1); dev.off()

#?ú†?†Ñ?ûê ?†ïÎ≥? Î≥ëÌï©
Rs03_T1_1Ay$ENSEMBL <- rownames(Rs03_T1_1Ay); Rs03_T1_1Ay <- merge(x=Rs03_T1_1Ay, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T1_1Ay) <- Rs03_T1_1Ay$ENSEMBL; Rs03_T1_1Ay <- Rs03_T1_1Ay[,-1]
Rs03_T1_1Cy$ENSEMBL <- rownames(Rs03_T1_1Cy); Rs03_T1_1Cy <- merge(x=Rs03_T1_1Cy, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T1_1Cy) <- Rs03_T1_1Cy$ENSEMBL; Rs03_T1_1Cy <- Rs03_T1_1Cy[,-1]
Rs03_T1_1Dy$ENSEMBL <- rownames(Rs03_T1_1Dy); Rs03_T1_1Dy <- merge(x=Rs03_T1_1Dy, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T1_1Dy) <- Rs03_T1_1Dy$ENSEMBL; Rs03_T1_1Dy <- Rs03_T1_1Dy[,-1]
Rs03_T1_1Ey$ENSEMBL <- rownames(Rs03_T1_1Ey); Rs03_T1_1Ey <- merge(x=Rs03_T1_1Ey, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T1_1Ey) <- Rs03_T1_1Ey$ENSEMBL; Rs03_T1_1Ey <- Rs03_T1_1Ey[,-1]
Rs03_T1_2Ay$ENSEMBL <- rownames(Rs03_T1_2Ay); Rs03_T1_2Ay <- merge(x=Rs03_T1_2Ay, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T1_2Ay) <- Rs03_T1_2Ay$ENSEMBL; Rs03_T1_2Ay <- Rs03_T1_2Ay[,-1]
Rs03_T1_2Cy$ENSEMBL <- rownames(Rs03_T1_2Cy); Rs03_T1_2Cy <- merge(x=Rs03_T1_2Cy, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T1_2Cy) <- Rs03_T1_2Cy$ENSEMBL; Rs03_T1_2Cy <- Rs03_T1_2Cy[,-1]
Rs03_T1_2Dy$ENSEMBL <- rownames(Rs03_T1_2Dy); Rs03_T1_2Dy <- merge(x=Rs03_T1_2Dy, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T1_2Dy) <- Rs03_T1_2Dy$ENSEMBL; Rs03_T1_2Dy <- Rs03_T1_2Dy[,-1]
Rs03_T1_2Ey$ENSEMBL <- rownames(Rs03_T1_2Ey); Rs03_T1_2Ey <- merge(x=Rs03_T1_2Ey, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T1_2Ey) <- Rs03_T1_2Ey$ENSEMBL; Rs03_T1_2Ey <- Rs03_T1_2Ey[,-1]
Rs03_T7_1Ay$ENSEMBL <- rownames(Rs03_T7_1Ay); Rs03_T7_1Ay <- merge(x=Rs03_T7_1Ay, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T7_1Ay) <- Rs03_T7_1Ay$ENSEMBL; Rs03_T7_1Ay <- Rs03_T7_1Ay[,-1]
Rs03_T7_1Cy$ENSEMBL <- rownames(Rs03_T7_1Cy); Rs03_T7_1Cy <- merge(x=Rs03_T7_1Cy, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T7_1Cy) <- Rs03_T7_1Cy$ENSEMBL; Rs03_T7_1Cy <- Rs03_T7_1Cy[,-1]
Rs03_T7_1Dy$ENSEMBL <- rownames(Rs03_T7_1Dy); Rs03_T7_1Dy <- merge(x=Rs03_T7_1Dy, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T7_1Dy) <- Rs03_T7_1Dy$ENSEMBL; Rs03_T7_1Dy <- Rs03_T7_1Dy[,-1]
Rs03_T7_1Ey$ENSEMBL <- rownames(Rs03_T7_1Ey); Rs03_T7_1Ey <- merge(x=Rs03_T7_1Ey, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T7_1Ey) <- Rs03_T7_1Ey$ENSEMBL; Rs03_T7_1Ey <- Rs03_T7_1Ey[,-1]
Rs03_T7_2Ay$ENSEMBL <- rownames(Rs03_T7_2Ay); Rs03_T7_2Ay <- merge(x=Rs03_T7_2Ay, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T7_2Ay) <- Rs03_T7_2Ay$ENSEMBL; Rs03_T7_2Ay <- Rs03_T7_2Ay[,-1]
Rs03_T7_2Cy$ENSEMBL <- rownames(Rs03_T7_2Cy); Rs03_T7_2Cy <- merge(x=Rs03_T7_2Cy, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T7_2Cy) <- Rs03_T7_2Cy$ENSEMBL; Rs03_T7_2Cy <- Rs03_T7_2Cy[,-1]
Rs03_T7_2Dy$ENSEMBL <- rownames(Rs03_T7_2Dy); Rs03_T7_2Dy <- merge(x=Rs03_T7_2Dy, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T7_2Dy) <- Rs03_T7_2Dy$ENSEMBL; Rs03_T7_2Dy <- Rs03_T7_2Dy[,-1]
Rs03_T7_2Ey$ENSEMBL <- rownames(Rs03_T7_2Ey); Rs03_T7_2Ey <- merge(x=Rs03_T7_2Ey, y=RNAseq_GENEs, by="ENSEMBL"); rownames(Rs03_T7_2Ey) <- Rs03_T7_2Ey$ENSEMBL; Rs03_T7_2Ey <- Rs03_T7_2Ey[,-1]

Rs03_T1_1Ayp <- pheatmap(dplyr::select(Rs03_T1_1Ay, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_1Ay$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1A, session), annotation_row = dplyr::select(Rs03_T1_1Ay, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T1_1Cyp <- pheatmap(dplyr::select(Rs03_T1_1Cy, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_1Cy$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1C, session), annotation_row = dplyr::select(Rs03_T1_1Cy, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T1_1Dyp <- pheatmap(dplyr::select(Rs03_T1_1Dy, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_1Dy$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1D, session), annotation_row = dplyr::select(Rs03_T1_1Dy, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T1_1Eyp <- pheatmap(dplyr::select(Rs03_T1_1Ey, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_1Ey$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1E, session), annotation_row = dplyr::select(Rs03_T1_1Ey, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T1_2Ayp <- pheatmap(dplyr::select(Rs03_T1_2Ay, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_2Ay$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C2A, session), annotation_row = dplyr::select(Rs03_T1_2Ay, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T1_2Cyp <- pheatmap(dplyr::select(Rs03_T1_2Cy, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_2Cy$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C2C, session), annotation_row = dplyr::select(Rs03_T1_2Cy, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T1_2Dyp <- pheatmap(dplyr::select(Rs03_T1_2Dy, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_2Dy$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C2D, session), annotation_row = dplyr::select(Rs03_T1_2Dy, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T1_2Eyp <- pheatmap(dplyr::select(Rs03_T1_2Ey, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_2Ey$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C2E, session), annotation_row = dplyr::select(Rs03_T1_2Ey, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T7_1Ayp <- pheatmap(dplyr::select(Rs03_T7_1Ay, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T7_1Ay$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T7, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1A, session), annotation_row = dplyr::select(Rs03_T7_1Ay, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T7_1Cyp <- pheatmap(dplyr::select(Rs03_T7_1Cy, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T7_1Cy$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T7, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1C, session), annotation_row = dplyr::select(Rs03_T7_1Cy, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T7_1Dyp <- pheatmap(dplyr::select(Rs03_T7_1Dy, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T7_1Dy$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T7, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1D, session), annotation_row = dplyr::select(Rs03_T7_1Dy, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T7_1Eyp <- pheatmap(dplyr::select(Rs03_T7_1Ey, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T7_1Ey$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T7, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1E, session), annotation_row = dplyr::select(Rs03_T7_1Ey, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T7_2Ayp <- pheatmap(dplyr::select(Rs03_T7_2Ay, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T7_2Ay$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T7, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C2A, session), annotation_row = dplyr::select(Rs03_T7_2Ay, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T7_2Cyp <- pheatmap(dplyr::select(Rs03_T7_2Cy, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T7_2Cy$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T7, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C2C, session), annotation_row = dplyr::select(Rs03_T7_2Cy, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T7_2Dyp <- pheatmap(dplyr::select(Rs03_T7_2Dy, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T7_2Dy$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T7, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C2D, session), annotation_row = dplyr::select(Rs03_T7_2Dy, logPV, FC), cutree_rows = 4, cutree_cols = 3)
Rs03_T7_2Eyp <- pheatmap(dplyr::select(Rs03_T7_2Ey, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T7_2Ey$SYMBOL, annotation_legend = FALSE, annotation_col = dplyr::select(PROPEL03T7, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C2E, session), annotation_row = dplyr::select(Rs03_T7_2Ey, logPV, FC), cutree_rows = 4, cutree_cols = 3)

ggsave("Output_20240519/04_Rs03_T1_1Ayp.pdf", plot = Rs03_T1_1Ayp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T1_1Cyp.pdf", plot = Rs03_T1_1Cyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T1_1Dyp.pdf", plot = Rs03_T1_1Dyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T1_1Eyp.pdf", plot = Rs03_T1_1Eyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T1_2Ayp.pdf", plot = Rs03_T1_2Ayp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T1_2Cyp.pdf", plot = Rs03_T1_2Cyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T1_2Dyp.pdf", plot = Rs03_T1_2Dyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T1_2Eyp.pdf", plot = Rs03_T1_2Eyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T7_1Ayp.pdf", plot = Rs03_T7_1Ayp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T7_1Cyp.pdf", plot = Rs03_T7_1Cyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T7_1Dyp.pdf", plot = Rs03_T7_1Dyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T7_1Eyp.pdf", plot = Rs03_T7_1Eyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T7_2Ayp.pdf", plot = Rs03_T7_2Ayp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T7_2Cyp.pdf", plot = Rs03_T7_2Cyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T7_2Dyp.pdf", plot = Rs03_T7_2Dyp$gtable, width = 12, height = 15)
ggsave("Output_20240519/04_Rs03_T7_2Eyp.pdf", plot = Rs03_T7_2Eyp$gtable, width = 12, height = 15)

#Heatmap?óê ?Çò?ò® ?ú†?†Ñ?ûê?ì§ boxplot?úºÎ°? ?ã§?ãú ?ôï?ù∏
LEGENDp <- pheatmap(dplyr::select(Rs03_T1_1Ay, contains("PROPEL")), border_color = NA, fontsize = 15, annotation_colors=annotation_colors, color = COLOR, labels_row = Rs03_T1_1Ay$SYMBOL, annotation_col = dplyr::select(PROPEL03T1, fatigue_TScore, cognition_TScore, mens_his, income, education_level, HLQ7, HLQ4, gender, age, C1A, session), annotation_row = dplyr::select(Rs03_T1_1Ay, logPV, FC), cutree_rows = 4, cutree_cols = 3)
ggsave("Output_20240519/04_LEGENDp.pdf", plot = LEGENDp$gtable, width = 12, height = 30)

Rs03_T1_2Dy <- Rs03_T1_2Dy[order(Rs03_T1_2Dy$PV), ]
head(Rs03_T1_2Dy)

PROPEL03T1x <- PROPEL03T1; PROPEL03T1x$PID <- rownames(PROPEL03T1x); PROPEL03T1x <- dplyr::select(PROPEL03T1x, PID, starts_with("C", ignore.case = FALSE))
PROPEL03T7x <- PROPEL03T7; PROPEL03T7x$PID <- rownames(PROPEL03T7x); PROPEL03T7x <- dplyr::select(PROPEL03T7x, PID, starts_with("C", ignore.case = FALSE))

Rs03_T1_2Dz <- as.data.frame(tidyr::pivot_longer(dplyr::select(Rs03_T1_2Dy, -PV,-logPV,-FC), cols = PROPEL_01:PROPEL_34, names_to = "PID", values_to = "Expression"))
Rs03_T1_2Dz <- merge(x=PROPEL03T1x,y=Rs03_T1_2Dz, by="PID",all.y=TRUE)
Rs03_T1_C2Dp <- ggplot(Rs03_T1_2Dz, aes(x=C2D, y=Expression, fill=C2D)) + geom_boxplot() + facet_wrap(~SYMBOL) + scale_fill_manual(values = c("E"="#4DD0E1", "N"="#E91E63")) + theme_light() + stat_compare_means(method = "t.test", label = "p.signif")
ggsave("Output_20240519/11_Boxplot_Rs03_T1_C2D.pdf", plot = Rs03_T1_C2Dp, width = 20, height = 20)

Rs03_T7_2Dy <- Rs03_T7_2Dy[order(Rs03_T7_2Dy$PV), ]
head(Rs03_T7_2Dy)

Rs03_T7_2Dz <- as.data.frame(tidyr::pivot_longer(dplyr::select(Rs03_T7_2Dy, -PV,-logPV,-FC), cols = PROPEL_01:PROPEL_34, names_to = "PID", values_to = "Expression"))
Rs03_T7_2Dz <- merge(x=PROPEL03T7x,y=Rs03_T7_2Dz, by="PID",all.y=TRUE)
Rs03_T7_C2Dp <- ggplot(Rs03_T7_2Dz, aes(x=C2D, y=Expression, fill=C2D)) + geom_boxplot() + facet_wrap(~SYMBOL) + scale_fill_manual(values = c("E"="#4DD0E1", "N"="#E91E63")) + theme_light() + stat_compare_means(method = "t.test", label = "p.signif")
ggsave("Output_20240519/11_Boxplot_Rs03_T7_C2D.pdf", plot = Rs03_T7_C2Dp, width = 20, height = 20)




#*p < 0.05, **p < 0.01, ***p < 0.001


library(org.Hs.eg.db); library(clusterProfiler)
head(Rs03_T1_1Az)
#?ú†?ùò?ïú ?ú†?†Ñ?ûêÍ∞Ä ?úÑÏ™ΩÏóê ?ò§?èÑÎ°? ?†ï?†¨, ?ú†?†Ñ?ûê ?ã¨Î≥? Ï§ëÎ≥µ ?†úÍ±?
Rs03_T1_1Az <- Rs03_T1_1Az[order(Rs03_T1_1Az$PV), ]; Rs03_T1_1Az <- Rs03_T1_1Az[-which(duplicated(Rs03_T1_1Az$SYMBOL)), ]; Rs03_T1_1Azx <- clusterProfiler::bitr(Rs03_T1_1Az$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T1_1Cz <- Rs03_T1_1Cz[order(Rs03_T1_1Cz$PV), ]; Rs03_T1_1Cz <- Rs03_T1_1Cz[-which(duplicated(Rs03_T1_1Cz$SYMBOL)), ]; Rs03_T1_1Czx <- clusterProfiler::bitr(Rs03_T1_1Cz$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T1_1Dz <- Rs03_T1_1Dz[order(Rs03_T1_1Dz$PV), ]; Rs03_T1_1Dz <- Rs03_T1_1Dz[-which(duplicated(Rs03_T1_1Dz$SYMBOL)), ]; Rs03_T1_1Dzx <- clusterProfiler::bitr(Rs03_T1_1Dz$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T1_1Ez <- Rs03_T1_1Ez[order(Rs03_T1_1Ez$PV), ]; Rs03_T1_1Ez <- Rs03_T1_1Ez[-which(duplicated(Rs03_T1_1Ez$SYMBOL)), ]; Rs03_T1_1Ezx <- clusterProfiler::bitr(Rs03_T1_1Ez$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T1_2Az <- Rs03_T1_2Az[order(Rs03_T1_2Az$PV), ]; Rs03_T1_2Az <- Rs03_T1_2Az[-which(duplicated(Rs03_T1_2Az$SYMBOL)), ]; Rs03_T1_2Azx <- clusterProfiler::bitr(Rs03_T1_2Az$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T1_2Cz <- Rs03_T1_2Cz[order(Rs03_T1_2Cz$PV), ]; Rs03_T1_2Cz <- Rs03_T1_2Cz[-which(duplicated(Rs03_T1_2Cz$SYMBOL)), ]; Rs03_T1_2Czx <- clusterProfiler::bitr(Rs03_T1_2Cz$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T1_2Dz <- Rs03_T1_2Dz[order(Rs03_T1_2Dz$PV), ]; Rs03_T1_2Dz <- Rs03_T1_2Dz[-which(duplicated(Rs03_T1_2Dz$SYMBOL)), ]; Rs03_T1_2Dzx <- clusterProfiler::bitr(Rs03_T1_2Dz$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T1_2Ez <- Rs03_T1_2Ez[order(Rs03_T1_2Ez$PV), ]; Rs03_T1_2Ez <- Rs03_T1_2Ez[-which(duplicated(Rs03_T1_2Ez$SYMBOL)), ]; Rs03_T1_2Ezx <- clusterProfiler::bitr(Rs03_T1_2Ez$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T7_1Az <- Rs03_T7_1Az[order(Rs03_T7_1Az$PV), ]; Rs03_T7_1Az <- Rs03_T7_1Az[-which(duplicated(Rs03_T7_1Az$SYMBOL)), ]; Rs03_T7_1Azx <- clusterProfiler::bitr(Rs03_T7_1Az$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T7_1Cz <- Rs03_T7_1Cz[order(Rs03_T7_1Cz$PV), ]; Rs03_T7_1Cz <- Rs03_T7_1Cz[-which(duplicated(Rs03_T7_1Cz$SYMBOL)), ]; Rs03_T7_1Czx <- clusterProfiler::bitr(Rs03_T7_1Cz$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T7_1Dz <- Rs03_T7_1Dz[order(Rs03_T7_1Dz$PV), ]; Rs03_T7_1Dz <- Rs03_T7_1Dz[-which(duplicated(Rs03_T7_1Dz$SYMBOL)), ]; Rs03_T7_1Dzx <- clusterProfiler::bitr(Rs03_T7_1Dz$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T7_1Ez <- Rs03_T7_1Ez[order(Rs03_T7_1Ez$PV), ]; Rs03_T7_1Ez <- Rs03_T7_1Ez[-which(duplicated(Rs03_T7_1Ez$SYMBOL)), ]; Rs03_T7_1Ezx <- clusterProfiler::bitr(Rs03_T7_1Ez$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T7_2Az <- Rs03_T7_2Az[order(Rs03_T7_2Az$PV), ]; Rs03_T7_2Az <- Rs03_T7_2Az[-which(duplicated(Rs03_T7_2Az$SYMBOL)), ]; Rs03_T7_2Azx <- clusterProfiler::bitr(Rs03_T7_2Az$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T7_2Cz <- Rs03_T7_2Cz[order(Rs03_T7_2Cz$PV), ]; Rs03_T7_2Cz <- Rs03_T7_2Cz[-which(duplicated(Rs03_T7_2Cz$SYMBOL)), ]; Rs03_T7_2Czx <- clusterProfiler::bitr(Rs03_T7_2Cz$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T7_2Dz <- Rs03_T7_2Dz[order(Rs03_T7_2Dz$PV), ]; Rs03_T7_2Dz <- Rs03_T7_2Dz[-which(duplicated(Rs03_T7_2Dz$SYMBOL)), ]; Rs03_T7_2Dzx <- clusterProfiler::bitr(Rs03_T7_2Dz$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Rs03_T7_2Ez <- Rs03_T7_2Ez[order(Rs03_T7_2Ez$PV), ]; Rs03_T7_2Ez <- Rs03_T7_2Ez[-which(duplicated(Rs03_T7_2Ez$SYMBOL)), ]; Rs03_T7_2Ezx <- clusterProfiler::bitr(Rs03_T7_2Ez$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#?ú†?†Ñ?ûê Î∂ÑÏÑù Í≤∞Í≥ºÎ•? KEGG??Ä Í≤∞Ìï©?
Rs03_T1_1Az <- merge(x=Rs03_T1_1Az, y=Rs03_T1_1Azx, by="SYMBOL"); Rs03_T1_1Azy <- dplyr::select(Rs03_T1_1Az, FC); rownames(Rs03_T1_1Azy) <- Rs03_T1_1Az$ENTREZID; Rs03_T1_1A_KG <- enrichKEGG(gene = rownames(Rs03_T1_1Azy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T1_1Cz <- merge(x=Rs03_T1_1Cz, y=Rs03_T1_1Czx, by="SYMBOL"); Rs03_T1_1Czy <- dplyr::select(Rs03_T1_1Cz, FC); rownames(Rs03_T1_1Czy) <- Rs03_T1_1Cz$ENTREZID; Rs03_T1_1C_KG <- enrichKEGG(gene = rownames(Rs03_T1_1Czy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T1_1Dz <- merge(x=Rs03_T1_1Dz, y=Rs03_T1_1Dzx, by="SYMBOL"); Rs03_T1_1Dzy <- dplyr::select(Rs03_T1_1Dz, FC); rownames(Rs03_T1_1Dzy) <- Rs03_T1_1Dz$ENTREZID; Rs03_T1_1D_KG <- enrichKEGG(gene = rownames(Rs03_T1_1Dzy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T1_1Ez <- merge(x=Rs03_T1_1Ez, y=Rs03_T1_1Ezx, by="SYMBOL"); Rs03_T1_1Ezy <- dplyr::select(Rs03_T1_1Ez, FC); rownames(Rs03_T1_1Ezy) <- Rs03_T1_1Ez$ENTREZID; Rs03_T1_1E_KG <- enrichKEGG(gene = rownames(Rs03_T1_1Ezy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T1_2Az <- merge(x=Rs03_T1_2Az, y=Rs03_T1_2Azx, by="SYMBOL"); Rs03_T1_2Azy <- dplyr::select(Rs03_T1_2Az, FC); rownames(Rs03_T1_2Azy) <- Rs03_T1_2Az$ENTREZID; Rs03_T1_2A_KG <- enrichKEGG(gene = rownames(Rs03_T1_2Azy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T1_2Cz <- merge(x=Rs03_T1_2Cz, y=Rs03_T1_2Czx, by="SYMBOL"); Rs03_T1_2Czy <- dplyr::select(Rs03_T1_2Cz, FC); rownames(Rs03_T1_2Czy) <- Rs03_T1_2Cz$ENTREZID; Rs03_T1_2C_KG <- enrichKEGG(gene = rownames(Rs03_T1_2Czy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T1_2Dz <- merge(x=Rs03_T1_2Dz, y=Rs03_T1_2Dzx, by="SYMBOL"); Rs03_T1_2Dzy <- dplyr::select(Rs03_T1_2Dz, FC); rownames(Rs03_T1_2Dzy) <- Rs03_T1_2Dz$ENTREZID; Rs03_T1_2D_KG <- enrichKEGG(gene = rownames(Rs03_T1_2Dzy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T1_2Ez <- merge(x=Rs03_T1_2Ez, y=Rs03_T1_2Ezx, by="SYMBOL"); Rs03_T1_2Ezy <- dplyr::select(Rs03_T1_2Ez, FC); rownames(Rs03_T1_2Ezy) <- Rs03_T1_2Ez$ENTREZID; Rs03_T1_2E_KG <- enrichKEGG(gene = rownames(Rs03_T1_2Ezy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T7_1Az <- merge(x=Rs03_T7_1Az, y=Rs03_T7_1Azx, by="SYMBOL"); Rs03_T7_1Azy <- dplyr::select(Rs03_T7_1Az, FC); rownames(Rs03_T7_1Azy) <- Rs03_T7_1Az$ENTREZID; Rs03_T7_1A_KG <- enrichKEGG(gene = rownames(Rs03_T7_1Azy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T7_1Cz <- merge(x=Rs03_T7_1Cz, y=Rs03_T7_1Czx, by="SYMBOL"); Rs03_T7_1Czy <- dplyr::select(Rs03_T7_1Cz, FC); rownames(Rs03_T7_1Czy) <- Rs03_T7_1Cz$ENTREZID; Rs03_T7_1C_KG <- enrichKEGG(gene = rownames(Rs03_T7_1Czy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T7_1Dz <- merge(x=Rs03_T7_1Dz, y=Rs03_T7_1Dzx, by="SYMBOL"); Rs03_T7_1Dzy <- dplyr::select(Rs03_T7_1Dz, FC); rownames(Rs03_T7_1Dzy) <- Rs03_T7_1Dz$ENTREZID; Rs03_T7_1D_KG <- enrichKEGG(gene = rownames(Rs03_T7_1Dzy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T7_1Ez <- merge(x=Rs03_T7_1Ez, y=Rs03_T7_1Ezx, by="SYMBOL"); Rs03_T7_1Ezy <- dplyr::select(Rs03_T7_1Ez, FC); rownames(Rs03_T7_1Ezy) <- Rs03_T7_1Ez$ENTREZID; Rs03_T7_1E_KG <- enrichKEGG(gene = rownames(Rs03_T7_1Ezy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T7_2Az <- merge(x=Rs03_T7_2Az, y=Rs03_T7_2Azx, by="SYMBOL"); Rs03_T7_2Azy <- dplyr::select(Rs03_T7_2Az, FC); rownames(Rs03_T7_2Azy) <- Rs03_T7_2Az$ENTREZID; Rs03_T7_2A_KG <- enrichKEGG(gene = rownames(Rs03_T7_2Azy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T7_2Cz <- merge(x=Rs03_T7_2Cz, y=Rs03_T7_2Czx, by="SYMBOL"); Rs03_T7_2Czy <- dplyr::select(Rs03_T7_2Cz, FC); rownames(Rs03_T7_2Czy) <- Rs03_T7_2Cz$ENTREZID; Rs03_T7_2C_KG <- enrichKEGG(gene = rownames(Rs03_T7_2Czy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T7_2Dz <- merge(x=Rs03_T7_2Dz, y=Rs03_T7_2Dzx, by="SYMBOL"); Rs03_T7_2Dzy <- dplyr::select(Rs03_T7_2Dz, FC); rownames(Rs03_T7_2Dzy) <- Rs03_T7_2Dz$ENTREZID; Rs03_T7_2D_KG <- enrichKEGG(gene = rownames(Rs03_T7_2Dzy), organism = 'hsa', pvalueCutoff = 0.2)
Rs03_T7_2Ez <- merge(x=Rs03_T7_2Ez, y=Rs03_T7_2Ezx, by="SYMBOL"); Rs03_T7_2Ezy <- dplyr::select(Rs03_T7_2Ez, FC); rownames(Rs03_T7_2Ezy) <- Rs03_T7_2Ez$ENTREZID; Rs03_T7_2E_KG <- enrichKEGG(gene = rownames(Rs03_T7_2Ezy), organism = 'hsa', pvalueCutoff = 0.2)

pdf(file = "Output_20240729/05_KEGG_Rs03_T1_1A.pdf",width=6,height=6); dotplot(Rs03_T1_1A_KG, showCategory=10, title = "PROPEL: T1_1A"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T1_1C.pdf",width=6,height=6); dotplot(Rs03_T1_1C_KG, showCategory=10, title = "PROPEL: T1_1C"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T1_1D.pdf",width=6,height=6); dotplot(Rs03_T1_1D_KG, showCategory=10, title = "PROPEL: T1_1D"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T1_1E.pdf",width=6,height=6); dotplot(Rs03_T1_1E_KG, showCategory=10, title = "PROPEL: T1_1E"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T1_2A.pdf",width=6,height=6); dotplot(Rs03_T1_2A_KG, showCategory=10, title = "PROPEL: T1_2A"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T1_2C.pdf",width=6,height=6); dotplot(Rs03_T1_2C_KG, showCategory=10, title = "PROPEL: T1_2C"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T1_2D.pdf",width=6,height=6); dotplot(Rs03_T1_2D_KG, showCategory=10, title = "PROPEL: T1_2D"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T1_2E.pdf",width=6,height=6); dotplot(Rs03_T1_2E_KG, showCategory=10, title = "PROPEL: T1_2E"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T7_1A.pdf",width=6,height=6); dotplot(Rs03_T7_1A_KG, showCategory=10, title = "PROPEL: T7_1A"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T7_1C.pdf",width=6,height=6); dotplot(Rs03_T7_1C_KG, showCategory=10, title = "PROPEL: T7_1C"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T7_1D.pdf",width=6,height=6); dotplot(Rs03_T7_1D_KG, showCategory=10, title = "PROPEL: T7_1D"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T7_1E.pdf",width=6,height=6); dotplot(Rs03_T7_1E_KG, showCategory=10, title = "PROPEL: T7_1E"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T7_2A.pdf",width=6,height=6); dotplot(Rs03_T7_2A_KG, showCategory=10, title = "PROPEL: T7_2A"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T7_2C.pdf",width=6,height=6); dotplot(Rs03_T7_2C_KG, showCategory=10, title = "PROPEL: T7_2C"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T7_2D.pdf",width=6,height=6); dotplot(Rs03_T7_2D_KG, showCategory=10, title = "PROPEL: T7_2D"); dev.off()
pdf(file = "Output_20240729/05_KEGG_Rs03_T7_2E.pdf",width=6,height=6); dotplot(Rs03_T7_2E_KG, showCategory=10, title = "PROPEL: T7_2E"); dev.off()

BiocManager::install("pathview")
library(pathview)
setwd("C:/Users/KU/OneDrive/Î∞îÌÉï ?ôîÎ©?/?Ç®?àò/MSN/?ó∞Íµ? ?åå?ùº/US PROPEL_pilot RNA seq/US.PROPEL/03_Rrunning/Output_20240519_KEGG")

summary(Rs03_T1_1Azy)
hist(Rs03_T1_1Azy$FC)
hist(Rs03_T1_1Czy$FC)
hist(Rs03_T1_1Dzy$FC)
hist(Rs03_T1_1Ezy$FC)
hist(Rs03_T1_2Azy$FC)
hist(Rs03_T1_2Czy$FC)
hist(Rs03_T1_2Dzy$FC)
hist(Rs03_T1_2Ezy$FC)
hist(Rs03_T7_1Azy$FC)
hist(Rs03_T7_1Czy$FC)
hist(Rs03_T7_1Dzy$FC)
hist(Rs03_T7_1Ezy$FC)
hist(Rs03_T7_2Azy$FC)
hist(Rs03_T7_2Czy$FC)
hist(Rs03_T7_2Dzy$FC)
hist(Rs03_T7_2Ezy$FC)



#	hsa05014; Amyotrophic lateral sclerosis
hsa05332 <- pathview(gene.data  = Rs03_T1_1Azy,
                     pathway.id = "hsa05332",
                     species    = "hsa",
                     out.suffix = "T1_1A", 
                     limit      = list(gene=0.5, cpd=1))

head(Rs03_T7_1A_KG)
hsa05310 <- pathview(gene.data  = Rs03_T1_1Azy,
                     pathway.id = "hsa05310",
                     species    = "hsa",
                     out.suffix = "T7_1A", 
                     limit      = list(gene=0.5, cpd=1))
head(Rs03_T7_1C_KG)
hsa05110 <- pathview(gene.data  = Rs03_T1_1Azy,
                     pathway.id = "hsa05110",
                     species    = "hsa",
                     out.suffix = "T7_1C", 
                     limit      = list(gene=0.5, cpd=1))
head(Rs03_T7_1D_KG)
hsa04012 <- pathview(gene.data  = Rs03_T1_1Azy,
                     pathway.id = "hsa04012",
                     species    = "hsa",
                     out.suffix = "T7_1D", 
                     limit      = list(gene=0.5, cpd=1))
head(Rs03_T7_1E_KG)
hsa04672 <- pathview(gene.data  = Rs03_T1_1Azy,
                     pathway.id = "hsa04672",
                     species    = "hsa",
                     out.suffix = "T7_1E", 
                     limit      = list(gene=0.5, cpd=1))
head(Rs03_T7_2A_KG)
hsa05330 <- pathview(gene.data  = Rs03_T1_1Azy,
                     pathway.id = "hsa05330",
                     species    = "hsa",
                     out.suffix = "T7_2A", 
                     limit      = list(gene=0.5, cpd=1))
head(Rs03_T7_2D_KG)
hsa05310 <- pathview(gene.data  = Rs03_T1_1Azy,
                     pathway.id = "hsa05310",
                     species    = "hsa",
                     out.suffix = "T7_2D", 
                     limit      = list(gene=0.5, cpd=1))
head(Rs03_T7_2E_KG)
hsa04662 <- pathview(gene.data  = Rs03_T1_1Azy,
                     pathway.id = "hsa04662",
                     species    = "hsa",
                     out.suffix = "T7_2E", 
                     limit      = list(gene=0.5, cpd=1))
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(DOSE)

PVFC_DEG02 <- PVFC_DEG02; PVFC_DEG02$ENSEMBL <- rownames(PVFC_DEG02)
PVFC_DEG02 <- merge(x=PVFC_DEG02, y=RNAseq_GENEs, by="ENSEMBL")
PVFC_DEG02 <- dplyr::filter(PVFC_DEG02, SYMBOL != "")
head(PVFC_DEG02); dim(PVFC_DEG02)
PVFC_DEG02x <- clusterProfiler::bitr(PVFC_DEG02$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(PVFC_DEG02x); dim(PVFC_DEG02x)
PVFC_DEG02 <- merge(x=PVFC_DEG02, y=PVFC_DEG02x, by="SYMBOL")
table(duplicated(PVFC_DEG02$ENSEMBL))
head(PVFC_DEG02)

#Time1?óê?Ñú ?ñ¥?ñ§ ÏßëÎã® Î≥? ÎπÑÍµêÎ•? ?ïú Í≤ÉÏù∏ÏßÄ?
PVFC_DEG02_T1G_T1P <- dplyr::select(PVFC_DEG02, contains("T1G_T1P"), ENTREZID)
PVFC_DEG02_T1G_T1P <- PVFC_DEG02_T1G_T1P[order(PVFC_DEG02_T1G_T1P$PV_T1G_T1P),]
PVFC_DEG02_T1G_T1P <- PVFC_DEG02_T1G_T1P[-which(duplicated(PVFC_DEG02_T1G_T1P$ENTREZID)), ]
PVFC_DEG02_T1G_T1P <- dplyr::filter(PVFC_DEG02_T1G_T1P, PV_T1G_T1P < 0.05)
head(PVFC_DEG02_T1G_T1P); dim(PVFC_DEG02_T1G_T1P)
PVFC_DEG02_T1G_T1P_kegg <- enrichKEGG(gene = PVFC_DEG02_T1G_T1P$ENTREZID, organism = 'hsa', pvalueCutoff = 0.2)
dotplot(PVFC_DEG02_T1G_T1P_kegg, showCategory=10, title = "KEGG: T1G and T1P DEGs")


#Time 7?óê?Ñú ?ñ¥?ñ§ ÏßëÎã® Î≥? ÎπÑÍµêÎ•??
PVFC_DEG02_T7G_T7P <- dplyr::select(PVFC_DEG02, contains("T7G_T7P"), ENTREZID)
PVFC_DEG02_T7G_T7P <- PVFC_DEG02_T7G_T7P[order(PVFC_DEG02_T7G_T7P$PV_T7G_T7P),]
head(PVFC_DEG02_T7G_T7P); dim(PVFC_DEG02_T7G_T7P)
PVFC_DEG02_T7G_T7P <- PVFC_DEG02_T7G_T7P[-which(duplicated(PVFC_DEG02_T7G_T7P$ENTREZID)), ]
PVFC_DEG02_T7G_T7P <- dplyr::filter(PVFC_DEG02_T7G_T7P, PV_T7G_T7P < 0.05)
head(PVFC_DEG02_T7G_T7P); dim(PVFC_DEG02_T7G_T7P)
colnames(PVFC_DEG02)  # T7G_T7P Í¥Ä?†® ?ó¥?ù¥ ?ûà?äîÏßÄ ?ôï?ù∏
str(PVFC_DEG02_T7G_T7P)
class(PVFC_DEG02_T7G_T7P)
PVFC_DEG02_T7G_T7P <- data.frame(ENTREZID = PVFC_DEG02_T7G_T7P)
str(PVFC_DEG02_T7G_T7P)


PVFC_DEG02_T7G_T7P_kegg <- enrichKEGG(gene = PVFC_DEG02_T7G_T7P$ENTREZID, organism = 'hsa', pvalueCutoff = 0.2)
dotplot(PVFC_DEG02_T7G_T7P_kegg, showCategory=10, title = "KEGG: T7G and T7P DEGs")



