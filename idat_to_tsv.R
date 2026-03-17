#For 450k (total 20mins)
setwd("D:/Met")
Sys.setenv(TMPDIR = "D:/TEMP")
library(data.table)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(FlowSorted.Blood.450k)

ASAS0400 <- read.metharray.exp(base = "ASAS0400", extended = TRUE) #5m
ASAS0400e <- getNBeads(ASAS0400)
data.table::fwrite(ASAS0400e, file="ASAS0400e.txt", sep="\t", quote=FALSE, row.names=TRUE)
ASAS0400p <- detectionP(ASAS0400)
data.table::fwrite(ASAS0400p, file="ASAS0400p.txt", sep="\t", quote=FALSE, row.names=TRUE)
ASAS0400c <- estimateCellCounts(ASAS0400)
data.table::fwrite(ASAS0400c, file="ASAS0400c.txt", sep="\t", quote=FALSE, row.names=TRUE)

ASAS0400n <- preprocessNoob(ASAS0400)
ASAS0400b <- getBeta(ASAS0400n)
data.table::fwrite(ASAS0400b, file="ASAS0400b.txt", sep="\t", quote=FALSE, row.names=TRUE)
ASAS0400m <- getM(ASAS0400n)
data.table::fwrite(ASAS0400m, file="ASAS0400m.txt", sep="\t", quote=FALSE, row.names=TRUE)
rm(ASAS0400); rm(ASAS0400e); rm(ASAS0400p); rm(ASAS0400c); rm(ASAS0400n); rm(ASAS0400b); rm(ASAS0400m)
gc()

#For EPIC (850k)(1528 4hours, 822 90mins)
setwd("D:/Met")
Sys.setenv(TMPDIR = "D:/TEMP")
library(data.table)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(FlowSorted.Blood.EPIC)

ASAS1528 <- read.metharray.exp(base = "ASAS1528", extended = TRUE) #50mins
ASAS1528e <- getNBeads(ASAS1528) #2mins
data.table::fwrite(ASAS1528e, file="ASAS1528e.txt", sep="\t", quote=FALSE, row.names=TRUE)
ASAS1528p <- detectionP(ASAS1528) #160mins
data.table::fwrite(ASAS1528p, file="ASAS1528p.txt", sep="\t", quote=FALSE, row.names=TRUE)
ASAS1528c <- estimateCellCounts(ASAS1528)
data.table::fwrite(ASAS1528c, file="ASAS1528c.txt", sep="\t", quote=FALSE, row.names=TRUE)
ASAS1528n <- preprocessNoob(ASAS1528)
ASAS1528b <- getBeta(ASAS1528n)
data.table::fwrite(ASAS1528b, file="ASAS1528b.txt", sep="\t", quote=FALSE, row.names=TRUE)
ASAS1528m <- getM(ASAS1528n)
data.table::fwrite(ASAS1528m, file="ASAS1528m.txt", sep="\t", quote=FALSE, row.names=TRUE)
rm(ASAS1528); rm(ASAS1528e); rm(ASAS1528p); rm(ASAS1528c); rm(ASAS1528n); rm(ASAS1528b); rm(ASAS1528m)
gc()

CITY0822 <- read.metharray.exp(base = "CITY0822", extended = TRUE, force = TRUE)
CITY0822e <- getNBeads(CITY0822)
data.table::fwrite(CITY0822e, file="CITY0822e.txt", sep="\t", quote=FALSE, row.names=TRUE)
CITY0822p <- detectionP(CITY0822)
data.table::fwrite(CITY0822p, file="CITY0822p.txt", sep="\t", quote=FALSE, row.names=TRUE)
CITY0822c <- estimateCellCounts(CITY0822)
data.table::fwrite(CITY0822c, file="CITY0822c.txt", sep="\t", quote=FALSE, row.names=TRUE)
CITY0822n <- preprocessNoob(CITY0822)
CITY0822b <- getBeta(CITY0822n)
data.table::fwrite(CITY0822b, file="CITY0822b.txt", sep="\t", quote=FALSE, row.names=TRUE)
CITY0822m <- getM(CITY0822n)
data.table::fwrite(CITY0822m, file="CITY0822m.txt", sep="\t", quote=FALSE, row.names=TRUE)
rm(CITY0822); rm(CITY0822e); rm(CITY0822p); rm(CITY0822c); rm(CITY0822n); rm(CITY0822b); rm(CITY0822m)
gc()
#############################################################################################################################################################################################################################
#Processing and batch effect: 1 hour
setwd("D:/Met")
Sys.setenv(TMPDIR = "D:/TEMP")
rm(list=ls())
gc()
set.seed(312)
library(minfi)
library(IlluminaHumanMethylation450kmanifest); library(FlowSorted.Blood.450k)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest); library(FlowSorted.Blood.EPIC)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19); library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno_450k <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
anno_EPIC <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))

library(dplyr); library(data.table); library(pheatmap); library(readxl)
library(ggplot2); library(ggvenn); library(gridExtra); library(viridis); library(ggpubr); library(grid)

ASAS0400b <- as.data.frame(fread("ASAS0400b.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS0400b) <- ASAS0400b$V1; ASAS0400b$V1 <- NULL
ASAS0400m <- as.data.frame(fread("ASAS0400m.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS0400m) <- ASAS0400m$V1; ASAS0400m$V1 <- NULL
ASAS1528b <- as.data.frame(fread("ASAS1528b.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS1528b) <- ASAS1528b$V1; ASAS1528b$V1 <- NULL
ASAS1528m <- as.data.frame(fread("ASAS1528m.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS1528m) <- ASAS1528m$V1; ASAS1528m$V1 <- NULL
CITY0822b <- as.data.frame(fread("CITY0822b.txt.gz",  header = TRUE, sep = "\t")); rownames(CITY0822b) <- CITY0822b$V1; CITY0822b$V1 <- NULL
CITY0822m <- as.data.frame(fread("CITY0822m.txt.gz",  header = TRUE, sep = "\t")); rownames(CITY0822m) <- CITY0822m$V1; CITY0822m$V1 <- NULL

ASAS0400p <- as.data.frame(fread("ASAS0400p.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS0400p) <- ASAS0400p$V1; ASAS0400p$V1 <- NULL
ASAS1528p <- as.data.frame(fread("ASAS1528p.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS1528p) <- ASAS1528p$V1; ASAS1528p$V1 <- NULL
CITY0822p <- as.data.frame(fread("CITY0822p.txt.gz",  header = TRUE, sep = "\t")); rownames(CITY0822p) <- CITY0822p$V1; CITY0822p$V1 <- NULL

ASAS0400e <- as.data.frame(fread("ASAS0400e.txt.gz",  header = TRUE, sep = "\t"))
ASAS1528e <- as.data.frame(fread("ASAS1528e.txt.gz",  header = TRUE, sep = "\t"))
CITY0822e <- as.data.frame(fread("CITY0822e.txt.gz",  header = TRUE, sep = "\t"))
ASAS0400e$V1 <- paste0("b",ASAS0400e$V1); rownames(ASAS0400e) <- ASAS0400e$V1; ASAS0400e$V1 <- NULL; ASAS0400e <- as.data.frame(t(ASAS0400e))
ASAS1528e$V1 <- paste0("b",ASAS1528e$V1); rownames(ASAS1528e) <- ASAS1528e$V1; ASAS1528e$V1 <- NULL; ASAS1528e <- as.data.frame(t(ASAS1528e))
CITY0822e$V1 <- paste0("b",CITY0822e$V1); rownames(CITY0822e) <- CITY0822e$V1; CITY0822e$V1 <- NULL; CITY0822e <- as.data.frame(t(CITY0822e))

ASAS0400c <- as.data.frame(fread("ASAS0400c.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS0400c) <- ASAS0400c$V1; ASAS0400c$V1 <- NULL
ASAS1528c <- as.data.frame(fread("ASAS1528c.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS1528c) <- ASAS1528c$V1; ASAS1528c$V1 <- NULL
CITY0822c <- as.data.frame(fread("CITY0822c.txt.gz",  header = TRUE, sep = "\t")); rownames(CITY0822c) <- CITY0822c$V1; CITY0822c$V1 <- NULL

ASAS0400b[1:5,1:20]; dim(ASAS0400b)
ASAS0400m[1:5,1:20]; dim(ASAS0400m)
ASAS0400e[1:5,1:20]; dim(ASAS0400e)
ASAS0400p[1:5,1:20]; dim(ASAS0400p)
ASAS1528b[1:5,1:20]; dim(ASAS1528b)
ASAS1528m[1:5,1:20]; dim(ASAS1528m)
ASAS1528e[1:5,1:20]; dim(ASAS1528e)
ASAS1528p[1:5,1:20]; dim(ASAS1528p)
CITY0822b[1:5,1:20]; dim(CITY0822b)
CITY0822m[1:5,1:20]; dim(CITY0822m)
CITY0822e[1:5,1:20]; dim(CITY0822e)
CITY0822p[1:5,1:20]; dim(CITY0822p)
head(ASAS0400c); dim(ASAS0400c)
head(ASAS1528c); dim(ASAS1528c)
head(CITY0822c); dim(CITY0822c)

ASAS0400_nCOL <- ncol(ASAS0400p)
ASAS1528_nCOL <- ncol(ASAS1528p)
CITY0822_nCOL <- ncol(CITY0822p)
ASAS0400p$MeanPV <- rowMeans(ASAS0400p[,1:ASAS0400_nCOL])
ASAS1528p$MeanPV <- rowMeans(ASAS1528p[,1:ASAS1528_nCOL])
CITY0822p$MeanPV <- rowMeans(CITY0822p[,1:CITY0822_nCOL])
ASAS0400p$n001PV <- rowSums(ASAS0400p[,1:ASAS0400_nCOL] > 0.01)
ASAS1528p$n001PV <- rowSums(ASAS1528p[,1:ASAS1528_nCOL] > 0.01)
CITY0822p$n001PV <- rowSums(CITY0822p[,1:CITY0822_nCOL] > 0.01)
ASAS0400p$cg <- substr(rownames(ASAS0400p),1,2); table(ASAS0400p$cg)
ASAS1528p$cg <- substr(rownames(ASAS1528p),1,2); table(ASAS1528p$cg)
CITY0822p$cg <- substr(rownames(CITY0822p),1,2); table(CITY0822p$cg)

#https://github.com/sirselim/illumina450k_filtering
blacklist_450k <- as.character(read.csv("https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors = FALSE)[, 1])
epic_cross_probes <- as.character(read.csv("https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/EPIC/13059_2016_1066_MOESM1_ESM.csv", stringsAsFactors = FALSE)[, 1])
epic_snp_probes <- as.character(read.csv("https://raw.githubusercontent.com/sirselim/illumina450k_filtering/master/EPIC/13059_2016_1066_MOESM4_ESM.csv", stringsAsFactors = FALSE)[, 1])
blacklist_EPIC <- unique(c(epic_cross_probes, epic_snp_probes))

ASAS0400px <- dplyr::filter(ASAS0400p, MeanPV < 0.01 & n001PV < ASAS0400_nCOL*0.01 & cg == "cg")[,1:ASAS0400_nCOL]
ASAS1528px <- dplyr::filter(ASAS1528p, MeanPV < 0.01 & n001PV < ASAS1528_nCOL*0.01 & cg == "cg")[,1:ASAS1528_nCOL]
CITY0822px <- dplyr::filter(CITY0822p, MeanPV < 0.01 & n001PV < CITY0822_nCOL*0.01 & cg == "cg")[,1:CITY0822_nCOL]
ASAS0400px <- ASAS0400px[!(rownames(ASAS0400px) %in% blacklist_450k), ]; dim(ASAS0400p); dim(ASAS0400px) #From 485512 to 451245
ASAS1528px <- ASAS1528px[!(rownames(ASAS1528px) %in% blacklist_EPIC), ]; dim(ASAS1528p); dim(ASAS1528px) #From 866091 to 805110
CITY0822px <- CITY0822px[!(rownames(CITY0822px) %in% blacklist_EPIC), ]; dim(CITY0822p); dim(CITY0822px) #From 865859 to 803393

ASAS0400bx <- ASAS0400b[(rownames(ASAS0400b) %in% rownames(ASAS0400px)), ]; table(nrow(ASAS0400bx)==nrow(ASAS0400px))
ASAS0400mx <- ASAS0400m[(rownames(ASAS0400m) %in% rownames(ASAS0400px)), ]; table(nrow(ASAS0400mx)==nrow(ASAS0400px))
ASAS1528bx <- ASAS1528b[(rownames(ASAS1528b) %in% rownames(ASAS1528px)), ]; table(nrow(ASAS1528bx)==nrow(ASAS1528px))
ASAS1528mx <- ASAS1528m[(rownames(ASAS1528m) %in% rownames(ASAS1528px)), ]; table(nrow(ASAS1528mx)==nrow(ASAS1528px))
CITY0822bx <- CITY0822b[(rownames(CITY0822b) %in% rownames(CITY0822px)), ]; table(nrow(CITY0822bx)==nrow(CITY0822px))
CITY0822mx <- CITY0822m[(rownames(CITY0822m) %in% rownames(CITY0822px)), ]; table(nrow(CITY0822mx)==nrow(CITY0822px))
ASAS0400bx[1:5,1:20]; dim(ASAS0400bx)
ASAS0400mx[1:5,1:20]; dim(ASAS0400mx)
ASAS1528bx[1:5,1:20]; dim(ASAS1528bx)
ASAS1528mx[1:5,1:20]; dim(ASAS1528mx)
CITY0822bx[1:5,1:20]; dim(CITY0822bx)
CITY0822mx[1:5,1:20]; dim(CITY0822mx)
rm(ASAS0400b); rm(ASAS0400m); rm(ASAS0400e); rm(ASAS0400p)
rm(ASAS1528b); rm(ASAS1528m); rm(ASAS1528e); rm(ASAS1528p)
rm(CITY0822b); rm(CITY0822m); rm(CITY0822e); rm(CITY0822p)
gc()

data.table::fwrite(ASAS0400bx, file="ASAS0400bx.txt", sep="\t", quote=FALSE, row.names=TRUE)
data.table::fwrite(ASAS0400mx, file="ASAS0400mx.txt", sep="\t", quote=FALSE, row.names=TRUE)
data.table::fwrite(ASAS1528bx, file="ASAS1528bx.txt", sep="\t", quote=FALSE, row.names=TRUE)
data.table::fwrite(ASAS1528mx, file="ASAS1528mx.txt", sep="\t", quote=FALSE, row.names=TRUE)
data.table::fwrite(CITY0822bx, file="CITY0822bx.txt", sep="\t", quote=FALSE, row.names=TRUE)
data.table::fwrite(CITY0822mx, file="CITY0822mx.txt", sep="\t", quote=FALSE, row.names=TRUE)

ASAS01to10 <- as.data.frame(fread("G:/내 드라이브/Rrunning/KoGES/Epid/ASAS01to10.txt.gz",  header = TRUE, sep = "\t"))
CITY01to02 <- as.data.frame(fread("G:/내 드라이브/Rrunning/KoGES/Epid/CITY01to02.txt.gz",  header = TRUE, sep = "\t"))
ASAS0400sx <- dplyr::filter(dplyr::select(ASAS01to10, DIST_ID,AS01_SEX), DIST_ID %in% colnames(ASAS0400bx)); dim(ASAS0400sx)
ASAS1528sx <- dplyr::filter(dplyr::select(ASAS01to10, DIST_ID,AS01_SEX), DIST_ID %in% colnames(ASAS1528bx)); dim(ASAS1528sx)
CITY0822sx <- dplyr::filter(dplyr::select(CITY01to02, DIST_ID,CT01_SEX), DIST_ID %in% colnames(CITY0822bx)); dim(CITY0822sx)

ASAS0400bxM <- dplyr::select(ASAS0400bx, dplyr::filter(ASAS0400sx, AS01_SEX==1)$DIST_ID); ASAS0400bxF <- dplyr::select(ASAS0400bx, dplyr::filter(ASAS0400sx, AS01_SEX==2)$DIST_ID)
ASAS0400mxM <- dplyr::select(ASAS0400mx, dplyr::filter(ASAS0400sx, AS01_SEX==1)$DIST_ID); ASAS0400mxF <- dplyr::select(ASAS0400mx, dplyr::filter(ASAS0400sx, AS01_SEX==2)$DIST_ID)
ASAS1528bxM <- dplyr::select(ASAS1528bx, dplyr::filter(ASAS1528sx, AS01_SEX==1)$DIST_ID); ASAS1528bxF <- dplyr::select(ASAS1528bx, dplyr::filter(ASAS1528sx, AS01_SEX==2)$DIST_ID)
ASAS1528mxM <- dplyr::select(ASAS1528mx, dplyr::filter(ASAS1528sx, AS01_SEX==1)$DIST_ID); ASAS1528mxF <- dplyr::select(ASAS1528mx, dplyr::filter(ASAS1528sx, AS01_SEX==2)$DIST_ID)
CITY0822bxM <- dplyr::select(CITY0822bx, dplyr::filter(CITY0822sx, CT01_SEX==1)$DIST_ID); CITY0822bxF <- dplyr::select(CITY0822bx, dplyr::filter(CITY0822sx, CT01_SEX==2)$DIST_ID)
CITY0822mxM <- dplyr::select(CITY0822mx, dplyr::filter(CITY0822sx, CT01_SEX==1)$DIST_ID); CITY0822mxF <- dplyr::select(CITY0822mx, dplyr::filter(CITY0822sx, CT01_SEX==2)$DIST_ID)
ASAS0400bxM_VAR <- apply(ASAS0400bxM, 1, var); ASAS0400bxM_VARtop10p <- ASAS0400bxM[names(sort(ASAS0400bxM_VAR, decreasing = TRUE))[1:round(nrow(ASAS0400bxM)*0.1,0)], ]; ASAS0400bxM_PCA <- prcomp(t(ASAS0400bxM_VARtop10p), scale. = TRUE); ASAS0400bxM_PCAdf <- data.frame(PC1 = ASAS0400bxM_PCA$x[, 1], PC2 = ASAS0400bxM_PCA$x[, 2])
ASAS0400mxM_VAR <- apply(ASAS0400mxM, 1, var); ASAS0400mxM_VARtop10p <- ASAS0400mxM[names(sort(ASAS0400mxM_VAR, decreasing = TRUE))[1:round(nrow(ASAS0400mxM)*0.1,0)], ]; ASAS0400mxM_PCA <- prcomp(t(ASAS0400mxM_VARtop10p), scale. = TRUE); ASAS0400mxM_PCAdf <- data.frame(PC1 = ASAS0400mxM_PCA$x[, 1], PC2 = ASAS0400mxM_PCA$x[, 2])
ASAS1528bxM_VAR <- apply(ASAS1528bxM, 1, var); ASAS1528bxM_VARtop10p <- ASAS1528bxM[names(sort(ASAS1528bxM_VAR, decreasing = TRUE))[1:round(nrow(ASAS1528bxM)*0.1,0)], ]; ASAS1528bxM_PCA <- prcomp(t(ASAS1528bxM_VARtop10p), scale. = TRUE); ASAS1528bxM_PCAdf <- data.frame(PC1 = ASAS1528bxM_PCA$x[, 1], PC2 = ASAS1528bxM_PCA$x[, 2])
ASAS1528mxM_VAR <- apply(ASAS1528mxM, 1, var); ASAS1528mxM_VARtop10p <- ASAS1528mxM[names(sort(ASAS1528mxM_VAR, decreasing = TRUE))[1:round(nrow(ASAS1528mxM)*0.1,0)], ]; ASAS1528mxM_PCA <- prcomp(t(ASAS1528mxM_VARtop10p), scale. = TRUE); ASAS1528mxM_PCAdf <- data.frame(PC1 = ASAS1528mxM_PCA$x[, 1], PC2 = ASAS1528mxM_PCA$x[, 2])
CITY0822bxM_VAR <- apply(CITY0822bxM, 1, var); CITY0822bxM_VARtop10p <- CITY0822bxM[names(sort(CITY0822bxM_VAR, decreasing = TRUE))[1:round(nrow(CITY0822bxM)*0.1,0)], ]; CITY0822bxM_PCA <- prcomp(t(CITY0822bxM_VARtop10p), scale. = TRUE); CITY0822bxM_PCAdf <- data.frame(PC1 = CITY0822bxM_PCA$x[, 1], PC2 = CITY0822bxM_PCA$x[, 2])
CITY0822mxM_VAR <- apply(CITY0822mxM, 1, var); CITY0822mxM_VARtop10p <- CITY0822mxM[names(sort(CITY0822mxM_VAR, decreasing = TRUE))[1:round(nrow(CITY0822mxM)*0.1,0)], ]; CITY0822mxM_PCA <- prcomp(t(CITY0822mxM_VARtop10p), scale. = TRUE); CITY0822mxM_PCAdf <- data.frame(PC1 = CITY0822mxM_PCA$x[, 1], PC2 = CITY0822mxM_PCA$x[, 2])
ASAS0400bxF_VAR <- apply(ASAS0400bxF, 1, var); ASAS0400bxF_VARtop10p <- ASAS0400bxF[names(sort(ASAS0400bxF_VAR, decreasing = TRUE))[1:round(nrow(ASAS0400bxF)*0.1,0)], ]; ASAS0400bxF_PCA <- prcomp(t(ASAS0400bxF_VARtop10p), scale. = TRUE); ASAS0400bxF_PCAdf <- data.frame(PC1 = ASAS0400bxF_PCA$x[, 1], PC2 = ASAS0400bxF_PCA$x[, 2])
ASAS0400mxF_VAR <- apply(ASAS0400mxF, 1, var); ASAS0400mxF_VARtop10p <- ASAS0400mxF[names(sort(ASAS0400mxF_VAR, decreasing = TRUE))[1:round(nrow(ASAS0400mxF)*0.1,0)], ]; ASAS0400mxF_PCA <- prcomp(t(ASAS0400mxF_VARtop10p), scale. = TRUE); ASAS0400mxF_PCAdf <- data.frame(PC1 = ASAS0400mxF_PCA$x[, 1], PC2 = ASAS0400mxF_PCA$x[, 2])
ASAS1528bxF_VAR <- apply(ASAS1528bxF, 1, var); ASAS1528bxF_VARtop10p <- ASAS1528bxF[names(sort(ASAS1528bxF_VAR, decreasing = TRUE))[1:round(nrow(ASAS1528bxF)*0.1,0)], ]; ASAS1528bxF_PCA <- prcomp(t(ASAS1528bxF_VARtop10p), scale. = TRUE); ASAS1528bxF_PCAdf <- data.frame(PC1 = ASAS1528bxF_PCA$x[, 1], PC2 = ASAS1528bxF_PCA$x[, 2])
ASAS1528mxF_VAR <- apply(ASAS1528mxF, 1, var); ASAS1528mxF_VARtop10p <- ASAS1528mxF[names(sort(ASAS1528mxF_VAR, decreasing = TRUE))[1:round(nrow(ASAS1528mxF)*0.1,0)], ]; ASAS1528mxF_PCA <- prcomp(t(ASAS1528mxF_VARtop10p), scale. = TRUE); ASAS1528mxF_PCAdf <- data.frame(PC1 = ASAS1528mxF_PCA$x[, 1], PC2 = ASAS1528mxF_PCA$x[, 2])
CITY0822bxF_VAR <- apply(CITY0822bxF, 1, var); CITY0822bxF_VARtop10p <- CITY0822bxF[names(sort(CITY0822bxF_VAR, decreasing = TRUE))[1:round(nrow(CITY0822bxF)*0.1,0)], ]; CITY0822bxF_PCA <- prcomp(t(CITY0822bxF_VARtop10p), scale. = TRUE); CITY0822bxF_PCAdf <- data.frame(PC1 = CITY0822bxF_PCA$x[, 1], PC2 = CITY0822bxF_PCA$x[, 2])
CITY0822mxF_VAR <- apply(CITY0822mxF, 1, var); CITY0822mxF_VARtop10p <- CITY0822mxF[names(sort(CITY0822mxF_VAR, decreasing = TRUE))[1:round(nrow(CITY0822mxF)*0.1,0)], ]; CITY0822mxF_PCA <- prcomp(t(CITY0822mxF_VARtop10p), scale. = TRUE); CITY0822mxF_PCAdf <- data.frame(PC1 = CITY0822mxF_PCA$x[, 1], PC2 = CITY0822mxF_PCA$x[, 2])

ASAS0400bxM_PCAdfp <- ggplot(ASAS0400bxM_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \nbeta value,   male (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(ASAS0400bxM_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(ASAS0400bxM_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
ASAS0400mxM_PCAdfp <- ggplot(ASAS0400mxM_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \n   M value,   male (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(ASAS0400mxM_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(ASAS0400mxM_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
ASAS0400bxF_PCAdfp <- ggplot(ASAS0400bxF_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \nbeta value, female (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(ASAS0400bxF_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(ASAS0400bxF_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
ASAS0400mxF_PCAdfp <- ggplot(ASAS0400mxF_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \n   M value, female (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(ASAS0400mxF_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(ASAS0400mxF_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
ASAS1528bxM_PCAdfp <- ggplot(ASAS1528bxM_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \nbeta value,   male (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(ASAS1528bxM_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(ASAS1528bxM_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
ASAS1528mxM_PCAdfp <- ggplot(ASAS1528mxM_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \n   M value,   male (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(ASAS1528mxM_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(ASAS1528mxM_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
ASAS1528bxF_PCAdfp <- ggplot(ASAS1528bxF_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \nbeta value, female (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(ASAS1528bxF_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(ASAS1528bxF_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
ASAS1528mxF_PCAdfp <- ggplot(ASAS1528mxF_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \n   M value, female (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(ASAS1528mxF_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(ASAS1528mxF_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
CITY0822bxM_PCAdfp <- ggplot(CITY0822bxM_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \nbeta value,   male (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(CITY0822bxM_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(CITY0822bxM_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
CITY0822mxM_PCAdfp <- ggplot(CITY0822mxM_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \n   M value,   male (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(CITY0822mxM_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(CITY0822mxM_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
CITY0822bxF_PCAdfp <- ggplot(CITY0822bxF_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \nbeta value, female (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(CITY0822bxF_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(CITY0822bxF_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
CITY0822mxF_PCAdfp <- ggplot(CITY0822mxF_PCAdf, aes(x = PC1, y = PC2)) + geom_point(alpha = 0.7, size = 2) + theme_minimal() + labs(title = "PCA of DNA Methylation, \n   M value, female (Top 10% Variable Probes)",x = paste0("PC1 (", round(summary(CITY0822mxF_PCA)$importance[2,1]*100, 1), "%)"),y = paste0("PC2 (", round(summary(CITY0822mxF_PCA)$importance[2,2]*100, 1), "%)")) + scale_color_brewer(palette = "Set1")
ggsave("01_PCA.pdf", plot = grid.arrange(ASAS0400bxM_PCAdfp,ASAS0400mxM_PCAdfp,ASAS0400bxF_PCAdfp,ASAS0400mxF_PCAdfp,ASAS1528bxM_PCAdfp,ASAS1528mxM_PCAdfp,ASAS1528bxF_PCAdfp,ASAS1528mxF_PCAdfp,CITY0822bxM_PCAdfp,CITY0822mxM_PCAdfp,CITY0822bxF_PCAdfp,CITY0822mxF_PCAdfp, ncol=4), width=20, height=15)

pdf("02_Boxplot.pdf", width = 20, height = 36)
par(mfrow = c(12, 1))
boxplot(ASAS0400bxM_VARtop10p,main = "Methylation Distribution,   male",xaxt = "n",xlab = "Samples",ylab = "beta value",outline = FALSE)
boxplot(ASAS0400mxM_VARtop10p,main = "Methylation Distribution,   male",xaxt = "n",xlab = "Samples",ylab = "   M value",outline = FALSE)
boxplot(ASAS0400bxF_VARtop10p,main = "Methylation Distribution, female",xaxt = "n",xlab = "Samples",ylab = "beta value",outline = FALSE)
boxplot(ASAS0400mxF_VARtop10p,main = "Methylation Distribution, female",xaxt = "n",xlab = "Samples",ylab = "   M value",outline = FALSE)
boxplot(ASAS1528bxM_VARtop10p,main = "Methylation Distribution,   male",xaxt = "n",xlab = "Samples",ylab = "beta value",outline = FALSE)
boxplot(ASAS1528mxM_VARtop10p,main = "Methylation Distribution,   male",xaxt = "n",xlab = "Samples",ylab = "   M value",outline = FALSE)
boxplot(ASAS1528bxF_VARtop10p,main = "Methylation Distribution, female",xaxt = "n",xlab = "Samples",ylab = "beta value",outline = FALSE)
boxplot(ASAS1528mxF_VARtop10p,main = "Methylation Distribution, female",xaxt = "n",xlab = "Samples",ylab = "   M value",outline = FALSE)
boxplot(CITY0822bxM_VARtop10p,main = "Methylation Distribution,   male",xaxt = "n",xlab = "Samples",ylab = "beta value",outline = FALSE)
boxplot(CITY0822mxM_VARtop10p,main = "Methylation Distribution,   male",xaxt = "n",xlab = "Samples",ylab = "   M value",outline = FALSE)
boxplot(CITY0822bxF_VARtop10p,main = "Methylation Distribution, female",xaxt = "n",xlab = "Samples",ylab = "beta value",outline = FALSE)
boxplot(CITY0822mxF_VARtop10p,main = "Methylation Distribution, female",xaxt = "n",xlab = "Samples",ylab = "   M value",outline = FALSE)
dev.off()

ASAS0400bxM_VARtop1p <- ASAS0400bxM[names(sort(ASAS0400bxM_VAR, decreasing = TRUE))[1:round(nrow(ASAS0400bxM)*0.01,0)], ]
ASAS0400mxM_VARtop1p <- ASAS0400mxM[names(sort(ASAS0400mxM_VAR, decreasing = TRUE))[1:round(nrow(ASAS0400mxM)*0.01,0)], ]
ASAS1528bxM_VARtop1p <- ASAS1528bxM[names(sort(ASAS1528bxM_VAR, decreasing = TRUE))[1:round(nrow(ASAS1528bxM)*0.01,0)], ]
ASAS1528mxM_VARtop1p <- ASAS1528mxM[names(sort(ASAS1528mxM_VAR, decreasing = TRUE))[1:round(nrow(ASAS1528mxM)*0.01,0)], ]
CITY0822bxM_VARtop1p <- CITY0822bxM[names(sort(CITY0822bxM_VAR, decreasing = TRUE))[1:round(nrow(CITY0822bxM)*0.01,0)], ]
CITY0822mxM_VARtop1p <- CITY0822mxM[names(sort(CITY0822mxM_VAR, decreasing = TRUE))[1:round(nrow(CITY0822mxM)*0.01,0)], ]
ASAS0400bxF_VARtop1p <- ASAS0400bxF[names(sort(ASAS0400bxF_VAR, decreasing = TRUE))[1:round(nrow(ASAS0400bxF)*0.01,0)], ]
ASAS0400mxF_VARtop1p <- ASAS0400mxF[names(sort(ASAS0400mxF_VAR, decreasing = TRUE))[1:round(nrow(ASAS0400mxF)*0.01,0)], ]
ASAS1528bxF_VARtop1p <- ASAS1528bxF[names(sort(ASAS1528bxF_VAR, decreasing = TRUE))[1:round(nrow(ASAS1528bxF)*0.01,0)], ]
ASAS1528mxF_VARtop1p <- ASAS1528mxF[names(sort(ASAS1528mxF_VAR, decreasing = TRUE))[1:round(nrow(ASAS1528mxF)*0.01,0)], ]
CITY0822bxF_VARtop1p <- CITY0822bxF[names(sort(CITY0822bxF_VAR, decreasing = TRUE))[1:round(nrow(CITY0822bxF)*0.01,0)], ]
CITY0822mxF_VARtop1p <- CITY0822mxF[names(sort(CITY0822mxF_VAR, decreasing = TRUE))[1:round(nrow(CITY0822mxF)*0.01,0)], ]

pdf("03_HM_ALL12.pdf", width = 8, height = 8)
pheatmap(ASAS0400bxM_VARtop1p, main="Landscape,   male, beta value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(ASAS0400mxM_VARtop1p, main="Landscape,   male,    M value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(ASAS0400bxF_VARtop1p, main="Landscape, female, beta value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(ASAS0400mxF_VARtop1p, main="Landscape, female,    M value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(ASAS1528bxM_VARtop1p, main="Landscape,   male, beta value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(ASAS1528mxM_VARtop1p, main="Landscape,   male,    M value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(ASAS1528bxF_VARtop1p, main="Landscape, female, beta value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(ASAS1528mxF_VARtop1p, main="Landscape, female,    M value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(CITY0822bxM_VARtop1p, main="Landscape,   male, beta value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(CITY0822mxM_VARtop1p, main="Landscape,   male,    M value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(CITY0822bxF_VARtop1p, main="Landscape, female, beta value", show_colnames=FALSE, show_rownames=FALSE)
pheatmap(CITY0822mxF_VARtop1p, main="Landscape, female,    M value", show_colnames=FALSE, show_rownames=FALSE)
dev.off()

fIQR <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE); Q3 <- quantile(x, 0.75, na.rm = TRUE); IQR <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR; upper_bound <- Q3 + 1.5 * IQR; return(x < lower_bound | x > upper_bound)}
fMAH <- function(data, p_cutoff = 0.001) {
  center <- colMeans(data); cov_mat <- cov(data); dist <- mahalanobis(data, center, cov_mat)
  p_values <- pchisq(dist, df = ncol(data), lower.tail = FALSE); return(p_values < p_cutoff)}

ASAS0400ci <- apply(ASAS0400c, 2, fIQR); ASAS0400cix <- ASAS0400c[rowSums(ASAS0400ci) > 0, ]; dim(ASAS0400cix)
ASAS1528ci <- apply(ASAS1528c, 2, fIQR); ASAS1528cix <- ASAS1528c[rowSums(ASAS1528ci) > 0, ]; dim(ASAS1528cix)
CITY0822ci <- apply(CITY0822c, 2, fIQR); CITY0822cix <- CITY0822c[rowSums(CITY0822ci) > 0, ]; dim(CITY0822cix)
ASAS0400cm <- fMAH(ASAS0400c, p_cutoff = 0.001); ASAS0400cmx <- ASAS0400c[ASAS0400cm, ]; dim(ASAS0400cmx)
ASAS1528cm <- fMAH(ASAS1528c, p_cutoff = 0.001); ASAS1528cmx <- ASAS1528c[ASAS1528cm, ]; dim(ASAS1528cmx)
CITY0822cm <- fMAH(CITY0822c, p_cutoff = 0.001); CITY0822cmx <- CITY0822c[CITY0822cm, ]; dim(CITY0822cmx)


