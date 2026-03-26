################################################################################################
#Module 00. Clear memory and Load library
################################################################################################
sTIME <- Sys.time()
rm(list=ls())
library(dplyr); library(tidyr); library(data.table); library(pheatmap); library(readxl)
library(ggvenn); library(gridExtra); library(viridis); library(ggpubr); library(grid)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
Rep_NA <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE); Q3 <- quantile(x, 0.75, na.rm = TRUE); IQR_val <- Q3 - Q1;
  lower <- Q1 - 1.5 * IQR_val;   upper <- Q3 + 1.5 * IQR_val; ifelse(x < lower | x > upper, NA, x)}
custom_col1 <- colorRampPalette(c("#60C5F1", "#FFFCCC", "#F47378"))(500)
custom_col2 <- colorRampPalette(c("#FFFCCC", "#D01C8B"))(1000)
custom_col3 <- colorRampPalette(c("#4DAC26", "#FFFCCC", "#CA0020"))(500)
custom_col4 <- colorRampPalette(c("#0571B0", "#FFFCCC", "#D01C8B"))(500)
custom_col5 <- colorRampPalette(c("#018571", "#FFFCCC", "#E66101"))(500)

annotation_colors = list(
  TYPE  = c("A"="#0D3692", "C"="#33A23D", "D"="#FE5B10", "R"="#32A1C8"),
  SEX  = c("M"="#0D3692", "F"="#CB333B"),
  AGE    = c("28"="#FFF3C7", "22"="#138535"),
  MSC    = c("28"="#FFF3C7", "22"="#B21016"),
  EXF    = c("28"="#FFF3C7", "22"="#232B99"),
  PV = c("28"="#FFF3C7", "22"="#F7418F"))

setwd("G:/내 드라이브/Rrunning/HW014/")
################################################################################################
#Module 01. Data loading
################################################################################################
#Loading KoGES methylation: Takes 15 mins
Sys.time()
source("Loading_KoGES_Met_ASAS0400_ASAS1528_CITY0822_Epid_QC.R")
#Loading TCGA methylation: Takes 15 mins
Sys.time()
source("Loading_TCGA_Met_All.R")
#Loading GSE223748 methylation: Takes 2 mins
Sys.time()
source("Loading_GSE223748_Module01.R")
Sys.time()

#Dataframe structure, check row and column number
Met_ASAS0400b[1:5,1:5]; dim(Met_ASAS0400b) #451245, 400
Met_ASAS1528b[1:5,1:5]; dim(Met_ASAS1528b) #805110,1528
Met_CITY0822b[1:5,1:5]; dim(Met_CITY0822b) #803393, 822
Met_ASAS0400m[1:5,1:5]; dim(Met_ASAS0400m) #451245, 400
Met_ASAS1528m[1:5,1:5]; dim(Met_ASAS1528m) #805110,1528
Met_CITY0822m[1:5,1:5]; dim(Met_CITY0822m) #803393, 822
GSE223748m[1:5, 1:5]; dim(GSE223748m) #37554,15043
GSE223748i[1:5, 1:5]; dim(GSE223748i) #15043,12

#INT_kte_CpG <- Reduce(intersect, list(rownames(Met_ASAS0400),rownames(Met_ASAS1528),rownames(Met_CITY0822),rownames(ACC__450m),rownames(BLCA_450m),rownames(BRCA_450m),rownames(CESC_450m),rownames(CHOL_450m),rownames(COAD_450m),rownames(DLBC_450m),rownames(ESCA_450m),rownames(GBM__450m),rownames(HNSC_450m),rownames(KICH_450m),rownames(KIRC_450m),rownames(KIRP_450m),rownames(LAML_450m),rownames(LGG__450m),rownames(LIHC_450m),rownames(LUAD_450m),rownames(LUSC_450m),rownames(MESO_450m),rownames(OV___450m),rownames(PAAD_450m),rownames(PCPG_450m),rownames(PRAD_450m),rownames(READ_450m),rownames(SKCM_450m),rownames(STAD_450m),rownames(TGCT_450m),rownames(THCA_450m),rownames(THYM_450m),rownames(UCEC_450m),rownames(UCS__450m),rownames(UVM__450m),rownames(GSE223748m)))
#INT_ktx_CpG <- Reduce(intersect, list(rownames(Met_ASAS0400),rownames(Met_ASAS1528),rownames(Met_CITY0822),rownames(ACC__450m),rownames(BLCA_450m),rownames(BRCA_450m),rownames(CESC_450m),rownames(CHOL_450m),rownames(COAD_450m),rownames(DLBC_450m),rownames(ESCA_450m),rownames(GBM__450m),rownames(HNSC_450m),rownames(KICH_450m),rownames(KIRC_450m),rownames(KIRP_450m),rownames(LAML_450m),rownames(LGG__450m),rownames(LIHC_450m),rownames(LUAD_450m),rownames(LUSC_450m),rownames(MESO_450m),rownames(OV___450m),rownames(PAAD_450m),rownames(PCPG_450m),rownames(PRAD_450m),rownames(READ_450m),rownames(SKCM_450m),rownames(STAD_450m),rownames(TGCT_450m),rownames(THCA_450m),rownames(THYM_450m),rownames(UCEC_450m),rownames(UCS__450m),rownames(UVM__450m)))
#INT_kex_CpG <- Reduce(intersect, list(rownames(Met_ASAS0400),rownames(Met_ASAS1528),rownames(Met_CITY0822),rownames(GSE223748m)))
#INT_etx_CpG <- Reduce(intersect, list(rownames(ACC__450m),rownames(BLCA_450m),rownames(BRCA_450m),rownames(CESC_450m),rownames(CHOL_450m),rownames(COAD_450m),rownames(DLBC_450m),rownames(ESCA_450m),rownames(GBM__450m),rownames(HNSC_450m),rownames(KICH_450m),rownames(KIRC_450m),rownames(KIRP_450m),rownames(LAML_450m),rownames(LGG__450m),rownames(LIHC_450m),rownames(LUAD_450m),rownames(LUSC_450m),rownames(MESO_450m),rownames(OV___450m),rownames(PAAD_450m),rownames(PCPG_450m),rownames(PRAD_450m),rownames(READ_450m),rownames(SKCM_450m),rownames(STAD_450m),rownames(TGCT_450m),rownames(THCA_450m),rownames(THYM_450m),rownames(UCEC_450m),rownames(UCS__450m),rownames(UVM__450m),rownames(GSE223748m)))
#write.table(INT_kte_CpG,file="INT_kte_CpG_QC.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
#write.table(INT_ktx_CpG,file="INT_ktx_CpG_QC.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
#write.table(INT_kex_CpG,file="INT_kex_CpG_QC.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
#write.table(INT_etx_CpG,file="INT_etx_CpG_QC.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
INT_kte_CpG <- as.vector(unlist(read.table("INT_kte_CpG_QC.txt")))
INT_ktx_CpG <- as.vector(unlist(read.table("INT_ktx_CpG_QC.txt")))
INT_kex_CpG <- as.vector(unlist(read.table("INT_kex_CpG_QC.txt")))
INT_etx_CpG <- as.vector(unlist(read.table("INT_etx_CpG_QC.txt")))
length(INT_kte_CpG) #5306 > 5213
length(INT_ktx_CpG) #452453 > 413411
length(INT_kex_CpG) #5306 > 5213
length(INT_etx_CpG) #5498 > 5498


ASAS0400c <- as.data.frame(fread("G:/내 드라이브/Rrunning/KoGES/Geno/Met/idat_to_tsv/ASAS0400c.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS0400c) <- ASAS0400c$V1; ASAS0400c$V1 <- NULL
ASAS1528c <- as.data.frame(fread("G:/내 드라이브/Rrunning/KoGES/Geno/Met/idat_to_tsv/ASAS1528c.txt.gz",  header = TRUE, sep = "\t")); rownames(ASAS1528c) <- ASAS1528c$V1; ASAS1528c$V1 <- NULL
CITY0822c <- as.data.frame(fread("G:/내 드라이브/Rrunning/KoGES/Geno/Met/idat_to_tsv/CITY0822c.txt.gz",  header = TRUE, sep = "\t")); rownames(CITY0822c) <- CITY0822c$V1; CITY0822c$V1 <- NULL

head(ASAS0400c); dim(ASAS0400c)
head(ASAS1528c); dim(ASAS1528c)
head(CITY0822c); dim(CITY0822c)

fIQR <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE); Q3 <- quantile(x, 0.75, na.rm = TRUE); IQR <- Q3 - Q1
  lower_bound <- Q1 - 2 * IQR; upper_bound <- Q3 + 2 * IQR; return(x < lower_bound | x > upper_bound)}
fMAH <- function(data, p_cutoff = 0.05) {
  center <- colMeans(data); cov_mat <- cov(data); dist <- mahalanobis(data, center, cov_mat)
  p_values <- pchisq(dist, df = ncol(data), lower.tail = FALSE); return(p_values < p_cutoff)}

ASAS0400ci <- apply(ASAS0400c, 2, fIQR); ASAS0400cix <- ASAS0400c[rowSums(ASAS0400ci) > 0, ]; dim(ASAS0400cix) #13
ASAS1528ci <- apply(ASAS1528c, 2, fIQR); ASAS1528cix <- ASAS1528c[rowSums(ASAS1528ci) > 0, ]; dim(ASAS1528cix) #34
CITY0822ci <- apply(CITY0822c, 2, fIQR); CITY0822cix <- CITY0822c[rowSums(CITY0822ci) > 0, ]; dim(CITY0822cix) #14
ASAS0400cm <- fMAH(ASAS0400c, p_cutoff = 0.01); ASAS0400cmx <- ASAS0400c[ASAS0400cm, ]; dim(ASAS0400cmx) #15
ASAS1528cm <- fMAH(ASAS1528c, p_cutoff = 0.01); ASAS1528cmx <- ASAS1528c[ASAS1528cm, ]; dim(ASAS1528cmx) #47
CITY0822cm <- fMAH(CITY0822c, p_cutoff = 0.01); CITY0822cmx <- CITY0822c[CITY0822cm, ]; dim(CITY0822cmx) #23

DROP_ASAS0400 <- intersect(rownames(ASAS0400cix),rownames(ASAS0400cmx)); length(DROP_ASAS0400) # 9
DROP_ASAS1528 <- intersect(rownames(ASAS1528cix),rownames(ASAS1528cmx)); length(DROP_ASAS1528) #20
DROP_CITY0822 <- intersect(rownames(CITY0822cix),rownames(CITY0822cmx)); length(DROP_CITY0822) #10

Met_ASAS0400b <- dplyr::select(Met_ASAS0400b, -DROP_ASAS0400)
Met_ASAS1528b <- dplyr::select(Met_ASAS1528b, -DROP_ASAS1528)
Met_CITY0822b <- dplyr::select(Met_CITY0822b, -DROP_CITY0822)
Met_ASAS0400m <- dplyr::select(Met_ASAS0400m, -DROP_ASAS0400)
Met_ASAS1528m <- dplyr::select(Met_ASAS1528m, -DROP_ASAS1528)
Met_CITY0822m <- dplyr::select(Met_CITY0822m, -DROP_CITY0822)

Met_ASAS0400b[1:5,1:5]; dim(Met_ASAS0400b) #451245, 391
Met_ASAS1528b[1:5,1:5]; dim(Met_ASAS1528b) #805110,1508
Met_CITY0822b[1:5,1:5]; dim(Met_CITY0822b) #803393, 812
Met_ASAS0400m[1:5,1:5]; dim(Met_ASAS0400m) #451245, 391
Met_ASAS1528m[1:5,1:5]; dim(Met_ASAS1528m) #805110,1508
Met_CITY0822m[1:5,1:5]; dim(Met_CITY0822m) #803393, 812

COMMON_A0105 <- intersect(colnames(Met_ASAS0400b),colnames(Met_ASAS1528b))

#hg19
GPL13534 <- dplyr::select(as.data.frame(fread("G:/내 드라이브/Rrunning/Annotations/21_Methyl/GPL13534-11288.txt.gz",  header = TRUE, sep = "\t")), ID, CHR,MAPINFO,UCSC_RefGene_Name,UCSC_RefGene_Group,Relation_to_UCSC_CpG_Island,Probe_SNPs)
GPL28271 <- dplyr::select(as.data.frame(fread("G:/내 드라이브/Rrunning/Annotations/21_Methyl/GPL28271-57075.txt.gz",  header = TRUE, sep = "\t")), ID, Human.CHR,Human.Hg19_CGstart,Human.Hg19_SYMBOL)
head(GPL13534)
head(GPL28271)

################################################################################################
#Module Ox11_OXI
library(ggplot2); library(reshape2); library(pheatmap); library(gridExtra); library(ape); library(dendextend); library(ggvenn)

#OBS DRSM
ASAS01to10$NA01_DR_DRK <- ASAS01to10$AS01_DRINK ; ASAS01to10$NA05_DR_DRK <- ASAS01to10$AS05_DRINK 
ASAS01to10$NA01_DR_DUR <- ASAS01to10$AS01_DRDUA ; ASAS01to10$NA05_DR_DUR <- ASAS01to10$AS05_DRDU  
ASAS01to10$NA01_DR_DUR <- ifelse(ASAS01to10$NA01_DR_DUR==4,3,ASAS01to10$NA01_DR_DUR)
ASAS01to10$NA01_DR_DUR <- ifelse(ASAS01to10$NA01_DR_DUR==5,4,ASAS01to10$NA01_DR_DUR)
ASAS01to10$NA01_DR_TOT <- ASAS01to10$AS01_TOTALC; ASAS01to10$NA05_DR_TOT <- ASAS01to10$AS05_TOTALC
ASAS01to10$NA01_SM_SMK <- ASAS01to10$AS01_SMOKEA; ASAS01to10$NA05_SM_SMK <- ASAS01to10$AS05_SMOKE 
ASAS01to10$NA01_SM_PYR <- ASAS01to10$AS01_PACKYR; ASAS01to10$NA05_SM_PYR <- ASAS01to10$AS05_PACKYR
ASAS01to10$NA01_SM_SMK <- ASAS01to10$NA01_SM_SMK + 1; ASAS01to10$NA01_SM_SMK <- ifelse(ASAS01to10$NA01_SM_SMK == 4, 3, ASAS01to10$NA01_SM_SMK)

CITY01to02$NC01_DR_DRK <- CITY01to02$CT01_DRINK
CITY01to02$NC01_DR_DUR <- CITY01to02$CT01_DRDU
CITY01to02$NC01_SM_SMK <- CITY01to02$CT01_SMOKE

ASAS01to10$NA01_DR_DUR <- ifelse(ASAS01to10$NA01_DR_DRK==1, 0, ASAS01to10$NA01_DR_DUR)
head(dplyr::select(dplyr::filter(ASAS01to10, NA01_DR_DRK==1), NA01_DR_DRK,NA01_DR_DUR))
ASAS01to10$NA05_DR_DUR <- ifelse(ASAS01to10$NA05_DR_DRK==1, 0, ASAS01to10$NA05_DR_DUR)
head(dplyr::select(dplyr::filter(ASAS01to10, NA05_DR_DRK==1), NA05_DR_DRK,NA05_DR_DUR))
CITY01to02$NC01_DR_DUR <- ifelse(CITY01to02$NC01_DR_DRK==1, 0, CITY01to02$NC01_DR_DUR)
head(dplyr::select(dplyr::filter(CITY01to02, NC01_DR_DRK==1), NC01_DR_DRK,NC01_DR_DUR))

ASAS01to10$NA01_DR_TOT <- ifelse(ASAS01to10$NA01_DR_DRK==1 | ASAS01to10$NA01_DR_DRK==2, 0, ASAS01to10$NA01_DR_TOT)
head(dplyr::select(dplyr::filter(ASAS01to10, NA01_DR_DRK==1), NA01_DR_DRK,NA01_DR_TOT))
head(dplyr::select(dplyr::filter(ASAS01to10, NA01_DR_DRK==2), NA01_DR_DRK,NA01_DR_TOT))
head(dplyr::select(dplyr::filter(ASAS01to10, NA01_DR_DRK==3), NA01_DR_DRK,NA01_DR_TOT))

ASAS01to10$NA01_SM_PYR <- ifelse(ASAS01to10$NA01_SM_SMK==1, 0, ASAS01to10$NA01_SM_PYR)
head(dplyr::select(dplyr::filter(ASAS01to10, NA01_SM_SMK==1), NA01_SM_SMK,NA01_SM_PYR))
head(dplyr::select(dplyr::filter(ASAS01to10, NA01_SM_SMK==2), NA01_SM_SMK,NA01_SM_PYR))
head(dplyr::select(dplyr::filter(ASAS01to10, NA01_SM_SMK==3), NA01_SM_SMK,NA01_SM_PYR))

ASAS01to10$NA05_DR_TOT <- ifelse(ASAS01to10$NA05_DR_DRK==1 | ASAS01to10$NA05_DR_DRK==2, 0, ASAS01to10$NA05_DR_TOT)
head(dplyr::select(dplyr::filter(ASAS01to10, NA05_DR_DRK==1), NA05_DR_DRK,NA05_DR_TOT))
head(dplyr::select(dplyr::filter(ASAS01to10, NA05_DR_DRK==2), NA05_DR_DRK,NA05_DR_TOT))
head(dplyr::select(dplyr::filter(ASAS01to10, NA05_DR_DRK==3), NA05_DR_DRK,NA05_DR_TOT))

ASAS01to10$NA05_SM_PYR <- ifelse(ASAS01to10$NA05_SM_SMK==1, 0, ASAS01to10$NA05_SM_PYR)
head(dplyr::select(dplyr::filter(ASAS01to10, NA05_SM_SMK==1), NA05_SM_SMK,NA05_SM_PYR))
head(dplyr::select(dplyr::filter(ASAS01to10, NA05_SM_SMK==2), NA05_SM_SMK,NA05_SM_PYR))
head(dplyr::select(dplyr::filter(ASAS01to10, NA05_SM_SMK==3), NA05_SM_SMK,NA05_SM_PYR))

table(ASAS01to10$NA01_DR_DRK); table(ASAS01to10$NA05_DR_DRK); table(CITY01to02$NC01_DR_DRK)
table(is.na(ASAS01to10$NA01_DR_DRK)); table(is.na(ASAS01to10$NA05_DR_DRK)); table(is.na(CITY01to02$NC01_DR_DRK))

table(ASAS01to10$NA01_DR_DUR); table(ASAS01to10$NA05_DR_DUR); table(CITY01to02$NC01_DR_DUR)
table(is.na(ASAS01to10$NA01_DR_DUR)); table(is.na(ASAS01to10$NA05_DR_DUR)); table(is.na(CITY01to02$NC01_DR_DUR))

table(ASAS01to10$NA01_SM_SMK); table(ASAS01to10$NA05_SM_SMK); table(CITY01to02$NC01_SM_SMK)
table(is.na(ASAS01to10$NA01_SM_SMK)); table(is.na(ASAS01to10$NA05_SM_SMK)); table(is.na(CITY01to02$NC01_SM_SMK))

summary(dplyr::select(ASAS01to10, NA01_DR_TOT,NA05_DR_TOT))
summary(dplyr::select(ASAS01to10, NA01_SM_PYR,NA05_SM_PYR))

#OBS Exercise
ASAS01to10$NA01_PHY_G1 <- ASAS01to10$AS01_PHYSTB
ASAS01to10$NA01_PHY_G2 <- ASAS01to10$AS01_PHYSIT
ASAS01to10$NA01_PHY_G3 <- ASAS01to10$AS01_PHYACTL
ASAS01to10$NA01_PHY_G4 <- ASAS01to10$AS01_PHYACTM
ASAS01to10$NA01_PHY_G5 <- ASAS01to10$AS01_PHYACTH
ASAS01to10$NC01_PHY_MT <- ASAS01to10$NA01_PHY_G1*1.0 + ASAS01to10$NA01_PHY_G2*2.5 + ASAS01to10$NA01_PHY_G3*4.0 + ASAS01to10$NA01_PHY_G4*8.0 + ASAS01to10$NA01_PHY_G5*12.0

CITY01to02$NC01_PHY_E1 <- CITY01to02$CT01_EXER
CITY01to02$NC01_PHY_E2 <- CITY01to02$CT01_EXERFQ
CITY01to02$NC01_PHY_E3 <- CITY01to02$CT01_EXERDU
CITY01to02$NC01_PHY_E4 <- CITY01to02$NC01_PHY_E2 * CITY01to02$NC01_PHY_E3

#Var selection
Ox01_OXI_ASAS01_LIST <- c("DIST_ID","AS01_SEX", "AS01_B21", "AS01_B19", "AS01_B12", "AS01_B13", "AS01_B16", "AS01_B17", "AS01_B14", "AS01_B23", "AS01_B05", "AS01_B15", "AS01_B03", "AS01_B07","AS01_BMI","NA01_DR_TOT","NA01_SM_SMK","NC01_PHY_MT")
Ox01_OXI_CITY01_LISx <- c("DIST_ID","CT01_SEX","CT01_SS21","CT01_SS19","CT01_SS12","CT01_SS13","CT01_SS16","CT01_SS17","CT01_SS14","CT01_SS23","CT01_SS05","CT01_SS15","CT01_SS03","CT01_SS07","CT01_BMI","NC01_DR_DRK")

A01_m04_OXI_M <- dplyr::filter(dplyr::select(ASAS01to10, DIST_ID, Ox01_OXI_ASAS01_LIST), DIST_ID %in% colnames(Met_ASAS0400b) & AS01_SEX==1); A01_m04_OXI_M <- A01_m04_OXI_M[complete.cases(A01_m04_OXI_M), ]; rownames(A01_m04_OXI_M) <- A01_m04_OXI_M$DIST_ID; A01_m04_OXI_M <- A01_m04_OXI_M[,-1:-2,drop=FALSE]
A01_m04_OXI_F <- dplyr::filter(dplyr::select(ASAS01to10, DIST_ID, Ox01_OXI_ASAS01_LIST), DIST_ID %in% colnames(Met_ASAS0400b) & AS01_SEX==2); A01_m04_OXI_F <- A01_m04_OXI_F[complete.cases(A01_m04_OXI_F), ]; rownames(A01_m04_OXI_F) <- A01_m04_OXI_F$DIST_ID; A01_m04_OXI_F <- A01_m04_OXI_F[,-1:-2,drop=FALSE]
C01_m08_OXI_M <- dplyr::filter(dplyr::select(CITY01to02, DIST_ID, Ox01_OXI_CITY01_LISx), DIST_ID %in% colnames(Met_CITY0822b) & CT01_SEX==1); C01_m08_OXI_M <- C01_m08_OXI_M[complete.cases(C01_m08_OXI_M), ]; rownames(C01_m08_OXI_M) <- C01_m08_OXI_M$DIST_ID; C01_m08_OXI_M <- C01_m08_OXI_M[,-1:-2,drop=FALSE]
C01_m08_OXI_F <- dplyr::filter(dplyr::select(CITY01to02, DIST_ID, Ox01_OXI_CITY01_LISx), DIST_ID %in% colnames(Met_CITY0822b) & CT01_SEX==2); C01_m08_OXI_F <- C01_m08_OXI_F[complete.cases(C01_m08_OXI_F), ]; rownames(C01_m08_OXI_F) <- C01_m08_OXI_F$DIST_ID; C01_m08_OXI_F <- C01_m08_OXI_F[,-1:-2,drop=FALSE]
head(A01_m04_OXI_M); dim(A01_m04_OXI_M); table(is.na(A01_m04_OXI_M)); sapply(A01_m04_OXI_M, function(x) length(unique(x))) #169,16
head(A01_m04_OXI_F); dim(A01_m04_OXI_F); table(is.na(A01_m04_OXI_F)); sapply(A01_m04_OXI_F, function(x) length(unique(x))) #173,16
head(C01_m08_OXI_M); dim(C01_m08_OXI_M); table(is.na(C01_m08_OXI_M)); sapply(C01_m08_OXI_M, function(x) length(unique(x))) #615,14
head(C01_m08_OXI_F); dim(C01_m08_OXI_F); table(is.na(C01_m08_OXI_F)); sapply(C01_m08_OXI_F, function(x) length(unique(x))) #194,14

m04bA01_OXI_M <- dplyr::select(Met_ASAS0400b[rownames(Met_ASAS0400b) %in% INT_kex_CpG, ], intersect(colnames(Met_ASAS0400b), rownames(A01_m04_OXI_M))); dim(m04bA01_OXI_M); table(colnames(m04bA01_OXI_M)==rownames(A01_m04_OXI_M)); table(is.na(m04bA01_OXI_M)); table(is.na(A01_m04_OXI_M))
m04bA01_OXI_F <- dplyr::select(Met_ASAS0400b[rownames(Met_ASAS0400b) %in% INT_kex_CpG, ], intersect(colnames(Met_ASAS0400b), rownames(A01_m04_OXI_F))); dim(m04bA01_OXI_F); table(colnames(m04bA01_OXI_F)==rownames(A01_m04_OXI_F)); table(is.na(m04bA01_OXI_F)); table(is.na(A01_m04_OXI_F))
m08bC01_OXI_M <- dplyr::select(Met_CITY0822b[rownames(Met_CITY0822b) %in% INT_kex_CpG, ], intersect(colnames(Met_CITY0822b), rownames(C01_m08_OXI_M))); dim(m08bC01_OXI_M); table(colnames(m08bC01_OXI_M)==rownames(C01_m08_OXI_M)); table(is.na(m08bC01_OXI_M)); table(is.na(C01_m08_OXI_M))
m08bC01_OXI_F <- dplyr::select(Met_CITY0822b[rownames(Met_CITY0822b) %in% INT_kex_CpG, ], intersect(colnames(Met_CITY0822b), rownames(C01_m08_OXI_F))); dim(m08bC01_OXI_F); table(colnames(m08bC01_OXI_F)==rownames(C01_m08_OXI_F)); table(is.na(m08bC01_OXI_F)); table(is.na(C01_m08_OXI_F))
m04mA01_OXI_M <- dplyr::select(Met_ASAS0400m[rownames(Met_ASAS0400m) %in% INT_kex_CpG, ], intersect(colnames(Met_ASAS0400m), rownames(A01_m04_OXI_M))); dim(m04mA01_OXI_M); table(colnames(m04mA01_OXI_M)==rownames(A01_m04_OXI_M)); table(is.na(m04mA01_OXI_M)); table(is.na(A01_m04_OXI_M))
m04mA01_OXI_F <- dplyr::select(Met_ASAS0400m[rownames(Met_ASAS0400m) %in% INT_kex_CpG, ], intersect(colnames(Met_ASAS0400m), rownames(A01_m04_OXI_F))); dim(m04mA01_OXI_F); table(colnames(m04mA01_OXI_F)==rownames(A01_m04_OXI_F)); table(is.na(m04mA01_OXI_F)); table(is.na(A01_m04_OXI_F))
m08mC01_OXI_M <- dplyr::select(Met_CITY0822m[rownames(Met_CITY0822m) %in% INT_kex_CpG, ], intersect(colnames(Met_CITY0822m), rownames(C01_m08_OXI_M))); dim(m08mC01_OXI_M); table(colnames(m08mC01_OXI_M)==rownames(C01_m08_OXI_M)); table(is.na(m08mC01_OXI_M)); table(is.na(C01_m08_OXI_M))
m08mC01_OXI_F <- dplyr::select(Met_CITY0822m[rownames(Met_CITY0822m) %in% INT_kex_CpG, ], intersect(colnames(Met_CITY0822m), rownames(C01_m08_OXI_F))); dim(m08mC01_OXI_F); table(colnames(m08mC01_OXI_F)==rownames(C01_m08_OXI_F)); table(is.na(m08mC01_OXI_F)); table(is.na(C01_m08_OXI_F))

#CRp <- cor(t(m04bA01_OXI_M), A01_m04_OXI_M, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m04bA01_OXI_M), A01_m04_OXI_M, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(A01_m04_OXI_M)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m04bA01_OXI_M_CRPV <- CRp; m04bA01_OXI_M_CRPV$CRs <- CRs$CRs; m04bA01_OXI_M_CRPV$PVp <- PVp$PVp; m04bA01_OXI_M_CRPV$PVs <- PVs$PVs; m04bA01_OXI_M_CRPV$CLASS <- "m04bA01_OXI_M"
#CRp <- cor(t(m04bA01_OXI_F), A01_m04_OXI_F, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m04bA01_OXI_F), A01_m04_OXI_F, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(A01_m04_OXI_F)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m04bA01_OXI_F_CRPV <- CRp; m04bA01_OXI_F_CRPV$CRs <- CRs$CRs; m04bA01_OXI_F_CRPV$PVp <- PVp$PVp; m04bA01_OXI_F_CRPV$PVs <- PVs$PVs; m04bA01_OXI_F_CRPV$CLASS <- "m04bA01_OXI_F"
#CRp <- cor(t(m08bC01_OXI_M), C01_m08_OXI_M, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m08bC01_OXI_M), C01_m08_OXI_M, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(C01_m08_OXI_M)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m08bC01_OXI_M_CRPV <- CRp; m08bC01_OXI_M_CRPV$CRs <- CRs$CRs; m08bC01_OXI_M_CRPV$PVp <- PVp$PVp; m08bC01_OXI_M_CRPV$PVs <- PVs$PVs; m08bC01_OXI_M_CRPV$CLASS <- "m08bC01_OXI_M"
#CRp <- cor(t(m08bC01_OXI_F), C01_m08_OXI_F, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m08bC01_OXI_F), C01_m08_OXI_F, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(C01_m08_OXI_F)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m08bC01_OXI_F_CRPV <- CRp; m08bC01_OXI_F_CRPV$CRs <- CRs$CRs; m08bC01_OXI_F_CRPV$PVp <- PVp$PVp; m08bC01_OXI_F_CRPV$PVs <- PVs$PVs; m08bC01_OXI_F_CRPV$CLASS <- "m08bC01_OXI_F"
#CRPVb_Ox01_OXI <- rbind(m04bA01_OXI_M_CRPV,m04bA01_OXI_F_CRPV,m08bC01_OXI_M_CRPV,m08bC01_OXI_F_CRPV); dim(CRPVb_Ox01_OXI); table(CRPVb_Ox01_OXI$CLASS)
#write.table(CRPVb_Ox01_OXI, file="Output_Ox11_OXI/CRPVb_Ox01_OXI.txt",sep="\t",quote=FALSE,row.names=FALSE)
#nrow(CRPVb_Ox01_OXI); CRPVb_Ox01_OXI <- dplyr::filter(CRPVb_Ox01_OXI, CRp * CRs > 0); nrow(CRPVb_Ox01_OXI); CRPVb_Ox01_OXI <- dplyr::filter(CRPVb_Ox01_OXI, PVp < 0.01 & PVs < 0.01); nrow(CRPVb_Ox01_OXI) #24804660 > 20824709 > 644500
#CRPVb_Ox01_OXI$HyPi <- (CRPVb_Ox01_OXI$CRp * -log10(CRPVb_Ox01_OXI$PVp)) + (CRPVb_Ox01_OXI$CRs * -log10(CRPVb_Ox01_OXI$PVs))
#write.table(CRPVb_Ox01_OXI, file="Output_Ox11_OXI/CRPVb_Ox01_OXI_HyPi.txt",sep="\t",quote=FALSE,row.names=FALSE)
#
#CRp <- cor(t(m04mA01_OXI_M), A01_m04_OXI_M, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m04mA01_OXI_M), A01_m04_OXI_M, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(A01_m04_OXI_M)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m04mA01_OXI_M_CRPV <- CRp; m04mA01_OXI_M_CRPV$CRs <- CRs$CRs; m04mA01_OXI_M_CRPV$PVp <- PVp$PVp; m04mA01_OXI_M_CRPV$PVs <- PVs$PVs; m04mA01_OXI_M_CRPV$CLASS <- "m04mA01_OXI_M"
#CRp <- cor(t(m04mA01_OXI_F), A01_m04_OXI_F, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m04mA01_OXI_F), A01_m04_OXI_F, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(A01_m04_OXI_F)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m04mA01_OXI_F_CRPV <- CRp; m04mA01_OXI_F_CRPV$CRs <- CRs$CRs; m04mA01_OXI_F_CRPV$PVp <- PVp$PVp; m04mA01_OXI_F_CRPV$PVs <- PVs$PVs; m04mA01_OXI_F_CRPV$CLASS <- "m04mA01_OXI_F"
#CRp <- cor(t(m08mC01_OXI_M), C01_m08_OXI_M, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m08mC01_OXI_M), C01_m08_OXI_M, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(C01_m08_OXI_M)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m08mC01_OXI_M_CRPV <- CRp; m08mC01_OXI_M_CRPV$CRs <- CRs$CRs; m08mC01_OXI_M_CRPV$PVp <- PVp$PVp; m08mC01_OXI_M_CRPV$PVs <- PVs$PVs; m08mC01_OXI_M_CRPV$CLASS <- "m08mC01_OXI_M"
#CRp <- cor(t(m08mC01_OXI_F), C01_m08_OXI_F, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m08mC01_OXI_F), C01_m08_OXI_F, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(C01_m08_OXI_F)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m08mC01_OXI_F_CRPV <- CRp; m08mC01_OXI_F_CRPV$CRs <- CRs$CRs; m08mC01_OXI_F_CRPV$PVp <- PVp$PVp; m08mC01_OXI_F_CRPV$PVs <- PVs$PVs; m08mC01_OXI_F_CRPV$CLASS <- "m08mC01_OXI_F"
#CRPVm_Ox01_OXI <- rbind(m04mA01_OXI_M_CRPV,m04mA01_OXI_F_CRPV,m08mC01_OXI_M_CRPV,m08mC01_OXI_F_CRPV); dim(CRPVm_Ox01_OXI); table(CRPVm_Ox01_OXI$CLASS)
#write.table(CRPVm_Ox01_OXI, file="Output_Ox11_OXI/CRPVm_Ox01_OXI.txt",sep="\t",quote=FALSE,row.names=FALSE)
#nrow(CRPVm_Ox01_OXI); CRPVm_Ox01_OXI <- dplyr::filter(CRPVm_Ox01_OXI, CRp * CRs > 0); nrow(CRPVm_Ox01_OXI); CRPVm_Ox01_OXI <- dplyr::filter(CRPVm_Ox01_OXI, PVp < 0.01 & PVs < 0.01); nrow(CRPVm_Ox01_OXI) #24804660 > 21026518 > 686098
#CRPVm_Ox01_OXI$HyPi <- (CRPVm_Ox01_OXI$CRp * -log10(CRPVm_Ox01_OXI$PVp)) + (CRPVm_Ox01_OXI$CRs * -log10(CRPVm_Ox01_OXI$PVs))
#write.table(CRPVm_Ox01_OXI, file="Output_Ox11_OXI/CRPVm_Ox01_OXI_HyPi.txt",sep="\t",quote=FALSE,row.names=FALSE)


head(A01_m04_OXI_M); dim(A01_m04_OXI_M)
head(A01_m04_OXI_F); dim(A01_m04_OXI_F)
head(C01_m08_OXI_M); dim(C01_m08_OXI_M)
head(C01_m08_OXI_F); dim(C01_m08_OXI_F)

SCR_U <- function(x) {cuts <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE); ifelse(x <= cuts[1], 0, ifelse(x <= cuts[2], 1, 2))}
SCR_D <- function(x) {cuts <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE); ifelse(x <= cuts[1], 2, ifelse(x <= cuts[2], 1, 0))}
SCR_0 <- function(x) {non_zero_values <- x[x != 0 & !is.na(x)]
if (length(non_zero_values) == 0) {return(ifelse(x == 0, 2, NA))}
cutoff <- median(non_zero_values)
ifelse(x == 0, 2, ifelse(x <= cutoff, 1, 0))}

A_OBS_LIST_U <- c("AS01_B21","AS01_B19","AS01_B12","AS01_B13","AS01_B16","AS01_B17","AS01_B14","AS01_B23","AS01_B05","AS01_B15","NC01_PHY_MT")
A_OBS_LIST_D <- c("AS01_B03","AS01_B07","AS01_BMI")
C_OBS_LIST_U <- c("CT01_SS21","CT01_SS19","CT01_SS12","CT01_SS13","CT01_SS16","CT01_SS17","CT01_SS14","CT01_SS23","CT01_SS05","CT01_SS15")
C_OBS_LIST_D <- c("CT01_SS03","CT01_SS07","CT01_BMI")

A01_m04_OXI_M$NA01_SM_SMK <- as.numeric(substr(ifelse(A01_m04_OXI_M$NA01_SM_SMK==1,"S2",ifelse(A01_m04_OXI_M$NA01_SM_SMK==2,"S1","S0")),2,2))
A01_m04_OXI_F$NA01_SM_SMK <- as.numeric(substr(ifelse(A01_m04_OXI_F$NA01_SM_SMK==1,"S2",ifelse(A01_m04_OXI_F$NA01_SM_SMK==2,"S1","S0")),2,2))
C01_m08_OXI_M$NC01_DR_DRK <- as.numeric(substr(ifelse(C01_m08_OXI_M$NC01_DR_DRK==3,"S0","S1"),2,2))
C01_m08_OXI_F$NC01_DR_DRK <- as.numeric(substr(ifelse(C01_m08_OXI_F$NC01_DR_DRK==3,"S0","S1"),2,2))
A01_m04_OXI_M[A_OBS_LIST_U] <- lapply(A01_m04_OXI_M[A_OBS_LIST_U], SCR_U)
A01_m04_OXI_F[A_OBS_LIST_U] <- lapply(A01_m04_OXI_F[A_OBS_LIST_U], SCR_U)
C01_m08_OXI_M[C_OBS_LIST_U] <- lapply(C01_m08_OXI_M[C_OBS_LIST_U], SCR_U)
C01_m08_OXI_F[C_OBS_LIST_U] <- lapply(C01_m08_OXI_F[C_OBS_LIST_U], SCR_U)
A01_m04_OXI_M[A_OBS_LIST_D] <- lapply(A01_m04_OXI_M[A_OBS_LIST_D], SCR_D)
A01_m04_OXI_F[A_OBS_LIST_D] <- lapply(A01_m04_OXI_F[A_OBS_LIST_D], SCR_D)
C01_m08_OXI_M[C_OBS_LIST_D] <- lapply(C01_m08_OXI_M[C_OBS_LIST_D], SCR_D)
C01_m08_OXI_F[C_OBS_LIST_D] <- lapply(C01_m08_OXI_F[C_OBS_LIST_D], SCR_D)
A01_m04_OXI_M$NA01_DR_TOT <- SCR_0(A01_m04_OXI_M$NA01_DR_TOT)
A01_m04_OXI_F$NA01_DR_TOT <- SCR_0(A01_m04_OXI_F$NA01_DR_TOT)

head(A01_m04_OXI_M); str(A01_m04_OXI_M)
head(A01_m04_OXI_F); str(A01_m04_OXI_F)
head(C01_m08_OXI_M); str(C01_m08_OXI_M)
head(C01_m08_OXI_F); str(C01_m08_OXI_F)

A01_m04_OXI_M$NA01_OBS <- rowSums(A01_m04_OXI_M); table(A01_m04_OXI_M$NA01_OBS)
A01_m04_OXI_F$NA01_OBS <- rowSums(A01_m04_OXI_F); table(A01_m04_OXI_F$NA01_OBS)
C01_m08_OXI_M$NC01_OBS <- rowSums(C01_m08_OXI_M); table(C01_m08_OXI_M$NC01_OBS)
C01_m08_OXI_F$NC01_OBS <- rowSums(C01_m08_OXI_F); table(C01_m08_OXI_F$NC01_OBS)
A01_m04_OXI_M <- dplyr::select(A01_m04_OXI_M,NA01_OBS)
A01_m04_OXI_F <- dplyr::select(A01_m04_OXI_F,NA01_OBS)
C01_m08_OXI_M <- dplyr::select(C01_m08_OXI_M,NC01_OBS)
C01_m08_OXI_F <- dplyr::select(C01_m08_OXI_F,NC01_OBS)

#CRp <- cor(t(m04bA01_OXI_M), A01_m04_OXI_M, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m04bA01_OXI_M), A01_m04_OXI_M, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(A01_m04_OXI_M)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m04bA01_OXI_M_CRPV <- CRp; m04bA01_OXI_M_CRPV$CRs <- CRs$CRs; m04bA01_OXI_M_CRPV$PVp <- PVp$PVp; m04bA01_OXI_M_CRPV$PVs <- PVs$PVs; m04bA01_OXI_M_CRPV$CLASS <- "m04bA01_OXI_M"
#CRp <- cor(t(m04bA01_OXI_F), A01_m04_OXI_F, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m04bA01_OXI_F), A01_m04_OXI_F, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(A01_m04_OXI_F)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m04bA01_OXI_F_CRPV <- CRp; m04bA01_OXI_F_CRPV$CRs <- CRs$CRs; m04bA01_OXI_F_CRPV$PVp <- PVp$PVp; m04bA01_OXI_F_CRPV$PVs <- PVs$PVs; m04bA01_OXI_F_CRPV$CLASS <- "m04bA01_OXI_F"
#CRp <- cor(t(m08bC01_OXI_M), C01_m08_OXI_M, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m08bC01_OXI_M), C01_m08_OXI_M, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(C01_m08_OXI_M)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m08bC01_OXI_M_CRPV <- CRp; m08bC01_OXI_M_CRPV$CRs <- CRs$CRs; m08bC01_OXI_M_CRPV$PVp <- PVp$PVp; m08bC01_OXI_M_CRPV$PVs <- PVs$PVs; m08bC01_OXI_M_CRPV$CLASS <- "m08bC01_OXI_M"
#CRp <- cor(t(m08bC01_OXI_F), C01_m08_OXI_F, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m08bC01_OXI_F), C01_m08_OXI_F, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(C01_m08_OXI_F)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m08bC01_OXI_F_CRPV <- CRp; m08bC01_OXI_F_CRPV$CRs <- CRs$CRs; m08bC01_OXI_F_CRPV$PVp <- PVp$PVp; m08bC01_OXI_F_CRPV$PVs <- PVs$PVs; m08bC01_OXI_F_CRPV$CLASS <- "m08bC01_OXI_F"
#CRPVb_Ox01_OXI <- rbind(m04bA01_OXI_M_CRPV,m04bA01_OXI_F_CRPV,m08bC01_OXI_M_CRPV,m08bC01_OXI_F_CRPV); dim(CRPVb_Ox01_OXI); table(CRPVb_Ox01_OXI$CLASS)
#write.table(CRPVb_Ox01_OXI, file="Output_Ox11_OXI/CRPVb_Ox01_OXI_OBS.txt",sep="\t",quote=FALSE,row.names=FALSE)
#nrow(CRPVb_Ox01_OXI); CRPVb_Ox01_OXI <- dplyr::filter(CRPVb_Ox01_OXI, CRp * CRs > 0); nrow(CRPVb_Ox01_OXI); CRPVb_Ox01_OXI <- dplyr::filter(CRPVb_Ox01_OXI, PVp < 0.01 & PVs < 0.01); nrow(CRPVb_Ox01_OXI) #1653644 > 1497298 > 20684
#CRPVb_Ox01_OXI$HyPi <- (CRPVb_Ox01_OXI$CRp * -log10(CRPVb_Ox01_OXI$PVp)) + (CRPVb_Ox01_OXI$CRs * -log10(CRPVb_Ox01_OXI$PVs))
#write.table(CRPVb_Ox01_OXI, file="Output_Ox11_OXI/CRPVb_Ox01_OXI_OBS_HyPi.txt",sep="\t",quote=FALSE,row.names=FALSE)
#
#CRp <- cor(t(m04mA01_OXI_M), A01_m04_OXI_M, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m04mA01_OXI_M), A01_m04_OXI_M, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(A01_m04_OXI_M)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m04mA01_OXI_M_CRPV <- CRp; m04mA01_OXI_M_CRPV$CRs <- CRs$CRs; m04mA01_OXI_M_CRPV$PVp <- PVp$PVp; m04mA01_OXI_M_CRPV$PVs <- PVs$PVs; m04mA01_OXI_M_CRPV$CLASS <- "m04mA01_OXI_M"
#CRp <- cor(t(m04mA01_OXI_F), A01_m04_OXI_F, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m04mA01_OXI_F), A01_m04_OXI_F, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(A01_m04_OXI_F)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m04mA01_OXI_F_CRPV <- CRp; m04mA01_OXI_F_CRPV$CRs <- CRs$CRs; m04mA01_OXI_F_CRPV$PVp <- PVp$PVp; m04mA01_OXI_F_CRPV$PVs <- PVs$PVs; m04mA01_OXI_F_CRPV$CLASS <- "m04mA01_OXI_F"
#CRp <- cor(t(m08mC01_OXI_M), C01_m08_OXI_M, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m08mC01_OXI_M), C01_m08_OXI_M, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(C01_m08_OXI_M)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m08mC01_OXI_M_CRPV <- CRp; m08mC01_OXI_M_CRPV$CRs <- CRs$CRs; m08mC01_OXI_M_CRPV$PVp <- PVp$PVp; m08mC01_OXI_M_CRPV$PVs <- PVs$PVs; m08mC01_OXI_M_CRPV$CLASS <- "m08mC01_OXI_M"
#CRp <- cor(t(m08mC01_OXI_F), C01_m08_OXI_F, method="pearson", use="pairwise.complete.obs"); CRs <- cor(t(m08mC01_OXI_F), C01_m08_OXI_F, method="spearman", use="pairwise.complete.obs")
#n_obs <- colSums(!is.na(C01_m08_OXI_F)); n_matrix <- matrix(rep(n_obs, each = nrow(CRp)), nrow = nrow(CRp)); t_pearson <- CRp * sqrt((n_matrix - 2) / (1 - CRp^2)); PVp <- 2 * pt(-abs(t_pearson), df = n_matrix - 2); t_spearman <- CRs * sqrt((n_matrix - 2) / (1 - CRs^2)); PVs <- 2 * pt(-abs(t_spearman), df = n_matrix - 2); CRp <- as.data.frame(CRp); CRp$ID <- rownames(CRp); CRp <- as.data.frame(tidyr::pivot_longer(data=CRp, cols=-ID, names_to="Var", values_to="CRp")); CRs <- as.data.frame(CRs); CRs$ID <- rownames(CRs); CRs <- as.data.frame(tidyr::pivot_longer(data=CRs, cols=-ID, names_to="Var", values_to="CRs")); PVp <- as.data.frame(PVp); PVp$ID <- rownames(PVp); PVp <- as.data.frame(tidyr::pivot_longer(data=PVp, cols=-ID, names_to="Var", values_to="PVp")); PVs <- as.data.frame(PVs); PVs$ID <- rownames(PVs); PVs <- as.data.frame(tidyr::pivot_longer(data=PVs, cols=-ID, names_to="Var", values_to="PVs"))
#m08mC01_OXI_F_CRPV <- CRp; m08mC01_OXI_F_CRPV$CRs <- CRs$CRs; m08mC01_OXI_F_CRPV$PVp <- PVp$PVp; m08mC01_OXI_F_CRPV$PVs <- PVs$PVs; m08mC01_OXI_F_CRPV$CLASS <- "m08mC01_OXI_F"
#CRPVm_Ox01_OXI <- rbind(m04mA01_OXI_M_CRPV,m04mA01_OXI_F_CRPV,m08mC01_OXI_M_CRPV,m08mC01_OXI_F_CRPV); dim(CRPVm_Ox01_OXI); table(CRPVm_Ox01_OXI$CLASS)
#write.table(CRPVm_Ox01_OXI, file="Output_Ox11_OXI/CRPVm_Ox01_OXI_OBS.txt",sep="\t",quote=FALSE,row.names=FALSE)
#nrow(CRPVm_Ox01_OXI); CRPVm_Ox01_OXI <- dplyr::filter(CRPVm_Ox01_OXI, CRp * CRs > 0); nrow(CRPVm_Ox01_OXI); CRPVm_Ox01_OXI <- dplyr::filter(CRPVm_Ox01_OXI, PVp < 0.01 & PVs < 0.01); nrow(CRPVm_Ox01_OXI) #1653644 > 1521765 > 25715
#CRPVm_Ox01_OXI$HyPi <- (CRPVm_Ox01_OXI$CRp * -log10(CRPVm_Ox01_OXI$PVp)) + (CRPVm_Ox01_OXI$CRs * -log10(CRPVm_Ox01_OXI$PVs))
#write.table(CRPVm_Ox01_OXI, file="Output_Ox11_OXI/CRPVm_Ox01_OXI_OBS_HyPi.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Ox01_OXI_CRPVb_Raw <- as.data.frame(fread("Output_Ox11_OXI/CRPVb_Ox01_OXI_OBS.txt.gz",  header = TRUE, sep = "\t"))
#Ox01_OXI_CRPVb_HPx <- as.data.frame(fread("Output_Ox11_OXI/CRPVb_Ox01_OXI_OBS_HyPi.txt.gz",  header = TRUE, sep = "\t"))
#Ox01_OXI_CRPVm_Raw <- as.data.frame(fread("Output_Ox11_OXI/CRPVm_Ox01_OXI_OBS.txt.gz",  header = TRUE, sep = "\t"))
#Ox01_OXI_CRPVm_HPx <- as.data.frame(fread("Output_Ox11_OXI/CRPVm_Ox01_OXI_OBS_HyPi.txt.gz",  header = TRUE, sep = "\t"))
#head(Ox01_OXI_CRPVb_Raw); dim(Ox01_OXI_CRPVb_Raw) #24804660
#head(Ox01_OXI_CRPVb_HPx); dim(Ox01_OXI_CRPVb_HPx) #  644500
#head(Ox01_OXI_CRPVm_Raw); dim(Ox01_OXI_CRPVm_Raw) #24804660
#head(Ox01_OXI_CRPVm_HPx); dim(Ox01_OXI_CRPVm_HPx) #  686098



#Select CpGs: intersect of KoGES & GSE223748
m04bA01_OXI_M <- dplyr::select(Met_ASAS0400b[rownames(Met_ASAS0400b) %in% INT_kex_CpG, ], intersect(colnames(Met_ASAS0400b), dplyr::filter(ASAS01to10, AS01_SEX==1)$DIST_ID)); m04bA01_OXI_M <- m04bA01_OXI_M[complete.cases(m04bA01_OXI_M), ]; m04bA01_OXI_M <- as.data.frame(as.data.table(m04bA01_OXI_M, keep.rownames = "ProbeID")); rownames(m04bA01_OXI_M) <- m04bA01_OXI_M$ProbeID; dim(m04bA01_OXI_M) #5306,201
m04bA01_OXI_F <- dplyr::select(Met_ASAS0400b[rownames(Met_ASAS0400b) %in% INT_kex_CpG, ], intersect(colnames(Met_ASAS0400b), dplyr::filter(ASAS01to10, AS01_SEX==2)$DIST_ID)); m04bA01_OXI_F <- m04bA01_OXI_F[complete.cases(m04bA01_OXI_F), ]; m04bA01_OXI_F <- as.data.frame(as.data.table(m04bA01_OXI_F, keep.rownames = "ProbeID")); rownames(m04bA01_OXI_F) <- m04bA01_OXI_F$ProbeID; dim(m04bA01_OXI_F) #5306,201
m08bC01_OXI_M <- dplyr::select(Met_CITY0822b[rownames(Met_CITY0822b) %in% INT_kex_CpG, ], intersect(colnames(Met_CITY0822b), dplyr::filter(CITY01to02, CT01_SEX==1)$DIST_ID)); m08bC01_OXI_M <- m08bC01_OXI_M[complete.cases(m08bC01_OXI_M), ]; m08bC01_OXI_M <- as.data.frame(as.data.table(m08bC01_OXI_M, keep.rownames = "ProbeID")); rownames(m08bC01_OXI_M) <- m08bC01_OXI_M$ProbeID; dim(m08bC01_OXI_M) #5306,623
m08bC01_OXI_F <- dplyr::select(Met_CITY0822b[rownames(Met_CITY0822b) %in% INT_kex_CpG, ], intersect(colnames(Met_CITY0822b), dplyr::filter(CITY01to02, CT01_SEX==2)$DIST_ID)); m08bC01_OXI_F <- m08bC01_OXI_F[complete.cases(m08bC01_OXI_F), ]; m08bC01_OXI_F <- as.data.frame(as.data.table(m08bC01_OXI_F, keep.rownames = "ProbeID")); rownames(m08bC01_OXI_F) <- m08bC01_OXI_F$ProbeID; dim(m08bC01_OXI_F) #5306,201
m04bA01_OXI_M <- m04bA01_OXI_M[order(m04bA01_OXI_M$ProbeID), ]
m04bA01_OXI_F <- m04bA01_OXI_F[order(m04bA01_OXI_F$ProbeID), ]; table(m04bA01_OXI_M$ProbeID==m04bA01_OXI_F$ProbeID) #5306
m08bC01_OXI_M <- m08bC01_OXI_M[order(m08bC01_OXI_M$ProbeID), ]; table(m04bA01_OXI_M$ProbeID==m08bC01_OXI_M$ProbeID) #5306
m08bC01_OXI_F <- m08bC01_OXI_F[order(m08bC01_OXI_F$ProbeID), ]; table(m04bA01_OXI_M$ProbeID==m08bC01_OXI_F$ProbeID) #5306
m04bA01_OXI_M[1:5, 1:5]; dim(m04bA01_OXI_M) #5213,196
m04bA01_OXI_F[1:5, 1:5]; dim(m04bA01_OXI_F) #5213,197
m08bC01_OXI_M[1:5, 1:5]; dim(m08bC01_OXI_M) #5213,617
m08bC01_OXI_F[1:5, 1:5]; dim(m08bC01_OXI_F) #5213,197

CRPV_OBS <- as.data.frame(fread("Output_Ox11_OXI/CRPVb_Ox01_OXI_OBS.txt.gz",  header = TRUE, sep = "\t")); CRPV_OBS <- dplyr::filter(CRPV_OBS, CRp * CRs > 0); dim(CRPV_OBS)
CRPV_OBS$adj_PVp <- p.adjust(CRPV_OBS$PVp, method = "fdr"); CRPV_OBS$adj_PVs <- p.adjust(CRPV_OBS$PVs, method = "fdr")
CRPV_OBSx <- dplyr::filter(CRPV_OBS, abs(CRp) > 0.24 & abs(CRs) > 0.24 & PVp < 0.002 & PVs < 0.002 & substr(CLASS,5,7)=="A01"); table(CRPV_OBSx$CLASS); nrow(CRPV_OBSx) #871
CRPV_OBSx$HyPi <- (CRPV_OBSx$CRp * -log10(CRPV_OBSx$PVp)) + (CRPV_OBSx$CRs * -log10(CRPV_OBSx$PVs))
CRPV_OBSm <- dplyr::filter(CRPV_OBSx, substr(CLASS,13,13)=="M"); CRPV_OBSf <- dplyr::filter(CRPV_OBSx, substr(CLASS,13,13)=="F"); intersect(CRPV_OBSm$ID,CRPV_OBSf$ID)
CRPV_OBSx$Sex <- substr(CRPV_OBSx$CLASS,13,13)

m04bA01_OXI_M_OBS <- as.data.frame(t(dplyr::select(dplyr::filter(m04bA01_OXI_M, ProbeID %in% CRPV_OBSm$ID), -ProbeID))); dim(m04bA01_OXI_M_OBS) #195,5
m04bA01_OXI_F_OBS <- as.data.frame(t(dplyr::select(dplyr::filter(m04bA01_OXI_F, ProbeID %in% CRPV_OBSf$ID), -ProbeID))); dim(m04bA01_OXI_F_OBS) #196,9
m04bA01_OXI_M_OBS <- merge(x=A01_m04_OXI_M, y=m04bA01_OXI_M_OBS,by=0); rownames(m04bA01_OXI_M_OBS) <- m04bA01_OXI_M_OBS$Row.names; m04bA01_OXI_M_OBS <- m04bA01_OXI_M_OBS[,-1]
m04bA01_OXI_F_OBS <- merge(x=A01_m04_OXI_F, y=m04bA01_OXI_F_OBS,by=0); rownames(m04bA01_OXI_F_OBS) <- m04bA01_OXI_F_OBS$Row.names; m04bA01_OXI_F_OBS <- m04bA01_OXI_F_OBS[,-1]

head(m04bA01_OXI_M_OBS); head(m04bA01_OXI_F_OBS)

head(GPL13534); dim(GPL13534) #485577,7
head(CRPV_OBS); dim(CRPV_OBS) #1497298,9
length(intersect(CRPV_OBS$ID,GPL13534$ID)) #413321
head(CRPV_OBSx); dim(CRPV_OBSx) #871,9
length(intersect(CRPV_OBSx$ID,GPL13534$ID)) #871
m04bA01_OXI_M[1:5,1:5]; dim(m04bA01_OXI_M)
m04bA01_OXI_F[1:5,1:5]; dim(m04bA01_OXI_F)


GPL13534x <- GPL13534
GPL13534x$SIG <- ifelse(GPL13534x$ID %in% colnames(m04bA01_OXI_M_OBS), "M", ifelse(GPL13534x$ID %in% colnames(m04bA01_OXI_F_OBS), "F","XXX"))
GPL13534x <- merge(x=GPL13534x,y=CRPV_OBSx[,c(1,10)])
GPL13534x <- dplyr::filter(GPL13534x, SIG != "XXX")
table(GPL13534x$SIG)
GPL13534x$UCSC_RefGene_Name <- sub(";.*", "", GPL13534x$UCSC_RefGene_Name); GPL13534x$UCSC_RefGene_Group <- sub(";.*", "", GPL13534x$UCSC_RefGene_Group); GPL13534x$HyPi <- round(as.numeric(GPL13534x$HyPi),2)
head(GPL13534x)
write.table(GPL13534x, file="Output_Ox11_OXI/01_CpG_Table.txt", sep="\t",quote=FALSE,row.names=FALSE)

#Plot 01: CpG-OBS correlation plots for 14 CpG sites
pS01 <- ggplot(m04bA01_OXI_M_OBS, aes(x = NA01_OBS, y = cg02132051)) + geom_point(color = "steelblue", size = 2) + geom_smooth(method = "lm", color =    "red", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg02132051\n(RGS6, Body)        ")
pS02 <- ggplot(m04bA01_OXI_M_OBS, aes(x = NA01_OBS, y = cg02499843)) + geom_point(color = "steelblue", size = 2) + geom_smooth(method = "lm", color =    "red", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg02499843\n(PLXNA4, 1stExon)   ")
pS03 <- ggplot(m04bA01_OXI_M_OBS, aes(x = NA01_OBS, y = cg04388989)) + geom_point(color = "steelblue", size = 2) + geom_smooth(method = "lm", color =    "red", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg04388989\n                    ")
pS04 <- ggplot(m04bA01_OXI_M_OBS, aes(x = NA01_OBS, y = cg05486213)) + geom_point(color = "steelblue", size = 2) + geom_smooth(method = "lm", color =    "red", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg05486213\n(MPPED2, TSS200)    ")
pS05 <- ggplot(m04bA01_OXI_M_OBS, aes(x = NA01_OBS, y = cg16802592)) + geom_point(color = "steelblue", size = 2) + geom_smooth(method = "lm", color =    "red", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg16802592\n(LOC145845, TSS1500)")
pS06 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg02650266)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg02650266\n                    ")
pS07 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg04218899)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg04218899\n(LAMC2, 1stExon)    ")
pS08 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg06458239)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg06458239\n(ZNF549, TSS200)    ")
pS09 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg06635150)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg06635150\n                    ")
pS10 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg11476211)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg11476211\n(PRKCE, 1stExon)    ")
pS11 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg16434657)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg16434657\n                    ")
pS12 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg20550118)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg20550118\n(CRABP1, Body)      ")
pS13 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg25925945)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg25925945\n(NCKAP1, 1stExon)   ")
pS14 <- ggplot(m04bA01_OXI_F_OBS, aes(x = NA01_OBS, y = cg26880891)) + geom_point(color = "darkgreen", size = 2) + geom_smooth(method = "lm", color = "orange", se = TRUE) + theme_minimal() + labs(title = "Scatter Plot: OBS vs cg26880891\n(ARIH1, 1stExon)    ")
ggsave("Output_Ox11_OXI/01_CpG_OBS_Corr.pdf", plot = grid.arrange(pS01,pS02,pS03,pS04,pS05,pS06,pS07,pS08,pS09,pS10,pS11,pS12,pS13,pS14, ncol=2), width = 10, height = 21)

CRPV_Suba <- CRPV_OBSx[, c("ID", "CRp", "CRs", "Sex")]; CRPV_Suba <- dplyr::filter(CRPV_Suba, ID %in% c(colnames(m04bA01_OXI_M_OBS),colnames(m04bA01_OXI_F_OBS))); dim(CRPV_Suba); table(CRPV_Suba$Sex)
Wide_CRPa <- reshape(CRPV_Suba, idvar = "ID", timevar = "Sex", direction = "wide"); rownames(Wide_CRPa) <- Wide_CRPa$ID; Wide_CRPa <- Wide_CRPa[, -1]; Wide_CRPa[is.na(Wide_CRPa)] <- 0
Wide_CRPa_M <- dplyr::select(Wide_CRPa, contains(".M")); ncol_Wide_CRPa_M <- ncol(Wide_CRPa_M)
Wide_CRPa_F <- dplyr::select(Wide_CRPa, contains(".F")); ncol_Wide_CRPa_F <- ncol(Wide_CRPa_F)
Wide_CRPa_M <- Wide_CRPa_M[rowSums(Wide_CRPa_M != 0) != 0, ]
Wide_CRPa_F <- Wide_CRPa_F[rowSums(Wide_CRPa_F != 0) != 0, ]
head(Wide_CRPa_M); dim(Wide_CRPa_M) #5,2
head(Wide_CRPa_F); dim(Wide_CRPa_F) #9,2

#Plot 02: Age-correlation CpG sites in five datasets, male and female
pdf("Output_Ox11_OXI/02_KoGES_Age_Corr_M.pdf", width = 3, height = 3)
pheatmap(Wide_CRPa_M, main = "OBS Correlation in Male",
         cluster_cols = FALSE, cluster_rows = TRUE, display_numbers = TRUE, number_color = "black", fontsize_number = 6,color = colorRampPalette(c("#1581BF", "white", "#BF124D"))(100), breaks = seq(-1, 1, length.out = 101)); dev.off()
pdf("Output_Ox11_OXI/02_KoGES_Age_Corr_F.pdf", width = 3, height = 3)
pheatmap(Wide_CRPa_F, main = "OBS Correlation in Female",
         cluster_cols = FALSE, cluster_rows = TRUE, display_numbers = TRUE, number_color = "black", fontsize_number = 6, color = colorRampPalette(c("#1581BF", "white", "#BF124D"))(100), breaks = seq(-1, 1, length.out = 101)); dev.off()


GSE223748m$ID <- rownames(GSE223748m); GSE223748m <- dplyr::filter(GSE223748m, ID %in% m04bA01_OXI_M$ProbeID); GSE223748m$ID <- NULL
GSE223748m[1:5, 1:5]; dim(GSE223748m) #37554->5213,15043
GSE223748i[1:5, 1:5]; dim(GSE223748i) #15043,12
GSE223748i$Sex <- as.numeric(GSE223748i$Sex); table(GSE223748i$Sex)
GSE223748i <- dplyr::filter(GSE223748i, Sex != 3)
head(GSE223748i); dim(GSE223748i); table(GSE223748i$Sex) #14667,12 / 1 6842, 2 7825

Blood_SciNameM <- table(dplyr::filter(GSE223748i, Tissue=="Blood"&Sex==1)$SciName); Blood_SciNameM <- Blood_SciNameM[Blood_SciNameM > 49]; Blood_SciNameM
Blood_SciNameF <- table(dplyr::filter(GSE223748i, Tissue=="Blood"&Sex==2)$SciName); Blood_SciNameF <- Blood_SciNameF[Blood_SciNameF > 49]; Blood_SciNameF
COMMON_Blood_SciName <- intersect(names(Blood_SciNameM),names(Blood_SciNameF)); length(COMMON_Blood_SciName) #9

GSE223748iM <- dplyr::filter(GSE223748i, Sex==1 & Tissue=="Blood" & SciName %in% COMMON_Blood_SciName)[,2:3]
GSE223748iF <- dplyr::filter(GSE223748i, Sex==2 & Tissue=="Blood" & SciName %in% COMMON_Blood_SciName)[,2:3]
GSE223748mM <- GSE223748m; GSE223748mM$ID <- rownames(GSE223748mM); GSE223748mM <- dplyr::select(dplyr::filter(GSE223748mM, ID %in% rownames(Wide_CRPa_M)), rownames(GSE223748iM))
GSE223748mF <- GSE223748m; GSE223748mF$ID <- rownames(GSE223748mF); GSE223748mF <- dplyr::select(dplyr::filter(GSE223748mF, ID %in% rownames(Wide_CRPa_F)), rownames(GSE223748iF))
GSE223748mM[1:5,1:5]; dim(GSE223748mM) #5,1257
GSE223748mF[1:5,1:5]; dim(GSE223748mF) #9,1252
head(GSE223748iM); dim(GSE223748iM) #1257,2
head(GSE223748iF); dim(GSE223748iF) #1252,2
table(colnames(GSE223748mM)==rownames(GSE223748iM)) #1257
table(colnames(GSE223748mF)==rownames(GSE223748iF)) #1252

GSE223748mM_Pr <- data.frame(row.names = rownames(GSE223748mM)); GSE223748mM_Sp <- data.frame(row.names = rownames(GSE223748mM))
GSE223748mF_Pr <- data.frame(row.names = rownames(GSE223748mF)); GSE223748mF_Sp <- data.frame(row.names = rownames(GSE223748mF))

for(spec in COMMON_Blood_SciName) {
  idx_Ma <- which(GSE223748iM$SciName == spec); idx_Fa <- which(GSE223748iF$SciName == spec)
  if(length(idx_Ma) > 2) {
    GSE223748mM_Pr[, spec] <- cor(t(GSE223748mM[, idx_Ma]), GSE223748iM$Age[idx_Ma], use="pairwise.complete.obs", method="pearson")
    GSE223748mM_Sp[, spec] <- cor(t(GSE223748mM[, idx_Ma]), GSE223748iM$Age[idx_Ma], use="pairwise.complete.obs", method="spearman")
  } else { GSE223748mM_Pr[, spec] <- NA; GSE223748mM_Sp[, spec] <- NA }
  if(length(idx_Fa) > 2) {
    GSE223748mF_Pr[, spec] <- cor(t(GSE223748mF[, idx_Fa]), GSE223748iF$Age[idx_Fa], use="pairwise.complete.obs", method="pearson")
    GSE223748mF_Sp[, spec] <- cor(t(GSE223748mF[, idx_Fa]), GSE223748iF$Age[idx_Fa], use="pairwise.complete.obs", method="spearman")
  } else { GSE223748mF_Pr[, spec] <- NA; GSE223748mF_Sp[, spec] <- NA }
}
GSE223748mM_Pr[is.na(GSE223748mM_Pr)] <- 0; GSE223748mM_Sp[is.na(GSE223748mM_Sp)] <- 0
GSE223748mF_Pr[is.na(GSE223748mF_Pr)] <- 0; GSE223748mF_Sp[is.na(GSE223748mF_Sp)] <- 0
dim(GSE223748mM_Pr); dim(GSE223748mM_Sp); dim(GSE223748mF_Pr); dim(GSE223748mF_Sp)

#Plot 03: Species Clustering Heatmap (Pearson vs Spearman)
pdf("Output_Ox11_OXI/02_Species_M_P.pdf", width = 4, height = 5)
pheatmap(GSE223748mM_Pr, main = "Male: Pearson Correlation, \nbetween OBS and beta value",    cluster_cols = TRUE, cluster_rows = TRUE, color = colorRampPalette(c("#4E61D3", "white", "#BF124D"))(100), breaks = seq(-1, 1, length.out = 101), fontsize_row = 8, fontsize_col = 9); dev.off()
pdf("Output_Ox11_OXI/02_Species_M_S.pdf", width = 4, height = 5)
pheatmap(GSE223748mM_Sp, main = "Male: Spearman Correlation, \nbetween OBS and beta value",   cluster_cols = TRUE, cluster_rows = TRUE, color = colorRampPalette(c("#4E61D3", "white", "#BF124D"))(100), breaks = seq(-1, 1, length.out = 101), fontsize_row = 8, fontsize_col = 9); dev.off()
pdf("Output_Ox11_OXI/02_Species_F_P.pdf", width = 4, height = 5)
pheatmap(GSE223748mF_Pr, main = "Female: Pearson Correlation, \nbetween OBS and beta value",  cluster_cols = TRUE, cluster_rows = TRUE, color = colorRampPalette(c("#4E61D3", "white", "#BF124D"))(100), breaks = seq(-1, 1, length.out = 101), fontsize_row = 8, fontsize_col = 9); dev.off()
pdf("Output_Ox11_OXI/02_Species_F_S.pdf", width = 4, height = 5)
pheatmap(GSE223748mF_Sp, main = "Female: Spearman Correlation, \nbetween OBS and beta value", cluster_cols = TRUE, cluster_rows = TRUE, color = colorRampPalette(c("#4E61D3", "white", "#BF124D"))(100), breaks = seq(-1, 1, length.out = 101), fontsize_row = 8, fontsize_col = 9); dev.off()

Long_Male_Pra <- melt(as.matrix(GSE223748mM_Pr)); colnames(Long_Male_Pra) <- c("CpG", "Species", "Pearson")
Long_Male_Spa <- melt(as.matrix(GSE223748mM_Sp)); colnames(Long_Male_Spa) <- c("CpG", "Species", "Spearman")
Combined_Malea <- merge(Long_Male_Pra, Long_Male_Spa, by = c("CpG", "Species")); rm(Long_Male_Pra); rm(Long_Male_Spa)
Combined_Malea$CLASS <- "M"

Long_Female_Pra <- melt(as.matrix(GSE223748mF_Pr)); colnames(Long_Female_Pra) <- c("CpG", "Species", "Pearson")
Long_Female_Spa <- melt(as.matrix(GSE223748mF_Sp)); colnames(Long_Female_Spa) <- c("CpG", "Species", "Spearman")
Combined_Femalea <- merge(Long_Female_Pra, Long_Female_Spa, by = c("CpG", "Species")); rm(Long_Female_Pra); rm(Long_Female_Spa)
Combined_Femalea$CLASS <- "F"


head(Combined_Malea); dim(Combined_Malea); table(Combined_Malea$CpG); length(unique(Combined_Malea$CpG)) #252,4 | 28*9=252
head(Combined_Femalea); dim(Combined_Femalea); table(Combined_Femalea$CpG) #252,4 | 28*9=252
length(intersect(unique(Combined_Malea$CpG),unique(Combined_Femalea$CpG)))
length(union(unique(Combined_Malea$CpG),unique(Combined_Femalea$CpG)))

Ox11_MF <- rbind(Combined_Malea,Combined_Femalea); colnames(Ox11_MF)[1] <- "ID"
Ox11_MF <- merge(x=Ox11_MF,y=GPL13534,by="ID",all.x=TRUE)
head(Ox11_MF); dim(Ox11_MF) #126,11
table(Ox11_MF[,c(2,1)])


Ox11_MFx <- unique(Ox11_MF[,-2:-4]); dim(Ox11_MFx)
Ox11_MF[] <- lapply(Ox11_MF, function(x) sub(";.*", "", x)); Ox11_MF$Pearson <- round(as.numeric(Ox11_MF$Pearson),2); Ox11_MF$Spearman <- round(as.numeric(Ox11_MF$Spearman),2)

head(Ox11_MFx); dim(Ox11_MFx); table(Ox11_MFx$ID)
write.table(Ox11_MF, file="Output_Ox11_OXI/03_STable.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(Ox11_MFx, file="Output_Ox11_OXI/03_Table.txt",sep="\t",quote=FALSE,row.names=FALSE)

#Plot 04: Linearity Check
Plot_Male_Linearity <- ggplot(Combined_Malea, aes(x = Pearson, y = Spearman)) +
  geom_point(alpha = 0.6, color = "dodgerblue") + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_wrap(~ Species, ncol = 3) + theme_bw() +
  labs(title = "Male Linearity Check: Pearson vs Spearman", x = "Pearson (r)", y = "Spearman (rho)")
Plot_Female_Linearity <- ggplot(Combined_Femalea, aes(x = Pearson, y = Spearman)) +
  geom_point(alpha = 0.6, color = "deeppink") + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  facet_wrap(~ Species, ncol = 3) + theme_bw() +
  labs(title = "Female Linearity Check: Pearson vs Spearman", x = "Pearson (r)", y = "Spearman (rho)")
ggsave("Output_Ox11_OXI/04_Linearity.pdf", plot=grid.arrange(Plot_Male_Linearity,Plot_Female_Linearity, ncol=2), width=12, height=6)


Ref_Order_a <- c("Homo sapiens", "Macaca mulatta", "Chlorocebus sabaeus", "Callithrix jacchus", "Mus musculus", "Rattus norvegicus", "Canis lupus familiaris", "Felis catus", "Bos taurus", "Sus scrofa", "Ovis aries", "Equus caballus", "Capreolus capreolus")
Evo_Order_a <- intersect(Ref_Order_a, colnames(GSE223748mM_Pr))

Long_M_Pr_a <- melt(as.matrix(GSE223748mM_Pr)); colnames(Long_M_Pr_a) <- c("CpG", "Species", "Cor"); Long_M_Pr_a$Method <- "Pearson"; Long_M_Pr_a$Sex <- "Male"
Long_M_Sp_a <- melt(as.matrix(GSE223748mM_Sp)); colnames(Long_M_Sp_a) <- c("CpG", "Species", "Cor"); Long_M_Sp_a$Method <- "Spearman"; Long_M_Sp_a$Sex <- "Male"
Long_F_Pr_a <- melt(as.matrix(GSE223748mF_Pr)); colnames(Long_F_Pr_a) <- c("CpG", "Species", "Cor"); Long_F_Pr_a$Method <- "Pearson"; Long_F_Pr_a$Sex <- "Female"
Long_F_Sp_a <- melt(as.matrix(GSE223748mF_Sp)); colnames(Long_F_Sp_a) <- c("CpG", "Species", "Cor"); Long_F_Sp_a$Method <- "Spearman"; Long_F_Sp_a$Sex <- "Female"
Long_M_Pr_a$Species <- factor(Long_M_Pr_a$Species, levels = Evo_Order_a); Long_M_Sp_a$Species <- factor(Long_M_Sp_a$Species, levels = Evo_Order_a)
Long_F_Pr_a$Species <- factor(Long_F_Pr_a$Species, levels = Evo_Order_a); Long_F_Sp_a$Species <- factor(Long_F_Sp_a$Species, levels = Evo_Order_a)
Dual_Dat_M_a <- rbind(Long_M_Pr_a, Long_M_Sp_a); Dual_Dat_F_a <- rbind(Long_F_Pr_a, Long_F_Sp_a)

#Plot 05: Dual-Phylo Decay Plot: Comparison of the decrease in absolute value of correlation coefficients according to evolutionary distance (Pearson vs Spearman)
Plot_Dual_M_a <- ggplot(Dual_Dat_M_a, aes(x = Species, y = abs(Cor), color = Method, group = Method)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) + stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_manual(values = c("Pearson" = "dodgerblue", "Spearman" = "firebrick")) +
  labs(title = "Dual-Phylo Decay: Male (Pearson vs Spearman)", y = "Mean |Correlation|", x = "") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
Plot_Dual_F_a <- ggplot(Dual_Dat_F_a, aes(x = Species, y = abs(Cor), color = Method, group = Method)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) + stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_manual(values = c("Pearson" = "dodgerblue", "Spearman" = "firebrick")) +
  labs(title = "Dual-Phylo Decay: Female (Pearson vs Spearman)", y = "Mean |Correlation|", x = "") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Output_Ox11_OXI/05_DecayPlot_Linegraph_twomethods.pdf", plot=grid.arrange(Plot_Dual_M_a, Plot_Dual_F_a, ncol=2), width=12, height=6)


#Plot 06: Phylo-Correlation Decay Plot (Boxplot with Jitter)
Plot_Phylo_M_Pr_a <- ggplot(Long_M_Pr_a, aes(x = Species, y = Cor)) +
  geom_boxplot(outlier.shape = NA, fill="aliceblue") + geom_jitter(width = 0.2, alpha = 0.5, color="navy") +
  geom_hline(yintercept = 0, linetype="dotted") +
  labs(title = "Phylo-Correlation Decay: Male (Pearson)", x = "", y = "Correlation") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
Plot_Phylo_M_Sp_a <- ggplot(Long_M_Sp_a, aes(x = Species, y = Cor)) +
  geom_boxplot(outlier.shape = NA, fill="aliceblue") + geom_jitter(width = 0.2, alpha = 0.5, color="navy") +
  geom_hline(yintercept = 0, linetype="dotted") +
  labs(title = "Phylo-Correlation Decay: Male (Spearman)", x = "", y = "Correlation") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
Plot_Phylo_F_Pr_a <- ggplot(Long_F_Pr_a, aes(x = Species, y = Cor)) +
  geom_boxplot(outlier.shape = NA, fill="mistyrose") + geom_jitter(width = 0.2, alpha = 0.5, color="darkred") +
  geom_hline(yintercept = 0, linetype="dotted") +
  labs(title = "Phylo-Correlation Decay: Female (Pearson)", x = "", y = "Correlation") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
Plot_Phylo_F_Sp_a <- ggplot(Long_F_Sp_a, aes(x = Species, y = Cor)) +
  geom_boxplot(outlier.shape = NA, fill="mistyrose") + geom_jitter(width = 0.2, alpha = 0.5, color="darkred") +
  geom_hline(yintercept = 0, linetype="dotted") +
  labs(title = "Phylo-Correlation Decay: Female (Spearman)", x = "", y = "Correlation") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Output_Ox11_OXI/06_DecayPlot_Boxplot.pdf", plot=grid.arrange(Plot_Phylo_M_Pr_a, Plot_Phylo_M_Sp_a, Plot_Phylo_F_Pr_a, Plot_Phylo_F_Sp_a, ncol=2), width=12, height=8)


#Plot 07: Consensus vs Divergence Plot (X=Homo sapiens, Y=Other Species)
Human_Vals_M_Pr <- GSE223748mM_Pr[, "Homo sapiens", drop=FALSE]; colnames(Human_Vals_M_Pr) <- "Human_Cor"; Human_Vals_M_Pr$CpG <- rownames(Human_Vals_M_Pr)
Human_Vals_M_Sp <- GSE223748mM_Sp[, "Homo sapiens", drop=FALSE]; colnames(Human_Vals_M_Sp) <- "Human_Cor"; Human_Vals_M_Sp$CpG <- rownames(Human_Vals_M_Sp)
Human_Vals_F_Pr <- GSE223748mF_Pr[, "Homo sapiens", drop=FALSE]; colnames(Human_Vals_F_Pr) <- "Human_Cor"; Human_Vals_F_Pr$CpG <- rownames(Human_Vals_F_Pr)
Human_Vals_F_Sp <- GSE223748mF_Sp[, "Homo sapiens", drop=FALSE]; colnames(Human_Vals_F_Sp) <- "Human_Cor"; Human_Vals_F_Sp$CpG <- rownames(Human_Vals_F_Sp)
Consensus_Dat_M_Pr <- merge(Long_M_Pr_a, Human_Vals_M_Pr, by="CpG"); Consensus_Dat_M_Pr <- Consensus_Dat_M_Pr[Consensus_Dat_M_Pr$Species != "Homo sapiens", ]
Consensus_Dat_M_Sp <- merge(Long_M_Sp_a, Human_Vals_M_Sp, by="CpG"); Consensus_Dat_M_Sp <- Consensus_Dat_M_Sp[Consensus_Dat_M_Sp$Species != "Homo sapiens", ]
Consensus_Dat_F_Pr <- merge(Long_F_Pr_a, Human_Vals_F_Pr, by="CpG"); Consensus_Dat_F_Pr <- Consensus_Dat_F_Pr[Consensus_Dat_F_Pr$Species != "Homo sapiens", ]
Consensus_Dat_F_Sp <- merge(Long_F_Sp_a, Human_Vals_F_Sp, by="CpG"); Consensus_Dat_F_Sp <- Consensus_Dat_F_Sp[Consensus_Dat_F_Sp$Species != "Homo sapiens", ]

Plot_Consens_M_Pr <- ggplot(Consensus_Dat_M_Pr, aes(x=Human_Cor, y=Cor)) + geom_point(alpha=0.6, color="darkgreen") + geom_smooth(method="lm", se=FALSE, color="grey", linetype="dashed") + geom_abline(intercept=0, slope=1, color="red", linetype="dotted") + facet_wrap(~Species, ncol=4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + labs(title="Consensus vs Divergence: Male Pearson (Ref: Homo sapiens)", x="Human Correlation", y="Species Correlation") + theme_bw()
Plot_Consens_M_Sp <- ggplot(Consensus_Dat_M_Sp, aes(x=Human_Cor, y=Cor)) + geom_point(alpha=0.6, color="darkgreen") + geom_smooth(method="lm", se=FALSE, color="grey", linetype="dashed") + geom_abline(intercept=0, slope=1, color="red", linetype="dotted") + facet_wrap(~Species, ncol=4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + labs(title="Consensus vs Divergence: Male Spearman (Ref: Homo sapiens)", x="Human Correlation", y="Species Correlation") + theme_bw()
Plot_Consens_F_Pr <- ggplot(Consensus_Dat_F_Pr, aes(x=Human_Cor, y=Cor)) + geom_point(alpha=0.6, color="darkgreen") + geom_smooth(method="lm", se=FALSE, color="grey", linetype="dashed") + geom_abline(intercept=0, slope=1, color="red", linetype="dotted") + facet_wrap(~Species, ncol=4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + labs(title="Consensus vs Divergence: Female Pearson (Ref: Homo sapiens)", x="Human Correlation", y="Species Correlation") + theme_bw()
Plot_Consens_F_Sp <- ggplot(Consensus_Dat_F_Sp, aes(x=Human_Cor, y=Cor)) + geom_point(alpha=0.6, color="darkgreen") + geom_smooth(method="lm", se=FALSE, color="grey", linetype="dashed") + geom_abline(intercept=0, slope=1, color="red", linetype="dotted") + facet_wrap(~Species, ncol=4) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + labs(title="Consensus vs Divergence: Female Spearman (Ref: Homo sapiens)", x="Human Correlation", y="Species Correlation") + theme_bw()
ggsave("Output_Ox11_OXI/07_ConsensPlot.pdf", plot=grid.arrange(Plot_Consens_M_Pr, Plot_Consens_M_Sp, Plot_Consens_F_Pr, Plot_Consens_F_Sp, ncol=2), width=12, height=8)




#Plot 08: Epigenetic Distance Matrix by Euclidean
table(colnames(GSE223748mM_Pr)==colnames(GSE223748mF_Pr))
table(colnames(GSE223748mM_Pr)==colnames(GSE223748mM_Sp))
table(colnames(GSE223748mM_Pr)==colnames(GSE223748mF_Sp))
Sp_COL_NAMES <- c("C. familiaris","C. capreolus","E. caballus","F. catus","H. sapiens","M. mulatta","M. musculus","O. aries","R. norvegicus")

colnames(GSE223748mM_Pr) <- Sp_COL_NAMES
colnames(GSE223748mF_Pr) <- Sp_COL_NAMES
colnames(GSE223748mM_Sp) <- Sp_COL_NAMES
colnames(GSE223748mF_Sp) <- Sp_COL_NAMES

Dist_M_Pra <- dist(t(GSE223748mM_Pr), method = "euclidean"); Dist_M_Spa <- dist(t(GSE223748mM_Sp), method = "euclidean")
Dist_F_Pra <- dist(t(GSE223748mF_Pr), method = "euclidean"); Dist_F_Spa <- dist(t(GSE223748mF_Sp), method = "euclidean")

#Neighbor-Joining (NJ): Unrooted Tree
Tree_M_Pra <- nj(Dist_M_Pra); Tree_M_Spa <- nj(Dist_M_Spa)
Tree_F_Pra <- nj(Dist_F_Pra); Tree_F_Spa <- nj(Dist_F_Spa)

#Epigenetic phylogeny
pdf("Output_Ox11_OXI/08_EpigeneticTree.pdf", width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 4, 4)) 
plot(Tree_M_Pra, main = "Epigenetic Tree: Male Pearson",    type = "unrooted", cex = 0.8, edge.width = 2, no.margin = FALSE)
plot(Tree_M_Spa, main = "Epigenetic Tree: Male Spearman",   type = "unrooted", cex = 0.8, edge.width = 2, no.margin = FALSE)
plot(Tree_F_Pra, main = "Epigenetic Tree: Female Pearson",  type = "unrooted", cex = 0.8, edge.width = 2, no.margin = FALSE)
plot(Tree_F_Spa, main = "Epigenetic Tree: Female Spearman", type = "unrooted", cex = 0.8, edge.width = 2, no.margin = FALSE)
dev.off()

#Plot 09: Tanglegram: Male vs Female, compare phylogeny structure by Entanglement
pdf("Output_Ox11_OXI/09_Tanglegram.pdf", width = 8, height = 4)
par(mfrow = c(1, 1))
Dend_M_Pra <- as.dendrogram(hclust(Dist_M_Pra)); Dend_F_Pra <- as.dendrogram(hclust(Dist_F_Pra))
labels(Dend_M_Pra) <- labels(Dend_F_Pra)[match(labels(Dend_M_Pra), labels(Dend_F_Pra))] 
tanglegram(Dend_M_Pra, Dend_F_Pra, main_left = "Male Pearson Tree", main_right = "Female Pearson Tree", lab.cex = 1, edge.lwd = 2, margin_inner = 7, columns_width = c(5, 2, 5))
dev.off()

#Plot 10: MDS Plot
MDS_M_Pra <- cmdscale(Dist_M_Pra); MDS_M_Pra <- data.frame(Species = rownames(MDS_M_Pra), Dim1 = MDS_M_Pra[, 1], Dim2 = MDS_M_Pra[, 2])
MDS_M_Spa <- cmdscale(Dist_M_Spa); MDS_M_Spa <- data.frame(Species = rownames(MDS_M_Spa), Dim1 = MDS_M_Spa[, 1], Dim2 = MDS_M_Spa[, 2])
MDS_F_Pra <- cmdscale(Dist_F_Pra); MDS_F_Pra <- data.frame(Species = rownames(MDS_F_Pra), Dim1 = MDS_F_Pra[, 1], Dim2 = MDS_F_Pra[, 2])
MDS_F_Spa <- cmdscale(Dist_F_Spa); MDS_F_Spa <- data.frame(Species = rownames(MDS_F_Spa), Dim1 = MDS_F_Spa[, 1], Dim2 = MDS_F_Spa[, 2])
Plot_M_Pra <- ggplot(MDS_M_Pra, aes(x = Dim1, y = Dim2, label = Species)) + geom_point(color = "firebrick", size = 4) + geom_text(vjust = -0.5, fontface = "italic") + labs(title = "Epigenetic Landscape: Male Pearson (MDS)   ", x = "Dimension 1", y = "Dimension 2") + theme_bw() + coord_cartesian(clip = "off")
Plot_M_Spa <- ggplot(MDS_M_Spa, aes(x = Dim1, y = Dim2, label = Species)) + geom_point(color = "firebrick", size = 4) + geom_text(vjust = -0.5, fontface = "italic") + labs(title = "Epigenetic Landscape: Male Spearman (MDS)  ", x = "Dimension 1", y = "Dimension 2") + theme_bw() + coord_cartesian(clip = "off")
Plot_F_Pra <- ggplot(MDS_F_Pra, aes(x = Dim1, y = Dim2, label = Species)) + geom_point(color = "firebrick", size = 4) + geom_text(vjust = -0.5, fontface = "italic") + labs(title = "Epigenetic Landscape: Female Pearson (MDS) ", x = "Dimension 1", y = "Dimension 2") + theme_bw() + coord_cartesian(clip = "off")
Plot_F_Spa <- ggplot(MDS_F_Spa, aes(x = Dim1, y = Dim2, label = Species)) + geom_point(color = "firebrick", size = 4) + geom_text(vjust = -0.5, fontface = "italic") + labs(title = "Epigenetic Landscape: Female Spearman (MDS)", x = "Dimension 1", y = "Dimension 2") + theme_bw() + coord_cartesian(clip = "off")
ggsave("Output_Ox11_OXI/10_MDS.pdf", plot=grid.arrange(Plot_M_Pra, Plot_M_Spa, Plot_F_Pra, Plot_F_Spa, ncol=2), width=12, height=8)



#Plot 11: PCA Plot
Run_PCA_Plot_a <- function(Mat_Data, Title_Text) {
  PCA_Res_a <- prcomp(t(Mat_Data))
  PCA_Dat_a <- data.frame(Species=rownames(PCA_Res_a$x), PC1=PCA_Res_a$x[,1], PC2=PCA_Res_a$x[,2])
  ggplot(PCA_Dat_a, aes(x=PC1, y=PC2, label=Species)) +
    geom_point(size=5, color="purple") + geom_text(vjust=-0.5, size=4, fontface="italic") +
    labs(title=Title_Text, x=paste0("PC1 (", round(summary(PCA_Res_a)$importance[2,1]*100,1), "%)"), y=paste0("PC2 (", round(summary(PCA_Res_a)$importance[2,2]*100,1), "%)")) + theme_bw() + coord_cartesian(clip = "off")}
P1_a <- Run_PCA_Plot_a(GSE223748mM_Pr, "PCA: Male Pearson"); P2_a <- Run_PCA_Plot_a(GSE223748mM_Sp, "PCA: Male Spearman")
P3_a <- Run_PCA_Plot_a(GSE223748mF_Pr, "PCA: Female Pearson"); P4_a <- Run_PCA_Plot_a(GSE223748mF_Sp, "PCA: Female Spearman")
ggsave("Output_Ox11_OXI/11_PCA.pdf", plot=grid.arrange(P1_a, P2_a, P3_a, P4_a, ncol=2), width=12, height=8)

