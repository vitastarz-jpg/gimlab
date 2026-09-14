################################################################################################
#Module 00. Clear memory and Load library
################################################################################################
rm(list=ls())
library(dplyr); library(tidyr); library(data.table); library(pheatmap); library(readxl)
library(ggvenn); library(gridExtra); library(viridis); library(ggpubr); library(grid)
library(ggplot2); library(reshape2); library(gridExtra); library(ape); library(dendextend); library(ggvenn)

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
#Module 11. eQTL
################################################################################################
source("Loading_TCGA_RNAseq_Met_Common_eQTL.R")

CpG_Table <- read.table("Output_Ox11_OXI/01_CpG_Table.txt",sep="\t",header=TRUE)
INTCpG01 <- CpG_Table$ID; INTCpG01 <- unique(INTCpG01[INTCpG01!=""])
INTExp01 <- CpG_Table$UCSC_RefGene_Name; INTExp01 <- unique(INTExp01[INTExp01!=""])
head(CpG_Table)

ACC__Exp01 <- ACC__Exp; ACC__Exp01$SYMBOL <- TCGA_GENE$SYMBOL; ACC__Exp01 <- dplyr::filter(ACC__Exp01, SYMBOL %in% INTExp01); rownames(ACC__Exp01) <- make.unique(ACC__Exp01$SYMBOL); ACC__Exp01$SYMBOL <- NULL; ACC__Met01 <- ACC__Met; ACC__Met01 <- na.omit(ACC__Met01[INTCpG01,])
BLCA_Exp01 <- BLCA_Exp; BLCA_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; BLCA_Exp01 <- dplyr::filter(BLCA_Exp01, SYMBOL %in% INTExp01); rownames(BLCA_Exp01) <- make.unique(BLCA_Exp01$SYMBOL); BLCA_Exp01$SYMBOL <- NULL; BLCA_Met01 <- BLCA_Met; BLCA_Met01 <- na.omit(BLCA_Met01[INTCpG01,])
BRCA_Exp01 <- BRCA_Exp; BRCA_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; BRCA_Exp01 <- dplyr::filter(BRCA_Exp01, SYMBOL %in% INTExp01); rownames(BRCA_Exp01) <- make.unique(BRCA_Exp01$SYMBOL); BRCA_Exp01$SYMBOL <- NULL; BRCA_Met01 <- BRCA_Met; BRCA_Met01 <- na.omit(BRCA_Met01[INTCpG01,])
CESC_Exp01 <- CESC_Exp; CESC_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; CESC_Exp01 <- dplyr::filter(CESC_Exp01, SYMBOL %in% INTExp01); rownames(CESC_Exp01) <- make.unique(CESC_Exp01$SYMBOL); CESC_Exp01$SYMBOL <- NULL; CESC_Met01 <- CESC_Met; CESC_Met01 <- na.omit(CESC_Met01[INTCpG01,])
CHOL_Exp01 <- CHOL_Exp; CHOL_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; CHOL_Exp01 <- dplyr::filter(CHOL_Exp01, SYMBOL %in% INTExp01); rownames(CHOL_Exp01) <- make.unique(CHOL_Exp01$SYMBOL); CHOL_Exp01$SYMBOL <- NULL; CHOL_Met01 <- CHOL_Met; CHOL_Met01 <- na.omit(CHOL_Met01[INTCpG01,])
COAD_Exp01 <- COAD_Exp; COAD_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; COAD_Exp01 <- dplyr::filter(COAD_Exp01, SYMBOL %in% INTExp01); rownames(COAD_Exp01) <- make.unique(COAD_Exp01$SYMBOL); COAD_Exp01$SYMBOL <- NULL; COAD_Met01 <- COAD_Met; COAD_Met01 <- na.omit(COAD_Met01[INTCpG01,])
DLBC_Exp01 <- DLBC_Exp; DLBC_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; DLBC_Exp01 <- dplyr::filter(DLBC_Exp01, SYMBOL %in% INTExp01); rownames(DLBC_Exp01) <- make.unique(DLBC_Exp01$SYMBOL); DLBC_Exp01$SYMBOL <- NULL; DLBC_Met01 <- DLBC_Met; DLBC_Met01 <- na.omit(DLBC_Met01[INTCpG01,])
ESCA_Exp01 <- ESCA_Exp; ESCA_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; ESCA_Exp01 <- dplyr::filter(ESCA_Exp01, SYMBOL %in% INTExp01); rownames(ESCA_Exp01) <- make.unique(ESCA_Exp01$SYMBOL); ESCA_Exp01$SYMBOL <- NULL; ESCA_Met01 <- ESCA_Met; ESCA_Met01 <- na.omit(ESCA_Met01[INTCpG01,])
GBM__Exp01 <- GBM__Exp; GBM__Exp01$SYMBOL <- TCGA_GENE$SYMBOL; GBM__Exp01 <- dplyr::filter(GBM__Exp01, SYMBOL %in% INTExp01); rownames(GBM__Exp01) <- make.unique(GBM__Exp01$SYMBOL); GBM__Exp01$SYMBOL <- NULL; GBM__Met01 <- GBM__Met; GBM__Met01 <- na.omit(GBM__Met01[INTCpG01,])
HNSC_Exp01 <- HNSC_Exp; HNSC_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; HNSC_Exp01 <- dplyr::filter(HNSC_Exp01, SYMBOL %in% INTExp01); rownames(HNSC_Exp01) <- make.unique(HNSC_Exp01$SYMBOL); HNSC_Exp01$SYMBOL <- NULL; HNSC_Met01 <- HNSC_Met; HNSC_Met01 <- na.omit(HNSC_Met01[INTCpG01,])
KICH_Exp01 <- KICH_Exp; KICH_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; KICH_Exp01 <- dplyr::filter(KICH_Exp01, SYMBOL %in% INTExp01); rownames(KICH_Exp01) <- make.unique(KICH_Exp01$SYMBOL); KICH_Exp01$SYMBOL <- NULL; KICH_Met01 <- KICH_Met; KICH_Met01 <- na.omit(KICH_Met01[INTCpG01,])
KIRC_Exp01 <- KIRC_Exp; KIRC_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; KIRC_Exp01 <- dplyr::filter(KIRC_Exp01, SYMBOL %in% INTExp01); rownames(KIRC_Exp01) <- make.unique(KIRC_Exp01$SYMBOL); KIRC_Exp01$SYMBOL <- NULL; KIRC_Met01 <- KIRC_Met; KIRC_Met01 <- na.omit(KIRC_Met01[INTCpG01,])
KIRP_Exp01 <- KIRP_Exp; KIRP_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; KIRP_Exp01 <- dplyr::filter(KIRP_Exp01, SYMBOL %in% INTExp01); rownames(KIRP_Exp01) <- make.unique(KIRP_Exp01$SYMBOL); KIRP_Exp01$SYMBOL <- NULL; KIRP_Met01 <- KIRP_Met; KIRP_Met01 <- na.omit(KIRP_Met01[INTCpG01,])
LAML_Exp01 <- LAML_Exp; LAML_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; LAML_Exp01 <- dplyr::filter(LAML_Exp01, SYMBOL %in% INTExp01); rownames(LAML_Exp01) <- make.unique(LAML_Exp01$SYMBOL); LAML_Exp01$SYMBOL <- NULL; LAML_Met01 <- LAML_Met; LAML_Met01 <- na.omit(LAML_Met01[INTCpG01,])
LGG__Exp01 <- LGG__Exp; LGG__Exp01$SYMBOL <- TCGA_GENE$SYMBOL; LGG__Exp01 <- dplyr::filter(LGG__Exp01, SYMBOL %in% INTExp01); rownames(LGG__Exp01) <- make.unique(LGG__Exp01$SYMBOL); LGG__Exp01$SYMBOL <- NULL; LGG__Met01 <- LGG__Met; LGG__Met01 <- na.omit(LGG__Met01[INTCpG01,])
LIHC_Exp01 <- LIHC_Exp; LIHC_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; LIHC_Exp01 <- dplyr::filter(LIHC_Exp01, SYMBOL %in% INTExp01); rownames(LIHC_Exp01) <- make.unique(LIHC_Exp01$SYMBOL); LIHC_Exp01$SYMBOL <- NULL; LIHC_Met01 <- LIHC_Met; LIHC_Met01 <- na.omit(LIHC_Met01[INTCpG01,])
LUAD_Exp01 <- LUAD_Exp; LUAD_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; LUAD_Exp01 <- dplyr::filter(LUAD_Exp01, SYMBOL %in% INTExp01); rownames(LUAD_Exp01) <- make.unique(LUAD_Exp01$SYMBOL); LUAD_Exp01$SYMBOL <- NULL; LUAD_Met01 <- LUAD_Met; LUAD_Met01 <- na.omit(LUAD_Met01[INTCpG01,])
LUSC_Exp01 <- LUSC_Exp; LUSC_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; LUSC_Exp01 <- dplyr::filter(LUSC_Exp01, SYMBOL %in% INTExp01); rownames(LUSC_Exp01) <- make.unique(LUSC_Exp01$SYMBOL); LUSC_Exp01$SYMBOL <- NULL; LUSC_Met01 <- LUSC_Met; LUSC_Met01 <- na.omit(LUSC_Met01[INTCpG01,])
MESO_Exp01 <- MESO_Exp; MESO_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; MESO_Exp01 <- dplyr::filter(MESO_Exp01, SYMBOL %in% INTExp01); rownames(MESO_Exp01) <- make.unique(MESO_Exp01$SYMBOL); MESO_Exp01$SYMBOL <- NULL; MESO_Met01 <- MESO_Met; MESO_Met01 <- na.omit(MESO_Met01[INTCpG01,])
PAAD_Exp01 <- PAAD_Exp; PAAD_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; PAAD_Exp01 <- dplyr::filter(PAAD_Exp01, SYMBOL %in% INTExp01); rownames(PAAD_Exp01) <- make.unique(PAAD_Exp01$SYMBOL); PAAD_Exp01$SYMBOL <- NULL; PAAD_Met01 <- PAAD_Met; PAAD_Met01 <- na.omit(PAAD_Met01[INTCpG01,])
PCPG_Exp01 <- PCPG_Exp; PCPG_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; PCPG_Exp01 <- dplyr::filter(PCPG_Exp01, SYMBOL %in% INTExp01); rownames(PCPG_Exp01) <- make.unique(PCPG_Exp01$SYMBOL); PCPG_Exp01$SYMBOL <- NULL; PCPG_Met01 <- PCPG_Met; PCPG_Met01 <- na.omit(PCPG_Met01[INTCpG01,])
PRAD_Exp01 <- PRAD_Exp; PRAD_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; PRAD_Exp01 <- dplyr::filter(PRAD_Exp01, SYMBOL %in% INTExp01); rownames(PRAD_Exp01) <- make.unique(PRAD_Exp01$SYMBOL); PRAD_Exp01$SYMBOL <- NULL; PRAD_Met01 <- PRAD_Met; PRAD_Met01 <- na.omit(PRAD_Met01[INTCpG01,])
READ_Exp01 <- READ_Exp; READ_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; READ_Exp01 <- dplyr::filter(READ_Exp01, SYMBOL %in% INTExp01); rownames(READ_Exp01) <- make.unique(READ_Exp01$SYMBOL); READ_Exp01$SYMBOL <- NULL; READ_Met01 <- READ_Met; READ_Met01 <- na.omit(READ_Met01[INTCpG01,])
SKCM_Exp01 <- SKCM_Exp; SKCM_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; SKCM_Exp01 <- dplyr::filter(SKCM_Exp01, SYMBOL %in% INTExp01); rownames(SKCM_Exp01) <- make.unique(SKCM_Exp01$SYMBOL); SKCM_Exp01$SYMBOL <- NULL; SKCM_Met01 <- SKCM_Met; SKCM_Met01 <- na.omit(SKCM_Met01[INTCpG01,])
STAD_Exp01 <- STAD_Exp; STAD_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; STAD_Exp01 <- dplyr::filter(STAD_Exp01, SYMBOL %in% INTExp01); rownames(STAD_Exp01) <- make.unique(STAD_Exp01$SYMBOL); STAD_Exp01$SYMBOL <- NULL; STAD_Met01 <- STAD_Met; STAD_Met01 <- na.omit(STAD_Met01[INTCpG01,])
TGCT_Exp01 <- TGCT_Exp; TGCT_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; TGCT_Exp01 <- dplyr::filter(TGCT_Exp01, SYMBOL %in% INTExp01); rownames(TGCT_Exp01) <- make.unique(TGCT_Exp01$SYMBOL); TGCT_Exp01$SYMBOL <- NULL; TGCT_Met01 <- TGCT_Met; TGCT_Met01 <- na.omit(TGCT_Met01[INTCpG01,])
THCA_Exp01 <- THCA_Exp; THCA_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; THCA_Exp01 <- dplyr::filter(THCA_Exp01, SYMBOL %in% INTExp01); rownames(THCA_Exp01) <- make.unique(THCA_Exp01$SYMBOL); THCA_Exp01$SYMBOL <- NULL; THCA_Met01 <- THCA_Met; THCA_Met01 <- na.omit(THCA_Met01[INTCpG01,])
THYM_Exp01 <- THYM_Exp; THYM_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; THYM_Exp01 <- dplyr::filter(THYM_Exp01, SYMBOL %in% INTExp01); rownames(THYM_Exp01) <- make.unique(THYM_Exp01$SYMBOL); THYM_Exp01$SYMBOL <- NULL; THYM_Met01 <- THYM_Met; THYM_Met01 <- na.omit(THYM_Met01[INTCpG01,])
UCEC_Exp01 <- UCEC_Exp; UCEC_Exp01$SYMBOL <- TCGA_GENE$SYMBOL; UCEC_Exp01 <- dplyr::filter(UCEC_Exp01, SYMBOL %in% INTExp01); rownames(UCEC_Exp01) <- make.unique(UCEC_Exp01$SYMBOL); UCEC_Exp01$SYMBOL <- NULL; UCEC_Met01 <- UCEC_Met; UCEC_Met01 <- na.omit(UCEC_Met01[INTCpG01,])
UCS__Exp01 <- UCS__Exp; UCS__Exp01$SYMBOL <- TCGA_GENE$SYMBOL; UCS__Exp01 <- dplyr::filter(UCS__Exp01, SYMBOL %in% INTExp01); rownames(UCS__Exp01) <- make.unique(UCS__Exp01$SYMBOL); UCS__Exp01$SYMBOL <- NULL; UCS__Met01 <- UCS__Met; UCS__Met01 <- na.omit(UCS__Met01[INTCpG01,])
UVM__Exp01 <- UVM__Exp; UVM__Exp01$SYMBOL <- TCGA_GENE$SYMBOL; UVM__Exp01 <- dplyr::filter(UVM__Exp01, SYMBOL %in% INTExp01); rownames(UVM__Exp01) <- make.unique(UVM__Exp01$SYMBOL); UVM__Exp01$SYMBOL <- NULL; UVM__Met01 <- UVM__Met; UVM__Met01 <- na.omit(UVM__Met01[INTCpG01,])

prefixes <- unlist(strsplit("ACC__ BLCA_ BRCA_ CESC_ CHOL_ COAD_ DLBC_ ESCA_ GBM__ HNSC_ KICH_ KIRC_ KIRP_ LAML_ LGG__ LIHC_ LUAD_ LUSC_ MESO_ PAAD_ PCPG_ PRAD_ READ_ SKCM_ STAD_ TGCT_ THCA_ THYM_ UCEC_ UCS__ UVM__", " ")); res_all <- NULL
for(ca in prefixes){ exp_df <- get(paste0(ca, "Exp01")); met_df <- get(paste0(ca, "Met01")); val <- CpG_Table[CpG_Table$ID %in% rownames(met_df) & CpG_Table$UCSC_RefGene_Name %in% rownames(exp_df), ]; if(nrow(val)>0){ res <- data.frame(Cancer=gsub("_", "", ca), ID=val$ID, Gene=val$UCSC_RefGene_Name, Pearson_R=0, Pearson_P=0, Spearman_R=0, Spearman_P=0); for(i in 1:nrow(res)){ met <- as.numeric(met_df[res$ID[i], ]); exp <- as.numeric(exp_df[res$Gene[i], ]); p_cor <- cor.test(met, exp, method="pearson"); s_cor <- cor.test(met, exp, method="spearman", exact=FALSE); res[i, 4:7] <- c(p_cor$estimate, p_cor$p.value, s_cor$estimate, s_cor$p.value) }; res_all <- rbind(res_all, res) } }
res_all$HyPi <- (res_all$Pearson_R * -log10(res_all$Pearson_P)) + (res_all$Spearman_R * -log10(res_all$Spearman_P)); res_all[,4:8] <- lapply(res_all[,4:8], function(x) round(as.numeric(x), 3))
write.table(res_all, "Output_Ox11_OXI/21_CpG_Gene_Correlation_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)

ACC__Exp01[1:5,1:5]; dim(ACC__Exp01); table(colnames(ACC__Exp01)==ACC__c$ID); ACC__Met01[1:5,1:5]; dim(ACC__Met01); table(colnames(ACC__Met01)==ACC__c$ID); head(ACC__c)
BLCA_Exp01[1:5,1:5]; dim(BLCA_Exp01); table(colnames(BLCA_Exp01)==BLCA_c$ID); BLCA_Met01[1:5,1:5]; dim(BLCA_Met01); table(colnames(BLCA_Met01)==BLCA_c$ID); head(BLCA_c)
BRCA_Exp01[1:5,1:5]; dim(BRCA_Exp01); table(colnames(BRCA_Exp01)==BRCA_c$ID); BRCA_Met01[1:5,1:5]; dim(BRCA_Met01); table(colnames(BRCA_Met01)==BRCA_c$ID); head(BRCA_c)
CESC_Exp01[1:5,1:5]; dim(CESC_Exp01); table(colnames(CESC_Exp01)==CESC_c$ID); CESC_Met01[1:5,1:5]; dim(CESC_Met01); table(colnames(CESC_Met01)==CESC_c$ID); head(CESC_c)
CHOL_Exp01[1:5,1:5]; dim(CHOL_Exp01); table(colnames(CHOL_Exp01)==CHOL_c$ID); CHOL_Met01[1:5,1:5]; dim(CHOL_Met01); table(colnames(CHOL_Met01)==CHOL_c$ID); head(CHOL_c)
COAD_Exp01[1:5,1:5]; dim(COAD_Exp01); table(colnames(COAD_Exp01)==COAD_c$ID); COAD_Met01[1:5,1:5]; dim(COAD_Met01); table(colnames(COAD_Met01)==COAD_c$ID); head(COAD_c)
DLBC_Exp01[1:5,1:5]; dim(DLBC_Exp01); table(colnames(DLBC_Exp01)==DLBC_c$ID); DLBC_Met01[1:5,1:5]; dim(DLBC_Met01); table(colnames(DLBC_Met01)==DLBC_c$ID); head(DLBC_c)
ESCA_Exp01[1:5,1:5]; dim(ESCA_Exp01); table(colnames(ESCA_Exp01)==ESCA_c$ID); ESCA_Met01[1:5,1:5]; dim(ESCA_Met01); table(colnames(ESCA_Met01)==ESCA_c$ID); head(ESCA_c)
GBM__Exp01[1:5,1:5]; dim(GBM__Exp01); table(colnames(GBM__Exp01)==GBM__c$ID); GBM__Met01[1:5,1:5]; dim(GBM__Met01); table(colnames(GBM__Met01)==GBM__c$ID); head(GBM__c)
HNSC_Exp01[1:5,1:5]; dim(HNSC_Exp01); table(colnames(HNSC_Exp01)==HNSC_c$ID); HNSC_Met01[1:5,1:5]; dim(HNSC_Met01); table(colnames(HNSC_Met01)==HNSC_c$ID); head(HNSC_c)
KICH_Exp01[1:5,1:5]; dim(KICH_Exp01); table(colnames(KICH_Exp01)==KICH_c$ID); KICH_Met01[1:5,1:5]; dim(KICH_Met01); table(colnames(KICH_Met01)==KICH_c$ID); head(KICH_c)
KIRC_Exp01[1:5,1:5]; dim(KIRC_Exp01); table(colnames(KIRC_Exp01)==KIRC_c$ID); KIRC_Met01[1:5,1:5]; dim(KIRC_Met01); table(colnames(KIRC_Met01)==KIRC_c$ID); head(KIRC_c)
KIRP_Exp01[1:5,1:5]; dim(KIRP_Exp01); table(colnames(KIRP_Exp01)==KIRP_c$ID); KIRP_Met01[1:5,1:5]; dim(KIRP_Met01); table(colnames(KIRP_Met01)==KIRP_c$ID); head(KIRP_c)
LAML_Exp01[1:5,1:5]; dim(LAML_Exp01); table(colnames(LAML_Exp01)==LAML_c$ID); LAML_Met01[1:5,1:5]; dim(LAML_Met01); table(colnames(LAML_Met01)==LAML_c$ID); head(LAML_c)
LGG__Exp01[1:5,1:5]; dim(LGG__Exp01); table(colnames(LGG__Exp01)==LGG__c$ID); LGG__Met01[1:5,1:5]; dim(LGG__Met01); table(colnames(LGG__Met01)==LGG__c$ID); head(LGG__c)
LIHC_Exp01[1:5,1:5]; dim(LIHC_Exp01); table(colnames(LIHC_Exp01)==LIHC_c$ID); LIHC_Met01[1:5,1:5]; dim(LIHC_Met01); table(colnames(LIHC_Met01)==LIHC_c$ID); head(LIHC_c)
LUAD_Exp01[1:5,1:5]; dim(LUAD_Exp01); table(colnames(LUAD_Exp01)==LUAD_c$ID); LUAD_Met01[1:5,1:5]; dim(LUAD_Met01); table(colnames(LUAD_Met01)==LUAD_c$ID); head(LUAD_c)
LUSC_Exp01[1:5,1:5]; dim(LUSC_Exp01); table(colnames(LUSC_Exp01)==LUSC_c$ID); LUSC_Met01[1:5,1:5]; dim(LUSC_Met01); table(colnames(LUSC_Met01)==LUSC_c$ID); head(LUSC_c)
MESO_Exp01[1:5,1:5]; dim(MESO_Exp01); table(colnames(MESO_Exp01)==MESO_c$ID); MESO_Met01[1:5,1:5]; dim(MESO_Met01); table(colnames(MESO_Met01)==MESO_c$ID); head(MESO_c)
PAAD_Exp01[1:5,1:5]; dim(PAAD_Exp01); table(colnames(PAAD_Exp01)==PAAD_c$ID); PAAD_Met01[1:5,1:5]; dim(PAAD_Met01); table(colnames(PAAD_Met01)==PAAD_c$ID); head(PAAD_c)
PCPG_Exp01[1:5,1:5]; dim(PCPG_Exp01); table(colnames(PCPG_Exp01)==PCPG_c$ID); PCPG_Met01[1:5,1:5]; dim(PCPG_Met01); table(colnames(PCPG_Met01)==PCPG_c$ID); head(PCPG_c)
PRAD_Exp01[1:5,1:5]; dim(PRAD_Exp01); table(colnames(PRAD_Exp01)==PRAD_c$ID); PRAD_Met01[1:5,1:5]; dim(PRAD_Met01); table(colnames(PRAD_Met01)==PRAD_c$ID); head(PRAD_c)
READ_Exp01[1:5,1:5]; dim(READ_Exp01); table(colnames(READ_Exp01)==READ_c$ID); READ_Met01[1:5,1:5]; dim(READ_Met01); table(colnames(READ_Met01)==READ_c$ID); head(READ_c)
SKCM_Exp01[1:5,1:5]; dim(SKCM_Exp01); table(colnames(SKCM_Exp01)==SKCM_c$ID); SKCM_Met01[1:5,1:5]; dim(SKCM_Met01); table(colnames(SKCM_Met01)==SKCM_c$ID); head(SKCM_c)
STAD_Exp01[1:5,1:5]; dim(STAD_Exp01); table(colnames(STAD_Exp01)==STAD_c$ID); STAD_Met01[1:5,1:5]; dim(STAD_Met01); table(colnames(STAD_Met01)==STAD_c$ID); head(STAD_c)
TGCT_Exp01[1:5,1:5]; dim(TGCT_Exp01); table(colnames(TGCT_Exp01)==TGCT_c$ID); TGCT_Met01[1:5,1:5]; dim(TGCT_Met01); table(colnames(TGCT_Met01)==TGCT_c$ID); head(TGCT_c)
THCA_Exp01[1:5,1:5]; dim(THCA_Exp01); table(colnames(THCA_Exp01)==THCA_c$ID); THCA_Met01[1:5,1:5]; dim(THCA_Met01); table(colnames(THCA_Met01)==THCA_c$ID); head(THCA_c)
THYM_Exp01[1:5,1:5]; dim(THYM_Exp01); table(colnames(THYM_Exp01)==THYM_c$ID); THYM_Met01[1:5,1:5]; dim(THYM_Met01); table(colnames(THYM_Met01)==THYM_c$ID); head(THYM_c)
UCEC_Exp01[1:5,1:5]; dim(UCEC_Exp01); table(colnames(UCEC_Exp01)==UCEC_c$ID); UCEC_Met01[1:5,1:5]; dim(UCEC_Met01); table(colnames(UCEC_Met01)==UCEC_c$ID); head(UCEC_c)
UCS__Exp01[1:5,1:5]; dim(UCS__Exp01); table(colnames(UCS__Exp01)==UCS__c$ID); UCS__Met01[1:5,1:5]; dim(UCS__Met01); table(colnames(UCS__Met01)==UCS__c$ID); head(UCS__c)
UVM__Exp01[1:5,1:5]; dim(UVM__Exp01); table(colnames(UVM__Exp01)==UVM__c$ID); UVM__Met01[1:5,1:5]; dim(UVM__Met01); table(colnames(UVM__Met01)==UVM__c$ID); head(UVM__c)

prefixes <- unlist(strsplit("ACC__ BLCA_ BRCA_ CESC_ CHOL_ COAD_ DLBC_ ESCA_ GBM__ HNSC_ KICH_ KIRC_ KIRP_ LAML_ LGG__ LIHC_ LUAD_ LUSC_ MESO_ PAAD_ PCPG_ PRAD_ READ_ SKCM_ STAD_ TGCT_ THCA_ THYM_ UCEC_ UCS__ UVM__", " ")); res_all <- NULL
for(ca in prefixes){ c_df <- get(paste0(ca, "c")); for(dt in c("Exp", "Met")){ dat <- get(paste0(ca, dt, "01")); g1 <- which(as.numeric(c_df$Surv)==1); g2 <- which(as.numeric(c_df$Surv)==2); res <- data.frame(ID=rownames(dat), FC=NA, P_value=NA, FDR=NA, Cancer=gsub("_", "", ca), Type=dt, stringsAsFactors=FALSE); for(i in 1:nrow(dat)){ v1 <- as.numeric(dat[i, g1]); v2 <- as.numeric(dat[i, g2]); if(sum(!is.na(v1))>1 & sum(!is.na(v2))>1){ if(var(v1, na.rm=TRUE)>0 | var(v2, na.rm=TRUE)>0){ res$FC[i] <- (mean(v2, na.rm=TRUE)+1e-6)/(mean(v1, na.rm=TRUE)+1e-6); res$P_value[i] <- t.test(v2, v1)$p.value } } }; res$FDR <- p.adjust(res$P_value, method="fdr"); res_all <- rbind(res_all, res) } }
write.table(res_all, "Output_Ox11_OXI/22_DEG_Survival_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)

prefixes <- unlist(strsplit("ACC__ BLCA_ BRCA_ CESC_ CHOL_ COAD_ DLBC_ ESCA_ GBM__ HNSC_ KICH_ KIRC_ KIRP_ LAML_ LGG__ LIHC_ LUAD_ LUSC_ MESO_ PAAD_ PCPG_ PRAD_ READ_ SKCM_ STAD_ TGCT_ THCA_ THYM_ UCEC_ UCS__ UVM__", " ")); res_all <- NULL
for(ca in prefixes){ c_df <- get(paste0(ca, "c")); for(dt in c("Exp", "Met")){ dat <- get(paste0(ca, dt, "01")); g1 <- which(as.numeric(c_df$Surv)==1); g2 <- which(as.numeric(c_df$Surv)==2); res <- data.frame(ID=rownames(dat), log2FC=NA, P_value=NA, FDR=NA, Cancer=gsub("_", "", ca), Type=dt, stringsAsFactors=FALSE); for(i in 1:nrow(dat)){ v1 <- as.numeric(dat[i, g1]); v2 <- as.numeric(dat[i, g2]); if(sum(!is.na(v1))>1 & sum(!is.na(v2))>1){ if(var(v1, na.rm=TRUE)>0 | var(v2, na.rm=TRUE)>0){ res$log2FC[i] <- log2((mean(v2, na.rm=TRUE)+1e-6)/(mean(v1, na.rm=TRUE)+1e-6)); res$P_value[i] <- t.test(v2, v1)$p.value } } }; res$FDR <- p.adjust(res$P_value, method="fdr"); res_all <- rbind(res_all, res) } }
write.table(res_all, "Output_Ox11_OXI/22_DEG_Survival_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)

head(res_all); dim(res_all); table(res_all[,c(6,5)])

pdf("Output_Ox11_OXI/22_DEG_Heatmaps.pdf", width=8, height=10); cancers <- unique(res_all$Cancer)
for(dt in c("Exp", "Met")){ sub <- res_all[res_all$Type==dt, ]; ids <- unique(sub$ID); m_fc <- matrix(NA, length(cancers), length(ids), dimnames=list(cancers, ids)); m_fdr <- m_fc; for(i in 1:nrow(sub)){ m_fc[sub$Cancer[i], sub$ID[i]] <- sub$log2FC[i]; m_fdr[sub$Cancer[i], sub$ID[i]] <- sub$FDR[i] }; m_sig <- matrix(ifelse(!is.na(m_fdr) & m_fdr < 0.05, "*", ""), nrow(m_fdr), dimnames=dimnames(m_fdr)); mx <- max(abs(m_fc), na.rm=TRUE); mx <- ifelse(is.infinite(mx) | mx==0, 1, mx); pheatmap(m_fc, color=colorRampPalette(c("blue", "white", "red"))(50), breaks=seq(-mx, mx, length.out=51), display_numbers=m_sig, cluster_rows=FALSE, cluster_cols=FALSE, na_col="gray90", main=paste(dt, "log2FC Heatmap"), fontsize_number=15, number_color="black") }
dev.off()

r_m <- res_all[res_all$Type=="Met" & !is.na(res_all$FDR) & res_all$FDR<0.05, c("Cancer", "ID", "log2FC", "FDR")]; r_e <- res_all[res_all$Type=="Exp" & !is.na(res_all$FDR) & res_all$FDR<0.05, c("Cancer", "ID", "log2FC", "FDR")]; cpg <- CpG_Table[, c("ID", "UCSC_RefGene_Name")]; colnames(r_m)[2] <- "CpG"; colnames(cpg) <- c("CpG", "Gene"); r_m <- merge(r_m, cpg, by="CpG"); colnames(r_e)[2] <- "Gene"; res_pair <- merge(r_m, r_e, by=c("Cancer", "Gene"), suffixes=c("_Met", "_Exp")); write.table(res_pair, "Output_Ox11_OXI/22_Significant_Pairs.txt", sep="\t", quote=FALSE, row.names=FALSE)
head(res_pair); dim(res_pair)

ACC__c_SG <- ACC__c; ACC__c_SG$Stage <- gsub("Stage ","",ACC__c_SG$Stage); s <- trimws(as.character(ACC__c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; ACC__c_SG$Stage <- g; ACC__c_SG <- ACC__c_SG[!is.na(ACC__c_SG$Stage), ]
BLCA_c_SG <- BLCA_c; BLCA_c_SG$Stage <- gsub("Stage ","",BLCA_c_SG$Stage); s <- trimws(as.character(BLCA_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; BLCA_c_SG$Stage <- g; BLCA_c_SG <- BLCA_c_SG[!is.na(BLCA_c_SG$Stage), ]
BRCA_c_SG <- BRCA_c; BRCA_c_SG$Stage <- gsub("Stage ","",BRCA_c_SG$Stage); s <- trimws(as.character(BRCA_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; BRCA_c_SG$Stage <- g; BRCA_c_SG <- BRCA_c_SG[!is.na(BRCA_c_SG$Stage), ]
CESC_c_SG <- CESC_c; CESC_c_SG$Stage <- gsub("Stage ","",CESC_c_SG$Stage); s <- trimws(as.character(CESC_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; CESC_c_SG$Stage <- g; CESC_c_SG <- CESC_c_SG[!is.na(CESC_c_SG$Stage), ]
CHOL_c_SG <- CHOL_c; CHOL_c_SG$Stage <- gsub("Stage ","",CHOL_c_SG$Stage); s <- trimws(as.character(CHOL_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; CHOL_c_SG$Stage <- g; CHOL_c_SG <- CHOL_c_SG[!is.na(CHOL_c_SG$Stage), ]
COAD_c_SG <- COAD_c; COAD_c_SG$Stage <- gsub("Stage ","",COAD_c_SG$Stage); s <- trimws(as.character(COAD_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; COAD_c_SG$Stage <- g; COAD_c_SG <- COAD_c_SG[!is.na(COAD_c_SG$Stage), ]
DLBC_c_SG <- DLBC_c; DLBC_c_SG$Stage <- gsub("Stage ","",DLBC_c_SG$Stage); s <- trimws(as.character(DLBC_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; DLBC_c_SG$Stage <- g; DLBC_c_SG <- DLBC_c_SG[!is.na(DLBC_c_SG$Stage), ]
ESCA_c_SG <- ESCA_c; ESCA_c_SG$Stage <- gsub("Stage ","",ESCA_c_SG$Stage); s <- trimws(as.character(ESCA_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; ESCA_c_SG$Stage <- g; ESCA_c_SG <- ESCA_c_SG[!is.na(ESCA_c_SG$Stage), ]
HNSC_c_SG <- HNSC_c; HNSC_c_SG$Stage <- gsub("Stage ","",HNSC_c_SG$Stage); s <- trimws(as.character(HNSC_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; HNSC_c_SG$Stage <- g; HNSC_c_SG <- HNSC_c_SG[!is.na(HNSC_c_SG$Stage), ]
KICH_c_SG <- KICH_c; KICH_c_SG$Stage <- gsub("Stage ","",KICH_c_SG$Stage); s <- trimws(as.character(KICH_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; KICH_c_SG$Stage <- g; KICH_c_SG <- KICH_c_SG[!is.na(KICH_c_SG$Stage), ]
KIRC_c_SG <- KIRC_c; KIRC_c_SG$Stage <- gsub("Stage ","",KIRC_c_SG$Stage); s <- trimws(as.character(KIRC_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; KIRC_c_SG$Stage <- g; KIRC_c_SG <- KIRC_c_SG[!is.na(KIRC_c_SG$Stage), ]
KIRP_c_SG <- KIRP_c; KIRP_c_SG$Stage <- gsub("Stage ","",KIRP_c_SG$Stage); s <- trimws(as.character(KIRP_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; KIRP_c_SG$Stage <- g; KIRP_c_SG <- KIRP_c_SG[!is.na(KIRP_c_SG$Stage), ]
LIHC_c_SG <- LIHC_c; LIHC_c_SG$Stage <- gsub("Stage ","",LIHC_c_SG$Stage); s <- trimws(as.character(LIHC_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; LIHC_c_SG$Stage <- g; LIHC_c_SG <- LIHC_c_SG[!is.na(LIHC_c_SG$Stage), ]
LUAD_c_SG <- LUAD_c; LUAD_c_SG$Stage <- gsub("Stage ","",LUAD_c_SG$Stage); s <- trimws(as.character(LUAD_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; LUAD_c_SG$Stage <- g; LUAD_c_SG <- LUAD_c_SG[!is.na(LUAD_c_SG$Stage), ]
LUSC_c_SG <- LUSC_c; LUSC_c_SG$Stage <- gsub("Stage ","",LUSC_c_SG$Stage); s <- trimws(as.character(LUSC_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; LUSC_c_SG$Stage <- g; LUSC_c_SG <- LUSC_c_SG[!is.na(LUSC_c_SG$Stage), ]
MESO_c_SG <- MESO_c; MESO_c_SG$Stage <- gsub("Stage ","",MESO_c_SG$Stage); s <- trimws(as.character(MESO_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; MESO_c_SG$Stage <- g; MESO_c_SG <- MESO_c_SG[!is.na(MESO_c_SG$Stage), ]
PAAD_c_SG <- PAAD_c; PAAD_c_SG$Stage <- gsub("Stage ","",PAAD_c_SG$Stage); s <- trimws(as.character(PAAD_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; PAAD_c_SG$Stage <- g; PAAD_c_SG <- PAAD_c_SG[!is.na(PAAD_c_SG$Stage), ]
READ_c_SG <- READ_c; READ_c_SG$Stage <- gsub("Stage ","",READ_c_SG$Stage); s <- trimws(as.character(READ_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; READ_c_SG$Stage <- g; READ_c_SG <- READ_c_SG[!is.na(READ_c_SG$Stage), ]
SKCM_c_SG <- SKCM_c; SKCM_c_SG$Stage <- gsub("Stage ","",SKCM_c_SG$Stage); s <- trimws(as.character(SKCM_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; SKCM_c_SG$Stage <- g; SKCM_c_SG <- SKCM_c_SG[!is.na(SKCM_c_SG$Stage), ]
STAD_c_SG <- STAD_c; STAD_c_SG$Stage <- gsub("Stage ","",STAD_c_SG$Stage); s <- trimws(as.character(STAD_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; STAD_c_SG$Stage <- g; STAD_c_SG <- STAD_c_SG[!is.na(STAD_c_SG$Stage), ]
TGCT_c_SG <- TGCT_c; TGCT_c_SG$Stage <- gsub("Stage ","",TGCT_c_SG$Stage); s <- trimws(as.character(TGCT_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; TGCT_c_SG$Stage <- g; TGCT_c_SG <- TGCT_c_SG[!is.na(TGCT_c_SG$Stage), ]
THCA_c_SG <- THCA_c; THCA_c_SG$Stage <- gsub("Stage ","",THCA_c_SG$Stage); s <- trimws(as.character(THCA_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; THCA_c_SG$Stage <- g; THCA_c_SG <- THCA_c_SG[!is.na(THCA_c_SG$Stage), ]
THYM_c_SG <- THYM_c; THYM_c_SG$Stage <- gsub("Stage ","",THYM_c_SG$Stage); s <- trimws(as.character(THYM_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; THYM_c_SG$Stage <- g; THYM_c_SG <- THYM_c_SG[!is.na(THYM_c_SG$Stage), ]
UCEC_c_SG <- UCEC_c; UCEC_c_SG$Stage <- gsub("Stage ","",UCEC_c_SG$Stage); s <- trimws(as.character(UCEC_c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; UCEC_c_SG$Stage <- g; UCEC_c_SG <- UCEC_c_SG[!is.na(UCEC_c_SG$Stage), ]
UCS__c_SG <- UCS__c; UCS__c_SG$Stage <- gsub("Stage ","",UCS__c_SG$Stage); s <- trimws(as.character(UCS__c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; UCS__c_SG$Stage <- g; UCS__c_SG <- UCS__c_SG[!is.na(UCS__c_SG$Stage), ]
UVM__c_SG <- UVM__c; UVM__c_SG$Stage <- gsub("Stage ","",UVM__c_SG$Stage); s <- trimws(as.character(UVM__c_SG$Stage)); g <- rep(NA, length(s)); g[grepl("^I$|^I[A-Ca-c]|^II$|^II[A-Ca-c]", s)] <- "L"; g[grepl("^III|^IV", s)] <- "M"; UVM__c_SG$Stage <- g; UVM__c_SG <- UVM__c_SG[!is.na(UVM__c_SG$Stage), ]

ACC__Exp01[1:5,1:5]; dim(ACC__Exp01); table(colnames(ACC__Exp01)==ACC__c$ID); ACC__Met01[1:5,1:5]; dim(ACC__Met01); table(colnames(ACC__Met01)==ACC__c$ID); head(ACC__c_SG)
BLCA_Exp01[1:5,1:5]; dim(BLCA_Exp01); table(colnames(BLCA_Exp01)==BLCA_c$ID); BLCA_Met01[1:5,1:5]; dim(BLCA_Met01); table(colnames(BLCA_Met01)==BLCA_c$ID); head(BLCA_c_SG)
BRCA_Exp01[1:5,1:5]; dim(BRCA_Exp01); table(colnames(BRCA_Exp01)==BRCA_c$ID); BRCA_Met01[1:5,1:5]; dim(BRCA_Met01); table(colnames(BRCA_Met01)==BRCA_c$ID); head(BRCA_c_SG)
CESC_Exp01[1:5,1:5]; dim(CESC_Exp01); table(colnames(CESC_Exp01)==CESC_c$ID); CESC_Met01[1:5,1:5]; dim(CESC_Met01); table(colnames(CESC_Met01)==CESC_c$ID); head(CESC_c_SG)
CHOL_Exp01[1:5,1:5]; dim(CHOL_Exp01); table(colnames(CHOL_Exp01)==CHOL_c$ID); CHOL_Met01[1:5,1:5]; dim(CHOL_Met01); table(colnames(CHOL_Met01)==CHOL_c$ID); head(CHOL_c_SG)
COAD_Exp01[1:5,1:5]; dim(COAD_Exp01); table(colnames(COAD_Exp01)==COAD_c$ID); COAD_Met01[1:5,1:5]; dim(COAD_Met01); table(colnames(COAD_Met01)==COAD_c$ID); head(COAD_c_SG)
DLBC_Exp01[1:5,1:5]; dim(DLBC_Exp01); table(colnames(DLBC_Exp01)==DLBC_c$ID); DLBC_Met01[1:5,1:5]; dim(DLBC_Met01); table(colnames(DLBC_Met01)==DLBC_c$ID); head(DLBC_c_SG)
ESCA_Exp01[1:5,1:5]; dim(ESCA_Exp01); table(colnames(ESCA_Exp01)==ESCA_c$ID); ESCA_Met01[1:5,1:5]; dim(ESCA_Met01); table(colnames(ESCA_Met01)==ESCA_c$ID); head(ESCA_c_SG)
HNSC_Exp01[1:5,1:5]; dim(HNSC_Exp01); table(colnames(HNSC_Exp01)==HNSC_c$ID); HNSC_Met01[1:5,1:5]; dim(HNSC_Met01); table(colnames(HNSC_Met01)==HNSC_c$ID); head(HNSC_c_SG)
KICH_Exp01[1:5,1:5]; dim(KICH_Exp01); table(colnames(KICH_Exp01)==KICH_c$ID); KICH_Met01[1:5,1:5]; dim(KICH_Met01); table(colnames(KICH_Met01)==KICH_c$ID); head(KICH_c_SG)
KIRC_Exp01[1:5,1:5]; dim(KIRC_Exp01); table(colnames(KIRC_Exp01)==KIRC_c$ID); KIRC_Met01[1:5,1:5]; dim(KIRC_Met01); table(colnames(KIRC_Met01)==KIRC_c$ID); head(KIRC_c_SG)
KIRP_Exp01[1:5,1:5]; dim(KIRP_Exp01); table(colnames(KIRP_Exp01)==KIRP_c$ID); KIRP_Met01[1:5,1:5]; dim(KIRP_Met01); table(colnames(KIRP_Met01)==KIRP_c$ID); head(KIRP_c_SG)
LIHC_Exp01[1:5,1:5]; dim(LIHC_Exp01); table(colnames(LIHC_Exp01)==LIHC_c$ID); LIHC_Met01[1:5,1:5]; dim(LIHC_Met01); table(colnames(LIHC_Met01)==LIHC_c$ID); head(LIHC_c_SG)
LUAD_Exp01[1:5,1:5]; dim(LUAD_Exp01); table(colnames(LUAD_Exp01)==LUAD_c$ID); LUAD_Met01[1:5,1:5]; dim(LUAD_Met01); table(colnames(LUAD_Met01)==LUAD_c$ID); head(LUAD_c_SG)
LUSC_Exp01[1:5,1:5]; dim(LUSC_Exp01); table(colnames(LUSC_Exp01)==LUSC_c$ID); LUSC_Met01[1:5,1:5]; dim(LUSC_Met01); table(colnames(LUSC_Met01)==LUSC_c$ID); head(LUSC_c_SG)
MESO_Exp01[1:5,1:5]; dim(MESO_Exp01); table(colnames(MESO_Exp01)==MESO_c$ID); MESO_Met01[1:5,1:5]; dim(MESO_Met01); table(colnames(MESO_Met01)==MESO_c$ID); head(MESO_c_SG)
PAAD_Exp01[1:5,1:5]; dim(PAAD_Exp01); table(colnames(PAAD_Exp01)==PAAD_c$ID); PAAD_Met01[1:5,1:5]; dim(PAAD_Met01); table(colnames(PAAD_Met01)==PAAD_c$ID); head(PAAD_c_SG)
READ_Exp01[1:5,1:5]; dim(READ_Exp01); table(colnames(READ_Exp01)==READ_c$ID); READ_Met01[1:5,1:5]; dim(READ_Met01); table(colnames(READ_Met01)==READ_c$ID); head(READ_c_SG)
SKCM_Exp01[1:5,1:5]; dim(SKCM_Exp01); table(colnames(SKCM_Exp01)==SKCM_c$ID); SKCM_Met01[1:5,1:5]; dim(SKCM_Met01); table(colnames(SKCM_Met01)==SKCM_c$ID); head(SKCM_c_SG)
STAD_Exp01[1:5,1:5]; dim(STAD_Exp01); table(colnames(STAD_Exp01)==STAD_c$ID); STAD_Met01[1:5,1:5]; dim(STAD_Met01); table(colnames(STAD_Met01)==STAD_c$ID); head(STAD_c_SG)
TGCT_Exp01[1:5,1:5]; dim(TGCT_Exp01); table(colnames(TGCT_Exp01)==TGCT_c$ID); TGCT_Met01[1:5,1:5]; dim(TGCT_Met01); table(colnames(TGCT_Met01)==TGCT_c$ID); head(TGCT_c_SG)
THCA_Exp01[1:5,1:5]; dim(THCA_Exp01); table(colnames(THCA_Exp01)==THCA_c$ID); THCA_Met01[1:5,1:5]; dim(THCA_Met01); table(colnames(THCA_Met01)==THCA_c$ID); head(THCA_c_SG)
THYM_Exp01[1:5,1:5]; dim(THYM_Exp01); table(colnames(THYM_Exp01)==THYM_c$ID); THYM_Met01[1:5,1:5]; dim(THYM_Met01); table(colnames(THYM_Met01)==THYM_c$ID); head(THYM_c_SG)
UCEC_Exp01[1:5,1:5]; dim(UCEC_Exp01); table(colnames(UCEC_Exp01)==UCEC_c$ID); UCEC_Met01[1:5,1:5]; dim(UCEC_Met01); table(colnames(UCEC_Met01)==UCEC_c$ID); head(UCEC_c_SG)
UCS__Exp01[1:5,1:5]; dim(UCS__Exp01); table(colnames(UCS__Exp01)==UCS__c$ID); UCS__Met01[1:5,1:5]; dim(UCS__Met01); table(colnames(UCS__Met01)==UCS__c$ID); head(UCS__c_SG)
UVM__Exp01[1:5,1:5]; dim(UVM__Exp01); table(colnames(UVM__Exp01)==UVM__c$ID); UVM__Met01[1:5,1:5]; dim(UVM__Met01); table(colnames(UVM__Met01)==UVM__c$ID); head(UVM__c_SG)

prefixes <- unlist(strsplit("ACC__ BLCA_ BRCA_ CESC_ CHOL_ COAD_ DLBC_ ESCA_ HNSC_ KICH_ KIRC_ KIRP_ LIHC_ LUAD_ LUSC_ MESO_ PAAD_ READ_ SKCM_ STAD_ TGCT_ THCA_ THYM_ UCEC_ UCS__ UVM__", " ")); res_all <- NULL
for(ca in prefixes){ c_df <- get(paste0(ca, "c_SG")); for(dt in c("Exp", "Met")){ dat <- get(paste0(ca, dt, "01")); g1 <- c_df$ID[which(c_df$Stage=="L")]; g1 <- g1[g1 %in% colnames(dat)]; g2 <- c_df$ID[which(c_df$Stage=="M")]; g2 <- g2[g2 %in% colnames(dat)]; res <- data.frame(ID=rownames(dat), log2FC=NA, P_value=NA, FDR=NA, Cancer=gsub("_", "", ca), Type=dt, stringsAsFactors=FALSE); for(i in 1:nrow(dat)){ v1 <- as.numeric(dat[i, g1]); v2 <- as.numeric(dat[i, g2]); if(sum(!is.na(v1))>1 & sum(!is.na(v2))>1){ if(var(v1, na.rm=TRUE)>0 | var(v2, na.rm=TRUE)>0){ res$log2FC[i] <- log2((mean(v2, na.rm=TRUE)+1e-6)/(mean(v1, na.rm=TRUE)+1e-6)); res$P_value[i] <- t.test(v2, v1)$p.value } } }; res$FDR <- p.adjust(res$P_value, method="fdr"); res_all <- rbind(res_all, res) } }
write.table(res_all, "Output_Ox11_OXI/23_DEG_Stage_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)

head(res_all); dim(res_all); table(res_all[,c(6,5)])

pdf("Output_Ox11_OXI/23_DEG_Heatmaps.pdf", width=8, height=10); cancers <- unique(res_all$Cancer)
for(dt in c("Exp", "Met")){ sub <- res_all[res_all$Type==dt, ]; ids <- unique(sub$ID); m_fc <- matrix(NA, length(cancers), length(ids), dimnames=list(cancers, ids)); m_fdr <- m_fc; for(i in 1:nrow(sub)){ m_fc[sub$Cancer[i], sub$ID[i]] <- sub$log2FC[i]; m_fdr[sub$Cancer[i], sub$ID[i]] <- sub$FDR[i] }; m_sig <- matrix(ifelse(!is.na(m_fdr) & m_fdr < 0.05, "*", ""), nrow(m_fdr), dimnames=dimnames(m_fdr)); mx <- max(abs(m_fc), na.rm=TRUE); mx <- ifelse(is.infinite(mx) | mx==0, 1, mx); pheatmap(m_fc, color=colorRampPalette(c("blue", "white", "red"))(50), breaks=seq(-mx, mx, length.out=51), display_numbers=m_sig, cluster_rows=FALSE, cluster_cols=FALSE, na_col="gray90", main=paste(dt, "log2FC Heatmap"), fontsize_number=15, number_color="black") }
dev.off()

r_m <- res_all[res_all$Type=="Met" & !is.na(res_all$FDR) & res_all$FDR<0.05, c("Cancer", "ID", "log2FC", "FDR")]; r_e <- res_all[res_all$Type=="Exp" & !is.na(res_all$FDR) & res_all$FDR<0.05, c("Cancer", "ID", "log2FC", "FDR")]; cpg <- CpG_Table[, c("ID", "UCSC_RefGene_Name")]; colnames(r_m)[2] <- "CpG"; colnames(cpg) <- c("CpG", "Gene"); r_m <- merge(r_m, cpg, by="CpG"); colnames(r_e)[2] <- "Gene"; res_pair <- merge(r_m, r_e, by=c("Cancer", "Gene"), suffixes=c("_Met", "_Exp")); write.table(res_pair, "Output_Ox11_OXI/23_Significant_Pairs.txt", sep="\t", quote=FALSE, row.names=FALSE)
head(res_pair); dim(res_pair)
################################################################################################
#Module 01. Data loading
################################################################################################
#Loading KoGES methylation: Takes 15 mins
Sys.time()
source("Loading_KoGES_Met_ASAS0400_ASAS1528_CITY0822_Epid_QC.R")
#Loading TCGA methylation: Takes 15 mins
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

#INT_kex_CpG <- Reduce(intersect, list(rownames(Met_ASAS0400),rownames(Met_ASAS1528),rownames(Met_CITY0822),rownames(GSE223748m)))
#write.table(INT_kex_CpG,file="INT_kex_CpG_QC.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
INT_kex_CpG <- as.vector(unlist(read.table("INT_kex_CpG_QC.txt")))
length(INT_kex_CpG) #5306 > 5213


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

COMMON_A0105 <- intersect(colnames(Met_ASAS0400b),colnames(Met_ASAS1528b)); length(COMMON_A0105)

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
A01_m04_OXI_M_Info <- A01_m04_OXI_M
A01_m04_OXI_F_Info <- A01_m04_OXI_F
C01_m08_OXI_M_Info <- C01_m08_OXI_M
C01_m08_OXI_F_Info <- C01_m08_OXI_F

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

A01_m04_OXI_M$NA01_OBS <- rowSums(A01_m04_OXI_M); table(A01_m04_OXI_M$NA01_OBS); A01_m04_ALL_OXI_M <- A01_m04_OXI_M
A01_m04_OXI_F$NA01_OBS <- rowSums(A01_m04_OXI_F); table(A01_m04_OXI_F$NA01_OBS); A01_m04_ALL_OXI_F <- A01_m04_OXI_F
C01_m08_OXI_M$NC01_OBS <- rowSums(C01_m08_OXI_M); table(C01_m08_OXI_M$NC01_OBS); C01_m08_ALL_OXI_M <- C01_m08_OXI_M
C01_m08_OXI_F$NC01_OBS <- rowSums(C01_m08_OXI_F); table(C01_m08_OXI_F$NC01_OBS); C01_m08_ALL_OXI_F <- C01_m08_OXI_F
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
m04mA01_OXI_M <- dplyr::select(Met_ASAS0400m[rownames(Met_ASAS0400m) %in% INT_kex_CpG, ], intersect(colnames(Met_ASAS0400m), dplyr::filter(ASAS01to10, AS01_SEX==1)$DIST_ID)); m04mA01_OXI_M <- m04mA01_OXI_M[complete.cases(m04mA01_OXI_M), ]; m04mA01_OXI_M <- as.data.frame(as.data.table(m04mA01_OXI_M, keep.rownames = "ProbeID")); rownames(m04mA01_OXI_M) <- m04mA01_OXI_M$ProbeID; dim(m04mA01_OXI_M) #5306,201
m04mA01_OXI_F <- dplyr::select(Met_ASAS0400m[rownames(Met_ASAS0400m) %in% INT_kex_CpG, ], intersect(colnames(Met_ASAS0400m), dplyr::filter(ASAS01to10, AS01_SEX==2)$DIST_ID)); m04mA01_OXI_F <- m04mA01_OXI_F[complete.cases(m04mA01_OXI_F), ]; m04mA01_OXI_F <- as.data.frame(as.data.table(m04mA01_OXI_F, keep.rownames = "ProbeID")); rownames(m04mA01_OXI_F) <- m04mA01_OXI_F$ProbeID; dim(m04mA01_OXI_F) #5306,201
m08mC01_OXI_M <- dplyr::select(Met_CITY0822m[rownames(Met_CITY0822m) %in% INT_kex_CpG, ], intersect(colnames(Met_CITY0822m), dplyr::filter(CITY01to02, CT01_SEX==1)$DIST_ID)); m08mC01_OXI_M <- m08mC01_OXI_M[complete.cases(m08mC01_OXI_M), ]; m08mC01_OXI_M <- as.data.frame(as.data.table(m08mC01_OXI_M, keep.rownames = "ProbeID")); rownames(m08mC01_OXI_M) <- m08mC01_OXI_M$ProbeID; dim(m08mC01_OXI_M) #5306,623
m08mC01_OXI_F <- dplyr::select(Met_CITY0822m[rownames(Met_CITY0822m) %in% INT_kex_CpG, ], intersect(colnames(Met_CITY0822m), dplyr::filter(CITY01to02, CT01_SEX==2)$DIST_ID)); m08mC01_OXI_F <- m08mC01_OXI_F[complete.cases(m08mC01_OXI_F), ]; m08mC01_OXI_F <- as.data.frame(as.data.table(m08mC01_OXI_F, keep.rownames = "ProbeID")); rownames(m08mC01_OXI_F) <- m08mC01_OXI_F$ProbeID; dim(m08mC01_OXI_F) #5306,201
m04mA01_OXI_M <- m04mA01_OXI_M[order(m04mA01_OXI_M$ProbeID), ]
m04mA01_OXI_F <- m04mA01_OXI_F[order(m04mA01_OXI_F$ProbeID), ]; table(m04mA01_OXI_M$ProbeID==m04mA01_OXI_F$ProbeID) #5306
m08mC01_OXI_M <- m08mC01_OXI_M[order(m08mC01_OXI_M$ProbeID), ]; table(m04mA01_OXI_M$ProbeID==m08mC01_OXI_M$ProbeID) #5306
m08mC01_OXI_F <- m08mC01_OXI_F[order(m08mC01_OXI_F$ProbeID), ]; table(m04mA01_OXI_M$ProbeID==m08mC01_OXI_F$ProbeID) #5306
m04mA01_OXI_M[1:5, 1:5]; dim(m04mA01_OXI_M) #5213,196
m04mA01_OXI_F[1:5, 1:5]; dim(m04mA01_OXI_F) #5213,197
m08mC01_OXI_M[1:5, 1:5]; dim(m08mC01_OXI_M) #5213,617
m08mC01_OXI_F[1:5, 1:5]; dim(m08mC01_OXI_F) #5213,197


A01_M_cID <- sort(intersect(colnames(m04mA01_OXI_M),rownames(A01_m04_ALL_OXI_M))); length(A01_M_cID) #169
A01_F_cID <- sort(intersect(colnames(m04mA01_OXI_F),rownames(A01_m04_ALL_OXI_F))); length(A01_F_cID) #173
C01_M_cID <- sort(intersect(colnames(m08mC01_OXI_M),rownames(C01_m08_ALL_OXI_M))); length(C01_M_cID) #615
C01_F_cID <- sort(intersect(colnames(m08mC01_OXI_F),rownames(C01_m08_ALL_OXI_F))); length(C01_F_cID) #194

A01_m04_ALL_OXI_M <- A01_m04_ALL_OXI_M[A01_M_cID,]
A01_m04_ALL_OXI_F <- A01_m04_ALL_OXI_F[A01_F_cID,]
C01_m08_ALL_OXI_M <- C01_m08_ALL_OXI_M[C01_M_cID,]
C01_m08_ALL_OXI_F <- C01_m08_ALL_OXI_F[C01_F_cID,]
m04mA01_OXI_M <- dplyr::select(m04mA01_OXI_M, A01_M_cID)
m04mA01_OXI_F <- dplyr::select(m04mA01_OXI_F, A01_F_cID)
m08mC01_OXI_M <- dplyr::select(m08mC01_OXI_M, C01_M_cID)
m08mC01_OXI_F <- dplyr::select(m08mC01_OXI_F, C01_F_cID)
A01_M_CT <- ASAS0400c[A01_M_cID,]
A01_F_CT <- ASAS0400c[A01_F_cID,]
C01_M_CT <- CITY0822c[C01_M_cID,]
C01_F_CT <- CITY0822c[C01_F_cID,]
A01_m04_OXI_M_Info <- A01_m04_OXI_M_Info[A01_M_cID,]
A01_m04_OXI_F_Info <- A01_m04_OXI_F_Info[A01_F_cID,]
C01_m08_OXI_M_Info <- C01_m08_OXI_M_Info[C01_M_cID,]
C01_m08_OXI_F_Info <- C01_m08_OXI_F_Info[C01_F_cID,]





head(A01_m04_ALL_OXI_M); dim(A01_m04_ALL_OXI_M)
head(A01_m04_ALL_OXI_F); dim(A01_m04_ALL_OXI_F)
head(C01_m08_ALL_OXI_M); dim(C01_m08_ALL_OXI_M)
head(C01_m08_ALL_OXI_F); dim(C01_m08_ALL_OXI_F)
m04mA01_OXI_M[1:5, 1:5]; dim(m04mA01_OXI_M) #5213,169
m04mA01_OXI_F[1:5, 1:5]; dim(m04mA01_OXI_F) #5213,173
m08mC01_OXI_M[1:5, 1:5]; dim(m08mC01_OXI_M) #5213,615
m08mC01_OXI_F[1:5, 1:5]; dim(m08mC01_OXI_F) #5213,194
head(A01_M_CT); dim(A01_M_CT)
head(A01_F_CT); dim(A01_F_CT)
head(C01_M_CT); dim(C01_M_CT)
head(C01_F_CT); dim(C01_F_CT)
table(rownames(A01_m04_ALL_OXI_M)==colnames(m04mA01_OXI_M))
table(rownames(A01_m04_ALL_OXI_F)==colnames(m04mA01_OXI_F))
table(rownames(C01_m08_ALL_OXI_M)==colnames(m08mC01_OXI_M))
table(rownames(C01_m08_ALL_OXI_F)==colnames(m08mC01_OXI_F))
table(rownames(A01_m04_ALL_OXI_M)==rownames(A01_M_CT))
table(rownames(A01_m04_ALL_OXI_F)==rownames(A01_F_CT))
table(rownames(C01_m08_ALL_OXI_M)==rownames(C01_M_CT))
table(rownames(C01_m08_ALL_OXI_F)==rownames(C01_F_CT))
table(rownames(A01_m04_ALL_OXI_M)==rownames(A01_m04_OXI_M_Info))
table(rownames(A01_m04_ALL_OXI_F)==rownames(A01_m04_OXI_F_Info))
table(rownames(C01_m08_ALL_OXI_M)==rownames(C01_m08_OXI_M_Info))
table(rownames(C01_m08_ALL_OXI_F)==rownames(C01_m08_OXI_F_Info))





CRPV_OBS <- as.data.frame(fread("Output_Ox11_OXI/CRPVm_Ox01_OXI_OBS.txt.gz",  header = TRUE, sep = "\t")); CRPV_OBS <- dplyr::filter(CRPV_OBS, CRp * CRs > 0); dim(CRPV_OBS)
CRPV_OBSx <- dplyr::filter(CRPV_OBS, abs(CRp) > 0.23 & abs(CRs) > 0.23 & PVp < 0.003 & PVs < 0.003 & substr(CLASS,5,7)=="A01"); table(CRPV_OBSx$CLASS); nrow(CRPV_OBSx)
CRPV_OBSy <- dplyr::filter(CRPV_OBS, abs(CRp) > 0.23 & abs(CRs) > 0.23 & PVp < 0.003 & PVs < 0.003 & substr(CLASS,5,7)=="C01"); table(CRPV_OBSy$CLASS); nrow(CRPV_OBSy)
intersect(CRPV_OBSx$ID,CRPV_OBSy$ID)
head(CRPV_OBSx); dim(CRPV_OBSx)
head(CRPV_OBSy); dim(CRPV_OBSy)

res <- data.frame()
dAM <- merge(A01_m04_ALL_OXI_M, A01_M_CT, by=0); rownames(dAM) <- dAM$Row.names; cAM <- CRPV_OBSx$ID[CRPV_OBSx$CLASS=="m04mA01_OXI_M"]
dAF <- merge(A01_m04_ALL_OXI_F, A01_F_CT, by=0); rownames(dAF) <- dAF$Row.names; cAF <- CRPV_OBSx$ID[CRPV_OBSx$CLASS=="m04mA01_OXI_F"]
dCF <- merge(C01_m08_ALL_OXI_F, C01_F_CT, by=0); rownames(dCF) <- dCF$Row.names; cCF <- CRPV_OBSy$ID[CRPV_OBSy$CLASS=="m08mC01_OXI_F"]

for(cg in cAM){ m <- lm(as.numeric(m04mA01_OXI_M[cg, rownames(dAM)]) ~ NA01_OBS + CD8T + CD4T + NK + Bcell + Mono, data=dAM); res <- rbind(res, data.frame(ID=cg, Cohort="ASAS_M", Coef=summary(m)$coef[2,1], PV=summary(m)$coef[2,4])) }
for(cg in cAF){ m <- lm(as.numeric(m04mA01_OXI_F[cg, rownames(dAF)]) ~ NA01_OBS + CD8T + CD4T + NK + Bcell + Mono, data=dAF); res <- rbind(res, data.frame(ID=cg, Cohort="ASAS_F", Coef=summary(m)$coef[2,1], PV=summary(m)$coef[2,4])) }
for(cg in cCF){ m <- lm(as.numeric(m08mC01_OXI_F[cg, rownames(dCF)]) ~ NC01_OBS + CD8T + CD4T + NK + Bcell + Mono, data=dCF); res <- rbind(res, data.frame(ID=cg, Cohort="CITY_F", Coef=summary(m)$coef[2,1], PV=summary(m)$coef[2,4])) }
res$FDR <- p.adjust(res$PV, method="fdr")
table(res$CLASS)
head(res); dim(res)

CRPV_OBSx <- merge(x=CRPV_OBSx,y=res,by="ID")
head(CRPV_OBSx); dim(CRPV_OBSx)


CRPV_OBSx$HyPi <- (CRPV_OBSx$CRp * -log10(CRPV_OBSx$PVp)) + (CRPV_OBSx$CRs * -log10(CRPV_OBSx$PVs))
CRPV_OBSm <- dplyr::filter(CRPV_OBSx, substr(CLASS,13,13)=="M"); CRPV_OBSf <- dplyr::filter(CRPV_OBSx, substr(CLASS,13,13)=="F"); intersect(CRPV_OBSm$ID,CRPV_OBSf$ID)
CRPV_OBSx$Sex <- substr(CRPV_OBSx$CLASS,13,13)
head(CRPV_OBSx); dim(CRPV_OBSx) #26,13
table(duplicated(CRPV_OBSx$ID))

m04mA01_OXI_M$ProbeID <- rownames(m04mA01_OXI_M)
m04mA01_OXI_F$ProbeID <- rownames(m04mA01_OXI_F)
m04mA01_OXI_M_OBS <- as.data.frame(t(dplyr::select(dplyr::filter(m04mA01_OXI_M, ProbeID %in% CRPV_OBSm$ID), -ProbeID))); dim(m04mA01_OXI_M_OBS) #169,12
m04mA01_OXI_F_OBS <- as.data.frame(t(dplyr::select(dplyr::filter(m04mA01_OXI_F, ProbeID %in% CRPV_OBSf$ID), -ProbeID))); dim(m04mA01_OXI_F_OBS) #173,14
m04mA01_OXI_M_OBS <- merge(x=A01_m04_OXI_M, y=m04mA01_OXI_M_OBS,by=0); rownames(m04mA01_OXI_M_OBS) <- m04mA01_OXI_M_OBS$Row.names; m04mA01_OXI_M_OBS <- m04mA01_OXI_M_OBS[,-1]
m04mA01_OXI_F_OBS <- merge(x=A01_m04_OXI_F, y=m04mA01_OXI_F_OBS,by=0); rownames(m04mA01_OXI_F_OBS) <- m04mA01_OXI_F_OBS$Row.names; m04mA01_OXI_F_OBS <- m04mA01_OXI_F_OBS[,-1]
m04mA01_OXI_M$ProbeID <- NULL
m04mA01_OXI_F$ProbeID <- NULL
head(m04mA01_OXI_M_OBS); dim(m04mA01_OXI_M_OBS)
head(m04mA01_OXI_F_OBS); dim(m04mA01_OXI_F_OBS)

hist(m04mA01_OXI_M_OBS$NA01_OBS)

head(GPL13534); dim(GPL13534) #485577,7
head(CRPV_OBS); dim(CRPV_OBS) #19196,9
length(intersect(CRPV_OBS$ID,GPL13534$ID)) #5212
head(CRPV_OBSx); dim(CRPV_OBSx) #26,13
length(intersect(CRPV_OBSx$ID,GPL13534$ID)) #26
m04mA01_OXI_M[1:5,1:5]; dim(m04mA01_OXI_M)
m04mA01_OXI_F[1:5,1:5]; dim(m04mA01_OXI_F)


GPL13534x <- GPL13534
GPL13534x$SEX <- ifelse(GPL13534x$ID %in% colnames(m04mA01_OXI_M_OBS), "M", ifelse(GPL13534x$ID %in% colnames(m04mA01_OXI_F_OBS), "F","XXX"))
GPL13534x <- merge(x=GPL13534x,y=CRPV_OBSx[,c(1,9:12)])
GPL13534x <- dplyr::filter(GPL13534x, SEX != "XXX")
table(GPL13534x$SEX)
GPL13534x$UCSC_RefGene_Name <- sub(";.*", "", GPL13534x$UCSC_RefGene_Name); GPL13534x$UCSC_RefGene_Group <- sub(";.*", "", GPL13534x$UCSC_RefGene_Group)
GPL13534x[, c(9,12)] <- lapply(GPL13534x[, c(9,12)], function(x) round(as.numeric(x), 3))
GPL13534x[,    9:11] <- lapply(GPL13534x[,    9:11], function(x) round(as.numeric(x), 5))

head(GPL13534x); dim(GPL13534x)
write.table(GPL13534x, file="Output_Ox11_OXI/01_CpG_Table.txt", sep="\t",quote=FALSE,row.names=FALSE)



A01_m04_OXI_M_Info$DIST_ID <- rownames(A01_m04_OXI_M_Info); A01_m04_OXI_M_Info <- merge(x=A01_m04_OXI_M_Info,y=dplyr::select(ASAS01to10,DIST_ID,AS01_AGE), by="DIST_ID"); rownames(A01_m04_OXI_M_Info) <- A01_m04_OXI_M_Info$DIST_ID; A01_m04_OXI_M_Info$DIST_ID <- NULL
A01_m04_OXI_F_Info$DIST_ID <- rownames(A01_m04_OXI_F_Info); A01_m04_OXI_F_Info <- merge(x=A01_m04_OXI_F_Info,y=dplyr::select(ASAS01to10,DIST_ID,AS01_AGE), by="DIST_ID"); rownames(A01_m04_OXI_F_Info) <- A01_m04_OXI_F_Info$DIST_ID; A01_m04_OXI_F_Info$DIST_ID <- NULL
head(A01_m04_OXI_M_Info); table(A01_m04_OXI_M_Info$NA01_SM_SMK)
head(A01_m04_OXI_F_Info); table(A01_m04_OXI_F_Info$NA01_SM_SMK)


#Plot 01: CpG-OBS correlation plots for 14 CpG sites
plist <- list(); GPL13534x <- GPL13534x[order(GPL13534x$SEX, decreasing=TRUE), ]
for(i in 1:nrow(GPL13534x)){ cg <- GPL13534x$ID[i]; gn <- GPL13534x$UCSC_RefGene_Name[i]; gp <- GPL13534x$UCSC_RefGene_Group[i]; sx <- GPL13534x$SEX[i]
tit <- paste0(ifelse(sx=="M", "[Male] ", "[Female] "), cg, ifelse(is.na(gn) | gn=="", "", paste0("\n(", gn, ", ", gp, ")")))
d <- if(sx=="M") m04mA01_OXI_M_OBS else m04mA01_OXI_F_OBS; pcol <- if(sx=="M") "steelblue" else "darkgreen"; lcol <- if(sx=="M") "red" else "orange"; d$Y <- d[,cg]
plist[[i]] <- ggplot(d, aes(x=NA01_OBS, y=Y)) + geom_point(color=pcol, size=1.5, alpha=0.7) + geom_smooth(method="lm", color=lcol, se=TRUE) + theme_bw() + labs(title=tit, x="OBS", y="M-value") + theme(plot.title=element_text(size=11, face="bold")) }
ggsave("Output_Ox11_OXI/01_CpG_OBS_Corr_26.pdf", plot=do.call(grid.arrange, c(plist, ncol=5)), width=20, height=18)

pM <- ggplot(m04mA01_OXI_M_OBS, aes(x=NA01_OBS)) + geom_histogram(binwidth=1, fill="steelblue", color="black", alpha=0.7) + theme_minimal() + labs(title="OBS Distribution (Male)", x="Oxidative Balance Score (OBS)", y="Frequency")
pF <- ggplot(m04mA01_OXI_F_OBS, aes(x=NA01_OBS)) + geom_histogram(binwidth=1, fill="darkgreen", color="black", alpha=0.7) + theme_minimal() + labs(title="OBS Distribution (Female)", x="Oxidative Balance Score (OBS)", y="Frequency")
ggsave("Output_Ox11_OXI/01_OBS_Histogram.pdf", plot=grid.arrange(pM, pF, ncol=2), width=10, height=5)

cv <- setdiff(colnames(A01_m04_OXI_M_Info), "NA01_SM_SMK"); out <- data.frame(Variable="N", Male=as.character(nrow(A01_m04_OXI_M_Info)), Female=as.character(nrow(A01_m04_OXI_F_Info)))
for(v in cv){ out <- rbind(out, data.frame(Variable=v, Male=paste0(round(mean(A01_m04_OXI_M_Info[,v], na.rm=TRUE),1), "\u00b1", round(sd(A01_m04_OXI_M_Info[,v], na.rm=TRUE),1)), Female=paste0(round(mean(A01_m04_OXI_F_Info[,v], na.rm=TRUE),1), "\u00b1", round(sd(A01_m04_OXI_F_Info[,v], na.rm=TRUE),1)))) }
for(i in 1:3){ out <- rbind(out, data.frame(Variable=paste0("NA01_SM_SMK_", i), Male=as.character(sum(A01_m04_OXI_M_Info$NA01_SM_SMK==i, na.rm=TRUE)), Female=as.character(sum(A01_m04_OXI_F_Info$NA01_SM_SMK==i, na.rm=TRUE)))) }
write.table(out, "Output_Ox11_OXI/01_Table1_Participants.txt", sep="\t", quote=FALSE, row.names=FALSE, fileEncoding="CP949")


CRPV_Suba <- CRPV_OBSx[, c("ID", "CRp", "CRs", "Sex")]; CRPV_Suba <- dplyr::filter(CRPV_Suba, ID %in% c(colnames(m04mA01_OXI_M_OBS),colnames(m04mA01_OXI_F_OBS))); dim(CRPV_Suba); table(CRPV_Suba$Sex)
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

table(rownames(m04mA01_OXI_M)==rownames(m04mA01_OXI_F))
GSE223748m$ID <- rownames(GSE223748m); GSE223748m <- dplyr::filter(GSE223748m, ID %in% rownames(m04mA01_OXI_M)); GSE223748m$ID <- NULL
GSE223748m[1:5, 1:5]; dim(GSE223748m) #37554->5213,15043
GSE223748i[1:5, 1:5]; dim(GSE223748i) #15043,12
GSE223748i$Sex <- as.numeric(GSE223748i$Sex); table(GSE223748i$Sex)
GSE223748i <- dplyr::filter(GSE223748i, Sex != 3)
head(GSE223748i); dim(GSE223748i); table(GSE223748i$Sex) #14667,12 / 1 6842, 2 7825

Blood_SciNameM <- table(dplyr::filter(GSE223748i, Tissue=="Blood"&Sex==1)$SciName); Blood_SciNameM <- Blood_SciNameM[Blood_SciNameM > 49]; Blood_SciNameM
Blood_SciNameF <- table(dplyr::filter(GSE223748i, Tissue=="Blood"&Sex==2)$SciName); Blood_SciNameF <- Blood_SciNameF[Blood_SciNameF > 49]; Blood_SciNameF
COMMON_Blood_SciName <- intersect(names(Blood_SciNameM),names(Blood_SciNameF)); length(COMMON_Blood_SciName) #9
COMMON_Blood_SciName

GSE223748iM <- dplyr::filter(GSE223748i, Sex==1 & Tissue=="Blood" & SciName %in% COMMON_Blood_SciName)[,2:3]
GSE223748iF <- dplyr::filter(GSE223748i, Sex==2 & Tissue=="Blood" & SciName %in% COMMON_Blood_SciName)[,2:3]
GSE223748mM <- GSE223748m; GSE223748mM$ID <- rownames(GSE223748mM); GSE223748mM <- dplyr::select(dplyr::filter(GSE223748mM, ID %in% rownames(Wide_CRPa_M)), rownames(GSE223748iM))
GSE223748mF <- GSE223748m; GSE223748mF$ID <- rownames(GSE223748mF); GSE223748mF <- dplyr::select(dplyr::filter(GSE223748mF, ID %in% rownames(Wide_CRPa_F)), rownames(GSE223748iF))
GSE223748mM[1:5,1:5]; dim(GSE223748mM) #12,1257
GSE223748mF[1:5,1:5]; dim(GSE223748mF) #14,1252
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







#Plot 08: Epigenetic Distance Matrix by Euclidean
Sp_COL_NAMES <- c("C. familiaris","C. capreolus","E. caballus","F. catus","H. sapiens","M. mulatta","M. musculus","O. aries","R. norvegicus")
colnames(GSE223748mM_Pr) <- Sp_COL_NAMES; colnames(GSE223748mF_Pr) <- Sp_COL_NAMES; colnames(GSE223748mM_Sp) <- Sp_COL_NAMES; colnames(GSE223748mF_Sp) <- Sp_COL_NAMES
Dist_M_Pra <- dist(t(GSE223748mM_Pr)); Dist_M_Spa <- dist(t(GSE223748mM_Sp)); Dist_F_Pra <- dist(t(GSE223748mF_Pr)); Dist_F_Spa <- dist(t(GSE223748mF_Sp))
Tree_M_Pra <- nj(Dist_M_Pra); Tree_M_Spa <- nj(Dist_M_Spa); Tree_F_Pra <- nj(Dist_F_Pra); Tree_F_Spa <- nj(Dist_F_Spa)
B <- 1000; bM_Pr <- list(); bM_Sp <- list(); bF_Pr <- list(); bF_Sp <- list()
for(i in 1:B){ idx <- sample(1:nrow(GSE223748mM_Pr), replace=TRUE); bM_Pr[[i]] <- nj(dist(t(GSE223748mM_Pr[idx, ]))); bM_Sp[[i]] <- nj(dist(t(GSE223748mM_Sp[idx, ]))); bF_Pr[[i]] <- nj(dist(t(GSE223748mF_Pr[idx, ]))); bF_Sp[[i]] <- nj(dist(t(GSE223748mF_Sp[idx, ]))) }
class(bM_Pr) <- "multiPhylo"; class(bM_Sp) <- "multiPhylo"; class(bF_Pr) <- "multiPhylo"; class(bF_Sp) <- "multiPhylo"
bp_M_Pra <- round(prop.clades(Tree_M_Pra, bM_Pr)/B*100); bp_M_Spa <- round(prop.clades(Tree_M_Spa, bM_Sp)/B*100); bp_F_Pra <- round(prop.clades(Tree_F_Pra, bF_Pr)/B*100); bp_F_Spa <- round(prop.clades(Tree_F_Spa, bF_Sp)/B*100)

#Plot 09: Tanglegram: Male vs Female, compare phylogeny structure by Entanglement
pdf("Output_Ox11_OXI/08_EpigeneticTree_Bootstrap_Fixed.pdf", width=12, height=12); par(mfrow=c(2,2), mar=c(2,2,4,2))
plot(Tree_M_Pra, main="Epigenetic Tree: Male Pearson", type="unrooted", lab4ut="axial", cex=0.9, edge.width=2); nodelabels(bp_M_Pra, cex=0.75, frame="none", col="royalblue", adj=c(1.2, -0.5))
plot(Tree_M_Spa, main="Epigenetic Tree: Male Spearman", type="unrooted", lab4ut="axial", cex=0.9, edge.width=2); nodelabels(bp_M_Spa, cex=0.75, frame="none", col="royalblue", adj=c(1.2, -0.5))
plot(Tree_F_Pra, main="Epigenetic Tree: Female Pearson", type="unrooted", lab4ut="axial", cex=0.9, edge.width=2); nodelabels(bp_F_Pra, cex=0.75, frame="none", col="royalblue", adj=c(1.2, -0.5))
plot(Tree_F_Spa, main="Epigenetic Tree: Female Spearman", type="unrooted", lab4ut="axial", cex=0.9, edge.width=2); nodelabels(bp_F_Spa, cex=0.75, frame="none", col="royalblue", adj=c(1.2, -0.5)); dev.off()


pdf("Output_Ox11_OXI/09_Tanglegram.pdf", width=8, height=4); par(mfrow=c(1,1))
Dend_M_Pra <- as.dendrogram(hclust(Dist_M_Pra)); Dend_F_Pra <- as.dendrogram(hclust(Dist_F_Pra)); labels(Dend_M_Pra) <- labels(Dend_F_Pra)[match(labels(Dend_M_Pra), labels(Dend_F_Pra))]
tanglegram(Dend_M_Pra, Dend_F_Pra, main_left="Male Pearson Tree", main_right="Female Pearson Tree", lab.cex=1, edge.lwd=2, margin_inner=7, columns_width=c(5,2,5)); dev.off()



cM <- matrix(0, nrow(GSE223748m), 9); rownames(cM) <- rownames(GSE223748m); colnames(cM) <- Sp_COL_NAMES; cF <- cM
FN <- c("Canis lupus familiaris", "Capreolus capreolus", "Equus caballus", "Felis catus", "Homo sapiens", "Macaca mulatta", "Mus musculus", "Ovis aries", "Rattus norvegicus")
for(j in 1:9){ iM <- intersect(rownames(GSE223748i)[GSE223748i$Sex==1 & GSE223748i$Tissue=="Blood" & GSE223748i$SciName==FN[j]], colnames(GSE223748m)); iF <- intersect(rownames(GSE223748i)[GSE223748i$Sex==2 & GSE223748i$Tissue=="Blood" & GSE223748i$SciName==FN[j]], colnames(GSE223748m)); if(length(iM)>2) cM[,j] <- cor(t(GSE223748m[, iM]), GSE223748i[iM, "Age"], use="pairwise.complete.obs"); if(length(iF)>2) cF[,j] <- cor(t(GSE223748m[, iF]), GSE223748i[iF, "Age"], use="pairwise.complete.obs") }
cM[is.na(cM)] <- 0; cF[is.na(cF)] <- 0; tM_T <- nj(dist(t(cM[cAM, ]))); tF_T <- nj(dist(t(cF[cAF, ]))); B <- 1000; rfM <- numeric(B); rfF <- numeric(B)
for(i in 1:B){ rfM[i] <- dist.topo(tM_T, nj(dist(t(cM[sample(nrow(cM), length(cAM)), ])))); rfF[i] <- dist.topo(tF_T, nj(dist(t(cF[sample(nrow(cF), length(cAF)), ])))) }
write.table(data.frame(Iter=1:B, RF_Male=rfM, RF_Female=rfF), "Output_Ox11_OXI/08_Sensitivity_NullDist.txt", sep="\t", row.names=FALSE, quote=FALSE)
pdf("Output_Ox11_OXI/08_Sensitivity_NullDist.pdf", width=10, height=5); par(mfrow=c(1,2)); hist(rfM, breaks=seq(-0.5,12.5,1), col="steelblue", main="Topological Dist (Male)\nTrue vs 1000 Random Sets", xlab="Robinson-Foulds Distance"); hist(rfF, breaks=seq(-0.5,12.5,1), col="darkgreen", main="Topological Dist (Female)\nTrue vs 1000 Random Sets", xlab="Robinson-Foulds Distance"); dev.off()




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

