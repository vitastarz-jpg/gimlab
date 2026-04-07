setwd("G:/내 드라이브/Rrunning/HW050")
rm(list=ls())
library(data.table); library(dplyr); library(matrixStats); library(pheatmap); library(corrplot); library(tidyr); library(ggpubr)
library(maps); library(plotrix)

TCGA_SNV  <- as.data.frame(fread("../TCGA_SNV.txt.gz", header = TRUE, sep = "\t"))
TCGA_Cln  <- as.data.frame(fread("../TCGA_Cln.txt.gz", header = TRUE, sep = "\t"))
LZp5_ALL  <- as.data.frame(fread("../LZALL_logPV5.txt.gz", header = TRUE, sep = "\t"))

KGAS_Af5  <- as.data.frame(fread("../KoGES/Geno/Affy5_ASAS_08840.txt.gz", header = TRUE, sep = "\t")) #1.5m
KGCT_Af6  <- as.data.frame(fread("../KoGES/Geno/Affy6_CITY_03693.txt.gz", header = TRUE, sep = "\t")) #10m
KGTW_Af6  <- as.data.frame(fread("../KoGES/Geno/Affy6_TWIN_01716.txt.gz", header = TRUE, sep = "\t")) #3m
KGAC_Exm  <- as.data.frame(fread("../KoGES/Geno/Exome_AACR_14025.txt.gz", header = TRUE, sep = "\t")) #3m


rownames(KGAS_Af5) <- KGAS_Af5$V1; KGAS_Af5$V1 <- NULL
rownames(KGCT_Af6) <- KGCT_Af6$V1; KGCT_Af6$V1 <- NULL
rownames(KGTW_Af6) <- KGTW_Af6$V1; KGTW_Af6$V1 <- NULL
rownames(KGAC_Exm) <- KGAC_Exm$V1; KGAC_Exm$V1 <- NULL

KGAS_Af5[1:5,1:5]; dim(KGAS_Af5) # 8840,352229
KGCT_Af6[1:5,1:5]; dim(KGCT_Af6) # 3693,627660
KGTW_Af6[1:5,1:5]; dim(KGTW_Af6) # 1716,516611
KGAC_Exm[1:5,1:5]; dim(KGAC_Exm) #14025,77473

GWS_annt <- as.data.frame(fread("../Annotations/31_SNV/gwas_catalog_v1.0-associations_e114_r2025-07-10.tsv.gz", header = TRUE, sep = "\t"))
Exm_annt <- as.data.frame(fread("../Annotations/31_SNV/HumanExome-12v1-1_A_Gene_Annotation.txt.gz", header = TRUE, sep = "\t"))
Af5_annt <- as.data.frame(fread("../Annotations/31_SNV/snpArrayAffy5.txt.gz", header = FALSE, sep = "\t"))[,-6]
Af6_annt <- as.data.frame(fread("../Annotations/31_SNV/snpArrayAffy6.txt.gz", header = FALSE, sep = "\t"))[,-6]
colnames(Exm_annt) <- c("OrID","Chr","Loc","Allele","Transcript","Symbol","Feature","Effect")
colnames(Af5_annt) <- c("No","Chr","Loc","Loc2","OrID","Dir","Allele","rsNo")
colnames(Af6_annt) <- c("No","Chr","Loc","Loc2","OrID","Dir","Allele","rsNo")
Exm_annt$rsNo <- Exm_annt$OrID; Exm_annt$rsNo <- gsub("exm\\-","",Exm_annt$rsNo)

GWS_annt[1:5, 1:5]; dim(GWS_annt)
head(Exm_annt); dim(Exm_annt)
head(Af5_annt); dim(Af5_annt)
head(Af6_annt); dim(Af6_annt)
write.table(GWS_annt[1:5,], file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)

NIBS_rsNo <- c("rs6265","rs4680","rs6277","rs1800497","rs25531","rs6295","rs7301328","rs1042098")
SLPT_rsNo <- c("rs113851554","rs13107325","rs62158211","rs4628869","rs1801260","rs121912617","rs11932595","rs228697","rs3923809","rs2300478","rs10499508","rs3104767","rs1154155","rs2858884","rs5770917")
ExAe_rsNo <- c("rs8192678","rs4253778","rs17602729","rs4340","rs4341","rs4343","rs1815739","rs1042713","rs1042714","rs2010963","rs6552828","rs3213537")
ExAn_rsNo <- c("rs1815739","rs4340","rs1805086","rs35767","rs1800795","rs1049434","rs1800169","rs2228570")
GERD_rsNo <- c("rs6766410","rs11789015","rs9257809","rs2687201","rs4965272","rs11557467")
CNDP_rsNo <- c("rs6566810","rs7489582","rs2346061","rs1187439","rs2294158")

GWS_annt_NIBS <- dplyr::filter(GWS_annt, SNPS %in% NIBS_rsNo); Exm_annt_NIBS <- dplyr::filter(Exm_annt, rsNo %in% NIBS_rsNo); Af5_annt_NIBS <- dplyr::filter(Af5_annt, rsNo %in% NIBS_rsNo); Af6_annt_NIBS <- dplyr::filter(Af6_annt, rsNo %in% NIBS_rsNo)
GWS_annt_SLPT <- dplyr::filter(GWS_annt, SNPS %in% SLPT_rsNo); Exm_annt_SLPT <- dplyr::filter(Exm_annt, rsNo %in% SLPT_rsNo); Af5_annt_SLPT <- dplyr::filter(Af5_annt, rsNo %in% SLPT_rsNo); Af6_annt_SLPT <- dplyr::filter(Af6_annt, rsNo %in% SLPT_rsNo)
GWS_annt_ExAe <- dplyr::filter(GWS_annt, SNPS %in% ExAe_rsNo); Exm_annt_ExAe <- dplyr::filter(Exm_annt, rsNo %in% ExAe_rsNo); Af5_annt_ExAe <- dplyr::filter(Af5_annt, rsNo %in% ExAe_rsNo); Af6_annt_ExAe <- dplyr::filter(Af6_annt, rsNo %in% ExAe_rsNo)
GWS_annt_ExAn <- dplyr::filter(GWS_annt, SNPS %in% ExAn_rsNo); Exm_annt_ExAn <- dplyr::filter(Exm_annt, rsNo %in% ExAn_rsNo); Af5_annt_ExAn <- dplyr::filter(Af5_annt, rsNo %in% ExAn_rsNo); Af6_annt_ExAn <- dplyr::filter(Af6_annt, rsNo %in% ExAn_rsNo)
GWS_annt_GERD <- dplyr::filter(GWS_annt, SNPS %in% GERD_rsNo); Exm_annt_GERD <- dplyr::filter(Exm_annt, rsNo %in% GERD_rsNo); Af5_annt_GERD <- dplyr::filter(Af5_annt, rsNo %in% GERD_rsNo); Af6_annt_GERD <- dplyr::filter(Af6_annt, rsNo %in% GERD_rsNo)
GWS_annt_CNDP <- dplyr::filter(GWS_annt, SNPS %in% CNDP_rsNo); Exm_annt_CNDP <- dplyr::filter(Exm_annt, rsNo %in% CNDP_rsNo); Af5_annt_CNDP <- dplyr::filter(Af5_annt, rsNo %in% CNDP_rsNo); Af6_annt_CNDP <- dplyr::filter(Af6_annt, rsNo %in% CNDP_rsNo)
length(NIBS_rsNo); dim(GWS_annt_NIBS); dim(Exm_annt_NIBS); dim(Af5_annt_NIBS); dim(Af6_annt_NIBS)
length(SLPT_rsNo); dim(GWS_annt_SLPT); dim(Exm_annt_SLPT); dim(Af5_annt_SLPT); dim(Af6_annt_SLPT)
length(ExAe_rsNo); dim(GWS_annt_ExAe); dim(Exm_annt_ExAe); dim(Af5_annt_ExAe); dim(Af6_annt_ExAe)
length(ExAn_rsNo); dim(GWS_annt_ExAn); dim(Exm_annt_ExAn); dim(Af5_annt_ExAn); dim(Af6_annt_ExAn)
length(GERD_rsNo); dim(GWS_annt_GERD); dim(Exm_annt_GERD); dim(Af5_annt_GERD); dim(Af6_annt_GERD)
length(CNDP_rsNo); dim(GWS_annt_CNDP); dim(Exm_annt_CNDP); dim(Af5_annt_CNDP); dim(Af6_annt_CNDP)

KGAS_Af5_NIBS <- dplyr::select(KGAS_Af5, any_of(Af5_annt_NIBS$OrID)); name_map <- setNames(Af5_annt_NIBS$rsNo, Af5_annt_NIBS$OrID); KGAS_Af5_NIBS <- rename_with(KGAS_Af5_NIBS, ~ name_map[.x], any_of(names(name_map)))
KGAS_Af5_SLPT <- dplyr::select(KGAS_Af5, any_of(Af5_annt_SLPT$OrID)); name_map <- setNames(Af5_annt_SLPT$rsNo, Af5_annt_SLPT$OrID); KGAS_Af5_SLPT <- rename_with(KGAS_Af5_SLPT, ~ name_map[.x], any_of(names(name_map)))
KGAS_Af5_ExAe <- dplyr::select(KGAS_Af5, any_of(Af5_annt_ExAe$OrID)); name_map <- setNames(Af5_annt_ExAe$rsNo, Af5_annt_ExAe$OrID); KGAS_Af5_ExAe <- rename_with(KGAS_Af5_ExAe, ~ name_map[.x], any_of(names(name_map)))
KGAS_Af5_ExAn <- dplyr::select(KGAS_Af5, any_of(Af5_annt_ExAn$OrID)); name_map <- setNames(Af5_annt_ExAn$rsNo, Af5_annt_ExAn$OrID); KGAS_Af5_ExAn <- rename_with(KGAS_Af5_ExAn, ~ name_map[.x], any_of(names(name_map)))
KGAS_Af5_GERD <- dplyr::select(KGAS_Af5, any_of(Af5_annt_GERD$OrID)); name_map <- setNames(Af5_annt_GERD$rsNo, Af5_annt_GERD$OrID); KGAS_Af5_GERD <- rename_with(KGAS_Af5_GERD, ~ name_map[.x], any_of(names(name_map)))
KGAS_Af5_CNDP <- dplyr::select(KGAS_Af5, any_of(Af5_annt_CNDP$OrID)); name_map <- setNames(Af5_annt_CNDP$rsNo, Af5_annt_CNDP$OrID); KGAS_Af5_CNDP <- rename_with(KGAS_Af5_CNDP, ~ name_map[.x], any_of(names(name_map)))
head(KGAS_Af5_NIBS); dim(KGAS_Af5_NIBS)
head(KGAS_Af5_SLPT); dim(KGAS_Af5_SLPT)
head(KGAS_Af5_ExAe); dim(KGAS_Af5_ExAe)
head(KGAS_Af5_ExAn); dim(KGAS_Af5_ExAn)
head(KGAS_Af5_GERD); dim(KGAS_Af5_GERD)
head(KGAS_Af5_CNDP); dim(KGAS_Af5_CNDP)

KGCT_Af6_NIBS <- dplyr::select(KGCT_Af6, any_of(Af5_annt_NIBS$OrID)); name_map <- setNames(Af5_annt_NIBS$rsNo, Af5_annt_NIBS$OrID); KGCT_Af6_NIBS <- rename_with(KGCT_Af6_NIBS, ~ name_map[.x], any_of(names(name_map)))
KGCT_Af6_SLPT <- dplyr::select(KGCT_Af6, any_of(Af5_annt_SLPT$OrID)); name_map <- setNames(Af5_annt_SLPT$rsNo, Af5_annt_SLPT$OrID); KGCT_Af6_SLPT <- rename_with(KGCT_Af6_SLPT, ~ name_map[.x], any_of(names(name_map)))
KGCT_Af6_ExAe <- dplyr::select(KGCT_Af6, any_of(Af5_annt_ExAe$OrID)); name_map <- setNames(Af5_annt_ExAe$rsNo, Af5_annt_ExAe$OrID); KGCT_Af6_ExAe <- rename_with(KGCT_Af6_ExAe, ~ name_map[.x], any_of(names(name_map)))
KGCT_Af6_ExAn <- dplyr::select(KGCT_Af6, any_of(Af5_annt_ExAn$OrID)); name_map <- setNames(Af5_annt_ExAn$rsNo, Af5_annt_ExAn$OrID); KGCT_Af6_ExAn <- rename_with(KGCT_Af6_ExAn, ~ name_map[.x], any_of(names(name_map)))
KGCT_Af6_GERD <- dplyr::select(KGCT_Af6, any_of(Af5_annt_GERD$OrID)); name_map <- setNames(Af5_annt_GERD$rsNo, Af5_annt_GERD$OrID); KGCT_Af6_GERD <- rename_with(KGCT_Af6_GERD, ~ name_map[.x], any_of(names(name_map)))
KGCT_Af6_CNDP <- dplyr::select(KGCT_Af6, any_of(Af5_annt_CNDP$OrID)); name_map <- setNames(Af5_annt_CNDP$rsNo, Af5_annt_CNDP$OrID); KGCT_Af6_CNDP <- rename_with(KGCT_Af6_CNDP, ~ name_map[.x], any_of(names(name_map)))
head(KGCT_Af6_NIBS); dim(KGCT_Af6_NIBS)
head(KGCT_Af6_SLPT); dim(KGCT_Af6_SLPT)
head(KGCT_Af6_ExAe); dim(KGCT_Af6_ExAe)
head(KGCT_Af6_ExAn); dim(KGCT_Af6_ExAn)
head(KGCT_Af6_GERD); dim(KGCT_Af6_GERD)
head(KGCT_Af6_CNDP); dim(KGCT_Af6_CNDP)

KGTW_Af6_NIBS <- dplyr::select(KGTW_Af6, any_of(Af5_annt_NIBS$OrID)); name_map <- setNames(Af5_annt_NIBS$rsNo, Af5_annt_NIBS$OrID); KGTW_Af6_NIBS <- rename_with(KGTW_Af6_NIBS, ~ name_map[.x], any_of(names(name_map)))
KGTW_Af6_SLPT <- dplyr::select(KGTW_Af6, any_of(Af5_annt_SLPT$OrID)); name_map <- setNames(Af5_annt_SLPT$rsNo, Af5_annt_SLPT$OrID); KGTW_Af6_SLPT <- rename_with(KGTW_Af6_SLPT, ~ name_map[.x], any_of(names(name_map)))
KGTW_Af6_ExAe <- dplyr::select(KGTW_Af6, any_of(Af5_annt_ExAe$OrID)); name_map <- setNames(Af5_annt_ExAe$rsNo, Af5_annt_ExAe$OrID); KGTW_Af6_ExAe <- rename_with(KGTW_Af6_ExAe, ~ name_map[.x], any_of(names(name_map)))
KGTW_Af6_ExAn <- dplyr::select(KGTW_Af6, any_of(Af5_annt_ExAn$OrID)); name_map <- setNames(Af5_annt_ExAn$rsNo, Af5_annt_ExAn$OrID); KGTW_Af6_ExAn <- rename_with(KGTW_Af6_ExAn, ~ name_map[.x], any_of(names(name_map)))
KGTW_Af6_GERD <- dplyr::select(KGTW_Af6, any_of(Af5_annt_GERD$OrID)); name_map <- setNames(Af5_annt_GERD$rsNo, Af5_annt_GERD$OrID); KGTW_Af6_GERD <- rename_with(KGTW_Af6_GERD, ~ name_map[.x], any_of(names(name_map)))
KGTW_Af6_CNDP <- dplyr::select(KGTW_Af6, any_of(Af5_annt_CNDP$OrID)); name_map <- setNames(Af5_annt_CNDP$rsNo, Af5_annt_CNDP$OrID); KGTW_Af6_CNDP <- rename_with(KGTW_Af6_CNDP, ~ name_map[.x], any_of(names(name_map)))
head(KGTW_Af6_NIBS); dim(KGTW_Af6_NIBS)
head(KGTW_Af6_SLPT); dim(KGTW_Af6_SLPT)
head(KGTW_Af6_ExAe); dim(KGTW_Af6_ExAe)
head(KGTW_Af6_ExAn); dim(KGTW_Af6_ExAn)
head(KGTW_Af6_GERD); dim(KGTW_Af6_GERD)
head(KGTW_Af6_CNDP); dim(KGTW_Af6_CNDP)

KGAC_Exm_NIBS <- dplyr::select(KGAC_Exm, any_of(Exm_annt_NIBS$OrID)); name_map <- setNames(Exm_annt_NIBS$rsNo, Exm_annt_NIBS$OrID); KGAC_Exm_NIBS <- rename_with(KGAC_Exm_NIBS, ~ name_map[.x], any_of(names(name_map)))
KGAC_Exm_SLPT <- dplyr::select(KGAC_Exm, any_of(Exm_annt_SLPT$OrID)); name_map <- setNames(Exm_annt_SLPT$rsNo, Exm_annt_SLPT$OrID); KGAC_Exm_SLPT <- rename_with(KGAC_Exm_SLPT, ~ name_map[.x], any_of(names(name_map)))
KGAC_Exm_ExAe <- dplyr::select(KGAC_Exm, any_of(Exm_annt_ExAe$OrID)); name_map <- setNames(Exm_annt_ExAe$rsNo, Exm_annt_ExAe$OrID); KGAC_Exm_ExAe <- rename_with(KGAC_Exm_ExAe, ~ name_map[.x], any_of(names(name_map)))
KGAC_Exm_ExAn <- dplyr::select(KGAC_Exm, any_of(Exm_annt_ExAn$OrID)); name_map <- setNames(Exm_annt_ExAn$rsNo, Exm_annt_ExAn$OrID); KGAC_Exm_ExAn <- rename_with(KGAC_Exm_ExAn, ~ name_map[.x], any_of(names(name_map)))
KGAC_Exm_GERD <- dplyr::select(KGAC_Exm, any_of(Exm_annt_GERD$OrID)); name_map <- setNames(Exm_annt_GERD$rsNo, Exm_annt_GERD$OrID); KGAC_Exm_GERD <- rename_with(KGAC_Exm_GERD, ~ name_map[.x], any_of(names(name_map)))
KGAC_Exm_CNDP <- dplyr::select(KGAC_Exm, any_of(Exm_annt_CNDP$OrID)); name_map <- setNames(Exm_annt_CNDP$rsNo, Exm_annt_CNDP$OrID); KGAC_Exm_CNDP <- rename_with(KGAC_Exm_CNDP, ~ name_map[.x], any_of(names(name_map)))
head(KGAC_Exm_NIBS); dim(KGAC_Exm_NIBS)
head(KGAC_Exm_SLPT); dim(KGAC_Exm_SLPT)
head(KGAC_Exm_ExAe); dim(KGAC_Exm_ExAe)
head(KGAC_Exm_ExAn); dim(KGAC_Exm_ExAn)
head(KGAC_Exm_GERD); dim(KGAC_Exm_GERD)
head(KGAC_Exm_CNDP); dim(KGAC_Exm_CNDP)


KG_LIST <- unique(c(colnames(KGAS_Af5_NIBS),colnames(KGAS_Af5_SLPT),colnames(KGAS_Af5_ExAe),colnames(KGAS_Af5_ExAn),colnames(KGAS_Af5_GERD),colnames(KGAS_Af5_CNDP),colnames(KGAC_Exm_NIBS),colnames(KGAC_Exm_SLPT),colnames(KGAC_Exm_ExAe),colnames(KGAC_Exm_ExAn),colnames(KGAC_Exm_GERD),colnames(KGAC_Exm_CNDP)))
length(KG_LIST)

dplyr::filter(TCGA_SNV, dbSNP_RS %in% KG_LIST)
dplyr::filter(LZp5_ALL, rsid %in% KG_LIST)

dplyr::filter(LZp5_ALL, rsid == "rs6265")
dplyr::filter(LZp5_ALL, rsid == "rs2346061")

LZp5_ALL_KG <- dplyr::filter(LZp5_ALL, rsid %in% KG_LIST)
head(LZp5_ALL_KG)
sort(table(LZp5_ALL_KG$rsid))
table(LZp5_ALL_KG$CLASS)


LZ_LIST <- c("LZ0000367","LZ0009147","LZ0011313","LZ0011854","LZ0019448","LZ0019554","LZ0023061","LZ0023704","LZ0032488","LZ0035319","LZ0037166","LZ0039839","LZ0041192","LZ0047372","LZ0048511","LZ0049363","LZ0049515","LZ0053259","LZ0054179","LZ0055950","LZ0065547","LZ0066622","LZ0070045","LZ0079035","LZ0081280","LZ0081663","LZ0086093","LZ0091228","LZ0096645","LZ0099623","LZ0111813","LZ0117023","LZ0122558","LZ0126765","LZ0131308","LZ0135281","LZ0135777","LZ0143313","LZ0154813","LZ0160474","LZ0161312","LZ0167711","LZ0170392","LZ0174066","LZ0179566","LZ0181648","LZ0184046","LZ0185011","LZ0191134","LZ0191954","LZ0194301","LZ0194438","LZ0205682","LZ0209873","LZ0213779","LZ0214561","LZ0216769","LZ0223886","LZ0224551","LZ0228917","LZ0230841","LZ0231949","LZ0234415","LZ0234452","LZ0234615","LZ0236887","LZ0238994","LZ0240498","LZ0244114","LZ0245558","LZ0252144","LZ0253931","LZ0255342","LZ0256459","LZ0256899","LZ0256934","LZ0259461","LZ0263813","LZ0264522","LZ0265129","LZ0265359","LZ0275893","LZ0276351","LZ0279362","LZ0283052","LZ0297110","LZ0298296","LZ0310896","LZ0311419","LZ0311578","LZ0312326","LZ0314329","LZ0321945","LZ0324250","LZ0325673","LZ0326230","LZ0329278","LZ0333349","LZ0345564","LZ0346909","LZ0349569","LZ0355203","LZ0361186","LZ0372565","LZ0376775","LZ0379969","LZ0397982","LZ0399768","LZ0401125","LZ0404523","LZ0404740","LZ0412689","LZ0413762","LZ0415536","LZ0417804","LZ0429447","LZ0434301","LZ0447498","LZ0447820","LZ0448067","LZ0450188","LZ0455484","LZ0459595","LZ0465704","LZ0468060","LZ0468791","LZ0473382","LZ0473490","LZ0483380","LZ0486254","LZ0489765","LZ0494698","LZ0499367","LZ0505304","LZ0512418","LZ0514818","LZ0517241","LZ0520569","LZ0521137","LZ0522538","LZ0528432","LZ0531301","LZ0533703","LZ0535896","LZ0544011","LZ0546873","LZ0549442","LZ0550539","LZ0552110","LZ0558637","LZ0575162","LZ0580032","LZ0581723","LZ0582023","LZ0584522","LZ0588053","LZ0589793","LZ0590514","LZ0596908","LZ0602739","LZ0607488","LZ0608812","LZ0614262","LZ0614511","LZ0614790","LZ0616502","LZ0619298","LZ0621701","LZ0621871","LZ0622873","LZ0627041","LZ0633434","LZ0633463","LZ0634596","LZ0635779","LZ0635934","LZ0642020","LZ0644492","LZ0645112","LZ0650563","LZ0651300","LZ0656755","LZ0658048","LZ0661311","LZ0664677","LZ0664785","LZ0669771","LZ0672170","LZ0673735","LZ0675137","LZ0675457","LZ0678270","LZ0679449","LZ0680289","LZ0680956","LZ0685531","LZ0687222","LZ0694421","LZ0707867","LZ0711549","LZ0716558","LZ0718604","LZ0722350","LZ0725414","LZ0726071","LZ0727541","LZ0728324","LZ0732622","LZ0733658","LZ0735967","LZ0736573","LZ0736976","LZ0738451","LZ0750543","LZ0751427","LZ0752157","LZ0754657","LZ0757370","LZ0761488","LZ0761777","LZ0765194","LZ0774274","LZ0776774","LZ0780570","LZ0781644","LZ0782059","LZ0789750","LZ0793872","LZ0794701","LZ0794777","LZ0796235","LZ0797689","LZ0803479","LZ0809171","LZ0810010","LZ0811240","LZ0812360","LZ0816934","LZ0819313","LZ0819478","LZ0826670","LZ0832235","LZ0832664","LZ0833892","LZ0834895","LZ0835235","LZ0841221","LZ0842685","LZ0844282","LZ0845276","LZ0845519","LZ0845782","LZ0847758","LZ0852108","LZ0862997","LZ0867416","LZ0876767","LZ0878271","LZ0879529","LZ0880496","LZ0880842","LZ0881707","LZ0890462","LZ0892646","LZ0894486","LZ0896371","LZ0898954","LZ0904725","LZ0914487","LZ0918951","LZ0919709","LZ0924322","LZ0929241","LZ0931496","LZ0934401","LZ0941195","LZ0945975","LZ0947479","LZ0953714","LZ0957567","LZ0957667","LZ0958269","LZ0965158","LZ0968359","LZ0977728","LZ0979922","LZ0980498","LZ0981834","LZ0986607","LZ0988545","LZ0990296","LZ0991276","LZ0991692","LZ0993531","LZ0995061","LZ0998663")
KoGES_COPD_rsNo <- c("rs2609264","rs7671167","rs8192575","rs2070600","rs2239688","rs16951883")

LZp5_ALLx <- dplyr::filter(LZp5_ALL, rsid %in% KoGES_COPD_rsNo); LZp5_ALLx <- LZp5_ALLx[order(LZp5_ALLx$neg_log_pvalue, decreasing = TRUE),]
TCGA_SNVx <- dplyr::filter(TCGA_SNV, dbSNP_RS %in% KoGES_COPD_rsNo)

head(LZp5_ALLx); dim(LZp5_ALLx); sort(table(LZp5_ALLx$CLASS)); sort(table(LZp5_ALLx$rsid))
head(TCGA_SNVx); dim(TCGA_SNVx); table(TCGA_SNVx$Project)
write.table(LZp5_ALLx, file="Output_COPD_MAP/LZp5_ALLx.txt",sep="\t",quote=FALSE,row.names=FALSE)

#setwd("//10.112.90.163/Omics/25_1000_Genomes_Project/")
#rsNo_ALL  <- as.data.frame(fread("rsNo_ALL.txt.gz", header = TRUE, sep = "\t"))
#head(rsNo_ALL); dim(rsNo_ALL)
#
#rs01 <- dplyr::filter(rsNo_ALL, V1 =="rs2609264" )
#rs02 <- dplyr::filter(rsNo_ALL, V1 =="rs7671167" )
#rs03 <- dplyr::filter(rsNo_ALL, V1 =="rs8192575" )
#rs04 <- dplyr::filter(rsNo_ALL, V1 =="rs2070600" )
#rs05 <- dplyr::filter(rsNo_ALL, V1 =="rs2239688" )
#rs06 <- dplyr::filter(rsNo_ALL, V1 =="rs16951883")
#rsALL <- rbind(rs01,rs02,rs03,rs04,rs05,rs06)
#rsALL <- rsALL[order(rsALL$Class),]
#rsALL; dim(rsALL)
#unique(rsALL$Class)
#rm(rsNo_ALL)
#
#CHR04_af  <- as.data.frame(fread("CHR04_af_VAR.txt.gz", header = TRUE, sep = "\t"))
#CHR06_ab  <- as.data.frame(fread("CHR06_ab_VAR.txt.gz", header = TRUE, sep = "\t"))
#CHR16_aa  <- as.data.frame(fread("CHR16_aa_VAR.txt.gz", header = TRUE, sep = "\t"))
#
#KG_rs2609264 <- as.data.frame(t(as.matrix(dplyr::select(dplyr::filter(CHR04_af, ID == "rs2609264"), -ID)))); names(KG_rs2609264) <- "rs2609264"
#KG_rs7671167 <- as.data.frame(t(as.matrix(dplyr::select(dplyr::filter(CHR04_af, ID == "rs7671167"), -ID)))); names(KG_rs7671167) <- "rs7671167"
#rm(CHR04_af)
#
#KG_rs8192575 <- as.data.frame(t(as.matrix(dplyr::select(dplyr::filter(CHR06_ab, ID == "rs8192575"), -ID)))); names(KG_rs8192575) <- "rs8192575"
#KG_rs2070600 <- as.data.frame(t(as.matrix(dplyr::select(dplyr::filter(CHR06_ab, ID == "rs2070600"), -ID)))); names(KG_rs2070600) <- "rs2070600"
#KG_rs2239688 <- as.data.frame(t(as.matrix(dplyr::select(dplyr::filter(CHR06_ab, ID == "rs2239688"), -ID)))); names(KG_rs2239688) <- "rs2239688"
#rm(CHR06_ab)
#
#KG_rs16951883 <- as.data.frame(t(as.matrix(dplyr::select(dplyr::filter(CHR16_aa, ID == "rs16951883"), -ID)))); names(KG_rs16951883) <- "rs16951883"
#rm(CHR16_aa)
#
#KG <- cbind(KG_rs2609264,KG_rs7671167,KG_rs8192575,KG_rs2070600,KG_rs2239688,KG_rs16951883)
#head(KG); dim(KG)
#setwd("G:/내 드라이브/Rrunning/KW25/")
#write.table(KG, "Output_COPD_MAP/KG.txt",sep="\t", quote=FALSE,col.names=NA)



KG <- read.csv("Output_COPD_MAP/KG.txt",sep="\t", row.names = 1)
for(i in 1:ncol(KG)){  KG[,i] <- ifelse(KG[,i]=="0|0",0,1)}
head(KG); dim(KG) #2504,6
table(KG[,1]); table(KG[,2]); table(KG[,3]); table(KG[,4]); table(KG[,5]); table(KG[,6])
KG$ID <- rownames(KG)

SI <- read.csv("1GP_MAP/igsr_samples.tsv",sep="\t")
names(SI)[1] <- "ID"
SI <- dplyr::filter(SI, ID %in% rownames(KG))
head(SI); dim(SI)
SI <- dplyr::select(SI, ID, Sex, Population.code, Superpopulation.code)
SI$Population.code <- gsub("IBS,MSL", "IBS", SI$Population.code)
table(SI$Sex)
table(SI$Population.code)
length(unique(SI$Population.code))

KG <- merge(x=KG,y=SI,by="ID")
head(KG); dim(KG) #2504,10
KGP <- data.frame(table(KG$Population.code))
names(KGP)[1] <- "POP"



KGdf01 <- data.frame(table(KG[,c(2,(ncol(KG)-1))]))
KGdf02 <- data.frame(table(KG[,c(3,(ncol(KG)-1))]))
KGdf03 <- data.frame(table(KG[,c(4,(ncol(KG)-1))]))
KGdf04 <- data.frame(table(KG[,c(5,(ncol(KG)-1))]))
KGdf05 <- data.frame(table(KG[,c(6,(ncol(KG)-1))]))
KGdf06 <- data.frame(table(KG[,c(7,(ncol(KG)-1))]))

KGdf01x <- cbind(KGdf01[KGdf01[,1] == 0, 3], KGdf01[KGdf01[,1] == 1, 3])
KGdf02x <- cbind(KGdf02[KGdf02[,1] == 0, 3], KGdf02[KGdf02[,1] == 1, 3])
KGdf03x <- cbind(KGdf03[KGdf03[,1] == 0, 3], KGdf03[KGdf03[,1] == 1, 3])
KGdf04x <- cbind(KGdf04[KGdf04[,1] == 0, 3], KGdf04[KGdf04[,1] == 1, 3])
KGdf05x <- cbind(KGdf05[KGdf05[,1] == 0, 3], KGdf05[KGdf05[,1] == 1, 3])
KGdf06x <- cbind(KGdf06[KGdf06[,1] == 0, 3], KGdf06[KGdf06[,1] == 1, 3])

rownames(KGdf01x) <- unique(KGdf01$Population.code); colnames(KGdf01x) <- c("VarX", "VarO")
rownames(KGdf02x) <- unique(KGdf02$Population.code); colnames(KGdf02x) <- c("VarX", "VarO")
rownames(KGdf03x) <- unique(KGdf03$Population.code); colnames(KGdf03x) <- c("VarX", "VarO")
rownames(KGdf04x) <- unique(KGdf04$Population.code); colnames(KGdf04x) <- c("VarX", "VarO")
rownames(KGdf05x) <- unique(KGdf05$Population.code); colnames(KGdf05x) <- c("VarX", "VarO")
rownames(KGdf06x) <- unique(KGdf06$Population.code); colnames(KGdf06x) <- c("VarX", "VarO")

KGdf01y <- KGdf01x; KGdf01y[,1] <- KGdf01x[,1]/rowSums(KGdf01x); KGdf01y[,2] <- KGdf01x[,2]/rowSums(KGdf01x)
KGdf02y <- KGdf02x; KGdf02y[,1] <- KGdf02x[,1]/rowSums(KGdf02x); KGdf02y[,2] <- KGdf02x[,2]/rowSums(KGdf02x)
KGdf03y <- KGdf03x; KGdf03y[,1] <- KGdf03x[,1]/rowSums(KGdf03x); KGdf03y[,2] <- KGdf03x[,2]/rowSums(KGdf03x)
KGdf04y <- KGdf04x; KGdf04y[,1] <- KGdf04x[,1]/rowSums(KGdf04x); KGdf04y[,2] <- KGdf04x[,2]/rowSums(KGdf04x)
KGdf05y <- KGdf05x; KGdf05y[,1] <- KGdf05x[,1]/rowSums(KGdf05x); KGdf05y[,2] <- KGdf05x[,2]/rowSums(KGdf05x)
KGdf06y <- KGdf06x; KGdf06y[,1] <- KGdf06x[,1]/rowSums(KGdf06x); KGdf06y[,2] <- KGdf06x[,2]/rowSums(KGdf06x)

KGdf01y <- as.data.frame(KGdf01y); rowSums(KGdf01y)
KGdf02y <- as.data.frame(KGdf02y); rowSums(KGdf02y)
KGdf03y <- as.data.frame(KGdf03y); rowSums(KGdf03y)
KGdf04y <- as.data.frame(KGdf04y); rowSums(KGdf04y)
KGdf05y <- as.data.frame(KGdf05y); rowSums(KGdf05y)
KGdf06y <- as.data.frame(KGdf06y); rowSums(KGdf06y)

GPS <- read.csv("1GP_MAP/sample_info_superpop_coord.tsv",sep="\t")
GPS <- dplyr::select(GPS, POP, lon, lat)
GPS <- merge(x=GPS, y=KGP, by="POP")
head(GPS); dim(GPS) #26,4

pdf("Output_COPD_MAP/map01_rs2609264x.pdf", width=32, height=16); map("world"); title(main = "rs2609264  Allele Frequency Map", cex.main = 3) ; for (i in 1:nrow(GPS)){  floating.pie(GPS$lon[i],GPS$lat[i], c(KGdf01y$VarX[i],KGdf01y$VarO[i]),radius=GPS$Freq[i]/30, col=c("#0571B0","#CA0020"))}; dev.off()
pdf("Output_COPD_MAP/map02_rs7671167x.pdf", width=32, height=16); map("world"); title(main = "rs7671167  Allele Frequency Map", cex.main = 3) ; for (i in 1:nrow(GPS)){  floating.pie(GPS$lon[i],GPS$lat[i], c(KGdf02y$VarX[i],KGdf02y$VarO[i]),radius=GPS$Freq[i]/30, col=c("#0571B0","#CA0020"))}; dev.off()
pdf("Output_COPD_MAP/map03_rs8192575x.pdf", width=32, height=16); map("world"); title(main = "rs8192575  Allele Frequency Map", cex.main = 3) ; for (i in 1:nrow(GPS)){  floating.pie(GPS$lon[i],GPS$lat[i], c(KGdf03y$VarX[i],KGdf03y$VarO[i]),radius=GPS$Freq[i]/30, col=c("#0571B0","#CA0020"))}; dev.off()
pdf("Output_COPD_MAP/map04_rs2070600x.pdf", width=32, height=16); map("world"); title(main = "rs2070600  Allele Frequency Map", cex.main = 3) ; for (i in 1:nrow(GPS)){  floating.pie(GPS$lon[i],GPS$lat[i], c(KGdf04y$VarX[i],KGdf04y$VarO[i]),radius=GPS$Freq[i]/30, col=c("#0571B0","#CA0020"))}; dev.off()
pdf("Output_COPD_MAP/map05_rs2239688x.pdf", width=32, height=16); map("world"); title(main = "rs2239688  Allele Frequency Map", cex.main = 3) ; for (i in 1:nrow(GPS)){  floating.pie(GPS$lon[i],GPS$lat[i], c(KGdf05y$VarX[i],KGdf05y$VarO[i]),radius=GPS$Freq[i]/30, col=c("#0571B0","#CA0020"))}; dev.off()
pdf("Output_COPD_MAP/map06_rs16951883.pdf", width=32, height=16); map("world"); title(main = "rs16951883 Allele Frequency Map", cex.main = 3) ; for (i in 1:nrow(GPS)){  floating.pie(GPS$lon[i],GPS$lat[i], c(KGdf06y$VarX[i],KGdf06y$VarO[i]),radius=GPS$Freq[i]/30, col=c("#0571B0","#CA0020"))}; dev.off()

KGdf01y$POP <- rownames(KGdf01y); KGdf01y$rsNo <- "rs2609264" 
KGdf02y$POP <- rownames(KGdf02y); KGdf02y$rsNo <- "rs7671167" 
KGdf03y$POP <- rownames(KGdf03y); KGdf03y$rsNo <- "rs8192575" 
KGdf04y$POP <- rownames(KGdf04y); KGdf04y$rsNo <- "rs2070600" 
KGdf05y$POP <- rownames(KGdf05y); KGdf05y$rsNo <- "rs2239688" 
KGdf06y$POP <- rownames(KGdf06y); KGdf06y$rsNo <- "rs16951883"
POPratio <- rbind(KGdf01y,KGdf02y,KGdf03y,KGdf04y,KGdf05y,KGdf06y)[,c(4,3,1:2)]
POPratio$VarX <- round(POPratio$VarX,2); POPratio$VarO <- round(POPratio$VarO,2)
write.table(POPratio, file="Output_COPD_MAP/POPratio.txt",sep="\t",quote=FALSE,row.names=FALSE)




