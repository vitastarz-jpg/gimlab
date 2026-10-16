################################################################################################
#Module 00. Clear memory and Load library
################################################################################################
Sys.time()
rm(list=ls())
library(dplyr); library(data.table); library(pheatmap); library(readxl)
library(ggvenn); library(gridExtra)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
Rep_NA <- function(x) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE); Q3 <- quantile(x, 0.75, na.rm = TRUE); IQR_val <- Q3 - Q1;
  lower <- Q1 - 1.5 * IQR_val;   upper <- Q3 + 1.5 * IQR_val; ifelse(x < lower | x > upper, NA, x)}
custom_col1 <- colorRampPalette(c("#60C5F1", "#FFFCCC", "#F47378"))(500)
custom_col2 <- colorRampPalette(c("#FFFCCC", "#D01C8B"))(6)
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


################################################################################################
#Module 01. Data loading
################################################################################################
Sys.time()
source("G:/내 드라이브/Rrunning/TCGA/Loading_TCGA_Met.R")
Met450pat <- ls(pattern="_450m"); Met450patdf <- Met450pat[sapply(mget(Met450pat), is.data.frame)]
COMMON_CpGs_TCGA <- Reduce(intersect, lapply(mget(Met450patdf), rownames))
Sys.time()


setwd("G:/내 드라이브/Rrunning/GEO/HW014/Primate")
#GSE103286
GSE103286m <- as.data.frame(fread("GSE103286_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=69))
rownames(GSE103286m) <- GSE103286m$ID_REF; GSE103286m <- GSE103286m[,-1]
GSE103286m <- GSE103286m[complete.cases(GSE103286m),]
GSE103286m[1:3,1:5]; dim(GSE103286m) #191954,112

GSE103286i <- as.data.frame(fread("GSE103286_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=28)); names(GSE103286i)[1] <- "ID"
GSE103286i <- dplyr::filter(GSE103286i, grepl("characteristics", ID))
GSE103286i <- as.data.frame(t(GSE103286i[,-1]))[,-5]
GSE103286i$V1 <- gsub("Sex: male",  "M",GSE103286i$V1)
GSE103286i$V1 <- gsub("Sex: female",  "F",GSE103286i$V1)
GSE103286i$V2 <- gsub("age \\(years\\): ",  "",GSE103286i$V2)
GSE103286i$V3 <- gsub("steady state weight \\(kg\\): ",  "",GSE103286i$V3)
GSE103286i$V4 <- gsub("tissue: ",  "",GSE103286i$V4)
GSE103286i$V6 <- gsub("disease state: ",  "",GSE103286i$V6)
colnames(GSE103286i) <- c("Sex","Age","Weight","Tissue","Disease")
head(GSE103286i); dim(GSE103286i) #69,9
table(colnames(GSE103286m)==rownames(GSE103286i))

#GSE103271
GSE103271m <- as.data.frame(fread("GSE103271_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=96))
rownames(GSE103271m) <- GSE103271m$ID_REF; GSE103271m <- GSE103271m[,-1]
GSE103271m <- GSE103271m[complete.cases(GSE103271m),]
GSE103271m[1:3,1:5]; dim(GSE103271m) #190898,10

GSE103271i <- as.data.frame(fread("GSE103271_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=29)); names(GSE103271i)[1] <- "ID"
GSE103271i <- dplyr::filter(GSE103271i, grepl("characteristics", ID))
GSE103271i <- as.data.frame(t(GSE103271i[,-1]))
GSE103271i$V1  <- gsub("Sex: "                                  ,  "",GSE103271i$V1 )
GSE103271i$V2  <- gsub("age \\(years\\): "                      ,  "",GSE103271i$V2 )
GSE103271i$V3  <- gsub("tissue: "                               ,  "",GSE103271i$V3 )
GSE103271i$V4  <- gsub("femur length \\(mm\\): "                ,  "",GSE103271i$V4 )
GSE103271i$V5  <- gsub("bicondylar femur length \\(mm\\): "     ,  "",GSE103271i$V5 )
GSE103271i$V6  <- gsub("maximum femur length \\(mm\\): "        ,  "",GSE103271i$V6 )
GSE103271i$V7  <- gsub("superior shaft width \\(mm\\): "        ,  "",GSE103271i$V7 )
GSE103271i$V8  <- gsub("superior shaft depth \\(mm\\): "        ,  "",GSE103271i$V8 )
GSE103271i$V9  <- gsub("middle shaft width \\(mm\\): "          ,  "",GSE103271i$V9 )
GSE103271i$V10 <- gsub("middle shaft depth \\(mm\\): "          ,  "",GSE103271i$V10)
GSE103271i$V11 <- gsub("inferior shaft width \\(mm\\): "        ,  "",GSE103271i$V11)
GSE103271i$V12 <- gsub("inferior shaft depth \\(mm\\): "        ,  "",GSE103271i$V12)
GSE103271i$V13 <- gsub("head height \\(mm\\): "                 ,  "",GSE103271i$V13)
GSE103271i$V14 <- gsub("head length \\(mm\\): "                 ,  "",GSE103271i$V14)
GSE103271i$V15 <- gsub("head width \\(mm\\): "                  ,  "",GSE103271i$V15)
GSE103271i$V16 <- gsub("anatomical neck length \\(mm\\): "      ,  "",GSE103271i$V16)
GSE103271i$V17 <- gsub("anatomical neck height \\(mm\\): "      ,  "",GSE103271i$V17)
GSE103271i$V18 <- gsub("anatomical neck depth \\(mm\\): "       ,  "",GSE103271i$V18)
GSE103271i$V19 <- gsub("biomechanical neck length \\(mm\\): "   ,  "",GSE103271i$V19)
GSE103271i$V20 <- gsub("proximal width \\(mm\\): "              ,  "",GSE103271i$V20)
GSE103271i$V21 <- gsub("lesser trochanter to head \\(mm\\): "   ,  "",GSE103271i$V21)
GSE103271i$V22 <- gsub("lesser trochanter to neck \\(mm\\): "   ,  "",GSE103271i$V22)
GSE103271i$V23 <- gsub("lesser to greater trochanter \\(mm\\): ",  "",GSE103271i$V23)
GSE103271i$V24 <- gsub("medial condyle height \\(mm\\): "       ,  "",GSE103271i$V24)
GSE103271i$V25 <- gsub("medial condyle depth \\(mm\\): "        ,  "",GSE103271i$V25)
GSE103271i$V26 <- gsub("medial condyle width \\(mm\\): "        ,  "",GSE103271i$V26)
GSE103271i$V27 <- gsub("lateral condyle height \\(mm\\): "      ,  "",GSE103271i$V27)
GSE103271i$V28 <- gsub("lateral condyle depth \\(mm\\): "       ,  "",GSE103271i$V28)
GSE103271i$V29 <- gsub("lateral condyle width \\(mm\\): "       ,  "",GSE103271i$V29)
GSE103271i$V30 <- gsub("intercondylar notch width \\(mm\\): "   ,  "",GSE103271i$V30)
GSE103271i$V31 <- gsub("intercondylar notch depth \\(mm\\): "   ,  "",GSE103271i$V31)
GSE103271i$V32 <- gsub("bicondylar width \\(mm\\): "            ,  "",GSE103271i$V32)
colnames(GSE103271i) <- c("Sex","Age","Tissue","FL","BFL","MFL","SSWd","SSDp","MSWd","MSDp","ISWd","ISDp","HH","HL","HW","ANL","ANH","AND","BNL","PW","LTH","LTN","LGT","MCH","MCDp","MCW","LCH","LCDp","LCW","INW","IND","BCW")
head(GSE103271i); dim(GSE103271i) #69,9
table(colnames(GSE103271m)==rownames(GSE103271i))

#GSE103279
GSE103279m <- as.data.frame(fread("GSE103279_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=99))
rownames(GSE103279m) <- GSE103279m$ID_REF; GSE103279m <- GSE103279m[,-1]
GSE103279m <- GSE103279m[complete.cases(GSE103279m),]
GSE103279m[1:3,1:5]; dim(GSE103279m) #189858,28

GSE103279i <- as.data.frame(fread("GSE103279_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=29)); names(GSE103279i)[1] <- "ID"
GSE103279i <- dplyr::filter(GSE103279i, grepl("characteristics", ID))
GSE103279i <- as.data.frame(t(GSE103279i[,-1]))
GSE103279i$V1  <- gsub("Sex: "                                  ,  "",GSE103279i$V1 )
GSE103279i$V2  <- gsub("age \\(years\\): "                      ,  "",GSE103279i$V2 )
GSE103279i$V3  <- gsub("steady state weight \\(kg\\): "         ,  "",GSE103279i$V3 )
GSE103279i$V4  <- gsub("tissue: "                               ,  "",GSE103279i$V4 )
GSE103279i$V5  <- gsub("age cohort: "                           ,  "",GSE103279i$V5 )
GSE103279i$V6  <- gsub("disease state: "                        ,  "",GSE103279i$V6 )
GSE103279i$V7  <- gsub("femur length \\(mm\\): "                ,  "",GSE103279i$V7 )
GSE103279i$V8  <- gsub("bicondylar femur length \\(mm\\): "     ,  "",GSE103279i$V8 )
GSE103279i$V9  <- gsub("maximum femur length \\(mm\\): "        ,  "",GSE103279i$V9 )
GSE103279i$V10 <- gsub("superior shaft width \\(mm\\): "        ,  "",GSE103279i$V10)
GSE103279i$V11 <- gsub("superior shaft depth \\(mm\\): "        ,  "",GSE103279i$V11)
GSE103279i$V12 <- gsub("middle shaft width \\(mm\\): "          ,  "",GSE103279i$V12)
GSE103279i$V13 <- gsub("middle shaft depth \\(mm\\): "          ,  "",GSE103279i$V13)
GSE103279i$V14 <- gsub("inferior shaft width \\(mm\\): "        ,  "",GSE103279i$V14)
GSE103279i$V15 <- gsub("inferior shaft depth \\(mm\\): "        ,  "",GSE103279i$V15)
GSE103279i$V16 <- gsub("head height \\(mm\\): "                 ,  "",GSE103279i$V16)
GSE103279i$V17 <- gsub("head length \\(mm\\): "                 ,  "",GSE103279i$V17)
GSE103279i$V18 <- gsub("head width \\(mm\\): "                  ,  "",GSE103279i$V18)
GSE103279i$V19 <- gsub("anatomical neck length \\(mm\\): "      ,  "",GSE103279i$V19)
GSE103279i$V20 <- gsub("anatomical neck height \\(mm\\): "      ,  "",GSE103279i$V20)
GSE103279i$V21 <- gsub("anatomical neck depth \\(mm\\): "       ,  "",GSE103279i$V21)
GSE103279i$V22 <- gsub("biomechanical neck length \\(mm\\): "   ,  "",GSE103279i$V22)
GSE103279i$V23 <- gsub("proximal width \\(mm\\): "              ,  "",GSE103279i$V23)
GSE103279i$V24 <- gsub("lesser trochanter to head \\(mm\\): "   ,  "",GSE103279i$V24)
GSE103279i$V25 <- gsub("lesser trochanter to neck \\(mm\\): "   ,  "",GSE103279i$V25)
GSE103279i$V26 <- gsub("lesser to greater trochanter \\(mm\\): ",  "",GSE103279i$V26)
GSE103279i$V27 <- gsub("medial condyle height \\(mm\\): "       ,  "",GSE103279i$V27)
GSE103279i$V28 <- gsub("medial condyle depth \\(mm\\): "        ,  "",GSE103279i$V28)
GSE103279i$V29 <- gsub("medial condyle width \\(mm\\): "        ,  "",GSE103279i$V29)
GSE103279i$V30 <- gsub("lateral condyle height \\(mm\\): "      ,  "",GSE103279i$V30)
GSE103279i$V31 <- gsub("lateral condyle depth \\(mm\\): "       ,  "",GSE103279i$V31)
GSE103279i$V32 <- gsub("lateral condyle width \\(mm\\): "       ,  "",GSE103279i$V32)
GSE103279i$V33 <- gsub("intercondylar notch width \\(mm\\): "   ,  "",GSE103279i$V33)
GSE103279i$V34 <- gsub("intercondylar notch depth \\(mm\\): "   ,  "",GSE103279i$V34)
GSE103279i$V35 <- gsub("bicondylar width \\(mm\\): "            ,  "",GSE103279i$V35)
colnames(GSE103279i) <- c("Sex","Age","SSW","Tissue","AgeC","Disease","FL","BFL","MFL","SSWd","SSDp","MSWd","MSDp","ISWd","ISDp","HH","HL","HW","ANL","ANH","AND","BNL","PW","LTH","LTN","LGT","MCH","MCDp","MCW","LCH","LCDp","LCW","INW","IND","BCW")
head(GSE103279i); dim(GSE103279i) #69,9
table(colnames(GSE103279m)==rownames(GSE103279i))

#GSE103280
GSE103280m <- as.data.frame(fread("GSE103280_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=96))
rownames(GSE103280m) <- GSE103280m$ID_REF; GSE103280m <- GSE103280m[,-1]
GSE103280m <- GSE103280m[complete.cases(GSE103280m),]
GSE103280m[1:3,1:5]; dim(GSE103280m) #191639,10

GSE103280i <- as.data.frame(fread("GSE103280_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=29)); names(GSE103280i)[1] <- "ID"
GSE103280i <- dplyr::filter(GSE103280i, grepl("characteristics", ID))
GSE103280i <- as.data.frame(t(GSE103280i[,-1]))
GSE103280i$V1  <- gsub("Sex: "                                  ,  "",GSE103280i$V1 )
GSE103280i$V2  <- gsub("age \\(years\\): "                      ,  "",GSE103280i$V2 )
GSE103280i$V3  <- gsub("tissue: "                               ,  "",GSE103280i$V3 )
GSE103280i$V4  <- gsub("femur length \\(mm\\): "                ,  "",GSE103280i$V4 )
GSE103280i$V5  <- gsub("bicondylar femur length \\(mm\\): "     ,  "",GSE103280i$V5 )
GSE103280i$V6  <- gsub("maximum femur length \\(mm\\): "        ,  "",GSE103280i$V6 )
GSE103280i$V7  <- gsub("superior shaft width \\(mm\\): "        ,  "",GSE103280i$V7 )
GSE103280i$V8  <- gsub("superior shaft depth \\(mm\\): "        ,  "",GSE103280i$V8 )
GSE103280i$V9  <- gsub("middle shaft width \\(mm\\): "          ,  "",GSE103280i$V9 )
GSE103280i$V10 <- gsub("middle shaft depth \\(mm\\): "          ,  "",GSE103280i$V10)
GSE103280i$V11 <- gsub("inferior shaft width \\(mm\\): "        ,  "",GSE103280i$V11)
GSE103280i$V12 <- gsub("inferior shaft depth \\(mm\\): "        ,  "",GSE103280i$V12)
GSE103280i$V13 <- gsub("head height \\(mm\\): "                 ,  "",GSE103280i$V13)
GSE103280i$V14 <- gsub("head length \\(mm\\): "                 ,  "",GSE103280i$V14)
GSE103280i$V15 <- gsub("head width \\(mm\\): "                  ,  "",GSE103280i$V15)
GSE103280i$V16 <- gsub("anatomical neck length \\(mm\\): "      ,  "",GSE103280i$V16)
GSE103280i$V17 <- gsub("anatomical neck height \\(mm\\): "      ,  "",GSE103280i$V17)
GSE103280i$V18 <- gsub("anatomical neck depth \\(mm\\): "       ,  "",GSE103280i$V18)
GSE103280i$V19 <- gsub("biomechanical neck length \\(mm\\): "   ,  "",GSE103280i$V19)
GSE103280i$V20 <- gsub("proximal width \\(mm\\): "              ,  "",GSE103280i$V20)
GSE103280i$V21 <- gsub("lesser trochanter to head \\(mm\\): "   ,  "",GSE103280i$V21)
GSE103280i$V22 <- gsub("lesser trochanter to neck \\(mm\\): "   ,  "",GSE103280i$V22)
GSE103280i$V23 <- gsub("lesser to greater trochanter \\(mm\\): ",  "",GSE103280i$V23)
GSE103280i$V24 <- gsub("medial condyle height \\(mm\\): "       ,  "",GSE103280i$V24)
GSE103280i$V25 <- gsub("medial condyle depth \\(mm\\): "        ,  "",GSE103280i$V25)
GSE103280i$V26 <- gsub("medial condyle width \\(mm\\): "        ,  "",GSE103280i$V26)
GSE103280i$V27 <- gsub("lateral condyle height \\(mm\\): "      ,  "",GSE103280i$V27)
GSE103280i$V28 <- gsub("lateral condyle depth \\(mm\\): "       ,  "",GSE103280i$V28)
GSE103280i$V29 <- gsub("lateral condyle width \\(mm\\): "       ,  "",GSE103280i$V29)
GSE103280i$V30 <- gsub("intercondylar notch width \\(mm\\): "   ,  "",GSE103280i$V30)
GSE103280i$V31 <- gsub("intercondylar notch depth \\(mm\\): "   ,  "",GSE103280i$V31)
GSE103280i$V32 <- gsub("bicondylar width \\(mm\\): "            ,  "",GSE103280i$V32)
colnames(GSE103280i) <- c("Sex","Age","Tissue","FL","BFL","MFL","SSWd","SSDp","MSWd","MSDp","ISWd","ISDp","HH","HL","HW","ANL","ANH","AND","BNL","PW","LTH","LTN","LGT","MCH","MCDp","MCW","LCH","LCDp","LCW","INW","IND","BCW")
head(GSE103280i); dim(GSE103280i) #69,9
table(colnames(GSE103280m)==rownames(GSE103280i))

#GSE103328
GSE103328m <- as.data.frame(fread("GSE103328_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=96))
rownames(GSE103328m) <- GSE103328m$ID_REF; GSE103328m <- GSE103328m[,-1]
GSE103328m <- GSE103328m[complete.cases(GSE103328m),]
GSE103328m[1:3,1:5]; dim(GSE103328m) #68709,545

GSE103328i <- as.data.frame(fread("GSE103328_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=29)); names(GSE103328i)[1] <- "ID"
GSE103328i <- dplyr::filter(GSE103328i, grepl("characteristics", ID))
GSE103328i <- as.data.frame(t(GSE103328i[,-1]))
GSE103328i$V1  <- gsub("Sex: "                                  ,  "",GSE103328i$V1 )
GSE103328i$V2  <- gsub("age \\(years\\): "                      ,  "",GSE103328i$V2 )
GSE103328i$V3  <- gsub("tissue: "                               ,  "",GSE103328i$V3 )
GSE103328i$V4  <- gsub("femur length \\(mm\\): "                ,  "",GSE103328i$V4 )
GSE103328i$V5  <- gsub("bicondylar femur length \\(mm\\): "     ,  "",GSE103328i$V5 )
GSE103328i$V6  <- gsub("maximum femur length \\(mm\\): "        ,  "",GSE103328i$V6 )
GSE103328i$V7  <- gsub("superior shaft width \\(mm\\): "        ,  "",GSE103328i$V7 )
GSE103328i$V8  <- gsub("superior shaft depth \\(mm\\): "        ,  "",GSE103328i$V8 )
GSE103328i$V9  <- gsub("middle shaft width \\(mm\\): "          ,  "",GSE103328i$V9 )
GSE103328i$V10 <- gsub("middle shaft depth \\(mm\\): "          ,  "",GSE103328i$V10)
GSE103328i$V11 <- gsub("inferior shaft width \\(mm\\): "        ,  "",GSE103328i$V11)
GSE103328i$V12 <- gsub("inferior shaft depth \\(mm\\): "        ,  "",GSE103328i$V12)
GSE103328i$V13 <- gsub("head height \\(mm\\): "                 ,  "",GSE103328i$V13)
GSE103328i$V14 <- gsub("head length \\(mm\\): "                 ,  "",GSE103328i$V14)
GSE103328i$V15 <- gsub("head width \\(mm\\): "                  ,  "",GSE103328i$V15)
GSE103328i$V16 <- gsub("anatomical neck length \\(mm\\): "      ,  "",GSE103328i$V16)
GSE103328i$V17 <- gsub("anatomical neck height \\(mm\\): "      ,  "",GSE103328i$V17)
GSE103328i$V18 <- gsub("anatomical neck depth \\(mm\\): "       ,  "",GSE103328i$V18)
GSE103328i$V19 <- gsub("biomechanical neck length \\(mm\\): "   ,  "",GSE103328i$V19)
GSE103328i$V20 <- gsub("proximal width \\(mm\\): "              ,  "",GSE103328i$V20)
GSE103328i$V21 <- gsub("lesser trochanter to head \\(mm\\): "   ,  "",GSE103328i$V21)
GSE103328i$V22 <- gsub("lesser trochanter to neck \\(mm\\): "   ,  "",GSE103328i$V22)
GSE103328i$V23 <- gsub("lesser to greater trochanter \\(mm\\): ",  "",GSE103328i$V23)
GSE103328i$V24 <- gsub("medial condyle height \\(mm\\): "       ,  "",GSE103328i$V24)
GSE103328i$V25 <- gsub("medial condyle depth \\(mm\\): "        ,  "",GSE103328i$V25)
GSE103328i$V26 <- gsub("medial condyle width \\(mm\\): "        ,  "",GSE103328i$V26)
GSE103328i$V27 <- gsub("lateral condyle height \\(mm\\): "      ,  "",GSE103328i$V27)
GSE103328i$V28 <- gsub("lateral condyle depth \\(mm\\): "       ,  "",GSE103328i$V28)
GSE103328i$V29 <- gsub("lateral condyle width \\(mm\\): "       ,  "",GSE103328i$V29)
GSE103328i$V30 <- gsub("intercondylar notch width \\(mm\\): "   ,  "",GSE103328i$V30)
GSE103328i$V31 <- gsub("intercondylar notch depth \\(mm\\): "   ,  "",GSE103328i$V31)
GSE103328i$V32 <- gsub("bicondylar width \\(mm\\): "            ,  "",GSE103328i$V32)
colnames(GSE103328i) <- c("Sex","Age","Tissue","FL","BFL","MFL","SSWd","SSDp","MSWd","MSDp","ISWd","ISDp","HH","HL","HW","ANL","ANH","AND","BNL","PW","LTH","LTN","LGT","MCH","MCDp","MCW","LCH","LCDp","LCW","INW","IND","BCW")
head(GSE103328i); dim(GSE103328i) #69,9
table(colnames(GSE103328m)==rownames(GSE103328i))

GSE103279i <- dplyr::select(GSE103279i, colnames(GSE103271i))
GSE103271i$GEO <- "GSE103271"
GSE103279i$GEO <- "GSE103279"
GSE103280i$GEO <- "GSE103280"
GSE103328i$GEO <- "GSE103328"

COMMON_CpGs <- sort(Reduce(intersect, list(rownames(GSE103271m),rownames(GSE103279m),rownames(GSE103280m),rownames(GSE103328m)))); length(COMMON_CpGs)
GSE103271m <- GSE103271m[COMMON_CpGs,]; dim(GSE103271m)
GSE103279m <- GSE103279m[COMMON_CpGs,]; dim(GSE103279m)
GSE103280m <- GSE103280m[COMMON_CpGs,]; dim(GSE103280m)
GSE103328m <- GSE103328m[COMMON_CpGs,]; dim(GSE103328m)
table(rownames(GSE103271m)==rownames(GSE103279m))
table(rownames(GSE103271m)==rownames(GSE103280m))
table(rownames(GSE103271m)==rownames(GSE103328m))

Primate2_Inf <- rbind(GSE103271i,GSE103279i,GSE103280i,GSE103328i)
Primate2_Met <- cbind(GSE103271m,GSE103279m,GSE103280m,GSE103328m)
table(rownames(Primate2_Inf)==colnames(Primate2_Met))
Primate2_Inf$Sex <- toupper(substr(Primate2_Inf$Sex,1,1))
head(Primate2_Inf); dim(Primate2_Inf); table(Primate2_Inf$GEO) #54,32
table(Primate2_Inf$Tissue)
Primate2_Inf$Tissue <- NULL
Primate2_Met[1:5,1:5]; dim(Primate2_Met) #42485,54
table(Primate2_Inf$GEO)
summary(Primate2_Inf)
Sys.time()

setwd("G:/내 드라이브/Rrunning/GEO/HW014/Bone/")

#GSE138307
GSE138307m <- as.data.frame(fread("GSE138307_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=65))
rownames(GSE138307m) <- GSE138307m$ID_REF; GSE138307m <- GSE138307m[,-1]
GSE138307m <- GSE138307m[complete.cases(GSE138307m),]
GSE138307m[1:3,1:5]; dim(GSE138307m) #234452,32

GSE138307i <- as.data.frame(fread("GSE138307_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=28)); names(GSE138307i)[1] <- "ID"
GSE138307i <- dplyr::filter(GSE138307i, grepl("Sample_characteristics_ch1", ID))
GSE138307i <- as.data.frame(t(GSE138307i[,-1]))[,1:2]
table(colnames(GSE138307m)==rownames(GSE138307i))
GSE138307i$V1 <- gsub("age: ", "",GSE138307i$V1); GSE138307i$V1  <- gsub(" years",  "",GSE138307i$V1)
GSE138307i$V2 <- substr(GSE138307i$V2,9,9)
colnames(GSE138307i) <- c("Age", "Sex")
head(GSE138307i); dim(GSE138307i) #32,2 | ethnicity: Korean-Asian, tissue: Femoral bone

#GSE64490
GSE64490_m <- as.data.frame(fread("GSE64490_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=70))
rownames(GSE64490_m) <- GSE64490_m$ID_REF; GSE64490_m <- GSE64490_m[,-1]
GSE64490_m <- GSE64490_m[complete.cases(GSE64490_m),]
GSE64490_m[1:3,1:5]; dim(GSE64490_m) #190898,10

GSE64490_i <- as.data.frame(fread("GSE64490_series_matrix.txt.gz",    header = TRUE, sep = "\t", skip=32)); names(GSE64490_i)[1] <- "ID"
GSE64490_i <- dplyr::filter(GSE64490_i, grepl("Sample_characteristics_ch1", ID))
GSE64490_i <- as.data.frame(t(GSE64490_i[,-1]))[,c(3,5)]
GSE64490_i$V3 <- toupper(substr(GSE64490_i$V3,6,6))
GSE64490_i$V5 <- gsub("age: ", "",GSE64490_i$V5)
colnames(GSE64490_i) <- c("Sex", "Age")
head(GSE64490_i); dim(GSE64490_i)
table(colnames(GSE64490_m)==rownames(GSE64490_i))

GSE103286i$Age <- as.numeric(GSE103286i$Age); GSE103286i$Weight <- as.numeric(GSE103286i$Weight)
GSE138307i$Age <- as.numeric(GSE138307i$Age); GSE64490_i$Age <- as.numeric(GSE64490_i$Age)
for(j in 2:31){Primate2_Inf[,j] <- suppressWarnings(as.numeric(as.character(Primate2_Inf[,j])))}

setwd("G:/내 드라이브/Rrunning/HW014/Output_Ev15_ALL")
library(sva); library(ggplot2)
mod <- model.matrix(~1, data=Primate2_Inf); pca_raw <- prcomp(t(na.omit(Primate2_Met)), center=TRUE, scale.=TRUE)
pca_df <- data.frame(PC1=pca_raw$x[,1], PC2=pca_raw$x[,2], GEO=Primate2_Inf$GEO); p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=GEO)) + geom_point(size=3, alpha=0.8) + theme_bw() + labs(title="PCA by GEO Batch Before Correction")
ggsave("03_Primate_Batch_Check_PCA.pdf", plot=p_pca, width=6, height=5); write.table(pca_df, "03_Primate_PCA_Stats.txt", sep="\t", quote=FALSE, row.names=TRUE)
Primate2_Met_Cor <- ComBat(dat=as.matrix(na.omit(Primate2_Met)), batch=Primate2_Inf$GEO, mod=mod, par.prior=TRUE, prior.plots=FALSE)
Primate2_Met_Cor[Primate2_Met_Cor < 0] <- 0; Primate2_Met_Cor[Primate2_Met_Cor > 1] <- 1
write.table(Primate2_Met_Cor, "03_Primate_Met_ComBat_Corrected.txt", sep="\t", quote=FALSE, col.names=NA)

pca_cor <- prcomp(t(Primate2_Met_Cor), center=TRUE, scale.=TRUE); pca_cor_df <- data.frame(PC1=pca_cor$x[,1], PC2=pca_cor$x[,2], GEO=Primate2_Inf$GEO)
p_pca_cor <- ggplot(pca_cor_df, aes(x=PC1, y=PC2, color=GEO)) + geom_point(size=3, alpha=0.8) + theme_bw() + labs(title="PCA by GEO Batch After Correction")
ggsave("04_Primate_Batch_Check_PCA_After.pdf", plot=p_pca_cor, width=6, height=5); write.table(pca_cor_df, "04_Primate_PCA_Stats_After.txt", sep="\t", quote=FALSE, row.names=TRUE)




head(Primate2_Met_Cor, 3); dim(Primate2_Met_Cor)
head(Primate2_Inf, 2); dim(Primate2_Inf); table(Primate2_Inf$Sex)
table(colnames(Primate2_Met_Cor)==rownames(Primate2_Inf))
table(Primate2_Inf$GEO)

GSE103286i <- dplyr::filter(GSE103286i, Tissue=="bone")
GSE103286m <- dplyr::select(GSE103286m, rownames(GSE103286i))

GSE103286m[1:3,1:5]; dim(GSE103286m)
head(GSE103286i); dim(GSE103286i); table(GSE103286i$Tissue); table(GSE103286i$Sex)
table(colnames(GSE103286m)==rownames(GSE103286i))

GSE103271m[1:3,1:5]; dim(GSE103271m)
head(GSE103271i); dim(GSE103271i); table(GSE103271i$Sex)
table(colnames(GSE103271m)==rownames(GSE103271i))

GSE103279m[1:3,1:5]; dim(GSE103279m)
head(GSE103279i); dim(GSE103279i); table(GSE103279i$Sex)
table(colnames(GSE103279m)==rownames(GSE103279i))

GSE103280m[1:3,1:5]; dim(GSE103280m)
head(GSE103280i); dim(GSE103280i); table(GSE103280i$Sex)
table(colnames(GSE103280m)==rownames(GSE103280i))

GSE103328m[1:3,1:5]; dim(GSE103328m)
head(GSE103328i); dim(GSE103328i); table(GSE103328i$Sex)
table(colnames(GSE103328m)==rownames(GSE103328i))

GSE138307i <- dplyr::filter(GSE138307i, Sex=="M")
GSE138307m <- dplyr::select(GSE138307m, rownames(GSE138307i))
GSE138307m[1:3,1:5]; dim(GSE138307m)
head(GSE138307i); dim(GSE138307i); table(GSE138307i$Sex)
table(colnames(GSE138307m)==rownames(GSE138307i))

GSE64490_i <- dplyr::filter(GSE64490_i, Sex=="F")
GSE64490_m <- dplyr::select(GSE64490_m, rownames(GSE64490_i))
GSE64490_m[1:3,1:5]; dim(GSE64490_m)
head(GSE64490_i); dim(GSE64490_i); table(GSE64490_i$Sex)
table(colnames(GSE64490_m)==rownames(GSE64490_i))

Sys.time()
m_list <- list(Primate2_c=Primate2_Met_Cor, GSE103286c=GSE103286m, GSE138307c=GSE138307m, GSE64490_c=GSE64490_m); i_list <- list(Primate2_c=Primate2_Inf, GSE103286c=GSE103286i, GSE138307c=GSE138307i, GSE64490_c=GSE64490_i); s_list <- list(Primate2_c=c("F","M"), GSE103286c=c("F","M"), GSE138307c="M", GSE64490_c="F")
for(k in 1:length(m_list)){ nm <- names(m_list)[k]; mat <- as.matrix(m_list[[k]]); storage.mode(mat) <- "numeric"; mat <- mat[complete.cases(mat), , drop=FALSE]; inf <- i_list[[k]]; inf$Sex <- toupper(substr(inf$Sex, 1, 1)); inf$Age <- suppressWarnings(as.numeric(as.character(inf$Age))); out <- data.frame(CpG=rownames(mat), stringsAsFactors=FALSE)
for(sx in s_list[[k]]){ idx <- which(inf$Sex == sx & !is.na(inf$Age)); n <- length(idx); age <- inf$Age[idx]; sub_t <- t(mat[, idx, drop=FALSE])
rp <- suppressWarnings(as.numeric(cor(age, sub_t, method="pearson"))); pp <- 2*pt(-abs(rp*sqrt((n-2)/pmax(1e-15, 1-rp^2))), df=n-2); rs <- suppressWarnings(as.numeric(cor(age, sub_t, method="spearman"))); ps <- 2*pt(-abs(rs*sqrt((n-2)/pmax(1e-15, 1-rs^2))), df=n-2)
out[, paste0(sx, c("_Pearson_R","_Pearson_P","_Spearman_R","_Spearman_P"))] <- cbind(rp, pp, rs, ps)}; assign(nm, out)}

Primate2_c$M_HyPi <- (Primate2_c$M_Pearson_R * -log10(Primate2_c$M_Pearson_P)) + (Primate2_c$M_Spearman_R * -log10(Primate2_c$M_Spearman_P)); Primate2_c$F_HyPi <- (Primate2_c$F_Pearson_R * -log10(Primate2_c$F_Pearson_P)) + (Primate2_c$F_Spearman_R * -log10(Primate2_c$F_Spearman_P))
GSE103286c$M_HyPi <- (GSE103286c$M_Pearson_R * -log10(GSE103286c$M_Pearson_P)) + (GSE103286c$M_Spearman_R * -log10(GSE103286c$M_Spearman_P)); GSE103286c$F_HyPi <- (GSE103286c$F_Pearson_R * -log10(GSE103286c$F_Pearson_P)) + (GSE103286c$F_Spearman_R * -log10(GSE103286c$F_Spearman_P))
GSE138307c$M_HyPi <- (GSE138307c$M_Pearson_R * -log10(GSE138307c$M_Pearson_P)) + (GSE138307c$M_Spearman_R * -log10(GSE138307c$M_Spearman_P))
GSE64490_c$F_HyPi <- (GSE64490_c$F_Pearson_R * -log10(GSE64490_c$F_Pearson_P)) + (GSE64490_c$F_Spearman_R * -log10(GSE64490_c$F_Spearman_P))
head(Primate2_c); dim(Primate2_c) # 42485,11 | M 21, F 33
head(GSE103286c); dim(GSE103286c) #191954,11 | M 20, F 36
head(GSE138307c); dim(GSE138307c) #234452, 6 | M 28
head(GSE64490_c); dim(GSE64490_c) #485159, 6 |       F 46

Sys.time()

setwd("G:/내 드라이브/Rrunning/HW014/Output_Ev15_ALL/ModuleD")
################################################################################################
# Module D01. 7개 골 데이터셋과 TCGA 10대 고형암 공통 CpG 교차 추출 및 메타데이터 통합 구조화
################################################################################################
tcga_cans <- c("BLCA","BRCA","COAD","HNSC","KIRC","LIHC","LUAD","LUSC","PRAD","THCA")
bone_cgs <- Reduce(intersect, list(rownames(GSE64490_m), rownames(GSE138307m), rownames(GSE103286m), rownames(Primate2_Met_Cor)))
tcga_cgs <- bone_cgs; for(cn in tcga_cans){ tn <- paste0(cn,"_450m_Tm"); if(exists(tn)) tcga_cgs <- intersect(tcga_cgs, rownames(get(tn))) }
COM_D_CpG <- sort(tcga_cgs); write.table(data.frame(CpG=COM_D_CpG), "01_Bone_TCGA_Common_CpGs.txt", sep="\t", quote=FALSE, row.names=FALSE)
meta_d01 <- data.frame(Sample=c(colnames(GSE64490_m), colnames(GSE138307m), colnames(Primate2_Met_Cor), colnames(GSE103286m)), Dataset=c(rep("GSE64490",ncol(GSE64490_m)), rep("GSE138307",ncol(GSE138307m)), Primate2_Inf$GEO, rep("GSE103286",ncol(GSE103286m))), Organ="Bone", Condition="Healthy_Bone", Sex=c(rep("F",ncol(GSE64490_m)), rep("M",ncol(GSE138307m)), Primate2_Inf$Sex, GSE103286i$Sex), Age=c(GSE64490_i$Age, GSE138307i$Age, Primate2_Inf$Age, GSE103286i$Age), Cancer="None")
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm"); cl_n <- paste0(cn,"_c")
if(exists(tm_n) && exists(cl_n)){ c_df <- get(cl_n); m_tm <- get(tm_n); s_tm <- intersect(colnames(m_tm), c_df$ID); if(length(s_tm)>0) meta_d01 <- rbind(meta_d01, data.frame(Sample=s_tm, Dataset=paste0("TCGA_",cn), Organ="Solid_Organ", Condition="Tumor", Sex=c_df$Sex[match(s_tm,c_df$ID)], Age=c_df$Age[match(s_tm,c_df$ID)], Cancer=cn)) }
if(exists(nm_n) && exists(cl_n)){ c_df <- get(cl_n); m_nm <- get(nm_n); s_nm <- intersect(colnames(m_nm), c_df$ID); if(length(s_nm)>0) meta_d01 <- rbind(meta_d01, data.frame(Sample=s_nm, Dataset=paste0("TCGA_",cn), Organ="Solid_Organ", Condition="Normal", Sex=c_df$Sex[match(s_nm,c_df$ID)], Age=c_df$Age[match(s_nm,c_df$ID)], Cancer=cn)) } }
meta_d01$Sex <- toupper(substr(meta_d01$Sex,1,1)); meta_d01$Age <- suppressWarnings(as.numeric(as.character(meta_d01$Age))); meta_d01 <- meta_d01[meta_d01$Sex %in% c("F","M") & is.finite(meta_d01$Age), ]
write.table(as.data.frame(table(meta_d01$Dataset, meta_d01$Condition)), "01_Cohort_Sample_Summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
p01 <- ggplot(meta_d01, aes(x=Dataset, fill=Condition)) + geom_bar(position="stack", color="black", width=0.7) + scale_fill_manual(values=c("Healthy_Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) + labs(title=paste0("Bone & TCGA Solid Cancer Integration Cohort Size (", length(COM_D_CpG), " Shared CpGs)"), x="Dataset / Cohort", y="Sample Count", fill="Condition"); ggsave("01_Bone_TCGA_Cohort_Distribution.pdf", plot=p01, width=11, height=5.5)

################################################################################################
# Module D02. 정상 골격계 기준선 기반 후성유전 시계 구축 및 고형암 정상 조직 간 노화 속도 비교
################################################################################################
mat_hbone <- cbind(as.matrix(GSE64490_m[COM_D_CpG, ]), as.matrix(GSE138307m[COM_D_CpG, ])); storage.mode(mat_hbone) <- "numeric"; age_hbone <- c(GSE64490_i$Age, GSE138307i$Age)
r_bone <- suppressWarnings(as.numeric(cor(age_hbone, t(mat_hbone)))); names(r_bone) <- COM_D_CpG; top_clk_cgs <- COM_D_CpG[order(abs(r_bone), decreasing=TRUE)[1:50]]; w_clk <- r_bone[top_clk_cgs]
clk_res02 <- data.frame(Sample=colnames(mat_hbone), Cohort="Human_Bone", Tissue="Healthy_Bone", Cancer="Bone", Age=age_hbone, Clock_Score=as.numeric(crossprod(mat_hbone[top_clk_cgs, ], w_clk)/sum(abs(w_clk))))
for(cn in tcga_cans){ nm_n <- paste0(cn,"_450m_Nm"); cl_n <- paste0(cn,"_c")
if(exists(nm_n) && exists(cl_n)){ c_df <- get(cl_n); mat_nm <- as.matrix(get(nm_n)[top_clk_cgs, , drop=FALSE]); storage.mode(mat_nm) <- "numeric"; s_ok <- intersect(colnames(mat_nm), c_df$ID)
if(length(s_ok)>=5){ a_v <- c_df$Age[match(s_ok, c_df$ID)]; sc_v <- as.numeric(crossprod(mat_nm[, s_ok], w_clk)/sum(abs(w_clk))); clk_res02 <- rbind(clk_res02, data.frame(Sample=s_ok, Cohort=paste0("TCGA_",cn,"_Norm"), Tissue="Solid_Normal", Cancer=cn, Age=a_v, Clock_Score=sc_v)) } } }
clk_res02 <- clk_res02[is.finite(clk_res02$Age) & is.finite(clk_res02$Clock_Score), ]; stat02 <- data.frame()
for(grp in unique(clk_res02$Cohort)){ sub_d <- clk_res02[clk_res02$Cohort==grp, ]; ft <- lm(Clock_Score ~ Age, data=sub_d); ct <- cor.test(sub_d$Age, sub_d$Clock_Score)
stat02 <- rbind(stat02, data.frame(Cohort=grp, N=nrow(sub_d), Slope=coef(ft)[2], Intercept=coef(ft)[1], R=ct$estimate, P_val=ct$p.value)) }
write.table(stat02, "02_Bone_vs_TCGA_Normal_Aging_Slope_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p02 <- ggplot(clk_res02, aes(x=Age, y=Clock_Score, color=Cohort)) + geom_point(alpha=0.6, size=1.5) + geom_smooth(method="lm", se=FALSE, size=1) + facet_wrap(~Cancer, scales="free_y", ncol=4) + theme_bw() + labs(title="Skeletal Aging Clock (50 Loci) Trajectory: \nHealthy Bone vs TCGA Adjacent Normals", x="Chronological Age (Years)", y="Predicted Epigenetic Aging Score"); ggsave("02_Bone_vs_TCGA_Normal_Aging_Slope_Plot.pdf", plot=p02, width=6, height=6)

################################################################################################
# Module D03. 정상 뼈 레퍼런스 대비 10대 고형암의 후성유전 노화 편차(EAD Z-score) 및 가속도 산출 (수정본)
################################################################################################
fit_bone_ref <- lm(Clock_Score ~ Age, data=clk_res02[clk_res02$Cohort=="Human_Bone", ]); sd_bone_res <- sd(resid(fit_bone_ref))
ead_all03 <- data.frame(); for(cn in tcga_cans){ for(st in c("Tm","Nm")){ m_name <- paste0(cn,"_450m_",st); cl_n <- paste0(cn,"_c")
if(!exists(m_name)){ for(sx in c("M","F")){ alt_n <- paste0(cn,"_450m_",st,sx); if(exists(alt_n)){ assign(m_name, get(alt_n)) } } }
if(exists(m_name) && exists(cl_n)){ c_df <- get(cl_n); mat_raw <- get(m_name); mat_s <- as.matrix(mat_raw[intersect(top_clk_cgs, rownames(mat_raw)), , drop=FALSE]); storage.mode(mat_s) <- "numeric"
m_pats <- toupper(substr(gsub("[_\\.]", "-", colnames(mat_s)), 1, 12)); c_pats <- toupper(substr(gsub("[_\\.]", "-", c_df$ID), 1, 12)); com_pats <- intersect(m_pats, c_pats)
if(length(com_pats)>=3){ idx_m <- match(com_pats, m_pats); idx_c <- match(com_pats, c_pats); sub_mat <- mat_s[, idx_m, drop=FALSE]; a_v <- suppressWarnings(as.numeric(as.character(c_df$Age[idx_c]))); sx_v <- toupper(substr(c_df$Sex[idx_c],1,1))
sc_v <- as.numeric(crossprod(sub_mat[top_clk_cgs, ], w_clk)/sum(abs(w_clk))); ok_s <- which(is.finite(a_v) & is.finite(sc_v) & sx_v %in% c("F","M"))
if(length(ok_s)>=3){ pred_v <- predict(fit_bone_ref, newdata=data.frame(Age=a_v[ok_s])); eaa_v <- sc_v[ok_s] - pred_v; ead_v <- eaa_v / sd_bone_res
ead_all03 <- rbind(ead_all03, data.frame(Sample=colnames(sub_mat)[ok_s], Cancer=cn, State=ifelse(st=="Tm","Tumor","Normal"), Sex=sx_v[ok_s], Age=a_v[ok_s], Clock_Score=sc_v[ok_s], EAA=eaa_v, EAD=ead_v)) } } } } }
if(nrow(ead_all03) > 0){ write.table(ead_all03, "03_Sample_Level_EAD_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  stat03 <- data.frame(); for(cn in unique(ead_all03$Cancer)){ for(st in unique(ead_all03$State)){ v <- ead_all03$EAD[ead_all03$Cancer==cn & ead_all03$State==st]; if(length(v)>0) stat03 <- rbind(stat03, data.frame(Cancer=cn, State=st, Mean=mean(v), SD=sd(v), Median=median(v), N=length(v))) } }
  write.table(stat03, "03_Cancer_State_EAD_Summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p03 <- ggplot(ead_all03, aes(x=reorder(Cancer, EAD, FUN=median), y=EAD, fill=State)) + geom_boxplot(outlier.size=0.5, alpha=0.85) + geom_hline(yintercept=0, linetype="dashed", color="darkgreen", size=1) + scale_fill_manual(values=c("Normal"="#377EB8","Tumor"="#E41A1C")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Epigenetic Age Deviation (EAD Z-score) Referenced to Healthy Human Bone Baseline", subtitle="Dashed line at Z = 0 indicates Healthy Bone Epigenetic Age Expectation", x="TCGA Solid Cancer Type", y="EAD (Z-score vs Bone Baseline)"); ggsave("03_Bone_Referenced_Cancer_EAD_Plot.pdf", plot=p03, width=10, height=5.5)
} else { cat("경고: TCGA 암 검체와 임상 데이터 간 매칭된 검체가 없어 ead_all03이 생성되지 않았습니다. 워크스페이스의 객체명을 확인하세요.\n") }

################################################################################################
# Module D04. 고형암 종양 조직의 골격 대비 후성유전 노화 가속도 성별 이형성 및 병기 연관성 분석 (수정본)
################################################################################################
if(!exists("ead_all03") || nrow(ead_all03) == 0){ cat("경고: ead_all03이 비어 있습니다. Module D03 수정본을 먼저 성공적으로 실행해 주세요.\n")
} else {
  ead_t04 <- ead_all03[ead_all03$State == "Tumor", ]; ead_t04$PatID <- toupper(substr(gsub("[_\\.]", "-", ead_t04$Sample), 1, 12))
  stage_info04 <- data.frame(); for(cn in tcga_cans){ cl_n <- paste0(cn,"_c")
  if(exists(cl_n)){ c_df <- get(cl_n); if("StageC" %in% colnames(c_df)){ pid <- toupper(substr(gsub("[_\\.]", "-", c_df$ID), 1, 12)); stage_info04 <- rbind(stage_info04, data.frame(PatID=pid, Cancer=cn, StageC=as.character(c_df$StageC))) } } }
  stage_info04 <- stage_info04[!duplicated(paste(stage_info04$PatID, stage_info04$Cancer, sep="_")), ]
  m_stage04 <- merge(ead_t04, stage_info04, by=c("PatID","Cancer")); m_stage04 <- m_stage04[m_stage04$StageC %in% c("L","M"), ]
  stat04_sex <- data.frame(); for(cn in unique(ead_t04$Cancer)){ sub_c <- ead_t04[ead_t04$Cancer==cn, ]; vm <- sub_c$EAD[sub_c$Sex=="M"]; vf <- sub_c$EAD[sub_c$Sex=="F"]
  if(length(vm)>=3 && length(vf)>=3 && (var(vm,na.rm=TRUE)>0 || var(vf,na.rm=TRUE)>0)){ tt <- t.test(vm, vf); stat04_sex <- rbind(stat04_sex, data.frame(Cancer=cn, Male_N=length(vm), Female_N=length(vf), Male_Mean_EAD=mean(vm,na.rm=TRUE), Female_Mean_EAD=mean(vf,na.rm=TRUE), Delta_M_minus_F=mean(vm,na.rm=TRUE)-mean(vf,na.rm=TRUE), P_value=tt$p.value)) } }
  stat04_stage <- data.frame(); for(cn in unique(m_stage04$Cancer)){ sub_s <- m_stage04[m_stage04$Cancer==cn, ]; vl <- sub_s$EAD[sub_s$StageC=="L"]; vm <- sub_s$EAD[sub_s$StageC=="M"]
  if(length(vl)>=3 && length(vm)>=3 && (var(vl,na.rm=TRUE)>0 || var(vm,na.rm=TRUE)>0)){ tt <- t.test(vm, vl); stat04_stage <- rbind(stat04_stage, data.frame(Cancer=cn, Early_Stage_N=length(vl), Late_Stage_N=length(vm), Early_Mean_EAD=mean(vl,na.rm=TRUE), Late_Mean_EAD=mean(vm,na.rm=TRUE), Delta_Late_minus_Early=mean(vm,na.rm=TRUE)-mean(vl,na.rm=TRUE), P_value=tt$p.value)) } }
  write.table(stat04_sex, "04_Tumor_EAD_Sex_Dimorphism_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE); write.table(stat04_stage, "04_Tumor_EAD_Clinical_Stage_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p04a <- ggplot(ead_t04, aes(x=Cancer, y=EAD, fill=Sex)) + geom_boxplot(outlier.size=0.4, alpha=0.8) + scale_fill_manual(values=c("M"="#0D3692","F"="#CB333B")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Sex Dimorphism in Skeletal-Referenced Tumor EAD", x="Cancer", y="Tumor EAD (Z)")
  p04b <- ggplot(m_stage04, aes(x=Cancer, y=EAD, fill=StageC)) + geom_boxplot(outlier.size=0.4, alpha=0.8) + scale_fill_manual(values=c("L"="#4DAC26","M"="#CA0020"), labels=c("L"="Early (I/II)", "M"="Advanced (III/IV)")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Tumor EAD Across Clinical Stages", x="Cancer", y="Tumor EAD (Z)", fill="Stage")
  ggsave("04_Tumor_EAD_Sex_Stage_Plot.pdf", plot=gridExtra::grid.arrange(p04a, p04b, ncol=2), width=13, height=5.5) }

################################################################################################
# Module D05. 정상 뼈 vs 암 인접 정상 조직(ANFE) vs 원발 고형암 간 점진적 메틸화 왜곡 및 현장 결손 분석
################################################################################################
mean_hbone <- rowMeans(mat_hbone, na.rm=TRUE); anfe_stat05 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[COM_D_CpG, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[COM_D_CpG, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=5 && ncol(m_nm)>=5){ mean_tm <- rowMeans(m_tm, na.rm=TRUE); mean_nm <- rowMeans(m_nm, na.rm=TRUE)
delta_field <- mean_nm - mean_hbone; delta_tumor <- mean_tm - mean_nm; conc_field <- (delta_field * delta_tumor) > 0
sig_anfe <- abs(delta_field) > 0.05 & abs(delta_tumor) > 0.1 & conc_field; ct <- cor.test(delta_field, delta_tumor)
anfe_stat05 <- rbind(anfe_stat05, data.frame(Cancer=cn, N_Normal=ncol(m_nm), N_Tumor=ncol(m_tm), Cor_Field_Tumor_R=ct$estimate, P_val=ct$p.value, Progressive_Field_CpGs=sum(sig_anfe, na.rm=TRUE), Progressive_Hyper=sum(sig_anfe & delta_field>0, na.rm=TRUE), Progressive_Hypo=sum(sig_anfe & delta_field<0, na.rm=TRUE))) } } }
write.table(anfe_stat05, "05_Adjacent_Normal_Field_Defect_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
top_can05 <- anfe_stat05$Cancer[which.max(anfe_stat05$Progressive_Field_CpGs)]; m_tm_top <- as.matrix(get(paste0(top_can05,"_450m_Tm"))[COM_D_CpG, ]); m_nm_top <- as.matrix(get(paste0(top_can05,"_450m_Nm"))[COM_D_CpG, ]); storage.mode(m_tm_top) <- "numeric"; storage.mode(m_nm_top) <- "numeric"
plot_df05 <- data.frame(CpG=COM_D_CpG, Delta_Field=rowMeans(m_nm_top,na.rm=TRUE)-mean_hbone, Delta_Tumor=rowMeans(m_tm_top,na.rm=TRUE)-rowMeans(m_nm_top,na.rm=TRUE))
plot_df05$Class <- "Background"; plot_df05$Class[plot_df05$Delta_Field > 0.05 & plot_df05$Delta_Tumor > 0.05] <- "Progressive Hypermethylation"; plot_df05$Class[plot_df05$Delta_Field < -0.05 & plot_df05$Delta_Tumor < -0.05] <- "Progressive Hypomethylation"
p05a <- ggplot(anfe_stat05, aes(x=reorder(Cancer, Progressive_Field_CpGs), y=Progressive_Field_CpGs, fill=Cor_Field_Tumor_R)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + scale_fill_gradient2(low="#377EB8", mid="white", high="#E41A1C", midpoint=0) + theme_bw() + labs(title="Progressive Field Defect Burden Across Solid Organs (vs Bone)", x="Cancer Type", y="Progressive Field CpGs (|Delta Beta| > 0.05)", fill="Correlation (R)")
p05b <- ggplot(plot_df05[seq(1,nrow(plot_df05),length.out=min(5000,nrow(plot_df05))),], aes(x=Delta_Field, y=Delta_Tumor, color=Class)) + geom_point(alpha=0.5, size=0.8) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + scale_color_manual(values=c("Progressive Hypermethylation"="#E41A1C","Progressive Hypomethylation"="#377EB8","Background"="grey80")) + theme_bw() + labs(title=paste0("Field Defect Vector in ", top_can05, " (Normal vs Bone vs Tumor)"), x="Delta Beta: Adjacent Normal - Healthy Bone", y="Delta Beta: Primary Tumor - Adjacent Normal")
ggsave("05_Bone_Normal_Tumor_FieldDefect_Plot.pdf", plot=gridExtra::grid.arrange(p05a, p05b, ncol=2), width=13, height=5.5)

################################################################################################
# Module D06. 골 전이 호발암(PRAD, BRCA, LUAD) vs 비호발암 간 골 지향성(Osteotropism) 메틸화 시그니처 발굴 (수정본)
################################################################################################
mat_hb <- cbind(as.matrix(GSE64490_m[COM_D_CpG, ]), as.matrix(GSE138307m[COM_D_CpG, ])); storage.mode(mat_hb) <- "numeric"; mean_hb <- rowMeans(mat_hb, na.rm=TRUE); names(mean_hb) <- COM_D_CpG
tropic_cans <- c("BRCA","PRAD","LUAD"); nontropic_cans <- c("COAD","KIRC","LIHC","BLCA")
dist06_df <- data.frame(); stat06 <- data.frame()
for(cn in c(tropic_cans, nontropic_cans)){ cand_list <- list()
for(sx in c("F","M")){ nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f")
if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); cand_list[[sx]] <- d[, f=="T", drop=FALSE] }
else if(exists(paste0(cn,"_450m_Tm",sx))){ cand_list[[sx]] <- get(paste0(cn,"_450m_Tm",sx)) }
else if(exists(paste0(cn,"_450m_Tm",sx,"1"))){ cand_list[[sx]] <- get(paste0(cn,"_450m_Tm",sx,"1")) } }
if(length(cand_list)==0 && exists(paste0(cn,"_450m_Tm"))){ cand_list[["all"]] <- get(paste0(cn,"_450m_Tm")) }
if(length(cand_list)>0){ cg_cn <- COM_D_CpG; for(k in 1:length(cand_list)){ cg_cn <- intersect(cg_cn, rownames(cand_list[[k]])) }
if(length(cg_cn)>=20){ m_tm <- NULL; for(k in 1:length(cand_list)){ sub_m <- as.matrix(cand_list[[k]][cg_cn, , drop=FALSE]); storage.mode(sub_m) <- "numeric"; m_tm <- if(is.null(m_tm)) sub_m else cbind(m_tm, sub_m) }
if(!is.null(m_tm) && ncol(m_tm)>=2){ hb_sub <- mean_hb[cg_cn]; d_vals <- sqrt(colMeans((m_tm - hb_sub)^2, na.rm=TRUE)); r_vals <- suppressWarnings(as.numeric(cor(hb_sub, m_tm)))
grp_t <- ifelse(cn %in% tropic_cans, "Bone-Tropic", "Non-Tropic")
dist06_df <- rbind(dist06_df, data.frame(Sample=colnames(m_tm), Cancer=cn, Tropism=grp_t, RMS_Distance=d_vals, Bone_Cor_R=r_vals, stringsAsFactors=FALSE)) } } } }
dist06_df <- dist06_df[is.finite(dist06_df$RMS_Distance) & is.finite(dist06_df$Bone_Cor_R), ]
for(cn in unique(dist06_df$Cancer)){ sub_c <- dist06_df[dist06_df$Cancer==cn, ]; stat06 <- rbind(stat06, data.frame(Cancer=cn, Tropism=sub_c$Tropism[1], N=nrow(sub_c), Mean_RMS_Dist=mean(sub_c$RMS_Distance), SD_RMS_Dist=sd(sub_c$RMS_Distance), Mean_Bone_R=mean(sub_c$Bone_Cor_R))) }
v_tr <- dist06_df$RMS_Distance[dist06_df$Tropism=="Bone-Tropic"]; v_nt <- dist06_df$RMS_Distance[dist06_df$Tropism=="Non-Tropic"]
tt06_d_p <- if(length(v_tr)>=2 && length(v_nt)>=2 && var(v_tr, na.rm=TRUE)>0 && var(v_nt, na.rm=TRUE)>0) t.test(v_tr, v_nt)$p.value else NA
r_tr <- dist06_df$Bone_Cor_R[dist06_df$Tropism=="Bone-Tropic"]; r_nt <- dist06_df$Bone_Cor_R[dist06_df$Tropism=="Non-Tropic"]
tt06_r_p <- if(length(r_tr)>=2 && length(r_nt)>=2 && var(r_tr, na.rm=TRUE)>0 && var(r_nt, na.rm=TRUE)>0) t.test(r_tr, r_nt)$p.value else NA
write.table(stat06, "06_Osteotropism_Signature_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p06a <- ggplot(dist06_df, aes(x=reorder(Cancer, RMS_Distance), y=RMS_Distance, fill=Tropism)) + geom_boxplot(alpha=0.8, outlier.size=0.4) + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title=paste0("Methylation Divergence from Healthy Bone (T-test P = ", ifelse(is.na(tt06_d_p), "NA", signif(tt06_d_p, 3)), ")"), x="Cancer Type", y="RMS Distance to Bone Baseline (Delta Beta)")
p06b <- ggplot(dist06_df, aes(x=reorder(Cancer, -Bone_Cor_R), y=Bone_Cor_R, fill=Tropism)) + geom_boxplot(alpha=0.8, outlier.size=0.4) + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title=paste0("Bone-Mimicry Correlation (T-test P = ", ifelse(is.na(tt06_r_p), "NA", signif(tt06_r_p, 3)), ")"), x="Cancer Type", y="Pearson R with Normal Bone Mean Profile")
ggsave("06_Osteotropism_Bone_Mimicry_Plot.pdf", plot=gridExtra::grid.arrange(p06a, p06b, ncol=2), width=13, height=5.5)

################################################################################################
# Module D07. 골격계 노화 유도성 CpG(HyPi)와 고형암 발암성 메틸화 변이(Delta Beta) 간 방향성 일치도 검증 (수정본)
################################################################################################
hp_m <- GSE138307c$M_HyPi[match(COM_D_CpG, GSE138307c$CpG)]; hp_f <- GSE64490_c$F_HyPi[match(COM_D_CpG, GSE64490_c$CpG)]
hm_hp <- (hp_m + hp_f) / 2; names(hm_hp) <- COM_D_CpG; hm_hp <- hm_hp[is.finite(hm_hp)]
stat07 <- data.frame(); plot07_df <- data.frame()
for(cn in tcga_cans){ cand_t <- list(); cand_n <- list()
for(sx in c("F","M")){ nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f")
if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); cand_t[[paste0(sx,"_T")]] <- d[, f=="T", drop=FALSE]; cand_n[[paste0(sx,"_N")]] <- d[, f=="N", drop=FALSE] }
tm_s <- paste0(cn,"_450m_Tm",sx); nm_s <- paste0(cn,"_450m_Nm",sx)
if(length(cand_t)==0 && exists(tm_s)) cand_t[[paste0(sx,"_Tm")]] <- get(tm_s)
if(length(cand_n)==0 && exists(nm_s)) cand_n[[paste0(sx,"_Nm")]] <- get(nm_s) }
if(length(cand_t)==0 && exists(paste0(cn,"_450m_Tm"))) cand_t[["all_T"]] <- get(paste0(cn,"_450m_Tm"))
if(length(cand_n)==0 && exists(paste0(cn,"_450m_Nm"))) cand_n[["all_N"]] <- get(paste0(cn,"_450m_Nm"))
if(length(cand_t)>0 && length(cand_n)>0){ all_m <- c(cand_t, cand_n); com_cg <- names(hm_hp)
for(k in 1:length(all_m)){ com_cg <- intersect(com_cg, rownames(all_m[[k]])) }
if(length(com_cg)>=20){ m_tm <- NULL; for(k in 1:length(cand_t)){ sub_m <- as.matrix(cand_t[[k]][com_cg, , drop=FALSE]); m_tm <- if(is.null(m_tm)) sub_m else cbind(m_tm, sub_m) }
m_nm <- NULL; for(k in 1:length(cand_n)){ sub_m <- as.matrix(cand_n[[k]][com_cg, , drop=FALSE]); m_nm <- if(is.null(m_nm)) sub_m else cbind(m_nm, sub_m) }
storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=2 && ncol(m_nm)>=2){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE); sub_hp <- hm_hp[com_cg]
sig_idx <- which(is.finite(sub_hp) & is.finite(db) & abs(sub_hp) > 1.0 & abs(db) > 0.02)
if(length(sig_idx)>=10){ ct <- cor.test(sub_hp[sig_idx], db[sig_idx]); conc_pct <- mean(sub_hp[sig_idx] * db[sig_idx] > 0) * 100
stat07 <- rbind(stat07, data.frame(Cancer=cn, Analyzed_CpGs=length(sig_idx), Cor_Aging_Cancer_R=ct$estimate, P_val=ct$p.value, Concordant_Pct=conc_pct, Hyper_Concordant=sum(sub_hp[sig_idx]>0 & db[sig_idx]>0), Hypo_Concordant=sum(sub_hp[sig_idx]<0 & db[sig_idx]<0), stringsAsFactors=FALSE))
sub_p <- data.frame(Cancer=cn, CpG=names(sub_hp)[sig_idx], Bone_Aging_HyPi=sub_hp[sig_idx], Cancer_Delta_Beta=db[sig_idx], Concordance=ifelse(sub_hp[sig_idx]*db[sig_idx]>0, "Concordant (Aging Mimicry)", "Discordant"), stringsAsFactors=FALSE)
plot07_df <- rbind(plot07_df, sub_p[seq(1, nrow(sub_p), length.out=min(800, nrow(sub_p))), ]) } } } } }
if(nrow(stat07) > 0){ write.table(stat07, "07_Aging_vs_Oncogenesis_Concordance_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p07a <- ggplot(stat07, aes(x=reorder(Cancer, Concordant_Pct), y=Concordant_Pct, fill=Cor_Aging_Cancer_R)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=50, linetype="dashed", color="red") + scale_fill_gradient2(low="#377EB8", mid="white", high="#E41A1C", midpoint=0) + theme_bw() + labs(title="Directional Concordance (% Same Sign): \nBone Aging HyPi vs Cancer Delta Beta", x="Solid Cancer", y="Concordance Rate (%)", fill="Correlation (R)")
  p07b <- ggplot(plot07_df, aes(x=Bone_Aging_HyPi, y=Cancer_Delta_Beta, color=Concordance)) + geom_point(alpha=0.5, size=0.7) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + facet_wrap(~Cancer, scales="free_y", ncol=min(5, length(unique(plot07_df$Cancer)))) + scale_color_manual(values=c("Concordant (Aging Mimicry)"="#E41A1C","Discordant"="#377EB8")) + theme_bw() + theme(legend.position="bottom") + labs(title="Bone Aging HyPi vs Primary Tumor Delta Beta Across Cancers", x="Human Bone Aging HyPi", y="Cancer Delta Beta (Tumor - Normal)")
#  ggsave("07_Aging_vs_Oncogenesis_Concordance_Plot.pdf", plot=gridExtra::grid.arrange(p07a, p07b, ncol=1, heights=c(1, 1.4)), width=13, height=8.5)
  ggsave("07a_Aging_vs_Oncogenesis_Concordance_Plot.pdf", plot=p07a, width= 5, height=5)
  ggsave("07b_Aging_vs_Oncogenesis_Concordance_Plot.pdf", plot=p07b, width=13, height=5)
} else { cat("경고: 유효한 정상/종양 쌍이 매칭된 암종이 없어 stat07이 비어 있습니다. 워크스페이스 내 객체를 확인하세요.\n") }

################################################################################################
# Module D08. 정상 뼈 메틸화 기저선 대비 10대 고형암 종양 조직의 후성유전학적 엔트로피(Shannon Entropy) 붕괴 측정
################################################################################################
ent_res08 <- data.frame(); b_hb <- pmin(pmax(mat_hb, 1e-6), 1-1e-6)
ent_hb <- colMeans(-(b_hb*log2(b_hb) + (1-b_hb)*log2(1-b_hb)), na.rm=TRUE)
ent_res08 <- rbind(ent_res08, data.frame(Sample=colnames(mat_hb), Cancer="Human_Bone", Group="Healthy_Bone", Entropy=ent_hb))
for(cn in tcga_cans){ for(st in c("Tm","Nm")){ m_name <- paste0(cn,"_450m_",st)
if(exists(m_name)){ mat_c <- as.matrix(get(m_name)[COM_D_CpG, , drop=FALSE]); storage.mode(mat_c) <- "numeric"
if(ncol(mat_c)>=3){ b_c <- pmin(pmax(mat_c, 1e-6), 1-1e-6); ent_c <- colMeans(-(b_c*log2(b_c) + (1-b_c)*log2(1-b_c)), na.rm=TRUE)
ent_res08 <- rbind(ent_res08, data.frame(Sample=colnames(mat_c), Cancer=cn, Group=ifelse(st=="Tm","Tumor","Normal"), Entropy=ent_c)) } } } }
ent_res08 <- ent_res08[is.finite(ent_res08$Entropy), ]; stat08 <- data.frame(); ref_ent <- ent_res08$Entropy[ent_res08$Group=="Healthy_Bone"]
for(cn in tcga_cans){ vt <- ent_res08$Entropy[ent_res08$Cancer==cn & ent_res08$Group=="Tumor"]; vn <- ent_res08$Entropy[ent_res08$Cancer==cn & ent_res08$Group=="Normal"]
tt_t <- if(length(vt)>=3) t.test(vt, ref_ent)$p.value else NA; tt_n <- if(length(vn)>=3) t.test(vn, ref_ent)$p.value else NA
stat08 <- rbind(stat08, data.frame(Cancer=cn, N_Tumor=length(vt), N_Normal=length(vn), Mean_Tumor_Entropy=ifelse(length(vt)>0,mean(vt),NA), Mean_Normal_Entropy=ifelse(length(vn)>0,mean(vn),NA), Bone_Baseline_Entropy=mean(ref_ent), Delta_Tumor_vs_Bone=ifelse(length(vt)>0,mean(vt)-mean(ref_ent),NA), P_Tumor_vs_Bone=tt_t, P_Normal_vs_Bone=tt_n)) }
write.table(stat08, "08_Epigenetic_Entropy_Collapse_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p08a <- ggplot(ent_res08, aes(x=reorder(Cancer, Entropy, FUN=median), y=Entropy, fill=Group)) + geom_boxplot(outlier.size=0.4, alpha=0.85) + geom_hline(yintercept=mean(ref_ent), linetype="dashed", color="darkgreen", size=1) + scale_fill_manual(values=c("Healthy_Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Genome-Wide Shannon Entropy Deregulation: Solid Cancers vs Healthy Bone Baseline", x="Cohort / Cancer Type", y="Mean Epigenetic Entropy (Bits)", fill="Tissue State")
p08b <- ggplot(stat08, aes(x=reorder(Cancer, Delta_Tumor_vs_Bone), y=Delta_Tumor_vs_Bone, fill=Delta_Tumor_vs_Bone>0)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Entropy Inflation","FALSE"="Entropy Loss")) + theme_bw() + labs(title="Tumor Entropy Deviation Relative to Bone Baseline", x="Cancer", y="Delta Shannon Entropy (Tumor - Bone)", fill="Direction")
ggsave("08_Epigenetic_Entropy_Collapse_Plot.pdf", plot=gridExtra::grid.arrange(p08a, p08b, ncol=2), width=13, height=5.5)

################################################################################################
# Module D09. 골격 마스터 전사인자(RUNX2, SP7, SOX9, BMP) 프로모터의 암종별 이상 과메틸화 및 억제 해제 분석
################################################################################################
if(!exists("GPL13534")){ GPL13534 <- as.data.frame(fread("G:/내 드라이브/Rrunning/Annotations/21_Methyl/GPL13534-11288.txt.gz", header=TRUE, sep="\t"))[, c("ID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")] }
skel_genes <- c("RUNX2","SP7","SOX9","COL1A1","ALPL","BGLAP","BMP2","BMP4")
skel_anno <- GPL13534[GPL13534$ID %in% COM_D_CpG & grepl(paste(skel_genes, collapse="|"), as.character(GPL13534$UCSC_RefGene_Name)), ]
skel_anno$Gene <- sub(";.*", "", as.character(skel_anno$UCSC_RefGene_Name)); skel_cgs <- unique(skel_anno$ID); stat09 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[skel_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[skel_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ for(cg in skel_cgs){ gn <- skel_anno$Gene[skel_anno$ID==cg][1]; v_hb <- mat_hb[cg, ]; v_tm <- m_tm[cg, ]; v_nm <- m_nm[cg, ]
tt_tb <- t.test(v_tm, v_hb); tt_tn <- t.test(v_tm, v_nm)
stat09 <- rbind(stat09, data.frame(Cancer=cn, CpG=cg, Gene=gn, Island=as.character(skel_anno$Relation_to_UCSC_CpG_Island[skel_anno$ID==cg][1]), Mean_Bone=mean(v_hb,na.rm=TRUE), Mean_Normal=mean(v_nm,na.rm=TRUE), Mean_Tumor=mean(v_tm,na.rm=TRUE), Delta_Tumor_vs_Bone=mean(v_tm,na.rm=TRUE)-mean(v_hb,na.rm=TRUE), Delta_Tumor_vs_Normal=mean(v_tm,na.rm=TRUE)-mean(v_nm,na.rm=TRUE), P_Tumor_vs_Bone=tt_tb$p.value, P_Tumor_vs_Normal=tt_tn$p.value)) } } } }
write.table(stat09, "09_Skeletal_Master_Regulators_Cancer_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
gene_summary09 <- aggregate(Delta_Tumor_vs_Bone ~ Cancer + Gene, data=stat09, FUN=mean)
p09a <- ggplot(gene_summary09, aes(x=Cancer, y=Gene, fill=Delta_Tumor_vs_Bone)) + geom_tile(color="white") + geom_text(aes(label=round(Delta_Tumor_vs_Bone,2)), size=3.2) + scale_fill_gradient2(low="#377EB8", mid="white", high="#E41A1C", midpoint=0) + theme_minimal() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Master Skeletal Transcription Factor Methylation Shift (Tumor vs Healthy Bone)", subtitle="Mean Delta Beta: Red = Repressed/Hypermethylated in Tumor, Blue = Hypomethylated/Activated", x="Solid Cancer", y="Osteogenic / Chondrogenic Master Gene", fill="Delta Beta")
p09b <- ggplot(stat09, aes(x=Gene, y=Delta_Tumor_vs_Normal, fill=Gene)) + geom_boxplot(alpha=0.8, outlier.size=0.5) + geom_hline(yintercept=0, linetype="dashed", color="red") + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none") + labs(title="Tumor-Specific Shift of Skeletal Regulators (Tumor vs Matched Normal)", x="Gene", y="Delta Beta (Tumor - Adjacent Normal)")
ggsave("09_Skeletal_Master_Regulators_Cancer_Plot.pdf", plot=gridExtra::grid.arrange(p09a, p09b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D10. 개코원숭이 골관절염(OA) 병적 메틸화 변이와 고형암 종양화 메틸화 변이의 수렴 진화 검증 (수정본)
################################################################################################
cgs_86 <- intersect(COM_D_CpG, rownames(GSE103286m)); mat86_d10 <- as.matrix(GSE103286m[cgs_86, , drop=FALSE]); storage.mode(mat86_d10) <- "numeric"; inf86_d10 <- GSE103286i
idx_h86 <- which(grepl("healthy", inf86_d10$Disease)); idx_o86 <- which(grepl("osteoarthritis", inf86_d10$Disease))
db_oa <- rowMeans(mat86_d10[, idx_o86, drop=FALSE], na.rm=TRUE) - rowMeans(mat86_d10[, idx_h86, drop=FALSE], na.rm=TRUE); names(db_oa) <- cgs_86
stat10 <- data.frame(); plot10_df <- data.frame()
for(cn in tcga_cans){ cand_t <- list(); cand_n <- list()
for(sx in c("F","M")){ nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f")
if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); cand_t[[paste0(sx,"_T")]] <- d[, f=="T", drop=FALSE]; cand_n[[paste0(sx,"_N")]] <- d[, f=="N", drop=FALSE] }
tm_s <- paste0(cn,"_450m_Tm",sx); nm_s <- paste0(cn,"_450m_Nm",sx)
if(length(cand_t)==0 && exists(tm_s)) cand_t[[paste0(sx,"_Tm")]] <- get(tm_s)
if(length(cand_n)==0 && exists(nm_s)) cand_n[[paste0(sx,"_Nm")]] <- get(nm_s) }
if(length(cand_t)==0 && exists(paste0(cn,"_450m_Tm"))) cand_t[["all_T"]] <- get(paste0(cn,"_450m_Tm"))
if(length(cand_n)==0 && exists(paste0(cn,"_450m_Nm"))) cand_n[["all_N"]] <- get(paste0(cn,"_450m_Nm"))
if(length(cand_t)>0 && length(cand_n)>0){ all_m <- c(cand_t, cand_n); com_cg <- names(db_oa)
for(k in 1:length(all_m)){ com_cg <- intersect(com_cg, rownames(all_m[[k]])) }
if(length(com_cg)>=20){ m_tm <- NULL; for(k in 1:length(cand_t)){ sub_m <- as.matrix(cand_t[[k]][com_cg, , drop=FALSE]); m_tm <- if(is.null(m_tm)) sub_m else cbind(m_tm, sub_m) }
m_nm <- NULL; for(k in 1:length(cand_n)){ sub_m <- as.matrix(cand_n[[k]][com_cg, , drop=FALSE]); m_nm <- if(is.null(m_nm)) sub_m else cbind(m_nm, sub_m) }
storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=2 && ncol(m_nm)>=2){ db_cn <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE); sub_oa <- db_oa[com_cg]
ok_cg <- which(is.finite(sub_oa) & is.finite(db_cn) & (abs(sub_oa) > 0.01 | abs(db_cn) > 0.03))
if(length(ok_cg)>=15){ ct <- cor.test(sub_oa[ok_cg], db_cn[ok_cg]); conc_pct <- mean((sub_oa[ok_cg] * db_cn[ok_cg]) > 0) * 100
stat10 <- rbind(stat10, data.frame(Cancer=cn, N_Loci=length(ok_cg), Cor_OA_Cancer_R=ct$estimate, P_val=ct$p.value, Convergence_Pct=conc_pct, Hyper_Shared=sum(sub_oa[ok_cg]>0 & db_cn[ok_cg]>0), Hypo_Shared=sum(sub_oa[ok_cg]<0 & db_cn[ok_cg]<0), stringsAsFactors=FALSE))
sub_10 <- data.frame(Cancer=cn, CpG=names(sub_oa)[ok_cg], Delta_OA=sub_oa[ok_cg], Delta_Cancer=db_cn[ok_cg], Mode=ifelse(sub_oa[ok_cg]*db_cn[ok_cg]>0, "Convergent Alteration", "Divergent"), stringsAsFactors=FALSE)
plot10_df <- rbind(plot10_df, sub_10[seq(1, nrow(sub_10), length.out=min(800, nrow(sub_10))), ]) } } } } }
if(nrow(stat10) > 0){ write.table(stat10, "10_Baboon_OA_vs_Solid_Tumor_Convergence_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p10a <- ggplot(stat10, aes(x=reorder(Cancer, Cor_OA_Cancer_R), y=Cor_OA_Cancer_R, fill=Cor_OA_Cancer_R>0)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Positive Convergence","FALSE"="Inverse Divergence")) + theme_bw() + labs(title="Convergence Between Degenerative Skeletal Disease (Baboon OA) & Cancer Epigenome", x="Solid Cancer", y="Correlation (Pearson R: Delta OA vs Delta Cancer)", fill="Direction")
  p10b <- ggplot(plot10_df, aes(x=Delta_OA, y=Delta_Cancer, color=Mode)) + geom_point(alpha=0.5, size=0.7) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + facet_wrap(~Cancer, scales="free_y", ncol=min(5, length(unique(plot10_df$Cancer)))) + scale_color_manual(values=c("Convergent Alteration"="#E41A1C","Divergent"="#377EB8")) + theme_bw() + theme(legend.position="bottom") + labs(title="Locus-by-Locus Delta Beta: Baboon OA vs TCGA Solid Tumors", x="Delta Beta: Baboon OA - Healthy Bone", y="Delta Beta: Primary Tumor - Matched Normal")
  ggsave("10_Baboon_OA_vs_Solid_Tumor_Convergence_Plot.pdf", plot=gridExtra::grid.arrange(p10a, p10b, ncol=1, heights=c(1, 1.4)), width=13, height=8.5)
} else { cat("경고: 유효한 정상/종양 쌍이 매칭된 암종이 없어 stat10이 비어 있습니다. 워크스페이스 내 객체를 확인하세요.\n") }

################################################################################################
# Module D11. CpG Island-Shore-Open Sea 구조에 따른 골격 유지 좌위 vs 종양 취약 좌위 분업 분석
################################################################################################
if(!exists("GPL13534")){ GPL13534 <- as.data.frame(fread("G:/내 드라이브/Rrunning/Annotations/21_Methyl/GPL13534-11288.txt.gz", header=TRUE, sep="\t"))[, c("ID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")] }
isl_v11 <- as.character(GPL13534$Relation_to_UCSC_CpG_Island[match(COM_D_CpG, GPL13534$ID)]); isl_v11[is.na(isl_v11) | isl_v11==""] <- "Open_Sea"; isl_v11 <- factor(isl_v11, levels=c("Island","N_Shore","S_Shore","N_Shelf","S_Shelf","Open_Sea"))
sd_hb11 <- sqrt(pmax(0, (rowSums(mat_hb^2) - ncol(mat_hb)*(rowMeans(mat_hb)^2))/(ncol(mat_hb)-1)))
can_db_mat11 <- matrix(0, nrow=length(COM_D_CpG), ncol=length(tcga_cans), dimnames=list(COM_D_CpG, tcga_cans))
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[COM_D_CpG, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[COM_D_CpG, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3) can_db_mat11[, cn] <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE) } }
mean_abs_db11 <- rowMeans(abs(can_db_mat11)); df11 <- data.frame(CpG=COM_D_CpG, Island=isl_v11, Bone_SD=sd_hb11, Mean_Cancer_Abs_DeltaBeta=mean_abs_db11)
stat11_isl <- aggregate(cbind(Bone_SD, Mean_Cancer_Abs_DeltaBeta) ~ Island, data=df11, FUN=mean)
stat11_cnt <- as.data.frame(table(df11$Island)); colnames(stat11_cnt) <- c("Island","Total_CpGs"); stat11 <- merge(stat11_cnt, stat11_isl, by="Island")
write.table(stat11, "11_CpG_Island_Partition_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p11a <- ggplot(stat11, aes(x=Island, y=Mean_Cancer_Abs_DeltaBeta, fill=Island)) + geom_bar(stat="identity", color="black", width=0.65) + scale_fill_brewer(palette="Set2") + theme_bw() + theme(legend.position="none") + labs(title="Solid Cancer Aberration Burden Across CpG Island Contexts", x="CpG Island Relation", y="Mean |Delta Beta| Across 10 Cancers")
p11b <- ggplot(df11[seq(1,nrow(df11),length.out=min(6000,nrow(df11))),], aes(x=Island, y=Bone_SD, fill=Island)) + geom_boxplot(alpha=0.8, outlier.size=0.4) + scale_fill_brewer(palette="Set2") + theme_bw() + theme(legend.position="none") + labs(title="Normal Bone Epigenetic Stability by CpG Island Context", x="CpG Island Relation", y="Normal Bone Intersample SD")
ggsave("11_CpG_Island_Partition_Plot.pdf", plot=gridExtra::grid.arrange(p11a, p11b, ncol=2), width=13, height=5.5)

################################################################################################
# Module D12. 영장류 종간 진화 가변성(Evo-SD) 상위 좌위와 고형암 특이적 발암 취약성 중첩 분석
################################################################################################
raw_p12 <- as.matrix(Primate2_Met[COM_D_CpG, ]); storage.mode(raw_p12) <- "numeric"
sp_m12 <- cbind(Mac=rowMeans(raw_p12[, Primate2_Inf$GEO=="GSE103271"], na.rm=TRUE), Bab=rowMeans(raw_p12[, Primate2_Inf$GEO=="GSE103279"], na.rm=TRUE), Ver=rowMeans(raw_p12[, Primate2_Inf$GEO=="GSE103280"], na.rm=TRUE), Mar=rowMeans(raw_p12[, Primate2_Inf$GEO=="GSE103328"], na.rm=TRUE), Hum=rowMeans(mat_hb, na.rm=TRUE))
evo_sd12 <- sqrt(pmax(0, (rowSums(sp_m12^2) - 5*(rowMeans(sp_m12)^2))/4)); cancer_vuln12 <- rowMeans(abs(can_db_mat11))
q_cuts12 <- quantile(evo_sd12, probs=seq(0,1,0.25), na.rm=TRUE); evo_q12 <- cut(evo_sd12, breaks=q_cuts12, labels=c("Q1 (Ultra-Conserved)","Q2 (Low Divergence)","Q3 (Moderate Divergence)","Q4 (Hyper-Divergent)"), include.lowest=TRUE)
df12 <- data.frame(CpG=COM_D_CpG, Evo_SD=evo_sd12, Evo_Quartile=evo_q12, Cancer_Aberration=cancer_vuln12); ct12 <- cor.test(df12$Evo_SD, df12$Cancer_Aberration)
stat12 <- aggregate(cbind(Evo_SD, Cancer_Aberration) ~ Evo_Quartile, data=df12, FUN=mean); write.table(stat12, "12_Evolutionary_SD_vs_Cancer_Aberration_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p12a <- ggplot(df12[!is.na(df12$Evo_Quartile),], aes(x=Evo_Quartile, y=Cancer_Aberration, fill=Evo_Quartile)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + scale_fill_brewer(palette="YlOrRd") + coord_cartesian(ylim=c(0, quantile(df12$Cancer_Aberration, 0.98, na.rm=TRUE))) + theme_bw() + theme(axis.text.x=element_text(angle=25, hjust=1), legend.position="none") + labs(title="Solid Tumor Methylation Plasticity Across Primate Evolutionary Quartiles", x="Interspecific Evolutionary Divergence (Evo-SD)", y="Mean |Delta Beta| Across 10 Solid Cancers")
p12b <- ggplot(df12[seq(1,nrow(df12),length.out=min(5000,nrow(df12))),], aes(x=Evo_SD, y=Cancer_Aberration)) + geom_point(alpha=0.35, size=0.8, color="#0571B0") + geom_smooth(method="lm", color="red", linetype="dashed", se=TRUE) + theme_bw() + labs(title=paste0("Evolutionary Epigenetic Divergence vs Cancer Vulnerability (r = ", round(ct12$estimate, 3), ", P = ", signif(ct12$p.value, 3), ")"), x="Primate Interspecific Evolutionary SD", y="Mean Solid Cancer |Delta Beta|")
ggsave("12_Evolutionary_SD_vs_Cancer_Aberration_Plot.pdf", plot=gridExtra::grid.arrange(p12a, p12b, ncol=2), width=13, height=5.5)

################################################################################################
# Module D13. 정상 뼈 vs 고형암 간 염색체 1~22번 상염색체별 메틸화 드리프트 편차 비교 (수정본)
################################################################################################
chr_v13 <- as.character(GPL13534$CHR[match(COM_D_CpG, GPL13534$ID)]); df13 <- data.frame(CpG=COM_D_CpG, CHR=chr_v13, Bone_SD=sd_hb11, Cancer_Shift=cancer_vuln12)
df13 <- df13[df13$CHR %in% as.character(1:22), ]; df13$CHR <- factor(df13$CHR, levels=as.character(1:22)); stat13 <- data.frame()
for(ch in levels(df13$CHR)){ sub_ch <- df13[df13$CHR==ch, ]
if(nrow(sub_ch)>=10){ ct <- cor.test(sub_ch$Bone_SD, sub_ch$Cancer_Shift)
stat13 <- rbind(stat13, data.frame(CHR=ch, N_CpGs=nrow(sub_ch), Mean_Bone_SD=mean(sub_ch$Bone_SD, na.rm=TRUE), Mean_Cancer_Shift=mean(sub_ch$Cancer_Shift, na.rm=TRUE), Bone_Cancer_Cor_R=ct$estimate, P_val=ct$p.value)) } }
stat13$CHR <- factor(stat13$CHR, levels=as.character(1:22)); write.table(stat13, "13_Chromosomal_Epigenetic_Drift_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p13a <- ggplot(stat13, aes(x=CHR, y=Mean_Cancer_Shift, fill=Mean_Bone_SD)) + geom_bar(stat="identity", color="black", width=0.7) + scale_fill_viridis_c(option="plasma") + theme_bw() + labs(title="Autosomal Cancer Aberration Burden Across Chromosomes 1-22", x="Human Autosome", y="Chromosomal Mean Cancer |Delta Beta|", fill="Bone SD")
p13b <- ggplot(stat13, aes(x=Mean_Bone_SD, y=Mean_Cancer_Shift)) + geom_point(size=3.5, color="#D95F02") + ggrepel::geom_text_repel(aes(label=paste0("Chr", CHR)), size=3.2) + geom_smooth(method="lm", formula=y~x, se=FALSE, linetype="dashed", color="black") + theme_bw() + labs(title="Chromosomal Epigenetic Drift: Normal Bone Variation vs Cancer Disturbance", x="Chromosomal Mean Normal Bone SD", y="Chromosomal Mean Cancer |Delta Beta|")
ggsave("13_Chromosomal_Epigenetic_Drift_Plot.pdf", plot=gridExtra::grid.arrange(p13a, p13b, ncol=2), width=13, height=5.5)

################################################################################################
# Module D14. 고형암 전이 억제 유전자군(Metastasis Suppressors)의 정상 뼈 메틸화 보호 패턴 프로파일링
################################################################################################
ms_genes <- c("CDH1","NME1","PEBP1","BRMS1","CD82","KISS1","TIMP3","CASP8","PTEN","TP53")
ms_anno <- GPL13534[GPL13534$ID %in% COM_D_CpG & grepl(paste(ms_genes, collapse="|"), as.character(GPL13534$UCSC_RefGene_Name)) & grepl("TSS1500|TSS200|1stExon", as.character(GPL13534$UCSC_RefGene_Group)), ]
ms_anno$Gene <- sub(";.*", "", as.character(ms_anno$UCSC_RefGene_Name)); ms_cgs <- unique(ms_anno$ID); stat14 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[ms_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[ms_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ for(cg in ms_cgs){ gn <- ms_anno$Gene[ms_anno$ID==cg][1]; v_hb <- mat_hb[cg, ]; v_tm <- m_tm[cg, ]; v_nm <- m_nm[cg, ]
stat14 <- rbind(stat14, data.frame(Cancer=cn, CpG=cg, Gene=gn, Island=as.character(ms_anno$Relation_to_UCSC_CpG_Island[ms_anno$ID==cg][1]), Mean_Bone=mean(v_hb,na.rm=TRUE), Mean_Normal=mean(v_nm,na.rm=TRUE), Mean_Tumor=mean(v_tm,na.rm=TRUE), Delta_Tumor_vs_Bone=mean(v_tm,na.rm=TRUE)-mean(v_hb,na.rm=TRUE), Delta_Tumor_vs_Normal=mean(v_tm,na.rm=TRUE)-mean(v_nm,na.rm=TRUE))) } } } }
write.table(stat14, "14_Metastasis_Suppressor_Promoter_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
ms_summary14 <- aggregate(cbind(Mean_Bone, Delta_Tumor_vs_Bone, Delta_Tumor_vs_Normal) ~ Cancer + Gene, data=stat14, FUN=mean)
p14a <- ggplot(ms_summary14, aes(x=Cancer, y=Gene, fill=Delta_Tumor_vs_Bone)) + geom_tile(color="white") + geom_text(aes(label=round(Delta_Tumor_vs_Bone,2)), size=3) + scale_fill_gradient2(low="#377EB8", mid="white", high="#E41A1C", midpoint=0) + theme_minimal() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Promoter Methylation Shift of Metastasis Suppressors (Tumor vs Healthy Bone)", subtitle="Red indicates Tumor Hypermethylation/Silencing; Blue indicates Bone Protection", x="Solid Cancer", y="Metastasis Suppressor Gene", fill="Delta Beta")
p14b <- ggplot(stat14, aes(x=Gene, y=Mean_Bone, fill=Gene)) + geom_boxplot(alpha=0.8, outlier.size=0.5) + geom_hline(yintercept=0.2, linetype="dashed", color="blue") + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none") + labs(title="Constitutive Hypomethylation (Protection) of Metastasis Suppressors in Healthy Bone", x="Gene", y="Healthy Bone Methylation Beta Value")
ggsave("14_Metastasis_Suppressor_Promoter_Plot.pdf", plot=gridExtra::grid.arrange(p14a, p14b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D15. 정상 뼈의 중간 메틸화 분율(IMF, beta in [0.3,0.7]) 쌍안정성과 종양 내 후성유전적 이질성(ITH) 연결
################################################################################################
bone_imf_freq15 <- rowMeans(mat_hb >= 0.3 & mat_hb <= 0.7, na.rm=TRUE)
imf_cat15 <- cut(bone_imf_freq15, breaks=c(-0.01, 0.1, 0.5, 1.01), labels=c("Fixed Un/Hyper (<10% IMF)","Dynamic IMF (10-50%)","Constitutive IMF (>50%)"))
ith_res15 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm")
if(exists(tm_n)){ m_tm <- as.matrix(get(tm_n)[COM_D_CpG, , drop=FALSE]); storage.mode(m_tm) <- "numeric"
if(ncol(m_tm)>=10){ t_var <- pmax(0, (rowSums(m_tm^2) - ncol(m_tm)*(rowMeans(m_tm)^2))/(ncol(m_tm)-1))
sub_ith <- data.frame(CpG=COM_D_CpG, Cancer=cn, Bone_IMF=bone_imf_freq15, IMF_Class=imf_cat15, Tumor_Variance=t_var)
ith_res15 <- rbind(ith_res15, sub_ith) } } }
stat15 <- aggregate(Tumor_Variance ~ Cancer + IMF_Class, data=ith_res15, FUN=mean)
write.table(stat15, "15_Bone_IMF_vs_Tumor_ITH_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p15a <- ggplot(ith_res15[!is.na(ith_res15$IMF_Class), ], aes(x=IMF_Class, y=Tumor_Variance, fill=IMF_Class)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + facet_wrap(~Cancer, scales="free_y", ncol=5) + scale_fill_brewer(palette="Set1") + coord_cartesian(ylim=c(0, 0.08)) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none") + labs(title="Tumor Epigenetic Heterogeneity (Variance) Stratified by Normal Bone IMF Plasticity", x="Healthy Bone IMF State", y="Intratumoral Epigenetic Variance (ITH)")
p15b <- ggplot(stat15, aes(x=Cancer, y=Tumor_Variance, fill=IMF_Class)) + geom_bar(stat="identity", position="dodge", width=0.7, color="black") + scale_fill_brewer(palette="Set1") + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Mean Tumor Epigenetic Variance Across Cancers by Bone IMF Category", x="Solid Cancer", y="Mean Tumor Methylation Variance", fill="Bone IMF Class")
ggsave("15_Bone_IMF_vs_Tumor_ITH_Plot.pdf", plot=gridExtra::grid.arrange(p15a, p15b, ncol=1, heights=c(1.3, 1)), width=13, height=8.5)

################################################################################################
# Module D16. 뼈-고형암 통합 다차원 차원 축소(MDS & PCA)를 통한 조직 분화 대 악성화 궤적 분리 (수정본)
################################################################################################
cand_mats16 <- list(Healthy_Bone=mat_hb); cand_meta16 <- data.frame(Sample=colnames(mat_hb), Group="Healthy_Bone", Cancer="Bone", stringsAsFactors=FALSE)
for(cn in tcga_cans){ for(sx in c("F","M")){ nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f")
if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); s_t <- head(colnames(d)[f=="T"], 10); s_n <- head(colnames(d)[f=="N"], 3)
if(length(s_t)>0){ cand_mats16[[paste0(cn,"_",sx,"_T")]] <- d[, s_t, drop=FALSE]; cand_meta16 <- rbind(cand_meta16, data.frame(Sample=s_t, Group="Tumor", Cancer=cn, stringsAsFactors=FALSE)) }
if(length(s_n)>0){ cand_mats16[[paste0(cn,"_",sx,"_N")]] <- d[, s_n, drop=FALSE]; cand_meta16 <- rbind(cand_meta16, data.frame(Sample=s_n, Group="Normal", Cancer=cn, stringsAsFactors=FALSE)) } }
else { tm_s <- paste0(cn,"_450m_Tm",sx); nm_s <- paste0(cn,"_450m_Nm",sx)
if(exists(tm_s)){ d_t <- get(tm_s); s_t <- head(colnames(d_t), 10); cand_mats16[[paste0(cn,"_",sx,"_Tm")]] <- d_t[, s_t, drop=FALSE]; cand_meta16 <- rbind(cand_meta16, data.frame(Sample=s_t, Group="Tumor", Cancer=cn, stringsAsFactors=FALSE)) }
if(exists(nm_s)){ d_n <- get(nm_s); s_n <- head(colnames(d_n), 3); cand_mats16[[paste0(cn,"_",sx,"_Nm")]] <- d_n[, s_n, drop=FALSE]; cand_meta16 <- rbind(cand_meta16, data.frame(Sample=s_n, Group="Normal", Cancer=cn, stringsAsFactors=FALSE)) } } } }
com_cg16 <- COM_D_CpG; for(k in 1:length(cand_mats16)){ com_cg16 <- intersect(com_cg16, rownames(cand_mats16[[k]])) }
m_comb16 <- NULL; for(k in 1:length(cand_mats16)){ sub_mat <- as.matrix(cand_mats16[[k]][com_cg16, , drop=FALSE]); storage.mode(sub_mat) <- "numeric"; m_comb16 <- if(is.null(m_comb16)) sub_mat else cbind(m_comb16, sub_mat) }
m_comb16 <- m_comb16[complete.cases(m_comb16), , drop=FALSE]
if(nrow(m_comb16)>=20 && ncol(m_comb16)>=5){ pca16 <- prcomp(t(m_comb16), center=TRUE, scale.=TRUE); mds16 <- cmdscale(dist(t(m_comb16)), k=2)
df16 <- data.frame(Sample=colnames(m_comb16), Group=cand_meta16$Group[match(colnames(m_comb16), cand_meta16$Sample)], Cancer=cand_meta16$Cancer[match(colnames(m_comb16), cand_meta16$Sample)], PC1=pca16$x[,1], PC2=pca16$x[,2], MDS1=mds16[,1], MDS2=mds16[,2], PC1_Var=round(summary(pca16)$importance[2,1]*100,1), PC2_Var=round(summary(pca16)$importance[2,2]*100,1), stringsAsFactors=FALSE)
write.table(df16, "16_Bone_TCGA_MDS_PCA_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p16a <- ggplot(df16, aes(x=PC1, y=PC2, color=Group, shape=Cancer)) + geom_point(size=2.8, alpha=0.85) + scale_color_manual(values=c("Healthy_Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + scale_shape_manual(values=c(16,17,15,3,7,8,9,10,12,13,14,1,2,4,5)[1:length(unique(df16$Cancer))]) + theme_bw() + labs(title=paste0("Integrated PCA: Normal Bone vs Solid Cancers (PC1: ",df16$PC1_Var[1],"%, PC2: ",df16$PC2_Var[1],"%)"), x="PC1 (Lineage / Tissue Axis)", y="PC2 (Malignancy Trajectory Axis)")
p16b <- ggplot(df16, aes(x=MDS1, y=MDS2, color=Group, shape=Cancer)) + geom_point(size=2.8, alpha=0.85) + scale_color_manual(values=c("Healthy_Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + scale_shape_manual(values=c(16,17,15,3,7,8,9,10,12,13,14,1,2,4,5)[1:length(unique(df16$Cancer))]) + theme_bw() + labs(title="Multidimensional Scaling (MDS) Map of Bone vs Solid Cancers", x="MDS Dimension 1", y="MDS Dimension 2")
ggsave("16_Bone_TCGA_MDS_PCA_Plot.pdf", plot=gridExtra::grid.arrange(p16a, p16b, ncol=2), width=14, height=6)
} else { cat("경고: 유효 공통 프로브 또는 검체 수가 부족하여 PCA/MDS를 수행할 수 없습니다 (행:", nrow(m_comb16), ", 열:", ncol(m_comb16), ")\n") }

################################################################################################
# Module D17. 골격계 공메틸화 네트워크(WGCNA) 모듈의 TCGA 고형암 투영 및 암종별 모듈 보존성 검정
################################################################################################
var_hb17 <- rowSums(mat_hb^2) - ncol(mat_hb)*(rowMeans(mat_hb)^2); top_v17 <- COM_D_CpG[order(var_hb17, decreasing=TRUE)[1:1200]]
c_mat17 <- cor(t(mat_hb[top_v17, ])); c_mat17[!is.finite(c_mat17)] <- 0; dis17 <- as.dist(1 - abs(c_mat17)); hc17 <- hclust(dis17, method="average"); mod_ids17 <- cutree(hc17, k=5); names(mod_ids17) <- top_v17
write.table(data.frame(CpG=top_v17, Bone_Module=paste0("Mod_", mod_ids17)), "17_Bone_WGCNA_Module_Assignments.txt", sep="\t", quote=FALSE, row.names=FALSE)
stat17 <- data.frame(); for(m in 1:5){ m_cgs <- names(mod_ids17)[mod_ids17 == m]; m_score_bone <- mean(mat_hb[m_cgs, ], na.rm=TRUE)
stat17 <- rbind(stat17, data.frame(Cancer="Healthy_Bone", Module=paste0("Mod_",m), N_CpGs=length(m_cgs), State="Bone", Mean_Score=m_score_bone))
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n)){ m_tm <- as.matrix(get(tm_n)[m_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; stat17 <- rbind(stat17, data.frame(Cancer=cn, Module=paste0("Mod_",m), N_CpGs=length(m_cgs), State="Tumor", Mean_Score=mean(m_tm, na.rm=TRUE))) }
if(exists(nm_n)){ m_nm <- as.matrix(get(nm_n)[m_cgs, , drop=FALSE]); storage.mode(m_nm) <- "numeric"; stat17 <- rbind(stat17, data.frame(Cancer=cn, Module=paste0("Mod_",m), N_CpGs=length(m_cgs), State="Normal", Mean_Score=mean(m_nm, na.rm=TRUE))) } } }
write.table(stat17, "17_Bone_Module_TCGA_Projection_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
stat17_w <- reshape(stat17[stat17$State %in% c("Tumor","Normal"), c("Cancer","Module","State","Mean_Score")], idvar=c("Cancer","Module"), timevar="State", direction="wide"); stat17_w$Delta_T_minus_N <- stat17_w$Mean_Score.Tumor - stat17_w$Mean_Score.Normal
p17a <- ggplot(stat17, aes(x=Module, y=Mean_Score, fill=State)) + geom_boxplot(alpha=0.8, outlier.size=0.5) + scale_fill_manual(values=c("Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + theme_bw() + labs(title="5 Bone Co-Methylation Modules Projected Across Solid Tissues", x="Bone Co-Methylation Module", y="Mean Module Methylation Beta Value")
p17b <- ggplot(stat17_w, aes(x=Cancer, y=Module, fill=Delta_T_minus_N)) + geom_tile(color="white") + geom_text(aes(label=round(Delta_T_minus_N, 2)), size=3.2) + scale_fill_gradient2(low="#377EB8", mid="white", high="#E41A1C", midpoint=0) + theme_minimal() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Tumor-Specific Deregulation of Bone Modules (Tumor - Matched Normal)", x="Solid Cancer", y="Bone Co-Methylation Module", fill="Delta Beta")
ggsave("17_Bone_WGCNA_TCGA_Projection_Plot.pdf", plot=gridExtra::grid.arrange(p17a, p17b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D18. 정상 뼈 vs 고형암 구분을 위한 최소 50개 핵심 후성유전 패널 선별(LASSO Penalized Regression)
################################################################################################
n_hb <- ncol(mat_hb); s_tm_pool <- c()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); if(exists(tm_n)){ s_cand <- intersect(colnames(get(tm_n)), meta_d01$Sample[meta_d01$Condition=="Tumor"]); s_tm_pool <- c(s_tm_pool, head(s_cand, 8)) } }
x_mat18 <- matrix(NA, nrow=length(COM_D_CpG), ncol=n_hb + length(s_tm_pool), dimnames=list(COM_D_CpG, c(colnames(mat_hb), s_tm_pool)))
x_mat18[, 1:n_hb] <- mat_hb
for(s in s_tm_pool){ for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); if(exists(tm_n) && s %in% colnames(get(tm_n))) x_mat18[, s] <- get(tm_n)[COM_D_CpG, s] } }
y_vec18 <- c(rep(0, n_hb), rep(1, length(s_tm_pool))); x_mat18 <- x_mat18[complete.cases(x_mat18), ]
set.seed(312); cv_lasso18 <- glmnet::cv.glmnet(x=t(x_mat18), y=y_vec18, family="binomial", alpha=1, nfolds=5)
cf18 <- coef(cv_lasso18, s="lambda.min"); cf_v18 <- as.numeric(cf18); names(cf_v18) <- rownames(cf18); cf_sub18 <- cf_v18[names(cf_v18) != "(Intercept)" & cf_v18 != 0]
top50_lasso18 <- names(sort(abs(cf_sub18), decreasing=TRUE))[1:min(50, length(cf_sub18))]
df18_cpg <- data.frame(CpG=top50_lasso18, Coef=cf_sub18[top50_lasso18], Odds_Ratio=exp(cf_sub18[top50_lasso18]), Direction=ifelse(cf_sub18[top50_lasso18]>0, "Tumor-Enriched (+)", "Bone-Enriched (-)"))
write.table(df18_cpg, "18_Bone_vs_Cancer_LASSO_Top50_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
cv_df18 <- data.frame(lambda=cv_lasso18$lambda, cvm=cv_lasso18$cvm, cvup=cv_lasso18$cvup, cvlo=cv_lasso18$cvlo)
p18a <- ggplot(cv_df18, aes(x=log(lambda), y=cvm)) + geom_point(color="red", size=1.5) + geom_errorbar(aes(ymin=cvlo, ymax=cvup), width=0.08, color="gray60") + geom_vline(xintercept=log(cv_lasso18$lambda.min), linetype="dashed", color="blue") + theme_bw() + labs(title="LASSO Cross-Validation: Distinguishing Healthy Bone from Solid Cancers", subtitle=paste0("Lambda Min: ", round(cv_lasso18$lambda.min, 4), " (Selected Loci: ", length(cf_sub18), ")"), x="log(Lambda)", y="Binomial Deviance")
p18b <- ggplot(head(df18_cpg, 20), aes(x=reorder(CpG, abs(Coef)), y=Coef, fill=Direction)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + scale_fill_manual(values=c("Tumor-Enriched (+)"="#E41A1C","Bone-Enriched (-)"="#1B9E77")) + theme_bw() + labs(title="Top 20 Discriminative LASSO Features", x="CpG Site", y="LASSO Penalized Coefficient")
ggsave("18_Bone_vs_Cancer_LASSO_Plot.pdf", plot=gridExtra::grid.arrange(p18a, p18b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D19. 잠재 공간 도메인 적응(Domain Adaptation PCA)을 통한 뼈-고형암 간 조직 불변 드라이버 좌위 추출
################################################################################################
pca_hb19 <- prcomp(t(mat_hb), center=TRUE, scale.=TRUE); rot_sub19 <- pca_hb19$rotation[, 1:3]; ctr19 <- pca_hb19$center; scl19 <- pca_hb19$scale
load_sum19 <- rowSums(abs(rot_sub19)); top_dr19 <- names(sort(load_sum19, decreasing=TRUE)[1:50])
g_dr19 <- sub(";.*", "", as.character(GPL13534$UCSC_RefGene_Name[match(top_dr19, GPL13534$ID)]))
write.table(data.frame(CpG=top_dr19, Gene=g_dr19, Loading_Sum_PC1_3=load_sum19[top_dr19]), "19_Domain_Adaptation_Driver_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
proj_df19 <- data.frame(Sample=colnames(mat_hb), Group="Healthy_Bone", Cancer="Bone", Latent_PC1=pca_hb19$x[,1], Latent_PC2=pca_hb19$x[,2])
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm")
if(exists(tm_n)){ m_tm <- as.matrix(get(tm_n)[COM_D_CpG, , drop=FALSE]); storage.mode(m_tm) <- "numeric"
s_cand <- intersect(colnames(m_tm), meta_d01$Sample[meta_d01$Condition=="Tumor"]); s_cand <- head(s_cand, 20)
if(length(s_cand)>=5){ m_scaled <- t((m_tm[, s_cand] - ctr19) / scl19); m_proj <- m_scaled %*% rot_sub19[, 1:2]
proj_df19 <- rbind(proj_df19, data.frame(Sample=s_cand, Group="Tumor", Cancer=cn, Latent_PC1=m_proj[,1], Latent_PC2=m_proj[,2])) } } }
write.table(proj_df19, "19_Domain_Adaptation_Sample_Projections.txt", sep="\t", quote=FALSE, row.names=FALSE)
p19a <- ggplot(proj_df19, aes(x=Latent_PC1, y=Latent_PC2, color=Group, shape=Cancer)) + geom_point(size=2.6, alpha=0.8) + scale_color_manual(values=c("Healthy_Bone"="#1B9E77","Tumor"="#E41A1C")) + scale_shape_manual(values=c(16,17,15,3,7,8,9,10,12,13,14)) + theme_bw() + labs(title="Domain Adaptation: Solid Cancers Projected onto Bone Latent Space", x="Bone Latent PC1", y="Bone Latent PC2")
df_dr_plot <- data.frame(Feature=paste0(g_dr19[1:20], " (", top_dr19[1:20], ")"), Loading=load_sum19[top_dr19[1:20]])
p19b <- ggplot(df_dr_plot, aes(x=reorder(Feature, Loading), y=Loading)) + geom_bar(stat="identity", fill="#8F87F1", color="black", width=0.7) + coord_flip() + theme_bw() + labs(title="Top 20 Latent Invariant Drivers", x="Gene (CpG)", y="Loading Sum (|PC1| + |PC2| + |PC3|)")
ggsave("19_Domain_Adaptation_Driver_Plot.pdf", plot=gridExtra::grid.arrange(p19a, p19b, ncol=2), width=14, height=5.5)

################################################################################################
# Module D20. 고형암 환자 생존율 및 사망 위험(Vital Status)과 골 모방성 점수(Bone Mimicry Score) 간 로지스틱 분석 (완전수정본)
################################################################################################
mean_hb20 <- rowMeans(mat_hb, na.rm=TRUE); names(mean_hb20) <- COM_D_CpG; stat20 <- data.frame(); surv_plot_df <- data.frame()
for(cn in tcga_cans){ cand_m <- list(); cand_c <- list()
for(sx in c("F","M")){ m1_n <- paste0(cn,"_450m_Tm",sx,"1"); c1_n <- paste0(cn,"_c450m",sx,"1")
if(exists(m1_n) && exists(c1_n)){ cand_m[[paste0(sx,"1")]] <- get(m1_n); cand_c[[paste0(sx,"1")]] <- get(c1_n) }
else { nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f"); cl_n <- paste0(cn,"_c")
if(exists(nt_n) && exists(nt_f) && exists(cl_n)){ d_nt <- get(nt_n); f_nt <- get(nt_f); cand_m[[sx]] <- d_nt[, f_nt=="T", drop=FALSE]; cand_c[[sx]] <- get(cl_n) } } }
if(length(cand_m)==0 && exists(paste0(cn,"_450m_Tm")) && exists(paste0(cn,"_c"))){ cand_m[["all"]] <- get(paste0(cn,"_450m_Tm")); cand_c[["all"]] <- get(paste0(cn,"_c")) }
if(length(cand_m)>0){ m_pool <- NULL; c_pool <- data.frame()
for(k in 1:length(cand_m)){ mat_k <- cand_m[[k]]; cln_k <- cand_c[[k]]; cg_k <- intersect(names(mean_hb20), rownames(mat_k))
pat_m <- toupper(substr(gsub("[_\\.]", "-", colnames(mat_k)), 1, 12)); pat_c <- toupper(substr(gsub("[_\\.]", "-", cln_k$ID), 1, 12)); com_p <- intersect(pat_m, pat_c)
if(length(com_p)>=5 && length(cg_k)>=20){ sub_mk <- as.matrix(mat_k[cg_k, match(com_p, pat_m), drop=FALSE]); storage.mode(sub_mk) <- "numeric"; sub_ck <- cln_k[match(com_p, pat_c), , drop=FALSE]
m_pool <- if(is.null(m_pool)) sub_mk else { c_com <- intersect(rownames(m_pool), rownames(sub_mk)); cbind(m_pool[c_com, , drop=FALSE], sub_mk[c_com, , drop=FALSE]) }
c_pool <- rbind(c_pool, sub_ck[, c("ID","Age","Sex","Surv")]) } }
if(!is.null(m_pool) && nrow(c_pool)>=15){ cg_use <- intersect(names(mean_hb20), rownames(m_pool)); mat_sub <- m_pool[cg_use, , drop=FALSE]
for(r in 1:nrow(mat_sub)){ na_idx <- which(is.na(mat_sub[r, ])); if(length(na_idx)>0){ rm_v <- mean(mat_sub[r, -na_idx], na.rm=TRUE); mat_sub[r, na_idx] <- ifelse(is.finite(rm_v), rm_v, 0.5) } }
r_bone_sc <- suppressWarnings(as.numeric(cor(mean_hb20[cg_use], mat_sub)))
surv_ch <- as.character(c_pool$Surv); ev_num <- ifelse(surv_ch=="D", 1, ifelse(surv_ch=="A", 0, NA)); ag_num <- suppressWarnings(as.numeric(as.character(c_pool$Age))); sx_ch <- toupper(substr(as.character(c_pool$Sex), 1, 1))
d_fit <- data.frame(Event=ev_num, Bone_Mimicry=as.numeric(scale(r_bone_sc)), Age=ag_num, Sex=sx_ch, stringsAsFactors=FALSE)
ok_fit <- which(is.finite(d_fit$Event) & is.finite(d_fit$Bone_Mimicry) & is.finite(d_fit$Age) & d_fit$Sex %in% c("F","M")); d_sub <- d_fit[ok_fit, ]
n_d <- sum(d_sub$Event==1); n_a <- sum(d_sub$Event==0)
if(n_d>=3 && n_a>=3 && nrow(d_sub)>=15 && is.finite(var(d_sub$Bone_Mimicry)) && var(d_sub$Bone_Mimicry)>0){
  fml <- as.formula(if(length(unique(d_sub$Sex))>1) "Event ~ Bone_Mimicry + Age + Sex" else "Event ~ Bone_Mimicry + Age")
  glm_fit <- suppressWarnings(glm(fml, data=d_sub, family=binomial)); cf_tab <- summary(glm_fit)$coefficients
  if("Bone_Mimicry" %in% rownames(cf_tab) && is.finite(cf_tab["Bone_Mimicry", 1]) && is.finite(cf_tab["Bone_Mimicry", 2])){ coef_bm <- cf_tab["Bone_Mimicry", ]
  or_v <- exp(coef_bm[1]); ci_l <- pmax(0.01, exp(coef_bm[1]-1.96*coef_bm[2])); ci_u <- pmin(100, exp(coef_bm[1]+1.96*coef_bm[2]))
  stat20 <- rbind(stat20, data.frame(Cancer=cn, N=nrow(d_sub), Events=n_d, Log_OR=coef_bm[1], SE=coef_bm[2], OR=or_v, CI_Lower=ci_l, CI_Upper=ci_u, P_value=coef_bm[4], stringsAsFactors=FALSE))
  surv_plot_df <- rbind(surv_plot_df, data.frame(Cancer=cn, Status=ifelse(d_sub$Event==1,"Deceased","Alive"), Bone_Mimicry=d_sub$Bone_Mimicry, stringsAsFactors=FALSE)) } } } } }
if(nrow(stat20)>0){ write.table(stat20, "20_TCGA_Survival_Bone_Mimicry_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p20a <- ggplot(stat20, aes(x=reorder(Cancer, OR), y=OR, color=P_value<0.05)) + geom_point(size=3) + geom_errorbar(aes(ymin=CI_Lower, ymax=CI_Upper), width=0.2, size=0.9) + geom_hline(yintercept=1, linetype="dashed", color="red") + coord_flip() + scale_color_manual(values=c("TRUE"="#E41A1C","FALSE"="gray50")) + theme_bw() + labs(title="Mortality Risk by Epigenetic Bone Mimicry Score (Adjusted for Age & Sex)", subtitle="Odds Ratio > 1 indicates bone-like methylation profile correlates with increased death risk", x="Solid Cancer", y="Odds Ratio (95% CI)")
  p20b <- ggplot(surv_plot_df, aes(x=Cancer, y=Bone_Mimicry, fill=Status)) + geom_boxplot(alpha=0.8, outlier.size=0.4) + scale_fill_manual(values=c("Alive"="#377EB8","Deceased"="#E41A1C")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Epigenetic Bone Mimicry Score Distribution by Vital Status", x="Solid Cancer", y="Standardized Bone Mimicry Score (Z)")
  ggsave("20_TCGA_Survival_Bone_Mimicry_Plot.pdf", plot=gridExtra::grid.arrange(p20a, p20b, ncol=2), width=14, height=6)
} else { cat("경고: 유효한 사망/생존 표본이 매칭된 암종이 없어 stat20이 비어 있습니다. 임상 객체(*_c450mF1/M1 등)를 확인하세요.\n") }

################################################################################################
# Module D21. 대퇴골 29개 형태학 지표 연관 핵심 후성유전 좌위의 10대 고형암 이상 메틸화(Delta Beta) 스크리닝
################################################################################################
mat_pr21 <- as.matrix(Primate2_Met_Cor[COM_D_CpG, ]); storage.mode(mat_pr21) <- "numeric"; inf_p21 <- Primate2_Inf
z_29_21 <- matrix(0, nrow(inf_p21), 29, dimnames=list(rownames(inf_p21), colnames(inf_p21)[3:31]))
for(g in unique(inf_p21$GEO)){ idx <- which(inf_p21$GEO==g); for(j in 1:29){ v <- inf_p21[idx, j+2]; sv <- as.numeric(scale(v)); sv[!is.finite(sv)] <- 0; z_29_21[idx, j] <- sv } }
r_m29_21 <- abs(cor(t(mat_pr21), z_29_21, use="pairwise.complete.obs")); r_m29_21[!is.finite(r_m29_21)] <- 0; top50_morph_cgs <- COM_D_CpG[order(rowMeans(r_m29_21), decreasing=TRUE)[1:50]]
stat21 <- data.frame(); for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[top50_morph_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[top50_morph_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE)
for(cg in top50_morph_cgs){ pv <- t.test(m_tm[cg, ], m_nm[cg, ])$p.value; stat21 <- rbind(stat21, data.frame(Cancer=cn, CpG=cg, Morpho_AbsR=rowMeans(r_m29_21)[cg], Delta_Beta=db[cg], P_value=pv, Sig=ifelse(pv<0.05,"*",""))) } } } }
write.table(stat21, "21_Femur_Morphology_CpG_Cancer_Shift_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
top20_plot_cgs <- head(top50_morph_cgs, 20); sub_stat21 <- stat21[stat21$CpG %in% top20_plot_cgs, ]
p21a <- ggplot(sub_stat21, aes(x=Cancer, y=CpG, fill=Delta_Beta)) + geom_tile(color="white") + geom_text(aes(label=Sig), size=4, vjust=0.75) + scale_fill_gradient2(low="#377EB8", mid="white", high="#E41A1C", midpoint=0) + theme_minimal() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Tumor-Specific Aberration of Top 20 Femur Morphology-Coupled CpGs (* P < 0.05)", x="Solid Cancer", y="Morphology-Coupled CpG", fill="Delta Beta")
p21b <- ggplot(stat21, aes(x=Cancer, y=abs(Delta_Beta), fill=Cancer)) + geom_boxplot(alpha=0.8, outlier.size=0.4) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none") + labs(title="Overall Tumor Dysregulation Burden Across 50 Skeletal Morpho-CpGs", x="Solid Cancer", y="|Delta Beta| (Tumor - Matched Normal)")
ggsave("21_Femur_Morphology_CpG_Cancer_Shift_Plot.pdf", plot=gridExtra::grid.arrange(p21a, p21b, ncol=2), width=14, height=6)

################################################################################################
# Module D22. 대퇴골 단면 관성모멘트(J, 극성 면적 모멘트) 강성 조절 CpG와 고형암 기질 교란 비교
################################################################################################
inf_p22 <- Primate2_Inf; L_m22 <- (inf_p22$FL + inf_p22$BFL + inf_p22$MFL)/3; J_mid22 <- (pi/64)*(inf_p22$MSWd*(inf_p22$MSDp^3) + inf_p22$MSDp*(inf_p22$MSWd^3))/(L_m22^4)*1e8; log_j22 <- log10(J_mid22)
z_j22 <- numeric(nrow(inf_p22)); for(g in unique(inf_p22$GEO)){ idx <- which(inf_p22$GEO==g); z_j22[idx] <- as.numeric(scale(log_j22[idx])) }
r_j22 <- suppressWarnings(as.numeric(cor(z_j22, t(mat_pr21)))); names(r_j22) <- COM_D_CpG; top30_j_cgs <- COM_D_CpG[order(abs(r_j22), decreasing=TRUE)[1:30]]
stat22 <- data.frame(); for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[top30_j_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[top30_j_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE)
stat22 <- rbind(stat22, data.frame(Cancer=cn, Mean_Rigidity_CpG_Shift=mean(abs(db)), Mean_Signed_Shift=mean(db), Hyper_CpGs=sum(db>0.05), Hypo_CpGs=sum(db < -0.05))) } } }
write.table(stat22, "22_Femur_Rigidity_J_Cancer_Shift_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
j_plot_df22 <- data.frame(CpG=rep(top30_j_cgs, length(tcga_cans)), Cancer=rep(tcga_cans, each=length(top30_j_cgs)), Bone_J_R=rep(r_j22[top30_j_cgs], length(tcga_cans)), Delta_Beta=NA)
for(i in 1:nrow(j_plot_df22)){ tm_n <- paste0(j_plot_df22$Cancer[i],"_450m_Tm"); nm_n <- paste0(j_plot_df22$Cancer[i],"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- get(tm_n); m_nm <- get(nm_n); cg <- j_plot_df22$CpG[i]; if(cg %in% rownames(m_tm) && cg %in% rownames(m_nm)) j_plot_df22$Delta_Beta[i] <- mean(as.numeric(m_tm[cg, ]), na.rm=TRUE) - mean(as.numeric(m_nm[cg, ]), na.rm=TRUE) } }
j_plot_df22 <- j_plot_df22[!is.na(j_plot_df22$Delta_Beta), ]
p22a <- ggplot(stat22, aes(x=reorder(Cancer, Mean_Rigidity_CpG_Shift), y=Mean_Rigidity_CpG_Shift, fill=Mean_Signed_Shift>0)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Net Hypermethylated","FALSE"="Net Hypomethylated")) + theme_bw() + labs(title="Tumor Disturbance of Diaphyseal Rigidity (J) Regulators", x="Solid Cancer", y="Mean |Delta Beta| across 30 Rigidity CpGs", fill="Net Shift")
p22b <- ggplot(j_plot_df22, aes(x=Bone_J_R, y=Delta_Beta, color=Cancer)) + geom_point(alpha=0.6, size=1.5) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + theme_bw() + labs(title="Bone Mechanical Rigidity (J) Correlation vs Cancer Delta Beta", x="Primate Bone Rigidity (J) Pearson R", y="Tumor Delta Beta (Tumor - Matched Normal)")
ggsave("22_Femur_Rigidity_J_Cancer_Shift_Plot.pdf", plot=gridExtra::grid.arrange(p22a, p22b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D23. 대퇴골 29차원 형태 가소성(Morphological Decanalization) 좌위와 고형암 진행도 연관성 (수정본)
################################################################################################
morph_decan23 <- sqrt(rowSums(z_29_21^2) / 29); r_decan23 <- suppressWarnings(as.numeric(cor(morph_decan23, t(mat_pr21))))
names(r_decan23) <- COM_D_CpG; top40_decan_cgs <- COM_D_CpG[order(abs(r_decan23), decreasing=TRUE)[1:40]]; w_decan23 <- r_decan23[top40_decan_cgs]
stat23 <- data.frame(); stage_score_df23 <- data.frame()
for(cn in tcga_cans){ cand_m <- list(); cand_c <- list(); cl_n <- paste0(cn,"_c")
for(sx in c("F","M")){ m2_n <- paste0(cn,"_450m_Tm",sx,"2"); c2_n <- paste0(cn,"_c450m",sx,"2")
if(exists(m2_n) && exists(c2_n)){ cand_m[[paste0(sx,"2")]] <- get(m2_n); cand_c[[paste0(sx,"2")]] <- get(c2_n) }
else { nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f")
if(exists(nt_n) && exists(nt_f) && exists(cl_n)){ d_nt <- get(nt_n); f_nt <- get(nt_f); cand_m[[sx]] <- d_nt[, f_nt=="T", drop=FALSE]; cand_c[[sx]] <- get(cl_n) } } }
if(length(cand_m)==0 && exists(paste0(cn,"_450m_Tm")) && exists(cl_n)){ cand_m[["all"]] <- get(paste0(cn,"_450m_Tm")); cand_c[["all"]] <- get(cl_n) }
if(length(cand_m)>0){ m_pool <- NULL; c_pool <- data.frame()
for(k in 1:length(cand_m)){ mat_k <- cand_m[[k]]; cln_k <- cand_c[[k]]; cg_k <- intersect(top40_decan_cgs, rownames(mat_k))
pat_m <- toupper(substr(gsub("[_\\.]", "-", colnames(mat_k)), 1, 12)); pat_c <- toupper(substr(gsub("[_\\.]", "-", cln_k$ID), 1, 12)); com_p <- intersect(pat_m, pat_c)
if(length(com_p)>=5 && length(cg_k)>=20){ sub_mk <- as.matrix(mat_k[cg_k, match(com_p, pat_m), drop=FALSE]); storage.mode(sub_mk) <- "numeric"; sub_ck <- cln_k[match(com_p, pat_c), , drop=FALSE]
m_pool <- if(is.null(m_pool)) sub_mk else { c_com <- intersect(rownames(m_pool), rownames(sub_mk)); cbind(m_pool[c_com, , drop=FALSE], sub_mk[c_com, , drop=FALSE]) }
c_pool <- rbind(c_pool, sub_ck[, c("ID","StageC")]) } }
if(!is.null(m_pool) && nrow(c_pool)>=10){ cg_use <- intersect(top40_decan_cgs, rownames(m_pool)); mat_sub <- m_pool[cg_use, , drop=FALSE]; w_sub <- w_decan23[cg_use]
for(r in 1:nrow(mat_sub)){ na_idx <- which(is.na(mat_sub[r, ])); if(length(na_idx)>0){ rm_v <- mean(mat_sub[r, -na_idx], na.rm=TRUE); mat_sub[r, na_idx] <- ifelse(is.finite(rm_v), rm_v, 0.5) } }
decan_sc <- as.numeric(crossprod(mat_sub, w_sub) / sum(abs(w_sub))); st_ch <- as.character(c_pool$StageC)
vl <- decan_sc[st_ch=="L"]; vm <- decan_sc[st_ch=="M"]; vl <- vl[!is.na(vl)]; vm <- vm[!is.na(vm)]
if(length(vl)>=3 && length(vm)>=3 && var(vl)>0 && var(vm)>0){ tt <- t.test(vm, vl)
stat23 <- rbind(stat23, data.frame(Cancer=cn, Early_N=length(vl), Late_N=length(vm), Mean_Early_Score=mean(vl), Mean_Late_Score=mean(vm), Delta_Late_minus_Early=mean(vm)-mean(vl), P_value=tt$p.value, stringsAsFactors=FALSE))
stage_score_df23 <- rbind(stage_score_df23, data.frame(Cancer=cn, Stage=st_ch[st_ch %in% c("L","M")], Decan_Score=as.numeric(scale(decan_sc[st_ch %in% c("L","M")])), stringsAsFactors=FALSE)) } } } }
if(nrow(stat23)>0){ write.table(stat23, "23_Morpho_Plasticity_Cancer_Stage_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p23a <- ggplot(stat23, aes(x=reorder(Cancer, Delta_Late_minus_Early), y=Delta_Late_minus_Early, fill=P_value<0.05)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="gray70")) + theme_bw() + labs(title="Activation of Skeletal Plasticity Loci in Advanced Cancers (Stage M - L)", x="Solid Cancer", y="Delta Epigenetic Plasticity Score (Advanced - Early Stage)")
  p23b <- ggplot(stage_score_df23, aes(x=Cancer, y=Decan_Score, fill=Stage)) + geom_boxplot(alpha=0.8, outlier.size=0.4) + scale_fill_manual(values=c("L"="#4DAC26","M"="#CA0020"), labels=c("L"="Early (I/II)", "M"="Advanced (III/IV)")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Skeletal Decanalization Epigenetic Score by Clinical Stage", x="Solid Cancer", y="Plasticity Decanalization Score (Z)", fill="Stage")
  ggsave("23_Morpho_Plasticity_Cancer_Stage_Plot.pdf", plot=gridExtra::grid.arrange(p23a, p23b, ncol=2), width=13.5, height=5.5)
} else { cat("경고: 유효한 병기(L/M) 표본이 매칭된 암종이 없어 stat23이 비어 있습니다. 임상 객체(*_c450mF2/M2 등)를 확인하세요.\n") }

################################################################################################
# Module D24. 연령 x 대퇴골 형태 상호작용(Mechano-Senescence) 마커의 암종별 기저 교란 검증
################################################################################################
ag_z24 <- numeric(nrow(inf_p21)); for(g in unique(inf_p21$GEO)){ idx <- which(inf_p21$GEO==g); ag_z24[idx] <- as.numeric(scale(inf_p21$Age[idx])) }
pca_m24 <- prcomp(z_29_21)$x[, 1]; int_term24 <- ag_z24 * pca_m24; X_int24 <- cbind(1, ag_z24, pca_m24, int_term24)
ft_int24 <- lm.fit(x=X_int24, y=t(mat_pr21)); df_r24 <- ft_int24$df.residual; se_int24 <- sqrt(chol2inv(ft_int24$qr$qr)[4,4] * colSums(ft_int24$residuals^2)/df_r24)
t_int24 <- ft_int24$coefficients[4, ] / se_int24; top40_ms_cgs <- COM_D_CpG[order(abs(t_int24), decreasing=TRUE)[1:40]]
stat24 <- data.frame(); ms_plot_df24 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[top40_ms_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[top40_ms_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE); ct <- cor.test(t_int24[top40_ms_cgs], db)
stat24 <- rbind(stat24, data.frame(Cancer=cn, N_Loci=length(top40_ms_cgs), Cor_IntT_Delta_R=ct$estimate, P_val=ct$p.value, Mean_Abs_Delta=mean(abs(db))))
ms_plot_df24 <- rbind(ms_plot_df24, data.frame(Cancer=cn, CpG=top40_ms_cgs, Bone_Int_T=t_int24[top40_ms_cgs], Tumor_Delta_Beta=db)) } } }
write.table(stat24, "24_Mechano_Senescence_Cancer_Shift_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p24a <- ggplot(stat24, aes(x=reorder(Cancer, Cor_IntT_Delta_R), y=Cor_IntT_Delta_R, fill=Cor_IntT_Delta_R>0)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Concordant Remodeling","FALSE"="Discordant")) + theme_bw() + labs(title="Mechano-Senescence (Age x Morphology) Concordance with Cancer Shift", x="Solid Cancer", y="Correlation (Bone Interaction T vs Cancer Delta Beta)", fill="Concordance")
p24b <- ggplot(ms_plot_df24, aes(x=Bone_Int_T, y=Tumor_Delta_Beta, color=Cancer)) + geom_point(alpha=0.6, size=1.5) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + facet_wrap(~Cancer, scales="free_y", ncol=5) + theme_bw() + theme(legend.position="none") + labs(title="Locus-by-Locus: Bone Mechano-Senescence T vs Tumor Delta Beta", x="Bone Age x Morphology Interaction T-stat", y="Tumor Delta Beta")
ggsave("24_Mechano_Senescence_Cancer_Shift_Plot.pdf", plot=gridExtra::grid.arrange(p24a, p24b, ncol=1, heights=c(1, 1.4)), width=13, height=8.5)

################################################################################################
# Module D25. 근위 고관절 추진 지렛대 대 원위 슬관절 하중부 역학 조절 좌위의 골 전이 호발암 특이성 분석 (수정본)
################################################################################################
prox_v25 <- c("ANL","ANH","AND","BNL","PW","LTH","LTN","LGT"); dist_v25 <- c("MCH","MCDp","MCW","LCH","LCDp","LCW","INW","IND","BCW")
z_prox25 <- rowMeans(z_29_21[, prox_v25]); z_dist25 <- rowMeans(z_29_21[, dist_v25]); pd_alloc25 <- z_prox25 - z_dist25
r_pd25 <- suppressWarnings(as.numeric(cor(pd_alloc25, t(mat_pr21)))); names(r_pd25) <- COM_D_CpG; top40_pd_cgs <- COM_D_CpG[order(abs(r_pd25), decreasing=TRUE)[1:40]]
trop_cans25 <- c("BRCA","PRAD","LUAD"); nontrop_cans25 <- c("COAD","KIRC","LIHC","BLCA"); stat25 <- data.frame(); trop_score_df25 <- data.frame()
for(cn in c(trop_cans25, nontrop_cans25)){ cand_t <- list(); cand_n <- list()
for(sx in c("F","M")){ nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f")
if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); cand_t[[paste0(sx,"_T")]] <- d[, f=="T", drop=FALSE]; cand_n[[paste0(sx,"_N")]] <- d[, f=="N", drop=FALSE] }
tm_s <- paste0(cn,"_450m_Tm",sx); nm_s <- paste0(cn,"_450m_Nm",sx)
if(length(cand_t)==0 && exists(tm_s)) cand_t[[paste0(sx,"_Tm")]] <- get(tm_s)
if(length(cand_n)==0 && exists(nm_s)) cand_n[[paste0(sx,"_Nm")]] <- get(nm_s) }
if(length(cand_t)==0 && exists(paste0(cn,"_450m_Tm"))) cand_t[["all_T"]] <- get(paste0(cn,"_450m_Tm"))
if(length(cand_n)==0 && exists(paste0(cn,"_450m_Nm"))) cand_n[["all_N"]] <- get(paste0(cn,"_450m_Nm"))
if(length(cand_t)>0 && length(cand_n)>0){ all_m <- c(cand_t, cand_n); com_cg <- top40_pd_cgs
for(k in 1:length(all_m)){ com_cg <- intersect(com_cg, rownames(all_m[[k]])) }
if(length(com_cg)>=10){ m_tm <- NULL; for(k in 1:length(cand_t)){ sub_m <- as.matrix(cand_t[[k]][com_cg, , drop=FALSE]); m_tm <- if(is.null(m_tm)) sub_m else cbind(m_tm, sub_m) }
m_nm <- NULL; for(k in 1:length(cand_n)){ sub_m <- as.matrix(cand_n[[k]][com_cg, , drop=FALSE]); m_nm <- if(is.null(m_nm)) sub_m else cbind(m_nm, sub_m) }
storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=2 && ncol(m_nm)>=2){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE); grp_t <- ifelse(cn %in% trop_cans25, "Bone-Tropic", "Non-Tropic"); w_sub <- r_pd25[com_cg]
p_sh <- if(sum(w_sub>0)>0) mean(db[w_sub>0]) else 0; d_sh <- if(sum(w_sub<0)>0) mean(db[w_sub<0]) else 0
stat25 <- rbind(stat25, data.frame(Cancer=cn, Tropism=grp_t, N_Loci=length(com_cg), Mean_Abs_Shift=mean(abs(db)), Mean_Signed_Shift=mean(db), Prox_Shift=p_sh, Dist_Shift=d_sh, stringsAsFactors=FALSE))
trop_score_df25 <- rbind(trop_score_df25, data.frame(Cancer=cn, Tropism=grp_t, CpG=com_cg, Delta_Beta=db, Bone_PD_Alloc_R=w_sub, stringsAsFactors=FALSE)) } } } }
if(nrow(stat25)>0){ write.table(stat25, "25_Proximal_Distal_Lever_Cancer_Tropism_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  v_tr <- stat25$Mean_Abs_Shift[stat25$Tropism=="Bone-Tropic"]; v_nt <- stat25$Mean_Abs_Shift[stat25$Tropism=="Non-Tropic"]
  tt25_p <- if(length(v_tr)>=2 && length(v_nt)>=2 && var(v_tr)>0 && var(v_nt)>0) t.test(v_tr, v_nt)$p.value else NA
  p25a <- ggplot(stat25, aes(x=reorder(Cancer, Mean_Abs_Shift), y=Mean_Abs_Shift, fill=Tropism)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title=paste0("Dysregulation of Proximal-Distal Lever Loci (T-test P = ", ifelse(is.na(tt25_p), "NA", signif(tt25_p, 3)), ")"), x="Solid Cancer", y="Mean |Delta Beta| Across Lever Allocation Loci")
  p25b <- ggplot(trop_score_df25, aes(x=Tropism, y=Delta_Beta, fill=Tropism)) + geom_boxplot(alpha=0.8, outlier.size=0.4) + facet_wrap(~ifelse(Bone_PD_Alloc_R>0, "Proximal Lever Enriched (+)", "Distal Condyle Enriched (-)")) + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + geom_hline(yintercept=0, linetype="dashed") + theme_bw() + labs(title="Directional Shift of Lever Allocation Loci by Cancer Bone Tropism", x="Cancer Tropism Class", y="Tumor Delta Beta (Tumor - Matched Normal)")
  ggsave("25_Proximal_Distal_Lever_Cancer_Tropism_Plot.pdf", plot=gridExtra::grid.arrange(p25a, p25b, ncol=2), width=13.5, height=5.5)
} else { cat("경고: 유효한 정상/종양 쌍이 매칭된 암종이 없어 stat25가 비어 있습니다.\n") }

################################################################################################
# Module D26. 파골/조골 세포 커플링 축(RANKL-OPG, TNFSF11 vs TNFRSF11B) 프로모터 메틸화와 고형암 기질 리모델링 비교
################################################################################################
rank_anno26 <- GPL13534[GPL13534$ID %in% COM_D_CpG & grepl("TNFSF11|TNFRSF11B", as.character(GPL13534$UCSC_RefGene_Name)), ]
rank_anno26$Gene <- ifelse(grepl("TNFSF11", rank_anno26$UCSC_RefGene_Name), "RANKL (TNFSF11)", "OPG (TNFRSF11B)"); rank_cgs26 <- unique(rank_anno26$ID)
stat26 <- data.frame(); rank_plot_df26 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[rank_cgs26, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[rank_cgs26, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ for(cg in rank_cgs26){ gn <- rank_anno26$Gene[rank_anno26$ID==cg][1]; v_hb <- mat_hb[cg, ]; v_tm <- m_tm[cg, ]; v_nm <- m_nm[cg, ]
tt_t <- t.test(v_tm, v_nm); tt_b <- t.test(v_tm, v_hb)
stat26 <- rbind(stat26, data.frame(Cancer=cn, CpG=cg, Gene=gn, Island=as.character(rank_anno26$Relation_to_UCSC_CpG_Island[rank_anno26$ID==cg][1]), Mean_Bone=mean(v_hb,na.rm=TRUE), Mean_Normal=mean(v_nm,na.rm=TRUE), Mean_Tumor=mean(v_tm,na.rm=TRUE), Delta_Tumor_vs_Normal=mean(v_tm,na.rm=TRUE)-mean(v_nm,na.rm=TRUE), Delta_Tumor_vs_Bone=mean(v_tm,na.rm=TRUE)-mean(v_hb,na.rm=TRUE), P_Tumor_vs_Normal=tt_t$p.value, P_Tumor_vs_Bone=tt_b$p.value)) } } } }
write.table(stat26, "26_RANKL_OPG_Axis_Cancer_Shift_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
rank_agg26 <- aggregate(cbind(Delta_Tumor_vs_Normal, Delta_Tumor_vs_Bone) ~ Cancer + Gene, data=stat26, FUN=mean)
p26a <- ggplot(rank_agg26, aes(x=Cancer, y=Gene, fill=Delta_Tumor_vs_Normal)) + geom_tile(color="white") + geom_text(aes(label=round(Delta_Tumor_vs_Normal, 2)), size=3.5, fontface="bold") + scale_fill_gradient2(low="#377EB8", mid="white", high="#E41A1C", midpoint=0) + theme_minimal() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Tumor-Specific Methylation Shift of RANKL vs OPG Promoters", subtitle="Red: Hypermethylation/Repression, Blue: Hypomethylation/Activation (Tumor - Matched Normal)", x="Solid Cancer", y="Osteoclastogenesis Regulator", fill="Delta Beta")
p26b <- ggplot(stat26, aes(x=Gene, y=Mean_Tumor - Mean_Bone, fill=Cancer %in% c("BRCA","PRAD","LUAD"))) + geom_boxplot(alpha=0.8, outlier.size=0.5) + geom_hline(yintercept=0, linetype="dashed", color="red") + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Bone-Tropic Cancers","FALSE"="Other Solid Cancers")) + theme_bw() + labs(title="RANKL/OPG Remodeling Divergence: Solid Tumors vs Healthy Bone Baseline", x="Coupling Axis", y="Delta Beta (Tumor - Healthy Bone)", fill="Cancer Type")
ggsave("26_RANKL_OPG_Axis_Cancer_Shift_Plot.pdf", plot=gridExtra::grid.arrange(p26a, p26b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D27. 골용해성(CTSK, ACP5, MMP9) 대 골형성성(BGLAP, SPP1, RUNX2) 이원적 마커의 암종별 탈조절 분석
################################################################################################
lyt_genes <- c("CTSK","ACP5","MMP9","CLCN7"); bls_genes <- c("BGLAP","SPP1","RUNX2","ALPL")
lyt_anno27 <- GPL13534[GPL13534$ID %in% COM_D_CpG & grepl(paste(lyt_genes, collapse="|"), as.character(GPL13534$UCSC_RefGene_Name)), ]; lyt_anno27$Type <- "Osteolytic"
bls_anno27 <- GPL13534[GPL13534$ID %in% COM_D_CpG & grepl(paste(bls_genes, collapse="|"), as.character(GPL13534$UCSC_RefGene_Name)), ]; bls_anno27$Type <- "Osteoblastic"
meta_bone_anno27 <- rbind(lyt_anno27, bls_anno27); meta_cgs27 <- unique(meta_bone_anno27$ID); stat27 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[meta_cgs27, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[meta_cgs27, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE)
for(cg in meta_cgs27){ tp <- meta_bone_anno27$Type[meta_bone_anno27$ID==cg][1]; gn <- sub(";.*","",as.character(meta_bone_anno27$UCSC_RefGene_Name[meta_bone_anno27$ID==cg][1]))
stat27 <- rbind(stat27, data.frame(Cancer=cn, CpG=cg, Gene=gn, Type=tp, Delta_Beta=db[cg], Tumor_Mean=mean(m_tm[cg, ], na.rm=TRUE), Normal_Mean=mean(m_nm[cg, ], na.rm=TRUE), Bone_Mean=mean(mat_hb[cg, ], na.rm=TRUE))) } } } }
write.table(stat27, "27_Osteolytic_vs_Osteoblastic_Shift_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
stat27_agg <- aggregate(Delta_Beta ~ Cancer + Type, data=stat27, FUN=mean); stat27_w <- reshape(stat27_agg, idvar="Cancer", timevar="Type", direction="wide")
stat27_w$Coupling_Index <- stat27_w$Delta_Beta.Osteolytic - stat27_w$Delta_Beta.Osteoblastic
p27a <- ggplot(stat27, aes(x=Cancer, y=Delta_Beta, fill=Type)) + geom_boxplot(alpha=0.8, outlier.size=0.4) + geom_hline(yintercept=0, linetype="dashed", color="red") + scale_fill_manual(values=c("Osteolytic"="#E41A1C","Osteoblastic"="#377EB8")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Tumor Dysregulation of Osteolytic vs Osteoblastic Machinery", x="Solid Cancer", y="Delta Beta (Tumor - Matched Normal)", fill="Functional Axis")
p27b <- ggplot(stat27_w, aes(x=reorder(Cancer, Coupling_Index), y=Coupling_Index, fill=Coupling_Index>0)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Osteolytic Skewed","FALSE"="Osteoblastic Skewed")) + theme_bw() + labs(title="Osteolytic vs Osteoblastic Skewing Index Across 10 Cancers", x="Solid Cancer", y="Skewing Index (Delta_Lytic - Delta_Blastic)", fill="Direction")
ggsave("27_Osteolytic_vs_Osteoblastic_Shift_Plot.pdf", plot=gridExtra::grid.arrange(p27a, p27b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D28. 정상 뼈 기저선 대비 고형암의 후성유전학적 정보 발산도(Manhattan / Epigenetic Divergence Distance)
################################################################################################
ref_bone28 <- rowMeans(mat_hb, na.rm=TRUE); div_stat28 <- data.frame(); div_samples28 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n)){ m_tm <- as.matrix(get(tm_n)[COM_D_CpG, , drop=FALSE]); storage.mode(m_tm) <- "numeric"
d_tm <- colMeans(abs(m_tm - ref_bone28), na.rm=TRUE); div_samples28 <- rbind(div_samples28, data.frame(Sample=colnames(m_tm), Cancer=cn, State="Tumor", Divergence=d_tm)) }
if(exists(nm_n)){ m_nm <- as.matrix(get(nm_n)[COM_D_CpG, , drop=FALSE]); storage.mode(m_nm) <- "numeric"
d_nm <- colMeans(abs(m_nm - ref_bone28), na.rm=TRUE); div_samples28 <- rbind(div_samples28, data.frame(Sample=colnames(m_nm), Cancer=cn, State="Normal", Divergence=d_nm)) } }
d_hb <- colMeans(abs(mat_hb - ref_bone28), na.rm=TRUE); div_samples28 <- rbind(div_samples28, data.frame(Sample=colnames(mat_hb), Cancer="Human_Bone", State="Bone", Divergence=d_hb))
for(cn in tcga_cans){ vt <- div_samples28$Divergence[div_samples28$Cancer==cn & div_samples28$State=="Tumor"]; vn <- div_samples28$Divergence[div_samples28$Cancer==cn & div_samples28$State=="Normal"]
if(length(vt)>=3){ div_stat28 <- rbind(div_stat28, data.frame(Cancer=cn, N_Tumor=length(vt), N_Normal=length(vn), Mean_Tumor_Div=mean(vt), Mean_Normal_Div=ifelse(length(vn)>0,mean(vn),NA), Bone_Baseline_Div=mean(d_hb), Delta_Div_Tumor_vs_Normal=ifelse(length(vn)>0,mean(vt)-mean(vn),NA), P_Tumor_vs_Normal=ifelse(length(vn)>=3, t.test(vt, vn)$p.value, NA))) } }
write.table(div_stat28, "28_Bone_Epigenetic_Divergence_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p28a <- ggplot(div_samples28, aes(x=reorder(Cancer, Divergence, FUN=median), y=Divergence, fill=State)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + geom_hline(yintercept=mean(d_hb), linetype="dashed", color="darkgreen", size=1) + scale_fill_manual(values=c("Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Genome-Wide Epigenetic Divergence Distance (Mean |Beta - Bone Ref|)", x="Cohort / Cancer", y="Mean Absolute Distance from Healthy Bone Baseline", fill="Tissue State")
p28b <- ggplot(div_stat28, aes(x=reorder(Cancer, Mean_Tumor_Div), y=Mean_Tumor_Div, fill=Cancer %in% c("BRCA","PRAD","LUAD"))) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Bone-Tropic","FALSE"="Other Solid Cancers")) + theme_bw() + labs(title="Tumor Divergence from Normal Bone by Tropism Class", x="Solid Cancer", y="Mean Tumor Divergence Distance", fill="Tropism")
ggsave("28_Bone_Epigenetic_Divergence_Plot.pdf", plot=gridExtra::grid.arrange(p28a, p28b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D29. 골 전이 호발암의 종양 후성유전학적 노화 가속도(EAD)와 임상 병기(Stage) 및 고령화 시너지 분석 (완전수정본)
################################################################################################
mat_hb29 <- cbind(as.matrix(GSE64490_m[COM_D_CpG, ]), as.matrix(GSE138307m[COM_D_CpG, ])); storage.mode(mat_hb29) <- "numeric"; age_hb29 <- c(GSE64490_i$Age, GSE138307i$Age)
r_b29 <- suppressWarnings(as.numeric(cor(age_hb29, t(mat_hb29)))); names(r_b29) <- COM_D_CpG; top_clk29 <- COM_D_CpG[order(abs(r_b29), decreasing=TRUE)[1:50]]; w_clk29 <- r_b29[top_clk29]
sc_hb29 <- as.numeric(crossprod(mat_hb29[top_clk29, ], w_clk29)/sum(abs(w_clk29))); ft_b29 <- lm(sc_hb29 ~ age_hb29); sd_b29 <- sd(resid(ft_b29))
m29_df <- data.frame(); stat29 <- data.frame(); trop_cans29 <- c("BRCA","PRAD","LUAD")
for(cn in tcga_cans){ cl_n <- paste0(cn,"_c"); if(exists(cl_n)){ c_df <- get(cl_n); if("StageC" %in% colnames(c_df)){ cand_m <- list()
for(sx in c("F","M")){ nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f"); if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); cand_m[[sx]] <- d[, f=="T", drop=FALSE] } else if(exists(paste0(cn,"_450m_Tm",sx,"2"))){ cand_m[[sx]] <- get(paste0(cn,"_450m_Tm",sx,"2")) } }
if(length(cand_m)>0){ com_cgs29 <- top_clk29; for(k in 1:length(cand_m)){ com_cgs29 <- intersect(com_cgs29, rownames(cand_m[[k]])) }
if(length(com_cgs29)>=20){ m_pool <- NULL; for(k in 1:length(cand_m)){ sub_m <- as.matrix(cand_m[[k]][com_cgs29, , drop=FALSE]); m_pool <- if(is.null(m_pool)) sub_m else cbind(m_pool, sub_m) }
pat_m <- toupper(substr(gsub("[_\\.]", "-", colnames(m_pool)), 1, 12)); pat_c <- toupper(substr(gsub("[_\\.]", "-", c_df$ID), 1, 12)); com_p <- intersect(pat_m, pat_c)
if(length(com_p)>=10){ idx_m <- match(com_p, pat_m); idx_c <- match(com_p, pat_c); sub_mat <- as.matrix(m_pool[, idx_m, drop=FALSE]); storage.mode(sub_mat) <- "numeric"
for(r in 1:nrow(sub_mat)){ na_idx <- which(is.na(sub_mat[r, ])); if(length(na_idx)>0){ rm_v <- mean(sub_mat[r, -na_idx], na.rm=TRUE); sub_mat[r, na_idx] <- ifelse(is.finite(rm_v), rm_v, 0.5) } }
sc_t <- as.numeric(crossprod(sub_mat, w_clk29[com_cgs29])/sum(abs(w_clk29[com_cgs29]))); a_t <- suppressWarnings(as.numeric(as.character(c_df$Age[idx_c]))); sx_t <- toupper(substr(as.character(c_df$Sex[idx_c]), 1, 1)); st_t <- as.character(c_df$StageC[idx_c])
ead_t <- (sc_t - predict(ft_b29, newdata=data.frame(age_hb29=a_t))) / sd_b29; ok_idx <- which(is.finite(ead_t) & is.finite(a_t) & st_t %in% c("L","M") & sx_t %in% c("F","M"))
if(length(ok_idx)>=10){ sub_df <- data.frame(Sample=colnames(sub_mat)[ok_idx], Cancer=cn, Tropism=ifelse(cn %in% trop_cans29, "Bone-Tropic", "Non-Tropic"), Age=a_t[ok_idx], Sex=sx_t[ok_idx], StageC=st_t[ok_idx], EAD=ead_t[ok_idx], stringsAsFactors=FALSE)
m29_df <- rbind(m29_df, sub_df); vl <- sub_df$EAD[sub_df$StageC=="L"]; vm <- sub_df$EAD[sub_df$StageC=="M"]
if(length(vl)>=3 && length(vm)>=3 && is.finite(var(vl)) && var(vl)>0 && is.finite(var(vm)) && var(vm)>0){ fml <- as.formula(if(length(unique(sub_df$Sex))>1) "EAD ~ Age + StageC + Sex" else "EAD ~ Age + StageC"); ft_reg <- suppressWarnings(lm(fml, data=sub_df)); cf <- summary(ft_reg)$coefficients; st_row <- grep("^StageC", rownames(cf))[1]
if(!is.na(st_row) && "Age" %in% rownames(cf)){ stat29 <- rbind(stat29, data.frame(Cancer=cn, Tropism=sub_df$Tropism[1], N=nrow(sub_df), Age_Beta=cf["Age",1], Age_P=cf["Age",4], Stage_M_Beta=cf[st_row,1], Stage_M_P=cf[st_row,4], stringsAsFactors=FALSE)) } } } } } } } } }
if(nrow(stat29)>0){ write.table(stat29, "29_Bone_EAD_Stage_Age_Synergy_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p29a <- ggplot(m29_df, aes(x=Age, y=EAD, color=StageC)) + geom_point(alpha=0.5, size=1.5) + geom_smooth(method="lm", formula=y~x, se=FALSE, size=1) + facet_wrap(~Tropism, ncol=2) + scale_color_manual(values=c("L"="#4DAC26","M"="#CA0020"), labels=c("L"="Early (I/II)", "M"="Advanced (III/IV)")) + theme_bw() + labs(title="Skeletal-Referenced Tumor EAD vs Chronological Age by Clinical Stage", x="Chronological Age (Years)", y="Tumor EAD (Z-score vs Bone)", color="Stage")
  p29b <- ggplot(stat29, aes(x=reorder(Cancer, Stage_M_Beta), y=Stage_M_Beta, fill=Tropism)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title="Independent Effect of Advanced Stage on Bone-EAD (Multivariate Beta)", x="Solid Cancer", y="Regression Beta (Stage Advanced vs Early)")
  ggsave("29_Bone_EAD_Stage_Age_Synergy_Plot.pdf", plot=gridExtra::grid.arrange(p29a, p29b, ncol=2), width=13.5, height=5.5)
} else { cat("경고: 유효 병기 매칭 검체가 부족하여 stat29가 비어 있습니다.\n") }

################################################################################################
# Module D30. 비선형 노화 임계 전이(Quadratic Age^2 Curvature) 좌위의 고형암 종양화 가속 메커니즘 검증
################################################################################################
age_hb30 <- c(GSE64490_i$Age, GSE138307i$Age); z_ag_hb30 <- as.numeric(scale(age_hb30))
X_hb30 <- cbind(1, z_ag_hb30, z_ag_hb30^2); ft_hb30 <- lm.fit(x=X_hb30, y=t(mat_hb)); df_hb30 <- ft_hb30$df.residual
se_hb30 <- sqrt(chol2inv(ft_hb30$qr$qr)[3,3] * colSums(ft_hb30$residuals^2)/df_hb30); t_quad_hb30 <- ft_hb30$coefficients[3, ] / se_hb30
top50_quad_cgs <- COM_D_CpG[order(abs(t_quad_hb30), decreasing=TRUE)[1:50]]; stat30 <- data.frame(); quad_plot_df30 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[top50_quad_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[top50_quad_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE); ct <- cor.test(t_quad_hb30[top50_quad_cgs], db)
stat30 <- rbind(stat30, data.frame(Cancer=cn, N_Loci=length(top50_quad_cgs), Cor_Quad_CancerShift_R=ct$estimate, P_val=ct$p.value, Mean_Delta=mean(db), Conc_Pct=mean(t_quad_hb30[top50_quad_cgs]*db > 0)*100))
quad_plot_df30 <- rbind(quad_plot_df30, data.frame(Cancer=cn, CpG=top50_quad_cgs, Bone_Quad_T=t_quad_hb30[top50_quad_cgs], Cancer_Delta_Beta=db)) } } }
write.table(stat30, "30_Nonlinear_Aging_Cancer_TippingPoint_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p30a <- ggplot(stat30, aes(x=reorder(Cancer, Cor_Quad_CancerShift_R), y=Cor_Quad_CancerShift_R, fill=Cor_Quad_CancerShift_R>0)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Late-Life Accelerated Shift","FALSE"="Decelerated Shift")) + theme_bw() + labs(title="Non-Linear Epigenetic Curvature: \nBone Age^2 T vs Solid Cancer Delta Beta", x="Solid Cancer", y="Correlation (Pearson R)", fill="Tipping-Point")
p30b <- ggplot(quad_plot_df30, aes(x=Bone_Quad_T, y=Cancer_Delta_Beta, color=Cancer)) + geom_point(alpha=0.6, size=1.5) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + facet_wrap(~Cancer, scales="free_y", ncol=5) + theme_bw() + theme(legend.position="none") + labs(title="Locus-by-Locus: Bone Curvature (Age^2 T) vs Tumor Delta Beta Across 10 Cancers", x="Bone Quadratic Aging Curvature T-statistic", y="Tumor Delta Beta (Tumor - Matched Normal)")
#ggsave("30_Nonlinear_Aging_Cancer_TippingPoint_Plot.pdf", plot=gridExtra::grid.arrange(p30a, p30b, ncol=1, heights=c(1, 1.4)), width=13, height=8.5)
ggsave("30_Nonlinear_Aging_Cancer_TippingPoint_Plot.pdf", plot=p30a, width=5 , height=5)
ggsave("30_Nonlinear_Aging_Cancer_TippingPoint_Plot.pdf", plot=p30b, width=13, height=5)

################################################################################################
# Module D31. 염색체 공간 자기상관(Lag-1 Autocorrelation) 및 국소 메틸화 거칠기(Roughness) 비교 (정상 뼈 vs 고형암)
################################################################################################
if(!exists("GPL13534")){ GPL13534 <- as.data.frame(fread("G:/내 드라이브/Rrunning/Annotations/21_Methyl/GPL13534-11288.txt.gz", header=TRUE, sep="\t"))[, c("ID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")] }
g_ord31 <- GPL13534[GPL13534$ID %in% COM_D_CpG & GPL13534$CHR %in% as.character(1:22), ]; g_ord31$CHR_N <- as.numeric(g_ord31$CHR); g_ord31$MAP_N <- as.numeric(g_ord31$MAPINFO)
g_ord31 <- g_ord31[order(g_ord31$CHR_N, g_ord31$MAP_N), ]; ord_cgs31 <- g_ord31$ID
same_pair31 <- which(g_ord31$CHR_N[-nrow(g_ord31)] == g_ord31$CHR_N[-1] & (g_ord31$MAP_N[-1] - g_ord31$MAP_N[-nrow(g_ord31)]) < 50000)
m_hb_ord <- mat_hb[ord_cgs31, ]; diff_hb <- m_hb_ord[-1, ] - m_hb_ord[-nrow(m_hb_ord), ]; rough_hb <- sqrt(colMeans(diff_hb[same_pair31, ]^2, na.rm=TRUE))
acf_hb <- numeric(ncol(m_hb_ord)); for(k in 1:ncol(m_hb_ord)){ acf_hb[k] <- cor(m_hb_ord[same_pair31, k], m_hb_ord[same_pair31+1, k], use="complete.obs") }
stat31 <- data.frame(Sample=colnames(mat_hb), Cancer="Human_Bone", State="Healthy_Bone", Lag1_Autocor=acf_hb, Roughness=rough_hb)
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n)){ m_tm <- as.matrix(get(tm_n)[ord_cgs31, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; s_t <- head(colnames(m_tm), 25); d_t <- m_tm[-1, s_t] - m_tm[-nrow(m_tm), s_t]; r_t <- sqrt(colMeans(d_t[same_pair31, ]^2, na.rm=TRUE)); a_t <- numeric(length(s_t)); for(k in 1:length(s_t)){ a_t[k] <- cor(m_tm[same_pair31, s_t[k]], m_tm[same_pair31+1, s_t[k]], use="complete.obs") }; stat31 <- rbind(stat31, data.frame(Sample=s_t, Cancer=cn, State="Tumor", Lag1_Autocor=a_t, Roughness=r_t)) }
if(exists(nm_n)){ m_nm <- as.matrix(get(nm_n)[ord_cgs31, , drop=FALSE]); storage.mode(m_nm) <- "numeric"; s_n <- head(colnames(m_nm), 10); if(length(s_n)>0){ d_n <- m_nm[-1, s_n] - m_nm[-nrow(m_nm), s_n]; r_n <- sqrt(colMeans(d_n[same_pair31, ]^2, na.rm=TRUE)); a_n <- numeric(length(s_n)); for(k in 1:length(s_n)){ a_n[k] <- cor(m_nm[same_pair31, s_n[k]], m_nm[same_pair31+1, s_n[k]], use="complete.obs") }; stat31 <- rbind(stat31, data.frame(Sample=s_n, Cancer=cn, State="Normal", Lag1_Autocor=a_n, Roughness=r_n)) } } }
stat31 <- stat31[is.finite(stat31$Lag1_Autocor) & is.finite(stat31$Roughness), ]; write.table(stat31, "31_Genomic_Spatial_Autocorrelation_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
sum31 <- aggregate(cbind(Lag1_Autocor, Roughness) ~ Cancer + State, data=stat31, FUN=mean); write.table(sum31, "31_Spatial_Autocor_Summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
p31a <- ggplot(stat31, aes(x=reorder(Cancer, Lag1_Autocor, FUN=median), y=Lag1_Autocor, fill=State)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + scale_fill_manual(values=c("Healthy_Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Genomic Spatial Coherence (Lag-1 Autocorrelation < 50kb)", x="Cancer Type", y="Lag-1 Spatial Autocorrelation (r)")
p31b <- ggplot(stat31, aes(x=reorder(Cancer, Roughness, FUN=median), y=Roughness, fill=State)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + scale_fill_manual(values=c("Healthy_Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Epigenomic Roughness (Spatial RMS Delta Beta)", x="Cancer Type", y="Local Methylation Roughness")
ggsave("31_Genomic_Spatial_Autocorrelation_Plot.pdf", plot=gridExtra::grid.arrange(p31a, p31b, ncol=2), width=14, height=5.5)

################################################################################################
# Module D32. 골격계 장수형 형태(Longevity Morphotype, -EAA) 연관 좌위의 고형암 억제(Tumor Suppressive) 활성 검증
################################################################################################
w_age32 <- (GSE138307c$M_HyPi[match(COM_D_CpG, GSE138307c$CpG)] + GSE64490_c$F_HyPi[match(COM_D_CpG, GSE64490_c$CpG)])/2
raw_clk32 <- as.numeric(crossprod(mat_hb, w_age32)/sum(abs(w_age32))); age_hb32 <- c(GSE64490_i$Age, GSE138307i$Age)
eaa_hb32 <- as.numeric(scale(resid(lm(raw_clk32 ~ age_hb32))))
r_eaa32 <- suppressWarnings(as.numeric(cor(eaa_hb32, t(mat_hb)))); names(r_eaa32) <- COM_D_CpG
top50_long_cgs <- COM_D_CpG[order(r_eaa32)[1:50]]; top50_acc_cgs <- COM_D_CpG[order(r_eaa32, decreasing=TRUE)[1:50]]
stat32 <- data.frame(); for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[COM_D_CpG, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[COM_D_CpG, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE)
stat32 <- rbind(stat32, data.frame(Cancer=cn, Type="Longevity Loci (-EAA)", Mean_Cancer_Shift=mean(db[top50_long_cgs]), Mean_Abs_Shift=mean(abs(db[top50_long_cgs])), Hyper_N=sum(db[top50_long_cgs]>0.05), Hypo_N=sum(db[top50_long_cgs]< -0.05)))
stat32 <- rbind(stat32, data.frame(Cancer=cn, Type="Accelerated Loci (+EAA)", Mean_Cancer_Shift=mean(db[top50_acc_cgs]), Mean_Abs_Shift=mean(abs(db[top50_acc_cgs])), Hyper_N=sum(db[top50_acc_cgs]>0.05), Hypo_N=sum(db[top50_acc_cgs]< -0.05))) } } }
write.table(stat32, "32_Longevity_Morphotype_Cancer_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p32a <- ggplot(stat32, aes(x=Cancer, y=Mean_Cancer_Shift, fill=Type)) + geom_bar(stat="identity", position="dodge", color="black", width=0.7) + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("Longevity Loci (-EAA)"="#1B9E77","Accelerated Loci (+EAA)"="#D95F02")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Tumor-Specific Methylation Shift of Skeletal Longevity vs Accelerated Aging Loci", x="Solid Cancer", y="Mean Delta Beta (Tumor - Matched Normal)", fill="Skeletal Aging Class")
p32b <- ggplot(stat32, aes(x=Type, y=Mean_Cancer_Shift, fill=Type)) + geom_boxplot(alpha=0.8, outlier.size=0.5) + geom_hline(yintercept=0, linetype="dashed", color="red") + scale_fill_manual(values=c("Longevity Loci (-EAA)"="#1B9E77","Accelerated Loci (+EAA)"="#D95F02")) + theme_bw() + labs(title="Overall Tumor Silencing (Hypermethylation) of Longevity Markers", x="Loci Category", y="Mean Tumor Shift Across 10 Cancers")
ggsave("32_Longevity_Morphotype_Cancer_Plot.pdf", plot=gridExtra::grid.arrange(p32a, p32b, ncol=2), width=13.5, height=5.5)

################################################################################################
# Module D33. 중간 메틸화(Intermediate Methylation, 0.3 <= beta <= 0.7) 쌍안정성 붕괴와 암종별 전이성 프로파일 매핑
################################################################################################
imf_hb33 <- rowMeans(mat_hb >= 0.3 & mat_hb <= 0.7, na.rm=TRUE)
top500_imf_cgs <- COM_D_CpG[order(imf_hb33, decreasing=TRUE)[1:500]]
stat33 <- data.frame(); imf_dist_df33 <- data.frame(Sample=colnames(mat_hb), Cancer="Human_Bone", State="Healthy_Bone", IMF_Fraction=colMeans(mat_hb[top500_imf_cgs, ] >= 0.3 & mat_hb[top500_imf_cgs, ] <= 0.7, na.rm=TRUE))
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n)){ m_tm <- as.matrix(get(tm_n)[top500_imf_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; s_t <- head(colnames(m_tm), 30)
imf_t <- colMeans(m_tm[, s_t] >= 0.3 & m_tm[, s_t] <= 0.7, na.rm=TRUE); imf_dist_df33 <- rbind(imf_dist_df33, data.frame(Sample=s_t, Cancer=cn, State="Tumor", IMF_Fraction=imf_t)) }
if(exists(nm_n)){ m_nm <- as.matrix(get(nm_n)[top500_imf_cgs, , drop=FALSE]); storage.mode(m_nm) <- "numeric"; s_n <- head(colnames(m_nm), 10)
if(length(s_n)>0){ imf_n <- colMeans(m_nm[, s_n] >= 0.3 & m_nm[, s_n] <= 0.7, na.rm=TRUE); imf_dist_df33 <- rbind(imf_dist_df33, data.frame(Sample=s_n, Cancer=cn, State="Normal", IMF_Fraction=imf_n)) } } }
for(cn in tcga_cans){ vt <- imf_dist_df33$IMF_Fraction[imf_dist_df33$Cancer==cn & imf_dist_df33$State=="Tumor"]; vn <- imf_dist_df33$IMF_Fraction[imf_dist_df33$Cancer==cn & imf_dist_df33$State=="Normal"]
stat33 <- rbind(stat33, data.frame(Cancer=cn, Tropism=ifelse(cn %in% c("BRCA","PRAD","LUAD"), "Bone-Tropic", "Non-Tropic"), Mean_Tumor_IMF=mean(vt), Mean_Normal_IMF=ifelse(length(vn)>0, mean(vn), NA), Bone_Ref_IMF=mean(imf_dist_df33$IMF_Fraction[imf_dist_df33$State=="Healthy_Bone"]), Loss_of_IMF=mean(imf_dist_df33$IMF_Fraction[imf_dist_df33$State=="Healthy_Bone"]) - mean(vt))) }
write.table(stat33, "33_IMF_Bistability_Breakdown_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p33a <- ggplot(imf_dist_df33, aes(x=reorder(Cancer, -IMF_Fraction, FUN=median), y=IMF_Fraction, fill=State)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + scale_fill_manual(values=c("Healthy_Bone"="#1B9E77","Normal"="#377EB8","Tumor"="#E41A1C")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Breakdown of Chromatin Metastability (IMF Fraction on Top 500 Bone IMF CpGs)", x="Cohort / Cancer", y="Proportion of CpGs in Dynamic Range (0.3 - 0.7)")
p33b <- ggplot(stat33, aes(x=reorder(Cancer, Loss_of_IMF), y=Loss_of_IMF, fill=Tropism)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title="Loss of Metastable Intermediate Methylation in Tumors vs Bone", x="Solid Cancer", y="IMF Depletion Index (Bone IMF - Tumor IMF)")
ggsave("33_IMF_Bistability_Breakdown_Plot.pdf", plot=gridExtra::grid.arrange(p33a, p33b, ncol=2), width=14, height=5.5)

################################################################################################
# Module D34. 영장류 대뇌화/생태 축(EQ & Day Range) 연결 골격 좌위와 고형암 발암 스트레스 수렴 분석
################################################################################################
eco_df34 <- data.frame(GEO=c("GSE103328","GSE103280","GSE103271","GSE103279"), EQ=c(1.74, 2.05, 2.18, 2.36), Day_Range_km=c(0.9, 1.8, 2.4, 5.5))
sp_cents34 <- matrix(0, nrow=length(COM_D_CpG), ncol=4, dimnames=list(COM_D_CpG, eco_df34$GEO))
for(g in eco_df34$GEO){ sp_cents34[, g] <- rowMeans(raw_p12[, Primate2_Inf$GEO==g, drop=FALSE], na.rm=TRUE) }
r_eq34 <- suppressWarnings(as.numeric(cor(eco_df34$EQ, t(sp_cents34)))); names(r_eq34) <- COM_D_CpG; top40_eq_cgs <- COM_D_CpG[order(abs(r_eq34), decreasing=TRUE)[1:40]]
stat34 <- data.frame(); eq_plot_df34 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[top40_eq_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[top40_eq_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE); ct <- cor.test(r_eq34[top40_eq_cgs], db)
stat34 <- rbind(stat34, data.frame(Cancer=cn, N_Loci=length(top40_eq_cgs), Cor_EQ_Delta_R=ct$estimate, P_val=ct$p.value, Mean_Abs_Delta=mean(abs(db))))
eq_plot_df34 <- rbind(eq_plot_df34, data.frame(Cancer=cn, CpG=top40_eq_cgs, Primate_EQ_R=r_eq34[top40_eq_cgs], Tumor_Delta_Beta=db)) } } }
write.table(stat34, "34_Encephalization_EQ_Bone_Cancer_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p34a <- ggplot(stat34, aes(x=reorder(Cancer, Cor_EQ_Delta_R), y=Cor_EQ_Delta_R, fill=Cor_EQ_Delta_R>0)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#377EB8"), labels=c("TRUE"="Positive Alignment","FALSE"="Negative Divergence")) + theme_bw() + labs(title="Primate Encephalization (EQ) Bone Axis vs Cancer Dysregulation", x="Solid Cancer", y="Correlation (EQ Pearson R vs Cancer Delta Beta)", fill="Direction")
p34b <- ggplot(eq_plot_df34, aes(x=Primate_EQ_R, y=Tumor_Delta_Beta, color=Cancer)) + geom_point(alpha=0.6, size=1.5) + geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") + facet_wrap(~Cancer, scales="free_y", ncol=5) + theme_bw() + theme(legend.position="none") + labs(title="Encephalization-Coupled Bone Loci: Primate EQ Axis vs Cancer Shift", x="Primate Encephalization Quotient (EQ) Correlation R", y="Tumor Delta Beta")
ggsave("34_Encephalization_EQ_Bone_Cancer_Plot.pdf", plot=gridExtra::grid.arrange(p34a, p34b, ncol=1, heights=c(1, 1.4)), width=13, height=8.5)

################################################################################################
# Module D35. 대퇴골 근위-원위 6단계 종축 비례 사슬(Proximo-Distal Scaling Chain) 좌위의 고형암 원격 전이 축 결합 검증
################################################################################################
z1_hd35 <- rowMeans(log(as.matrix(inf_p21[,c("HH","HL","HW")]))); z2_nk35 <- rowMeans(log(as.matrix(inf_p21[,c("ANL","ANH","AND","BNL","PW","LTH","LTN","LGT")])))
z3_ss35 <- rowMeans(log(as.matrix(inf_p21[,c("SSWd","SSDp")]))); z4_ms35 <- rowMeans(log(as.matrix(inf_p21[,c("MSWd","MSDp")])))
z5_is35 <- rowMeans(log(as.matrix(inf_p21[,c("ISWd","ISDp")]))); z6_cd35 <- rowMeans(log(as.matrix(inf_p21[,c("MCH","MCDp","MCW","LCH","LCDp","LCW","INW","IND","BCW")])))
tr_chain35 <- cbind(T1=z2_nk35-z1_hd35, T2=z3_ss35-z2_nk35, T3=z4_ms35-z3_ss35, T4=z5_is35-z4_ms35, T5=z6_cd35-z5_is35)
z_chain35 <- tr_chain35; for(g in unique(inf_p21$GEO)){ idx <- which(inf_p21$GEO==g); for(j in 1:5){ z_chain35[idx, j] <- as.numeric(scale(tr_chain35[idx, j])) } }
r_ch35 <- cor(t(mat_pr21), z_chain35, use="pairwise.complete.obs"); r_ch35[!is.finite(r_ch35)] <- 0
grad_slope35 <- as.numeric((r_ch35 %*% c(-2,-1,0,1,2)) / 10); names(grad_slope35) <- COM_D_CpG
top40_grad_cgs <- COM_D_CpG[order(abs(grad_slope35), decreasing=TRUE)[1:40]]
stat35 <- data.frame(); pd_plot_df35 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[top40_grad_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[top40_grad_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE)
ct <- cor.test(grad_slope35[top40_grad_cgs], db); grp_t <- ifelse(cn %in% c("BRCA","PRAD","LUAD"), "Bone-Tropic", "Non-Tropic")
stat35 <- rbind(stat35, data.frame(Cancer=cn, Tropism=grp_t, N_Loci=length(top40_grad_cgs), Cor_Gradient_Cancer_R=ct$estimate, P_val=ct$p.value, Mean_Abs_Shift=mean(abs(db))))
pd_plot_df35 <- rbind(pd_plot_df35, data.frame(Cancer=cn, Tropism=grp_t, CpG=top40_grad_cgs, Bone_PD_Slope=grad_slope35[top40_grad_cgs], Tumor_Delta_Beta=db)) } } }
write.table(stat35, "35_ProximoDistal_Gradient_Cancer_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p35a <- ggplot(stat35, aes(x=reorder(Cancer, Cor_Gradient_Cancer_R), y=Cor_Gradient_Cancer_R, fill=Tropism)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title="Proximo-Distal Femur Spatial Gradient vs Tumor Dysregulation", x="Solid Cancer", y="Correlation (P-D Gradient Slope vs Cancer Delta Beta)", fill="Tropism")
p35b <- ggplot(pd_plot_df35, aes(x=Bone_PD_Slope, y=Tumor_Delta_Beta, color=Tropism)) + geom_point(alpha=0.6, size=1.5) + geom_smooth(method="lm", se=FALSE) + facet_wrap(~Tropism) + scale_color_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title="Skeletal Spatial Gradient Sensitivity: Bone-Tropic vs Non-Tropic Cancers", x="Femur Proximo-Distal Scaling Gradient Slope", y="Tumor Delta Beta")
ggsave("35_ProximoDistal_Gradient_Cancer_Plot.pdf", plot=gridExtra::grid.arrange(p35a, p35b, ncol=2), width=13.5, height=5.5)


################################################################################################
# Module D36. 29차원 가우시안 그래피컬 모델(GGM) 편상관 네트워크 허브의 암종별 선택적 침식 분석
################################################################################################
z_29_36 <- matrix(0, nrow(Primate2_Inf), 29, dimnames=list(rownames(Primate2_Inf), colnames(Primate2_Inf)[3:31]))
for(g in unique(Primate2_Inf$GEO)){ idx <- which(Primate2_Inf$GEO==g); for(j in 1:29){ v <- Primate2_Inf[idx, j+2]; sv <- as.numeric(scale(v)); sv[!is.finite(sv)] <- 0; z_29_36[idx, j] <- sv } }
S_cov36 <- cor(z_29_36, use="pairwise.complete.obs"); S_cov36[!is.finite(S_cov36)] <- 0; Omega36 <- solve(S_cov36 + diag(0.25, 29)); pcor36 <- -Omega36 / sqrt(outer(diag(Omega36), diag(Omega36))); diag(pcor36) <- 0
edge_cgs36 <- c(); for(i in 1:28){ for(j in (i+1):29){ if(abs(pcor36[i,j]) > 0.15){ ep <- z_29_36[, i] * z_29_36[, j]; rv <- suppressWarnings(as.numeric(cor(ep, t(mat_pr21)))); edge_cgs36 <- c(edge_cgs36, COM_D_CpG[which.max(abs(rv))]) } } }
top30_ggm_cgs <- unique(edge_cgs36)[1:min(30, length(unique(edge_cgs36)))]; stat36 <- data.frame(); ggm_plot_df36 <- data.frame()
for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm"); nm_n <- paste0(cn,"_450m_Nm")
if(exists(tm_n) && exists(nm_n)){ m_tm <- as.matrix(get(tm_n)[top30_ggm_cgs, , drop=FALSE]); m_nm <- as.matrix(get(nm_n)[top30_ggm_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=3 && ncol(m_nm)>=3){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE); grp_t <- ifelse(cn %in% c("BRCA","PRAD","LUAD"), "Bone-Tropic", "Non-Tropic")
stat36 <- rbind(stat36, data.frame(Cancer=cn, Tropism=grp_t, N_Loci=length(top30_ggm_cgs), Mean_Abs_Shift=mean(abs(db)), Mean_Signed_Shift=mean(db), Hyper_Shift_N=sum(db>0.05), Hypo_Shift_N=sum(db< -0.05)))
ggm_plot_df36 <- rbind(ggm_plot_df36, data.frame(Cancer=cn, Tropism=grp_t, CpG=top30_ggm_cgs, Delta_Beta=db)) } } }
write.table(stat36, "36_GGM_Bone_Network_Drivers_Cancer_Shift_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p36a <- ggplot(stat36, aes(x=reorder(Cancer, Mean_Abs_Shift), y=Mean_Abs_Shift, fill=Tropism)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title="Tumor Disturbance of GGM Structural Network Drivers", x="Solid Cancer", y="Mean |Delta Beta| across 30 GGM Edge Drivers")
p36b <- ggplot(ggm_plot_df36, aes(x=reorder(CpG, abs(Delta_Beta), FUN=median), y=Delta_Beta, fill=Tropism)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + geom_hline(yintercept=0, linetype="dashed", color="red") + coord_flip() + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title="Locus-by-Locus Dysregulation of Bone GGM Hubs Across Cancers", x="GGM Driver CpG", y="Tumor Delta Beta")
ggsave("36_GGM_Bone_Network_Drivers_Cancer_Shift_Plot.pdf", plot=gridExtra::grid.arrange(p36a, p36b, ncol=2), width=14, height=6)

################################################################################################
# Module D37. 골격계 의사시간(Pseudotime tau) 다양체 기반 TCGA 암종 악성화 축 매핑 (수정본)
################################################################################################
hp_m37 <- GSE138307c$M_HyPi[match(COM_D_CpG, GSE138307c$CpG)]; hp_f37 <- GSE64490_c$F_HyPi[match(COM_D_CpG, GSE64490_c$CpG)]
w_pt37 <- (hp_m37 + hp_f37) / 2; names(w_pt37) <- COM_D_CpG; w_pt37 <- w_pt37[is.finite(w_pt37)]
top50_pt_cgs <- names(sort(abs(w_pt37), decreasing=TRUE))[1:min(50, length(w_pt37))]; w50_pt <- w_pt37[top50_pt_cgs]
pt_raw_hb <- as.numeric(crossprod(mat_hb[top50_pt_cgs, ], w50_pt) / sum(abs(w50_pt))); tau_hb <- (pt_raw_hb - min(pt_raw_hb)) / max(1e-6, max(pt_raw_hb) - min(pt_raw_hb))
df37_pt <- data.frame(Sample=colnames(mat_hb), Cancer="Human_Bone", State="Healthy_Bone", StageC="Normal", Pseudotime_Tau=tau_hb, stringsAsFactors=FALSE); stat37 <- data.frame()
for(cn in tcga_cans){ cl_n <- paste0(cn, "_c")
if(exists(cl_n)){ c_df <- get(cl_n); cand_m <- list()
for(sx in c("F", "M")){ nt_n <- paste0(cn, "_450m_NT", sx); nt_f <- paste0(cn, "_450m_NT", sx, "f")
if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); cand_m[[sx]] <- d[, f == "T", drop=FALSE] } else if(exists(paste0(cn, "_450m_Tm", sx, "2"))){ cand_m[[sx]] <- get(paste0(cn, "_450m_Tm", sx, "2")) } else if(exists(paste0(cn, "_450m_Tm", sx))){ cand_m[[sx]] <- get(paste0(cn, "_450m_Tm", sx)) } }
if(length(cand_m) == 0 && exists(paste0(cn, "_450m_Tm"))){ cand_m[["all"]] <- get(paste0(cn, "_450m_Tm")) }
if(length(cand_m) > 0){ com_cg37 <- top50_pt_cgs; for(k in 1:length(cand_m)){ com_cg37 <- intersect(com_cg37, rownames(cand_m[[k]])) }
if(length(com_cg37) >= 20){ m_tm <- NULL; for(k in 1:length(cand_m)){ sub_m <- as.matrix(cand_m[[k]][com_cg37, , drop=FALSE]); m_tm <- if(is.null(m_tm)) sub_m else cbind(m_tm, sub_m) }
pat_m <- toupper(substr(gsub("[_\\.]", "-", colnames(m_tm)), 1, 12)); pat_c <- toupper(substr(gsub("[_\\.]", "-", c_df$ID), 1, 12)); com_p <- intersect(pat_m, pat_c)
if(length(com_p) >= 5){ idx_m <- match(com_p, pat_m); idx_c <- match(com_p, pat_c); sub_mat <- as.matrix(m_tm[, idx_m, drop=FALSE]); storage.mode(sub_mat) <- "numeric"
for(r in 1:nrow(sub_mat)){ na_idx <- which(is.na(sub_mat[r, ])); if(length(na_idx) > 0){ rm_v <- mean(sub_mat[r, -na_idx], na.rm=TRUE); sub_mat[r, na_idx] <- ifelse(is.finite(rm_v), rm_v, 0.5) } }
w_sub <- w50_pt[com_cg37]; sc_t <- as.numeric(crossprod(sub_mat, w_sub) / sum(abs(w_sub))); t_tau <- (sc_t - min(pt_raw_hb)) / max(1e-6, max(pt_raw_hb) - min(pt_raw_hb))
sub_c <- c_df[idx_c, ]; st_vec <- if("StageC" %in% colnames(sub_c)) as.character(sub_c$StageC) else rep("X", length(idx_c))
df37_pt <- rbind(df37_pt, data.frame(Sample=colnames(sub_mat), Cancer=cn, State="Tumor", StageC=st_vec, Pseudotime_Tau=t_tau, stringsAsFactors=FALSE))
vl <- t_tau[st_vec == "L"]; vm <- t_tau[st_vec == "M"]; vl <- vl[is.finite(vl)]; vm <- vm[is.finite(vm)]
tt_st <- if(length(vl) >= 3 && length(vm) >= 3 && is.finite(var(vl)) && var(vl) > 0 && is.finite(var(vm)) && var(vm) > 0) suppressWarnings(t.test(vm, vl)$p.value) else NA
m_vl <- if(length(vl) > 0) mean(vl) else NA; m_vm <- if(length(vm) > 0) mean(vm) else NA
stat37 <- rbind(stat37, data.frame(Cancer=cn, N_Tumor=length(t_tau), Mean_Tumor_Tau=mean(t_tau, na.rm=TRUE), Mean_Early_Tau=m_vl, Mean_Late_Tau=m_vm, Delta_Late_vs_Early=ifelse(is.finite(m_vl) && is.finite(m_vm), m_vm - m_vl, NA), Stage_P=tt_st, stringsAsFactors=FALSE)) } } } } }
if(nrow(stat37) > 0){ write.table(stat37, "37_Bone_Pseudotime_TCGA_Projection_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p37a <- ggplot(df37_pt[df37_pt$State=="Tumor", ], aes(x=reorder(Cancer, Pseudotime_Tau, FUN=median), y=Pseudotime_Tau, fill=Cancer)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + geom_hline(yintercept=mean(tau_hb), linetype="dashed", color="darkgreen", size=1) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1), legend.position="none") + labs(title="Projection of Bone Epigenetic Pseudotime (tau) Across TCGA Tumors", subtitle="Dashed line indicates Healthy Bone Pseudotime Reference Centroid", x="Solid Cancer", y="Bone Pseudotime (tau)")
  p37b <- ggplot(df37_pt[df37_pt$StageC %in% c("L","M"), ], aes(x=StageC, y=Pseudotime_Tau, fill=StageC)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + facet_wrap(~Cancer, ncol=5) + scale_fill_manual(values=c("L"="#4DAC26","M"="#CA0020"), labels=c("L"="Early (I/II)","M"="Late (III/IV)")) + theme_bw() + theme(legend.position="bottom") + labs(title="Pseudotime Progression Across Tumor Stages", x="Clinical Stage", y="Pseudotime (tau)", fill="Stage")
  ggsave("37_Bone_Pseudotime_TCGA_Projection_Plot.pdf", plot=gridExtra::grid.arrange(p37a, p37b, ncol=1, heights=c(1, 1.3)), width=13.5, height=8.5)
} else { cat("경고: 유효한 암종 표본이 매칭되지 않아 stat37이 비어 있습니다.\n") }

################################################################################################
# Module D38. 골-내외측 비대칭 과간 하중비(Medial vs Lateral Condyle Dominance) 조절 좌위의 고형암 교란 분석 (수정본)
################################################################################################
c_dom38 <- (log(inf_p21$MCH/inf_p21$LCH) + log(inf_p21$MCDp/inf_p21$LCDp) + log(inf_p21$MCW/inf_p21$LCW))/3
z_cdom38 <- numeric(nrow(inf_p21)); for(g in unique(inf_p21$GEO)){ idx <- which(inf_p21$GEO==g); v <- c_dom38[idx]; if(sum(is.finite(v))>2) z_cdom38[idx] <- as.numeric(scale(v)) }
r_cd38 <- suppressWarnings(as.numeric(cor(z_cdom38, t(mat_pr21)))); names(r_cd38) <- COM_D_CpG; r_cd38[!is.finite(r_cd38)] <- 0
top40_cdom_cgs <- COM_D_CpG[order(abs(r_cd38), decreasing=TRUE)[1:40]]
stat38 <- data.frame(); cdom_plot_df38 <- data.frame()
for(cn in tcga_cans){ cand_t <- list(); cand_n <- list()
for(sx in c("F","M")){ nt_n <- paste0(cn,"_450m_NT",sx); nt_f <- paste0(cn,"_450m_NT",sx,"f")
if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); cand_t[[paste0(sx,"_T")]] <- d[, f=="T", drop=FALSE]; cand_n[[paste0(sx,"_N")]] <- d[, f=="N", drop=FALSE] }
else if(exists(paste0(cn,"_450m_Tm",sx))){ cand_t[[paste0(sx,"_Tm")]] <- get(paste0(cn,"_450m_Tm",sx)) }
else if(exists(paste0(cn,"_450m_Nm",sx))){ cand_n[[paste0(sx,"_Nm")]] <- get(paste0(cn,"_450m_Nm",sx)) } }
if(length(cand_t)==0 && exists(paste0(cn,"_450m_Tm"))) cand_t[["all_T"]] <- get(paste0(cn,"_450m_Tm"))
if(length(cand_n)==0 && exists(paste0(cn,"_450m_Nm"))) cand_n[["all_N"]] <- get(paste0(cn,"_450m_Nm"))
if(length(cand_t)>0 && length(cand_n)>0){ all_m <- c(cand_t, cand_n); com_cg <- top40_cdom_cgs
for(k in 1:length(all_m)){ com_cg <- intersect(com_cg, rownames(all_m[[k]])) }
if(length(com_cg)>=10){ m_tm <- NULL; for(k in 1:length(cand_t)){ sub_m <- as.matrix(cand_t[[k]][com_cg, , drop=FALSE]); m_tm <- if(is.null(m_tm)) sub_m else cbind(m_tm, sub_m) }
m_nm <- NULL; for(k in 1:length(cand_n)){ sub_m <- as.matrix(cand_n[[k]][com_cg, , drop=FALSE]); m_nm <- if(is.null(m_nm)) sub_m else cbind(m_nm, sub_m) }
storage.mode(m_tm) <- "numeric"; storage.mode(m_nm) <- "numeric"
if(ncol(m_tm)>=2 && ncol(m_nm)>=2){ db <- rowMeans(m_tm, na.rm=TRUE) - rowMeans(m_nm, na.rm=TRUE); grp_t <- ifelse(cn %in% c("BRCA","PRAD","LUAD"), "Bone-Tropic", "Non-Tropic"); w_sub <- r_cd38[com_cg]
ok_p <- which(is.finite(w_sub) & is.finite(db))
if(length(ok_p)>=5 && is.finite(var(w_sub[ok_p])) && var(w_sub[ok_p])>0 && is.finite(var(db[ok_p])) && var(db[ok_p])>0){
  ct <- suppressWarnings(cor.test(w_sub[ok_p], db[ok_p]))
  stat38 <- rbind(stat38, data.frame(Cancer=cn, Tropism=grp_t, N_Loci=length(ok_p), Cor_Condyle_Delta_R=ct$estimate, P_val=ct$p.value, Mean_Abs_Delta=mean(abs(db[ok_p])), stringsAsFactors=FALSE))
  cdom_plot_df38 <- rbind(cdom_plot_df38, data.frame(Cancer=cn, Tropism=grp_t, CpG=com_cg[ok_p], Condyle_Dominance_R=w_sub[ok_p], Tumor_Delta_Beta=db[ok_p], stringsAsFactors=FALSE)) } } } } }
if(nrow(stat38)>0){ write.table(stat38, "38_Medial_Lateral_Condyle_Cancer_Shift_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p38a <- ggplot(stat38, aes(x=reorder(Cancer, Cor_Condyle_Delta_R), y=Cor_Condyle_Delta_R, fill=Tropism)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title="Medial-vs-Lateral Condyle Dominance Coupling with Cancer Shift", x="Solid Cancer", y="Correlation (Condyle Asymmetry R vs Cancer Delta Beta)", fill="Tropism")
  p38b <- ggplot(cdom_plot_df38, aes(x=Condyle_Dominance_R, y=Tumor_Delta_Beta, color=Tropism)) + geom_point(alpha=0.6, size=1.5) + geom_smooth(method="lm", formula=y~x, se=FALSE) + facet_wrap(~Tropism) + scale_color_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + labs(title="Condyle Asymmetry Mechanics: Bone-Tropic vs Non-Tropic Shift", x="Primate Condylar Asymmetry Correlation R", y="Tumor Delta Beta")
  ggsave("38_Medial_Lateral_Condyle_Cancer_Shift_Plot.pdf", plot=gridExtra::grid.arrange(p38a, p38b, ncol=2), width=13.5, height=5.5)
} else { cat("경고: 유효한 정상/종양 쌍이 매칭된 암종이 없어 stat38이 비어 있습니다.\n") }

################################################################################################
# Module D39. 개체별 진화적 후성유전 이상치(IEA Score)의 TCGA 고형암 전이율 및 악성도 연결 (수정본)
################################################################################################
iea_res39 <- mat_pr21 - rowMeans(mat_pr21); iea_score_bone39 <- sqrt(colMeans(iea_res39^2, na.rm=TRUE))
r_iea39 <- suppressWarnings(as.numeric(cor(iea_score_bone39, t(mat_pr21)))); names(r_iea39) <- COM_D_CpG; r_iea39[!is.finite(r_iea39)] <- 0
top40_iea_cgs <- COM_D_CpG[order(abs(r_iea39), decreasing=TRUE)[1:40]]; w_iea39 <- r_iea39[top40_iea_cgs]
stat39 <- data.frame(); iea_tumor_df39 <- data.frame(); trop_cans39 <- c("BRCA","PRAD","LUAD")
for(cn in tcga_cans){ cl_n <- paste0(cn, "_c")
if(exists(cl_n)){ c_df <- get(cl_n); cand_m <- list()
for(sx in c("F", "M")){ nt_n <- paste0(cn, "_450m_NT", sx); nt_f <- paste0(cn, "_450m_NT", sx, "f")
if(exists(nt_n) && exists(nt_f)){ d <- get(nt_n); f <- get(nt_f); cand_m[[sx]] <- d[, f == "T", drop=FALSE] } else if(exists(paste0(cn, "_450m_Tm", sx, "2"))){ cand_m[[sx]] <- get(paste0(cn, "_450m_Tm", sx, "2")) } else if(exists(paste0(cn, "_450m_Tm", sx))){ cand_m[[sx]] <- get(paste0(cn, "_450m_Tm", sx)) } }
if(length(cand_m) == 0 && exists(paste0(cn, "_450m_Tm"))){ cand_m[["all"]] <- get(paste0(cn, "_450m_Tm")) }
if(length(cand_m) > 0){ com_cg39 <- top40_iea_cgs; for(k in 1:length(cand_m)){ com_cg39 <- intersect(com_cg39, rownames(cand_m[[k]])) }
if(length(com_cg39) >= 20){ m_tm <- NULL; for(k in 1:length(cand_m)){ sub_m <- as.matrix(cand_m[[k]][com_cg39, , drop=FALSE]); m_tm <- if(is.null(m_tm)) sub_m else cbind(m_tm, sub_m) }
pat_m <- toupper(substr(gsub("[_\\.]", "-", colnames(m_tm)), 1, 12)); pat_c <- toupper(substr(gsub("[_\\.]", "-", c_df$ID), 1, 12)); com_p <- intersect(pat_m, pat_c)
if(length(com_p) >= 5){ idx_m <- match(com_p, pat_m); idx_c <- match(com_p, pat_c); sub_mat <- as.matrix(m_tm[, idx_m, drop=FALSE]); storage.mode(sub_mat) <- "numeric"
for(r in 1:nrow(sub_mat)){ na_idx <- which(is.na(sub_mat[r, ])); if(length(na_idx) > 0){ rm_v <- mean(sub_mat[r, -na_idx], na.rm=TRUE); sub_mat[r, na_idx] <- ifelse(is.finite(rm_v), rm_v, 0.5) } }
w_sub <- w_iea39[com_cg39]; iea_sc <- as.numeric(crossprod(sub_mat, w_sub) / sum(abs(w_sub)))
sub_c <- c_df[idx_c, ]; st_vec <- if("StageC" %in% colnames(sub_c)) as.character(sub_c$StageC) else rep("X", length(idx_c))
surv_vec <- if("Surv" %in% colnames(sub_c)) as.character(sub_c$Surv) else rep("X", length(idx_c))
grp_t <- ifelse(cn %in% trop_cans39, "Bone-Tropic", "Non-Tropic")
iea_tumor_df39 <- rbind(iea_tumor_df39, data.frame(Sample=colnames(sub_mat), Cancer=cn, Tropism=grp_t, StageC=st_vec, Surv=surv_vec, IEA_Score=as.numeric(scale(iea_sc)), stringsAsFactors=FALSE))
vl <- iea_sc[st_vec == "L"]; vm <- iea_sc[st_vec == "M"]; vl <- vl[is.finite(vl)]; vm <- vm[is.finite(vm)]
tt_st <- if(length(vl) >= 3 && length(vm) >= 3 && is.finite(var(vl)) && var(vl) > 0 && is.finite(var(vm)) && var(vm) > 0) suppressWarnings(t.test(vm, vl)$p.value) else NA
m_vl <- if(length(vl) > 0) mean(vl) else NA; m_vm <- if(length(vm) > 0) mean(vm) else NA
stat39 <- rbind(stat39, data.frame(Cancer=cn, Tropism=grp_t, N=length(iea_sc), Mean_IEA=mean(iea_sc, na.rm=TRUE), Delta_Stage_M_L=ifelse(is.finite(m_vl) && is.finite(m_vm), m_vm - m_vl, NA), Stage_P=tt_st, stringsAsFactors=FALSE)) } } } } }
if(nrow(stat39) > 0 && nrow(iea_tumor_df39) > 0){ write.table(stat39, "39_Tumor_Evolutionary_Anomaly_IEA_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
  p39a <- ggplot(iea_tumor_df39, aes(x=reorder(Cancer, IEA_Score, FUN=median), y=IEA_Score, fill=Tropism)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title="Tumor Projection of Evolutionary Epigenetic Anomaly (IEA)", x="Solid Cancer", y="Tumor IEA Score (Z-score)")
  sub_stat39 <- stat39[!is.na(stat39$Delta_Stage_M_L), ]
  p39b <- if(nrow(sub_stat39) > 0){ ggplot(sub_stat39, aes(x=reorder(Cancer, Delta_Stage_M_L), y=Delta_Stage_M_L, fill=!is.na(Stage_P) & Stage_P < 0.05)) + geom_bar(stat="identity", color="black", width=0.7) + coord_flip() + geom_hline(yintercept=0, linetype="dashed") + scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="gray70"), labels=c("TRUE"="P < 0.05","FALSE"="Not Sig")) + theme_bw() + labs(title="IEA Activation in Advanced Stage Tumors (Stage Advanced - Early)", x="Solid Cancer", y="Delta IEA Score", fill="Significance") } else { ggplot() + theme_void() }
  ggsave("39_Tumor_Evolutionary_Anomaly_IEA_Plot.pdf", plot=gridExtra::grid.arrange(p39a, p39b, ncol=2), width=13.5, height=5.5)
} else { cat("경고: 유효한 암종 표본이 매칭되지 않아 stat39가 비어 있습니다.\n") }

################################################################################################
# Module D40. 다단계 인과 경로 아틀라스(진화 기저선 -> 뼈 리모델링/기하학 -> 고형암 발암성 왜곡)
################################################################################################
mods5_40 <- list(M1_Length=c("FL","BFL","MFL"), M2_Diaphysis=c("SSWd","SSDp","MSWd","MSDp","ISWd","ISDp"), M3_Head=c("HH","HL","HW"), M4_Neck=c("ANL","ANH","AND","BNL","PW","LTH","LTN","LGT"), M5_Condyle=c("MCH","MCDp","MCW","LCH","LCDp","LCW","INW","IND","BCW"))
evo_base40 <- as.numeric(scale(rowMeans(sp_m12))); cancer_shift40 <- as.numeric(scale(rowMeans(abs(can_db_mat11))))
stat40 <- data.frame()
for(mn in names(mods5_40)){ m_cgs_bone <- names(sort(rowMeans(abs(cor(t(mat_pr21), z_29_36[, mods5_40[[mn]], drop=FALSE]))), decreasing=TRUE)[1:300])
idx_sub <- match(m_cgs_bone, COM_D_CpG); x_e <- evo_base40[idx_sub]; y_m <- rowMeans(abs(cor(t(mat_pr21[m_cgs_bone, ]), z_29_36[, mods5_40[[mn]], drop=FALSE]))); z_c <- cancer_shift40[idx_sub]
f1 <- summary(lm(y_m ~ x_e))$coef[2, ]; f2 <- summary(lm(z_c ~ y_m + x_e))$coef[2, ]; ind_eff <- f1[1] * f2[1]
stat40 <- rbind(stat40, data.frame(Bone_Module=mn, Step1_Evo_to_Bone_Beta=f1[1], Step1_P=f1[4], Step2_Bone_to_Cancer_Beta=f2[1], Step2_P=f2[4], Indirect_Cascade_Beta=ind_eff, Total_Coupling=abs(ind_eff))) }
write.table(stat40, "40_MultiStage_Causal_Cascade_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
m40_l <- reshape2::melt(stat40[, c("Bone_Module","Step1_Evo_to_Bone_Beta","Step2_Bone_to_Cancer_Beta","Indirect_Cascade_Beta")], id.vars="Bone_Module", variable.name="Causal_Step", value.name="Beta_Value")
p40 <- ggplot(m40_l, aes(x=Causal_Step, y=Bone_Module, fill=Beta_Value)) + geom_tile(color="white") + geom_text(aes(label=round(Beta_Value, 3)), size=3.8, fontface="bold") + scale_fill_gradient2(low="#0571B0", mid="white", high="#CA0020", midpoint=0) + theme_minimal() + theme(axis.text.x=element_text(angle=20, hjust=1, face="bold")) + labs(title="Multi-Stage Causal Path: Primate Evolution -> Bone Morphological Coupling -> Cancer Aberration", subtitle="Step 1: Evolutionary Baseline to Bone Mechanics | Step 2: Bone Mechanics to Cancer Dysregulation | Indirect Cascade", x="Sequential Path Cascade Step", y="5 Femur Anatomical Modules", fill="Path Beta"); ggsave("40_MultiStage_Causal_Cascade_Plot.pdf", plot=p40, width=11, height=6)

################################################################################################
# Module D41. 범영장류 뼈-고형암 통합 시료 간 전유전체 메타-HyPi 가중치 산출 및 종합 맨해튼 플롯
################################################################################################
n_bone41 <- ncol(mat_hb); n_tcga41 <- ncol(can_db_mat11)*15; w_n41 <- c(sqrt(n_bone41), sqrt(n_tcga41))
meta_sign_pi41 <- as.numeric((cbind(hm_hp, rowMeans(can_db_mat11)*5) %*% w_n41) / sqrt(sum(w_n41^2)))
man_df41 <- data.frame(CpG=COM_D_CpG, Meta_CancerBone_Pi=meta_sign_pi41, Abs_Meta_Pi=abs(meta_sign_pi41))
man_df41 <- merge(man_df41, GPL13534[, c("ID","CHR","MAPINFO","UCSC_RefGene_Name","Relation_to_UCSC_CpG_Island")], by.x="CpG", by.y="ID", all.x=TRUE)
man_df41 <- man_df41[man_df41$CHR %in% as.character(1:22), ]; man_df41$CHR <- factor(man_df41$CHR, levels=as.character(1:22)); man_df41$MAPINFO <- as.numeric(as.character(man_df41$MAPINFO))
man_df41 <- man_df41[order(as.numeric(as.character(man_df41$CHR)), man_df41$MAPINFO), ]; man_df41$Pos_Index <- 1:nrow(man_df41)
write.table(head(man_df41[order(man_df41$Abs_Meta_Pi, decreasing=TRUE), ], 100), "41_Grand_Meta_CancerBone_Top100_CpGs.txt", sep="\t", quote=FALSE, row.names=FALSE)
man_sub41 <- man_df41[unique(c(which(man_df41$Abs_Meta_Pi > 2.5), seq(1, nrow(man_df41), length.out=min(8000, nrow(man_df41))))), ]
p41 <- ggplot(man_sub41, aes(x=Pos_Index, y=Meta_CancerBone_Pi, color=CHR)) + geom_point(alpha=0.65, size=0.8) + geom_hline(yintercept=c(-3, 3), linetype="dashed", color="red") + geom_hline(yintercept=0, color="black") + theme_bw() + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(title="Grand Pan-Bone & TCGA Solid Cancer Integrated Manhattan Plot (Signed Meta-HyPi)", subtitle="Autosomes 1 to 22; Red dashed lines indicate |Meta-HyPi| = 3", x="Genomic Position Across Autosomes (Chr 1 - 22)", y="Sample-Size Weighted Meta-HyPi"); ggsave("41_Grand_Meta_Manhattan_Plot.pdf", plot=p41, width=12.5, height=5.5)

################################################################################################
# Module D42. 골 전이 호발암(BRCA, PRAD, LUAD) 특이적 Osteomimicry 점수 기반 다중 코호트 성과 평가
################################################################################################
w_osteo42 <- (hm_hp * sign(rowMeans(can_db_mat11))); top50_osteo_cgs <- COM_D_CpG[order(abs(w_osteo42), decreasing=TRUE)[1:50]]
osteo_perf42 <- data.frame(); for(cn in tcga_cans){ tm_n <- paste0(cn,"_450m_Tm")
if(exists(tm_n)){ m_tm <- as.matrix(get(tm_n)[top50_osteo_cgs, , drop=FALSE]); storage.mode(m_tm) <- "numeric"
sc_o <- as.numeric(crossprod(m_tm, w_osteo42[top50_osteo_cgs])/sum(abs(w_osteo42[top50_osteo_cgs])))
osteo_perf42 <- rbind(osteo_perf42, data.frame(Sample=colnames(m_tm), Cancer=cn, Tropism=ifelse(cn %in% c("BRCA","PRAD","LUAD"),"Bone-Tropic","Non-Tropic"), Osteomimicry_Score=as.numeric(scale(sc_o)))) } }
stat42 <- aggregate(Osteomimicry_Score ~ Cancer + Tropism, data=osteo_perf42, FUN=function(x){ c(Mean=mean(x), SD=sd(x), Median=median(x)) }); stat42 <- data.frame(stat42[,1:2], stat42$Osteomimicry_Score)
tt42 <- t.test(osteo_perf42$Osteomimicry_Score[osteo_perf42$Tropism=="Bone-Tropic"], osteo_perf42$Osteomimicry_Score[osteo_perf42$Tropism=="Non-Tropic"])
write.table(stat42, "42_Osteomimicry_Scorecard_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
p42 <- ggplot(osteo_perf42, aes(x=reorder(Cancer, Osteomimicry_Score, FUN=median), y=Osteomimicry_Score, fill=Tropism)) + geom_boxplot(alpha=0.8, outlier.size=0.3) + geom_hline(yintercept=0, linetype="dashed", color="black") + scale_fill_manual(values=c("Bone-Tropic"="#E41A1C","Non-Tropic"="#377EB8")) + theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1)) + labs(title=paste0("Multi-Cancer Epigenetic Osteomimicry Scorecard (T-test P = ", signif(tt42$p.value, 3), ")"), x="Solid Cancer", y="Standardized Osteomimicry Score (Z)", fill="Bone Tropism"); ggsave("42_Osteomimicry_Scorecard_Plot.pdf", plot=p42, width=10, height=5.5)

################################################################################################
# Module D43. 5대 핵심 축(진화·노화·역학기하·관절병증·고형암 발암성) 통합 35대 마스터 좌위 아틀라스
################################################################################################
ax1_evo43 <- as.numeric(scale(evo_sd12)); ax2_age43 <- as.numeric(scale(abs(hm_hp))); ax3_bio43 <- as.numeric(scale(abs(r_j22)))
ax4_oa43 <- as.numeric(scale(abs(db_oa))); ax5_can43 <- as.numeric(scale(rowMeans(abs(can_db_mat11))))
mat5_43 <- cbind(Evo_Divergence=ax1_evo43, Bone_Aging=ax2_age43, Mechanical_Rigidity=ax3_bio43, Arthropathy_OA=ax4_oa43, Solid_Tumor_Shift=ax5_can43); rownames(mat5_43) <- COM_D_CpG; mat5_43[!is.finite(mat5_43)] <- 0
master_sc43 <- rowMeans(mat5_43) + 0.4*rowSums(mat5_43 > 1.0); top35_m_cgs <- COM_D_CpG[order(master_sc43, decreasing=TRUE)[1:35]]
g_sym43 <- sub(";.*","",as.character(GPL13534$UCSC_RefGene_Name[match(top35_m_cgs, GPL13534$ID)])); lbl35_43 <- ifelse(is.na(g_sym43)|g_sym43=="", top35_m_cgs, paste0(g_sym43," (",top35_m_cgs,")"))
df43_out <- data.frame(CpG=top35_m_cgs, Gene=g_sym43, Master_Score=master_sc43[top35_m_cgs], mat5_43[top35_m_cgs, ]); write.table(df43_out, "43_Grand_Master_5Axes_Synthesis_Stats.txt", sep="\t", quote=FALSE, row.names=FALSE)
m43_l <- data.frame(Locus=rep(factor(lbl35_43, levels=rev(lbl35_43)), 5), Axis=rep(factor(colnames(mat5_43), levels=colnames(mat5_43)), each=35), Z_Score=as.vector(pmin(4, pmax(-2, mat5_43[top35_m_cgs, ]))))
p43 <- ggplot(m43_l, aes(x=Axis, y=Locus, fill=Z_Score)) + geom_tile(color="white") + geom_text(aes(label=round(Z_Score, 2)), size=2.8, fontface="bold") + scale_fill_viridis_c(option="turbo") + theme_minimal() + theme(axis.text.x=element_text(angle=25, hjust=1, face="bold")) + labs(title="Grand Master 5-Axis Synthesis (Primate Evolution, Bone Aging, Biomechanics, Arthropathy & Cancer)", subtitle="Top 35 Pleiotropic Master Epigenetic Loci Across All 7 GEO Bone Datasets and TCGA 10 Solid Cancers", x="Core Biological Axis", y="Master Gene (CpG Probe)", fill="Z-Score"); ggsave("43_Grand_Master_5Axes_Synthesis_Plot.pdf", plot=p43, width=11, height=9)

################################################################################################
# Module D44. 전체 통합 코호트 특성 요약표(Table 1) 및 핵심 마스터 유전자 아틀라스(Table 2) 반출
################################################################################################
tb1_44 <- data.frame(); for(coh in unique(meta_d01$Dataset)){ sub_m <- meta_d01[meta_d01$Dataset==coh, ]; ag <- sub_m$Age[!is.na(sub_m$Age)]
tb1_44 <- rbind(tb1_44, data.frame(Dataset=coh, Organ=sub_m$Organ[1], Condition=sub_m$Condition[1], Total_N=nrow(sub_m), Male_N=sum(sub_m$Sex=="M"), Female_N=sum(sub_m$Sex=="F"), Mean_Age=round(mean(ag),1), SD_Age=round(sd(ag),1))) }
write.table(tb1_44, "44_Table1_Integrated_Cohorts_Summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
tb2_44 <- merge(df43_out[1:30, ], GPL13534[, c("ID","CHR","MAPINFO","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")], by.x="CpG", by.y="ID", all.x=TRUE)
write.table(tb2_44[order(tb2_44$Master_Score, decreasing=TRUE), ], "44_Table2_Top30_Master_Genes.txt", sep="\t", quote=FALSE, row.names=FALSE)
p44 <- ggplot(tb2_44, aes(x=reorder(Gene, Master_Score), y=Master_Score, size=Solid_Tumor_Shift, color=Bone_Aging)) + geom_point(alpha=0.85) + coord_flip() + scale_color_gradient2(low="#377EB8", mid="white", high="#E41A1C", midpoint=1) + theme_bw() + labs(title="Top 30 Master Genes: Overall Synthesis Score vs Cancer & Bone Alteration", x="Master Gene Symbol", y="Composite Master Score", size="Cancer Shift (Z)", color="Bone Aging (Z)"); ggsave("44_Master_Gene_Bubble_Plot.pdf", plot=p44, width=9, height=7.5)

################################################################################################
# Module D45. 세션 정보 기록, 출력 파일 무결성 점검 및 워크스페이스 정리 (Final Integrity Check)
################################################################################################
sink("45_Session_Info.txt"); print(sessionInfo()); sink()
out_files45 <- list.files(pattern="^([0-9]{2}|D[0-9]{2}).*\\.(txt|pdf)$")
chk45 <- data.frame(File_Name=out_files45, Size_Bytes=file.info(out_files45)$size, Status=ifelse(file.info(out_files45)$size > 0, "Valid (Non-Empty)", "Corrupt (Zero-Byte)"))
write.table(chk45, "45_Final_Pipeline_Integrity_Report.txt", sep="\t", quote=FALSE, row.names=FALSE)
cat("=== Module D 통합 분석 완료: 총", nrow(chk45), "개 통계 및 그래픽 파일이 정상 생성되었습니다. ===\n")
gc()


