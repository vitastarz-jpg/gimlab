library(data.table); library(dplyr); library(readxl)

GSE223748m <- as.data.frame(fread("G:/내 드라이브/Rrunning/GEO/HW014/GSE223748/GSE223748_datBetaNormalized.csv.gz", header = TRUE, sep = ","))
rownames(GSE223748m) <- GSE223748m$V1; GSE223748m <- GSE223748m[,-1]

GPL28271 <- as.data.frame(fread("G:/내 드라이브/Rrunning/Annotations/21_Methyl/GPL28271-57075.txt.gz", header=TRUE, skip=65, sep="\t"))[,-c(4:8)]
head(GPL28271); dim(GPL28271)

GSE223748i <- as.data.frame(fread("G:/내 드라이브/Rrunning/GEO/HW014/GSE223748/GSE223748_series_matrix.txt.gz", header = TRUE, sep = "\t", skip=381)); colnames(GSE223748i)[1] <- "ID"
GSE223748i <- dplyr::filter(GSE223748i, grepl("Sample_source_name", ID) | grepl("organism", ID) | grepl("characteristics", ID))
GSE223748i <- as.data.frame(t(as.matrix(GSE223748i[,-1])))
colnames(GSE223748i)[1:9] <- c("V01","V02","V03","V04","V05","V06","V07","V08","V09")
head(GSE223748i,2); dim(GSE223748i)
colnames(GSE223748m) <- rownames(GSE223748i)

table(substr(GSE223748i$V01,1, 5))
table(substr(GSE223748i$V02,1, 5))
table(substr(GSE223748i$V03,1, 3))
table(substr(GSE223748i$V04,1,22)) #2 rows were pushed to 4th column: needs to be corrected

dplyr::filter(GSE223748i, substr(V04,1,6)=="female")
FILTER_C01 <- which(grepl("female", GSE223748i$V04))
GSE223748i[FILTER_C01, ]
GSE223748i[FILTER_C01, 5:11] <- GSE223748i[FILTER_C01, 4:10]
GSE223748i[FILTER_C01, 4] <- "confidenceinageestima: "
GSE223748i[FILTER_C01, ]
table(substr(GSE223748i$V04,1,23)) #4th column is corrected

table(substr(GSE223748i$V05,1, 6)) #96 rows were pushed to 5th column: needs to be corrected

dplyr::filter(GSE223748i, substr(V05,1,6)=="tissue")[1:2,]
table(dplyr::filter(GSE223748i, substr(V05,1,6)=="tissue")$V01)
table(dplyr::filter(GSE223748i, substr(V05,1,6)=="tissue")$V05) #Maybe blood extracted from ear?
FILTER_C02 <- which(grepl("tissue", GSE223748i$V05))
table(GSE223748i[FILTER_C02, ]$V05)
GSE223748i[FILTER_C02, 5:11] <- GSE223748i[FILTER_C02, c(6,8,9,10,11,14,15)]
GSE223748i[FILTER_C02, 12:15] <- ""
table(substr(GSE223748i$V05,1, 6)) #5th column is corrected

table(substr(GSE223748i$V06,1,18))
table(substr(GSE223748i$V07,1,16))
table(substr(GSE223748i$V08,1,19))
table(substr(GSE223748i$V09,1,19))
table(substr(GSE223748i$V10,1,19))
table(substr(GSE223748i$V11,1,19))

GSE223748i <- GSE223748i[,1:11]
head(GSE223748i,2)
GSE223748i$V03 <- round(abs(as.numeric(gsub("age\\: ", "", GSE223748i$V03))),2)
GSE223748i$V04 <- gsub("confidenceinageestimate\\: ",      "", GSE223748i$V04)
GSE223748i$V04 <- gsub("confidenceinageestima\\: ",        "", GSE223748i$V04)
GSE223748i$V05 <- gsub("female\\: ",                       "", GSE223748i$V05)
GSE223748i$V06 <- gsub("speciescommonname\\: ",            "", GSE223748i$V06)
GSE223748i$V07 <- gsub("maximumlifespanyears\\: ",         "", GSE223748i$V07)
GSE223748i$V07 <- gsub("maximum_age\\: ",                  "", GSE223748i$V07)
GSE223748i$V08 <- gsub("average_adultweight\\: ",          "", GSE223748i$V08)
GSE223748i$V08 <- gsub("average_weight\\: ",               "", GSE223748i$V08)
GSE223748i$V09 <- gsub("ageatsexualmaturityyears\\: ",     "", GSE223748i$V09)
GSE223748i$V09 <- gsub("age\\.sexualmaturity\\: ",         "", GSE223748i$V09)
GSE223748i$V10 <- gsub("networktrainingset\\: ",           "", GSE223748i$V10)
GSE223748i$V11 <- gsub("panmammalianclocktrainingset\\: ", "", GSE223748i$V11)
head(GSE223748i)

table(substr(GSE223748i$V01,1,3))
table(substr(GSE223748i$V02,1,3))
table(substr(GSE223748i$V03,1,3))
table(substr(GSE223748i$V04,1,3))
table(substr(GSE223748i$V05,1,3))
table(substr(GSE223748i$V06,1,3))
table(substr(GSE223748i$V07,1,3))
table(substr(GSE223748i$V08,1,3))
table(substr(GSE223748i$V09,1,3))
table(substr(GSE223748i$V10,1,3))
table(substr(GSE223748i$V11,1,3))

GSE223748i[,c(3:4,7:9)] <- lapply(GSE223748i[,c(3:4,7:9)], as.numeric)
GSE223748i[,c(1:2,5:6,10:11)] <- lapply(GSE223748i[,c(1:2,5:6,10:11)], as.factor)
GSE223748i[is.na(GSE223748i)] <- ""
colnames(GSE223748i) <- c("Tissue","SciName","Age","Age_Estimate","Sex","ComName","MaxSpan","Weight","Age_Sexual","Net_train","Clock_train")
head(GSE223748i)
str(GSE223748i)

table(rownames(GSE223748i) == colnames(GSE223748m))

GSE223748m[1:5, 1:5]; dim(GSE223748m) #37554,15043
GSE223748i[1:5, 1:5]; dim(GSE223748i) #15043,11
GSE223748i$Age <- as.numeric(GSE223748i$Age)
GSE223748i$ID <- rownames(GSE223748i)

head(GSE223748i); dim(GSE223748i); length(unique(GSE223748i$SciName)) #15043,12 / 347 species
hist(GSE223748i$Age)
table(GSE223748i$Sex)
table(GSE223748i$SciName)
table(unique(GSE223748i$SciName) %in% "Taurotragus oryx")

Species_DF <- as.data.frame(table(GSE223748i$SciName))
write.table(Species_DF, file="clipboard",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

GSE223748i_add <- as.data.frame(fread("Input_20250415/GSE223748i_add.txt", header = TRUE, sep = "\t"))
head(GSE223748i_add); dim(GSE223748i_add) #15043,22

GSE223748i_add <- merge(x=GSE223748i_add, y=GSE223748i, by="SciName", all.y=TRUE); rownames(GSE223748i_add) <- GSE223748i_add$ID
GSE223748i_add$Sex <- ifelse(GSE223748i_add$Sex == 0, "M", ifelse(GSE223748i_add$Sex == 1, "F", "X"))
head(GSE223748i_add); dim(GSE223748i_add)
table(duplicated(GSE223748i_add$ID))
sort(table(GSE223748i_add$Tissue))
write.table(sort(table(GSE223748i_add$Tissue)), file="clipboard",sep="\t",quote=FALSE,row.names=FALSE)
