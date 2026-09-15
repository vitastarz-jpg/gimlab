rm(list=ls())
setwd("G:/내 드라이브/Rrunning/KW21")

#Load library
library(dplyr); library(data.table); library(matrixStats); library(pheatmap); library(progress); library(readxl); library(ggplot2); library(ggpubr); library(rlang); library(gridExtra); library(ggvenn); library(Hmisc)

annotation_colors = list(
  Cohort  = c("Guro"="#1DCD9F", "Anam"="#FE4F2D"),
  Group     = c("G1"="#E9A5F1", "G2"="#C68EFD", "G3"="#8F87F1", "G4"="#F8F8E1", "G5"="#F8F8E1"),
  CDH     = c("N"="#A7D477", "Y"="#F72C5B", "X"="#4E1F00"),
  Morphology     = c("Lobular"="#A7D477", "Mucinous"="#F72C5B", "Ductal"="#FEBA17", "Normal"="#F8F8E1"),
  IHC     = c("Positive"="#A6D6D6", "Membranous"="#F7CFD8", "Loss"="#CA0020", "X"="#4E1F00"),
  FC     = c("28"="#A7D477", "22"="#F72C5B"),
  logPV  = c("28"="#FCF596", "22"="#FF4545"))

my_palette1 <- colorRampPalette(c("#60C5F1", "#FFFCCC", "#F47378"))(500) #Cellline
my_palette2 <- colorRampPalette(c("#4DAC26", "#FFFCCC", "#CA0020"))(500) #Blood
my_palette3 <- colorRampPalette(c("#0571B0", "#FFFCCC", "#D01C8B"))(500) #Tissue
my_palette4 <- colorRampPalette(c("#018571", "#FFFCCC", "#E66101"))(500) #PDx
my_palette5 <- colorRampPalette(c("#FFFCCC", "#CA0020"))(500) #PDx
make.unique.2 = function(x, sep='.'){ ave(x, x, FUN=function(a){if(length(a) > 1){paste(a, 1:length(a), sep=sep)} else {a}}) }
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
##########################################################################################################################################################

library(gwasrapidd)
Sys.time()
asc_rs <- get_associations(variant_id=c("rs1815739", "rs7412"))
trt_rs <- asc_rs@traits
write.table(trt_rs, "Output_Evol/06_Traits_from_RS.txt", sep="\t", row.names=FALSE, quote=FALSE)

Sys.time()
var_gn <- get_variants(gene_name="TP53"); asc_gn <- get_associations(variant_id=var_gn@variants$variant_id[1:10]); trt_gn <- asc_gn@traits
write.table(trt_gn, "Output_Evol/06_Traits_from_Gene.txt", sep="\t", row.names=FALSE, quote=FALSE)

Sys.time()
asc_tr <- get_associations(efo_trait="breast cancer"); rs_tr <- asc_tr@risk_alleles; gn_tr <- asc_tr@genes
write.table(merge(rs_tr, gn_tr, by="association_id", all=TRUE), "Output_Evol/06_RS_Gene_from_Trait.txt", sep="\t", row.names=FALSE, quote=FALSE)

Sys.time()
trt_cnt <- as.data.frame(table(trt_rs$disease_trait))
pdf("Output_Evol/06_GWAS_Summary.pdf", width=10, height=8); print(ggplot(trt_cnt[order(-trt_cnt$Freq)[1:min(15, nrow(trt_cnt))], ], aes(x=reorder(Var1, Freq), y=Freq)) + geom_bar(stat="identity", fill="darkred", alpha=0.8) + coord_flip() + labs(title="Top GWAS Traits Associated with Target SNPs", x="Disease / Trait", y="Association Count") + theme_minimal(base_size=12) + theme(axis.text.y=element_text(size=10, face="bold"))); dev.off()

Sys.time()



##########################################################################################################################################################
#Evolution
library(bigsnpr)
bed_file <- "D:/KoGES/Geno/Aff/2024-045_KARE_affy5_8840_QCed.bed"; rds_file <- sub("\\.bed$", ".rds", bed_file) #5s
if (file.exists(rds_file)) {KoGES_Af5 <- snp_attach(rds_file)} else {rds_path_new <- snp_readBed(bed_file); KoGES_Af5 <- snp_attach(rds_path_new)}

bed_file <- "D:/KoGES/Geno/Aff/2024-045_HEXA_affy6_3693_QCed.bed"; rds_file <- sub("\\.bed$", ".rds", bed_file) #5s
if (file.exists(rds_file)) {KoGES_Af6 <- snp_attach(rds_file)} else {rds_path_new <- snp_readBed(bed_file); KoGES_Af6 <- snp_attach(rds_path_new)}

bed_file <- "D:/KoGES/Geno/Exo/2024-045_ALL_Exomechip_14025_QCed_Rare.bed"; rds_file <- sub("\\.bed$", ".rds", bed_file) #5s
if (file.exists(rds_file)) {KoGES_Exo <- snp_attach(rds_file)} else {rds_path_new <- snp_readBed(bed_file); KoGES_Exo <- snp_attach(rds_path_new)}

bed_file <- "D:/KoGES/Geno/Kch/2024-045_ALL_KCHIP_72291_QCed.bed"; rds_file <- sub("\\.bed$", ".rds", bed_file) #5m
if (file.exists(rds_file)) {KoGES_Kch <- snp_attach(rds_file)} else {rds_path_new <- snp_readBed(bed_file); KoGES_Kch <- snp_attach(rds_path_new)}

head(KoGES_Af5$map); head(KoGES_Af5$fam); KoGES_Af5$genotypes[1:5,1:5]
head(KoGES_Af6$map); head(KoGES_Af6$fam); KoGES_Af6$genotypes[1:5,1:5]
head(KoGES_Exo$map); head(KoGES_Exo$fam); KoGES_Exo$genotypes[1:5,1:5]
head(KoGES_Kch$map); head(KoGES_Kch$fam); KoGES_Kch$genotypes[1:5,1:5]
dim(KoGES_Af5$map); dim(KoGES_Af5$fam); dim(KoGES_Af5$genotypes)
dim(KoGES_Af6$map); dim(KoGES_Af6$fam); dim(KoGES_Af6$genotypes)
dim(KoGES_Exo$map); dim(KoGES_Exo$fam); dim(KoGES_Exo$genotypes)
dim(KoGES_Kch$map); dim(KoGES_Kch$fam); dim(KoGES_Kch$genotypes)
KoGES_Af5map <- KoGES_Af5$map
KoGES_Af6map <- KoGES_Af6$map
KoGES_Exomap <- KoGES_Exo$map
KoGES_Kchmap <- KoGES_Kch$map

CNT_0 <- colSums(KoGES_Af5$genotypes[] == 0, na.rm=TRUE)
CNT_1 <- colSums(KoGES_Af5$genotypes[] == 1, na.rm=TRUE)
CNT_2 <- colSums(KoGES_Af5$genotypes[] == 2, na.rm=TRUE)
REVERSED <- sum(CNT_2 > CNT_0, na.rm=TRUE)
print(paste("2가 0보다 많은(Alt가 Major인) SNP 개수:", REVERSED))


Sys.time()
library(BSgenome.Hsapiens.UCSC.hg19)



MAP_ENH <- as.data.frame(data.table::fread("Output_COM_rsNo/KoGES_Af5map.txt.gz", header=TRUE, sep="\t")); CUR_RS <- MAP_ENH$rs_number[match(KoGES_Af5map$marker.ID, MAP_ENH$marker.ID)]
WGT <- KoGES_Af5mapx$PhastCons_Score[match(CUR_RS, KoGES_Af5mapx$refsnp_id)]; W_C <- ifelse(!is.na(WGT) & WGT >= 0.8, WGT, 0); W_H <- ifelse(!is.na(WGT) & WGT <= 0.2, 1 - WGT, 0)
P_C <- numeric(nrow(KoGES_Af5$genotypes)); P_H <- numeric(nrow(KoGES_Af5$genotypes))
for(I in 1:nrow(KoGES_Af5$genotypes)) { R <- KoGES_Af5$genotypes[I, ]; R[FLP] <- 2 - R[FLP]; P_C[I] <- sum(R * W_C, na.rm=TRUE); P_H[I] <- sum(R * W_H, na.rm=TRUE) }
KoGES_Af5_PRS <- data.frame(DIST_ID=KoGES_Af5$fam$sample.ID, PRS_CON=P_C, PRS_HUM=P_H); write.table(KoGES_Af5_PRS, "Output_Evol/07_PRS_ALIGNED.txt", sep="\t", row.names=FALSE, quote=FALSE)



KoGES_Af5mapx <- as.data.frame(data.table::fread("Output_COM_rsNo/KoGES_Af5mapx.txt.gz", header=TRUE, sep="\t"))
BSG <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19; CHR <- paste0("chr", KoGES_Af5map$chromosome)
VAL <- CHR %in% GenomeInfoDb::seqnames(BSG) & !is.na(KoGES_Af5map$physical.pos) & KoGES_Af5map$physical.pos > 0; VAL[VAL] <- KoGES_Af5map$physical.pos[VAL] <= GenomeInfoDb::seqlengths(BSG)[CHR[VAL]]
GR <- GenomicRanges::GRanges(seqnames=CHR[VAL], ranges=IRanges::IRanges(start=KoGES_Af5map$physical.pos[VAL], width=1)); REF <- as.character(Biostrings::getSeq(BSG, GR))
MAP <- KoGES_Af5map; MAP$REF <- NA; MAP$REF[VAL] <- REF; FLP <- which(MAP$allele1 == MAP$REF); print(paste("Alignment Correction: Flipping", length(FLP), "SNPs"))
WGT <- KoGES_Af5mapx$PhastCons_Score[match(KoGES_Af5map$rs_number, KoGES_Af5mapx$refsnp_id)]; P_C <- numeric(nrow(KoGES_Af5$genotypes)); P_H <- numeric(nrow(KoGES_Af5$genotypes))
W_C <- ifelse(!is.na(WGT) & WGT >= 0.8, WGT, 0); W_H <- ifelse(!is.na(WGT) & WGT <= 0.2, 1 - WGT, 0)
for(I in 1:nrow(KoGES_Af5$genotypes)) { R <- KoGES_Af5$genotypes[I, ]; R[FLP] <- 2 - R[FLP]; P_C[I] <- sum(R * W_C, na.rm=TRUE); P_H[I] <- sum(R * W_H, na.rm=TRUE) }
write.table(data.frame(DIST_ID=KoGES_Af5$fam$sample.ID, PRS_CON=P_C, PRS_HUM=P_H), "Output_Evol/07_PRS_ALIGNED.txt", sep="\t", row.names=FALSE, quote=FALSE)

KoGES_Af5_PRS <- data.frame(DIST_ID=KoGES_Af5$fam$sample.ID, PRS_CON=P_C, PRS_HUM=P_H)
head(KoGES_Af5_PRS)
summary(KoGES_Af5_PRS)

hist(KoGES_Af5_PRS$PRS_HUM)

Sys.time()


table(KoGES_Af5$genotypes[,1])


#Takes 10mins
#dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE); .libPaths(new = c(Sys.getenv("R_LIBS_USER"), .libPaths()))
#install.packages(pkgs = "G:/내 드라이브/Downloads/SNPlocs.Hsapiens.dbSNP144.GRCh37_0.99.20.tar.gz", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))
library(data.table); library(GenomicRanges); library(rtracklayer); library(GenomicFeatures); library(org.Hs.eg.db); library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

#download.file(url = "https://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz", destfile = paste0(tempdir(), "/hg18ToHg19.over.chain.gz"), mode = "wb")
#writeLines(text = readLines(con = gzfile(paste0(tempdir(), "/hg18ToHg19.over.chain.gz"))), con = paste0(tempdir(), "/hg18ToHg19.over.chain")); chain_hg18_19 <- rtracklayer::import.chain(paste0(tempdir(), "/hg18ToHg19.over.chain"))
#dt_lo <- data.table::as.data.table(unlist(rtracklayer::liftOver(x = GenomicRanges::GRanges(seqnames = paste0("chr", KoGES_Af5map$chromosome), ranges = IRanges::IRanges(start = KoGES_Af5map$physical.pos, end = KoGES_Af5map$physical.pos), marker.ID = KoGES_Af5map$marker.ID), chain = chain_hg18_19)))
#KoGES_Af5map <- as.data.frame(data.table::as.data.table(KoGES_Af5map)[dt_lo, `:=`(chromosome = gsub("chr", "", as.character(i.seqnames)), physical.pos = i.start), on = .(marker.ID == marker.ID)])
#dt_lo <- data.table::as.data.table(unlist(rtracklayer::liftOver(x = GenomicRanges::GRanges(seqnames = paste0("chr", KoGES_Af6map$chromosome), ranges = IRanges::IRanges(start = KoGES_Af6map$physical.pos, end = KoGES_Af6map$physical.pos), marker.ID = KoGES_Af6map$marker.ID), chain = chain_hg18_19)))
#KoGES_Af6map <- as.data.frame(data.table::as.data.table(KoGES_Af6map)[dt_lo, `:=`(chromosome = gsub("chr", "", as.character(i.seqnames)), physical.pos = i.start), on = .(marker.ID == marker.ID)])
#
#tx_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene; g_hg19 <- genes(tx_hg19)
#sym_hg19 <- mapIds(org.Hs.eg.db, keys = g_hg19$gene_id, column = "SYMBOL", keytype = "ENTREZID")
#dt_g19 <- data.table(chromosome = gsub("chr", "", as.character(seqnames(g_hg19))), start_position = start(g_hg19), end_position = end(g_hg19), SYMBOL = sym_hg19)
#snps_hg19 <- snpsBySeqname(SNPlocs.Hsapiens.dbSNP144.GRCh37::SNPlocs.Hsapiens.dbSNP144.GRCh37, as.character(c(1:22, "X", "Y"))); dt_snp19 <- data.table(chromosome = gsub("chr", "", as.character(seqnames(snps_hg19))), physical.pos = pos(snps_hg19), rs_number = snps_hg19$RefSNP_id)
#KoGES_Af5map$chromosome <- as.character(KoGES_Af5map$chromosome); KoGES_Af5map <- as.data.frame(data.table::as.data.table(KoGES_Af5map)[dt_g19, SYMBOL := i.SYMBOL, on = .(chromosome == chromosome, physical.pos >= start_position, physical.pos <= end_position)][dt_snp19, rs_number := i.rs_number, on = .(chromosome == chromosome, physical.pos == physical.pos)])
#KoGES_Af6map$chromosome <- as.character(KoGES_Af6map$chromosome); KoGES_Af6map <- as.data.frame(data.table::as.data.table(KoGES_Af6map)[dt_g19, SYMBOL := i.SYMBOL, on = .(chromosome == chromosome, physical.pos >= start_position, physical.pos <= end_position)][dt_snp19, rs_number := i.rs_number, on = .(chromosome == chromosome, physical.pos == physical.pos)])
#KoGES_Exomap$chromosome <- as.character(KoGES_Exomap$chromosome); KoGES_Exomap <- as.data.frame(data.table::as.data.table(KoGES_Exomap)[dt_g19, SYMBOL := i.SYMBOL, on = .(chromosome == chromosome, physical.pos >= start_position, physical.pos <= end_position)][dt_snp19, rs_number := i.rs_number, on = .(chromosome == chromosome, physical.pos == physical.pos)])
#KoGES_Kchmap$chromosome <- as.character(KoGES_Kchmap$chromosome); KoGES_Kchmap <- as.data.frame(data.table::as.data.table(KoGES_Kchmap)[dt_g19, SYMBOL := i.SYMBOL, on = .(chromosome == chromosome, physical.pos >= start_position, physical.pos <= end_position)][dt_snp19, rs_number := i.rs_number, on = .(chromosome == chromosome, physical.pos == physical.pos)])
#KoGES_Af5map$gencol <- ifelse(is.na(KoGES_Af5map$rs_number),KoGES_Af5map$marker.ID,KoGES_Af5map$rs_number)
#KoGES_Af6map$gencol <- ifelse(is.na(KoGES_Af6map$rs_number),KoGES_Af6map$marker.ID,KoGES_Af6map$rs_number)
#KoGES_Exomap$gencol <- ifelse(is.na(KoGES_Exomap$rs_number),KoGES_Exomap$marker.ID,KoGES_Exomap$rs_number)
#KoGES_Kchmap$gencol <- ifelse(is.na(KoGES_Kchmap$rs_number),KoGES_Kchmap$marker.ID,KoGES_Kchmap$rs_number)
#write.table(KoGES_Af5map, file="Output_COM_rsNo/KoGES_Af5map.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(KoGES_Af6map, file="Output_COM_rsNo/KoGES_Af6map.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(KoGES_Exomap, file="Output_COM_rsNo/KoGES_Exomap.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(KoGES_Kchmap, file="Output_COM_rsNo/KoGES_Kchmap.txt",sep="\t",quote=FALSE,row.names=FALSE)
KoGES_Af5map <- as.data.frame(fread("Output_COM_rsNo/KoGES_Af5map.txt.gz", header = TRUE, sep = "\t"))
KoGES_Af6map <- as.data.frame(fread("Output_COM_rsNo/KoGES_Af6map.txt.gz", header = TRUE, sep = "\t"))
KoGES_Exomap <- as.data.frame(fread("Output_COM_rsNo/KoGES_Exomap.txt.gz", header = TRUE, sep = "\t"))
KoGES_Kchmap <- as.data.frame(fread("Output_COM_rsNo/KoGES_Kchmap.txt.gz", header = TRUE, sep = "\t"))
head(KoGES_Af5map); str(KoGES_Af5map); dim(KoGES_Af5map) #352211,9
head(KoGES_Af6map); str(KoGES_Af6map); dim(KoGES_Af6map) #626716,9
head(KoGES_Exomap); str(KoGES_Exomap); dim(KoGES_Exomap) # 42730,9
head(KoGES_Kchmap); str(KoGES_Kchmap); dim(KoGES_Kchmap) #447926,9

library(biomaRt); library(GenomicRanges); library(GenomicScores); library(phastCons100way.UCSC.hg38); library(rtracklayer)
#snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
#snp_mart <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp")

#20m / 20m / 3m / 42m
#KoGES_Af5mapx <- KoGES_Af5map[!is.na(KoGES_Af5map$rs_number), ]; nrow(KoGES_Af5map); nrow(KoGES_Af5mapx)
#KoGES_Af6mapx <- KoGES_Af6map[!is.na(KoGES_Af6map$rs_number), ]; nrow(KoGES_Af6map); nrow(KoGES_Af6mapx)
#KoGES_Exomapx <- KoGES_Exomap[!is.na(KoGES_Exomap$rs_number), ]; nrow(KoGES_Exomap); nrow(KoGES_Exomapx)
#KoGES_Kchmapx <- KoGES_Kchmap[!is.na(KoGES_Kchmap$rs_number), ]; nrow(KoGES_Kchmap); nrow(KoGES_Kchmapx)
#snp_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"), filters = "snp_filter", values = KoGES_Af5mapx$rs_number, mart = snp_mart); snp_pos <- snp_pos[snp_pos$chr_name %in% as.character(1:22), ]; gr <- GRanges(seqnames = paste0("chr", snp_pos$chr_name), ranges = IRanges(start = snp_pos$chrom_start, width = 1)); scores <- gscores(phastCons100way.UCSC.hg38, gr); snp_pos$PhastCons_Score <- scores$default; snp_pos$Level <- ifelse(is.na(snp_pos$PhastCons_Score), "Unknown", ifelse(snp_pos$PhastCons_Score >= 0.9, "Highly Conserved (Vertebrates)", ifelse(snp_pos$PhastCons_Score >= 0.3, "Moderately Conserved (Mammals)", "Poorly Conserved (Human-specific)"))); KoGES_Af5mapx <- snp_pos
#snp_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"), filters = "snp_filter", values = KoGES_Af6mapx$rs_number, mart = snp_mart); snp_pos <- snp_pos[snp_pos$chr_name %in% as.character(1:22), ]; gr <- GRanges(seqnames = paste0("chr", snp_pos$chr_name), ranges = IRanges(start = snp_pos$chrom_start, width = 1)); scores <- gscores(phastCons100way.UCSC.hg38, gr); snp_pos$PhastCons_Score <- scores$default; snp_pos$Level <- ifelse(is.na(snp_pos$PhastCons_Score), "Unknown", ifelse(snp_pos$PhastCons_Score >= 0.9, "Highly Conserved (Vertebrates)", ifelse(snp_pos$PhastCons_Score >= 0.3, "Moderately Conserved (Mammals)", "Poorly Conserved (Human-specific)"))); KoGES_Af6mapx <- snp_pos
#snp_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"), filters = "snp_filter", values = KoGES_Exomapx$rs_number, mart = snp_mart); snp_pos <- snp_pos[snp_pos$chr_name %in% as.character(1:22), ]; gr <- GRanges(seqnames = paste0("chr", snp_pos$chr_name), ranges = IRanges(start = snp_pos$chrom_start, width = 1)); scores <- gscores(phastCons100way.UCSC.hg38, gr); snp_pos$PhastCons_Score <- scores$default; snp_pos$Level <- ifelse(is.na(snp_pos$PhastCons_Score), "Unknown", ifelse(snp_pos$PhastCons_Score >= 0.9, "Highly Conserved (Vertebrates)", ifelse(snp_pos$PhastCons_Score >= 0.3, "Moderately Conserved (Mammals)", "Poorly Conserved (Human-specific)"))); KoGES_Exomapx <- snp_pos
#snp_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"), filters = "snp_filter", values = KoGES_Kchmapx$rs_number, mart = snp_mart); snp_pos <- snp_pos[snp_pos$chr_name %in% as.character(1:22), ]; gr <- GRanges(seqnames = paste0("chr", snp_pos$chr_name), ranges = IRanges(start = snp_pos$chrom_start, width = 1)); scores <- gscores(phastCons100way.UCSC.hg38, gr); snp_pos$PhastCons_Score <- scores$default; snp_pos$Level <- ifelse(is.na(snp_pos$PhastCons_Score), "Unknown", ifelse(snp_pos$PhastCons_Score >= 0.9, "Highly Conserved (Vertebrates)", ifelse(snp_pos$PhastCons_Score >= 0.3, "Moderately Conserved (Mammals)", "Poorly Conserved (Human-specific)"))); KoGES_Kchmapx <- snp_pos
#write.table(KoGES_Af5mapx, file="Output_COM_rsNo/KoGES_Af5mapx.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(KoGES_Af6mapx, file="Output_COM_rsNo/KoGES_Af6mapx.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(KoGES_Exomapx, file="Output_COM_rsNo/KoGES_Exomapx.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(KoGES_Kchmapx, file="Output_COM_rsNo/KoGES_Kchmapx.txt",sep="\t",quote=FALSE,row.names=FALSE)
KoGES_Af5mapx <- as.data.frame(fread("Output_COM_rsNo/KoGES_Af5mapx.txt.gz", header = TRUE, sep = "\t"))
KoGES_Af6mapx <- as.data.frame(fread("Output_COM_rsNo/KoGES_Af6mapx.txt.gz", header = TRUE, sep = "\t"))
KoGES_Exomapx <- as.data.frame(fread("Output_COM_rsNo/KoGES_Exomapx.txt.gz", header = TRUE, sep = "\t"))
KoGES_Kchmapx <- as.data.frame(fread("Output_COM_rsNo/KoGES_Kchmapx.txt.gz", header = TRUE, sep = "\t"))
head(KoGES_Af5mapx); dim(KoGES_Af5mapx) #343315,5
head(KoGES_Af6mapx); dim(KoGES_Af6mapx) #597551,5
head(KoGES_Exomapx); dim(KoGES_Exomapx) # 40560,5
head(KoGES_Kchmapx); dim(KoGES_Kchmapx) #440670,5


source("G:/내 드라이브/Rrunning/KoGES/Loading_KoGES_Epid.R") #1.5m
ASAS01to10[1:5,1:5]; dim(ASAS01to10)
CITY01to02[1:5,1:5]; dim(CITY01to02)
RURL01to05[1:5,1:5]; dim(RURL01to05)

KoGES_Af5gen <- KoGES_Af5$genotypes; KoGES_Af5mapy <- merge(KoGES_Af5map, KoGES_Af5mapx, by.x = "rs_number", by.y = "refsnp_id"); WEIGHT_PC_Af5 <- KoGES_Af5mapy$PhastCons_Score; all(KoGES_Af5mapy$marker.ID == colnames(KoGES_Af5gen)); PRS_Af5 <- numeric(nrow(KoGES_Af5gen))
KoGES_Af6gen <- KoGES_Af6$genotypes; KoGES_Af6mapy <- merge(KoGES_Af6map, KoGES_Af6mapx, by.x = "rs_number", by.y = "refsnp_id"); WEIGHT_PC_Af6 <- KoGES_Af6mapy$PhastCons_Score; all(KoGES_Af6mapy$marker.ID == colnames(KoGES_Af6gen)); PRS_Af6 <- numeric(nrow(KoGES_Af6gen))
KoGES_Exogen <- KoGES_Exo$genotypes; KoGES_Exomapy <- merge(KoGES_Exomap, KoGES_Exomapx, by.x = "rs_number", by.y = "refsnp_id"); WEIGHT_PC_Exo <- KoGES_Exomapy$PhastCons_Score; all(KoGES_Exomapy$marker.ID == colnames(KoGES_Exogen)); PRS_Exo <- numeric(nrow(KoGES_Exogen))
KoGES_Kchgen <- KoGES_Kch$genotypes; KoGES_Kchmapy <- merge(KoGES_Kchmap, KoGES_Kchmapx, by.x = "rs_number", by.y = "refsnp_id"); WEIGHT_PC_Kch <- KoGES_Kchmapy$PhastCons_Score; all(KoGES_Kchmapy$marker.ID == colnames(KoGES_Kchgen)); PRS_Kch <- numeric(nrow(KoGES_Kchgen))
#for (i in 1:nrow(KoGES_Af5gen)) {row_genotypes <- KoGES_Af5gen[i, ];   PRS_Af5[i] <- sum(row_genotypes * WEIGHT_PC_Af5, na.rm = TRUE)}
#for (i in 1:nrow(KoGES_Af6gen)) {row_genotypes <- KoGES_Af6gen[i, ];   PRS_Af6[i] <- sum(row_genotypes * WEIGHT_PC_Af6, na.rm = TRUE)}
#for (i in 1:nrow(KoGES_Exogen)) {row_genotypes <- KoGES_Exogen[i, ];   PRS_Exo[i] <- sum(row_genotypes * WEIGHT_PC_Exo, na.rm = TRUE)}
#for (i in 1:nrow(KoGES_Kchgen)) {row_genotypes <- KoGES_Kchgen[i, ];   PRS_Kch[i] <- sum(row_genotypes * WEIGHT_PC_Kch, na.rm = TRUE)}
#PRSdf_Af5 <- data.frame(DIST_ID = KoGES_Af5$fam$sample.ID, PRS = PRS_Af5)
#PRSdf_Af6 <- data.frame(DIST_ID = KoGES_Af6$fam$sample.ID, PRS = PRS_Af6)
#PRSdf_Exo <- data.frame(DIST_ID = KoGES_Exo$fam$sample.ID, PRS = PRS_Exo)
#PRSdf_Kch <- data.frame(DIST_ID = KoGES_Kch$fam$sample.ID, PRS = PRS_Kch)
#write.table(PRSdf_Af5, file="Output_COM_rsNo/PRSdf_Af5.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(PRSdf_Af6, file="Output_COM_rsNo/PRSdf_Af6.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(PRSdf_Exo, file="Output_COM_rsNo/PRSdf_Exo.txt",sep="\t",quote=FALSE,row.names=FALSE)
#write.table(PRSdf_Kch, file="Output_COM_rsNo/PRSdf_Kch.txt",sep="\t",quote=FALSE,row.names=FALSE)
PRSdf_Af5 <- as.data.frame(fread("Output_COM_rsNo/PRSdf_Af5.txt.gz", header = TRUE, sep = "\t"))
PRSdf_Af6 <- as.data.frame(fread("Output_COM_rsNo/PRSdf_Af6.txt.gz", header = TRUE, sep = "\t"))
PRSdf_Exo <- as.data.frame(fread("Output_COM_rsNo/PRSdf_Exo.txt.gz", header = TRUE, sep = "\t"))
PRSdf_Kch <- as.data.frame(fread("Output_COM_rsNo/PRSdf_Kch.txt.gz", header = TRUE, sep = "\t"))
head(PRSdf_Af5); dim(PRSdf_Af5) # 8840,2
head(PRSdf_Af6); dim(PRSdf_Af6) # 3693,2
head(PRSdf_Exo); dim(PRSdf_Exo) #14025,2
head(PRSdf_Kch); dim(PRSdf_Kch) #72291,2

ASAS01to10_BMI <- dplyr::select(ASAS01to10, DIST_ID, contains("BS_BMI")); head(ASAS01to10_BMI); dim(ASAS01to10_BMI)
ASAS01to10_WHR <- dplyr::select(ASAS01to10, contains("BS_WIT")) / dplyr::select(ASAS01to10, contains("BS_HIP")); colnames(ASAS01to10_WHR) <- gsub("_WIT","_WHR",colnames(ASAS01to10_WHR))
ASAS01to10_BMI <- cbind(ASAS01to10_BMI,ASAS01to10_WHR)
head(ASAS01to10_BMI); dim(ASAS01to10_BMI); length(unique(substr(colnames(ASAS01to10_BMI)[-1],6,11)))

ASAS01to10_LAB <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_LB_"));   ASAS01to10_LAB[1:5,1:5]; dim(ASAS01to10_LAB); length(unique(substr(colnames(ASAS01to10_LAB)[-1],6,11)))
ASAS01to10_MMI <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_M_MMI")); ASAS01to10_MMI[1:5,1:5]; dim(ASAS01to10_MMI); length(unique(substr(colnames(ASAS01to10_MMI)[-1],6,11)))
ASAS01to10_IBW <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_IB_"));   ASAS01to10_IBW[1:5,1:5]; dim(ASAS01to10_IBW); length(unique(substr(colnames(ASAS01to10_IBW)[-1],6,11)))
ASAS01to10_FAT <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_FT_"));   ASAS01to10_FAT[1:5,1:5]; dim(ASAS01to10_FAT); length(unique(substr(colnames(ASAS01to10_FAT)[-1],6,11)))
ASAS01to10_BDS <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_BDS_"));  ASAS01to10_BDS[1:5,1:5]; dim(ASAS01to10_BDS); length(unique(substr(colnames(ASAS01to10_BDS)[-1],6,11)))
ASAS01to10_SPR <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_SPR_"));  ASAS01to10_SPR[1:5,1:5]; dim(ASAS01to10_SPR); length(unique(substr(colnames(ASAS01to10_SPR)[-1],6,11)))
ASAS01to10_ALL <- cbind(ASAS01to10_LAB,ASAS01to10_MMI[,-1],ASAS01to10_IBW[,-1],ASAS01to10_FAT[,-1],ASAS01to10_BDS[,-1],ASAS01to10_SPR[,-1]); ASAS01to10_ALL[1:5,1:5]; dim(ASAS01to10_ALL); length(unique(substr(colnames(ASAS01to10_ALL)[-1],6,11)))

VARS <- unique(substr(colnames(ASAS01to10_ALL)[-1], 6, 11)); DAT <- merge(ASAS01to10_ALL, PRSdf_Af5, by="DIST_ID"); ALL_COLS <- as.vector(outer(paste0("NA", sprintf("%02d", 1:10), "_"), VARS, paste0)); MISS_COLS <- ALL_COLS[!(ALL_COLS %in% colnames(DAT))]; if(length(MISS_COLS)>0) DAT[, MISS_COLS] <- NA
SLP_MAT <- matrix(NA, nrow=nrow(DAT), ncol=length(VARS)); colnames(SLP_MAT) <- VARS

for(V in 1:length(VARS)) { COLS <- paste0("NA", sprintf("%02d", 1:10), "_", VARS[V]); VALID <- rowSums(!is.na(DAT[, COLS])) >= 6; for(I in 1:nrow(DAT)) { if(!VALID[I]) next; VALS <- as.numeric(DAT[I, COLS]); TMP <- numeric(0); for(J in 2:10) { for(K in 1:(J-1)) { if(!is.na(VALS[J]) && !is.na(VALS[K])) TMP <- c(TMP, (VALS[J]-VALS[K])/(J-K)) } }; if(length(TMP)>0) SLP_MAT[I,V] <- median(TMP) } }
DAT <- cbind(DAT, as.data.frame(SLP_MAT)); QNT <- quantile(DAT$PRS, probs=c(0.33, 0.66), na.rm=TRUE); DAT$GRP <- ifelse(DAT$PRS<=QNT[1], "LOW", ifelse(DAT$PRS<=QNT[2], "MID", "HIGH"))
DAT[1:3, 1:19]; dim(DAT)

BETA_V <- numeric(length(VARS)); PVAL_V <- numeric(length(VARS)); RSQ_V <- numeric(length(VARS))
pdf("Output_Evol/01_ALL_EVOLUTIONARY_ANALYSIS.pdf", width=10, height=8); par(mfrow=c(2, 2))
for(V in 1:length(VARS)) { if(sum(!is.na(DAT[, VARS[V]])) < 10) next; COLS <- paste0("NA", sprintf("%02d", 1:10), "_", VARS[V]); MDL <- lm(DAT[, VARS[V]] ~ DAT$PRS); SMR <- summary(MDL); if(nrow(SMR$coefficients)>1) { BETA_V[V] <- SMR$coefficients[2,1]; PVAL_V[V] <- SMR$coefficients[2,4]; RSQ_V[V] <- SMR$r.squared }; plot(DAT$PRS, DAT[, VARS[V]], main=paste("PRS vs", VARS[V], "Slope"), xlab="PRS Score", ylab="Sen's Slope", col="gray"); abline(MDL, col="red", lwd=2); PLOTa <- colMeans(DAT[DAT$GRP=="LOW", COLS], na.rm=TRUE); PLOTb <- colMeans(DAT[DAT$GRP=="MID", COLS], na.rm=TRUE); PLOTc <- colMeans(DAT[DAT$GRP=="HIGH", COLS], na.rm=TRUE); Y_MIN <- suppressWarnings(min(c(PLOTa, PLOTb, PLOTc), na.rm=TRUE)); Y_MAX <- suppressWarnings(max(c(PLOTa, PLOTb, PLOTc), na.rm=TRUE)); if(is.infinite(Y_MIN)) { Y_MIN <- 0; Y_MAX <- 1 }; plot(1:10, PLOTa, type="b", col="blue", ylim=c(Y_MIN, Y_MAX), xlab="Time Point", ylab=paste("Mean", VARS[V]), main=paste(VARS[V], "Change by PRS")); lines(1:10, PLOTb, type="b", col="green"); lines(1:10, PLOTc, type="b", col="red"); legend("topright", legend=c("Low PRS", "Mid", "High PRS"), col=c("blue", "green", "red"), lty=1); boxplot(DAT[, VARS[V]] ~ DAT$GRP, main=paste("Slope Dist by PRS Group (", VARS[V], ")"), xlab="PRS Group", ylab="Slope", col=c("lightblue", "lightgreen", "pink")); plot(DAT$PRS, jitter(as.numeric(DAT[, VARS[V]])), main=paste("PRS vs", VARS[V], "Jitter"), xlab="PRS", ylab="Jittered Slope", pch=16, col=rgb(0,0,0,0.2)) }; dev.off()

RESULT_TABLE <- data.frame(MODEL=VARS, BETA=BETA_V, PVAL=PVAL_V, R_SQ=RSQ_V)
write.table(RESULT_TABLE, "Output_Evol/01_ALL_PRS_ANALYSIS_RESULT.txt", sep="\t", row.names=FALSE, quote=FALSE)


TARGETS <- c("LB_GLU", "LB_G0h", "LB_CR1", "LB_CR4", "LB_EPC", "LB_NTR", "LB_TGY"); SUB_RES <- RESULT_TABLE[RESULT_TABLE$MODEL %in% TARGETS, ]; write.table(SUB_RES, "Output_Evol/05_STORY_TOP5_STATS.txt", sep="\t", row.names=FALSE, quote=FALSE)
library(ggplot2); pdf("Output_Evol/05_STORY_TOP5_PLOTS.pdf", width=10, height=8)
for(TGT in TARGETS) { if(TGT %in% colnames(DAT)) { TMP_DF <- data.frame(PRS=DAT$PRS, VAL=DAT[, TGT], GRP=DAT$GRP); PL <- ggplot(TMP_DF[!is.na(TMP_DF$VAL), ], aes(x=PRS, y=VAL)) + geom_point(aes(color=GRP), alpha=0.5, size=2, position=position_jitter(width=0.015)) + geom_smooth(method="lm", color="darkred", fill="lightcoral", se=TRUE, linetype="dashed") + labs(title=paste("Evolutionary Mismatch:", TGT, "10-Year Progression"), subtitle=sprintf("Effect Size (BETA): %.3e | P-Value: %.4f", SUB_RES$BETA[SUB_RES$MODEL==TGT], SUB_RES$PVAL[SUB_RES$MODEL==TGT]), x="Evolutionary Polygenic Risk Score (PRS)", y="10-Year Change Rate (Sen's Slope)") + scale_color_manual(values=c("LOW"="steelblue", "MID"="gray50", "HIGH"="firebrick")) + theme_minimal(base_size=14) + theme(plot.title=element_text(face="bold", size=16), plot.subtitle=element_text(color="midnightblue", face="italic"), legend.position="bottom", panel.grid.minor=element_blank()); print(PL) } }
dev.off()


# 1. 표현형(LB_GLU 기울기) 산출 및 유전형 데이터와 1:1 매칭
FAM <- KoGES_Af5$fam$sample.ID; PHE <- rep(NA, length(FAM)); COL <- paste0("NA", sprintf("%02d", 1:10), "_LB_GLU"); for(IDX in 1:length(FAM)) { MAT <- which(ASAS01to10_LAB$DIST_ID == FAM[IDX]); if(length(MAT)==0) next; VAL <- as.numeric(ASAS01to10_LAB[MAT, COL]); if(sum(!is.na(VAL))<6) next; TMP <- numeric(0); for(J in 2:10) { for(K in 1:(J-1)) { if(!is.na(VAL[J]) && !is.na(VAL[K])) TMP <- c(TMP, (VAL[J]-VAL[K])/(J-K)) } }; if(length(TMP)>0) PHE[IDX] <- median(TMP) }

# 2. 결측치 대치(Imputation) 및 GWAS 수행
GEN <- KoGES_Af5$genotypes; VLD <- which(!is.na(PHE)); GEN_IMP <- snp_fastImputeSimple(GEN, method="mean0"); GWS <- big_univLinReg(GEN_IMP, y.train=PHE[VLD], ind.train=VLD)

# 3. P-value 산출, 상위 20개 SNP 추출, SYMBOL 및 진화적 보존도(Level) 병합
PVL <- predict(GWS, log10=FALSE); TOP <- order(PVL)[1:20]; MAP <- KoGES_Af5$map[TOP, ]; MAP$PVAL <- PVL[TOP]; OUT <- merge(MAP, KoGES_Af5map, by="marker.ID", all.x=TRUE); OUT <- merge(OUT, KoGES_Af5mapx, by.x="rs_number", by.y="refsnp_id", all.x=TRUE); OUT <- OUT[order(OUT$PVAL), c("marker.ID", "chromosome.x", "physical.pos.x", "rs_number", "SYMBOL", "PVAL", "PhastCons_Score", "Level")]

# 4. 최종 결과 출력 및 TXT 파일 저장
print(head(OUT, 20)); write.table(OUT, "Output_Evol/01_TOP_SNPS_LB_GLU_EVOLUTION.txt", sep="\t", row.names=FALSE, quote=FALSE)



# 1. 진화 수준별 가중치 분리 및 2개의 독립적 PRS 산출 (매칭 오류 방지를 위해 열 이름 재정렬)
# 1. 진화 수준별 가중치 분리 (순서 꼬임 및 FBM 열 이름 의존성 완벽 차단)
MAP_IDX <- match(KoGES_Af5map$rs_number, KoGES_Af5mapx$refsnp_id); WGT_ALL <- KoGES_Af5mapx$PhastCons_Score[MAP_IDX]
WGT_CON <- ifelse(!is.na(WGT_ALL) & WGT_ALL >= 0.8, WGT_ALL, 0); WGT_HUM <- ifelse(!is.na(WGT_ALL) & WGT_ALL <= 0.2, 1 - WGT_ALL, 0)
PRS_CON <- numeric(nrow(KoGES_Af5gen)); PRS_HUM <- numeric(nrow(KoGES_Af5gen))
for(IDX in 1:nrow(KoGES_Af5gen)) { ROW_GEN <- KoGES_Af5gen[IDX, ]; PRS_CON[IDX] <- sum(ROW_GEN * WGT_CON, na.rm=TRUE); PRS_HUM[IDX] <- sum(ROW_GEN * WGT_HUM, na.rm=TRUE) }
PRS_DF_DUAL <- data.frame(DIST_ID=KoGES_Af5$fam$sample.ID, PRS_CON=PRS_CON, PRS_HUM=PRS_HUM)

# (이후 기존 코드의 2단계 DAT_MRG 생성 부분부터 그대로 실행)
# 2. 역학 데이터 병합 및 43개 LAB 지표의 Sen's Slope 일괄 산출 (결측 열 자동 생성 포함)
VAR_NMS <- unique(substr(colnames(ASAS01to10_LAB)[-1], 6, 11)); DAT_MRG <- merge(ASAS01to10_LAB, PRS_DF_DUAL, by="DIST_ID")
ALL_COL <- as.vector(outer(paste0("NA", sprintf("%02d", 1:10), "_"), VAR_NMS, paste0)); MIS_COL <- ALL_COL[!(ALL_COL %in% colnames(DAT_MRG))]; if(length(MIS_COL)>0) DAT_MRG[, MIS_COL] <- NA
SLP_MAT <- matrix(NA, nrow=nrow(DAT_MRG), ncol=length(VAR_NMS)); colnames(SLP_MAT) <- VAR_NMS
for(VDX in 1:length(VAR_NMS)) { COL_NMS <- paste0("NA", sprintf("%02d", 1:10), "_", VAR_NMS[VDX]); VAL_IDX <- rowSums(!is.na(DAT_MRG[, COL_NMS])) >= 6; for(IDX in 1:nrow(DAT_MRG)) { if(!VAL_IDX[IDX]) next; VAL_ARR <- as.numeric(DAT_MRG[IDX, COL_NMS]); TMP_ARR <- numeric(0); for(JDX in 2:10) { for(KDX in 1:(JDX-1)) { if(!is.na(VAL_ARR[JDX]) && !is.na(VAL_ARR[KDX])) TMP_ARR <- c(TMP_ARR, (VAL_ARR[JDX]-VAL_ARR[KDX])/(JDX-KDX)) } }; if(length(TMP_ARR)>0) SLP_MAT[IDX, VDX] <- median(TMP_ARR) } }
DAT_MRG <- cbind(DAT_MRG, as.data.frame(SLP_MAT))

# 3. 경쟁적 다중 회귀 분석(Competitive Regression) 수행 및 결과 정리
BET_CON <- numeric(length(VAR_NMS)); BET_HUM <- numeric(length(VAR_NMS)); PVS_CON <- numeric(length(VAR_NMS)); PVS_HUM <- numeric(length(VAR_NMS))
for(VDX in 1:length(VAR_NMS)) { if(sum(!is.na(DAT_MRG[, VAR_NMS[VDX]])) < 10) next; MDL_FIT <- lm(DAT_MRG[, VAR_NMS[VDX]] ~ DAT_MRG$PRS_CON + DAT_MRG$PRS_HUM); SMR_FIT <- summary(MDL_FIT); if(nrow(SMR_FIT$coefficients)>2) { BET_CON[VDX] <- SMR_FIT$coefficients[2,1]; BET_HUM[VDX] <- SMR_FIT$coefficients[3,1]; PVS_CON[VDX] <- SMR_FIT$coefficients[2,4]; PVS_HUM[VDX] <- SMR_FIT$coefficients[3,4] } }
RES_OUT <- data.frame(LAB_VAR=VAR_NMS, BETA_CONSERVED=BET_CON, BETA_HUMAN=BET_HUM, PVAL_CONSERVED=PVS_CON, PVAL_HUMAN=PVS_HUM)
write.table(RES_OUT, "Output_Evol/02_DUAL_PRS_EVOLUTION_RESULT.txt", sep="\t", row.names=FALSE, quote=FALSE)

# 4. 진화적 맵핑 2D Scatter Plot 시각화 및 PDF 반출
pdf("Output_Evol/02_EVOLUTIONARY_MAPPING_43LAB.pdf", width=12, height=10)
SIG_IDX <- ifelse(RES_OUT$PVAL_CONSERVED < 0.05 | RES_OUT$PVAL_HUMAN < 0.05, "red", "gray70")
X_LIM <- c(min(RES_OUT$BETA_CONSERVED, na.rm=TRUE) * 1.5, max(RES_OUT$BETA_CONSERVED, na.rm=TRUE) * 1.5); Y_LIM <- c(min(RES_OUT$BETA_HUMAN, na.rm=TRUE) * 1.5, max(RES_OUT$BETA_HUMAN, na.rm=TRUE) * 1.5)
plot(RES_OUT$BETA_CONSERVED, RES_OUT$BETA_HUMAN, main="Evolutionary Landscape of 43 Clinical Markers", xlab="Effect of Conserved PRS (BETA)", ylab="Effect of Human-Specific PRS (BETA)", pch=16, col=SIG_IDX, cex=1.5, xlim=X_LIM, ylim=Y_LIM)
abline(h=0, v=0, lty=2, col="black", lwd=1.5)
text(RES_OUT$BETA_CONSERVED, RES_OUT$BETA_HUMAN, labels=RES_OUT$LAB_VAR, pos=3, cex=0.85, col=ifelse(SIG_IDX=="red", "darkred", "gray40"), font=2)
legend("topleft", legend=c("Significant (p < 0.05 in either model)", "Not Significant"), pch=16, col=c("red", "gray70"), pt.cex=1.5)
dev.off()



ASAS01to10_cXr <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_CH_"));   ASAS01to10_cXr[1:5,1:5]; dim(ASAS01to10_cXr); length(unique(substr(colnames(ASAS01to10_cXr)[-1],6,11)))
ASAS01to10_DRG <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_DRG_"));  ASAS01to10_DRG[1:5,1:5]; dim(ASAS01to10_DRG); length(unique(substr(colnames(ASAS01to10_DRG)[-1],6,11)))
ASAS01to10_ChD <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("_ChD_"));  ASAS01to10_ChD[1:5,1:5]; dim(ASAS01to10_ChD); length(unique(substr(colnames(ASAS01to10_ChD)[-1],6,11)))
ASAS01to10_EKG <- dplyr::select(ASAS01to10, DIST_ID, starts_with("NA") & contains("EKG"));    ASAS01to10_EKG[1:5,1:5]; dim(ASAS01to10_EKG); length(unique(substr(colnames(ASAS01to10_EKG)[-1],6,11)))

ASAS01to10_cXr[1:5,1:5]; dim(ASAS01to10_cXr); table(ASAS01to10_cXr[,2]); unique(substr(colnames(ASAS01to10_cXr)[-1],6,11))
ASAS01to10_DRG[1:5,1:5]; dim(ASAS01to10_DRG); table(ASAS01to10_DRG[,2]); unique(substr(colnames(ASAS01to10_DRG)[-1],6,11))
ASAS01to10_ChD[1:5,1:5]; dim(ASAS01to10_ChD); table(ASAS01to10_ChD[,2]); unique(substr(colnames(ASAS01to10_ChD)[-1],6,11))
ASAS01to10_EKG[1:5,1:5]; dim(ASAS01to10_EKG); table(ASAS01to10_EKG[,2]); unique(substr(colnames(ASAS01to10_EKG)[-1],6,11))



# 1. 4개 범주형 데이터셋 리스트화 및 코딩 통일용 변환 벡터(0/1 차감 기준) 설정
LST_DAT <- list(ASAS01to10_cXr, ASAS01to10_DRG, ASAS01to10_ChD, ASAS01to10_EKG); LST_TYP <- c(0, 1, 1, 0); LST_NMS <- c("CXR", "DRG", "CHD", "EKG"); RES_ALL <- data.frame(CATEGORY=character(), VAR_NAME=character(), BETA_CON=numeric(), BETA_HUM=numeric(), PVAL_CON=numeric(), PVAL_HUM=numeric())

# 2. 4개 데이터셋을 순회하며 누적 발현 부담률(Burden Rate) 산출 및 다중 회귀분석 일괄 수행 (결측 변수 자동 스킵 방어 코드 적용)
Sys.time() #"2026-09-15 00:38:10 KST"
for(DDX in 1:4) { CUR_DAT <- LST_DAT[[DDX]]; CUR_VAR <- unique(substr(colnames(CUR_DAT)[-1], 6, 11)); MRG_DAT <- merge(CUR_DAT, PRS_DF_DUAL, by="DIST_ID"); ALL_COL <- as.vector(outer(paste0("NA", sprintf("%02d", 1:10), "_"), CUR_VAR, paste0)); MIS_COL <- ALL_COL[!(ALL_COL %in% colnames(MRG_DAT))]; if(length(MIS_COL)>0) MRG_DAT[, MIS_COL] <- NA; BRD_MAT <- matrix(NA, nrow=nrow(MRG_DAT), ncol=length(CUR_VAR)); colnames(BRD_MAT) <- CUR_VAR; for(VDX in 1:length(CUR_VAR)) { COL_NMS <- paste0("NA", sprintf("%02d", 1:10), "_", CUR_VAR[VDX]); VAL_IDX <- rowSums(!is.na(MRG_DAT[, COL_NMS])) >= 6; for(IDX in 1:nrow(MRG_DAT)) { if(!VAL_IDX[IDX]) next; VAL_ARR <- as.numeric(MRG_DAT[IDX, COL_NMS]); VAL_ARR <- VAL_ARR - LST_TYP[DDX]; BRD_MAT[IDX, VDX] <- mean(VAL_ARR, na.rm=TRUE) } }; MRG_DAT <- cbind(MRG_DAT, as.data.frame(BRD_MAT)); for(VDX in 1:length(CUR_VAR)) { if(sum(!is.na(MRG_DAT[, CUR_VAR[VDX]])) < 10) next; MDL_FIT <- lm(MRG_DAT[, CUR_VAR[VDX]] ~ MRG_DAT$PRS_CON + MRG_DAT$PRS_HUM); SMR_FIT <- summary(MDL_FIT); if(nrow(SMR_FIT$coefficients)>2) { RES_ALL <- rbind(RES_ALL, data.frame(CATEGORY=LST_NMS[DDX], VAR_NAME=CUR_VAR[VDX], BETA_CON=SMR_FIT$coefficients[2,1], BETA_HUM=SMR_FIT$coefficients[3,1], PVAL_CON=SMR_FIT$coefficients[2,4], PVAL_HUM=SMR_FIT$coefficients[3,4])) } } }
Sys.time() #"2026-09-15 01:06:19 KST"

# 3. 전체 범주형 지표 분석 결과 텍스트 반출 및 진화적 Landscape 시각화
write.table(RES_ALL, "Output_Evol/03_CATEGORICAL_DUAL_PRS_RESULT.txt", sep="\t", row.names=FALSE, quote=FALSE)
pdf("Output_Evol/03_CATEGORICAL_EVOLUTIONARY_MAPPING.pdf", width=14, height=12)
SIG_IDX <- ifelse(RES_ALL$PVAL_CON < 0.05 | RES_ALL$PVAL_HUM < 0.05, "red", "gray80")
X_LIM <- c(min(RES_ALL$BETA_CON, na.rm=TRUE)*1.5, max(RES_ALL$BETA_CON, na.rm=TRUE)*1.5); Y_LIM <- c(min(RES_ALL$BETA_HUM, na.rm=TRUE)*1.5, max(RES_ALL$BETA_HUM, na.rm=TRUE)*1.5)
plot(RES_ALL$BETA_CON, RES_ALL$BETA_HUM, main="Evolutionary Landscape of Disease/Drug Phenotypes", xlab="Effect of Conserved PRS (BETA)", ylab="Effect of Human-Specific PRS (BETA)", pch=16, col=SIG_IDX, cex=1.2, xlim=X_LIM, ylim=Y_LIM)
abline(h=0, v=0, lty=2, col="black", lwd=1.5)
SIG_DAT <- RES_ALL[SIG_IDX=="red", ]; if(nrow(SIG_DAT)>0) text(SIG_DAT$BETA_CON, SIG_DAT$BETA_HUM, labels=SIG_DAT$VAR_NAME, pos=3, cex=0.8, col="darkred", font=2)
legend("topleft", legend=c("Significant (p < 0.05)", "Not Significant"), pch=16, col=c("red", "gray80"), pt.cex=1.2)
dev.off()

x02 <- as.data.frame(fread("Output_Evol/02.txt",sep="\t"))[,1:5]; dim(x02)
x03 <- as.data.frame(fread("Output_Evol/03.txt",sep="\t"))[,1:5]; dim(x03)
colnames(x02) <- colnames(x03)
x01 <- rbind(x02,x03)
head(x01); dim(x01)

# 1. 데이터 로드 (x01 데이터 프레임 사용)
DAT_RAW <- x01

# 2. P-value 기반 색상 그라데이션 매핑 (0.05 이상은 일괄 회색 gray80 처리)
PAL_GRN <- colorRampPalette(c("#A1D99B", "#00441B"))(100); PAL_YLW <- colorRampPalette(c("#FFF200", "#B36B00"))(100)
IDX_CON <- pmin(100, pmax(1, round((-log10(DAT_RAW$PVAL_CON) - 1.3) / 2 * 100))); COL_CON <- ifelse(DAT_RAW$PVAL_CON >= 0.05, "gray80", PAL_GRN[IDX_CON])
IDX_HUM <- pmin(100, pmax(1, round((-log10(DAT_RAW$PVAL_HUM) - 1.3) / 2 * 100))); COL_HUM <- ifelse(DAT_RAW$PVAL_HUM >= 0.05, "gray80", PAL_YLW[IDX_HUM])

# 3. 막대 그래프용 데이터 교차 배열(Interleaving) 및 위아래 출력 순서 반전(rev)
VEC_BETA <- rev(as.vector(rbind(DAT_RAW$BETA_CON, DAT_RAW$BETA_HUM))); VEC_COL <- rev(as.vector(rbind(COL_CON, COL_HUM)))
VEC_LBL <- rev(as.vector(rbind(paste0(DAT_RAW$VAR_NAME, " +CON"), paste0(DAT_RAW$VAR_NAME, " +HUM"))))

# 4. 극단값(0.0005 이상)을 제외한 Zoom-in X축 스케일 설정 및 범례 오른쪽 아래(bottomright) 이동
VAL_SUB <- VEC_BETA[abs(VEC_BETA) < 0.0005]; X_LIM <- c(min(VAL_SUB) * 1.2, max(VAL_SUB) * 1.4)
pdf("Output_Evol/04_EVOLUTIONARY_BETA_CON_HUM_BARPLOT_ZOOMED.pdf", width=10, height=7); par(mar=c(5, 8, 4, 2))
barplot(VEC_BETA, horiz=TRUE, names.arg=VEC_LBL, col=VEC_COL, las=1, cex.names=0.9, xlab="BETA Coefficient (Zoomed Scale)", main="Evolutionary Effect Sizes (Zoomed in on Micro-Effects)", border="gray50", xlim=X_LIM, xpd=FALSE)
abline(v=0, col="black", lwd=1.5)
legend("bottomright", legend=c("CON (p < 0.05)", "HUM (p < 0.05)", "Not Significant (p >= 0.05)"), fill=c(PAL_GRN[75], PAL_YLW[75], "gray80"), border="gray50", bg="white", inset=c(0.02, 0.02))
dev.off()






