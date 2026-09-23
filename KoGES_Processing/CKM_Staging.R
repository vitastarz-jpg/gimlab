# Module 0
library(data.table); library(openxlsx); library(PooledCohort)
ASAS_CKM <- fread("G:/내 드라이브/Rrunning/KoGES/Epid/ASAS01to10.txt.gz")
ASAS_CKM[1:5, 1:5]; dim(ASAS_CKM) #10030,17802
orig_cols <- colnames(ASAS_CKM)

# Module 1: UACR
ASAS_CKM$UACR_ALBU <- ASAS_CKM$AS01_ALBU_U/ASAS_CKM$AS01_CREATINE_U*1000
ASAS_CKM$UACR_MICALBU <- ASAS_CKM$AS01_MICALBU_U/ASAS_CKM$AS01_CREATINE_U*100
ASAS_CKM$UACR_combined <- ifelse(!is.na(ASAS_CKM$UACR_ALBU), ASAS_CKM$UACR_ALBU, ASAS_CKM$UACR_MICALBU)

# Module 2: eGFR
ASAS_CKM$kappa <- ifelse(ASAS_CKM$AS01_SEX==2, 0.7, 0.9)
ASAS_CKM$alpha <- ifelse(ASAS_CKM$AS01_SEX==2, -0.241, -0.302)
ASAS_CKM$sex_factor <- ifelse(ASAS_CKM$AS01_SEX==2, 1.012, 1)
ASAS_CKM$AS01_eGFR <- 142*pmin(ASAS_CKM$AS01_CREATININE_TR1/ASAS_CKM$kappa,1)^ASAS_CKM$alpha * pmax(ASAS_CKM$AS01_CREATININE_TR1/ASAS_CKM$kappa,1)^(-1.2) * 0.9938^ASAS_CKM$AS01_AGE * ASAS_CKM$sex_factor
ASAS_CKM$CKM_eligible <- !is.na(ASAS_CKM$AS01_eGFR) & !is.na(ASAS_CKM$UACR_combined);

# Module 3: GFR category + Albuminuria category
ASAS_CKM$GFR_category <- cut(ASAS_CKM$AS01_eGFR, breaks=c(-Inf,15,30,45,60,90,Inf), labels=c("G5 <15","G4 15-29","G3b 30-44","G3a 45-59","G2 60-89","G1 >=90"), right=FALSE)
ASAS_CKM$Albuminuria_category <- cut(ASAS_CKM$UACR_combined, breaks=c(-Inf,30,300,Inf), labels=c("A1 <30","A2 30-299","A3 >=300"), right=FALSE)

# Module 4: KDIGO risk
ASAS_CKM$KDIGO_risk <- NA_character_
K_av <- !is.na(ASAS_CKM$GFR_category) & !is.na(ASAS_CKM$Albuminuria_category)
ASAS_CKM$KDIGO_risk[K_av & ASAS_CKM$GFR_category %in% c("G1 >=90","G2 60-89") & ASAS_CKM$Albuminuria_category=="A1 <30"] <- "Low"
ASAS_CKM$KDIGO_risk[K_av & ((ASAS_CKM$GFR_category %in% c("G1 >=90","G2 60-89") & ASAS_CKM$Albuminuria_category=="A2 30-299") | (ASAS_CKM$GFR_category=="G3a 45-59" & ASAS_CKM$Albuminuria_category=="A1 <30"))] <- "Moderate"
ASAS_CKM$KDIGO_risk[K_av & ((ASAS_CKM$GFR_category %in% c("G1 >=90","G2 60-89") & ASAS_CKM$Albuminuria_category=="A3 >=300") | (ASAS_CKM$GFR_category=="G3a 45-59" & ASAS_CKM$Albuminuria_category=="A2 30-299") | (ASAS_CKM$GFR_category=="G3b 30-44" & ASAS_CKM$Albuminuria_category=="A1 <30"))] <- "High"
ASAS_CKM$KDIGO_risk[K_av & ((ASAS_CKM$GFR_category=="G3a 45-59" & ASAS_CKM$Albuminuria_category=="A3 >=300") | (ASAS_CKM$GFR_category=="G3b 30-44" & ASAS_CKM$Albuminuria_category %in% c("A2 30-299","A3 >=300")) | (ASAS_CKM$GFR_category %in% c("G4 15-29","G5 <15")))] <- "Very high"

# Module 5: Metabolic variables
ASAS_CKM$WAIST <- rowMeans(cbind(ASAS_CKM$AS01_WAIST1, ASAS_CKM$AS01_WAIST2, ASAS_CKM$AS01_WAIST3), na.rm=TRUE); ASAS_CKM$WAIST[is.nan(ASAS_CKM$WAIST)] <- NA
ASAS_CKM$SBP <- rowMeans(cbind(ASAS_CKM$AS01_BPSIT1LS, ASAS_CKM$AS01_BPSIT1RS), na.rm=TRUE); ASAS_CKM$SBP[is.nan(ASAS_CKM$SBP)] <- NA
ASAS_CKM$DBP <- rowMeans(cbind(ASAS_CKM$AS01_BPSIT1LD, ASAS_CKM$AS01_BPSIT1RD), na.rm=TRUE); ASAS_CKM$DBP[is.nan(ASAS_CKM$DBP)] <- NA
ASAS_CKM$FPG <- ASAS_CKM$AS01_GLU0_TR
ASAS_CKM$TG <- ASAS_CKM$AS01_TG_TR
ASAS_CKM$HDL <- ASAS_CKM$AS01_HDL_TR
ASAS_CKM$TC <- ASAS_CKM$AS01_TCHL_TR
ASAS_CKM$HbA1c <- ASAS_CKM$AS01_HBA1C

# Module 6-1: Stage 1
ASAS_CKM$CKM_adiposity <- (!is.na(ASAS_CKM$AS01_BMI) & ASAS_CKM$AS01_BMI>=23) | (ASAS_CKM$AS01_SEX==1 & !is.na(ASAS_CKM$WAIST) & ASAS_CKM$WAIST>=90) | (ASAS_CKM$AS01_SEX==2 & !is.na(ASAS_CKM$WAIST) & ASAS_CKM$WAIST>=80)
ASAS_CKM$CKM_prediabetes <- (!is.na(ASAS_CKM$FPG) & ASAS_CKM$FPG>=100 & ASAS_CKM$FPG<126) | (!is.na(ASAS_CKM$HbA1c) & ASAS_CKM$HbA1c>=5.7 & ASAS_CKM$HbA1c<6.5)
ASAS_CKM$CKM_stage1_condition <- ASAS_CKM$CKM_adiposity | ASAS_CKM$CKM_prediabetes

# Module 6-2: Stage 2
ASAS_CKM$CKM_HTN <- (!is.na(ASAS_CKM$SBP) & ASAS_CKM$SBP>=130) | (!is.na(ASAS_CKM$DBP) & ASAS_CKM$DBP>=80) | (!is.na(ASAS_CKM$AS01_DRUGHTCU) & ASAS_CKM$AS01_DRUGHTCU==2)
ASAS_CKM$CKM_high_TG <- !is.na(ASAS_CKM$TG) & ASAS_CKM$TG>=150
ASAS_CKM$CKM_DM <- (!is.na(ASAS_CKM$FPG) & ASAS_CKM$FPG>=126) | (!is.na(ASAS_CKM$HbA1c) & ASAS_CKM$HbA1c>=6.5) | (!is.na(ASAS_CKM$AS01_PDDM) & ASAS_CKM$AS01_PDDM==2) | (!is.na(ASAS_CKM$AS01_DRUGDMCU) & ASAS_CKM$AS01_DRUGDMCU==2) | (!is.na(ASAS_CKM$AS01_DRUGINSCU) & ASAS_CKM$AS01_DRUGINSCU==2)
ASAS_CKM$MetS_count <- rowSums(cbind((ASAS_CKM$AS01_SEX==1 & !is.na(ASAS_CKM$WAIST) & ASAS_CKM$WAIST>=90)|(ASAS_CKM$AS01_SEX==2 & !is.na(ASAS_CKM$WAIST) & ASAS_CKM$WAIST>=80), !is.na(ASAS_CKM$TG) & ASAS_CKM$TG>=150, (ASAS_CKM$AS01_SEX==1 & !is.na(ASAS_CKM$HDL) & ASAS_CKM$HDL<40)|(ASAS_CKM$AS01_SEX==2 & !is.na(ASAS_CKM$HDL) & ASAS_CKM$HDL<50), ASAS_CKM$CKM_HTN, !is.na(ASAS_CKM$FPG) & ASAS_CKM$FPG>=100), na.rm=TRUE)
ASAS_CKM$CKM_MetS <- ASAS_CKM$MetS_count>=3
ASAS_CKM$CKM_CKD_stage2 <- ASAS_CKM$KDIGO_risk %in% c("Moderate","High")
ASAS_CKM$CKM_stage2_condition <- ASAS_CKM$CKM_HTN | ASAS_CKM$CKM_high_TG | ASAS_CKM$CKM_MetS | ASAS_CKM$CKM_DM | ASAS_CKM$CKM_CKD_stage2
ASAS_CKM$CKM_clinical_CVD <- (!is.na(ASAS_CKM$AS01_PDMI) & ASAS_CKM$AS01_PDMI==2) | (!is.na(ASAS_CKM$AS01_PDCD) & ASAS_CKM$AS01_PDCD==2) | (!is.na(ASAS_CKM$AS01_PDCH) & ASAS_CKM$AS01_PDCH==2) | (!is.na(ASAS_CKM$AS01_PDCV) & ASAS_CKM$AS01_PDCV==2);
ASAS_CKM$PREVENT_sex <- ifelse(ASAS_CKM$AS01_SEX==1, "male", ifelse(ASAS_CKM$AS01_SEX==2, "female", NA))
ASAS_CKM$PREVENT_smoke <- ifelse(ASAS_CKM$AS01_SMOKEA %in% c(2,3), "yes", ifelse(ASAS_CKM$AS01_SMOKEA %in% c(0,1), "no", NA))
ASAS_CKM$PREVENT_bpmed <- NA_character_
ASAS_CKM$PREVENT_bpmed[ASAS_CKM$AS01_DRUGHTCU==2] <- "yes"
ASAS_CKM$PREVENT_bpmed[ASAS_CKM$AS01_DRUGHT==1 | ASAS_CKM$AS01_DRUGHTCU==1] <- "no"
ASAS_CKM$PREVENT_statin <- NA_character_
ASAS_CKM$PREVENT_statin[ASAS_CKM$AS01_DRUGLPCU==2] <- "yes"
ASAS_CKM$PREVENT_statin[ASAS_CKM$AS01_DRUGLP==1 | ASAS_CKM$AS01_DRUGLPCU==1] <- "no"
ASAS_CKM$PREVENT_diabetes <- NA_character_
ASAS_CKM$PREVENT_diabetes[ASAS_CKM$CKM_DM==TRUE] <- "yes"
ASAS_CKM$PREVENT_diabetes[is.na(ASAS_CKM$PREVENT_diabetes) & ASAS_CKM$AS01_PDDM==1 & ASAS_CKM$FPG<126 & ASAS_CKM$HbA1c<6.5] <- "no"
P_av <- ASAS_CKM$CKM_eligible & ASAS_CKM$AS01_AGE>=30 & ASAS_CKM$AS01_AGE<=79 & !ASAS_CKM$CKM_clinical_CVD & !is.na(ASAS_CKM$PREVENT_sex) & !is.na(ASAS_CKM$PREVENT_smoke) & !is.na(ASAS_CKM$TC) & !is.na(ASAS_CKM$HDL) & !is.na(ASAS_CKM$SBP) & !is.na(ASAS_CKM$PREVENT_bpmed) & !is.na(ASAS_CKM$PREVENT_statin) & !is.na(ASAS_CKM$PREVENT_diabetes) & !is.na(ASAS_CKM$AS01_BMI) & !is.na(ASAS_CKM$AS01_eGFR)
ASAS_CKM$PREVENT_CVD_10yr <- NA_real_
ASAS_CKM$PREVENT_CVD_10yr[P_av] <- predict_10yr_cvd_risk(age_years=ASAS_CKM$AS01_AGE[P_av], sex=ASAS_CKM$PREVENT_sex[P_av], smoke_current=ASAS_CKM$PREVENT_smoke[P_av], chol_total_mgdl=ASAS_CKM$TC[P_av], chol_hdl_mgdl=ASAS_CKM$HDL[P_av], bp_sys_mmhg=ASAS_CKM$SBP[P_av], bp_meds=ASAS_CKM$PREVENT_bpmed[P_av], statin_meds=ASAS_CKM$PREVENT_statin[P_av], diabetes=ASAS_CKM$PREVENT_diabetes[P_av], bmi=ASAS_CKM$AS01_BMI[P_av], egfr_mlminm2=ASAS_CKM$AS01_eGFR[P_av], equation_version="Khan_2023", prevent_type="base", override_boundary_errors=TRUE)

# Module 6-3: Stage 3
ASAS_CKM$CKM_PREVENT_stage3 <- !is.na(ASAS_CKM$PREVENT_CVD_10yr) & ASAS_CKM$PREVENT_CVD_10yr>=0.20;
ASAS_CKM$CKM_CKD_stage3 <- ASAS_CKM$KDIGO_risk=="Very high"
ASAS_CKM$CKM_subclinical_CVD <- FALSE
ASAS_CKM$CKM_stage3_condition <- ASAS_CKM$CKM_CKD_stage3 | (ASAS_CKM$CKM_stage2_condition & ASAS_CKM$CKM_PREVENT_stage3)

# Module 6-4: Stage 4
ASAS_CKM$CKM_any_risk <- ASAS_CKM$CKM_stage1_condition | ASAS_CKM$CKM_stage2_condition | ASAS_CKM$CKM_CKD_stage3
ASAS_CKM$CKM_stage4_condition <- ASAS_CKM$CKM_clinical_CVD & ASAS_CKM$CKM_any_risk

# Module 7: Final CKM stage
ASAS_CKM$CKM_stage_final <- NA_integer_
ASAS_CKM$CKM_stage_final[ASAS_CKM$CKM_eligible] <- 0
ASAS_CKM$CKM_stage_final[ASAS_CKM$CKM_eligible & ASAS_CKM$CKM_stage1_condition] <- 1
ASAS_CKM$CKM_stage_final[ASAS_CKM$CKM_eligible & ASAS_CKM$CKM_stage2_condition] <- 2
ASAS_CKM$CKM_stage_final[ASAS_CKM$CKM_eligible & !ASAS_CKM$CKM_clinical_CVD & ASAS_CKM$CKM_stage3_condition] <- 3
ASAS_CKM$CKM_stage_final[ASAS_CKM$CKM_eligible & ASAS_CKM$CKM_stage4_condition] <- 4
ASAS_CKM$CKM_stage4_subtype <- NA_character_
ASAS_CKM$CKM_stage4_subtype[ASAS_CKM$CKM_stage_final==4 & ASAS_CKM$AS01_eGFR>=15] <- "Stage 4a"
ASAS_CKM$CKM_stage4_subtype[ASAS_CKM$CKM_stage_final==4 & ASAS_CKM$AS01_eGFR<15] <- "Stage 4b"

# Module 8: Select new variables
ASAS_CKM_add <- as.data.frame(ASAS_CKM)[, c("DIST_ID",setdiff(names(ASAS_CKM), orig_cols))]
tb_GA <- as.data.frame.matrix(table(ASAS_CKM$GFR_category, ASAS_CKM$Albuminuria_category, useNA="no"))
tb_GA <- cbind(GFR_category=rownames(tb_GA), tb_GA)
rownames(tb_GA) <- NULL