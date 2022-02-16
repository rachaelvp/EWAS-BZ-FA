################################################################################
# In this Import script we consider the sample level filtering performed by the
# CGR. On top of that, we also consider further recommended filtering on both
# the probes and the samples.
################################################################################

# Set up workspace and read in IDAT files if necessary
library(here)
library(dplyr)
library(ChAMP)

# setting project and data directories
raw_data <- ("raw_data.RData")

# check existence of raw data file 
if (!file.exists(normalizePath(raw_data))) {

  print("raw_data object not found; attempting to read raw files...")

  load("rgsetAll.Rdata")

  RGset_raw <- rgsetAll
  RGset_raw <- BiocGenerics::updateObject(RGset_raw)
  
  # ADD extra step to make the object valid 
  e <- new.env()
  load("rgsetAll.Rdata", envir = e)
  valid_RGset <- updateObject(rgsetAll)
  # check validity 
  validObject(valid_RGset)
  RGset_raw <- valid_RGset
  # check dimension 
  dim(RGset_raw)
  
  Mset_raw <- minfi::preprocessRaw(RGset_raw)
  betas_raw <- minfi::getBeta(Mset_raw)
# 1 - 4 are the NA07057_1 to NA07057_14 Internal_Controls

  ## load phenotype date 
  ## START HERE
  pheno <- read.csv("pheno_raw.csv")

  ##### only include the samples for which we have phenotype information #######
  # keep NA rows 
  samples <- pheno %>%
  filter(Internal_CGR_control != 1 | is.na(Internal_CGR_control)) 
  # remove failed CGR_QC
  samples <- samples %>% filter(Failed_CGR_QC != 1 | is.na(Failed_CGR_QC))
  
  # remove replicate
  # row 50 and 51
  # (1) keep row 50 for now
  #samples <- samples[-51,]
  # (2) keep row 51 for now: replicate 2
  samples <- samples[-50,]
  # convert factor ID to string ID
  rownames(samples) <- as.character(samples$TargetID)
  # remove redundant ID column 
  #samples <- samples[,-1]

  # adding the X in front of the name
  names <- sub("^", "X", colnames(betas_raw))
  colnames(betas_raw) <- names
  colnames(betas_raw) %in% rownames(samples)

# remove FALSE columns 
# ====== CHANGE column index =======
  # r1
  #betas <- betas_raw[,-c(1,2,3,4,55,c(77:88))]
  # r2
  betas <- betas_raw[,-c(1,2,3,4,54,c(77:88))]
  dim(betas)
 # [1] 485512     71

  #### only include the phenotype info for samples that CGR didn't filter out ##
  all_samples <- as.character(samples$TargetID)
  unfiltered_samples <- colnames(betas)
  filtered <- all_samples[!all_samples %in% unfiltered_samples]
  filtered

  
  ## ===========Replace pheno with samples ================ ##
  pheno <- samples
  
  #which(rownames(samples) %in% filtered)
  #pheno <- samples[-c(26,43,74,75,79,80,97,101,105,111,112,114,115,127,128,129,
  #130,132,134,142,143,144,145),]


  # we need to make sure that colnames(betas) == rownames(pheno)
  # we have two places where they are switched!

  #betas <- betas[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
  #19,18,20:ncol(betas))]

  #betas <- betas[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
  #24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,43,42,44:ncol(betas))]

  # numeric variables
  #CD8_Central <- as.numeric(as.character(pheno$GZ17_CD45RAN_CCR7P_CD8))
  #CD8_Effector <- as.numeric(as.character(pheno$GZ21_CD45RAN_CCR7N_CD8))
  #CD8_Naive <- as.numeric(as.character(pheno$GZ19_CD45RAP_CCR7P_CD8))
  #CD4_Effector <- as.numeric(as.character(pheno$GZ15_CD45RAN_CCR7N_CD4))
  #CD4_Central <- as.numeric(as.character(pheno$GZ11_CD45RAN_CCR7P_CD4))
  #CD4_Naive <- as.numeric(as.character(pheno$GZ13_CD45RAP_CCR7P_CD4P))
  #CD4_FoxP3 <- as.numeric(as.character(pheno$GZ26_CD4P_FOXP3P))
  #CD25_FoxP3 <- as.numeric(as.character(pheno$CD25P_FOXP3P))
  #CD4_CD25 <- as.numeric(as.character(pheno$GZ25_CD4P_CD25P))
  #NKcell <- as.numeric(as.character(pheno$GZ08LY_CD3N_CD56P))
  #Bcell <- as.numeric(as.character(pheno$GZ07_CD19P))
  #CD8_Tcell <- as.numeric(as.character(pheno$GZ03LY_CD8P_CD3P))
  #CD4_Tcell <- as.numeric(as.character(pheno$GZ02LY_CD4P_CD3P))
  #Tcell <- as.numeric(as.character(pheno$GZ01LY_CD3P))
  Granulocyte <- as.numeric(as.character(pheno$GR_CBC))
  Monocyte <- as.numeric(as.character(pheno$MO))
  Lymphocyte <- as.numeric(as.character(pheno$LY))
  #Platelet <- as.numeric(as.character(pheno$PLT))
  WBC <- as.numeric(as.character(pheno$WBC))
  #pack_yrs <- as.numeric(as.character(pheno$pack_yrs))
  #smk_duration <- as.numeric(as.character(pheno$smk_duration))
  #smk_intensity <- as.numeric(as.character(pheno$smk_intensity))
  #TCA <- as.numeric(as.character(pheno$TCA))
  
  fa_ppm <- as.numeric(as.character(pheno$fa_ppm_))
  BMI <- as.numeric(as.character(pheno$BMI))
  age <- as.numeric(as.character(pheno$age))
  
  # factor variables
  #TCE_M_C <- as.factor(as.character(pheno$TCE_M_C))
  #dich_tce <- as.factor(as.character(pheno$dich_tce))
  fa <- as.factor(as.character(pheno$fa))
  Subject_ID <- as.factor(as.character(pheno$Subject_ID))
  Sample_ID <- as.factor(as.character(pheno$Sample_ID))
  Sample_Plate <- as.factor(as.character(pheno$Sample_Plate))
  Target_ID <- as.factor(as.character(pheno$TargetID))
  Replicate <- as.factor(as.character(pheno$Replicate))
  sex <- as.factor(as.character(pheno$sex))
  smk <- as.factor(as.character(pheno$smk))
  alcol <- as.factor(as.character(pheno$alcol))
  infection <- as.factor(as.character(pheno$infection))

  
  pheno <- data.frame(Subject_ID, Sample_ID, Target_ID, Replicate, Sample_Plate,fa,
  fa_ppm, sex, age, BMI, smk, alcol, infection,
  WBC, Lymphocyte, Monocyte, Granulocyte)
  # rownames = targetID
  rownames(pheno) <- pheno[,3]

  ################################################################################
  # Examining the QC plots
  # Considering an additional layer of filtering
  ################################################################################

 #Mset <- Mset_raw[,-c(1,2,3,4,55,c(77:88))]
  Mset <- Mset_raw[,-c(1,2,3,4,54,c(77:88))]
  # we need to make sure that colnames(Mset) == rownames(pheno)
  # we have two places where they are switched!
  #Mset <- Mset[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
  #19,18,20:ncol(Mset))]
  #Mset <- Mset[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
  #24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,43,42,44:ncol(Mset))]
  names <- sub("^", "X", colnames(Mset))
  colnames(Mset) <- names
  colnames(Mset) == rownames(pheno)
  # TRUE

  #RGset <- RGset_raw[,-c(1,2,3,4,55,c(77:88))]
  RGset <- RGset_raw[,-c(1,2,3,4,54,c(77:88))]
  # we need to make sure that colnames(RGset) == rownames(pheno)
  # we have two places where they are switched!
  #RGset <- RGset[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
  #19,18,20:ncol(RGset))]
  #RGset <- RGset[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
  #24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,43,42,44:ncol(RGset))]
  names <- sub("^", "X", colnames(RGset))
  colnames(RGset) <- names
  colnames(RGset) == rownames(pheno)
  # TRUE

  champ.QC(beta = getBeta(Mset), pheno=pheno$fa,
  mdsPlot=TRUE, densityPlot=TRUE, dendrogram=TRUE, PDFplot=TRUE, Rplot= FALSE,
  Feature.sel="None", resultsDir=here::here("SuperFund/replicate2", "QC",
  "preFiltering"))

  ## ERROR msg: 
  # << PhenoTypes.lv generated successfully. >>
  # Error in while (step == 0) { : missing value where TRUE/FALSE needed
  
  ## Failed to generate any plots 
  
  ## =================== Debug ===================== #
  champ.SVD(beta = betas, rgSet=RGset, pd=pheno,
  RGEffect=FALSE, Rplot=FALSE, resultsDir=here::here("QC",
  "preFiltering", "SVD"))

  ######################### Internal Control Probes ###########################


  # minfi::qcReport(RGset, sampGroups = pheno$fa, sampNames = pheno$Subject_ID,
  # pdf = here::here("SuperFund/graphs","QC", "ControlProbes.pdf"))
  
  minfi::qcReport(RGset, sampGroups = pheno$fa, sampNames = pheno$Subject_ID,
                  pdf = here::here("replicate2","QC", "ControlProbes.pdf"))
  
  # pdf(file = here::here("SuperFund/graphs","QC",
  # "Internal_ControlProbes.pdf"))
  
  pdf(file = here::here("replicate2","QC",
  "Internal_ControlProbes.pdf"))
  
  
  ENmix::plotCtrl(RGset)
  dev.off()

  ############################ Probe/CpG Level QC ############################

  ##################### without the sample filtering before
  detP <- minfi::detectionP(RGset)

  filtered_P <- champ.filter(beta=getBeta(Mset), M=getM(Mset),
              pd = pheno, beadcount = NULL, detP = detP,
              arraytype = "450K")

  ######  replicate 1 ################
  # Filtering probes with a detection p-value above 0.01.
  # Removing 13988 probes.
  # If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
  # 
  # Filtering NoCG Start
  # Only Keep CpGs, removing 2686 probes from the analysis.
  # 
  # Filtering SNPs Start
  # Using general 450K SNP list for filtering.
  # Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
  # Removing 56931 probes from the analysis.
  # 
  # Filtering MultiHit Start
  # Filtering probes that align to multiple locations as identified in Nordlund et al
  # Removing 11 probes from the analysis.
  # 
  # Filtering XY Start
  # Filtering probes located on X,Y chromosome, removing 9429 probes from the analysis.
  # 
  # Updating PD file
  # 
  # Fixing Outliers Start
  # Replacing all value smaller/equal to 0 with smallest positive value.
  # Replacing all value greater/equal to 1 with largest value below 1..
  # [ Section 2: Filtering Done ]
  # 
  # All filterings are Done, now you have 402467 probes and 71 samples.

  
  
  
################ replicate 2 #################
  
  # Filtering probes with a detection p-value above 0.01.
  # Removing 14173 probes.
  # If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
  # 
  # Filtering NoCG Start
  # Only Keep CpGs, removing 2662 probes from the analysis.
  # 
  # Filtering SNPs Start
  # Using general 450K SNP list for filtering.
  # Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
  # Removing 56909 probes from the analysis.
  # 
  # Filtering MultiHit Start
  # Filtering probes that align to multiple locations as identified in Nordlund et al
  # Removing 11 probes from the analysis.
  # 
  # Filtering XY Start
  # Filtering probes located on X,Y chromosome, removing 9419 probes from the analysis.
  # 
  # Updating PD file
  # 
  # Fixing Outliers Start
  # Replacing all value smaller/equal to 0 with smallest positive value.
  # Replacing all value greater/equal to 1 with largest value below 1..
  # [ Section 2: Filtering Done ]
  # 
  # All filterings are Done, now you have 402338 probes and 71 samples.
  # 
  # [<<<<< ChAMP.FILTER END >>>>>>]
  # [===========================]
  # [You may want to process champ.QC() next.]
  
  pheno_P <- filtered_P$pd
  betas_P <- filtered_P$beta
  mvals_P <- filtered_P$M

  filtered_CpG <- rownames(betas_P)
  unfiltered_CpG <- rownames(betas)
  outCpG_P <- unfiltered_CpG[!(unfiltered_CpG %in% filtered_CpG)]

  champ.QC(beta = betas_P, pheno=pheno_P$fa, mdsPlot=TRUE,
  densityPlot=TRUE, dendrogram=TRUE, PDFplot=TRUE, Rplot= FALSE,
  Feature.sel="None", resultsDir=here::here("replicate2", "QC",
  "postProbeFiltering"))

  champ.SVD(beta = betas_P, rgSet=RGset, pd=pheno_P,
  RGEffect=TRUE, Rplot=FALSE, resultsDir=here::here("replicate2", "QC",
  "postProbeFiltering", "SVD"))

  ############################ Sample Level QC ################################

  QC_minfi <- minfiQC(Mset, fixOutliers = TRUE)

  pdf(file = here::here("replicate2","QC", "PotentialOutliers.pdf"))
  plotQC(QC_minfi$qc)
  dev.off()

  qc <- QC_minfi$qc
  
  ## Question: Also use 10 for FA??
  badSampleCutoff <- 10
  meds <- (qc$mMed + qc$uMed)/2
  whichBad <- which((meds < badSampleCutoff))
  whichBad
  # [1]  1  2  3 26 32
  RGset_filtered <- RGset[,-c(1,2,3,26,32)]
  Mset_filtered <- Mset[,-c(1,2,3,26,32)]
  pheno_filtered <- pheno[-c(1,2,3,26,32),]

  pheno$fa[c(1,2,3,26,32)]
  # [1] 1 1 1 0 1

  pheno$fa_ppm[c(1,2,3,26,32)]
  # [1] 0.318102981 1.999647696 0.877723577 0.008485281 1.081327913

  pheno$age[c(1,2,3,26,32)]
  # [1] 36.22177 36.12047 32.39699 33.55236 24.56674

  pheno$Subject_ID[c(1,2,3,26,32)]
  # [1] 1001 1003 1004 2046 1065

  ########################### with the sample filtering before
  detP_filtered <- detectionP(RGset_filtered)

  filtered_PS <- champ.filter(beta=getBeta(Mset_filtered),
  M=getM(Mset_filtered), pd = pheno_filtered, beadcount = NULL,
  detP = detP_filtered, arraytype = "450K")

  ############## replicate 1################
  # Filtering probes with a detection p-value above 0.01.
  # Removing 12446 probes.
  # If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
  # 
  # Filtering NoCG Start
  # Only Keep CpGs, removing 2688 probes from the analysis.
  # 
  # Filtering SNPs Start
  # Using general 450K SNP list for filtering.
  # Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
  # Removing 57129 probes from the analysis.
  # 
  # Filtering MultiHit Start
  # Filtering probes that align to multiple locations as identified in Nordlund et al
  # Removing 11 probes from the analysis.
  # 
  # Filtering XY Start
  # Filtering probes located on X,Y chromosome, removing 9479 probes from the analysis.
  # 
  # Updating PD file
  # 
  # Fixing Outliers Start
  # Replacing all value smaller/equal to 0 with smallest positive value.
  # Replacing all value greater/equal to 1 with largest value below 1..
  # [ Section 2: Filtering Done ]
  # 
  # All filterings are Done, now you have 403759 probes and 66 samples.

  ############ replicate 2 ###############
  # Filtering probes with a detection p-value above 0.01.
  # Removing 12640 probes.
  # If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
  # 
  # Filtering NoCG Start
  # Only Keep CpGs, removing 2664 probes from the analysis.
  # 
  # Filtering SNPs Start
  # Using general 450K SNP list for filtering.
  # Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
  # Removing 57107 probes from the analysis.
  # 
  # Filtering MultiHit Start
  # Filtering probes that align to multiple locations as identified in Nordlund et al
  # Removing 11 probes from the analysis.
  # 
  # Filtering XY Start
  # Filtering probes located on X,Y chromosome, removing 9469 probes from the analysis.
  # 
  # Updating PD file
  # 
  # Fixing Outliers Start
  # Replacing all value smaller/equal to 0 with smallest positive value.
  # Replacing all value greater/equal to 1 with largest value below 1..
  # [ Section 2: Filtering Done ]
  # 
  # All filterings are Done, now you have 403621 probes and 66 samples.
  # 
  # [<<<<< ChAMP.FILTER END >>>>>>]
  # [===========================]
  # [You may want to process champ.QC() next.]

  pheno_PS <- filtered_PS$pd
  betas_PS <- filtered_PS$beta
  mvals_PS <- filtered_PS$M

  filtered_CpG <- rownames(betas_PS)
  unfiltered_CpG <- rownames(betas)
  outCpG_PS <- unfiltered_CpG[!(unfiltered_CpG %in% filtered_CpG)]

  champ.QC(beta = betas_PS, pheno=pheno_PS$fa, mdsPlot=TRUE,
  densityPlot=TRUE, dendrogram=TRUE, PDFplot=TRUE, Rplot= FALSE,
  Feature.sel="None", resultsDir=here::here("replicate2", "QC",
  "postProbeandSampleFiltering"))

  champ.SVD(beta = betas_PS, rgSet=RGset_filtered, pd=pheno_PS,
  RGEffect=TRUE, Rplot=FALSE, resultsDir=here::here("replicate2", "QC",
  "postProbeandSampleFiltering", "SVD"))

  ######################### save the relevant data ############################

  # save(pheno, betas, RGset, Mset,
  #   pheno_filtered, RGset_filtered, Mset_filtered,
  #   pheno_PS, betas_PS, mvals_PS,
  #   mvals_P, betas_P, pheno_P, outCpG_P, outCpG_PS,
  # file = here("SuperFund/data", "raw_data.RData"))
  
  save(pheno, betas, RGset, Mset,
    pheno_filtered, RGset_filtered, Mset_filtered,
    pheno_PS, betas_PS, mvals_PS,
    mvals_P, betas_P, pheno_P, outCpG_P, outCpG_PS,
  file = here("replicate2", "raw_data_r2.RData"))

  write.csv(pheno, file=here("replicate2", "pheno_CGR_r2.csv"))
  write.csv(pheno_P, file=here("replicate2", "pheno_P_r2.csv"))
  write.csv(pheno_PS, file=here("replicate2", "pheno_PS_r2.csv"))

} else {
   load(raw_data)
}
