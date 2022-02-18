library(here)
library(data.table)
library(tidyverse)
library(limma)
library(minfi)
library(missMethyl)
library(doParallel)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(DMRcate)
library(tidyverse)
library(org.Hs.eg.db)
library(omicsPrint)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(qvalue)
library(qqman)
library(ggplot2)
library(reshape2)
source(here("BZ-FA-450K", "benzene", "utils.R"))

load(here("BZ-FA-450K", "benzene", "data", "mvals.RData"))
pheno <- read.csv(
  here("BZ-FA-450K", "benzene", "data", "pheno_noduplicates.csv")
)
rownames(pheno) <- pheno$target_ID
all(rownames(pheno) == colnames(mvals))

0.05/nrow(mvals)
# [1] 1.227774e-07

design <- model.matrix(
  ~benz1mo+sex+smk+age+bmi+Gran_est+Mono_est+Bcell_est+NK_est+CD4T_est+CD8T_est, 
  data = pheno
)
################################################################################
# DMP 
################################################################################
# lambda < 1 == deflated, under-powered, over adjustment (larger p-values)
# more of uniform distribution
DMP <- run_DMP(mvals = mvals, design = design)
save(DMP, here("BZ-FA-450K", "benzene", "results", "DMP.RData"), compress = TRUE)
lambda(DMP[,5]) # 1.048371
DMP$Bonferroni <- p.adjust(DMP$P.Value, method = "bonferroni")
nrow(dplyr::filter(DMP, P.Value < .05)) # 23,609
DMP_sig <- dplyr::filter(DMP, adj.P.Val < .05)
length(unique(DMP_sig$ProbeID)) # 61
nrow(dplyr::filter(DMP, Bonferroni < .05)) # 26
DMP_annotated <- annotate_and_filterEAS(DMP_sig)
# [1] "Number probes did not annotate to genome: 0"
# [1] "Number probes did not annotate to gene: 11"
# [1] "Total BH-sig EAS probes: 6"
# [1] "Total Bonferroni-sig EAS probes: 4"
length(unique(DMP_annotated$ProbeID)) # 55

# Manhattan plot
results <- data.frame(DMP)
rownames(results) <- results$ProbeID
results$pvalue <- results$P.Value
plot_manhattan(
  results, here("BZ-FA-450K", "benzene", "results", "manhattanplot_DMP.pdf")
)

################################################################################
# DVP 
################################################################################
DVP <- run_DVP(mvals = mvals, design = design)
save(DVP, here("BZ-FA-450K", "benzene", "results", "DVP.RData"), compress = TRUE)
DVP <- DVP[-is.na(DVP$LogVarRatio), ]
lambda(DVP[,6]) # 0.7003848
DVP$Bonferroni <- p.adjust(DVP$P.Value, method = "bonferroni")
nrow(dplyr::filter(DVP, P.Value < .05)) # 12,477
DVP_sig <- dplyr::filter(DVP, Adj.P.Value < .05)
length(unique(DVP_sig$Probe_ID)) # 1,688
nrow(dplyr::filter(DVP, Bonferroni < .05)) # 324
DVP_annotated <- annotate_and_filterEAS(DVP_sig, TRUE)
# [1] "Number probes did not annotate to genome: 0"
# [1] "Number probes did not annotate to gene: 343"
# [1] "Total BH-sig EAS probes: 32"
# [1] "Total Bonferroni-sig EAS probes: 7"
length(unique(DVP_annotated$ProbeID)) # 1,656

# Manhattan plot
results <- data.frame(DVP)
rownames(results) <- results$Probe_ID
results$pvalue <- results$P.Value
plot_manhattan(
  results, here("BZ-FA-450K", "benzene", "results", "manhattanplot_DVP.pdf")
)

keep <- c(
  "ProbeID","SampleVar","LogVarRatio","DiffLevene","t","P.Value", "Adj.P.Value",
  "Bonferroni","chr","pos","strand","Type","UCSC_RefGene_Name","DMR","Enhancer",
  "DHS"
)
DVP_annotated <- data.table(DVP_annotated)
to_keep <- !duplicated(DVP_annotated[,keep,with=F])
DVP_annotated <- DVP_annotated[to_keep,]

################################################################################
# DVP+DMP overlap
################################################################################
DMP_DVP <- DMP_annotated[which(DMP_annotated$ProbeID %in% DVP_annotated$ProbeID),]
unique(DMP_DVP$ProbeID)
# [1] "cg17670477" "cg14075413" "cg17905084" "cg21394778" "cg20536921"
# [6] "cg13459303" "cg18552413" "cg16619049" "cg14051111" "cg06169961"
# [11] "cg20159193" "cg04759112" "cg25608626" "cg06713830" "cg06530725"
# [16] "cg08175635" "cg16211055" "cg07832006" "cg25114611" "cg19539385"
# [21] "cg11822739" "cg05668674" "cg02938413" "cg08781710" "cg14254433"
# [26] "cg03992168" "cg07987250" "cg05308829" "cg18380742" "cg23673397"
# [31] "cg15924831" "cg13869366"
length(unique(DMP_DVP$ProbeID)) # 32
DMP_annotated$DVP_Overlap <- ifelse(
  DMP_annotated$ProbeID %in% DMP_DVP$ProbeID, "Yes", ""
)
keep <- c(
  "ProbeID", "logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B", "Bonferroni",
  "chr","pos","strand","Type","UCSC_RefGene_Name", "DMR","Enhancer", "DHS", 
  "DVP_Overlap"
)
DMP_annotated <- data.table(DMP_annotated)
to_keep <- !duplicated(DMP_annotated[, keep, with=F])
DMP_annotated <- DMP_annotated[to_keep,]

################################################################################
# DMR 
################################################################################
DMR <- run_DMR_bumphunter(
  mvals = mvals, design = design, pheno = pheno, num_cores = 20
)
# cutoff: 0.007
# Found 1,202 bumps.
nrow(dplyr::filter(DMR, p.value < .05)) # 168
nrow(dplyr::filter(DMR, fwer < .05)) # 0
DMR_annotated <- annotate_regions(DMR)
DMR_annot <- DMR_annotated[,c(1:3,17,4:16,18:20)]
write.csv(
  DMR_annot, row.names = FALSE,
  file = here("BZ-FA-450K", "benzene", "results", "DMR.csv"), 
)
DMR_annot <- dplyr::filter(DMR_annot, p.value < .05)

################################################################################
# DMB 
################################################################################
load(here("BZ-FA-450K", "benzene", "data", "normalized_assays.RData"))
CN <- assays(normalized_assays)$CN
filtered_CpG <- rownames(mvals)
unfiltered_CpG <- rownames(CN)
outCpG <- unfiltered_CpG[!(unfiltered_CpG %in% filtered_CpG)]
CN <- CN[!(rownames(CN) %in% outCpG), (colnames(CN) %in% colnames(mvals))]

DMB <- run_DMB(
  mvals = mvals, design = design, pheno = pheno, CN = CN, num_cores = 20
)
# cutoff: 0.003
# Found 440 bumps.
nrow(dplyr::filter(DMB, p.value < .05)) # 70
nrow(dplyr::filter(DMB, fwer < .05)) # 1
DMB_annotated <- annotate_regions(DMB)
DMB_annot <- DMB_annotated[,c(1:3,17,4:16,18:20)]
write.csv(
  DMB_annot, row.names = FALSE,
  file = here("BZ-FA-450K", "benzene", "results", "DMB.csv")
)
DMB_annot <- dplyr::filter(DMB_annot, p.value < .05)
write.csv(
  DMB_annot, row.names = FALSE,
  here("BZ-FA-450K", "benzene", "results", "DMB_significant.csv")
)
################################ gene overlap ##################################
DMB_genes <- DMB_annot$gene_symbol
DMR_genes <- DMR_annot$gene_symbol
DVP_genes <- DVP_annotated$UCSC_RefGene_Name
DMP_genes <- DMP_annotated$UCSC_RefGene_Name

DMP_DVP_genes <- unique(DMP_genes[which(DMP_genes %in% DVP_genes)])
DMP_DVP_genes
# [1] "TBCD"        ""            "SERPINA4"    "FCRLB"       "PRDM16"
# [6] "TMEM176A"    "TMEM176B"    "DARC"        "FAM41C"      "PTPRN2"
# [11] "C10orf26"    "NUDT3"       "RGS12"       "CMIP"        "DSCAML1"
# [16] "PPAN"        "PPAN-P2RY11" "SNORD105B"   "DIO1"        "ZBTB46"
# [21] "TRIM11"      "SYCP1"       "LOC285847"   "FKBP5"       "EIF2AK1"
# [26] "STX16"       "GPR56"       "ZAN"         "MEFV"        "PACSIN1"
# [31] "GNA12"       "NFIA"        "NXN"         "FLJ23867"    "QSOX1"
# [36] "C1orf159"    "PDGFRB"      "OR5E1P"
DMP_DVP_genes <- DMP_DVP_genes[-2]

DMP_DMR_genes <- unique(DMP_genes[which(DMP_genes %in% DMR_genes)])
DMP_DMR_genes
# "PTPRN2"

DMP_DMB_genes <- unique(DMP_genes[which(DMP_genes %in% DMB_genes)])
DMP_DMB_genes
# character(0)

DMR_DVP_genes <- unique(DMR_genes[which(DMR_genes %in% DVP_genes)])
DMR_DVP_genes
# [1] "UNC45A" "MUC4"   "CIDEB"  "LTB4R2" "BRF1"   "BTBD6"  "F11R"   "EHHADH"
# [9] "TXNRD1" "PTPRN2" "GTDC1"  "ACTA1"

DMB_DVP_genes <- unique(DMB_genes[which(DMB_genes %in% DVP_genes)])
DMB_DVP_genes
# [1] "CCDC130"  "FUT10"    "STK35"    "TMEM155"  "CCNA2"    "SLC23A2"
# [7] "INSR"     "YTHDF1"   "SNORD50A" "SNHG5"    "SNORD50B"
DMR_DMB_genes <- unique(DMR_genes[which(DMR_genes %in% DMB_genes)])
DMR_DMB_genes
# [1] NA  "ARPP21" "TRPC3"
DMR_DMB_genes <- DMR_DMB_genes[-1]

DMP_annotated$DVPgene_Overlap <- ifelse(
  DMP_annotated$UCSC_RefGene_Name %in% DMP_DVP_genes, "Yes", ""
)
DMP_annotated$DMR_Overlap <- ifelse(
  DMP_annotated$UCSC_RefGene_Name %in% DMP_DMR_genes, "Yes", ""
)
DMP_annotated$DMB_Overlap <- ifelse(
  DMP_annotated$UCSC_RefGene_Name %in% DMP_DMB_genes, "Yes", ""
)
write.csv(
  DMP_annotated, 
  file = here("BZ-FA-450K", "benzene", "results", "DMP_BHsignificant.csv"), 
  row.names = FALSE
)

DVP_annotated$DMR_Overlap <- ifelse(
  DVP_annotated$UCSC_RefGene_Name %in% DMR_DVP_genes, "Yes", "")
DVP_annotated$DMB_Overlap <- ifelse(
  DVP_annotated$UCSC_RefGene_Name %in% DMB_DVP_genes, "Yes", "")
write.csv(
  DVP_annotated, 
  file = here("BZ-FA-450K", "benzene", "results", "DVP_BHsignificant.csv"), 
  row.names = FALSE
)

DMR_annot$DMB_Overlap <- ifelse(DMR_annot$gene_symbol %in% DMR_DMB_genes, "Yes", "")
write.csv(
  DMR_annot, 
  file = here("BZ-FA-450K", "benzene", "results", "DMR_significant.csv"), 
  row.names = FALSE
)
