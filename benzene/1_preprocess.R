library(here)
library(tidyverse)
library(minfi)
library(ENmix)
library(ChAMP)
library(sva)

################################################################################
# Import RG Channel Set to match with phenotype information
################################################################################
pheno <- read.csv(here("BZ-FA-450K", "benzene", "data", "pheno.csv"))
dim(pheno)
# [1] 107  30

load(here("BZ-FA-450K", "benzene", "data", "RGset.RData"))
dim(RGset)
# [1] 622399    107

# crucial that the methylation data columns and phenotype rows correspond
rownames(pheno) <- pheno$target_ID
all(colnames(RGset) == rownames(pheno))

# see if we need to remove any samples
QC_minfi <- minfiQC(preprocessRaw(RGset), fixOutliers = TRUE)
pdf(file = here("BZ-FA-450K", "benzene", "QC", "PotentialOutliers.pdf"))
plotQC(QC_minfi$qc)
dev.off()
# no bad samples

################################################################################
# Normalize the data and estimate the counts for each cell type
################################################################################

cellcounts <- estimateCellCounts(
  RGset, referencePlatform = "IlluminaHumanMethylation450k", 
  returnAll = TRUE, sex = pheno$sex
)

colnames(cellcounts$counts) <- c(
  "CD8T_est","CD4T_est","NK_est","Bcell_est", "Mono_est","Gran_est"
)
counts_df <- data.frame(cellcounts$counts)
counts_df$target_ID <- rownames(counts_df)
pheno <- left_join(pheno, counts_df)
write.csv(
  pheno, file = here("BZ-FA-450K", "benzene", "data", "pheno_all.csv"), 
  row.names = FALSE
)

normalized_assays <- cellcounts$normalizedData
normalized_betas <- M2B(assays(normalized_assays)$M)
normalized_betas <- Harman::shiftBetas(normalized_betas, shiftBy = 1e-4)
normalized_mvals <- B2M(newBetas)
assays(normalized_assays)$Beta <- normalized_betas
assays(normalized_assays)$M <- normalized_mvals

all(colnames(RGset) == rownames(pheno))
all(colnames(assays(normalized_assays)$M) == rownames(pheno))
all(colnames(assays(normalized_assays)$Beta) == rownames(pheno))

save(
  normalized_assays, 
  file = here("BZ-FA-450K", "benzene", "data", "normalized_assays.RData"),
  compress = TRUE
)
################################################################################
# Filter probes
################################################################################

# detP - so we can filter the CpGs that didn't successfully hybridize to array
detP <- detectionP(RGset)
detP <- detP[match(rownames(assays(normalized_assays)$Beta), rownames(detP)),]
all(rownames(assays(normalized_assays)$Beta) == rownames(detP))

Mset <- preprocessRaw(RGset)
Mset <- Mset[match(rownames(assays(normalized_assays)$Beta), rownames(Mset)),]
all(rownames(assays(normalized_assays)$Beta) == rownames(M))

filtered <- champ.filter(
  beta = assays(normalized_assays)$Beta, 
  M = assays(normalized_assays)$M, 
  intensity = Mset, pd = pheno, detP = detP
)

# Filtering probes with a detection p-value above 0.01.
# Removing 8646 probes.
# If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples
#
# Filtering NoCG Start
# Only Keep CpGs, removing 2797 probes from the analysis.
#
# Filtering SNPs Start
# Using general 450K SNP list for filtering.
# Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
# Removing  57190 probes from the analysis.
#
# Filtering MultiHit Start
# Filtering probes that align to multiple locations as identified in Nordlund et al
# Removing 11 probes from the analysis.
#
# Filtering XY Start
# Filtering probes located on X,Y chromosome, removing 9627 probes from the analysis.
#
# Updating PD file
#
# Fixing Outliers Start
# Replacing all value smaller/equal to 0 with smallest positive value.
# Replacing all value greater/equal to 1 with largest value below 1..
#
# All filterings are Done, now you have 407241 probes and 107 samples.

betas <- filtered$beta
mvals <- filtered$M
Mset <- filtered$intensity

champ.QC(
  beta = betas, pheno = pheno$benz_dich, mdsPlot = TRUE,
  densityPlot = TRUE, dendrogram = TRUE, PDFplot = TRUE, Rplot = FALSE,
  Feature.sel = "None", resultsDir = here("BZ-FA-450K", "benzene", "QC")
)
champ.SVD(
  beta = betas, rgSet = rgsetAll, pd = pheno, RGEffect = TRUE,
  Rplot = FALSE, resultsDir = here("BZ-FA-450K", "benzene", "QC")
)

##############################################################################
# Remove duplicated samples
# Decision based on probes quality, as measured by detectionP
##############################################################################

dups <- which(duplicated(pheno[,9:20]))
dups2 <- which(duplicated(pheno[,9:20], fromLast = TRUE))

detP <- detectionP(RGset)
detP_relevant <- detP[which(rownames(detP) %in% rownames(mvals)), ]
nrow(detP_relevant) == nrow(mvals)
detP_sum <- colSums(detP_relevant)
sum_detP_tbl <- data.frame(sum_duplicate1 = detP_sum[dups], 
                           sum_duplicate2 = detP_sum[dups2])
dups_tbl <- data.frame(dups, dups2)
dups_tbl$min_detP_dup <- apply(sum_detP_tbl, 1, which.min)
worse_detP <- ifelse(dups_tbl$min_detP_dup == 2, dups_tbl$dups, dups_tbl$dups2)
worse_detP <- c(5,13,26,28,36,42,50,84,90)
pheno_nodups <- pheno[-worse_detP,]  
table(pheno_nodups$benz_dich)
# 0  1
# 48 50                                          
mvals_nodups <- mvals[, -worse_detP]

##############################################################################
# ComBat
# The ComBat function adjusts for known batches using an empirical Bayesian
# framework. ComBat allows users to adjust for batch effects in datasets where
# the batch cov is known, using methodology described in Johnson et al. 2007.
##############################################################################
batch <- pheno_nodups$sample_plate
mod <- model.matrix(~benz1mo, data = pheno_nodups)
mvals <- ComBat(dat = mvals_nodups, batch, mod)
betas <- M2B(mvals)
betas <- Harman::shiftBetas(betas, shiftBy = 1e-4)
mvals <- B2M(betas)

save(
  mvals, compress = TRUE,
  file = here("BZ-FA-450K", "benzene", "data", "mvals.RData")
)
write.csv(
  pheno_nodups, row.names = FALSE,
  file = here("BZ-FA-450K", "benzene", "data", "pheno_noduplicates.csv")
)

# examine SVD plot after batch effect correction
champ.SVD(
  beta = betas, rgSet=rgsetAll[,-worse_detP], pd = pheno_nodups, 
  RGEffect = TRUE, Rplot = FALSE, 
  resultsDir = here("BZ-FA-450K", "benzene", "QC", "postCombat_")
)