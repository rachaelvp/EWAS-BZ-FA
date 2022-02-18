library(here)
source(here("utils.R"))

#pheno <- read.csv(here("Desktop/SuperFund/data", "pheno_all.csv"))
pheno <- read.csv(here("replicate2", "pheno_all_r2.csv"))
sapply(pheno, class)
pheno <- run_class(pheno, fac = c(5,6,8,12,13), num = c(9,14:23))
rownames(pheno) <- pheno$Target_ID

# remove infection
design <- model.matrix(~fa+sex+age+Gran_est+Mono_est+Bcell_est+
                        NK_est+CD4T_est+CD8T_est, data = pheno)

# unadjusted design matrix
design_u <-  model.matrix(~fa, data = pheno)

load(here("replicate2", "combat_r2.RData"))

# bonferroni threshold
0.05/nrow(mvals_combat)


#################################### DMP #######################################
# Differential Mean DNA methylation 
DMP <- run_DMP(mvals = mvals_combat, design = design)
lambda(DMP[,5])



DMP_u <- run_DMP(mvals = mvals_combat, design = design_u)
lambda(DMP_u[,5])



DMP$Bonferroni <- p.adjust(DMP$P.Value, method = "bonferroni")


DMP_sig <- DMP %>%
filter(P.Value < .05)
dim(DMP_sig)



############ NOT after bonferrroni #############
########### ONLY after BH method ###############
DMP_annot <- annotate_probes(probeID_column = 1, sig_dat = DMP_sig)
DMP_annot <- DMP_annot[,-1]
DMP_annot <- dplyr::arrange(DMP_annot, P.Value)
# order by p.value
#DMP_annot <- DMP_annot[order(DMP_annot$P.Value),]
#write.csv(DMP_annot, here("Desktop/SuperFund/results","DMP.csv"),row.names = FALSE)
#DMP_annot <- read.csv(here("Desktop/SuperFund/results","DMP.csv"))

# EAS filtering
library(omicsPrint)
data(hm450.manifest.pop.GoNL)
dat <- data.frame(hm450.manifest.pop.GoNL)
rownames(dat) <- names(hm450.manifest.pop.GoNL)
maskname <- rownames(dat)[which(dat[,paste("MASK.general.","EAS",sep="")]==TRUE)]

DMP_CpG <- DMP_annot$ProbeID
length(unique(DMP_CpG[(DMP_CpG %in% maskname)]))

DMP_annot <- DMP_annot[!(DMP_annot$ProbeID %in% maskname),]
bonf <- dplyr::filter(DMP_annot, Bonferroni < 0.05)
DMP_bonf <- bonf$ProbeID
length(unique(DMP_bonf[(DMP_bonf %in% maskname)]))

# unadjusted
DMP_sig_u <- DMP_u %>%
  filter(P.Value < .05)
dim(DMP_sig_u)


DMP_sig <- DMP %>%
filter(adj.P.Val < .05)
dim(DMP_sig)
# [1] 0 7

# unadjusted 
DMP_sig_u <- DMP_u %>%
  filter(adj.P.Val < .05)
dim(DMP_sig_u)
# [1] 0 7



#################################### DVP #######################################
# Differential Variable DNA methylation 
DVP <- run_DVP(mvals = mvals_combat, design = design)
lambda(DVP[,6])


# unadjusted
DVP_u <- run_DVP(mvals = mvals_combat, design = design_u)
lambda(DVP_u[,6])

DVP$Bonferroni <- p.adjust(DVP$P.Value, method = "bonferroni")


DVP_sig <- DVP %>%
filter(P.Value < .05)
dim(DVP_sig)


# unadjusted 
DVP_sig_u <- DVP_u %>%
  filter(P.Value < .05)
dim(DVP_sig_u)

DVP_annot <- annotate_probes(probeID_column = 1, sig_dat = DVP_sig)
# EAS
library(omicsPrint)
data(hm450.manifest.pop.GoNL)
dat <- data.frame(hm450.manifest.pop.GoNL)
rownames(dat) <- names(hm450.manifest.pop.GoNL)
maskname <- rownames(dat)[which(dat[,paste("MASK.general.","EAS",sep="")]==TRUE)]

DVP_CpG <- DVP_annot$Probe_ID
length(unique(DVP_CpG[(DVP_CpG %in% maskname)]))


DVP_annot <- DVP_annot[!(DVP_annot$Probe_ID %in% maskname),]
bonf <- dplyr::filter(DVP_annot, Bonferroni < 0.05)
DVP_bonf <- bonf$Probe_ID
length(unique(DVP_bonf[(DVP_bonf %in% maskname)]))
# 0
DVP_annot <- dplyr::arrange(DVP_annot, P.Value)

DVP_sig <- DVP %>%
filter(Adj.P.Value < .05)
dim(DVP_sig)

# unadjusted
DVP_sig_u <- DVP_u %>%
  filter(Adj.P.Value < .05)
dim(DVP_sig_u)


############ 
DVP_annot <- annotate_probes(probeID_column = 1, sig_dat = DVP_sig)
DVP_annot <- DVP[order(DVP$Adj.P.Value),]

# unadjusted 
DVP_annot_u <- annotate_probes(probeID_column = 1, sig_dat = DVP_sig_u)

#################################### DMR #######################################
# DNA methylation across small regions
## 1000 bootstraps
# run on terminal 
DMR <- run_DMR(mvals = mvals_combat, design = design, pheno = pheno)

DMR_sig <- DMR %>%
filter(p.value < .05)
dim(DMR_sig)

DMR_annot <- annotate_regions(DMR_sig)
DMR_annot <- DMR_annot[order(DMR_annot$p.value),]


# family-wise error rate 
DMR_sig <- DMR %>%
filter(fwer < .05)
dim(DMR_sig)


#################################### DMB #######################################
# large, cluster regions of CpG probes, blocks 
load(here("data","normalized.RData"))
CN <- assays(normalized)$CN
filtered_CpG <- rownames(mvals_combat)
unfiltered_CpG <- rownames(CN)
outCpG <- unfiltered_CpG[!(unfiltered_CpG %in% filtered_CpG)]
CN <- CN[!rownames(CN) %in% outCpG,]
DMB <- run_DMB(mvals = mvals_combat, design = design, pheno = pheno,CN = CN)


DMB_sig <- DMB %>%
filter(p.value < .05)
dim(DMB_sig)

#### Annotate DMB ##########
DMB_annot <- annotate_regions(DMB_sig)
DMB_annot <- DMB_annot[order(DMB_annot$p.value),]

# family-wise error rate 
DMB_sig <- DMB %>%
filter(fwer < .05)
dim(DMB_sig)
