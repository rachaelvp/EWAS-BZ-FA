library(limma)
library(minfi)
library(missMethyl)
library(doParallel)
library(GenomicFeatures)
library(org.Hs.eg.db)
#library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

run_DMP <- function(mvals, design){
  lFit <- limma::lmFit(object = mvals, design = design)
  eFit <- limma::eBayes(lFit, robust = TRUE)
  pTop <- limma::topTable(eFit, coef = 2, num = Inf, sort.by = "p")
  pTop <- pTop[order(pTop$P.Value), , drop = FALSE]
  ProbeID <- rownames(pTop)
  DMP <- data.frame(ProbeID, pTop)
  return(DMP)
}

run_DVP <- function(mvals, design){
  fitvar <- varFit(mvals, design = design, coef = c(1,2))
  topDV <- topVar(fitvar, coef=2, number = 399439)
  topDV <- data.frame(Probe_ID = rownames(topDV), topDV)
  return(topDV)
}

# change bootstrap iter to 200
run_DMR <- function(mvals, design, pheno){
  registerDoParallel(cores = 1)
  GRset <- makeGenomicRatioSetFromMatrix(mat = mvals,
    pData = pheno, array = "IlluminaHumanMethylation450k",
    annotation = "ilmn12.hg19", what = "M")
  Bumphunter <- minfi::bumphunter(object = GRset, design = design, coef = 2,
    B = 300, type = "M", pickCutoff=TRUE, nullMethod="bootstrap",
    maxGap=1000, smooth = TRUE, smoothFunction = loessByCluster)
  dmr_results <- Bumphunter$table
  return(dmr_results)
}

# run B = 500 (test)
run_DMB <- function(mvals, design, pheno, CN){
  registerDoParallel(cores = 1)
  GRset <- makeGenomicRatioSetFromMatrix(mat = mvals,
    pData = pheno, array = "IlluminaHumanMethylation450k",
    annotation = "ilmn12.hg19", what = "M")
  assays(GRset)$CN <- CN
  cluster <- cpgCollapse(GRset, what = "M")
  blocks <- minfi::blockFinder(cluster$object, design, coef = 2,
  what = "M", nullMethod = "bootstrap", smooth = TRUE, cluster = NULL,
  pickCutoff=TRUE, B = 500, smoothFunction = loessByCluster)
  dmb_results <- blocks$table
  return(dmb_results)
}

lambda <- function(p) median(qchisq(p, df=1, lower.tail=FALSE),
na.rm=TRUE) / qchisq(0.5, df=1)

## ============== Where to get annotation data ============ ##

annotate_probes <- function(probeID_column, sig_dat){
  load(here("SuperFund/data", "450k_annot.RData"))
  sig_probes <- as.character(sig_dat[,probeID_column])
  sig_dat$Probe_ID <- sig_dat[,probeID_column]
  annot <- annotation[(annotation$Probe_ID %in% sig_probes),]
  dat <- merge(sig_dat, annot, by = "Probe_ID")
  return(dat)
}

# annote DMRs/DMBs with genes
annotate_regions <- function(regions_data) {
  GR <- makeGRangesFromDataFrame(regions_data)
  GR1 <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), GR)
  sym <- as.data.frame(org.Hs.egSYMBOL)
  merged <- merge(GR1, sym, by = "gene_id")
  GR2 <- mergeByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), GR)
  m <- DataFrame(geneid = GR2[,2], GR_region = GR2[,3])
  m$gene_id <- make.names(m$geneid, unique=TRUE)
  m$gene_id <- substring(m$gene_id, 2)
  #m1 <- names(m$GR_region)
  #names(m$GR_region) <- make.names(m1, unique=TRUE)
  merged2 <- merge(m, merged,by = "gene_id")
  annot <- data.frame(merged2)
  #annot <- data.frame(merged)
  annot1 <- annot[,c(2:5,13,6,9:11)]
  colnames(annot1) <- c("gene_id","chr", "start","end","gene_symbol","width",
                        "gene_start","gene_end","gene_width")
  join <- left_join(regions_data, annot1, by = c("chr","start","end"))
  join$width <- join$end-join$start
  join1 <- filter(join, width != 0)
  return(join1)
}


######## run_class ########

library(tidyverse)
# smk col = 11, rm smk 
run_class <- function(pheno,fac = c(5,6,8,12,13), num = c(9,14:23)) {
  cols_all <- c(colnames(pheno))
  pheno[cols_all] <- lapply(pheno[cols_all], as.character)
  cols_num <- c(colnames(pheno[num]))
  pheno[cols_num] <- lapply(pheno[cols_num], as.numeric)
  cols_fac <- c(colnames(pheno[fac]))
  pheno[cols_fac] <- lapply(pheno[cols_fac], as.factor)
  return(pheno)
}
