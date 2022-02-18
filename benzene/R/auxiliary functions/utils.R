################################################################################
# analyses across individual probes (aka CpG site)
################################################################################

# identify PROBES exhibiting differential mean methylation
run_DMP <- function(mvals, design) {
  lFit <- limma::lmFit(object = mvals, design = design)
  eFit <- limma::eBayes(lFit, robust = TRUE)
  pTop <- limma::topTable(eFit, coef = 2, num = Inf, sort.by = "p")
  pTop <- pTop[order(pTop$P.Value), , drop = FALSE]
  ProbeID <- rownames(pTop)
  return(data.frame(ProbeID, pTop))
}

# identify PROBES exhibiting differential variance of methylation
run_DVP <- function(mvals, design) {
  fitvar <- varFit(mvals, design = design, coef=NULL, type = "SQ")
  topDV <- topVar(fitvar, coef=2, number = nrow(mvals))
  return(data.frame(Probe_ID = rownames(topDV), topDV))
}

# calculate the enrichment score (hopefully ~ 1) for DVP and DMP p-values
lambda <- function(p){
  median(qchisq(p, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, df = 1)
}

# annote DVPs/DMPs with genes and other information (i.e., CpG islands, SNPs)
load(here("BZ-FA-450K", "benzene", "data", "450Kannotation.RData"))
annotate_probes <- function(probeID_column, sig_dat){
  sig_probes <- as.character(sig_dat[,probeID_column])
  sig_dat$Probe_ID <- sig_dat[,probeID_column]
  annot <- annotation[(annotation$Probe_ID %in% sig_probes),]
  missing <- length(sig_probes[!(sig_probes %in% annot$Probe_ID)])
  print(paste0("Number probes did not annotate to genome: ", missing))
  dat <- merge(sig_dat, annot, by = "Probe_ID", all.x=T)
  CpGs_nogenes <- unique(dat[which(dat$UCSC_RefGene_Name == ""),1])
  print(paste0("Number probes did not annotate to gene: ", length(CpGs_nogenes)))
  return(dat)
}

# EAS filtering
data(hm450.manifest.pop.GoNL)
dat <- data.frame(hm450.manifest.pop.GoNL)
rownames(dat) <- names(hm450.manifest.pop.GoNL)
maskname <- rownames(dat)[which(dat[,paste("MASK.general.","EAS",sep="")]==TRUE)]
annotate_and_filterEAS <- function(BHsig_probe_res, DVP = FALSE){
  if(DVP){
    col_to_rename <- which(colnames(BHsig_probe_res) == "Probe_ID")
    colnames(BHsig_probe_res)[col_to_rename] <- "ProbeID"
  }
  
  DMP_annotated <- annotate_probes(probeID_column = 1, sig_dat = BHsig_probe_res)
  if(!DVP){
    DMP_annotated <- DMP_annotated[,-1]
  }
  DMP_annotated <- dplyr::arrange(DMP_annotated, P.Value)
  
  DMP_CpG <- DMP_annotated$ProbeID
  BH_EAS_probes <- length(unique(DMP_CpG[(DMP_CpG %in% maskname)]))
  print(paste0("Total BH-sig EAS probes: ", BH_EAS_probes))
  
  bonf <- dplyr::filter(DMP_annotated, Bonferroni < 0.05)
  DMP_bonf <- bonf$ProbeID
  Bonf_EAS_probes <- length(unique(DMP_bonf[(DMP_bonf %in% maskname)]))
  print(paste0("Total Bonferroni-sig EAS probes: ", Bonf_EAS_probes))
  
  return(DMP_annotated[!(DMP_annotated$ProbeID %in% maskname),])
}

################################################################################
# analyses across regions (i.e. across many probes)
################################################################################

# identify REGIONS exhibiting differential mean methylation with bumphunter
run_DMR_bumphunter <- function(mvals, design, pheno, num_cores = 1, 
                               num_permutations = 1000) {
  GRset <- makeGenomicRatioSetFromMatrix(mat = mvals,
    pData = pheno, array = "IlluminaHumanMethylation450k",
    annotation = "ilmn12.hg19", what = "M")
  registerDoParallel(cores = num_cores)
  Bumphunter <- minfi::bumphunter(object = GRset, design = design, coef = 2,
    B = num_permutations, type = "M", pickCutoff=TRUE, nullMethod="bootstrap",
    maxGap=1000, smooth = TRUE, smoothFunction = loessByCluster)
  dmr_results <- Bumphunter$table
  return(dmr_results)
}

# identify BLOCKS (aka huge regions) exhibiting differential mean methylation
run_DMB <- function(mvals, design, pheno, CN, num_cores = 1, 
                    num_permutations = 1000) {
  GRset <- makeGenomicRatioSetFromMatrix(
    mat = mvals, pData = pheno, array = "IlluminaHumanMethylation450k",
    annotation = "ilmn12.hg19", what = "M"
  )
  assays(GRset)$CN <- CN
  cluster <- cpgCollapse(GRset, what = "M", maxGap = 1000)
  registerDoParallel(cores = num_cores)
  blocks <- minfi::blockFinder(
    cluster$object, design, coef = 2, what = "M", nullMethod = "bootstrap", 
    smooth = TRUE, pickCutoff=TRUE, B = num_permutations, 
    smoothFunction = loessByCluster
  )
  dmb_results <- blocks$table
  return(dmb_results)
}

# annote DMRs/DMBs with genes
annotate_regions <- function(regions_data) {
  GR <- makeGRangesFromDataFrame(regions_data)
  GR1 <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), GR)
  sym <- as.data.frame(org.Hs.egSYMBOL)
  merged <- merge(GR1, sym, by = "gene_id")
  GR2 <- mergeByOverlaps(
    genes(TxDb.Hsapiens.UCSC.hg19.knownGene,single.strand.genes.only = FALSE), 
    GR
  )
  m <- DataFrame(geneid = GR2[,2], GR_region = GR2[,3])
  m$gene_id <- make.names(m$geneid, unique=TRUE)
  m$gene_id <- substring(m$gene_id, 2)
  m1 <- names(m$GR_region)
  names(m$GR_region) <- make.names(m1, unique=TRUE)
  merged2 <- merge(m, merged,by = "gene_id")
  annot <- data.frame(merged2)
  annot1 <- annot[,c(2:5,13,6,9:11)]
  colnames(annot1) <- c(
    "gene_id","chr", "start","end","gene_symbol","width","gene_start","gene_end",
    "gene_width"
  )
  join <- left_join(regions_data, annot1, by = c("chr","start","end"))
  join$width <- join$end-join$start
  join1 <- filter(join, width != 0)
  return(join1)
}


################################################################################
# EAS SNP utility functions
################################################################################
get_population_SNPs <- function(population = "EAS"){
  require(omicsPrint)
  data(hm450.manifest.pop.GoNL)
  d <- data.frame(hm450.manifest.pop.GoNL)
  rownames(d) <- names(hm450.manifest.pop.GoNL)
  rownames(d)[which(d[,paste("MASK.general.",population,sep="")]==TRUE)]
}

filter_SNPs <- function(methylation_data, SNPs){
  SNP_probes <- which(rownames(methylation_data) %in% SNPs)
  methylation_data[-SNP_probes,]
}

################################################################################
# 450K annotation 
################################################################################
# data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# annotation <- data.frame(annotation)
# annotation <- data.table(annotation)
# annotation$Probe_ID <- as.character(annotation$Name)
# keep_cols <- c("Probe_ID","chr", "pos", "strand", "Type", "Color","Probe_rs",
#                "Probe_maf", "CpG_rs", "CpG_maf", "Islands_Name",
#                "Relation_to_Island", "Random_Loci", "UCSC_RefGene_Name", 
#                "UCSC_RefGene_Group", "Islands_Name", "Relation_to_Island", 
#                "Phantom", "DMR", "Enhancer", "DHS", "Regulatory_Feature_Name", 
#                "Regulatory_Feature_Group")
# annotation <- annotation[, keep_cols, with=F]
# annotation <- data.frame(annotation)
# annotation <- annotation %>% tidyr::separate_rows(UCSC_RefGene_Name, sep = ";")
# save(annotation, 
#      file = here("BZ-FA-450K", "benzene", "data", "450Kannotation.RData")
# )

# https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf

################################################################################
# miscellaneous 
################################################################################
run_class <- function(pheno, fac, num) {
  cols_all <- c(colnames(pheno))
  pheno[cols_all] <- lapply(pheno[cols_all], as.character)
  cols_num <- c(colnames(pheno[num]))
  pheno[cols_num] <- lapply(pheno[cols_num], as.numeric)
  cols_fac <- c(colnames(pheno[fac]))
  pheno[cols_fac] <- lapply(pheno[cols_fac], as.factor)
  return(pheno)
}

################################################################################
# plotting
################################################################################

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
IlluminaAnnot450K <-getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
IlluminaAnnot450K <- as.data.frame(IlluminaAnnot450K)

#'# Manhattan Plot
# results is a data.frame with Probe IDs as rownames and a column titled "pvalue"
plot_manhattan <- function(results, save_path){
  annot <- IlluminaAnnot450K[rownames(IlluminaAnnot450K) %in% rownames(results),]
  results <- results[rownames(annot) %in% rownames(results),]
  results <- results[sort(rownames(results)),]
  annot <- annot[sort(rownames(annot)),]
  
  #' Needs to be TRUE
  stopifnot(identical(rownames(annot),rownames(results)))
  stopifnot(all(rownames(annot) == rownames(results)))
  
  #' q-value and reduce annotation
  qval <- qvalue::qvalue(results$pvalue)
  qThresh <- max(results$pvalue[qval$qv <= 0.05])
  
  #' Only these elements are needed from annotation
  AnnotMan <- annot[,c("Name","chr","pos")]
  stopifnot(identical(rownames(AnnotMan),rownames(results))) 
  stopifnot(all(AnnotMan$Name==rownames(results))) 
  stopifnot(all(AnnotMan$Name==results$Probe_ID)) 
  results$Name <- rownames(results)
  stopifnot(all(results$Name==AnnotMan$Name))
  resultsMan_AnnotMan <- cbind(AnnotMan,results)
  stopifnot(all(resultsMan_AnnotMan$Name==resultsMan_AnnotMan$Probe_ID))
  
  #'  CHR should be numeric
  resultsMan_AnnotMan$CHR <- sub('.*chr', '',resultsMan_AnnotMan$chr)
  resultsMan_AnnotMan$CHR <- as.numeric(resultsMan_AnnotMan$CHR)
  # Colors
  variable <- as.factor(resultsMan_AnnotMan$CHR)
  myColors <- rainbow(length(table(variable)))
  names(myColors) <- levels(variable)
  
  nCpG <- length(results$LogFC)
  
  pdf(save_path)
  par(mfrow=c(1,1), mar = c(5,5,1,2), oma = c(0,0,0,0))
  qqman::manhattan(resultsMan_AnnotMan, chr = "CHR", bp = "pos", 
                   p = "pvalue", snp = "Name", genomewideline=-log10(0.05/nCpG),
                   suggestiveline=-log10(qThresh), ylim=range(-log10(resultsMan_AnnotMan$pvalue)), 
                   yaxt='n')
  axis(2, cex=1)
  dev.off()
}


# colnames of betas needs to be rownames pheno
make_DMR_tbl <- function(betas, pheno, gene_name, DMR_start, DMR_end,
                         exposure_colname = "benz_dich",
                         exposure_levels = c(0,1)){
  match_idx <- unlist(lapply(gene_name, function(x) grep(x, IlluminaAnnot450K$UCSC_RefGene_Name)))
  matches <- unique(IlluminaAnnot450K$UCSC_RefGene_Name[unique(match_idx)])
  gene <- dplyr::filter(IlluminaAnnot450K, UCSC_RefGene_Name %in% matches)
  gene <- dplyr::select(gene, "pos")
  betas_gene <- betas[rownames(betas) %in% rownames(gene),]
  betas_annot <- merge(betas_gene, gene, by = "row.names")
  rownames(betas_annot) <- betas_annot$Row.names
  betas_annot <- betas_annot[,-which(colnames(betas_annot) == "Row.names")]
  
  DMR_annot <- betas_annot[which(betas_annot$pos >= DMR_start & betas_annot$pos <= DMR_end),]
  DMR_annot_melted <- melt(DMR_annot, id.var = "pos")
  
  pheno_relevant <- dplyr::select(pheno, all_of(exposure_colname))
  pheno_relevant$variable <- rownames(pheno_relevant)
  
  merged <- merge(DMR_annot_melted, pheno_relevant, by = "variable", all = TRUE)
  colnames(merged) <- c("TargetID", "Position", "Value", "Exposure")
  
  site_pos <- data.frame("ProbeID" = rownames(DMR_annot), "Position" = DMR_annot$pos)
  plot_tbl <- merge(merged, site_pos, by = "Position", all = TRUE)
  plot_tbl <- dplyr::arrange(plot_tbl, Position)
  plot_tbl$Position_ProbeID <- paste(plot_tbl$Position, plot_tbl$ProbeID, sep="_")
  plot_tbl$Exposure <- factor(plot_tbl$Exposure, levels = exposure_levels)
  return(plot_tbl)
}

plot_DMR <- function(betas, pheno, gene_name, gene_name_label, DMR_start, DMR_end,
                     exposure_colname = "benz_dich",
                     exposure_levels = c(0,1), 
                     exposure_colors = c("#009E73", "#D55E00"),
                     exposure_labels = c("Benzene Controls", "Benzene Exposed"),
                     probe_ids = NULL){
  
  plottbl <- make_DMR_tbl(betas=betas, pheno=pheno, gene_name=gene_name, 
                          DMR_start=DMR_start, DMR_end=DMR_end,
                          exposure_colname = exposure_colname,
                          exposure_levels = exposure_levels)
  
  if(!is.null(probe_ids)){
    plottbl <- dplyr::filter(plottbl, ProbeID %in% probe_ids)
  }
  
  p <- plottbl %>%
    dplyr::group_by(Exposure, Position_ProbeID) %>%
    ggplot(., aes(x = Position_ProbeID, y = Value, fill = Exposure)) +
    geom_boxplot(outlier.size = 1) +
    scale_x_discrete("Position_ProbeID") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = exposure_colors, name = "", labels = exposure_labels) +
    labs(x = "CpG Probe ID and Genomic Location",
         y = "Beta Value (% DNA Methylation)",
         title = substitute(paste(italic(gene_name_label), " CpG Probes within DMR"))) +
    theme(legend.position="bottom")
  return(p)
}
