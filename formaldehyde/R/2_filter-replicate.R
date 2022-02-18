library(data.table)
library(here)
library(ChAMP)


################### load RGset, mvals, and phenotype table ######################
# note, I am assuming that the bad samples that appeared in red in the 
# PotentialOutliers plot have already been removed. 

# load phenotype table
pheno <- read.csv(here("data/pheno_all.csv")) 
rownames(pheno) <- pheno$Target_ID

# load RGset, the one that removed the bad samples in PotentialOutliers
# assuming that the name of this object is RGset
e <- new.env()
load(here("data", "rgsetAll.Rdata"), envir = e)
RGset <- updateObject(rgsetAll)
validObject(RGset)
names <- sub("^", "X", colnames(RGset))
colnames(RGset) <- names
# filter rgset
idx <- names %in% rownames(pheno) 
RGset <- RGset[,idx]
all.equal(colnames(RGset), rownames(pheno)) # make sure this is TRUE

# load mvalues, the ones produced after all preprocessing (no bad probes or samples)
# assuming that the name of this object is mvals
load(here("data/combat.RData"))
all.equal(colnames(mvals_combat), rownames(pheno))  # make sure this is TRUE
all.equal(colnames(mvals_combat), colnames(RGset))  # make sure this is TRUE

######################## identify duplicated samples 3##########################
#cols_with_duplicates <- ("age","sex","alcol", "BMI", "infection", "fa")
# duplicate <- which(duplicated(pheno[,c("age","sex","alcol", "BMI", 
#                                        "infection", "fa"), with=F]))
duplicate <- which(duplicated(pheno[,c("age","sex","alcol","BMI","infection","fa")]))
other_duplicate <- which(duplicated(pheno[,c("age","sex","alcol","BMI","infection","fa")], 
                         fromLast = T))
duplicates_tbl <- data.frame(duplicate, other_duplicate)

################## look at probe quality for the duplicates ####################
detP <- minfi::detectionP(RGset)
# we only care about the probes that remained after filtering
detP_relevant <- detP[which(rownames(detP) %in% rownames(mvals_combat)), ]
nrow(detP_relevant) == nrow(mvals_combat) #  make sure this is TRUE

detP_sum <- apply(detP_relevant, 2, sum)
detP_tbl <- data.frame(summary_duplicate = detP_sum[duplicate], 
                       summary_other_duplicate = detP_sum[other_duplicate])
detP_tbl



duplicates_tbl$min_detP <- apply(detP_tbl, 1, which.min)
duplicates_tbl                       



good_duplicates <- ifelse(duplicates_tbl$min_detP == 1, 
                          duplicates_tbl$duplicate, 
                          duplicates_tbl$other_duplicate)     
                          
worse_duplicates <- ifelse(duplicates_tbl$min_detP == 2, 
                           duplicates_tbl$duplicate, 
                           duplicates_tbl$other_duplicate)   

###### subset phenotype table and mvals to include good duplicate ##############
pheno_good <- pheno[-worse_duplicates, ]
design_good <- model.matrix(..., data = pheno_good)
mvals_good <- mvals_filtered[, -worse]

# proceed with analysis (DMP,DVP,DMR,DMB) using pheno_good, design_good, mvals_good
   
write.csv(pheno_good, file = here("data/pheno_good.csv"))  

# row 50 good: choose replicate 1
