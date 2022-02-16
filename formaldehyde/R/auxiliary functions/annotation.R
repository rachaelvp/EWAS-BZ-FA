library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

data(Locations)
data(Other)
annot <- data.frame(Other)
loc <- data.frame(Locations)
probes <- merge(loc, annot, by = "row.names")
rownames(probes) <- probes$Row.names
probes <- probes[,c(2:4,9:18)]

data(SNPs.Illumina)
data(Islands.UCSC)
snps <- data.frame(SNPs.Illumina)
islands <- data.frame(Islands.UCSC)
probes2 <- merge(snps, islands, by = "row.names")
rownames(probes2) <- probes2$Row.names
probes2 <- probes2[,c(2:5)]


annotation <- merge(probes, probes2, by = "row.names")
colnames(annotation)[1] <- "Probe_ID"

save(annotation, file = here("SuperFund/data","450k_annot.RData"))


# https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf
