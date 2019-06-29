#!/applications/R/R-3.3.2/bin/Rscript

# Calculate mean library-size-normalised wt and kss SPO11oligos
# coverage between the start and end of each DMR
# Order DMRs by decreasing SPO11oligos loss in kss

# Usage:
# ./DMRs_wt_kss_SPO11oligos_diff_ranked.R suvh456_hypoCHG_DMR kss_hypoCHG_DMRs wt_kss_SPO11oligos_diff

args <- commandArgs(trailingOnly = T)
featureNameIn <- args[1]
featureNameOut <- args[2]
diffName <- args[3]

library(segmentSeq)

# Chromosome definitions
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")

# Import DMRs as GRanges object
DMRs <- read.table(paste0("../", featureNameIn, "_vs3reps_min4filter_mg200.bed"),
                   header = F)
DMRsGR <- GRanges(seqnames = DMRs[,1],
                  ranges = IRanges(start = DMRs[,2], end = DMRs[,3]),
                  strand = "*")
seqlevels(DMRsGR) <- sub("chr", "Chr", seqlevels(DMRsGR))
print("***********DMRs***********")
print(DMRsGR)

wt_norm <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/WT/coverage/log2ChIPinput/log2wtSPO11oligoRPI1NakedDNA_norm_allchrs_coverage_coord_tab.bed",
                     colClasses = c(NA, NA, "NULL", NA))
kss_norm <- read.table("/projects/ajt200/BAM_masters/SPO11-oligo/suvh456/coverage/log2ChIPinput/log2suvh456SPO11oligoRPI34NakedDNA_norm_allchrs_coverage_coord_tab.bed",
                       colClasses = c(NA, NA, "NULL", NA))

allCovDMRs <- data.frame()
for(i in 1:5) {
  print(i)
  wt_chrCov <- wt_norm[wt_norm[,1] == chrs[i],]
  wt_chrCov <- wt_chrCov[,3]
  kss_chrCov <- kss_norm[kss_norm[,1] == chrs[i],]
  kss_chrCov <- kss_chrCov[,3]
  print(i)
  covCoords <- seq(1, length(wt_chrCov), by = 1)
  covGR <- GRanges(seqnames = chrs[i],
                   ranges = IRanges(start = covCoords,
                                    width = 1),
                   strand = "*")
  print(i)
  chrDMRsGR <- DMRsGR[seqnames(DMRsGR) == chrs[i]]
  DMRsOverlaps <- getOverlaps(chrDMRsGR,
                              covGR,
                              whichOverlaps = T)
  wt_norm_DMRCov <- sapply(DMRsOverlaps,
                           function(x) mean(wt_chrCov[x]))
  kss_norm_DMRCov <- sapply(DMRsOverlaps,
                            function(x) mean(kss_chrCov[x]))
  print(i)
  chrDMRs <- data.frame(chr = rep(chrs[i], length(chrDMRsGR)),
                        start = as.integer(start(chrDMRsGR)),
                        end = as.integer(end(chrDMRsGR)),
                        wt_SPO11oligos_mean = as.numeric(wt_norm_DMRCov),
                        kss_SPO11oligos_mean = as.numeric(kss_norm_DMRCov),
                        wt_kss_SPO11oligos_diff = as.numeric(wt_norm_DMRCov-kss_norm_DMRCov))
  allCovDMRs <- rbind(allCovDMRs, chrDMRs)
}
allCovDMRs_ordered <- allCovDMRs[order(as.numeric(allCovDMRs$wt_kss_SPO11oligos_diff),
                                       decreasing = T),]
write.table(allCovDMRs_ordered,
            file = paste0("./", featureNameOut, "_", diffName, "_ordered.txt"),
            sep = "\t", row.names = F, col.names = T, quote = F)
