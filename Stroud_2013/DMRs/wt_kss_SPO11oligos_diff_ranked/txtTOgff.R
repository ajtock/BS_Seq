#!/applications/R/R-3.5.0/bin/Rscript

# Convert txt file containing peaks to gff file

# Usage:
# ./txtTOgff.R kss_hypoCHG_DMRs wt_kss_SPO11oligos_diff

args <- commandArgs(trailingOnly = T)
featureName <- args[1]
diffName <- args[2]

peaks <- read.table(paste0("./", featureName, "_",
                           diffName, "_ordered.txt"),
                    header = T)

peaks_gff <- cbind(peaks$chr,
                   rep("."),
                   rep(featureName),
                   peaks$start,
                   peaks$end,
                   peaks$wt_kss_SPO11oligos_diff,
                   rep("."),
                   rep("."),
                   rep("."))
colnames(peaks_gff) <- c("chr",
                         "source",
                         "feature",
                         "start",
                         "end",
                         "score",
                         "strand",
                         "frame",
                         "attribute")
write.table(peaks_gff,
            file = paste0("./", featureName, "_",
                          diffName, "_ordered.gff"),
            sep = "\t", quote = F, row.names = F, col.names = F)
