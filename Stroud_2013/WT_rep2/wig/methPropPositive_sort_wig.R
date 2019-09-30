#!/applications/R/R-3.5.0/bin/Rscript

# Usage:
# ./methPropPositive_sort_wig.R GSM980986_WT_rep2 

args <- commandArgs(trailingOnly = T)
libName <- args[1]

library(rtracklayer)

# Import wig files
CG <- import.wig(paste0(libName, "_CG_original.wig"))
CHG <- import.wig(paste0(libName, "_CHG_original.wig"))
CHH <- import.wig(paste0(libName, "_CHH_original.wig"))

# Convert 'negative' DNA methylation proportion values into positive
CG$score <- abs(CG$score)
CHG$score <- abs(CHG$score)
CHH$score <- abs(CHH$score)

# Ensure seqlevels are sorted
CG <- sortSeqlevels(CG)
CHG <- sortSeqlevels(CHG)
CHH <- sortSeqlevels(CHH)

# Sort each UCSCData object
CG <- sort(CG)
CHG <- sort(CHG)
CHH <- sort(CHH)

# Export sorted wig files
export.wig(object = CG,
           con = paste0(libName, "_CG_methPropPositive_sorted.wig"))
export.wig(object = CHG,
           con = paste0(libName, "_CHG_methPropPositive_sorted.wig"))
export.wig(object = CHH,
           con = paste0(libName, "_CHH_methPropPositive_sorted.wig"))

# Convert into TDF format for loading into IGV
system(paste0("igvtools toTDF ",
              libName, "_CG_methPropPositive_sorted.wig ",
              libName, "_CG_methPropPositive_sorted.wig.tdf /projects/ajt200/TAIR10/tair10.chrom.sizes"))
system(paste0("igvtools toTDF ",
              libName, "_CHG_methPropPositive_sorted.wig ",
              libName, "_CHG_methPropPositive_sorted.wig.tdf /projects/ajt200/TAIR10/tair10.chrom.sizes"))
system(paste0("igvtools toTDF ",
              libName, "_CHH_methPropPositive_sorted.wig ",
              libName, "_CHH_methPropPositive_sorted.wig.tdf /projects/ajt200/TAIR10/tair10.chrom.sizes"))
