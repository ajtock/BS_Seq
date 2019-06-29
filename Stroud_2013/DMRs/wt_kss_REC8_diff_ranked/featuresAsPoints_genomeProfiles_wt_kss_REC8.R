#!/applications/R/R-3.5.0/bin/Rscript

# Plot chromosome profiles of features presented as points along the genome (x-axis)
# and log2(ChIP/input) on the y-axis in wild type and mutant

# Usage:
# ./featuresAsPoints_genomeProfiles_wt_kss_REC8.R kss_hypoCHG_DMRs_wt_kss_REC8_diff_ordered.txt kss_hypoCHG_DMRs 'HypoCHG DMRs in kyp suvh5 suvh6'

args <- commandArgs(trailingOnly = T)
featureFile <- args[1]
featureName <- args[2]
featureNamePlot <- args[3]

featureTab <- read.table(paste0("./", featureFile),
                         header = T)
featureTab <- featureTab[order(featureTab$chr,
                               featureTab$start,
                               featureTab$end),]

# Genomic definitions
chrs <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
# Pericentromeric regions are as defined in Supplemental Table S26
# of Ziolkowski et al. (2017) Genes Dev. 31
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)

# Make cumulative genomes
sumchr <- cumsum(c(0, chrLens))
print(sumchr)
sumchr_tot <- sumchr[length(sumchr)]
print(sumchr_tot)

centromeres <- sapply(seq_along(centromeres), function(x) {
  centromeres[x] + sumchr[x]
})
print(centromeres)
pericenStart <- sapply(seq_along(pericenStart), function(x) {
  pericenStart[x] + sumchr[x]
})
print(pericenStart)
pericenEnd <- sapply(seq_along(pericenEnd), function(x) {
  pericenEnd[x] + sumchr[x]
})
print(pericenEnd)

plotDir <- "./plots/"
system(paste0("[ -d ", plotDir, " ] || mkdir ", plotDir))

# Calculate cumulative coordinates and midpoints
featureTabCumStart <- NULL
featureTabCumEnd <- NULL
featureTabCumMidpoint <- NULL
for(i in 1:length(chrs)) {
  featureTabChr <- featureTab[featureTab$chr == chrs[i],]
  featureTabChrCumStart <- featureTabChr$start+sumchr[i]
  featureTabChrCumEnd <- featureTabChr$end+sumchr[i]
  featureTabChrCumMidpoint <- round(featureTabChrCumStart+
                                    ((featureTabChrCumEnd-featureTabChrCumStart)/2))
  featureTabCumStart <- c(featureTabCumStart, featureTabChrCumStart)
  featureTabCumEnd <- c(featureTabCumEnd, featureTabChrCumEnd)
  featureTabCumMidpoint <- c(featureTabCumMidpoint, featureTabChrCumMidpoint)
}

featureTabCum <- data.frame(featureTab,
                            midpoint = as.integer(round(featureTab$start+
                                                        ((featureTab$end-featureTab$start)/2))),
                            cum_start = as.integer(featureTabCumStart),
                            cum_end = as.integer(featureTabCumEnd),
                            cum_midpoint = as.integer(featureTabCumMidpoint))

# Function to plot features as points along the genome (x-axis)
# and log2(ChIP/input) on the y-axis in wild type and mutant
dat1Vdat2GenomePointPlot <- function(xplot,
                                     dat1,
                                     dat2,
                                     legendNames,
                                     plotColours) {
  plot(x = xplot, y = dat1,
       type = "p", pch = 16, cex = 0.25, col = plotColours[1],
       ylim = c(min(dat1, dat2),
                max(dat1, dat2)),
       xlim = c(0, sumchr_tot),
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       main = featureNamePlot)
  points(x = xplot, y = dat2,
         type = "p", pch = 16, cex = 0.25, col = plotColours[2])
  axis(side = 1, lwd.tick = 2, cex.axis = 1.25,
       at = c(0, 2e7, 4e7, 6e7, 8e7, 1e8, 1.2e8),
       labels = c("0", "20", "40", "60", "80", "100", "120"))
  mtext(side = 1, line = 2.30, cex = 1.5,
        text = "Coordinates (Mb)")
  axis(side = 2, lwd.tick = 2, cex.axis = 1.25)
  mtext(side = 2, line = 2.30, cex = 1.5,
        text = expression("Log"[2]*"(REC8 ChIP/input)"))
  abline(v = sumchr, lty = 1, lwd = 2)
  abline(v = centromeres, lty = 5, lwd = 2)
  box(lwd = 2)
  legend("topleft",
         legend = legendNames,
         col = "white",
         text.col = plotColours,
         text.font = c(1, 3),
         ncol = 1, cex = 0.7, lwd = 2, bty = "n")
} 

# Plot
pdf(paste0(plotDir, featureName, "_as_points_wt_kss_REC8_genomeProfiles_v120619.pdf"),
    height = 4.5, width = 12)
par(mfcol = c(1, 1))
par(mar = c(4.1, 5.1, 2.1, 6.1))
par(mgp = c(3, 1, 0))
dat1Vdat2GenomePointPlot(xplot = featureTabCum$cum_midpoint,
                         dat1 = featureTabCum$wt_REC8_mean,
                         dat2 = featureTabCum$kss_REC8_mean,
                         legendNames = c("Wild type", "kyp suvh5 suvh6"),
                         plotColours = c("red", "navy"))
dev.off()
