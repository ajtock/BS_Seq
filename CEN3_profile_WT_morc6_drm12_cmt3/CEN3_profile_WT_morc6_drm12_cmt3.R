######################################################################################################
# Plot DNA methylation in CG, CHG and CHH contexts in wild-type Columbia (WT) (Stroud et al., 2013), # 
# and morc6 (Moissiard et al., 2014), drm1/2 (Stroud et al., 2013)                                   #
######################################################################################################

library(dplyr)
library(segmentSeq)
library(parallel)

stroud.bed.dir <- "/home/meiosis/ajt200/BS_Seq/Stroud_2013/wig/bed/"
moissiard.bed.dir <- "/home/meiosis/ajt200/BS_Seq/Moissiard_2014/wig/bed/"
out.dir <- "/home/meiosis/ajt200/BS_Seq/CEN3_profile_WT_morc6_drm12_cmt3/"
plot.dir <- "/home/meiosis/ajt200/BS_Seq/CEN3_profile_WT_morc6_drm12_cmt3/plots/"

WT_rep2_CG <- read.table(file = paste0(stroud.bed.dir, "GSM980986_WT_rep2_CG.wig.bed.gr.tab.bed"))
WT_rep2_CHG <- read.table(file = paste0(stroud.bed.dir, "GSM980986_WT_rep2_CHG.wig.bed.gr.tab.bed"))
WT_rep2_CHH <- read.table(file = paste0(stroud.bed.dir, "GSM980986_WT_rep2_CHH.wig.bed.gr.tab.bed"))

drm12_CG <- read.table(file = paste0(stroud.bed.dir, "GSM981015_drm12_CG.wig.bed.gr.tab.bed"))
drm12_CHG <- read.table(file = paste0(stroud.bed.dir, "GSM981015_drm12_CHG.wig.bed.gr.tab.bed"))
drm12_CHH <- read.table(file = paste0(stroud.bed.dir, "GSM981015_drm12_CHH.wig.bed.gr.tab.bed"))

morc6_CG <- read.table(file = paste0(moissiard.bed.dir, "morc6_CG.wig.bed.gr.tab.bed"))
morc6_CHG <- read.table(file = paste0(moissiard.bed.dir, "morc6_CHG.wig.bed.gr.tab.bed"))
morc6_CHH <- read.table(file = paste0(moissiard.bed.dir, "morc6_CHH.wig.bed.gr.tab.bed"))

cmt3_CG <- read.table(file = paste0(stroud.bed.dir, "GSM981003_cmt3_CG.wig.bed.gr.tab.bed"))
cmt3_CHG <- read.table(file = paste0(stroud.bed.dir, "GSM981003_cmt3_CHG.wig.bed.gr.tab.bed"))
cmt3_CHH <- read.table(file = paste0(stroud.bed.dir, "GSM981003_cmt3_CHH.wig.bed.gr.tab.bed"))

CEN3.start <- c(11115724)
CEN3.end <- c(16520560)

WT_rep2_CG_CEN3 <- WT_rep2_CG %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)
WT_rep2_CHG_CEN3 <- WT_rep2_CHG %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)
WT_rep2_CHH_CEN3 <- WT_rep2_CHH %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)

drm12_CG_CEN3 <- drm12_CG %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)
drm12_CHG_CEN3 <- drm12_CHG %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)
drm12_CHH_CEN3 <- drm12_CHH %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)

morc6_CG_CEN3 <- morc6_CG %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)
morc6_CHG_CEN3 <- morc6_CHG %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)
morc6_CHH_CEN3 <- morc6_CHH %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)

cmt3_CG_CEN3 <- cmt3_CG %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)
cmt3_CHG_CEN3 <- cmt3_CHG %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)
cmt3_CHH_CEN3 <- cmt3_CHH %>%
  filter(V1 == "chr3") %>%
  filter(V2 >= CEN3.start & V2 <= CEN3.end)

CEN3.list <- list(WT_rep2_CG_CEN3, drm12_CG_CEN3, morc6_CG_CEN3, cmt3_CG_CEN3,
                  WT_rep2_CHG_CEN3, drm12_CHG_CEN3, morc6_CHG_CEN3, cmt3_CHG_CEN3,
                  WT_rep2_CHH_CEN3, drm12_CHH_CEN3, morc6_CHH_CEN3, cmt3_CHH_CEN3)
CEN3.names <- c("WT_rep2_CG_CEN3", "drm12_CG_CEN3", "morc6_CG_CEN3", "cmt3_CG_CEN3",
                "WT_rep2_CHG_CEN3", "drm12_CHG_CEN3", "morc6_CHG_CEN3", "cmt3_CHG_CEN3",
                "WT_rep2_CHH_CEN3", "drm12_CHH_CEN3", "morc6_CHH_CEN3", "cmt3_CHH_CEN3")

irl <- lapply(seq_along(CEN3.list), function(x) {
	IRanges(start = CEN3.list[[x]]$V2, width = 1)
})

grl <- lapply(seq_along(irl), function(x) {
	GRanges(seqnames = "chr3", strand = "+", ranges = irl[[x]])
})

windows <- seq(CEN3.start, CEN3.end, by = 20000)
windows <- c(windows, CEN3.end)
windows.iranges <- IRanges(start = windows, width = 20000)
windows.granges <- GRanges(seqnames = "chr3", strand = "+", ranges = windows.iranges)

overlaps <- mclapply(seq_along(grl), function(x) {
	getOverlaps(windows.granges, grl[[x]], whichOverlaps = TRUE)
}, mc.cores = 12)

win.vals <- mclapply(seq_along(overlaps), function(x) {
	sapply(overlaps[[x]], function(y) mean(as.numeric(CEN3.list[[x]]$V4[y])))
}, mc.cores = 12)
j = 50
test <- seq(1, 200, by = 1)
ma <- rep(1, test[j])/test[j]
filt <- mclapply(seq_along(win.vals), function(x) {
	stats::filter(win.vals[[x]], ma)
}, mc.cores = 12)
filtn <- mclapply(seq_along(filt), function(x) {
	filt[[x]][!is.na(filt[[x]])]
}, mc.cores = 12)

win.dat <- mclapply(seq_along(win.vals), function(x) {
	cbind(windows, win.vals[[x]])
}, mc.cores = 12)
lapply(seq_along(win.dat), function(x) {
	write.table(win.dat[[x]], file = paste0(out.dir, CEN3.names[[x]], ".20kb.txt"))
})

ymin.CG <- min(filtn[[1]], filtn[[2]], filtn[[3]], filtn[[4]])
ymax.CG <- max(filtn[[1]], filtn[[2]], filtn[[3]], filtn[[4]])
ymin.CHG <- min(filtn[[5]], filtn[[6]], filtn[[7]], filtn[[8]])
ymax.CHG <- max(filtn[[5]], filtn[[6]], filtn[[7]], filtn[[8]])
ymin.CHH <- min(filtn[[9]], filtn[[10]], filtn[[11]], filtn[[12]])
ymax.CHH <- max(filtn[[9]], filtn[[10]], filtn[[11]], filtn[[12]])

yminl <- list(ymin.CG, ymin.CHG, ymin.CHH)
ymaxl <- list(ymax.CG, ymax.CHG, ymax.CHH)

plotCEN3 <- function(filtData1, filtData2, filtData3, filtData4, ymin, ymax, context) {
	plot(windows, filtData1, type = "l", col = "black",
             ylim = c(ymin, ymax),
             main = "",
             xlab = "Chromosome 3 coordinates (bp)",
             ylab = paste(context, " methylation proportion", sep = ""))
	title(expression(bolditalic("CEN3") ~ bold(" (20-kb windows)")), line = 1)
	lines(windows, filtData2, col = "red")
	lines(windows, filtData3, col = "blue")
	lines(windows, filtData4, col = "purple")
	legend("topright", legend = c("wt", expression(italic("cmt3")), expression(italic("drm1/2")), expression(italic("morc6"))),
               col = c("black", "purple", "red", "blue"),
               text.col = c("black", "purple", "red", "blue"),
               ncol = 1, cex = 0.6, lwd = 1.2, bty = "n")
}

pdf(paste0(plot.dir, "CEN3_profile_wt_cmt3_drm12_morc6_20kb.pdf"))
par(mfrow = c(1, 3))
par(mar = c(20, 4.1, 20, 0.5))
par(mgp = c(2, 0.75, 0))
plotCEN3(filt[[1]], filt[[2]], filt[[3]], filt[[4]], yminl[[1]], ymaxl[[1]], "CG")
plotCEN3(filt[[5]], filt[[6]], filt[[7]], filt[[8]], yminl[[2]], ymaxl[[2]], "CHG")
plotCEN3(filt[[9]], filt[[10]], filt[[11]], filt[[12]], yminl[[3]], ymaxl[[3]], "CHH")
dev.off()


sessionInfo()
