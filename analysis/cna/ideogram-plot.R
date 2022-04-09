# Load example data
require(data.table)
require(pbapply)
require(GenomicRanges)
require(karyoploteR)
require(naturalsort)

total = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO205+206_MS18+NZ24_1e5_gatk-cnv.rds")

#sample
sample = "BICRO205+206_MS18+NZ24"
binsize = "100kb"

#set cores  
cl = 30

total = pblapply(total, function(x) {
  # set proper names for granges object
  setnames(x, "CONTIG", "seqnames")
  
  #set 'chr' and add column for HER2
  x[, seqnames := paste0("chr", seqnames)]
  x[FEATURE == "17:37800001-37900000", HER2 := MEAN_LOG2_COPY_RATIO]
  makeGRangesFromDataFrame(x, keep.extra.columns = T)
}, cl = cl)

if(!dir.exists(paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ideograms/", sample))) {
  dir.create(paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ideograms/", sample))
}

pblapply(1:length(total), function(x) {
  
  #save path
  cairo_ps(file = paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ideograms/", sample, "/",
                     names(total[x]), "_", binsize, "_cairo_ps.eps"), height = 6, width = 12)
  
  #plotting
  kp <- plotKaryotype(chromosomes = "chr17", cex = 2)
  kpPoints(kp, data = total[[x]], y = total[[x]]$LOG2_COPY_RATIO, cex = 0.6, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "grey50")
  kpPoints(kp, data = total[[x]], y = total[[x]]$MEAN_LOG2_COPY_RATIO, cex = 0.8, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "black")
  if(total[[x]]$CALL[total[[x]]$FEATURE == "17:37800001-37900000"] == "+") {
    kpPoints(kp, data = total[[x]], y = total[[x]]$HER2, cex = 1.5, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "red")
  }
  if(total[[x]]$CALL[total[[x]]$FEATURE == "17:37800001-37900000"] == "-") {
    kpPoints(kp, data = total[[x]], y = total[[x]]$HER2, cex = 1.5, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "blue")
  }
  if(total[[x]]$CALL[total[[x]]$FEATURE == "17:37800001-37900000"] == "0") {
    kpPoints(kp, data = total[[x]], y = total[[x]]$HER2, cex = 1.5, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "black")
  }
  
  kpAxis(kp, ymin=-2, ymax=4, col="grey50", cex=2, numticks = 4)
  #kpAbline(kp, h = log2(2.5/2), lty = 2, r0 = 0, r1 = 1, ymin = -2, ymax = 4)
  #kpAbline(kp, h = log2(1.5/2), lty = 2, r0 = 0, r1 = 1, ymin = -2, ymax = 4)
  kpAddMainTitle(kp, main = names(total[x]), cex = 4)
  dev.off()

}, cl = cl)


