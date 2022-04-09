## Author: Luuk Harbers
## Date: 2020-07-22
## Script to create ideograms per chromosome to highlight genes in genelist

## Load/install packages
packages = c("data.table", "pbapply", "GenomicRanges", "karyoploteR", "naturalsort")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load in data
total = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO205+206_MS18+NZ24_1e5_gatk-cnv.rds")
genes = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/genes/genelist_nik-zainal-etal.tsv")

# Sample characteristics
sample = "BICRO205+206_MS18+NZ24"
binsize = "100kb"

# Set threads
nthreads = 32

# Subset genelist
genes = genes[Gene %in% c("ERBB2", "CCND1", "FGFR1", "ZNF217", "IGF1R", "SMAD4", "CCNE1", "MDM2", "FGFR2", "PTEN")]

# Set key and make chr character
genes[, chr := paste0("chr", chr)]
setkey(genes, chr, start, end)

# Add column names with gene name in bins
total = pblapply(total, function(dt) {

  # Set chr prefix
  dt[, CONTIG := paste0("chr", CONTIG)]
  setkey(dt, CONTIG, START, END)
  
  # Get overlaps with genes
  overlap = foverlaps(genes, dt)
  overlap = unique(overlap, by = "Gene")
  dt = merge(dt, overlap[, .(FEATURE, Gene)], by = "FEATURE", all.x = T)
  
  # Add value columns
  dt[!is.na(Gene), gene_value := MEAN_LOG2_COPY_RATIO]
  
  # Change NAs 
  dt[is.na(Gene), Gene := ""]

  
  # Make granges
  setnames(dt, "CONTIG", "seqnames")
  return(makeGRangesFromDataFrame(dt, keep.extra.columns = T))
}, cl = nthreads)

sapply(genes$Gene, function(gene) {
  dir.create(paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ideograms/", sample, "/", gene), showWarnings = F)
})

lapply(1:nrow(genes), function(row){
  
  pblapply(1:length(total), function(x) {
    
    cairo_ps(file = paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ideograms/", sample, "/", genes[row]$Gene, "/",
                           names(total[x]), "_", binsize, "_cairo_ps.eps"), height = 6, width = 12)
    # Plotting
    kp <- plotKaryotype(chromosomes = genes[row]$chr, cex = 2)
    kpPoints(kp, data = total[[x]], y = total[[x]]$LOG2_COPY_RATIO, cex = 0.6, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "grey50")
    kpPoints(kp, data = total[[x]], y = total[[x]]$MEAN_LOG2_COPY_RATIO, cex = 0.8, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "black")
    
    if(total[[x]][total[[x]]$Gene == genes[row]$Gene]$CALL == "+") {
      kpPoints(kp, data = total[[x]], y = total[[x]]$gene_value, cex = 1.5, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "red")
    }
    if(total[[x]][total[[x]]$Gene == genes[row]$Gene]$CALL == "-") {
      kpPoints(kp, data = total[[x]], y = total[[x]]$gene_value, cex = 1.5, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "blue")
    }
    if(total[[x]][total[[x]]$Gene == genes[row]$Gene]$CALL == "0") {
      kpPoints(kp, data = total[[x]], y = total[[x]]$gene_value, cex = 1.5, r0 = 0, r1 = 1, ymax = 4, ymin = -2, col = "black")
    }
    kpText(kp, ymin = -1, ymax = 4, col = "black", data = total[[x]], labels = total[[x]]$Gene, 
           y = total[[x]][total[[x]]$Gene == genes[row]$Gene]$LOG2_COPY_RATIO + 1)
    kpAxis(kp, ymin=-2, ymax=4, col="grey50", cex=2, numticks = 4)
    kpAddMainTitle(kp, main = names(total[x]), cex = 4)
  
    dev.off()
    
  }, cl = nthreads)
})
  


