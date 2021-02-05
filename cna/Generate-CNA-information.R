## Author: Luuk Harbers
## Date: 2020-07-16
## Script for generation of CNA information (overlaps and lengths)

## Load/install packages
packages = c("data.table", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Variables
sample = "BICRO205+206_MS18+NZ24"
binsize = "100kb"

# list files
files = list.files("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/gatk-cnv/BICRO205+206_MS18+NZ24/1e5/model_tumor/",
                   pattern = "called.seg", full.names = T)
annotation = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/barcodes/BICRO205_206/ffpe_barcodes.csv")
genelist = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/genes/genelist_nik-zainal-etal.tsv")

# Extract barcodes
barcodes = data.table(BARCODE = gsub(".*\\/|.dedup.*", "", files))
samples = merge(barcodes, annotation, all.x = T)

# Extract CNA regions
cnas = pblapply(1:length(files), function(i){
  dt = fread(files[i], skip = "CONTIG")
  dt[, sample := samples[i]$ID]
  dt[(CALL == "-" | CALL == "+") & !grepl("GL|MT|Y", CONTIG), c(1:3, 6:7)]
})
total = rbindlist(cnas)

# Set keys and find overlaps
setnames(genelist, c("CONTIG", "START", "END", "GENE"))
genelist[, CONTIG := factor(CONTIG)]
setkey(total, CONTIG, START, END)
setkey(genelist, CONTIG, START, END)

overlaps = foverlaps(genelist, total)
overlaps = overlaps[, c(1, 6:7, 4:5, 8)]
setnames(overlaps, c("CONTIG", "START", "END", "CALL", "sample", "gene"))

# Write results
write.table(total, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/cna-information/", 
                          sample, "_", binsize, "_CNAinformation.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
write.table(overlaps, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/cna-information/", 
                             sample, "_", binsize, "_Gene-overlaps_genelist.tsv"), quote = F, row.names = F, col.names = T, sep = "\t")
