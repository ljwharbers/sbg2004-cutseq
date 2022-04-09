# get amplified/deleted genes

# LIBRARIES and functions
require(data.table)

# SET PARAMETERS

sample = "BICRO205_206_HQ"
binsize = "100kb"

#load rds file
total <- readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Scripts/rds-files/BICRO205_206_HQ_100kb_CNAinformation.rds")

#load gene list
genes <- fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/genes/gene-list.tsv")

# RUN

#make character contig names
genes[, CONTIG := as.character(CONTIG)]

#setkeys
setkey(total, CONTIG, START, END)
setkey(genes, CONTIG, START, END)

#get overlaps
combined = foverlaps(genes, total)
combined = combined[, c(1, 6:7, 4:5, 8)]

#setnames
setnames(combined, c("CONTIG", "START", "END", "CALL", "patient_id", "gene"))

#save rds
saveRDS(combined, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Scripts/rds-files/", 
                         sample, "-", binsize, "_CNA-geneoverlaps.rds"))

#save table of with all samples and gene status
tbl = dcast(combined, patient_id ~ gene, value.var = "CALL")
tbl[is.na(tbl)] = "0"

#write table
write.table(tbl, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/alteration-tables/", 
                        sample, "-", binsize, "altered-genes.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
