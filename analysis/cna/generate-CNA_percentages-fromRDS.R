#Create full info RDS file

# LIBRARIES
require(data.table)
require(naturalsort)
require(pbapply)

# PARAMETERS
#sample name
sample <- "BICRO205_206HQ"
binsize <- "100kb"

#read rds
output <- readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Scripts/rds-files/BICRO205_206_HQ_100kb.rds")

#set cores to use
cores = 20

# RUN

#calculate percentage of amp/del/aneuploidy
totalLength <- nrow(output[[1]])
percAmp <- unlist(pblapply(output, function(x) nrow(x[CALL == "+"]) / totalLength * 100, cl = cores))
percDel <- unlist(pblapply(output, function(x) nrow(x[CALL == "-"]) / totalLength * 100, cl = cores))
percAneu <- percAmp + percDel

#combine all data
total <- as.data.table(cbind(names(percAmp), percAmp, percDel, percAneu))
setnames(total, "V1", "Patient_ID")

#save data
write.table(total, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/CNA-percentage/CNA-percentage_", sample, "_", binsize, ".tsv"), sep = "\t", quote = F, col.names = T, row.names = F)

