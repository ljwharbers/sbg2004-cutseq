## Author: Luuk Harbers
## Date: 2020-07-16
## Script for stratification of CNA lengths

## Load/install packages
packages = c("data.table", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Variables
sample = "BICRO205+206_MS18+NZ24"
binsize = "100kb"

## Load data
cna = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/cna-information/BICRO205+206_MS18+NZ24_100kb_CNAinformation.tsv")
cna[, length_mb := (END - START) / 1e6]

cna[length_mb < 1, stratification := "focal"]
cna[length_mb >= 1 & length_mb < 10, stratification := "medium"]
cna[length_mb >= 10, stratification := "large"]

# Get counts per sample, stratification and type of cna
strat = cna[, .N, by = c("sample", "stratification", "CALL")]
strat = dcast(strat, sample ~ stratification + CALL, value.var = "N")

# Rename
setnames(strat, c("sample", "focal_amp", "focal_del", "large_amp", "large_del", "medium_amp", "medium_del"))

# Get total
strat[is.na(strat)] = 0L
strat[, total_amp := focal_amp + large_amp + medium_amp]
strat[, total_del := focal_del + large_del + medium_del]

# Get percentage
strat[, focal_amp_perc := focal_amp / total_amp * 100]
strat[, medium_amp_perc := medium_amp / total_amp * 100]
strat[, large_amp_perc := large_amp / total_amp * 100]
strat[, focal_del_perc := focal_del / total_del * 100]
strat[, medium_del_perc := medium_del / total_del * 100]
strat[, large_del_perc := large_del / total_del * 100]

# Split stratification and cna
freq_melt = melt(strat[, grepl("sample|perc", colnames(strat)), with = F], id.vars = "sample")
freq_melt[, variable := gsub("_perc", "", variable)]
freq_melt[, c("stratification", "cna") := tstrsplit(variable, "_")]
freq_melt[, stratification := factor(stratification, levels = c("focal", "medium", "large"))]

# Plot
plt = ggplot(freq_melt, aes(x = stratification, y = value, color = cna)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15)) +
  scale_color_brewer(palette = "Set1") +
  scale_x_discrete(labels = c("< 1Mb", "1-10Mb", "> 10Mb")) +
  labs(y = "Proportion of alteration type in sample (%)",
       x = "",
       color = "")

save_and_plot(plt, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/cna-length/", sample, "_", binsize, 
                          "_cna-length_stratification"), height = 7, width = 7)
