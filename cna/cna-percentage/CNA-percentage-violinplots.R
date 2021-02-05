## Author: Luuk Harbers
## Date: 2020-07-16
## Script for analysis of percentage of AMP/DEL for 1 cohort (including subtypes)

## Load/install packages
packages = c("data.table", "ggplot2", "pbapply")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Variables
sample = "BICRO205+206_MS18+NZ24"
binsize = "100kb"
threads = 20

# Load in data
data = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO205+206_MS18+NZ24_1e5_gatk-cnv.rds")
#annotation = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/")

# Extract AMP/DEL percentage
total_rows = nrow(data[[1]])
altered = pblapply(1:length(data), function(i){
  data.table(sample = names(data[i]),
             Amplified = (sum(data[[i]]$CALL == "+") / total_rows) * 100,
             Deleted = (sum(data[[i]]$CALL == "-") / total_rows) * 100)
})
total = rbindlist(altered)
total_melt = melt(total, id.vars = "sample")

# Plot
plt = ggplot(total_melt, aes(x = variable, y = value, color = variable)) +
  #geom_violin(aes(fill = variable, color = variable)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(y = "Percentage of genome altered (%)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

save_and_plot(plt, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/cna_burden/", sample, "_",
                          binsize, "_CNAburden_boxplot"), height = 7, width = 6)

