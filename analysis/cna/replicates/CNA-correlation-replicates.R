## Author: Luuk Harbers
## Date: 2020-07-15
## Script for analysis of percentage of amplification/deletion in the replicates

## Load/install packages
packages = c("data.table", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

## Variables
threads = 20

## Load replicates
rep1 = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO205+206_MS18_1e5_gatk-cnv.rds")
rep2 = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO206_NZ24_1e5_gatk-cnv.rds")

## Get samples in both replicates
samples = intersect(names(rep1), names(rep2))
total_rows = nrow(rep1[[1]])

## Apply through files and extract percentage amplified/deleted
cna_perc = pblapply(samples, function(sample){
  rep1 = (nrow(rep1[[sample]][CALL == "+" | CALL == "-"]) / total_rows * 100)
  rep2 = (nrow(rep2[[sample]][CALL == "+" | CALL == "-"]) / total_rows * 100)
  data.table(sample = sample, rep1 = rep1, rep2 = rep2)
}, cl = threads)

total = rbindlist(cna_perc)

## Plot
plt = ggplot(total, aes(x = rep1, y = rep2)) +
  geom_point() +
  #geom_text(aes(label = sample))
  geom_smooth(method = "lm", se = F, color = "red", linetype = 2) +
  stat_cor() +
  labs(x = "Replicate 1",
       y = "Replicate 2")

save_and_plot(plt, 
              "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/replicates/CNA-correlation-replicates",
              width = 7, height = 7)
