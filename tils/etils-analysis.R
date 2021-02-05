## Author: Luuk Harbers
## Date: 2020-12-24
## Script for analysing eTILs

## Load/install packages
packages = c("data.table")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

wsi_etils = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/etils_wsi.tsv", select = c(1, 30:32, 34))
tma_etils = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/etils_tma.tsv", select = c(1:4, 6))
vis_tils = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/TILs.tsv")
ihc = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/total_ratios.tsv")

# Set new names to merge
setnames(wsi_etils, c("pid", "wsi_etils", "wsi_ettils", "wsi_estils", "wsi_eastils"))
setnames(tma_etils, c("pid", "tma_etils", "tma_ettils", "tma_estils", "tma_eastils"))
setnames(vis_tils, c("pid", "til"))

total = Reduce(function(x, y) merge(x, y, by = "pid", all = T), list(wsi_etils, tma_etils, vis_tils))
total = total[rowSums(is.na(total)) < 9]

# Plot wsi versus visual tils
wsi_subset = total[, .(pid, til, wsi_etils, wsi_eastils)]
wsi_subset = wsi_subset[complete.cases(wsi_subset)]

wsi_subset[, lbpc := ifelse(til >= 50, "LBPC" , "non-LBPC")]
wsi_subset_m = melt(wsi_subset, id.vars = c("pid", "lbpc", "til"))

wsi_subset_m[, lbpc := factor(lbpc, levels = c("non-LBPC", "LBPC"))]
plt1 = ggplot(wsi_subset_m[variable == "wsi_etils"], aes(x = lbpc, y = value, color = lbpc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .2) +
  stat_compare_means(paired = F, comparisons = list(c("non-LBPC", "LBPC"))) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(y = "WSI_eTILs (%)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/tils/wsi_etils-vis_tils-boxplot",
              width = 7, height = 7)

plt2 = ggplot(wsi_subset_m[variable == "wsi_eastils"], aes(x = lbpc, y = value, color = lbpc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .2) +
  stat_compare_means(paired = F, comparisons = list(c("non-LBPC", "LBPC"))) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(y = "WSI_easTILs (%)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

save_and_plot(plt2, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/tils/wsi_eastils-vis_tils-boxplot",
              width = 7, height = 7)

# Plot tma versus visual tils
tma_subset = total[, .(pid, til, tma_etils, tma_eastils)]
tma_subset = tma_subset[complete.cases(tma_subset)]

tma_subset[, lbpc := ifelse(til >= 50, "LBPC" , "non-LBPC")]
tma_subset_m = melt(tma_subset, id.vars = c("pid", "lbpc", "til"))

tma_subset_m[, lbpc := factor(lbpc, levels = c("non-LBPC", "LBPC"))]
plt3 = ggplot(tma_subset_m[variable == "tma_etils"], aes(x = lbpc, y = value, color = lbpc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .2) +
  stat_compare_means(paired = F, comparisons = list(c("non-LBPC", "LBPC"))) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(y = "TMA_eTILs (%)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

save_and_plot(plt3, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/tils/tma_etils-vis_tils-boxplot",
              width = 7, height = 7)

plt4 = ggplot(tma_subset_m[variable == "tma_eastils"], aes(x = lbpc, y = value, color = lbpc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .2) +
  stat_compare_means(paired = F, comparisons = list(c("non-LBPC", "LBPC"))) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(y = "TMA_easTILs (%)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

save_and_plot(plt4, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/tils/tma_eastils-vis_tils-boxplot",
              width = 7, height = 7)

# Plot IHC density versus visual tils
setnames(ihc, "patient_id", "pid")
ihc_til = merge(ihc[, .(pid, CD4, CD8, location)], vis_tils)
ihc_til[, `CD4 + CD8` := CD4 + CD8]
ihc_til[, lbpc := ifelse(til >= 50, "LBPC" , "non-LBPC")]
ihc_til = ihc_til[complete.cases(ihc_til)]

ihc_til_m = melt(ihc_til, id.vars = c("pid", "til", "lbpc", "location"))
ihc_til_m[, lbpc := factor(lbpc, levels = c("non-LBPC", "LBPC"))]

plt5 = ggplot(ihc_til_m[location == "stroma"], aes(x = lbpc, y = value, color = lbpc)) +
  facet_wrap(~variable) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .2) +
  stat_compare_means(paired = F, comparisons = list(c("non-LBPC", "LBPC"))) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(y = "IHC cell density (cells/mm2)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

save_and_plot(plt5, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/tils/ihc_density-lbpc-boxplot-stroma",
              width = 7, height = 11)

plt6 = ggplot(ihc_til_m[location == "tumor"], aes(x = lbpc, y = value, color = lbpc)) +
  facet_wrap(~variable) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .2) +
  stat_compare_means(paired = F, comparisons = list(c("non-LBPC", "LBPC"))) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(y = "IHC cell density (cells/mm2)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")

save_and_plot(plt6, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/tils/ihc_density-lbpc-boxplot-tumor",
              width = 7, height = 11)
