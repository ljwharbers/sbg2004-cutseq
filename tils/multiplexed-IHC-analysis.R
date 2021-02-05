## Author: Luuk Harbers
## Date: 2020-08-19
## Script for analysis of multiplexed IHC data

## Load/install packages
packages = c("data.table", "ComplexHeatmap")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load in data
total = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/total_ratios.tsv")
tils = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/TILs.tsv")
foxp3 = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/foxp3_ratios.tsv")
pd1 = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/pd1_ratios.tsv")
pdl1 = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/pdl1_ratios.tsv")
pdl1_ck = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/pdl1_ck.tsv")

# # Add missing samples in pdl1_ck
# pdl1_ck = rbind(pdl1_ck, data.table(pid = c(203, 638, 705), dens = NA), use.names = F)
# setorder(pdl1_ck, PATIENT.ID)

# Remove pid 203 638 and 705
total = total[!grepl("203|638|705", patient_id)]
tils = tils[!grepl("203|638|705", patient_id)]
foxp3 = foxp3[!grepl("203|638|705", patient_id)]
pd1 = pd1[!grepl("203|638|705", patient_id)]
pdl1 = pdl1[!grepl("203|638|705", patient_id)]

# Plot total ratios
total_m = melt.data.table(total, id.vars = c("patient_id", "location"))
total_m[, location := factor(location, levels = c("stroma", "tumor"))]

# Prepare values for log transformation
total_m[, value_log := value + 1]

plt_all_box = ggplot(total_m, aes(x = variable, y = value_log, fill = location)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(y = "Cell density +1",
       x = "")

total_median = total_m[!is.na(value), median(value), by = .(location, variable)]
total_mean = total_m[!is.na(value), mean(value), by = .(location, variable)]

save_and_plot(plt_all_box, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ihc-plots/density_singles-boxplot",
              height = 7, width = 7)

plt_all_heat = ggplot(total_mean, aes(x = variable, y = location, fill = V1)) +
  geom_tile() +
  coord_flip() +
  labs(y = "",
       x = "", 
       fill = "Mean cell density") +
  scale_fill_gradient2(high = "red")

save_and_plot(plt_all_heat, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ihc-plots/density_singles-heatmap",
              height = 7, width = 7)

# Combine all co-expression
coexp = Reduce(merge, list(foxp3, pd1, pdl1))
coexp = coexp[, grepl("PD|patient|location", colnames(coexp)), with = F]
coexp_m = melt.data.table(coexp, id.vars = c("patient_id", "location"))

# Prepare values for log transformation
coexp_m[, value_log := value + 1]

plt_coexp_box = ggplot(coexp_m, aes(x = variable, y = value_log, fill = location)) +
  geom_boxplot() +
  scale_y_log10() +
  #coord_flip() +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(y = "Cell density +1",
       x = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

save_and_plot(plt_coexp_box, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ihc-plots/density_coexp-boxplot",
              height = 7, width = 7)

coexp_median = coexp_m[!is.na(value), median(value), by = .(location, variable)]
coexp_mean = coexp_m[!is.na(value), mean(value), by = .(location, variable)]

plt_coexp_heat = ggplot(coexp_mean, aes(x = variable, y = location, fill = V1)) +
  geom_tile() +
  coord_flip() +
  labs(y = "",
       x = "", 
       fill = "Mean cell density") +
  scale_fill_gradient2(high = "red")

save_and_plot(plt_coexp_heat, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ihc-plots/density_coexp-heatmap",
              height = 7, width = 7)
# Plot sTILs
stils = merge(tils, total[location == "stroma", total[location == "stroma", list(total_immune = sum(CD4, CD8)), by = patient_id]])
stils_scatter = ggplot(stils, aes(x = TIL_percentage, y = total_immune)) +
  geom_point() +
  labs(y = "Total stromal immune cells by IF (density)", x = "s-TILs by H&E (%)")

# Cluster analysis
total_matrix_stroma = cbind(total[location == "stroma", 1:5], 
                            pd1[location == "stroma", 2:5,],
                            pdl1[location == "stroma", 2:5])
total_matrix_tumor = cbind(total[location == "tumor", 1:5], 
                           pd1[location == "tumor", 2:5,],
                           pdl1[location == "tumor", 2:5],
                           pdl1_ck[, 2])

# Set as DF to make rownames into columns
setDF(total_matrix_stroma)
setDF(total_matrix_tumor)

rownames(total_matrix_stroma) = as.character(total_matrix_stroma$patient_id)
rownames(total_matrix_tumor) = as.character(total_matrix_tumor$patient_id)

total_matrix_stroma$patient_id = NULL
total_matrix_tumor$patient_id = NULL

total_matrix_tumor = total_matrix_tumor[rowSums(is.na(total_matrix_tumor)) != ncol(total_matrix_tumor),]

# Plot and cluster
total_heatmap_stroma = Heatmap(t(as.matrix(total_matrix_stroma)))
total_heatmap_tumor = Heatmap(t(as.matrix(total_matrix_tumor)))

save_and_plot(total_heatmap_stroma, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ihc-plots/clustered-heatmap_stroma",
              height = 7, width = 14)
save_and_plot(total_heatmap_tumor, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ihc-plots/clustered-heatmap_tumor",
              height = 7, width = 14)

# sTILs correlation versus CD4+CD8 coexp
tilcd4cd8 = merge(tils, total[location == "stroma", c(1, 2, 4)])
tilcd4cd8[, CD4CD8 := CD4+CD8]
tilcd4cd8[, state := ifelse(TIL_percentage >= 50, "lymphocytic predominant", "non-predominant")]

tilcd4cd8_m = melt(tilcd4cd8, id.vars = c("patient_id", "state", "TIL_percentage"))

tilPlot = ggplot(tilcd4cd8_m[complete.cases(tilcd4cd8_m)], aes(x = variable, y = value, color = state)) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge()) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) + 
  labs(y = "Density", x = "", color = "") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

save_and_plot(tilPlot, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/ihc-plots/til_cd4cd8_boxplots",
              height = 7, width = 7)
