## Author: Luuk Harbers
## Date: 2020-12-24
## Script for analysing eTILs (correlation matrices)

## Load/install packages
packages = c("data.table", "Hmisc")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

wsi_etils = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/etils_wsi.tsv", select = c(1, 30:32, 34))
tma_etils = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/etils_tma.tsv", select = c(1:4, 6))
vis_tils = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/TILs.tsv")
ihc = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/total_ratios.tsv")
foxp3 = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/foxp3_ratios.tsv")
pd1 = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/pd1_ratios.tsv")
pdl1 = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/pdl1_ratios.tsv")
pdl1_ck = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/multiplexed-IHC/pdl1_ck.tsv")

# Set new names to merge
setnames(wsi_etils, c("pid", "eTILs (WSI)", "etTILs (WSI)", "esTILs (WSI)", "easTILs (WSI)"))
setnames(tma_etils, c("pid", "eTILs (TMA)", "etTILs (TMA)", "esTILs (TMA)", "easTILs (TMA)"))
setnames(vis_tils, c("pid", "TIL"))

total = Reduce(function(x, y) merge(x, y, by = "pid", all = T), list(wsi_etils, tma_etils, vis_tils))
total = total[rowSums(is.na(total)) < 9]

# Merge coexpressions
coexp_stromal = Reduce(function(x, y) merge(x, y, by = "patient_id", all = T), 
                       list(ihc[location == "stroma", 1:5],
                            foxp3[location == "stroma", 1:5],
                            pd1[location == "stroma", 1:5],
                            pdl1[location == "stroma", 1:5]))
setnames(coexp_stromal, c("pid", paste0(names(coexp_stromal)[2:ncol(coexp_stromal)], " (Stroma)")))

coexp_tumor = Reduce(function(x, y) merge(x, y, by = "patient_id", all = T), 
                       list(ihc[location == "tumor", 1:5],
                            foxp3[location == "tumor", 1:5],
                            pd1[location == "tumor", 1:5],
                            pdl1[location == "tumor", 1:5]))
setnames(coexp_tumor, c("pid", paste0(names(coexp_tumor)[2:ncol(coexp_tumor)], " (Tumor)")))

# Merge everything
total_mat = Reduce(function(x, y) merge(x, y, by = "pid", all = T), list(total, coexp_stromal, coexp_tumor))

# Get correlation and melt
cor_mat = cor(total_mat[, 2:ncol(total_mat)], use = "complete.obs", method = "spearman")
mat_compl = as.matrix(total_mat[, 2:ncol(total_mat)][complete.cases(total_mat[, 2:ncol(total_mat)])])
cor_mat_p = rcorr(mat_compl, type = "spearman")$P
write.table(cor_mat_p, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/correlation-matrix-pvalues.tsv",
            quote = F, col.names = NA, row.names = T, sep = "\t")

cor_melt = melt(cor_mat)
setDT(cor_melt)

# Set factor
cor_melt[, Var2 := factor(Var2, levels = rev(levels(Var2)))]

# Plot matrix
plt = ggplot(cor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis("Spearman's\ncorrelation", direction = 1) +
  scale_y_discrete("") + 
  scale_x_discrete("", position = "top") +
  theme_minimal() +
  theme(text = element_text(family = "Helvetica"),
        axis.text.x.top = element_text(angle = 45, hjust = 0))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/correlation-matrices/ihc+tils+etils-correlation_matrix",
               width = 9, height = 7)  

# Plot heatmap with CNA
clusters = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/sample-clustering.txt")
total_clustered = merge(total_mat, clusters, by.x = "pid", by.y = "sample", all.x = T)

total_cluster_m = melt(total_clustered[!is.na(cluster)], id.vars = c("pid", "cluster"))

setorder(total_cluster_m, cluster, pid)
total_cluster_m[, pid := factor(pid, levels = unique(pid))]
ggplot(total_cluster_m, aes(x = pid, y = variable, fill = value)) +
  geom_tile() 

# Get mean values per cluster
datacols = names(total_mat)[2:ncol(total_mat)]
cluster_means = total_clustered[!is.na(cluster), lapply(.SD, mean, na.rm = T), by = "cluster", .SDcols = datacols]

# Melt and plot
cluster_means = melt(cluster_means, id.vars = "cluster")
cluster_means[, cluster := factor(cluster)]

# Split into sections
cluster_means[grepl("TIL", variable), section := "TILs"]
cluster_means[grepl("Stroma", variable), section := "Stromal IHC"]
cluster_means[grepl("Tumor", variable), section := "Tumoral IHC"]

plt1 = ggplot(cluster_means[grepl("TIL", variable),], aes(x = cluster, y = variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient2(high = "red") +
  labs(y = "", x = "Cluster", fill = "Mean Percentage (%)")

plt2 = ggplot(cluster_means[grepl("Stroma", variable),], aes(x = cluster, y = variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient2(high = "red") +
  labs(y = "", x = "Cluster", fill = "Mean IHC cell density\n(cells/mm2)")


plt3 = ggplot(cluster_means[grepl("Tumor", variable),], aes(x = cluster, y = variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient2(high = "red") +
  labs(y = "", x = "Cluster", fill = "Mean IHC cell density\n(cells/mm2)")

plt = plot_grid(plt1, plt2, plt3, ncol = 1, align = "v")
save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/clusters/mean-ihc_til-values-per_cluster",
              height = 10, width = 7)
