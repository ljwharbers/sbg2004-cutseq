# Author: Luuk Harbers
# Date: 2020-07-08
# Script for generation of OncoPrint plot of phase II breast cancer samples

# Load/install packages
packages = c("data.table", "ComplexHeatmap", "pbapply", "RColorBrewer")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")


# Read data
genes = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/cna-information/BICRO205+206_MS18+NZ24_100kb_Gene-overlaps.tsv")

#parameters
sample = "BICRO205+206_MS18+NZ24"
binsize = "100kb"

#cores for parallel processing
cl = 20

# List all samples
samples = names(readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO205+206_MS18+NZ24_1e5_gatk-cnv.rds"))
samples = paste0("PID ", samples)
annot = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/ER-status.tsv")

#get gene amplifications/deletions
plotting = genes[, c(6, 5, 4)]
setnames(plotting, c("gene", "sample", "variant_class"))
plotting[, sample := paste0("PID ", sample)]
plotting[variant_class == "+", variant_class := 1]
plotting[variant_class == "-", variant_class := -1]

#bind together
plotting[variant_class == -1, variant_class := "DEL"]
plotting[variant_class == 1, variant_class := "AMP"]
plotting[variant_class == 0, variant_class := ""]

# dcast to matrix
plotting[, sample := factor(sample)]

plotting_df = dcast(plotting, gene ~ sample, value.var = "variant_class", 
                                  fun.aggregate = function(x) paste(x, collapse = ";"))

# Add missing samples
add_samples = samples[!samples %in% colnames(plotting_df)]
plotting_df[, eval(add_samples) := ""]

# To dataframe
setDF(plotting_df)

rownames(plotting_df) = plotting_df$gene
plotting_df$gene = NULL
plotting_df$`PID NA` = NULL

# # # Subset for top n genes
# freq = apply(plotting_df, 1, function(x) sum(grepl("AMP|DEL", x)))
# freq = sort(freq, decreasing = T)
# plotting_df = plotting_df[rownames(plotting_df) %in% names(freq[1:100]),]

plotting_mat = as.matrix(plotting_df)

# Set colors
cols = brewer.pal(3, "Set1")[1:2]
names(cols) = c("AMP", "DEL")

# plot oncoprint
plt = oncoPrint(plotting_mat,
                alter_fun = list(
                  background = alter_graphic("rect", width = 0.9, height = 0.9, fill = "#FFFFFF", col = "black"),
                  AMP = alter_graphic("rect", width = 0.85, height = 0.85, fill = cols["AMP"]),
                  DEL = alter_graphic("rect", width = 0.85, height = 0.85, fill = cols["DEL"])),
                col = cols,
                border = "black",
                show_column_names = T)

# Save
save_and_plot(plt, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/oncoprint/", 
                                     sample, "_", binsize, "_", "_oncoprint"), 
              height = 10, width = 12)

# Seperate for ER- and ER+
annot[, pid := paste0("PID ", pid)]
plt2 = oncoPrint(plotting_mat[, colnames(plotting_mat) %in% annot[er_status == 0, pid]],
                alter_fun = list(
                  background = alter_graphic("rect", width = 0.9, height = 0.9, fill = "#FFFFFF", col = "black"),
                  AMP = alter_graphic("rect", width = 0.85, height = 0.85, fill = cols["AMP"]),
                  DEL = alter_graphic("rect", width = 0.85, height = 0.85, fill = cols["DEL"])),
                col = cols,
                border = "black",
                show_column_names = T)

# Save
save_and_plot(plt2, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/oncoprint/", 
                          sample, "_", binsize, "_", "_oncoprint-erNeg"),
              height = 10, width = 12)

plt3 = oncoPrint(plotting_mat[, colnames(plotting_mat) %in% annot[er_status == 1, pid]],
                alter_fun = list(
                  background = alter_graphic("rect", width = 0.9, height = 0.9, fill = "#FFFFFF", col = "black"),
                  AMP = alter_graphic("rect", width = 0.85, height = 0.85, fill = cols["AMP"]),
                  DEL = alter_graphic("rect", width = 0.85, height = 0.85, fill = cols["DEL"])),
                col = cols,
                border = "black",
                show_column_names = T)

# Save
save_and_plot(plt3, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/oncoprint/", 
                          sample, "_", binsize, "_", "_oncoprint-erPos"),
              height = 10, width = 12)