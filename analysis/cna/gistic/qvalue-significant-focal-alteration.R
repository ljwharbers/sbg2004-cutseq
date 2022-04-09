## Author: Luuk Harbers
## Date: 2020-07-17
## Script for extraction of GISTIC qvalues while subsetting for genes in a specified genelist

## Load/install packages
packages = c("data.table", "RColorBrewer")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load in data
amp_genes = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/GISTIC/GISTIC2-output/BICRO205+206_MS18+NZ24/amp_genes.conf_99.txt",
                  drop = 1)
del_genes = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/GISTIC/GISTIC2-output/BICRO205+206_MS18+NZ24/del_genes.conf_99.txt",
                  drop = 1)
genelist = fread("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/genes/genelist_nik-zainal-etal.tsv")

# Notate datasets for looping (only note the gene level ones)
datasets = c("amp_genes", "del_genes")

# Apply through datasets and over columns to get all genes with their respective q values
tmp = lapply(datasets, function(y){
  lapply(colnames(get(y)), function(x) {
    tmp = data.table(get(y)[5:nrow(get(y)), ..x], get(y)[2, ..x], y)
    setnames(tmp, c("gene", "qvalue", "variable"))
    tmp = tmp[grepl(".", tmp$gene),]
    return(tmp)
  })
})

# Bind the nested list of output
tmp = lapply(1:length(tmp), function(x){
  rbindlist(tmp[[x]])
})

total = rbindlist(tmp)
total[, qvalue := as.numeric(qvalue)]
total_subset = unique(total[gene %in% genelist$Gene])

# Transform q values to -log10
total_subset[, q_transf := -log10(qvalue)]
setorder(total_subset, -q_transf)
total_subset[, gene := factor(gene, levels = gene)]


plt = ggplot(total_subset, aes(x = gene, y = q_transf, fill = variable)) +
  geom_col() +
  geom_hline(yintercept = -log10(0.05), linetype = 2, size = 1) +
  labs(y = "-log10(q-value)",
       x = "") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/gistic/significantly-focal-alterations",
              height = 6, width = 7)
