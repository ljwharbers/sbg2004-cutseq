## Author: Luuk Harbers
## Date: 2020-07-15
## Script for analysis of frequency of genome alterations genome wide in replicates

## Load/install packages
packages = c("data.table", "pbapply", "ggplot2", "naturalsort")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

## Variables
threads = 20

## Load replicates
rep1 = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO205+206_MS18_1e5_gatk-cnv.rds")
rep2 = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO206_NZ24_1e5_gatk-cnv.rds")

## specify binsize
binsize = "100kb"

## Load bins
bins = fread(paste0("/mnt/AchTeraD/Documents/bins/fullbins/hg19_", binsize, "Bins_nochr.bed"))
setnames(bins, 1:3, c("CONTIG", "START", "END"))
bins = bins[CONTIG != "Y"]
bins[, START := START + 1]
bins[, FEATURE := paste0(CONTIG, ":", START, "-", END)]

## Reorder
bins = bins[naturalorder(FEATURE, )]

#get x-axis tick locations, switch locations
ticks = numeric(0)
for(k in c(1:22, "X")){
  ticks[k] = nrow(bins[bins$CONTIG == k]) / 2
}

#get switch point locations
chrSwitch = c(0, which(bins$CONTIG != dplyr::lag(bins$CONTIG)))
labels = c(1:22, "X")

## Get samples in both replicates
samples = intersect(names(rep1), names(rep2))
total_rows = nrow(rep1[[1]])

## Get contig information
contigs = rep1[[1]][, c(1:3, 8)]

## lapply through lists and extract 'CALL' column
rep1_dt = pblapply(samples, function(sample){
  dt = rep1[[sample]][, "CALL"]
  setnames(dt, "CALL", sample)
  return(dt)
})
## Bind together and add contigs/features
rep1_total = dplyr::bind_cols(rep1_dt)

## lapply through lists and extract 'CALL' column
rep2_dt = pblapply(samples, function(sample){
  dt = rep2[[sample]][, "CALL"]
  setnames(dt, "CALL", sample)
  return(dt)
})
rep2_total = dplyr::bind_cols(rep2_dt)

## Get frequencies
freq = function(f){
  cna = sum(f != "0")
  data.table(cna = cna)
}

rep1_freq = Reduce(rbind, apply(rep1_total, 1, freq))
rep2_freq = Reduce(rbind, apply(rep2_total, 1, freq))
setnames(rep1_freq, "replicate 1")
setnames(rep2_freq, "replicate 2")

frequency = cbind(contigs, rep1_freq, rep2_freq)

## Add missing bins
add_bins = bins[!FEATURE %in% frequency$FEATURE,]
frequency = rbind(frequency, add_bins, fill = T)

## Reorder and set factor levels
frequency = frequency[naturalorder(FEATURE),]
frequency[, FEATURE := factor(FEATURE, levels = rev(FEATURE))]

## Melt
freq_melt = melt(frequency, id.vars = "FEATURE", measure.vars = c("replicate 1", "replicate 2"))

## Prepare for plotting
freq_melt[, value := (value / length(samples)) * 100]
freq_melt[variable == "replicate 1", value := value * -1]

## Plot
plt = ggplot(freq_melt, aes(x = FEATURE, y = value, color = variable)) +
  geom_col() +
  geom_vline(xintercept = (nrow(bins) - chrSwitch[2:length(chrSwitch)]), linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.6) +
  labs(y = "Percentage of samples with alteration", fill = "") +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks]) +
  scale_y_continuous(limits = c(-75, 75), breaks = seq(-75, 75, by = 25), 
                     labels = paste0(abs(seq(-75, 75, by = 25)), "%")) + 
  coord_flip() +
  scale_color_viridis_d(begin = 0.2, end = 0.5) +
  theme(legend.position = "top",
        line = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/replicates/Genome-alteration-frequency",
              height = 12, width = 8)

## Save empty plot for axes
plt_empty = ggplot(freq_melt, aes(x = FEATURE, y = value, color = variable)) +
  geom_path() +
  geom_vline(xintercept = (nrow(bins) - chrSwitch[2:length(chrSwitch)]), linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.6) +
  labs(y = "Percentage of samples with alteration", fill = "") +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks]) +
  scale_y_continuous(limits = c(-75, 75), breaks = seq(-75, 75, by = 25), 
                     labels = paste0(abs(seq(-75, 75, by = 25)), "%")) + 
  coord_flip() +
  scale_color_viridis_d(begin = 0.2, end = 0.5) +
  theme(legend.position = "top",
        line = element_blank(),
        axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"))

save_and_plot(plt_empty, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/replicates/Genome-alteration-frequency_EMPTY",
              height = 12, width = 8)


## Save only bars
plt_graph = ggplot(freq_melt, aes(x = FEATURE, y = value, color = variable)) +
  geom_col() +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks]) +
  scale_y_continuous(limits = c(-75, 75), breaks = seq(-75, 75, by = 25), 
                     labels = paste0(abs(seq(-75, 75, by = 25)), "%")) + 
  coord_flip() +
  scale_color_viridis_d(begin = 0.2, end = 0.5) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

save_and_plot(plt_graph, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/replicates/Genome-alteration-frequency_GRAPHONLY",
              height = 12, width = 8)
