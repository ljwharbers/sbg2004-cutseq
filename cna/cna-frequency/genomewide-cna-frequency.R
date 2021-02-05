## Author: Luuk Harbers
## Date: 2020-07-15
## Script for analysis of frequency of genome alterations genome wide in one cohort (combined replicates)

## Load/install packages
packages = c("data.table", "pbapply", "ggplot2", "naturalsort")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

## Variables
sample= "BICRO205+206_MS18+NZ24"
binsize = "100kb"
threads = 20

## Load replicates
data = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO205+206_MS18+NZ24_1e5_gatk-cnv.rds")

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
total_rows = nrow(data[[1]])

## Get contig information
contigs = data[[1]][, c(1:3, 8)]

## lapply through lists and extract 'CALL' column
dt = pblapply(1:length(data), function(i){
  dt = data[[i]][, "CALL"]
  setnames(dt, "CALL", names(data[i]))
  return(dt)
}, cl = threads)

## Bind together and add contigs/features
total = dplyr::bind_cols(dt)

## Get frequencies
getFreq = function(f){
  amp = sum(f == "+")
  del = sum(f == "-")
  data.table(amp = amp, del = del)
}

freq = Reduce(rbind, apply(total, 1, getFreq))
freq_c1 = Reduce(rbind, apply(total[, colnames(total) %in% clusters[cluster == 1]$sample, with = F], 1, getFreq))
freq_c2 = Reduce(rbind, apply(total[, colnames(total) %in% clusters[cluster == 2]$sample, with = F], 1, getFreq))
frequency = cbind(contigs, freq)
frequency_c1 = cbind(contigs, freq_c1)
frequency_c2 = cbind(contigs, freq_c2)

## Add missing bins
add_bins = bins[!FEATURE %in% frequency$FEATURE,]
frequency = rbind(frequency, add_bins, fill = T)
frequency_c1 = rbind(frequency_c1, add_bins, fill = T)
frequency_c2 = rbind(frequency_c2, add_bins, fill = T)

## Reorder and set factor levels
frequency = frequency[naturalorder(FEATURE),]
frequency[, FEATURE := factor(FEATURE, levels = rev(FEATURE))]

frequency_c1 = frequency_c1[naturalorder(FEATURE),]
frequency_c1[, FEATURE := factor(FEATURE, levels = rev(FEATURE))]

frequency_c2 = frequency_c2[naturalorder(FEATURE),]
frequency_c2[, FEATURE := factor(FEATURE, levels = rev(FEATURE))]

## Melt
freq_melt = melt(frequency, id.vars = "FEATURE", measure.vars = c("amp", "del"))
freq_melt_c1 = melt(frequency_c1, id.vars = "FEATURE", measure.vars = c("amp", "del"))
freq_melt_c2 = melt(frequency_c2, id.vars = "FEATURE", measure.vars = c("amp", "del"))

## Prepare for plotting
freq_melt[, value := (value / length(data)) * 100]
freq_melt[variable == "del", value := value * -1]

freq_melt_c1[, value := (value / sum(clusters$cluster == 1)) * 100]
freq_melt_c1[variable == "del", value := value * -1]

freq_melt_c2[, value := (value / sum(clusters$cluster == 2)) * 100]
freq_melt_c2[variable == "del", value := value * -1]

# Plot
plt = ggplot(freq_melt, aes(x = FEATURE, y = value, color = variable)) +
  geom_col() +
  geom_vline(xintercept = (nrow(bins) - chrSwitch[2:length(chrSwitch)]), linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.6) +
  labs(y = "Percentage of samples with alteration", fill = "") +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "top",
        line = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequency",
              height = 12, width = 8)

plt_clean = ggplot(freq_melt, aes(x = FEATURE, y = value, color = variable)) +
  geom_path() +
  geom_vline(xintercept = (nrow(bins) - chrSwitch[2:length(chrSwitch)]), linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.6) +
  labs(y = "Percentage of samples with alteration", fill = "") +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "top",
        line = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"))

save_and_plot(plt_clean, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequencyEMPTY",
              height = 12, width = 8)

plt_graph = ggplot(freq_melt, aes(x = FEATURE, y = value, color = variable)) +
  geom_col() +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

save_and_plot(plt_graph, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequencyGRAPH",
              height = 12, width = 8)

# Plot
plt = ggplot(freq_melt_c1, aes(x = FEATURE, y = value, color = variable)) +
  geom_col() +
  geom_vline(xintercept = (nrow(bins) - chrSwitch[2:length(chrSwitch)]), linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.6) +
  labs(y = "Percentage of samples with alteration", fill = "") +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "top",
        line = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequency_c1",
              height = 12, width = 8)

plt_clean = ggplot(freq_melt_c1, aes(x = FEATURE, y = value, color = variable)) +
  geom_path() +
  geom_vline(xintercept = (nrow(bins) - chrSwitch[2:length(chrSwitch)]), linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.6) +
  labs(y = "Percentage of samples with alteration", fill = "") +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "top",
        line = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"))

save_and_plot(plt_clean, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequency_c1EMPTY",
              height = 12, width = 8)

plt_graph = ggplot(freq_melt_c1, aes(x = FEATURE, y = value, color = variable)) +
  geom_col() +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

save_and_plot(plt_graph, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequency_c1GRAPH",
              height = 12, width = 8)

# Plot
plt = ggplot(freq_melt_c2, aes(x = FEATURE, y = value, color = variable)) +
  geom_col() +
  geom_vline(xintercept = (nrow(bins) - chrSwitch[2:length(chrSwitch)]), linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.6) +
  labs(y = "Percentage of samples with alteration", fill = "") +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "top",
        line = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"))

save_and_plot(plt, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequency_c2",
              height = 12, width = 8)

plt_clean = ggplot(freq_melt_c2, aes(x = FEATURE, y = value, color = variable)) +
  geom_path() +
  geom_vline(xintercept = (nrow(bins) - chrSwitch[2:length(chrSwitch)]), linetype = 2, size = 0.3) +
  geom_hline(yintercept = 0, size = 0.6) +
  labs(y = "Percentage of samples with alteration", fill = "") +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "top",
        line = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title.y = element_blank(),
        text = element_text(family = "Helvetica"))

save_and_plot(plt_clean, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequency_c2EMPTY",
              height = 12, width = 8)

plt_graph = ggplot(freq_melt_c2, aes(x = FEATURE, y = value, color = variable)) +
  geom_col() +
  scale_x_discrete(drop = F, breaks = frequency$FEATURE[chrSwitch + ticks], label = labels) +
  scale_y_continuous(limits = c(-100, 100), breaks = seq(-100, 100, by = 25), 
                     labels = paste0(abs(seq(-100, 100, by = 25)), "%")) + 
  coord_flip() +
  scale_color_brewer(palette = "Set1") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        legend.position = "none")

save_and_plot(plt_graph, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/alteration-frequency/Genome-alteration-frequency_c2GRAPH",
              height = 12, width = 8)