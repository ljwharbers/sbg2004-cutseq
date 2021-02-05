## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for hierarchical clustering of CN segments

## Load/install packages
packages = c("data.table", "ggdendro", "ggrastr")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

data = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/BICRO205+206_MS18+NZ24_1e5_gatk-cnv.rds")

# Extract columns
data_cols = lapply(1:length(data), function(x) setnames(data[[x]], "MEAN_LOG2_COPY_RATIO", names(data[x])))
segments = do.call(cbind, lapply(data_cols, function(x) x[, 6]))

# Transpose
segments_t = t(segments)

#create distance and cluster
dd = dist(scale(segments_t), method = "euclidean")
hc = hclust(dd, method = "ward.D2")
clust = cutree(hc, k=2)
clust_dt = data.table(label = names(clust), cluster = factor(clust))

#make dendrogram
dendr = dendro_data(hc, type = "rectangle") 
dendr[["labels"]] = merge(dendr[["labels"]], clust_dt, by = "label")

#plot
plt1 = ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color = cluster), 
            size=4) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
  ylab("Height") +
  scale_color_brewer(palette = "Set1") +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

save_and_plot(plt1, "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/hierarchical-clustering/CN-segments_clustering",
              height = 12, width = 7)

# Calculate CNA percentage
total_rows = nrow(data[[1]])
altered = lapply(1:length(data), function(i){
  data.table(sample = names(data[i]),
             Amplified = (sum(data[[i]]$CALL == "+") / total_rows) * 100,
             Deleted = (sum(data[[i]]$CALL == "-") / total_rows) * 100)
})
total = rbindlist(altered)

# Plot CNA percentage for different clusters
total_merged = merge(total, clust_dt, by.x = "sample", by.y = "label")
#total_merged[, Aneuploid := Amplified + Deleted]

total_m = melt(total_merged)

# Plot Percentage
plt2 = ggplot(total_m, aes(x = cluster, y = value, color = variable)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2)) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(y = "Percentage of genome altered (%)", x = "Cluster", color = "") +
  theme(legend.position = "top")

save_and_plot(plt2, 
              "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/cna_burden/BICRO205+206_MS18+NZ24_100kb_CNAburden_clustered",
              height = 7, width = 7)

# Plot dendrogram with genomewide heatmap
# Load in bins
bins = fread(paste0("/mnt/AchTeraD/Documents/bins/fullbins/hg19_100kbBins_nochr.bed"))
setnames(bins, 1:3, c("chr", "start", "end"))

# Add 1 for consistency between datasets
bins[, start := start+1]

#remove X+Y
bins = bins[chr != "Y",]

#setorder
bins = bins[naturalsort::naturalorder(chr)]
setkey(bins, chr, start, end)

#get x-axis tick locations, switch locations
ticks = numeric(0)
for(k in c(as.character(1:22), "X")){
  ticks[k] = nrow(bins[bins$chr == k]) / 2
}

#get switch point locations
chrSwitch = c(0, which(bins$chr != dplyr::lag(bins$chr)))
labels = c("1", rep("", 3), "5", rep("", 4), "10", rep("", 4), "15", rep("", 4), "20", "", "", "X")

setnames(bins, 1:3, c("chr", "start", "end"))
bins[, feature := paste0(chr, ":", start, "-", end)]

# Add feature in segments dt
segments[, feature := data[[1]][, FEATURE]]

#add missing rows
to.add = bins[!bins$feature %in% segments$feature]
segments_total = rbind(segments, data.table(feature = to.add$feature), fill = T, use.names=TRUE)

#set factor levels and melt
segments_total = segments_total[naturalsort::naturalorder(feature)]
segments_total[, feature := factor(feature, levels = feature)]

# Melt and set factor levels
segments_total = melt(segments_total, id.vars = "feature")

# Set values above/below (-)2 to (-)2
segments_total[value > 1.5, value := 1.5]
segments_total[value < -1.5, value := -1.5]

segments_total[, variable := factor(variable, levels = dendr$labels$label[order(dendr$labels$x)])]

# Plot
plt3 = ggplot(segments_total, aes(x = feature, y = variable, fill = value)) +
  geom_tile_rast() +
  geom_vline(xintercept = chrSwitch[2:length(chrSwitch)], linetype = 2, size = 0.6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1.5, 1.5), na.value = "white",
                       name = "log2 ratio", labels = c("< -1.5", -1, 0, 1, "> 1.5"), breaks = c(-1.5, -1, 0, 1, 1.5)) +
  scale_x_discrete(drop = F, breaks = bins$feature[chrSwitch + ticks], labels = labels) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())

plt4 = ggplot(segment(dendr)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0.013, 0.013)) +
  theme_dendro()

cowplot::plot_grid(plt4, plt3, align = "h", rel_widths = c(0.1, 1), ncol = 2)

