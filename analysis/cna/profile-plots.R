require(data.table)
require(ggplot2)
require(naturalsort)
require(pbapply)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

#load in data 
combined = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Scripts/rds-files/BICRO205_206_HQ_100kb.rds")

#outdir
outdir = "/mnt/AchTeraD/data/BICRO250/MS76/gatk-cnv/profile-plots/"
sample = "MS76"
binsize = "100kb"

#set cores
cl = 20

# RUN

#get bins for chrom separation in plots
#load in full bins 
bins <- fread(paste0("/mnt/AchTeraD/Documents/bins/fullbins/hg19_", binsize, "Bins_nochr.bed"))
setnames(bins, c("CONTIG", "START", "END"))

#add 1 for consistency between datasets
bins[, START := START + 1]

#remove Y
bins <- bins[!grepl("Y", CONTIG)]

#add feature column and set keys and levels
bins[, FEATURE := paste0(CONTIG, ":", START, "-", END)]
bins[, FEATURE := factor(FEATURE, levels = FEATURE)]

#add missing rows in dataset
add = lapply(combined, function(x) {
  bins[!bins$FEATURE %in% x$FEATURE, ]
})

#rbind missing rows
combined = pblapply(1:length(combined), function(x) {
  rbind(combined[[x]], add[[x]], fill = T)
}, cl = cl)
names(combined) = names(add)

#remove 'add' list
rm(add)

#reorder
combined = pblapply(combined, function(x) {
  setkey(x, CONTIG, START, END)
  x[naturalorder(CONTIG)]
}, cl = cl)

combined = pblapply(combined, function(x) {
  x[, FEATURE := factor(FEATURE, levels = FEATURE)]
}, cl = cl)

combined = pblapply(combined, function(x) {
  x[CALL == "+", color := "#E41A1C"]
  x[CALL == "-", color := "#377EB8"]
  x[CALL == "0", color := "black"]
})

#get chromosome switch point locations
chrSwitch <- c(0, which(combined[[1]]$CONTIG != dplyr::lag(combined[[1]]$CONTIG)))

#get x-axis tick locations, switch locations
ticks = sapply(c(as.character(1:22), "X"), function(x) nrow(combined[[1]][CONTIG == x]) / 2)

#generate plots
plots_full = pblapply(combined, function(x) {
  ggplot() +
    geom_point(data = x, aes(x = FEATURE, y = LOG2_COPY_RATIO), color = "grey", size = 0.5) +
    geom_point(data = x, aes(x = FEATURE, y = MEAN_LOG2_COPY_RATIO, color = color), size = 0.5) +
    geom_vline(xintercept = chrSwitch[2:length(chrSwitch)], linetype = 2, size = 1) +
    scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, by = 2)) + 
    scale_x_discrete(drop = F, breaks = x$FEATURE[chrSwitch + ticks], labels = c(as.character(1:22), "X")) +
    scale_color_identity() +
    theme(text = element_text(family = "Helvetica"),
          legend.position = "",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.y = element_text(size = 10),
          #axis.text = element_blank(),
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          axis.line.y = element_blank())
}, cl = cl)

plots_clean = pblapply(combined, function(x) {
  ggplot() +
    geom_point(data = x, aes(x = FEATURE, y = LOG2_COPY_RATIO), color = "grey", size = 0.5) +
    geom_point(data = x, aes(x = FEATURE, y = MEAN_LOG2_COPY_RATIO, color = color), size = 0.5) +
    scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, by = 2)) + 
    scale_x_discrete(drop = F, breaks = x$FEATURE[chrSwitch + ticks], labels = c(as.character(1:22), "X")) +
    scale_color_identity() +
    theme(text = element_text(family = "Helvetica"),
          legend.position = "",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          line = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank())
}, cl = cl)

#generate plots
plots_full_black = pblapply(combined, function(x) {
  ggplot() +
    geom_point(data = x, aes(x = FEATURE, y = LOG2_COPY_RATIO), color = "grey", size = 0.5) +
    geom_point(data = x, aes(x = FEATURE, y = MEAN_LOG2_COPY_RATIO, color = "black"), size = 0.5) +
    geom_vline(xintercept = chrSwitch[2:length(chrSwitch)], linetype = 2, size = 1) +
    scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, by = 2)) + 
    scale_x_discrete(drop = F, breaks = x$FEATURE[chrSwitch + ticks], labels = c(as.character(1:22), "X")) +
    scale_color_identity() +
    theme(text = element_text(family = "Helvetica"),
          legend.position = "",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          #axis.text.y = element_text(size = 10),
          #axis.text = element_blank(),
          axis.ticks.y = element_line(),
          axis.ticks.x = element_line(),
          axis.line.y = element_blank())
}, cl = cl)

plots_clean_black = pblapply(combined, function(x) {
  ggplot() +
    geom_point(data = x, aes(x = FEATURE, y = LOG2_COPY_RATIO), color = "grey", size = 0.5) +
    geom_point(data = x, aes(x = FEATURE, y = MEAN_LOG2_COPY_RATIO, color = "black"), size = 0.5) +
    scale_y_continuous(limits = c(-2, 4), breaks = seq(-2, 4, by = 2)) + 
    scale_x_discrete(drop = F, breaks = x$FEATURE[chrSwitch + ticks], labels = c(as.character(1:22), "X")) +
    scale_color_identity() +
    theme(text = element_text(family = "Helvetica"),
          legend.position = "",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          line = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank())
}, cl = cl)


#create directories
if(!dir.exists(paste0(outdir, sample, "/full_color"))){
  dir.create(paste0(outdir, sample, "/full_color"))
}
if(!dir.exists(paste0(outdir, sample, "/clean_color"))){
  dir.create(paste0(outdir, sample, "/clean_color"))
}
if(!dir.exists(paste0(outdir, sample, "/full_black"))){
  dir.create(paste0(outdir, sample, "/full_black"))
}
if(!dir.exists(paste0(outdir, sample, "/clean_black"))){
  dir.create(paste0(outdir, sample, "/clean_black"))
}

#save plots
invisible(pblapply(names(plots_full), function(x) {
  ggsave(filename = paste0(outdir, sample, "/full_color/", x, "_color_full.png"), plot = plots_full[[x]], dpi = 300, units = "in", height = 10, width = 16)
}, cl = cl))

invisible(pblapply(names(plots_clean), function(x) {
  ggsave(filename = paste0(outdir, sample, "/clean_color/", x, "_color_clean.png"), plot = plots_clean[[x]], dpi = 300, units = "in", height = 10, width = 16)
}, cl = cl))

invisible(pblapply(names(plots_full_black), function(x) {
  ggsave(filename = paste0(outdir, sample, "/full_black/", x, "_black_full.png"), plot = plots_full_black[[x]], dpi = 300, units = "in", height = 10, width = 16)
}, cl = cl))

invisible(pblapply(names(plots_clean_black), function(x) {
  ggsave(filename = paste0(outdir, sample, "/clean_black/", x, "_black_clean.png"), plot = plots_clean_black[[x]], dpi = 300, units = "in", height = 10, width = 16)
}, cl = cl))