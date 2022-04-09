#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: 2020-07-14
## Script for creating of rds files after gatk-cnv calling

## Load/install packages
packages = c("data.table", "argparser", "pbapply", "RColorBrewer", "naturalsort")
if(FALSE %in% sapply(packages, require, character.only = T)) stop("Failed while loading packages")

invisible(source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R"))

## Parse arguments
parser = arg_parser("Generation of full rds files after gatk-cnv calling")
parser = add_argument(parser, "--dots", help = "Path to raw datapoints of gatk-cnv analysis, ends with '*denoisedCR.tsv'")
parser = add_argument(parser, "--calls", help = "Path to cnv calls of gatk-cnv analysis, ends with '*called.seg'")
parser = add_argument(parser, "--binsize", help = "Binsize used for analysis in scientific notation (i.e. 1e5)")
parser = add_argument(parser, "--annotation", help = "List of sample IDs and matching barcodes, if not present uses barcode as sample id,
                      first column should be ID and second column the matching barcode")
parser = add_argument(parser, "--rdspath", help = "Output path for rds files (i.e. path/to/outputrds/samplename.rds)")
parser = add_argument(parser, "--threads", help = "Number of threads to use")
parser = add_argument(parser, "--plot", flag = TRUE, short = "-p", help = "Plot profiles")
parser = add_argument(parser, "--plotdir", help = "Output directory for plots")
paser = add_argument(parser, "--binpath", default = "/mnt/AchTeraD/Documents/bins/fullbins/", help = "Path to bin files")


argv = parse_args(parser)

## Check if rds files exist, if not create
if(!file.exists(argv$rdspath)) {
  ## Check for args
  if(is.na(argv$dots)) stop("Directory path of --dots not specified")
  if(is.na(argv$calls)) stop("Directory path of --calls not specified")
  if(is.na(argv$annotation)) stop("Filepath of --annotation not specified")
  
  print("Generating rds file ...")
  ## List files and load in annotation
  dots = list.files(argv$dots, pattern = "denoisedCR.tsv", full.names = T)
  calls = list.files(argv$calls, pattern = "called.seg", full.names = T)
  annotation = fread(argv$annotation)
  
  ## Create output directory
  dir.create(gsub("\\/[^\\/]+$", "", argv$rdspath), showWarnings = F)
  
  ## Get ordered names/barcodes
  names_calls = as.data.table(gsub(".*\\/|.dedup.*", "", calls))
  ordered_names = merge(names_calls, annotation, all.x = T, by.y = eval(names(annotation)[2]), by.x = "V1")
  
  ## Read in data
  calls_dt = lapply(calls, function(x) {
    fread(x, skip = "CONTIG")
  })
  names(calls_dt) = ordered_names[[2]]
  
  dots_dt = lapply(dots, function(x) {
    fread(x, skip = "CONTIG")
  })
  names(dots_dt) = ordered_names[[2]]
  
  ## Subset data and setkeys
  calls_dt = pblapply(calls_dt, function(x) {
    setkey(x, CONTIG, START, END)
    x[!grepl("GL|MT|Y", x$CONTIG)]
  }, cl = argv$threads)
  
  dots_dt = pblapply(dots_dt, function(x) {
    setkey(x, CONTIG, START, END)
    x[!grepl("GL|MT|Y", CONTIG)]
  }, cl = argv$threads)
  
  ## Check if reorder necessary
  if(!all.equal(names(dots_dt), names(calls_dt))) {
    calls_dt[names(dots_dt)]
  }
  
  ## Get overlaps and combine into one list
  combined = pblapply(1:length(calls_dt), function(x) {
    foverlaps(calls_dt[[x]], dots_dt[[x]])[, c(1:4, 7:9)]
  }, cl = argv$threads)
  names(combined) = names(calls_dt)
    
  ## Remove previous lists
  rm(calls_dt, dots_dt)
  
  ## Set levels and set NA logratios for rows that a ratio of 1 or less bins
  combined = lapply(combined, function(x) {
    x[, FEATURE := paste0(CONTIG, ":", START, "-", END)]
    x[NUM_POINTS_COPY_RATIO <= 1, MEAN_LOG2_COPY_RATIO := NA]
  })
  
  ## Save RDS file with all information and clean up
  saveRDS(combined, argv$rdspath)
  invisible(gc())
} else print("rds file already exist, skipping...")

## Generate plots if plot = TRUE
if(argv$plot) {
  ## Check for args
  if(is.na(argv$binsize)) stop("Binsize argument not specified")
  if(is.na(argv$plotdir)) stop("Plotdir argument not specified")
  if(is.na(argv$binsize)) stop("Binsize argument not specified")
  
  print("Generating plots ...")
  dir.create(argv$plotdir, recursive = T, showWarnings = F)
  
  ## Load in rds if not in environment
  if(!exists("combined")) combined = readRDS(argv$rdspath)

  ## Load in full bins 
  options(scipen = 999)
  binsize = paste0(1e5 / 1000, "kb")
  bins = fread(paste0("/mnt/AchTeraD/Documents/bins/fullbins/hg19_", binsize, "Bins_nochr.bed"))
  setnames(bins, c("CONTIG", "START", "END"))
  
  ## Add 1 for consistency between datasets
  bins[, START := START + 1]
  
  ## Remove Y
  bins = bins[!grepl("Y", CONTIG)]
  
  ## Add feature column and set keys and levels
  bins[, FEATURE := paste0(CONTIG, ":", START, "-", END)]
  bins[, FEATURE := factor(FEATURE, levels = FEATURE)]
  
  ## Add missing rows in dataset
  add = lapply(combined, function(x) {
    bins[!bins$FEATURE %in% x$FEATURE, ]
  })
  
  ## rbind missing rows
  combined = pblapply(1:length(combined), function(x) {
    rbind(combined[[x]], add[[x]], fill = T)
  }, cl = argv$threads)
  names(combined) = names(add)
  
  #remove 'add' list
  rm(add)
  
  #reorder
  combined = pblapply(combined, function(dt) {
    setkey(dt, CONTIG, START, END)
    dt[naturalorder(CONTIG)]
  }, cl = argv$threads)
  
  combined = pblapply(combined, function(dt) {
    dt[, FEATURE := factor(FEATURE, levels = FEATURE)]
  }, cl = argv$threads)
  
  combined = pblapply(combined, function(x) {
    x[CALL == "+", color := "#E41A1C"]
    x[CALL == "-", color := "#377EB8"]
    x[CALL == "0", color := "black"]
  })
  
  ## Get chromosome switch point locations
  chrSwitch = c(0, which(combined[[1]]$CONTIG != dplyr::lag(combined[[1]]$CONTIG)))
  
  ## Get x-axis tick locations, switch locations
  ticks = sapply(c(as.character(1:22), "X"), function(x) nrow(combined[[1]][CONTIG == x]) / 2)
  
  ## Generate plots
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
  }, cl = argv$threads)
  
  ## Save plots
  #save plots
  pblapply(names(plots_full), function(x) {
    # png(filename = paste0(argv$plotdir, x, ".png"), units="in", height=10, width=16, res=300)
    # plots_full[[x]]
    # dev.off()
    ggsave(filename = paste0(argv$plotdir, x, ".png"), plot = plots_full[[x]], dpi = 300, units = "in", height = 10, width = 16)
  }, cl = argv$threads)
}