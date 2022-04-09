## Author: Luuk Harbers
## Date: 2020-07-15
## Script for pearson correlation analysis of replicates at different resolution

## Load/install packages
packages = c("data.table", "pbapply", "ggplot2")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

## Variables
rep1= "BICRO205\\+206_MS18_"
rep2= "BICRO206_NZ24"
threads = 20

## list directories to compare
files_rep1 = list.files("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/", 
                        recursive = F, full.names = T, pattern = rep1)
files_rep2 = list.files("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/data/rds-files/", 
                        recursive = F, full.names = T, pattern = rep2)

## Check if length is equal
if(length(files_rep1) != length(files_rep2)) stop("Replicate 1 and 2 not of equal length")

## Apply through files and extract mean copy number ratio of all samples for each resolution
## Rep1
data_rep1 = lapply(files_rep1, function(x){
  dt = readRDS(x)
  ratios = pblapply(1:length(dt), function(i){
    cn = dt[[i]][, "MEAN_LOG2_COPY_RATIO"]
    setnames(cn, "MEAN_LOG2_COPY_RATIO", names(dt[i]))
    return(cn)
  }, cl = threads)
  dplyr::bind_cols(ratios)
})
names(data_rep1) = gsub(paste0(".*", rep1, "|_gatk.*|_"), "", files_rep1)

## Rep2
data_rep2 = lapply(files_rep2, function(x){
  dt = readRDS(x)
  ratios = pblapply(1:length(dt), function(i){
    cn = dt[[i]][, "MEAN_LOG2_COPY_RATIO"]
    setnames(cn, "MEAN_LOG2_COPY_RATIO", names(dt[i]))
    return(cn)
  }, cl = threads)
  dplyr::bind_cols(ratios)
})
names(data_rep2) = gsub(paste0(".*", rep2, "|_gatk.*|_"), "", files_rep2)

## Check if names are the same
if(names(data_rep1) != names(data_rep2)) stop("Replicate differ in (order of) column names")

## Get colnames from both replicates
colnames_rep1 = names(data_rep1[[1]])
colnames_rep2 = names(data_rep2[[1]])
shared_samples = intersect(colnames_rep1, colnames_rep2)

## Vector of colnames to use for Pearson Correlation
pcc = lapply(1:length(data_rep1), function(i){
  sapply(shared_samples, function(sample){
    cor(data_rep1[[i]][, as.character(sample), with = F], data_rep2[[i]][, as.character(sample), with = F], use = "complete.obs")
  })
})
names(pcc) = names(data_rep1)
total = melt(as.data.table(dplyr::bind_cols(pcc)))


## Set factor levels
total[, variable := factor(variable, levels = c("5e5", "2e5", "1e5", "5e4"))]

## Plot
plt = ggplot(total, aes(x = variable, y = value)) +
  geom_violin(aes(fill = variable, color = variable)) +
  geom_boxplot(width = 0.04, outlier.size = 0.95) +
  scale_fill_viridis_d(begin = 0.3, end = 0.98) +
  scale_color_viridis_d(begin = 0.3, end = 0.98) +
  labs(y = "Pearson Correlation Coefficient",
       x = "Binsize") +
  theme(legend.position = "None")

save_and_plot(plt, 
              "/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/replicates/pearson-correlation-resolutions_rep1vsrep2",
              width = 7, height = 7)
