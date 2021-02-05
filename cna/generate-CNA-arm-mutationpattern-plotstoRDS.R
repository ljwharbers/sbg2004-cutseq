require(data.table)
require(pbapply)
require(dplyr)
require(grid)
require(GenVisR)
require(RColorBrewer)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

#load data
total = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Scripts/rds-files/BICRO205_206_HQ_100kb.rds")
genes = readRDS("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Scripts/rds-files/BICRO205_206_HQ-100kb_CNA-geneoverlaps.rds")

splits = fread("/mnt/AchTeraD/Documents/Projects/LungCancer/Data/hg19_chromosome_splits.bed")


#parameters
sample = "BICRO205_206_HQ"
binsize = "100kb"

#cores for parallel processing
cl = 20

#run

#add arms to all data.tables
total = pblapply(total, function(x){
  for(i in c(1:22, "X")){
    x[CONTIG == i, ARM := ifelse(END < splits$V2[splits$V1 == i], paste0(i, "p"), paste0(i, "q"))]
  };
  return(x)
})

#get list of medians per chromosomal arm for each data.table
output = data.table(ARM = unique(total[[1]]$ARM))
output = output[-c(42,41)] #remove x rows, can't distinguish male/female for this analysis
output[, names(total) := 0]

# get amplified/deleted arms
lapply(1:length(total), function(x){
  for(i in output$ARM){
    if(nrow(total[[x]][ARM == i & CALL == "+"]) >= nrow(total[[x]][ARM == i]) / 2) output[ARM == i, names(total[x]) := 1]
    if(nrow(total[[x]][ARM == i & CALL == "-"]) >= nrow(total[[x]][ARM == i]) / 2) output[ARM == i, names(total[x]) := -1]
  }
})

#write output to table
write.table(output, paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/alteration-tables/", sample, "_", binsize,
                           "_arm-level-alterations.tsv"), sep = "\t", quote = F, col.names = T, row.names = F)

#melt data
plotting = melt.data.table(output[, 1:(ncol(output))])
setnames(plotting, c("gene", "sample", "variant_class"))

#get gene amplifications/deletions
genes = genes[, c(6, 5, 4)]
setnames(genes, c("gene", "sample", "variant_class"))
genes[variant_class == "+", variant_class := 1]
genes[variant_class == "-", variant_class := -1]


#bind together
plotting = rbind(plotting, genes)

#set names and set variant classes
plotting[variant_class == 1, variant_class := "Amplified"]
plotting[variant_class == -1, variant_class := "Deleted"]
plotting[variant_class == 0, variant_class := NA]

#remove NA rows
plotting = plotting[!is.na(plotting$variant_class)]
plotting = as.data.frame(plotting)


#calculate 'mutation burden' per sample, in this case just the total number of altered arms/genes
mutburden = as.data.frame(table(plotting$sample))
setnames(mutburden, c("sample", "mut_burden"))

#create ggplot2 layer
mut_burden_layer = theme(axis.title.y = element_blank())

#plot
# plt = waterfall(plotting, fileType = "custom", variant_class_order = c("AMP", "DEL"), plotMutBurden = T,
#                  mutBurden = mutburden, mainXlabel = F, mainPalette = brewer.pal(6, "Set1")[c(5, 2)], mainGrid = T,
#                  maxGenes = 328,
#                  mainDropMut = T, mutBurdenLayer = mut_burden_layer, out = "grob")

plt = waterfall(plotting, fileType = "custom", variant_class_order = c("Amplified", "Deleted"),
                plotMutBurden = T, mutBurden = mutburden, mutBurdenLayer = mut_burden_layer, mainXlabel = T, 
                mainGrid = T, mainPalette = brewer.pal(5, "Set1")[c(5, 2)], out = "grob")

plt$grobs[[3]] <- grid.rect(gp = gpar(col = "white"))

grid.draw(plt)

#save
save_and_plot(grid.draw(plt), paste0("/mnt/AchTeraD/Documents/Projects/phaseII-clinicaltrial/Plots/mutation_burden/", 
                                     sample, "_", binsize, "_", "_CNA-arm-mutationplot"), 
              height = 10, width = 20)


