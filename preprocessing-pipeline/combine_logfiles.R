#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script to concatenate log files

## Load/install packages
packages = c("data.table", "argparser")
invisible(lapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

## Parse arguments
parser = arg_parser("Sum read numbers of log files from different lanes and output single log file")
parser = add_argument(parser, "--files", help = "Path to logfiles", nargs = Inf)
parser = add_argument(parser, "--output", help = "Path to logfiles", nargs = Inf)

argv = parse_args(parser)

log = lapply(argv$files, function(file) { 
  fread(file, sep = ":")
  })

full_log = rbindlist(log)

combined = full_log[, sum(V2), by = V1]

write.table(combined, argv$output, sep = ": ", col.names = F, row.names = F, quote = F)