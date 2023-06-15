#!/usr/bin/env Rscript
# R CLI script to concatenate all cruise filter params into one CSV file.

csv.list <- list.files(path=".", pattern="filterparams.csv", recursive=T, full.names=T)
csv.list <- csv.list[-grep("ALL-filterparams.csv", csv.list)]
DF <- do.call(rbind, lapply(csv.list, function(x) read.csv(x)))

write.csv(DF,"ALL-filterparams.csv", quote=F, row.names=F)

