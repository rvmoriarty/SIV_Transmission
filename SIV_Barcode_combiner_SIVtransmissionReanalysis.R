options(stringsAsFactors = FALSE)

if("vegan" %in% installed.packages() == F){install.packages("vegan")}
if("stringr" %in% installed.packages() == F){install.packages("stringr")}

library(vegan)
library(stringr)

to.df = function(longlist){
  df = data.frame()
  for(item in 1:length(longlist)){
    df = rbind(df, longlist[[item]])
  }
  return(df)
}

setwd("/Users/ryanmoriarty/Documents/SIV_transmission_reanalysis_info/SIV_transmission_reanalysis/combined_cutoff_csvs")

ir_animals = c("r10001","r04103")
iv_animals = c("cy0428", "cy0575")

animals = c(ir_animals, iv_animals)

barcode_csvs = list.files(pattern = "Barcode ")
barcode_csvs_list = lapply(barcode_csvs, read.csv)
barcode_csvs_names = gsub(".fastq_Barcode Analysis_Cutoff.csv", "", barcode_csvs)
names(barcode_csvs_list) = barcode_csvs_names

number_barcodes = data.frame()
for(b in 1:length(barcode_csvs_list)){
  number_barcodes[b, "sample"] = barcode_csvs_names[b]
  number_barcodes[b, "number_sequences"] = barcode_csvs_list[[b]][2, 2]
  number_barcodes[b, "number_barcodes"] = barcode_csvs_list[[b]][2, 3]
  number_barcodes[b, "proportion_unique"] = barcode_csvs_list[[b]][2, 4]
}

write.csv(number_barcodes, "SIV_transmission_reanalysis_num_barcodes.csv")

barcodes_found = data.frame()
sample_barcodes_found = list()
for(b in 1:length(barcode_csvs_list)){
  data = barcode_csvs_list[[b]][grep("SIVmac239M", barcode_csvs_list[[b]][,1]), c(1:4)]
  data = cbind(barcode_csvs_names[b], data)
  colnames(data) = c("sample", "barcode", "sequence", "count", "proportion")
  barcodes_found = rbind(barcodes_found, data)
  sample_barcodes_found[[b]] = data
}

by_animal = list()
for(a in 1:length(animals)){
  by_animal[[a]] = barcodes_found[grep(animals[a], barcodes_found$sample), ]
}
names(by_animal) = animals

samples = sort(unique(barcodes_found$sample))
sample_list = list()
for(i in 1:length(samples)){
  sample_list[[i]] = barcodes_found[which(barcodes_found$sample == samples[i]), ]
}
names(sample_list) = samples

for(i in 1:length(sample_list)){
  for(j in 1:nrow(sample_list[[i]])){
    sample_list[[i]][j, "proportion"] = as.numeric(sample_list[[i]][j, "count"])/sum(as.numeric(sample_list[[i]][,"count"]))
  }
}

sample_sdis = data.frame()
for(i in 1:length(sample_list)){
  sample_sdis[i, "sample"] = samples[i]
  sample_sdis[i, "sdi"] = diversity(as.numeric(sample_list[[i]][,"count"]), index = "simpson")
  sample_sdis[i, "total_sequences"] = number_barcodes[i, "number_sequences"]
  sample_sdis[i, "animal"] = animals[!is.na(str_extract(samples[i], animals))]
  sample_sdis[i, "rep"] = str_extract(samples[i], "rep1|rep2|Rep1|Rep2")
}
print(sample_sdis)
write.csv(sample_sdis, "SIV_transmission_reanalysis_sdi.csv")

barcode_3pc_cutoff = list()
for(s in 1:length(sample_list)){
  sample_list[[s]]$proportion = as.numeric(sample_list[[s]]$proportion)
  data = sample_list[[s]][which(sample_list[[s]]$proportion > 0.03), ]
  otherrow = data.frame(names(sample_list)[s], "other", "other", sum(as.numeric(sample_list[[s]]$count)) - sum(as.numeric(data$count)), 1-sum(data$proportion))
  colnames(otherrow) = colnames(data)
  barcode_3pc_cutoff[[s]] = rbind(data, otherrow)
}
names(barcode_3pc_cutoff) = names(sample_list)

barcode_3pc_cutoff_df = to.df(barcode_3pc_cutoff)

barcode_3pc_cutoff_df$animal = str_extract(barcode_3pc_cutoff_df$sample, "cy0428|cy0575|r10001|r04103")
barcode_3pc_cutoff_df$rep = str_extract(barcode_3pc_cutoff_df$sample, "rep1|rep2")

dpi = c()
for(d in 1:nrow(barcode_3pc_cutoff_df)){
  dpi[d] = str_split(barcode_3pc_cutoff_df[d, "sample"], "_")[[1]][2]
}
dpi = as.numeric(gsub("day", "", dpi))
barcode_3pc_cutoff_df$dpi = dpi
write.csv(barcode_3pc_cutoff_df, "barcode_3pc_cutoff_SIV_transmission_reanalysis.csv")


