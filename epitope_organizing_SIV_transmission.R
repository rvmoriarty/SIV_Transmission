# Written by Ryan V. Moriarty in R 4.1.0

options(stringsAsFactors = F)
if("vegan" %in% installed.packages() == F){install.packages("vegan")}
if("stringr" %in% installed.packages() == F){install.packages("stringr")}

library(vegan)
library(stringr)

setwd("/Users/ryanmoriarty/Documents/SIV_transmission_reanalysis_info/SIV_transmission_reanalysis/Epitope_analysis/outputs")

# Grab the .csv files you made in txt_to_modified_csv.ipynb 
csvs = list.files(pattern = "_.csv")
csvs_list = lapply(csvs, read.csv)
names(csvs_list) = gsub(".csv", "", csvs)

# Add a column telling the frequencies of each epitope in a given sample 
epitope_percentage = data.frame()
for(c in 1:length(csvs_list)){
  tmp = csvs_list[[c]]
  for(t in 1:nrow(tmp)){
    tmp[t, "freq"] = as.numeric(tmp[t, "barcode_count"])/sum(tmp$barcode_count)
  }
  epitope_percentage = rbind(epitope_percentage, tmp)
}
# Write the results to a .csv file 
write.csv(epitope_percentage, "epitope_percentages.csv")

# Organize by animal and epitope 
epitope_percentage$rep = "rep1" # there is only one replicate for each sample here
by_animal_epitope_freqs = list()
animals = unique(epitope_percentage$animal)
l = 1
for(a in animals){
  tmpdf = epitope_percentage[which(epitope_percentage$animal == a), ]
  dpi = unique(tmpdf$dpi)
  reps = unique(tmpdf$rep)
  epitopes = unique(tmpdf$epitope)
  animal_list = list()
  i = 1
  for(e in epitopes){
    seqs = unique(tmpdf[which(tmpdf$epitope == e), "barcode_sequence"])
    t = tmpdf[which(tmpdf$epitope == e & tmpdf$animal == a),]
    epitopedf = data.frame()
    for(s in seqs){
      for(d in dpi){
        for(r in reps){
          epitopedf[s, "Sequence"] = s
          epitopedf[s, "animal"] = a
          
          if(length(t[which(t$barcode_sequence == s & t$dpi == d & t$rep == r), "freq"]) > 0){
            epitopedf[s, paste(d,r,sep = "_")] = as.numeric(t[which(t$barcode_sequence == s & t$dpi == d & t$rep == r), "freq"])
          }else{
            epitopedf[s, paste(d,r,sep = "_")] = NA
          }
        }
      }
    }
    epitopedf$rowsum = rowSums(epitopedf[,3:ncol(epitopedf)], na.rm = T)
    epitopedf = epitopedf[-which(epitopedf$rowsum < 0.01), ]
    other_epitope = c("OTHER", unique(epitopedf$animal), 1- colSums(epitopedf[,grep("dpi", colnames(epitopedf)) ], na.rm = T), 0)
    epitopedf = rbind(epitopedf, other_epitope)
    # For each animal, write a .csv for each organized epitope
    write.csv(epitopedf, paste(a, e, "organized_epitopes_all.csv", sep = "_"))
    animal_list[[i]] = epitopedf
    i = i + 1
  }
  names(animal_list) = epitopes
  by_animal_epitope_freqs[[l]] = animal_list
  l = l + 1
}
