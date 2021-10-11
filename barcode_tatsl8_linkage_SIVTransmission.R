# examine linkage data set 

options(stringsAsFactors = FALSE)

#### Install the packages needed to get this all to work ####
if("vegan" %in% installed.packages() == F){install.packages("vegan")}
if("stringr" %in% installed.packages() == F){install.packages("stringr")}
if("stringdist" %in% installed.packages() == F){install.packages("stringdist")}

library(vegan)
library(stringr)
library(stringdist)

## Because I often switch from lists to data frames, I made this handy dandy function to convert a list to a data frame with one command 
to.df = function(longlist){
  df = data.frame()
  for(item in 1:length(longlist)){
    df = rbind(df, longlist[[item]])
  }
  return(df)
}

#### Set the working directory so we can find these files #### 
setwd("/Users/ryanmoriarty/Documents/SIV_transmission_reanalysis_info/SIV_transmission_reanalysis/6196_files_linkage")

#### Import the CSVs that we generated in the Jupyter notebook ####
r10001_linkage_df = read.csv("/Users/ryanmoriarty/Documents/SIV_transmission_reanalysis_info/SIV_transmission_reanalysis/6196_files_linkage/r10001_sep_mapped/r10001_TatSL8toBarcodeLinkage_sep.csv")
r04103_linkage_df = read.csv("/Users/ryanmoriarty/Documents/SIV_transmission_reanalysis_info/SIV_transmission_reanalysis/6196_files_linkage/r04103_sep_mapped/r04103_TatSL8toBarcodeLinkage_sep.csv")
linkage_df = rbind(r10001_linkage_df, r04103_linkage_df)

r10001_IR_linkage_df = read.csv("/Users/ryanmoriarty/Documents/SIV_transmission_reanalysis_info/SIV_transmission_reanalysis/6196_files_linkage/r10001_IR_sep_mapped/r10001_IR_TatSL8toBarcodeLinkage_sep.csv")
r04103_IR_linkage_df = read.csv("/Users/ryanmoriarty/Documents/SIV_transmission_reanalysis_info/SIV_transmission_reanalysis/6196_files_linkage/r04103_IR_sep_mapped/r04103_IR_TatSL8toBarcodeLinkage_sep.csv")
IR_linkage_df = rbind(r10001_IR_linkage_df, r04103_IR_linkage_df)

#### Because some time points were listed as "week" rather than "day", convert them to days, then all values to numeric
linkage_df$Time.Point = gsub("week12", "day84", linkage_df$Time.Point)
linkage_df$Time.Point = gsub("week10", "day70", linkage_df$Time.Point)
linkage_df$Time.Point = gsub("week8", "day56", linkage_df$Time.Point)
linkage_df$Time.Point = gsub("week6", "day42", linkage_df$Time.Point)
linkage_df$Time.Point = as.numeric(gsub("day", "", linkage_df$Time.Point))

IR_linkage_df$Time.Point = as.numeric(gsub("day", "", IR_linkage_df$Time.Point))

#### The barcodes identified are currently listed as nucleotide strings, so we will convert them to known barcode IDs. #### 
unique_barcodes = unique(linkage_df$SIVmac239M.Barcode)
unique_barcodes_IR = unique(IR_linkage_df$SIVmac239M.Barcode)

barcode_list = read.delim("239M reference 071619.fasta", header = F, sep = ">")
barcode_list_names = barcode_list[grep("SIV", barcode_list$V2), "V2"]
barcode_list_seqs = barcode_list[grep("SIV", barcode_list$V1, invert = T),"V1" ]
barcode_list_seqs = barcode_list_seqs[which(barcode_list_seqs > 1)]
barcode_list_2 = as.data.frame(cbind(barcode_list_names, barcode_list_seqs))

# Because we merged, there may be some issues with the barcode region, so factoring a hamming distance of 1 can account for barcodes that are only 1nt away from a known barcode
hamming_dist_1_from_known = data.frame()
for(u in 1:length(unique_barcodes)){
  hamming_dist_matrix = t(as.data.frame(stringdistmatrix(unique_barcodes[u], barcode_list_2$barcode_list_seqs)))
  hamming_dist_matrix = as.data.frame(cbind(hamming_dist_matrix, barcode_list_names))
  min_dist = as.numeric(min(hamming_dist_matrix$V1))
  if(min_dist == 0){ 
    hamming_dist_1_from_known[u, "sequence"] = unique_barcodes[u]
    hamming_dist_1_from_known[u, "known_barcode"] = hamming_dist_matrix[which(hamming_dist_matrix$V1 == min_dist), "barcode_list_names"]
    hamming_dist_1_from_known[u, "hamming_dist"] = min_dist
  }else if(min_dist < 2){
    hamming_dist_1_from_known[u, "sequence"] = unique_barcodes[u]
    hamming_dist_1_from_known[u, "known_barcode"] = paste(hamming_dist_matrix[which(hamming_dist_matrix$V1 == min_dist), "barcode_list_names"], "hd", collapse = " ")
    hamming_dist_1_from_known[u, "hamming_dist"] = min_dist
  }else{ # These barcodes do not match a known barcode with a max mismatch of 2, so they are likely PCR errors and should not be considered further. 
    hamming_dist_1_from_known[u, "sequence"] = unique_barcodes[u]
    hamming_dist_1_from_known[u, "known_barcode"] = "HammingDist over 2"
    hamming_dist_1_from_known[u, "hamming_dist"] = min_dist
  }
}

hamming_dist_1_from_known_IR = data.frame()
for(u in 1:length(unique_barcodes_IR)){
  hamming_dist_matrix = t(as.data.frame(stringdistmatrix(unique_barcodes_IR[u], barcode_list_2$barcode_list_seqs)))
  hamming_dist_matrix = as.data.frame(cbind(hamming_dist_matrix, barcode_list_names))
  min_dist = as.numeric(min(hamming_dist_matrix$V1))
  if(min_dist == 0){ 
    hamming_dist_1_from_known_IR[u, "sequence"] = unique_barcodes_IR[u]
    hamming_dist_1_from_known_IR[u, "known_barcode"] = hamming_dist_matrix[which(hamming_dist_matrix$V1 == min_dist), "barcode_list_names"]
    hamming_dist_1_from_known_IR[u, "hamming_dist"] = min_dist
  }else if(min_dist < 2){
    hamming_dist_1_from_known_IR[u, "sequence"] = unique_barcodes_IR[u]
    hamming_dist_1_from_known_IR[u, "known_barcode"] = paste(hamming_dist_matrix[which(hamming_dist_matrix$V1 == min_dist), "barcode_list_names"], "hd", collapse = " ")
    hamming_dist_1_from_known_IR[u, "hamming_dist"] = min_dist
  }else{ # These barcodes do not match a known barcode with a max mismatch of 2, so they are likely PCR errors and should not be considered further. 
    hamming_dist_1_from_known_IR[u, "sequence"] = unique_barcodes_IR[u]
    hamming_dist_1_from_known_IR[u, "known_barcode"] = "HammingDist over 2"
    hamming_dist_1_from_known_IR[u, "hamming_dist"] = min_dist
  }
}

#### Now that we have the list of found nucleotide sequences and their respective barcode name, let's include that into our linkage df ####
for(f in 1:nrow(linkage_df)){
  linkage_df[f, "SIVmac239M.Barcode"] = hamming_dist_1_from_known[which(hamming_dist_1_from_known$sequence == linkage_df[f, "SIVmac239M.Barcode"]), "known_barcode"]
}

for(f in 1:nrow(IR_linkage_df)){
  IR_linkage_df[f, "SIVmac239M.Barcode"] = hamming_dist_1_from_known_IR[which(hamming_dist_1_from_known_IR$sequence == IR_linkage_df[f, "SIVmac239M.Barcode"]), "known_barcode"]
}

write.csv(linkage_df, "linkage_df_edited.csv")
write.csv(IR_linkage_df, "IR_linkage_df_edited.csv")

#### We now will organize the data so it is easier to plot #### 
## For each day, first we will determine the proportion of barcode-epitope combinations present as a whole ##

animals = unique(linkage_df$Animal.ID)
animal_list_all_lineage_prop = list()
animal_list_all_lineage_prop_filtered = list()
animal_hd_comb = list()
animal_hd_comb_names = c()
h = 1
for(a in 1:length(animals)){
  animal_df = linkage_df[which(linkage_df$Animal.ID == animals[a]), 2:ncol(linkage_df)]
  dpi = sort(unique(animal_df$Time.Point))
  day_list = list()
  day_list_filtered = list()
  for(d in 1:length(dpi)){
    day_df = animal_df[which(animal_df$Time.Point == dpi[d]), ]
    erroneous_bcs = day_df[which(day_df$SIVmac239M.Barcode == "HammingDist over 2"), ]
    day_df_real = day_df[which(day_df$SIVmac239M.Barcode != "HammingDist over 2"), ]
    if(nrow(day_df_real) > 0){
      reps = unique(day_df_real$Rep.Number)
      rep_list = list()
      rep_list_filtered = list()
      for(r in 1:length(reps)){
        dpi_rep_df = day_df_real[which(day_df_real$Rep.Number == reps[r]), ]
        bcs_found = unique(dpi_rep_df$SIVmac239M.Barcode)
        bc_epitope_comb = paste0(dpi_rep_df$SIVmac239M.Barcode, "_", dpi_rep_df$Tat.SL8.Amino.Acid.Sequence)
        dpi_rep_df = cbind(dpi_rep_df, bc_epitope_comb)
        unique_combs = unique(bc_epitope_comb)
        modified_df = data.frame()
        for(u in 1:length(unique_combs)){
          modified_df[u, "animal"] = animals[a]
          modified_df[u, "dpi"] = dpi[d]
          modified_df[u, "rep"] = reps[r]
          modified_df[u, "barcode_epitope"] = unique_combs[u]
          modified_df[u, "barcode"] = strsplit(unique_combs[u], "_")[[1]][1]
          modified_df[u, "epitope"] = strsplit(unique_combs[u], "_")[[1]][2]
          modified_df[u, "count"] = sum(dpi_rep_df[which(dpi_rep_df$bc_epitope_comb == unique_combs[u]), "count"])
        }
        animal_hd_comb[[h]] = modified_df
        animal_hd_comb_names[h] = paste(animals[a], dpi[d], reps[r], sep = "_")
        h = h + 1 
        modified_df$proportion = modified_df$count/sum(modified_df$count)
        rep_list[[r]] = modified_df
        
        modified_df_filtered = modified_df[which(modified_df$proportion >= 0.03), ]
        newrrow = c(animals[a], dpi[d], reps[r], "other", "other", "other", sum(modified_df$count)-sum(modified_df_filtered$count), 1-sum(modified_df_filtered$proportion))
        modified_df_filtered = rbind(modified_df_filtered, newrrow)
        rep_list_filtered[[r]] = modified_df_filtered
      }
      day_list[[d]] = rep_list
      day_list_filtered[[d]] = rep_list_filtered
    }
  }
  names(day_list) = paste0("day", dpi)
  names(day_list_filtered) = paste0("day", dpi)
  animal_list_all_lineage_prop[[a]] = day_list
  animal_list_all_lineage_prop_filtered[[a]] = day_list_filtered
}
names(animal_hd_comb) = animal_hd_comb_names

animal_comb = to.df(animal_hd_comb)
animal_comb = animal_comb[grep("hd", animal_comb$barcode, invert = T), ]

# This file contains barcodes that are a hamming distance of less than 2 from a known barcode 
write.csv(animal_comb, "animal_comb.csv")

# Only keep the barcodes that are an exact match to our list
write.csv(to.df(animal_hd_comb), "animal_hd_comb.csv")
