# SIV Transmission Study Analysis 

The set of scripts here were used to generate and manipulate SIVmac239M barcoded data from intrarectally (IR) and intravenously (IV) challenged rhesus macaques (RMs) and cynomolgus macaques (CM). 

These scripts use either Jupyter Notebook or R. Jupyter Notebook is free to download, and we recommend downloading Anaconda Navigator
and using Anaconda to access Jupyter Notebook. Anaconda Navigator is available at: https://www.anaconda.com/products/individual. R and RStudio are also free to download. The R version used for these scripts is R 4.1.0 (aside from the Barcoded Virus Analysis App Modified Cutoff.R, which was provided to us courtesy of Dr. Brandon Keele), and the RStudio version is 1.4.1103. All R scripts were run using RStudio. 

The description of each script is listed below. 

## Barcoded Virus Analysis App Modified Cutoff.R

This R script must be run in R-3.3.0 and was provided to us courtesy of Dr. Brandon Keele. This produces a Shiny App that identifies the number and count of unique barcodes present in a .fastq file. This outputs a .csv file containing this information that was then further analyzed using SIV_Barcode_combiner_SIVtransmissionReanalysis.R. 

## SIV_Barcode_combiner_SIVtransmissionReanalysis.R

This R script takes the output of Barcoded Virus Analysis App Modified Cutoff.R and returns the number of unique barcodes present in each day, replicate, and animal, as well as the Simpson's Diversity Index (SDI) of the population. This also filters the barcodes present at a frequency of less than 3% in the population on a given day. 

## TatSL8_to_barcode_RVM.ipynb

This Jupyter Notebook takes an input of paired FASTQ files, merges the reads, maps them to the SIVmac239M reference file, identifies the unique SIV barcode present, and extracts the TatSL8 epitope sequence. Because the barcode region fand the Tat SL8 epitope are present on the same amplicon, we can use this data to determine linkage. This data is output as a .csv file. 

## barcode_tatsl8_linkage_SIVTransmission.R

Once the .csv file is generated reporting linakge data, this R script is used to further manipulate the data to make plotting and filtering data more simple. It also converts the barcode sequence to the barcode names, as well as removes barcodes that do not match the barcode list exactly. A filter to bin links that are present at a frequency of less than 3% is also implemented, and data is reorganized to make for easier plotting in Prism. 

## epitope_analysis.ipynb

This Jupyter Notebook takes an input of paired Fastq files from the cynomolgus macaques, merges the reads, and maps them to the SIVmac239M reference. This script can identify both Gag GW9 and Nef RM9 epitopes depending on the upstream and downstream flanking sequence you select. Following the analysis, a series of .txt files are generated with data of epitope sequence and counts. This script was modified from the Zequencer Zika virus sequence analysis software developed by Dr. Dave O'Connor's lab at the University of Wisconsin-Madison. The original script is documented and provided at https://bitbucket.org/dhoconno/zequencer/src/master/. 

## txt_to_modified_csv_SIV_transmission.ipynb 

This Jupyter Notebook reads the .txt files made by epitope_analysis.ipynb and identifies the animal, epitope sequence, day post infection, and frequency of each epitope sequence, then combines the data into a .csv file, which is then manipulated using epitope_organzing_SIV_transmission.R. 

## epitope_organizing_SIV_transmission.R 

This R script takes the animal .csv from txt_to_modified_csv_SIV_transmission.ipynb and converts counts to frequencies, filters by frequency, and reorganizes the data to make it easier to plot in Prism. 
