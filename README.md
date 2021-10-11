# SIV Transmission Study Analysis 

The set of scripts here were used to generate and manipulate SIVmac239M barcoded data from intrarectally (IR) and intravenously (IV) challenged rhesus macaques (RMs) and cynomolgus macaques (CM). 

These scripts use either Jupyter Notebook or R. Jupyter Notebook is free to download, and we recommend downloading Anaconda Navigator
and using Anaconda to access Jupyter Notebook. Anaconda Navigator is available at: https://www.anaconda.com/products/individual. R and RStudio are also free to download. The R version used for these scripts is R 4.1.0 (aside from the Barcoded Virus Analysis App Modified Cutoff.R, which was provided to us courtesy of Dr. Brandon Keele), and the RStudio version is 1.4.1103. All R scripts were run using RStudio. 

The description of each script is listed below. 

# Barcoded Virus Analysis App Modified Cutoff.R

This R script must be run in R-3.3.0 and was provided to us courtesy of Dr. Brandon Keele. This produces a Shiny App that identifies the number and count of unique barcodes present in a .fastq file. This outputs a .csv file containing this information that was then further analyzed using SIV_Barcode_combiner_SIVtransmissionReanalysis.R. 

# SIV_Barcode_combiner_SIVtransmissionReanalysis.R

This R script takes the output of Barcoded Virus Analysis App Modified Cutoff.R and returns the number of unique barcodes present in each day, replicate, and animal, as well as the Simpson's Diversity Index (SDI) of the population. This also filters the barcodes present at a frequency of less than 3% in the population on a given day. 


