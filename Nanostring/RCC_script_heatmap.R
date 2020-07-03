# # # # # # # # # # # # # # # # # #
# # #  Monika Balint. 02 Jul 2020
# # # R version 3.6.3 (2020-02-29) -- "Holding the Windsock"
# # # # # # # # # # # # # # # # #
    
####################################################################### Required R packages ######################################

# Install NanoStringNorm, Bioconductor and other necessary packages
#install.packages("nlme")
#install.packages("spatial")
#install.packages("survival")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("vsn")
#install.packages("NanoStringNorm")
#install.packages("ggplot2", dependencies = TRUE)
#install.packages("gplots")
#install.packages("dplyr")
#install.packages("ggpubr")

# Import packages to the environment
library(NanoStringNorm)
library(ggplot2)
library(gplots)
library(devtools)
library(ggpubr)
library(dplyr)
library(tidyr)

# Set the environment with the location to the input files
setwd("/home/moni/data_science_kaggle/ccbr1022-nanostring-master/RawData/Cytopenia_nanostring/")

##################################################################  Reading sample files #########################################
# Read in RCC files, with read.markup.RCC function. 
# Generate 3 lists: 
#   mRNA_tot - contains results from all visits
#   mRNA_b   - contains results from baseline visits
#   mRNA_pt  - contains results from post treatment visits

mRNA_tot  <- read.markup.RCC(rcc.path = '.') 

# default values are not explicitly written: rcc.pattern = "*.RCC|*.rcc",exclude = NULL,include = NULL,nprobes = -1

# create 2 subset datasets for the 2 visit types: baseline and post-treatment.
PFM_annotation <- read.csv2("case_study_annotations.csv",header = TRUE,sep=",")
baseline <- subset(PFM_annotation, visit == "Baseline")
post_treatment <- subset(PFM_annotation, visit =="Post-Treatment")
mRNA_b  <- read.markup.RCC(rcc.path = 'baseline')
mRNA_pt  <- read.markup.RCC(rcc.path = 'posttreatment')

##################################################################  Normalization ################################################

mRNA_norm <- NanoStringNorm(
    x = mRNA_tot,
    anno = NA,
    CodeCount = 'geo.mean',
    Background = 'mean',
    SampleContent = 'housekeeping.geo.mean',
    round.values = TRUE,
    take.log = TRUE,
    return.matrix.of.endogenous.probes = TRUE
)

head(mRNA_norm, 130)

mRNA_tot[1]
mRNA_tot[2]

################################################################## 3.1 Data Quality ##############################################

# Generate heatmap. 
#To generate heatmap, and transpose the mRNA matrix, with positive and negative genes in the column, 
# and the samples are in row, we use the transpose function (t).

library("RColorBrewer")
col_heat <- grDevices::hcl.colors(20, "viridis")

# create function RCC_heatmap, to avoid copy-pasting code
RCC_heatmap <- function(output_file, input_file){
    pdf(output_file, width=10, height=10)    
    heatmap.2(
        t(input_file),  
        col = col_heat, 
        trace = "none",  
        cexRow = 0.8, 
        scale="none",
        margins=c(4,13)
    )
    dev.off()
}

# create heatmap for mRNA results.
RCC_heatmap("mRNA_tot_heatmap.pdf", mRNA_norm)

