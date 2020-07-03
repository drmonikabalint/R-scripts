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
PFM_annotation <- read.csv2("../ccbr1022_metadata.csv",header = TRUE,sep=",")
PFM_annotation
baseline <- subset(PFM_annotation, treatment == "pre")
post_treatment <- subset(PFM_annotation, treatment =="post")
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

################################################################## 3.2 Data Analysis ############################################
# Subset dataframe generated for the 2 genes of interest (CXCL1 and MCL1) 
# mRNA_norm contains data from both baseline and post treatment visits

table_CXCL1_TFRC <- data.frame(
    patient = PFM_annotation$patient,
    visit = PFM_annotation$treatment,
    CXCL1 = mRNA_norm["CXCL1",],
    TFRC = mRNA_norm["TFRC",])

# Gather results from both genes under the same column
table_tot_unite <- gather(table_CXCL1_TFRC, key = gene, value = 'results', CXCL1:TFRC)
table_tot_unite
# Factorize the datapoints under the column called "gene". This will be used in grouping values for ggplot 
table_tot_unite$gene <- factor(table_tot_unite$gene)
summary(table_tot_unite)

#Generate boxplots for each gene 
boxplot1 <- ggplot(subset(table_tot_unite, gene == "TFRC"),
                   aes(x = visit, y = results, group = visit)) +
    scale_fill_manual(values=c('#007A96', '#d6d925')) + 
    scale_color_manual(values=c('#007A96', '#DFE22B')) + 
    geom_boxplot(aes(fill=visit, color = visit, group = visit), alpha = 0.7) +
    geom_point(aes(fill=visit, color = visit)) +
    ggrepel::geom_text_repel(aes(label = results), size = 2.5) +
    ggtitle("Results for CXCL1 gene by visits") +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 4/4,  legend.position = 'none') +
    scale_x_discrete(name = "")+
    scale_y_continuous(name = "Results")

boxplot2 <- ggplot(subset(table_tot_unite, gene == "CXCL1"),
                   aes(x = visit, y = results, group = visit)) +
    scale_fill_manual(values=c('#007A96', '#d6d925')) + 
    scale_color_manual(values=c('#007A96', '#DFE22B')) + 
    geom_boxplot(aes(fill=visit, color = visit, group = visit), alpha = 0.7) +
    geom_point(aes(fill=visit, color = visit)) +
    ggrepel::geom_text_repel(aes(label = results), size = 2.5) +
    ggtitle("Results for TFRC gene by visits") +
    theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 4/4, legend.position = 'none') +
    scale_x_discrete(name = "")+
    scale_y_continuous(name = "")

ggarrange(boxplot1, boxplot2, widths = c(2,2), align = 'h')

#Summary for the subset tables used in each boxplot

summary(subset(table_tot_unite, gene == "CXCL1"))
summary(subset(table_tot_unite, gene == "TFRC"))

################################## 3.3 Data Reporting ###########################################
# Additional data explorations with ggplot-barplot 
# CXCL1 and MCL1 values by patient

barplot1 <- ggplot(subset(table_tot_unite, gene == "TFRC"), 
                   aes(x=patient, y = results, fill = visit)) + 
    geom_bar(stat="identity",position=position_dodge(), colour="lightgrey") + 
    labs (x= "Patient", y = "Results for TFRC") + 
    scale_fill_manual(values=c('#007A96', '#DFE22B')) +
    ggtitle("TFRC - results by patient") +
    theme(plot.title = element_text(hjust = 0.5),  legend.position = 'bottom', aspect.ratio = 4/4)

barplot2 <- ggplot(subset(table_tot_unite, gene == "CXCL1"), 
                   aes(x=patient, y = results, fill = visit)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="lightgrey") + 
    labs (x= "Patient", y = "Results for CXCL1") + 
    scale_fill_manual(values=c('#007A96', '#DFE22B')) +
    ggtitle("CXCL1 - results by patient") +
    theme(plot.title = element_text(hjust = 0.5), legend.position = 'none', aspect.ratio = 4/4)

ggarrange(barplot1, barplot2, common.legend = TRUE, legend="bottom")

