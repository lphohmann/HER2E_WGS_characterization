# Script: Driver mutations in the HER2E subtype (WGS based)

# TODO:
# - 

# empty environment
rm(list=ls())

# set working directory to where the data is
setwd("~/Desktop/MTP_project")

#packages
library(ggplot2)
library(tidyverse)
library(readxl)

################################################################################
# Driver waterfall plot etc.
################################################################################

# load the cn data foe one sample
mut.data <- as.data.frame((read_excel("Data/SCAN_B/WGS/CodingAndDriversForJohan_21July22.xlsx")))
View(mut.data)
str(mut.data)
driver.data <- mut.data %>% filter(grepl("yes",mut.data$driver)) %>% 
    select(Sample,VD_Gene,VC) %>% 
    dplyr::rename(sample=Sample,gene=VD_Gene,variant_class=VC)
View(driver.data)

# plot

# Create a vector to save mutation priority order for plotting: CHANGE THIS
mutation_priority <- as.character(unique(driver.data$variant_class))

mutation_priority <- c("nonsense","missense")
# plotting parameters
# 1. mainRecurCutoff accepts a numeric value between 0 and 1, and will only plot genes with mutations in x proportion of samples.
# 2. if there are specific genes of interest those can be specified directly via the plotGenes parameter. Input to plotGenes should be a character vector of a list of genes that are desireable to be shown and is case sensitive. 
# 3. plot only specific samples. This can be achieved via the parameter plotSamples
# 4. the maxGenes parameter will only plot the top x genes and takes an integer value. This is usefull for example if when using the mainRecurCutoff parameter a vector of genes have values at x cutoff and all of them are not desired. 
# 5. the rmvSilent parameter will remove all silent mutations from the data.


# plot for 
# her2 driver
pdf(file = "~/Desktop/MTP_project/Data/SCAN_B/WGS/waterfall_driver_her2.pdf",width=20, height=10)
waterfall(driver.data,fileType = "Custom", variant_class_order=mutation_priority, mainRecurCutoff = 0, maxGenes = 27, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15)
dev.off()

################################################################################
# her2 all wf plot
################################################################################
mut.all <- mut.data %>%
    select(Sample,VD_Gene,VC) %>% 
    dplyr::rename(sample=Sample,gene=VD_Gene,variant_class=VC)

# Create a vector to save mutation priority order for plotting: CHANGE THIS
mutation_priority <- as.character(unique(mut.all$variant_class))

mutation_priority <- c("nonsense","start_lost","stop_lost","missense","ess_splice","5prime_UTR_ess_splice","splice_region","silent")
# make a custom colour pallete
custom_pallete <- c("#0e0421", "#d4136d", "#12e0dd", "#c70c0c", "#2a18cc", "#c7c41e",
                    "#37e019","#a903fc")
# plot for each pam50 subtype
# her2 driver
pdf(file = "~/Desktop/MTP_project/Data/SCAN_B/WGS/waterfall_all_her2.pdf",width=20, height=10)
waterfall(mut.all,fileType = "Custom", variant_class_order=mutation_priority, mainRecurCutoff = 0, maxGenes = 20, plotMutBurden = FALSE, mainGrid = TRUE,main_geneLabSize=15,rmvSilent = TRUE,mainPalette = custom_pallete)
dev.off()
