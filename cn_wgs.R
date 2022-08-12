# Script: Copy number alterations in the HER2E subtype (WGS based)

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
# loading all required data
################################################################################

cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/CN_state_matrix.tsv', sep = '\t', header = TRUE))) 

str(cn.scanb)
#View(head(cn.scanb))

# have to check if position is genome position or chromosome position
# also deal with X chromosome
max(cn.scanb$Position)
hist(cn.scanb$Position) # looks like it is a chromosome position

# correct position or add genome position column



# rename X chromosome to 23
cn.scanb$Chr[cn.scanb$Chr == "X"] <- "23" 
cn.scanb$Chr <- as.numeric(cn.scanb$Chr)
# get chromosome lengths (i simply take the last pos on each)
chr.lengths <- cn.scanb %>% group_by(Chr) %>% summarise(length = max(Position)) %>% as.data.frame()
chr.lengths <- chr.lengths[with(chr.lengths, order(Chr)),] %>% mutate(genome = cumsum(as.numeric(length)))

# add chr lengths to df
cn.scanb <- cn.scanb %>% group_by(Chr) %>% 
    mutate(chr_genomelengths = 
               case_when(Chr==1 ~ chr.lengths$genome[1],
                         Chr==2 ~ chr.lengths$genome[2],
                         Chr==3 ~ chr.lengths$genome[3],
                         Chr==4 ~ chr.lengths$genome[4],
                         Chr==5 ~ chr.lengths$genome[5],
                         Chr==6 ~ chr.lengths$genome[6],
                         Chr==7 ~ chr.lengths$genome[7],
                         Chr==8 ~ chr.lengths$genome[8],
                         Chr==9 ~ chr.lengths$genome[9],
                         Chr==10 ~ chr.lengths$genome[10],
                         Chr==11 ~ chr.lengths$genome[11],
                         Chr==12 ~ chr.lengths$genome[12],
                         Chr==13 ~ chr.lengths$genome[13],
                         Chr==14 ~ chr.lengths$genome[14],
                         Chr==15 ~ chr.lengths$genome[15],
                         Chr==16 ~ chr.lengths$genome[16],
                         Chr==17 ~ chr.lengths$genome[17],
                         Chr==18 ~ chr.lengths$genome[18],
                         Chr==19 ~ chr.lengths$genome[19],
                         Chr==20 ~ chr.lengths$genome[20],
                         Chr==21 ~ chr.lengths$genome[21],
                         Chr==22 ~ chr.lengths$genome[22],
                         Chr==23 ~ chr.lengths$genome[23])) %>% 
    relocate(chr_genomelengths, .after=Position) %>% ungroup()

# add up position + chr length to get genome position
cn.scanb$genome_pos <- cn.scanb$Position + cn.scanb$chr_genomelengths 
cn.scanb <- cn.scanb %>% relocate(genome_pos, .after=chr_genomelengths)

#View(cn.scanb)
str(cn.scanb)
#View(head(cn.scanb))
cn.scanb$chr_genomelengths <- NULL
cn.scanb$genome <- NULL

#save(cn.scanb, file = "~/Desktop/MTP_project/Data/SCAN_B/WGS/cn_scanb.RData")
load(file = "~/Desktop/MTP_project/Data/SCAN_B/WGS/cn_scanb.RData")
################################################################################
# plot each sample
################################################################################

pdf(file = paste("~/Desktop/MTP_project/Data/SCAN_B/WGS/CN_sampleplots.pdf", sep =""),onefile = TRUE, height = 14.8, width = 21)

pb = txtProgressBar(min = 0, max = ncol(cn.scanb)-4, initial = 0, style = 3) 
# for each sample
for(i in c(1:(ncol(cn.scanb)-4))) { 
    setTxtProgressBar(pb,i)
    # loop iteration example
    sampleID <- colnames(cn.scanb)[i+4] #gsub("^.{0,3}", "", colnames(cn.scanb)[i+4])
    #print(sampleID)
    sample.data <- as.data.frame(cn.scanb[,c(1:4,i+4)]) 
    #print(head(sample.data))
    plot <- ggplot(sample.data, aes(x=sample.data$genome, y=sample.data[,sampleID])) + #, group=1 #scale_y_continuous(limits=c(-2, 2), breaks=c(-2, -1, 0, 1, 2)) +
        geom_step() + #geom_path()
        geom_point() +
        ggtitle(paste("Copy number profile: ",gsub("^.{0,3}", "", sampleID),sep="")) +
        geom_vline(xintercept = chr.lengths$genome, linetype="dotted") +
        ylab("copy number state") +
        xlab("genome position (chromosome)") +
        scale_x_continuous(breaks=chr.lengths$genome,
                           labels=chr.lengths$Chr)+
        theme(plot.title = element_text(size = 30),
              axis.text.x = element_text(size = 20),
              axis.title.x = element_text(size = 25),
              axis.text.y = element_text(size = 20),
              axis.title.y = element_text(size = 25),
              legend.position = "none")
    print(plot)
}
dev.off()


################################################################################
# plot for her2e group
################################################################################

# calc loss/gain freqs per group
cn.scanb$freqloss.H2 <- apply(cn.scanb[,5:16], 1, function(x) (length(which(x<2))/ncol(cn.scanb[,5:16]))*-100) # i add a minus to make it easer for plotting
cn.scanb$freqgain.H2 <- apply(cn.scanb[,5:16], 1, function(x) (length(which(x>2))/ncol(cn.scanb[,5:16]))*-100)


plot <- ggplot(cn.scanb, aes(x=genome)) + 
    ggtitle("Frequency of gain/loss CN alterations in luminal BC") +
    geom_line(aes(y = freqgain.H2, color = "red"),size=4) + 
    geom_line(aes(y = freqloss.H2, color = "red"),size=4) + 

    scale_colour_manual(name="PAM50", values = c("#d72e2b", "#07a109", "#2176d5"),
                        labels = c("HER2", "LUMA", "LUMB")) + # get corect colors
    geom_vline(xintercept = chr.lengths$genome, linetype="dotted",size=2) +
    scale_x_continuous(name="Genome position (chromosome)",
                       breaks=chr.lengths$genome,
                       labels=c(as.character(1:22),"X"),
                       limits = c(0,max(cn.data.meta$genome)+50000000),
                       expand = c(0, 0)) +
    scale_y_continuous(name="Alteration frequency (%)", # \n Loss          Gain", # \n Loss & G
                       breaks=c(seq(-100,100,25)),
                       labels=c(100,75,50,25,0,25,50,75,100),
                       expand = c(0, 0),
                       limits = c(-80,80)) +
    theme(text=element_text(size=30),
          legend.title = element_blank(),
          axis.title.y = element_text(vjust = 0.5),
          legend.position = c(0.97, 0.90)) +
    annotate(x=45000000,y=c(-50,50), label=c("Loss","Gain"), 
             geom="text", angle=90, hjust=0.5, size=9, colour=c("black","black")) 


print(plot)

################################################################################
# loading all required data
################################################################################

#del
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S000763_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S001143_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S001347_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S002097_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S002236_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S002369_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S002605_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S003100_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S003516_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S004149_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S004725_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S005529_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))

# # load the cn data foe one sample#
#cn.scanb <- as.data.frame((read.table(file = 'Data/SCAN_B/WGS/ASCAT/S000763_l_d_a.copynumber.caveman.csv', sep = ',', header = FALSE)))
#View(cn.scanb)
# # readme for colnames
#readme <- as.data.frame(read.table(file = 'Data/SCAN_B/WGS/ASCATreadme.tsv', sep = '\t', header = FALSE))
# # downloaded file for probe positions at:
# # ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+.tar.gz
# probe.file <- as.data.frame(read.table(file = 'Data/SCAN_B/WGS/CNV_SV_ref_GRCh38_hla_decoy_ebv_brass6+/ascat/SnpGcCorrections.tsv', sep = '\t', header = TRUE))[,1:3] 
# names(probe.file)[1] <- "ProbeID" # QUestion: This is probe right? what does the probe column express int he original file?
# # set colnames
#names(cn.scanb) <- readme$V1[2:9]
# 
#View(cn.scanb)
# View(head(probe.file))
# str(probe.file)
