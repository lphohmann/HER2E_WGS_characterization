# Script: Mutational and rearrangement signatures in the HER2E subtype (WGS based)

# TODO:
# - exclude non-shared signatures and then recalc. proportions (how to do in BASIS as i only have proportions and not raw counts?)

# empty environment
rm(list=ls())

# set working directory to where the data is
setwd("~/Desktop/MTP_project")

#packages
library(ggplot2)
library(tidyverse)
library(readxl)

################################################################################
# Mutational Signatures
################################################################################

######### SCANB (ERpHER2n HER2E) #########

# load the raw signature counts
mut.scanb <- as.data.frame(t(read.table(file = 'Data/SCAN_B/WGS/Exposures_subs.tsv', sep = '\t', header = TRUE)))
#names <- colnames(mut.scanb)
names(mut.scanb) <- gsub('SBS', 'S', colnames(mut.scanb))

# present sigs: "SBS1"   "SBS2"   "SBS3"   "SBS5"   "SBS8"   "SBS13"  "SBS18"  "SBS127" "SBS126"

######### BASIS (ERpHER2n LumA; ERpHER2n LumA; ERpHER2p) #########

# I have to define 3 groups: ERpHER2n LumA; ERpHER2n LumA; ERpHER2p

# load the raw signature counts
load("Data/BASIS/Summarized_Annotations_BASIS.RData")

mut.basis <- as.data.frame(anno) %>% 
    select(c(sample_name,final.ER,final.HER2,
             ClinicalGroup,PAM50_AIMS,MutSigProportions)) %>%
    mutate(Subtype = 
               case_when(final.ER=="positive" & final.HER2=="positive" ~ "ERpHER2p",
                         final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumA" ~ "LumA",
                         final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumB" ~ "LumB")) %>% 
    filter(Subtype %in% c("ERpHER2p","LumA","LumB")) %>% 
    select(Subtype,MutSigProportions)
mut.basis$MutSigProportions <- as.data.frame(mut.basis$MutSigProportions)
# present sigs: "S1"  "S2"  "S3"  "S5"  "S6"  "S8"  "S13" "S17" "S18" "S20" "S26" "S30"
table(mut.basis$Subtype)
# check which signatures are present in both cohorts
common.mut.sigs <- intersect(colnames(mut.basis$MutSigProportions),colnames(mut.scanb))

# exclude non-shared ones
mut.scanb <- mut.scanb[,which(names(mut.scanb) %in% common.mut.sigs)]
mut.basis$MutSigProportions <- mut.basis$MutSigProportions[,which(names(mut.basis$MutSigProportions) %in% common.mut.sigs)]

# convert SCANB to proportion per sample
prop.mut.scanb <- as.data.frame(t(apply(t(mut.scanb),2,function(x) (x/sum(x, na.rm=TRUE))))) %>% 
    mutate(Subtype = "HER2E")

# get basis data in right format
mut.basis <- as.data.frame(do.call(cbind, mut.basis))
names(mut.basis) <- gsub("MutSigProportions.","", names(mut.basis))

# put together in an object ready for plotting 
mut.data <- rbind(mut.basis,prop.mut.scanb)
#View(mut.data)#
# plot
group.colors <- c(LumA = "#0f1bbd", LumB = "#09d3e6", HER2E ="#b609e6", ERpHER2p = "#e67b09")
pdf(file = "~/Desktop/MTP_project/Output/Plots/WGS/mut_signature_plots.pdf", onefile= TRUE)
for(i in c(1:7)) {
    plot <- ggplot(mut.data, aes(x=as.factor(Subtype),y=mut.data[,i+1],fill=as.factor(Subtype))) +
                geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
                xlab("Subtype") +
                ylab("Proportion") +
                ylim(c(0,1)) +
                ggtitle(paste("Mutational signature ",colnames(mut.data)[i+1],sep="")) +
                theme(plot.title = element_text(size = 30),
                      axis.text.x = element_text(size = 20),
                      axis.title.x = element_text(size = 25),
                      axis.text.y = element_text(size = 20),
                      axis.title.y = element_text(size = 25),
                      legend.position = "none") +
                scale_fill_manual(values=group.colors)
    
    print(plot)
}
dev.off()

################################################################################
# Statistics
################################################################################

mut.data$Subtype <- as.factor(mut.data$Subtype)

sig_ttest <- function(mut.data,sign){
    # vars
    hdata <- subset(mut.data, Subtype == "HER2E", select=c(sign)) 
    adata <- subset(mut.data, Subtype == "LumA", select=c(sign)) 
    bdata <- subset(mut.data, Subtype == "LumB", select=c(sign)) 
    edata <- subset(mut.data, Subtype == "ERpHER2p", select=c(sign)) 
    # for Her2 vs LumA
    # equal variance check
    if (var.test(unlist(hdata),unlist(adata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = FALSE)
    } else {
        H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = TRUE)
    }
    # save results
    print(paste("HER2E vs. LumA: ", H2vsLA_ttest_result$p.value,sep = "")) # 0.01478637
    
    # for Her2 vs LumB
    # equal variance check
    if (var.test(unlist(hdata),unlist(bdata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = FALSE)
    } else {
        H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = TRUE)
    }
    # save results
    print(paste("HER2E vs. LumB: ", H2vsLB_ttest_result$p.value,sep = ""))  # 0.07441933
    
    # for Her2 vs ERpHER2p
    # equal variance check
    if (var.test(unlist(hdata),unlist(edata), alternative = "two.sided")$p.value <= 0.05) {
        H2vsE_ttest_result <- t.test(hdata,edata, var.equal = FALSE)
    } else {
        H2vsE_ttest_result <- t.test(hdata,edata, var.equal = TRUE)
    }
    # save results
    print(paste("HER2E vs. ERpHER2p: ", H2vsE_ttest_result$p.value,sep = ""))  #0.01534438
    
}
sig_ttest(mut.data,"S1") #[1] "HER2E vs. LumA: * 
sig_ttest(mut.data,"S2") #
sig_ttest(mut.data,"S3") #[1] "HER2E vs. LumA: **  "HER2E vs. ERpHER2p: **
sig_ttest(mut.data,"S5") #[1] "HER2E vs. LumA: *** "HER2E vs. LumB: ***
sig_ttest(mut.data,"S8") #
sig_ttest(mut.data,"S13") #
sig_ttest(mut.data,"S18") #

################################################################################
# Rearrangement Signatures
################################################################################
rm(list=ls())

######### SCANB (ERpHER2n HER2E) #########

# load the raw signature counts
rearr.scanb <- as.data.frame(t(read.table(file = 'Data/SCAN_B/WGS/Exposures_rearr.tsv', sep = '\t', header = TRUE)))
#View(rearr.scanb)
#names <- colnames(mut.scanb)
names(rearr.scanb) <- gsub('RefSig R', 'RS', colnames(rearr.scanb))

# present sigs: "RefSig R2"  "RefSig R4"  "RefSig R6a" "RefSig R1"  "RefSig R5"  "RefSig R6b" "RefSig R3"  "RefSig R11"

######### BASIS (ERpHER2n LumA; ERpHER2n LumA; ERpHER2p) #########

# load the raw signature counts
load("Data/BASIS/Summarized_Annotations_BASIS.RData")
rearr.basis <- as.data.frame(anno) %>% 
    select(c(sample_name,final.ER,final.HER2,
             ClinicalGroup,PAM50_AIMS,RSproportions)) %>%
    mutate(Subtype = 
               case_when(final.ER=="positive" & final.HER2=="positive" ~ "ERpHER2p",
                         final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumA" ~ "LumA",
                         final.ER=="positive" & final.HER2=="negative" & PAM50_AIMS=="LumB" ~ "LumB")) %>% 
    filter(Subtype %in% c("ERpHER2p","LumA","LumB")) %>% 
    select(Subtype,RSproportions)
rearr.basis$RSproportions <- as.data.frame(rearr.basis$RSproportions)
# present sigs: "RS1" "RS2" "RS3" "RS4" "RS5" "RS6"

# FIND SOLUTION FOR RS6 vs RS6a RS6b

six.method <- "pool" # a; b; pool

if (six.method=="a") {
    rearr.scanb <- rearr.scanb %>% rename(RS6 = RS6a) %>% select(-c(RS6b))
} else if (six.method=="b") {
    rearr.scanb <- rearr.scanb %>% rename(RS6 = RS6b) %>% select(-c(RS6a))
} else if (six.method=="pool") {
    rearr.scanb <- rearr.scanb %>% mutate(RS6 = RS6a+RS6b) %>% select(-c(RS6a,RS6b))
} else { stop("check input") }

str(rearr.scanb)
# check which signatures are present in both cohorts
common.rearr.sigs <- intersect(colnames(rearr.basis$RSproportions),colnames(rearr.scanb))

# exclude non-shared ones
rearr.scanb <- rearr.scanb[,which(names(rearr.scanb) %in% common.rearr.sigs)]
rearr.basis$RSproportions <- rearr.basis$RSproportions[,which(names(rearr.basis$RSproportions) %in% common.rearr.sigs)]

# convert SCANB to proportion per sample
prop.rearr.scanb <- as.data.frame(t(apply(t(rearr.scanb),2,function(x) (x/sum(x, na.rm=TRUE))))) %>% 
    mutate(Subtype = "HER2E")

# get basis data in right format
rearr.basis <- as.data.frame(do.call(cbind, rearr.basis))
names(rearr.basis) <- gsub("RSproportions.","", names(rearr.basis))

# put together in an object ready for plotting 
rearr.data <- rbind(rearr.basis,prop.rearr.scanb)

# plot
group.colors <- c(LumA = "#0f1bbd", LumB = "#09d3e6", HER2E ="#b609e6", ERpHER2p = "#e67b09")
pdf(file = "~/Desktop/MTP_project/Output/Plots/WGS/rearr_signature_plots.pdf", onefile= TRUE)
for(i in c(1:6)) {
    plot <- ggplot(rearr.data, aes(x=as.factor(Subtype),y=rearr.data[,i+1],fill=as.factor(Subtype))) +
        geom_boxplot(alpha=0.7, size=1.5, outlier.size = 5) +
        xlab("Subtype") +
        ylab("Proportion") +
        ylim(c(0,1)) +
        ggtitle(paste("Rearrangement signature ",colnames(rearr.data)[i+1],sep="")) +
        theme(plot.title = element_text(size = 30),
              axis.text.x = element_text(size = 20),
              axis.title.x = element_text(size = 25),
              axis.text.y = element_text(size = 20),
              axis.title.y = element_text(size = 25),
              legend.position = "none") +
        scale_fill_manual(values=group.colors)
    
    print(plot)
}
dev.off()


################################################################################
# Statistics
################################################################################
sig_ttest(rearr.data,"RS1") 
sig_ttest(rearr.data,"RS2") 
sig_ttest(rearr.data,"RS3") 
sig_ttest(rearr.data,"RS4") #[1] "HER2E vs. LumA: **" lumb **** erpher2p ****
sig_ttest(rearr.data,"RS5") 
sig_ttest(rearr.data,"RS6") #[1] "HER2E vs. LumA: **




# for RS4

rearr.data$Subtype <- as.factor(rearr.data$Subtype)

# vars
hdata <- subset(rearr.data, Subtype == "HER2E", select=c(RS4)) 
adata <- subset(rearr.data, Subtype == "LumA", select=c(RS4)) 
bdata <- subset(rearr.data, Subtype == "LumB", select=c(RS4)) 
edata <- subset(rearr.data, Subtype == "ERpHER2p", select=c(RS4)) 
# for Her2 vs LumA
# equal variance check
if (var.test(unlist(hdata),unlist(adata), alternative = "two.sided")$p.value <= 0.05) {
    H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = FALSE)
} else {
    H2vsLA_ttest_result <- t.test(hdata,adata, var.equal = TRUE)
}
# save results
print(H2vsLA_ttest_result$p.value) # 0.001971769

# for Her2 vs LumB
# equal variance check
if (var.test(unlist(hdata),unlist(bdata), alternative = "two.sided")$p.value <= 0.05) {
    H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = FALSE)
} else {
    H2vsLB_ttest_result <- t.test(hdata,bdata, var.equal = TRUE)
}
# save results
print(H2vsLB_ttest_result$p.value) # 0.00003 . 3.333821e-05

# for Her2 vs ERpHER2p
# equal variance check
if (var.test(unlist(hdata),unlist(edata), alternative = "two.sided")$p.value <= 0.05) {
    H2vsE_ttest_result <- t.test(hdata,edata, var.equal = FALSE)
} else {
    H2vsE_ttest_result <- t.test(hdata,edata, var.equal = TRUE)
}
# save results
print(H2vsE_ttest_result$p.value) #4.62071e-09
