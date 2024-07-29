# Plotting and testing diatom single-cell richness and shannon index (SC-seq)
# as a response to strain identity (A1, A2, A5), culture medium condition (full medium,
# nitrate limited, vitamin limited), DNA processing method (extracted DNA vs. direct PCR).
#_______________________________________________________________________________

rm(list = ls())

library(tidyr)
library(ggplot2)
library(viridis)
library(SRS)
library(plyr)
library(phyloseq)
library(vegan)
library(readr)
library(tibble)
library(gridExtra)
library(ggthemes)
library(zCompositions)
library(propr)
library(corpcor)
library(dendextend)
library(microbiome)

setwd("YOUR/DIRECTORY")

options(getClass.msg=FALSE)

#_______________________________________________________________________________

tax <- read.delim('exp_pcr_taxa.txt')
info <- read.delim('exp_pcr_info.txt')
count <- read.delim('exp_pcr_counts.txt') 
count[is.na(count)] <- 0

tax$Genus[tax$ASV_ID == "ASV_363"] <- "Yoonia-Loktanella"
tax$Genus[tax$ASV_ID == "ASV_10"] <- "Roseobacter clade NAC11-7 lineage"
tax$Genus[tax$ASV_ID == "ASV_431"] <- "Klebsiella"
tax$Genus[tax$ASV_ID == "ASV_15"] <- "Colwellia"

#cut-off as recommended in the paper
count_new <- count
rownames(count_new) <- count_new$ASV_ID
count_new$ASV_ID <- NULL
head(count_new)
count_prop <- apply(count_new,2,function(x){(100/sum(x))*x})
colSums(count_prop)
count_prop[count_prop<=0.25]<-0
count_m <- count_prop
count_m[count_m>0]<-1

#some data wrangling
count <- count_new*count_m
count$ASV_ID <- rownames(count)
count <- count[,c(145,1:144)]
count <- count[rowSums(count[,2:145])>0,]

#remove common contaminants
con <- read.delim('contaminants.txt') #import list with common contaminates
int <- tax$ASV_ID[tax$Genus %in% intersect(tax$Genus,con$Contaminat_taxa)] #intersect 
count <- count[!count$ASV_ID %in% int, ] #remove
count <- count[rowSums(count[,2:145])>0,]
head(count)
head(tax)

# remove chloroplasts and mitochondria
tax_new <- tax[tax$ASV_ID %in% rownames(count),]
tax_new[is.na(tax_new)] <- 'x1'
tax_new <- subset(tax_new, tax_new$Order!='Chloroplast')
tax_new <- subset(tax_new, tax_new$Family!='Mitochondria')
count <- count[rownames(count) %in% tax_new$ASV_ID,]

#Cutoff at low 10% quantile or at 10.000 reads sequencing depth
depths <- colSums(count[,2:145]) #prepare df containing the sequencing depth
plot(depths[order(depths)])
abline(h=5000, col="blue")
abline(h=10000, col="red")#look for outliers of sequencing depth to cut off
quantile(depths) #we use 10000 as cut-off
quantile(depths, probs = seq(0,1,.1)) #Show all sequencing depth quantiles 
rownames(count) <- count$ASV_ID
count <- count[,2:145]
count <- count[,colSums(count)>=2337] #cut off below 10% quantile or 10.000
info_new <- info[info$sample_id %in% colnames(count),]

#Scaling via ranked subsampling
Cmin <- min(colSums(count))
df <- SRS(count, Cmin)
count[,1:129] <- df

#Normalization with power transformation 
count_new <- count^0.25
count_new <- count_new[!rowSums(count_new)==0,]

#Preparations for phyloseq object
rownames(tax_new) <- tax_new$ASV_ID #set ASVs as rownames
tax_new <- tax_new[rownames(tax_new) %in% rownames(count_new),]
tax_new$ASV_ID <- NULL
tax_new <- as.matrix(tax_new)

rownames(info_new) <- info_new$sample_id 
info_new <- info_new[rownames(info_new) %in% colnames(count_new),]
rownames(info_new) <- info_new$sample_id #set sample IDs as row names
#info_new$sample_id <- NULL #delete "old" Sample ID column

# Check and adjust that rownames of samples_df match colnames of otu_mat
all(rownames(info_new) == colnames(count_new)) #should be TRUE

# name elements for subsequent phyloseq object
OTU = otu_table(count_new, taxa_are_rows = TRUE)
TAX = tax_table(tax_new)
samples = sample_data(info_new)

# create phyloseq object
phy <- phyloseq(OTU, TAX, samples)
phy

#_______________________________________________________________________________

#Calculate Diversity Measures
rich <- richness(phy)
shann <- diversity(phy, index = "shannon")
div <- cbind(info_new, rich, shann)
div$chao1 <- NULL
#Richness

p1 <- div %>% ggplot(aes(method, observed))+
  geom_boxplot(width=0.3, outlier.shape = NA, fill = "lightgrey")+
  geom_jitter(width = 0.1)+
  ylab("Richness")+
  theme_few()
p1
#ggsave("rich_methods_2024.eps",p1, device = "eps", units = "cm", width = 10, height = 10, dpi = 300)

p2 <- div %>% ggplot(aes(strain, observed))+
  geom_boxplot(width=0.3, outlier.shape = NA, fill = "lightgrey")+
  geom_jitter(width = 0.1)+
  ylab("Richness")+
  facet_grid(~condition)+
  theme_light()
p2
#ggsave("rich_con_2024.eps",p2, device = "eps", units = "cm", width = 14, height = 6, dpi = 300)

#anova model (used in paper)
m1 <- aov(observed~strain*condition+method, data = div)
#plot(m1)
summary(m1)
TukeyHSD(m1)

#Shannon
p3 <- div %>% ggplot(aes(method, shannon))+
  geom_boxplot(width=0.3, outlier.shape = NA, fill="lightgrey")+
  geom_jitter(width = 0.1)+
  ylab("Shannon Diversity")+
  theme_few()
p3
#ggsave("shann_methods_2024.eps",p3, device = "eps", units = "cm", width = 10, height = 10, dpi = 300)

p4 <- div %>% ggplot(aes(strain, shannon))+
  geom_boxplot(width=0.3, outlier.shape = NA, fill="lightgrey")+
  geom_jitter(width = 0.1)+
  facet_grid(~condition)+
  ylab("Shannon Diversity")+
  theme_light()
p4
#ggsave("shann_con_2024.eps",p4, device = "eps", units = "cm", width = 14, height = 6, dpi = 300)

#anova model
m2 <- aov(shannon~strain*condition+method, data = div)
#plot(m2)
summary(m2)
TukeyHSD(m2)

