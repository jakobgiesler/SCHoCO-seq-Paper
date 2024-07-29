# Plotting single-cell community heatmap, analyze core microbiome community
# shared between diatom genotypes, test effects for each bacterial Genus separately.
#_______________________________________________________________________________

rm(list=ls())

library(tidyr)
library(ggplot2)
library(viridis)
library(SRS)
library(plyr)
library(phyloseq)
library(viridis)
library(dplyr)
library(vegan)
library(readr)
library(tibble)
library(gridExtra)
library(ggthemes)
library(zCompositions)
library(propr)
library(corpcor)
library(dendextend)
library(stringr)
library(gplots)
library(plotrix)

setwd('YOUR/DIRECTORY')

#_______________________________________________________________________________

#import files
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
count <- count_new*count_m
count$ASV_ID <- rownames(count)
count <- count[,c(145,1:144)]
count <- count[rowSums(count[,2:145])>0,]

#remove common contaminants
con <- read.delim('contaminants.txt')
int <- tax$ASV_ID[tax$Genus %in% intersect(tax$Genus,con$Contaminat_taxa)] #intersect and then remove the contaminats
count <- count[!count$ASV_ID %in% int, ]
count <- count[rowSums(count[,2:145])>0,]
tax <- tax[tax$ASV_ID %in% count$ASV_ID,]
count <- count[rownames(count) %in% tax$ASV_ID,]

#summarize on the genus level
tax <- tax[,1:7]#
asvs_tax <- merge(count,tax,by='ASV_ID') #bind to tax info
asvs_tax$ASV_ID <- NULL #remove ASV_ID
asvs_tax <- aggregate(.~Domain+Phylum+Class+Order+Family+Genus,data=asvs_tax,FUN=sum) #aggregate based on family
asvs_tax$id <- paste('new_ASV',seq(1,51,1),sep='_')#make new ID column

#remove chloroplasts, mitochondria, eukaryote ASVs
asvs_tax <- subset(asvs_tax, asvs_tax$Order!='Chloroplast')
asvs_tax <- subset(asvs_tax, asvs_tax$Family!='Mitochondria')
rownames(asvs_tax) <- asvs_tax[,151]
tax_backup <- asvs_tax[,1:6]
asvs_tax <- asvs_tax[,7:150]

#Scaling via ranked subsampling (SRS)
depths <- colSums(asvs_tax) #prepare df containing the sequencing depth
plot(depths[order(depths)])
abline(h=10000, col="blue")
abline(h=5000, col="red")#look for outliers of sequencing depth to cut off
quantile(depths) #we use 5000 as cut-off
quantile(depths, probs = seq(0,1,.1)) #Show all sequencing depth quantiles 

#cut off below 10% quantile
asvs_tax <- asvs_tax[,colSums(asvs_tax)>2304.0] #remove samples with depth < 10000 (none)

#running the SRS function
Cmin <- min(colSums(asvs_tax))
asvs_tax <- SRS(asvs_tax, Cmin)
plot(colSums(asvs_tax))

#power transformation
asvs_tax <- cbind(tax_backup,asvs_tax)
asvs_tax <- asvs_tax[rowSums(asvs_tax[,7:135])>0,]
tax_backup <- asvs_tax[,1:6]
asvs_tax <-asvs_tax[,7:135]
asvs_tax <- asvs_tax^0.25
asvs_tax <- cbind(tax_backup, asvs_tax)

#_______________________________________________________________________________

#Creating Heatmap

#create tidy data format table
df <-  gather(asvs_tax, sample_id, trans_counts, E.A1.MFa:P.A5.vit.h)

#add sample info
df <- merge(df, info, by = "sample_id")

#store list with desired genera
fam_list <- unique(df$Genus)

#filter the desired families from the original data frame to keep it tidy
df2 <- asvs_tax[asvs_tax$Genus %in% fam_list,]
df2 <-  gather(df2, sample_id, trans_counts, E.A1.MFa:P.A5.vit.h)
df2 <- merge(df2, info, by = "sample_id")
unique(df2$Genus)


#Plotting
p1 <- df2 %>% ggplot(aes(x=sample_id, y=Genus, fill=trans_counts))+
  geom_tile(colour="white", linewidth=0.25)+
  labs(x="", y="")+
  theme_grey(base_size=8)+
  facet_grid(Phylum~strain+condition+method,scales='free', space = "free")+
  scale_fill_viridis_c()+
  theme(
    legend.text=element_text(face="bold"),
    plot.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_blank(),
    #axis.text.x = element_text(angle = 90),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
p1
#ggsave("heatmap.eps",p1, device = "eps", units = "cm", width = 19, height = 15)

#Calculate some numbers that are mentioned in the manuscript results/discussion
df_for_text <- df2
df_for_text <- df_for_text[,2:11]
df_for_text <- dplyr::group_by(df_for_text, strain, Genus)
df_for_text <- dplyr::summarize(df_for_text, summedgenus = sum(trans_counts))

# Wich genera occurr in all strains?
species_in_all_sites <- df_for_text %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(sites_with_species = sum(summedgenus > 0)) %>%
  dplyr::filter(sites_with_species == length(unique(df_for_text$strain)))
print(species_in_all_sites)

#most abundant of the genera shared across all strains?
most_df <- df_for_text
most_df <- most_df[most_df$Genus %in% species_in_all_sites$Genus,]
most_df <- dplyr::ungroup(most_df)
most_df <- most_df %>% dplyr::group_by(Genus) %>% dplyr::mutate(summedstrain = sum(summedgenus))

#which ones are shared between strain A1 and A2?
shareA1A2 <- df_for_text
shareA1A2 <- shareA1A2 %>%
  pivot_wider(names_from = strain, values_from = summedgenus, values_fill = list(summedgenus = 0))

# Find genera present in strain A1 and A5, but absent in A2
species_in_A_and_B_not_C <- shareA1A2 %>%
  filter(A1 > 0, A5 > 0, A2 == 0) %>%
  nrow()

# Identify genera that are unique to a specific strain
unique_to_A1 <- shareA1A2 %>%
  filter(A1 > 0, A2 == 0, A5 == 0) 

# Identify genera that are unique to a specific strain
unique_to_A2 <- shareA1A2 %>%
  filter(A2 > 0, A1 == 0, A5 == 0) 

# Identify genera that are unique to a specific strain
unique_to_A5 <- shareA1A2 %>%
  filter(A5 > 0, A2 == 0, A1 == 0) 

#_______________________________________________________________________________

#Testing for differences for each bacterial Genus seperately using multiple ANOVAs

#delete prior groupings
df2 <- ungroup(df2)
df2 <- as.data.frame(df2)
is.na(df2)
df2[is.na(df2)] <- 0
count(is.na(df2)==TRUE)

#define as factors
df2$method <- as.factor(df2$method)
df2$strain <- as.factor(df2$strain)
df2$condition <- as.factor(df2$condition)
df_fams <- split(df2, df2$Genus) #create list with subsetted genera
str(df_fams)

#apply 
anova_results_list <- lapply(df_fams, function(x){as.data.frame(anova(aov(trans_counts~strain*condition,data=x)))[["Pr(>F)"]]})

#extract p-values for each effect for each genus and store as dataframe
anova_results_df <- do.call(rbind, anova_results_list)
colnames(anova_results_df) <- c('strain','condition','strain:condition','residuals')
anova_results_df <- data.frame(anova_results_df[,1:4]) 

#write function to replace p-values with asterisks
replace_values <- function(value) {
  if (value < 0.001) {
    return('***')
  } else if (value > 0.001 && value < 0.01) {
    return('**')
  } else if (value > 0.01 && value < 0.05) {
    return('*')
  } else {
    return(NA)
  }
}

#exclude residuals
anova_results_df <- dplyr::select(anova_results_df, strain, condition, strain.condition)

#Apply the function to each element of the data frame
replaced_df <- apply(anova_results_df, c(1, 2), replace_values)

# Convert the result back to a data frame
replaced_df <- as.data.frame(replaced_df)
