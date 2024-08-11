### load library ####
library(ggplot2)

### load files ####
data <- read.delim('final_count_data.txt')
info <- read.delim('final_info_data.txt')
tax <- read.delim('final_tax_data.txt')


### subset data (per strain) ####
rownames(data) <- data$ASV_ID
data$ASV_ID <- NULL
a1 <- data[,grepl('_A1_',colnames(data))]
a2 <- data[,grepl('_A2_',colnames(data))]
h5 <- data[,grepl('_H5_',colnames(data))]

### remove zero count only ASVs ####
a1 <- a1[rowSums(a1)>0,]
a2 <- a2[rowSums(a2)>0,]
h5 <- h5[rowSums(h5)>0,]

### subset metadata file ####
info_a1 <- info[info$newName %in% colnames(a1),]
info_a2 <- info[info$newName %in% colnames(a2),]
info_h5 <- info[info$newName %in% colnames(h5),]

### add taxonomy to ASV table and group ASVs for each strain data set ####
a1$ASV_ID <- rownames(a1)
a1 <- merge(a1,tax,by='ASV_ID',all.x = T)
a1[is.na(a1)] <- 'NULL'
a1$ID <- ifelse(a1$Order=='Chloroplast',"CHL",ifelse(a1$Family=='Mitochondria'| is.na(a1$Order),"OTHER","ASV"))

a2$ASV_ID <- rownames(a2)
a2 <- merge(a2,tax,by='ASV_ID',all.x = T)
a2[is.na(a2)] <- 'NULL'
a2$ID <- ifelse(a2$Order=='Chloroplast',"CHL",ifelse(a2$Family=='Mitochondria'| is.na(a2$Order),"OTHER","ASV"))

h5$ASV_ID <- rownames(h5)
h5 <- merge(h5,tax,by='ASV_ID',all.x = T)
h5[is.na(h5)] <- 'NULL'
h5$ID <- ifelse(h5$Order=='Chloroplast',"CHL",ifelse(h5$Family=='Mitochondria'| is.na(h5$Order),"OTHER","ASV"))

### select bacterial ASVs for each strain data set and plot richness (CAS vs. noCAS) ####
a1_asv <- subset(a1,a1$ID %in% 'ASV')
a1_asv <- a1_asv[,c(1:45)]
a1_asv[,c(2:45)][a1_asv[,c(2:45)] > 0] <- 1
a1_asv$cas9 <- rowSums(a1_asv[,grepl('_Cas_',colnames(a1_asv))])
a1_asv$no_cas9 <- rowSums(a1_asv[,grepl('_NoCas_',colnames(a1_asv))])
a1_asv <- merge(a1_asv,tax,by='ASV_ID',all.x = T)
a1_rich <- as.data.frame(cbind(colnames(a1_asv)[2:45],colSums(a1_asv[2:45])),)
a1_rich$V2 <- as.numeric(a1_rich$V2)
a1x <- a1_rich$V1[grepl('_Cas_',a1_rich$V1)]
a1_rich$V3 <- ifelse(a1_rich$V1 %in% a1x, 'xcas9', 'no_cas9')
a1_rich$V4 <- info_a1$pair[match(a1_rich$V1,info_a1$newName)]
a1_rich$V5 <- 'A1'

a2_asv <- subset(a2,a2$ID=='ASV')
a2_asv <- a2_asv[,c(1:49)]
a2_asv[,c(2:49)][a2_asv[,c(2:49)] > 0] <- 1
a2_asv$cas9 <- rowSums(a2_asv[,grepl('_Cas_',colnames(a2_asv))])
a2_asv$no_cas9 <- rowSums(a2_asv[,grepl('_NoCas_',colnames(a2_asv))])
a2_asv <- merge(a2_asv,tax,by='ASV_ID',all.x = T)
a2_rich <- as.data.frame(cbind(colnames(a2_asv)[2:49],colSums(a2_asv[2:49])),)
a2_rich$V2 <- as.numeric(a2_rich$V2)
a2x <- a2_rich$V1[grepl('_Cas_',a2_rich$V1)]
a2_rich$V3 <- ifelse(a2_rich$V1 %in% a2x, 'xcas9', 'no_cas9')
a2_rich$V4 <- info_a2$pair[match(a2_rich$V1,info_a2$newName)]
a2_rich$V5 <- 'A2'

h5_asv <- subset(h5,h5$ID=='ASV')
h5_asv <- h5_asv[,c(1:33)]
h5_asv[,c(2:33)][h5_asv[,c(2:33)] > 0] <- 1
h5_asv$cas9 <- rowSums(h5_asv[,grepl('_Cas_',colnames(h5_asv))])
h5_asv$no_cas9 <- rowSums(h5_asv[,grepl('_NoCas_',colnames(h5_asv))])
h5_asv <- merge(h5_asv,tax,by='ASV_ID',all.x = T)
h5_rich <- as.data.frame(cbind(colnames(h5_asv)[2:33],colSums(h5_asv[2:33])),)
h5_rich$V2 <- as.numeric(h5_rich$V2)
h5x <- h5_rich$V1[grepl('_Cas_',h5_rich$V1)]
h5_rich$V3 <- ifelse(h5_rich$V1 %in% h5x, 'xcas9', 'no_cas9')
h5_rich$V4 <- info_h5$pair[match(h5_rich$V1,info_h5$newName)]
h5_rich$V5 <- 'H5'

### combined data and richness plot ####
all_rich <- as.data.frame(rbind(a1_rich,a2_rich,h5_rich))
ggplot(all_rich, aes(V3,V2,shape=V5,col=V5))+geom_boxplot(outlier.shape = NA,col="black")+geom_jitter(width = 0.2,height = 0.2)+geom_line(mapping = aes(group = V4),col='blue')+facet_grid(.~V5,scales='free')+theme_light()


### calculate the indicene rate ratios per strain (at the family level) ####

## prepare data for strain a1 ####
fa1 <- a1[,c(2:45,50)]
fa1 <- aggregate(.~Family,fa1,FUN=sum)
fa1l <- split(fa1,fa1$Family)
a1ld <- lapply(fa1l, function(x) {as.data.frame(cbind(counts <- unlist(c(x[colnames(x) %in% colnames(x)[grepl('_Cas_', colnames(x))]],x[colnames(x) %in% colnames(x)[grepl('_NoCas_', colnames(x))]])),group <- rep(c("xCas", "NoCas"), each = 22)))})
a1ld <- lapply(a1ld,function(x){x$V1 <- as.numeric(x$V1); return(x)})
info_mean <- lapply(a1ld, function(x){aggregate(V1~V2,x,FUN=sum)})
info_mean <- Map(cbind, info_mean, SampleID = names(info_mean))
info_mean <- do.call(rbind,info_mean)

## Fit Poisson regression model for families (strain a1) and extract results ####
list_poisson_model <- lapply(a1ld, function(x){glm(V1 ~ V2,data=x, family = poisson(link="log"))})
list_pvalues <- lapply(list_poisson_model,function(x){summary(x)[["coefficients"]][2,4]})
list_pvalues <- as.data.frame(do.call(rbind,list_pvalues))
info_mean$p_poi <- list_pvalues$V1[match(info_mean$SampleID,rownames(list_pvalues))]
info_mean$sign <- ifelse(info_mean$p_poi<0.05,'yes','no')
# Extract coefficients for the grouping variable ####
l_IRR <- lapply(list_poisson_model, function(x){exp(coef(x)["V2xCas"])})
l_coef <- lapply(list_poisson_model, function(x){coef(x)["V2xCas"]})
L_sde <- lapply(list_poisson_model, function(x){summary(x)$coefficients[2, "Std. Error"]})
## combine data for visualization ####
IRR <- as.data.frame(do.call(rbind,l_IRR))
coef <- as.data.frame(do.call(rbind,l_coef))
sde <- as.data.frame(do.call(rbind,L_sde))
info_mean$IRR <- IRR$V2xCas[match(info_mean$SampleID,rownames(IRR))]
info_mean$coef <- coef$V2xCas[match(info_mean$SampleID,rownames(coef))]
info_mean$sde <- sde$V1[match(info_mean$SampleID,rownames(sde))]
# Calculate confidence intervals of IRRs ####
info_mean$CI_lower <- exp(info_mean$coef-1.96*info_mean$sde)
info_mean$CI_upper <- exp(info_mean$coef+1.96*info_mean$sde)
# store final data frame ####
a1_irr_family_all <- info_mean
a1_irr_family_all$strain <- 'a1'

## prepare data for strain a2 ####
fa2 <- a2[,c(2:49,54)]
fa2 <- aggregate(.~Family,fa2,FUN=sum)
fa2l <- split(fa2,fa2$Family)
a2ld <- lapply(fa2l, function(x) {as.data.frame(cbind(counts <- unlist(c(x[colnames(x) %in% colnames(x)[grepl('_Cas_', colnames(x))]],x[colnames(x) %in% colnames(x)[grepl('_NoCas_', colnames(x))]])),group <- rep(c("xCas", "NoCas"), each = 24)))})
a2ld <- lapply(a2ld, function(x) {x$V1 <- as.numeric(x$V1); return(x)})
info_mean <- lapply(a2ld, function(x){aggregate(V1~V2,x,FUN=sum)})
info_mean <- Map(cbind, info_mean, SampleID = names(info_mean))
info_mean <- do.call(rbind,info_mean)
## Fit Poisson regression model for families (strain a2) and extract results ####
list_poisson_model <- lapply(a2ld, function(x){glm(V1 ~ V2,data=x, family = poisson(link='log'))})
list_pvalues <- lapply(list_poisson_model,function(x){summary(x)[["coefficients"]][2,4]})
list_pvalues <- as.data.frame(do.call(rbind,list_pvalues))
info_mean$p_poi <- list_pvalues$V1[match(info_mean$SampleID,rownames(list_pvalues))]
info_mean$sign <- ifelse(info_mean$p_poi<0.05,'yes','no')
# Extract coefficients for the grouping variable ####
l_IRR <- lapply(list_poisson_model, function(x){exp(coef(x)["V2xCas"])})
l_coef <- lapply(list_poisson_model, function(x){coef(x)["V2xCas"]})
L_sde <- lapply(list_poisson_model, function(x){summary(x)$coefficients[2, "Std. Error"]})
## combine data for visualization ####
IRR <- as.data.frame(do.call(rbind,l_IRR))
coef <- as.data.frame(do.call(rbind,l_coef))
sde <- as.data.frame(do.call(rbind,L_sde))
info_mean$IRR <- IRR$V2xCas[match(info_mean$SampleID,rownames(IRR))]
info_mean$coef <- coef$V2xCas[match(info_mean$SampleID,rownames(coef))]
info_mean$sde <- sde$V1[match(info_mean$SampleID,rownames(sde))]
# Calculate confidence intervals of IRRs ####
info_mean$CI_lower <- exp(info_mean$coef-1.96*info_mean$sde)
info_mean$CI_upper <- exp(info_mean$coef+1.96*info_mean$sde)
# store final data frame ####
a2_irr_family_all <- info_mean
a2_irr_family_all$strain <- 'a2'

## prepare data for strain h5 ####
fh5 <- h5[,c(2:33,38)]
fh5 <- aggregate(.~Family,fh5,FUN=sum)
fh5l <- split(fh5,fh5$Family)
h5ld <- lapply(fh5l, function(x) {as.data.frame(cbind(counts <- unlist(c(x[colnames(x) %in% colnames(x)[grepl('_Cas_', colnames(x))]],x[colnames(x) %in% colnames(x)[grepl('_NoCas_', colnames(x))]])),group <- rep(c("xCas", "NoCas"), each = 16)))})
h5ld <- lapply(h5ld, function(x){x$V1 <- as.numeric(x$V1); return(x)})
info_mean <- lapply(h5ld, function(x){aggregate(V1~V2,x,FUN=sum)})
info_mean <- Map(cbind, info_mean, SampleID = names(info_mean))
info_mean <- do.call(rbind,info_mean)
## Fit Poisson regression model for families (strain h5) and extract results ####
list_poisson_model <- lapply(h5ld, function(x){glm(V1 ~ V2,data=x, family = poisson(link='log'))})
list_pvalues <- lapply(list_poisson_model,function(x){summary(x)[["coefficients"]][2,4]})
list_pvalues <- as.data.frame(do.call(rbind,list_pvalues))
info_mean$p_poi <- list_pvalues$V1[match(info_mean$SampleID,rownames(list_pvalues))]
info_mean$sign <- ifelse(info_mean$p_poi<0.05,'yes','no')
# Extract coefficients for the grouping variable ####
l_IRR <- lapply(list_poisson_model, function(x){exp(coef(x)["V2xCas"])})
l_coef <- lapply(list_poisson_model, function(x){coef(x)["V2xCas"]})
L_sde <- lapply(list_poisson_model, function(x){summary(x)$coefficients[2, "Std. Error"]})
## combine data for visualization ####
IRR <- as.data.frame(do.call(rbind,l_IRR))
coef <- as.data.frame(do.call(rbind,l_coef))
sde <- as.data.frame(do.call(rbind,L_sde))
info_mean$IRR <- IRR$V2xCas[match(info_mean$SampleID,rownames(IRR))]
info_mean$coef <- coef$V2xCas[match(info_mean$SampleID,rownames(coef))]
info_mean$sde <- sde$V1[match(info_mean$SampleID,rownames(sde))]
# Calculate confidence intervals of IRRs ####
info_mean$CI_lower <- exp(info_mean$coef-1.96*info_mean$sde)
info_mean$CI_upper <- exp(info_mean$coef+1.96*info_mean$sde)
# store final data frame ####
h5_irr_family_all <- info_mean
h5_irr_family_all$strain <- 'h5'

### combine data for final graph ####
final <- rbind(h5_irr_family_all,a1_irr_family_all,a2_irr_family_all)
final <- subset(final, final$sign=='yes')
final <- final[!final$SampleID %in% c('Mitochondria','NULL'),]
final$Phylum <- tax$Phylum[match(final$SampleID,tax$Family)]

ggplot(final, aes(y = SampleID, x = IRR)) +
  geom_point(size = 2, aes(shape = strain)) +
  geom_errorbar(aes(xmin = CI_lower,xmax = CI_upper, color = strain), width = 0.2) + scale_color_viridis_d(option='turbo',direction=-1)+geom_vline(xintercept=1)+
  labs(y = "Family", x = "Incidence Rate Ratio (IRR)", title = "Incidence Rate Ratios with Confidence Intervals") +
  theme_light()+facet_grid(Phylum~.,scales="free_y", space = 'free') 

#### End ####