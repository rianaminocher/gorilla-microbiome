# initial OTU analysis 

rm(list=ls())

library(phyloseq)
library(ggplot2)
library(vegan)
library(reshape2)
library(ape)
library(plyr)
library(dplyr)

setwd("/Volumes/riana_minocher/gorilla-microbiome-gh/") # change wd to the directory with data & scripts

# setwd("/Users/riana/ownCloud/riana_minocher/gorilla-microbiome/")

      
      ### load data raw ###

otu_raw <- import_biom(file.path("data/from_uc_w_tax.biom"))
tree <- read_tree_greengenes("data/rep_set.tre")
taxa_names(tree) <- gsub("><-><", ";", taxa_names(tree)) # change the taxa names to match otu tab
load("data/meta.robj")
otu_raw <- merge_phyloseq(meta, otu_raw, tree)



      ### sample summary ###

# populations
nrow(meta[meta$Population == "KB", ])
nrow(meta[meta$Population == "NB", ])
length(unique(meta[meta$Population == "KB", ]$Composite_Indidvidual_ID))
length(unique(meta[meta$Population == "NB", ]$Composite_Indidvidual_ID))
unique(meta[duplicated(meta$Composite_Indidvidual_ID), ]$Composite_Indidvidual_ID)

# sampling date
setdiff(meta$Date.of.coll., meta$Collection_Date)
length(grep("/06/|/07/", meta$Collection_Date)) # 95 samples sampled in June & July 
meta[!grepl("/06/|/07/", meta$Collection_Date), ]$Social_Group_G # only chims not collected in june july
levels(meta[!grepl("/06/|/07/", meta$Collection_Date), ]$Collection_Date)

# repeated samples
rep_samples <- unique(meta[duplicated(meta$Composite_Indidvidual_ID), ]$Composite_Indidvidual_ID)[2:24]
sampling_dates <- lapply(rep_samples,
                         function(x) {
                           levels(meta[meta$Composite_Indidvidual_ID == x]$Collection_Date)
                           })
sum(lengths(sampling_dates) == 1)
sum(lengths(sampling_dates) == 2)
names(sampling_dates) <- rep_samples

# no. of unique individuals per group
unique_inds <- as.data.frame(levels(meta$Social_Group_G))
unique_inds$unique_inds <- c(length(unique(meta[meta$Social_Group_G == "Chimanuka"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "Mpungwe"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "Mufanzala 1"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "Mufanzala2"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "Mugaruka"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "Namadiriri"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "NB"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "Nouvelle Famille"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "Part of Mankoto"]$Composite_Indidvidual_ID)) - 1,
                           length(unique(meta[meta$Social_Group_G == "Unknown"]$Composite_Indidvidual_ID)) - 1
                           )



      ### taxonomy summary ###

tax_table <- as.data.frame(otu_raw@tax_table)
levels(tax_table$Rank1)
(nrow(tax_table[grep("Eukaryota|Archaea", tax_table$Rank1), ]) / 
    nrow(tax_table)) * 100 # 0.36% Archaea
(nrow(tax_table[grep("Unclassified", tax_table$Rank1), ]) / 
    nrow(tax_table)) * 100 # 0.01% Unclassified



      ### read count summary ###

total_reads <- sum(meta$Read_Count)
total_reads_quality <- sum(meta$Read_Count_quality)
pct_reads_discarded <- (sum(meta$Read_Count) - sum(meta$Read_Count_quality)) / sum(meta$Read_Count) * 100
summary(meta[grep("Blank", meta$Sample_ID), ]$Read_Count_quality)
summary(meta[grep("KB|NB", meta$Sample_ID), ]$Read_Count_quality)
sort(meta[grep("KB|NB", meta$Sample_ID), ]$Read_Count_quality, decreasing = TRUE)



      ### load data bacteria only ###

otu_bact <- import_biom(file.path("data/from_uc_w_tax_bact.biom"))
otu_bact <- merge_phyloseq(meta, otu_bact, tree)


      ### read depth summary ###

head(sort(sample_sums(otu_bact)))
tail(sort(sample_sums(otu_bact)))
readdepth <- data.frame(readdepth = sample_sums(otu_bact))
readdepth$sampletype <- "sample"
readdepth[grepl("Blank",rownames(readdepth)), ]$sampletype <- "blank"
p_readdepth <- ggplot(readdepth, aes(x = readdepth, fill = factor(sampletype)))
p_readdepth <- p_readdepth + geom_histogram(color = "black", binwidth = 3000) + theme_bw()
p_readdepth <- p_readdepth + xlab("Sequencing read depth") + ylab("")
p_readdepth <- p_readdepth + scale_fill_discrete(name = "Sample type", labels = c("Extraction blank", "Faecal sample"))
p_readdepth <- p_readdepth + scale_x_continuous(breaks = pretty_breaks(10))
p_readdepth



      ### analysis of extraction blanks ###

# check distribution across samples 

otu_tab_bact <- as.matrix(otu_bact@otu_table)      # make otu tab
otu_tab_bact <- prop.table(otu_tab_bact, 2)        # transform counts to proportions
colSums(otu_tab_bact)                              # check that all proportions add to 1
blank_otu_tab_bact <- otu_tab_bact[ , grepl("Blank", colnames(otu_tab_bact))]    # subset only blanks
blank_otu_tab_bact <- blank_otu_tab_bact[!rowSums(blank_otu_tab_bact) == 0, ]    # select rows(OTUs) for which any blank has an observation
blank_otu_tab_bact_all <- merge(blank_otu_tab_bact,
                                (otu_tab_bact[ , !grepl("Blank", colnames(otu_tab_bact))]), by = 0) # merge with samples
rownames(blank_otu_tab_bact_all) <- blank_otu_tab_bact_all$Row.names
blank_otu_tab_bact_all <- blank_otu_tab_bact_all[ , -1]   # dataframe of all "Blank" OTUs and their presence in samples
colSums(blank_otu_tab_bact_all) # now only blanks have complete proportions
summary(colSums(blank_otu_tab_bact_all[ , !grepl("Blank", colnames(blank_otu_tab_bact_all))]))* 100 # rep. on average 64.9% of all sample OTUs
hist(colSums(blank_otu_tab_bact_all))

# the mean proportional presence of all "Blank" OTUs in samples is 65%

blank_otu_tab_bact[blank_otu_tab_bact > 0] <- 1 # make it binary
max(rowSums(blank_otu_tab_bact)) # max otus in all blanks
table(rowSums(blank_otu_tab_bact))
which(rowSums(blank_otu_tab_bact) == 5)
bin_blank_otu_tab_bact_all <- as.data.frame(blank_otu_tab_bact_all)
bin_blank_otu_tab_bact_all[bin_blank_otu_tab_bact_all > 0] <- 1
rowSums(bin_blank_otu_tab_bact_all)
hist(rowSums(bin_blank_otu_tab_bact_all)) # blank OTUs present in 100-120 samples mostly
table(rowSums(bin_blank_otu_tab_bact_all))

# make a heatplot of the OTUs present in more than 3 of 9 blanks, to see proportional abundance in samples

blankotus_phyl <- transform_sample_counts(subset_samples(otu_bact, !sample_names(otu_bact) 
                                                         %in% meta[meta$Read_Count_quality < 50000, ]$Sample_ID),
                                          fun = function(x) {x / sum(x)})
otus_to_keep <- names(which(rowSums(blank_otu_tab_bact) > 3))
blankotus_phyl <- prune_taxa(otus_to_keep, blankotus_phyl)
mean(colSums(blankotus_phyl@otu_table))
ntaxa(blankotus_phyl) # 20 OTUs present in more than 3 blanks
ntaxa(otu_bact) # 23684 total OTUs
heatmap <- plot_heatmap(blankotus_phyl)
heatmap

# OTUs with high average abundance across samples

sample_otu_tab_bact <- otu_tab_bact[ , !grepl("Blank",colnames(otu_tab_bact))]# otu tab samples
mean_abund_sample_otus <- rowMeans(sample_otu_tab_bact) * 100
max(mean_abund_sample_otus) # maximum is 2.45%
ggplot(melt(mean_abund_sample_otus), aes(x = value)) + 
  geom_histogram(color = "black", fill = "goldenrod3", bins = 100) +
  ggtitle("Distribution of mean OTU abundances") + 
  xlab("Mean proportional abundance across samples") +
  theme(axis.title.y = element_blank())
mean_abund_sample_otus <- sort(mean_abund_sample_otus, decreasing = TRUE) # order by abundance
top30 <- mean_abund_sample_otus[1 : (23684 / 3)]
setdiff(rownames(blank_otu_tab_bact), names(top30))
length(which(mean_abund_sample_otus > 1)) # 9 OTUs have a mean proportional abundance > 1%
length(which(mean_abund_sample_otus > 0.5)) # 42 OTUs have an abundance greater than 0.5%
length(which(mean_abund_sample_otus > 0.2)) # 102 OTUs have an abundance greater than 0.2%
setdiff(names(mean_abund_sample_otus[mean_abund_sample_otus > 1]), rownames(blank_otu_tab_bact))
setdiff(names(mean_abund_sample_otus[mean_abund_sample_otus > 0.5]), rownames(blank_otu_tab_bact))
setdiff(names(mean_abund_sample_otus[mean_abund_sample_otus > 0.2]), rownames(blank_otu_tab_bact))

## All but 9 of 356 "Blank OTUs" are present in the top 1/3 most abundant sample OTUs
## All OTUs with a mean abundance greater than 1% are present in the "Blank OTUs"
## All but 1 of 42 OTUs with a mean abundance greater than 0.5% are present in the "Blank OTUs"
## 9 of 102 OTUs with a mean abundance greater than 0.2% are present in the "Blank OTUs"

# Blank OTUs from vsearch output show similar results to QIIME output



      ### taxonomy summary  ###


# make the otu table with only samples
otu_bact_samples <- subset_samples(otu_bact, !sample_names(otu_bact) 
                                   %in% meta[meta$Read_Count_quality < 50000, ]$Sample_ID) # by selecting readcount > 50,000
str(otu_bact_samples@tax_table)
tax_table_bact_samples <- as.data.frame(otu_bact_samples@tax_table)
tax_table_bact_samples$Rank1 <- gsub("D_0__","",tax_table_bact_samples$Rank1)
tax_table_bact_samples$Rank2 <- gsub("D_1__","",tax_table_bact_samples$Rank2)
tax_table_bact_samples$Rank3 <- gsub("D_2__","",tax_table_bact_samples$Rank3)
tax_table_bact_samples$Rank4 <- gsub("D_3__","",tax_table_bact_samples$Rank4)
tax_table_bact_samples$Rank5 <- gsub("D_4__","",tax_table_bact_samples$Rank5)
tax_table_bact_samples$Rank6 <- gsub("D_5__","",tax_table_bact_samples$Rank6)
tax_table_bact_samples$Rank7 <- gsub("D_6__","",tax_table_bact_samples$Rank7)

# bacterial phyla
bact_phyla <- read.delim("data/map_R_L2.txt", sep = "\t", h = TRUE)      # read bacterial phyla tab
rownames(bact_phyla) <- bact_phyla$X.SampleID
bact_phyla <- select(bact_phyla, grep("Bacteria", colnames(bact_phyla))) # select only bacteria data (remove meta)
bact_phyla <- as.data.frame(t(bact_phyla)) # transposeee for plot

bact_phyla <- bact_phyla[, !grepl("Blank", colnames(bact_phyla))]        # get rid of blanks
bact_phyla <- bact_phyla[ , !grepl("KBMu21611|KBNa1921", colnames(bact_phyla))] 
colSums(bact_phyla)
avg_bact_phyla_comp <- sort(rowMeans(bact_phyla) * 100, decreasing = TRUE)


bact_low <- colSums(subset(bact_phyla, (rowMeans(bact_phyla) < 0.01) &     # bact present in < 0.01% mean proportion across samples
                     (!grepl("Other", rownames(bact_phyla)))))             # but not the unclassified bact
bact_high <- subset(bact_phyla, rowMeans(bact_phyla) >= 0.01)              # bact present in > 0.01% mean proportion
bact_other <- subset(bact_phyla, grepl("Other", rownames(bact_phyla)))     # group all unclassified 
bact_phyla_summary <- as.matrix(rbind(bact_high, bact_low, bact_other) * 100) # stitch together
rownames(bact_phyla_summary) <- gsub("D_0__Bacteria.D_1__|D_0__Bacteria.", "", rownames(bact_phyla_summary))
rownames(bact_phyla_summary)[rownames(bact_phyla_summary) == "10"] <- "Bacterial phyla < 1%"
rownames(bact_phyla_summary)[rownames(bact_phyla_summary) == "Other"] <- "Unclassified phyla"
colSums(bact_phyla_summary)

# phylum stacked barplot

meta$Social_Group <- factor(meta$Social_Group_G, levels = c("Ch","Mp","Mu1","Mu2","Mg","NB","Na","NF","Ma","Unk"))
meta1 <- meta[order(meta$Social_Group, na.last = FALSE), ] 
# meta1 <- meta1[!grepl("Blank|KBMu21611|KBNa1921", rownames(meta1)), ]
meta1 <- meta1[!grepl("Blank", rownames(meta1)), ]
bact_phyla_summary_plot <- t(t(bact_phyla_summary)[match(rownames(meta1), 
                                                 rownames(t(bact_phyla_summary))), ])

phylum_order <- c("Actinobacteria", "Bacteroidetes", "Chloroflexi", 
                  "Cyanobacteria", "Firmicutes", "Proteobacteria", "Spirochaetae",
                  "Tenericutes", "Verrucomicrobia", "Bacterial phyla < 1%", "Unclassified phyla")

png()
layout(matrix(c(1,2), 1, 2), widths = c(10, 3))
par(xpd = NA)
color_p <- c("goldenrod2",
             "navy",
             "seashell3",
             "indianred3",
             "pink3",
             "antiquewhite4",
             "green3",
             "dodgerblue",
             "firebrick3",
             "orange4",
             "darkviolet")
barplot(bact_phyla_summary_plot, 
        axes = FALSE, 
        space = 0.3, 
        axisnames = FALSE, 
        col=color_p)
axis(2, at = c(0, 20, 40, 60, 80, 100), 
     pos = -0.8, 
     cex.axis = 0.8, 
     ylab = "Relative phylum abundance (%)") 
mtext("Relative phylum abundance (%)", 
      side = 2, 
      line = 1.5, 
      cex = 1) 
legend(x = 155, 
       y = 100, 
       legend = rev(rownames(bact_phyla_summary_plot)), 
       fill = rev(color_p), 
       cex = 0.7, 
       bty = "n")
dev.off()


# bacterial genera 

genus <- read.delim("data/map_R_L6.txt", 
                    sep="\t", h=T)
rownames(genus) <- genus[,1]
genus <- genus[ , -c(1:7)]
genus <- as.data.frame(t(genus))
# remove blanks and poor qual samples KBMu21611 and KBNa1921
genus <- genus[ ,!grepl("Blank", colnames(genus))]
genus <- genus[ ,!grepl("KBMu21611|KBNa1921", colnames(genus))]

# genera present in all samples are 
genus_prev <- as.data.frame(rowSums(genus != 0))
colnames(genus_prev) <- "prevalence_n"
genus_prev100 <- subset(genus_prev, prevalence_n == 117 & !(rownames(genus_prev)
                                                         %in% 
                                                         grep("Other|uncultured|uncultured_bacterium$", rownames(genus_prev),
                                                         value=TRUE)))
genus100 <- subset(genus, rownames(genus) %in% rownames(genus_prev100))
genus100 <- genus100[grepl("D_5__",rownames(genus100)),]
rownames(genus100) <- gsub("D_1__|D_0__|D_2__|D_3__|D_4__|D_5__","", rownames(genus100))
# Calculate genus prevalence and % relative abundance (mean, SD, minimum, and maximum)
genus_summary_tab <- cbind(as.data.frame(rowSums(genus100 != 0)),
                   round(as.data.frame(rowSums(genus100 != 0))/ncol(genus100)*100, 2),
                   format(as.data.frame(rowMeans(genus100)*100), digits=3, scientific=TRUE),
                   format(as.data.frame(apply(genus100, 1, sd, na.rm = TRUE)*100), digits = 3, scientific = TRUE),
                   format(as.data.frame(apply(genus100, 1, min, na.rm = TRUE)*100), digits = 3, scientific = TRUE),
                   format(as.data.frame(apply(genus100, 1, max, na.rm = TRUE)*100), digits = 3, scientific = TRUE))
colnames(genus_summary_tab) <- c("Prevalence (n)", "Prevalence (%)", "Abundance (mean)", "Abundance (sd)", "Abundance (min)", "Abundance (max)")
rownames(genus_summary_tab)

#### LEfSe 

lefse <- read.table("data/Galaxy22-[B)_LDA_Effect_Size_(LEfSe)_on_data_15].txt",
                    sep="\t", h=T, row.names = 1)
colnames(lefse) <- c("log.class.avg","class","effectsize","p.value")
lefse <- lefse[complete.cases(lefse),]
rownames(genus) <- gsub("D_1__|D_2__|D_0__|D_3__|D_4__|D_5__","",rownames(genus))
rownames(genus) <- gsub("\\.","_",rownames(genus))
str(lefse)
lefse$p.value <- as.numeric(as.character(lefse$p.value))
lefse <- lefse[lefse$p.value < 0.01,] # keep only p < 0.01
lefse <- lefse[!grepl("uncultured|Other",rownames(lefse)),] # remove uncultured bacteria
setdiff(rownames(lefse), rownames(genus))

# change names for non matching genera
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_FamilyXIII_FamilyXIIIUCG_001"] <-
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Family_XIII_Family_XIII_UCG_001"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_LachnospiraceaeNK4A136group"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Lachnospiraceae_NK4A136_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_PrevotellaceaeUCG_003"] <- 
  "Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotellaceae_UCG_003"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_ClostridialesvadinBB60group"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiales_vadinBB60_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae__Eubacterium_oxidoreducensgroup"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae__Eubacterium__oxidoreducens_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_ClostridialesvadinBB60group_ClostridialesbacteriumJN18_A56_K"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiales_vadinBB60_group_Clostridiales_bacterium_JN18_A56_K"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae1_Sarcina"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_1_Sarcina"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_RuminococcaceaeNK4A214group"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcaceae_NK4A214_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_RuminococcaceaeUCG_008"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcaceae_UCG_008"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_RuminococcaceaeUCG_010"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcaceae_UCG_010"
rownames(lefse)[rownames(lefse)=="Bacteria_Cyanobacteria_Chloroplast_Cercisgigantea_Cercisgigantea"] <- 
  "Bacteria_Cyanobacteria_Chloroplast_Cercis_gigantea_Cercis_gigantea"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminiclostridium5"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminiclostridium_5"
rownames(lefse)[rownames(lefse)=="Bacteria_Verrucomicrobia_Opitutae_OpitutaevadinHA64"] <- 
  "Bacteria_Verrucomicrobia_Opitutae_Opitutae_vadinHA64"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcus1"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcus_1" 
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae__Eubacterium_halliigroup"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae__Eubacterium__hallii_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_BacteroidalesS24_7group"] <- 
  "Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidales_S24_7_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Cyanobacteria_Chloroplast_Cercisgigantea"] <- 
  "Bacteria_Cyanobacteria_Chloroplast_Cercis_gigantea"
rownames(lefse)[rownames(lefse)=="Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_BacteroidalesBS11gutgroup"] <- 
  "Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidales_BS11_gut_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_LachnospiraceaeND3007group"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Lachnospiraceae_ND3007_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella7"] <- 
  "Bacteria_Bacteroidetes_Bacteroidia_Bacteroidales_Prevotellaceae_Prevotella_7"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae1"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Clostridiaceae_1" 
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_RuminococcaceaeUCG_003"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Ruminococcaceae_Ruminococcaceae_UCG_003"
rownames(lefse)[rownames(lefse)=="Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_LachnospiraceaeND3007group"] <- 
  "Bacteria_Firmicutes_Clostridia_Clostridiales_Lachnospiraceae_Lachnospiraceae_ND3007_group"
rownames(lefse)[rownames(lefse)=="Bacteria_Cyanobacteria_Chloroplast_Cercisgigantea_Cercisgigantea_Cercisgigantea"] <- 
  "Bacteria_Cyanobacteria_Chloroplast_Cercis_gigantea_Cercis_gigantea_Cercis_gigantea" 

# pick out taxa based on lefse list
lefse_mat_full <- genus[rownames(genus) %in% rownames(lefse), ]
# pick out genus level
lefse_mat <- lefse_mat_full[grepl("vadinHA64|BS11|S24|Cercis_gigantea_Cercis_gigantea_Cercis_gigantea|Clostridiaceae_1|vadinBB60|Olsenella|Senegalimassilia|Alloprevotella|Prevotella_7|UCG_003|Bacillus|Kurthia|Solibacillus|Lactobacillus|Sarcina|JN18|UCG_001|ND3007|NK4A136|Oribacterium|Eubacterium|Faecalibacterium|Ruminiclostridium|NK4A214|UCG_003|UCG_008|UCG_010|Ruminococcus_1|Subdoliaranulum|Phascolarctobacterium|Sutterella|Desulfovibrionaceae_Desulfovibrio|Sphaerochaeta|Anaeroplasmataceae_Anaeroplasma|Catenisphaera|Solobacterium|Subdoligranulum",
     rownames(lefse_mat_full)),]
rownames(lefse_mat) <- c("Verrucomicrobia_Opitutae vadinHA64","Bacteroidetes_Bacteroidales BS11 gut group","Bacteroidetes_Bacteroidales S24 7 group",
                         "Firmicutes_Clostridiaceae 1", "Firmicutes_Clostridiales vadinBB60 group", "Actinobacteria_Olsenella", "Actinobacteria_Senegalimassilia",
                         "Bacteroidetes_Alloprevotella", "Bacteroidetes_Prevotella 7", "Bacteroidetes_Prevotellaceae UCG-003", "Cyanobacteria_Chloroplast Cercis gigantea",
                         "Firmicutes_Bacillus", "Firmicutes_Kurthia", "Firmicutes_Solibacillus", "Firmicutes_Lactobacillus", "Firmicutes_Sarcina", "Firmicutes_Clostridiales bacterium JN18 A56 K",
                         "Firmicutes_Clostridiales Family XIII UCG-001", "Firmicutes_Lachnospiraceae ND3007 group", "Firmicutes_Lachnospiraceae NK4A136 group",
                         "Firmicutes_Oribacterium", "Firmicutes_Eubacterium hallii group", "Firmicutes_Eubacterium oxidoreducens group", "Firmicutes_Faecalibacterium",
                         "Firmicutes_Ruminiclostridium 5", "Firmicutes_Ruminococcaceae NK4A214 group", "Firmicutes_Ruminococcaceae UCG-003",
                         "Firmicutes_Ruminococcaceae UCG-008", "Firmicutes_Ruminococcaceae UCG-010", "Firmicutes_Ruminococcus 1", "Firmicutes_Subdoligranulum",
                         "Firmicutes_Catenisphaera", "Firmicutes_Solobacterium", "Firmicutes_Phascolarctobacterium",
                         "Proteobacteria_Sutterella","Proteobacteria_Desulfovibrio","Spirochaetae_Sphaerochaeta", "Tenericutes_Anaeroplasma")

# aggregate by KB or NB 
pop <- data.frame(pop=meta$Population, row.names=row.names(meta))
lefse_mat <- as.data.frame(t(lefse_mat))
lefse_mat_2 <- merge(pop, lefse_mat, by="row.names")
rownames(lefse_mat_2) <- lefse_mat_2[,1]
lefse_mat_2 <- lefse_mat_2[,-1]
lefse_agg <- aggregate(lefse_mat_2[2:ncol(lefse_mat_2)], by=list(lefse_mat_2$pop), FUN=mean, data=lefse_mat_2)
# melt dfs 
melt_lefse_agg <- melt(lefse_agg)
colnames(melt_lefse_agg) <- c("group", "taxa", "average")
melt_lefse_agg <- melt_lefse_agg[order(as.character(melt_lefse_agg$taxa)),]
# add a column with phylum (for coloring circles later)
library(splitstackshape)
melt_lefse_agg <- data.frame(cSplit(melt_lefse_agg, "taxa", sep="_",drop=FALSE))
melt_lefse_agg <- subset(melt_lefse_agg, select = -c(taxa_2))
names(melt_lefse_agg)[4] <- c("phylum")
melt_lefse_agg$phylum <- with(melt_lefse_agg, factor(phylum, levels = rev(levels(phylum))))
melt_lefse_agg$taxa <- as.factor(melt_lefse_agg$taxa)
melt_lefse_agg$taxa <- with(melt_lefse_agg, factor(taxa, levels = rev(levels(taxa))))
melt_lefse_agg$average <- melt_lefse_agg$average*100
lefse_genus <- ggplot(melt_lefse_agg, aes(group,taxa)) +
  geom_point(aes(fill=factor(phylum), size=average), shape=21) +
  theme_bw() + ylab("") + xlab("") + guides(fill = FALSE) + 
  scale_size_continuous(range = c(0.1,9), breaks = c(0.01,0.1,1,5), labels = c('\u2264 0.01%',"0.1%","1%",'\u2265 5%'), name="Mean\nRelative\nAbundance") +
  theme(legend.position="right", legend.text.align=0, legend.title.align=0.5, 
legend.background = element_rect(fill="white", size=.5, linetype="solid", color="black"), 
axis.text.x=element_text(size=15), axis.text.y = element_text(size=12)) + 
  scale_fill_manual(values=c("firebrick3","dodgerblue","forestgreen","seashell3","pink3","indianred3","navy","goldenrod2"))
lefse_genus
