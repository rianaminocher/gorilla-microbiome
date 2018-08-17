# diversity analyses

setwd("/Volumes/riana_minocher/gorilla-microbiome-gh/") # set wd to repo folder

rm(list=ls())

load("data/vsearch_phyloseq_css.RData") ### can load workspace to start or run code below (may have issues with version metagenomeSeq)


      ### load data and cleanup ###

otu_bact <- import_biom(file.path("data/from_uc_w_tax_bact.biom"))
tree <- read_tree_greengenes("data/rep_set.tre")
taxa_names(tree) <- gsub("><-><", ";", taxa_names(tree)) # change the taxa names to match otu tab
load("data/meta.robj")

otu_bact <- merge_phyloseq(otu_bact, meta, tree) # make phylo object

# normalise using CSS
library(metagenomeSeq)
otu_bact_CSS <- otu_bact # copy a phyloseq object
otu_table_CSS <- as.data.frame(otu_bact_CSS@otu_table) # extract the otu table
obj <- newMRexperiment(as.matrix(otu_table_CSS[1:128])) # make the metagenome seq obj
p <- cumNormStat(obj) # compute the percentile to scale by
obj_sf <- cumNorm(obj, p = p) # calculates column's quantiles and calculates sum
otu_norm_CSS <- MRcounts(obj_sf, norm=T) # normalise counts
otu_bact_CSS@otu_table <- otu_table(otu_norm_CSS, taxa_are_rows = T) # replace otu table in phyloseq object

# load relatedness matrix 
relate <- read.table("data/ML-Relate_KBNP_OUT_Matrix.txt", 
                     h = TRUE, 
                     fill = TRUE, 
                     sep = "\t")
# load parentage tab
parent <- read.csv("data/mother-offspring_pairs.csv")

# edit relatedness matrix
rownames(relate) <- relate[ , 1]
relate <- relate[ , -1]
relate <- as.matrix(relate)
# change names of Chim15+Chim16
rownames(relate)[1] <- "Chim_15_Chim_16"
colnames(relate)[1] <- "Chim_15_Chim_16"
setdiff(colnames(relate), rownames(relate))
relate <- setNames(melt(relate), c("Ind1", "Ind2", "Relate"))
length(which(is.na(relate)))
# for relatedness stuff, we can only use KB
meta_relate <- meta[!grepl("NB|Blank", meta$Sample_ID), ]
length(unique(meta_relate$Composite_Indidvidual_ID))

# now all the names should fit
# trim relate list based on names in composite ID
relate <- relate[relate$Ind1 %in% meta_relate$Composite_Indidvidual_ID, ]
relate <- relate[relate$Ind2 %in% meta_relate$Composite_Indidvidual_ID, ] 
# 6400 observations for 80 individuals, remove duplicated "pairs" (with "NA")
relate <- relate[!is.na(relate$Relate), ] # now we have 3200 obs, each pair appears only once
table(relate$Ind1) + table(relate$Ind2)

# edit parentage 

parent <- parent[ , -1]
parent$Offspring.ID <- as.character(parent$Offspring.ID)
parent$Candidate.mother.ID <- as.character(parent$Candidate.mother.ID)
# check overlap of composite ids with our list
diff <- setdiff(c(parent$Offspring.ID, parent$Candidate.mother.ID), meta_relate$Composite_Indidvidual_ID)
# remove rows with nonmatching names
parent <- parent[!parent$Offspring.ID %in% diff, ]
parent <- parent[!parent$Candidate.mother.ID %in% diff, ]
setdiff(c(parent$Offspring.ID, parent$Candidate.mother.ID), meta_relate$Composite_Indidvidual_ID)
# now there should only be individuals in my list
parent$mother_offspring_pair <- "Y"

# check that blanks and samples do not cluster
# remove the poor quality samples KBMu21611 and KBNa1921

rownames(meta)[(meta$Read_Count_quality < 50000)] 
poor_qual_samples <- c("KBMu21611", "KBNa1921") # remove low read count samples and blanks
otu_bact_CSS_blanks <- subset_samples(otu_bact_CSS, 
                                      !sample_names(otu_bact_CSS) %in% poor_qual_samples) 

dist_methods <- unlist(distanceMethodList)
print(dist_methods) # we select only the ones we want 
plist <- vector("list", length(dist_methods[c(8, 10)])) # bray, unifrac, wunifrac, jaccard
dist_methods <- dist_methods[c(8, 10)]
names(plist) <- dist_methods
for(i in dist_methods){
  iDist <- phyloseq::distance(otu_bact_CSS_blanks, method = i) # build dm
  iPCoA <- ordinate(otu_bact_CSS_blanks, "PCoA", distance = iDist) # ordinate
  p <- NULL # don't carry prev plot, plot a temp "p" for each
  p <- plot_ordination(otu_bact_CSS_blanks, iPCoA) 
  p <- p + ggtitle(paste("PCoA using distance method ", i, sep=""))
  plist[[i]] <- p # save to file 
}
df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

df[grep("bray", df$distance), ]$distance <- "Bray-Curtis"
df[grep("jaccard", df$distance), ]$distance <- "Jaccard"

p_pcoa_all <- ggplot(df, aes(Axis.1, Axis.2, color=Population, shape=Population))
p_pcoa_all <- p_pcoa_all + geom_point(size=2.5) + scale_shape_manual(values=c(15,16,17))
p_pcoa_all <- p_pcoa_all + facet_wrap(~distance, scales="free") 
p_pcoa_all <- p_pcoa_all + theme_bw() + theme(strip.background=element_rect(fill=NA, colour="black"),
                                                      strip.text = element_text(size=12, face="bold"),
                                                      legend.text = element_text(size=10, face="bold"),
                                                      legend.title = element_text(size=12, face="bold"),
                                                      axis.title = element_text(size=10),
                                                      legend.key.size = unit(1, "line"))
p_pcoa_all <- p_pcoa_all + ggtitle("") + xlab("PC.1") + ylab("PC.2")
p_pcoa_all <- p_pcoa_all + scale_color_manual(name = "Population", values = c("coral","goldenrod2","midnightblue"))
p_pcoa_all

# first PC separates samples and blanks, as expected
# two samples fall between samples and blanks
# ... Chim individuals (check sampling date, taxa, read depth)

# trim to samples
samples_to_keep <- rownames(meta)[!(meta$Read_Count_quality < 50000)] 
# remove low read count samples and blanks
otu_bact_CSS <- prune_samples(samples_to_keep, otu_bact_CSS) 

# we have CSS normalised data, which is not normalising for distance-based metrics (uw unifrac, jaccard)
# we may need to use a rarefied table for that, proceed with bray

#### PCOA 
dist_methods <- unlist(distanceMethodList)
print(dist_methods) # we select only the ones we want
plist <- vector("list", length(dist_methods[c(1, 2, 8, 10)])) # bray, unifrac, wunifrac, jaccard
dist_methods <- dist_methods[c(1, 2, 8, 10)]
names(plist) <- dist_methods
for( i in dist_methods ){
  iDist <- distance(otu_bact_CSS, method=i) # build dm
  iPCoA <- ordinate(otu_bact_CSS, "PCoA", distance=iDist) # ordinate
  p <- NULL # don't carry prev plot, plot a temp "p" for each
  p <- plot_ordination(otu_bact_CSS, iPCoA) 
  p <- p + ggtitle(paste("PCoA using distance method ", i, sep=""))
  plist[[i]] <- p # save to file 
}
df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"

## 1. POPULATION
df[grep("bray", df$distance) ,]$distance <- "Bray-Curtis"
df[grep("jaccard", df$distance) ,]$distance <- "Jaccard"
df[grep("wunifrac", df$distance) ,]$distance <- "Weighted Unifrac"
df[grep("unifrac", df$distance) ,]$distance <- "Unifrac"

colours <- c("goldenrod2", "midnightblue")
p_pcoa_CSS_pop <- ggplot(df, aes(Axis.1, Axis.2, color=Population, shape=Population))
p_pcoa_CSS_pop <- p_pcoa_CSS_pop + geom_point(size=2, alpha=0.8)
p_pcoa_CSS_pop <- p_pcoa_CSS_pop + facet_wrap(~distance, scales="free")
p_pcoa_CSS_pop <- p_pcoa_CSS_pop + theme_bw() + theme(strip.background=element_rect(fill=NA, colour="black"),
                                                      strip.text = element_text(size=12, face="bold"),
                                                        legend.text = element_text(size=10, face="bold"),
                                                        legend.title = element_text(size=12, face="bold"),
                                                        axis.title = element_text(size=10),
                                                        legend.key.size = unit(1, "line"))
p_pcoa_CSS_pop <- p_pcoa_CSS_pop + ggtitle("") + xlab("PC.1") + ylab("PC.2")
p_pcoa_CSS_pop <- p_pcoa_CSS_pop + scale_color_manual(name = "Population", values = colours)
p_pcoa_CSS_pop 


      ### Dissimilarity tests ###

## 1. INDIVIDUAL

dist <- (distance(otu_bact_CSS, method = "bray"))
dist_melt <- melt(as.matrix(dist))
colnames(dist_melt) <- c("SampleID_1", 
                         "SampleID_2", 
                         "value")
library(cluster)
df <- as(sample_data(otu_bact_CSS), "data.frame")
bid <- data.frame(sname = df$Individual_ID)
bid_gower <- daisy(bid, metric = "gower") # makes the dissimilarity matrix
sid <- data.frame(id = df$Individual_ID, 
                  row.names = rownames(df)) # df of sampleIDs and indIDs
dyads <- data.frame(t(combn(row.names(sid), 2)), 
                    as.factor(as.numeric(bid_gower))) # generates all combos of samples and adds dissimilarity
colnames(dyads) <- c("SampleID_1", 
                     "SampleID_2", 
                     "Same_Ind")
dyads <- merge(dyads, 
               dist_melt, 
               by = c("SampleID_1", "SampleID_2")) # only unique dyads
wilcox.test(dyads$value ~ dyads$Same_Ind)
dyads_melt <- melt(dyads)
p_diss_ind <- ggplot(data = dyads_melt, aes(x = variable,y = value))
p_diss_ind <- p_diss_ind + 
  geom_violin(aes(fill = Same_Ind), 
              width = 0.5, 
              position = position_dodge(0.6)) + 
  stat_summary(fun.y = "median", 
               geom = "point", 
               aes(fill = Same_Ind), 
               position = position_dodge(0.6))
p_diss_ind  <- p_diss_ind + 
  ylab("Bray-Curtis dissimilarity") + 
  xlab("") + 
  theme_bw()
p_diss_ind <- p_diss_ind + 
  theme(axis.text.x = element_blank(), 
        strip.background = element_rect(size = 0), 
        strip.text = element_text(size = 10,face = "bold"), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.key.size = unit(1, "cm"))
p_diss_ind <- p_diss_ind + 
  scale_fill_manual(values = c("pink3", "goldenrod2"),
                    labels = c("Same Individual", "Different Individual"))
p_diss_ind

## 2. SOCIAL GROUP

# we need only the unique individuals here
# select one sample of each individual

meta_unique <- meta[!grepl("Blank", meta$Sample_ID_M), ] # remove blanks
meta_unique <- arrange(meta_unique, date, Composite_Indidvidual_ID)
meta_unique <- meta_unique[!duplicated(meta_unique$Individual_ID), ] # remove repeated samples
meta_unique <- meta_unique[!grepl("KBNa1921|KBMu21611", meta_unique$Sample_ID), ] # remove bad samples
unique_samples <- meta_unique$Sample_ID
rownames(meta_unique) <- meta_unique$Sample_ID
# make a unique phyloseq object, n=89
otu_bact_CSS_unique <- subset_samples(otu_bact_CSS, 
                                      sample_names(otu_bact_CSS) %in% unique_samples)

dist_unique <- (distance(otu_bact_CSS_unique, method = "bray"))
dist_unique_melt <- melt(as.matrix(dist_unique)) # n*n 89*89 = 7921
colnames(dist_unique_melt) <- c("SampleID_1", "SampleID_2", "value")
df <- as(sample_data(otu_bact_CSS_unique), "data.frame")
bsoc <- data.frame(sname = meta_unique$Social_Group_G)
bsoc_gower <- daisy(bsoc, metric = "gower") # makes the dissimilarity matrix
ssoc <- data.frame(soc = meta_unique$Social_Group_G, 
                   row.names = rownames(meta_unique)) # df of sampleIDs and indIDs
dyads_2 <- data.frame(t(combn(row.names(ssoc), 2)), 
                      as.factor(as.numeric(bsoc_gower))) # generates all combos of samples and adds dissimilarity
# correct n*(n-1)/2 (89*88)/2 = 3916
colnames(dyads_2) <- c("SampleID_1", 
                       "SampleID_2", 
                       "Same_Group")
dyads_full <- merge(dyads_2, 
                    dist_unique_melt, 
                    by = c("SampleID_1","SampleID_2"))
wilcox.test(dyads_full$value ~ dyads_full$Same_Group)

dyads_full_melt <- melt(dyads_full)

p_diss_group <- ggplot(data=dyads_full_melt, aes(x=variable,y=value))
p_diss_group <- p_diss_group + geom_violin(aes(fill=Same_Group), width=0.5, position=position_dodge(0.6)) + stat_summary(fun.y="median", geom="point", aes(fill=Same_Group), position=position_dodge(0.6))
p_diss_group  <- p_diss_group + ylab("Bray-Curtis dissimilarity") + xlab("") + theme_bw()
p_diss_group <- p_diss_group + theme(axis.text.x=element_blank(), strip.background=element_rect(size=0), strip.text=element_text(size=10,face="bold"), 
                                     legend.title=element_blank(), legend.text = element_text(size=10, face="bold"), axis.title.y = element_text(size=12, face="bold"),
                                     legend.key.size = unit(1, "cm"))
p_diss_group <- p_diss_group + scale_fill_manual(values=c("pink3","goldenrod2"), labels=c("Same Group","Different Group"))
p_diss_group

## 3. MOTHER-INFANT PAIRS

setdiff(parent$Offspring.ID, meta_full$Composite_Indidvidual_ID) # names match
# merge parent with meta df names to get Sample IDs
merge_parent <- merge(parent, meta_unique[,c("Composite_Indidvidual_ID","Sample_ID")], by.x="Offspring.ID", by.y="Composite_Indidvidual_ID")
colnames(merge_parent) <- c("Offspring.ID","Candidate.mother.ID","mother_offspring_pair","Offspring_Sample_ID")
merge_parent <- merge(merge_parent, meta_unique[,c("Composite_Indidvidual_ID","Sample_ID")], by.x="Candidate.mother.ID", by.y="Composite_Indidvidual_ID")
merge_parent <- merge_parent[,c(5,4,3)]
colnames(merge_parent) <- c("Mother_Sample_ID","Offspring_Sample_ID","mother_offspring_pair")
# merge with df based on age and sex
merge_parent <- merge(merge_parent, meta_unique[,c("Sample_ID", "Age_recode", "Sex")], 
                      by.x = c("Mother_Sample_ID"), by.y = c("Sample_ID"))
colnames(merge_parent) <- c(colnames(merge_parent)[1:3], "Mother_Age", "Mother_Sex")
merge_parent <- merge(merge_parent, meta_unique[,c("Sample_ID", "Age_recode", "Sex")], 
                      by.x = c("Offspring_Sample_ID"), by.y = c("Sample_ID"))
colnames(merge_parent) <- c(colnames(merge_parent)[1:5], "Offspring_Age", "Offspring_Sex")

# we use the unique individuals dissimilarity mat
# we need just the unrelated adult female, unrelated juvenile pairs
# need to control for group?

# merge merge_parent to dyads_full
dyads_full_off <- NULL
dyads_full <- merge(dyads_2, dist_unique_melt, by=c("SampleID_1","SampleID_2"))
# need to add age and sex data to dyads full
dyads_full_off <- merge(dyads_full, meta_unique[,c("Sample_ID", "Age_recode", "Sex")],
                        by.x= c("SampleID_1"), by.y = c("Sample_ID"))
colnames(dyads_full_off) <- c(colnames(dyads_full_off)[1:4],"SampleID_1_Age","SampleID_1_Sex")
dyads_full_off <- merge(dyads_full_off, meta_unique[,c("Sample_ID", "Age_recode", "Sex")],
                        by.x= c("SampleID_2"), by.y = c("Sample_ID"))  
colnames(dyads_full_off) <- c(colnames(dyads_full_off)[1:6],"SampleID_2_Age","SampleID_2_Sex")

# now dyads_full has two ids and age and sex per id
dyads_full_off <- merge(dyads_full_off, merge_parent, by.x=c("SampleID_1","SampleID_2"), by.y=c("Mother_Sample_ID","Offspring_Sample_ID"), all.x=T)
# assign N to non m-off pairs
dyads_full_off[is.na(dyads_full_off$mother_offspring_pair),]$mother_offspring_pair <- "N"
# trim the dataset to adult female, !adult any pairs
# is sampleid 1 or 2 the mother? sample id 1 is the mother 
identical(dyads_full_off[complete.cases(dyads_full_off),]$SampleID_1_Age, 
          dyads_full_off[complete.cases(dyads_full_off),]$Mother_Age)
identical(dyads_full_off[complete.cases(dyads_full_off),]$SampleID_1_Sex, 
          dyads_full_off[complete.cases(dyads_full_off),]$Mother_Sex)
# subset out any combos id1 that are adult & female in id1 and not adult in id2 
sub1 <- subset(dyads_full_off, dyads_full_off$SampleID_1_Age=="adult" & dyads_full_off$SampleID_1_Sex=="F" & !dyads_full_off$SampleID_2_Age=="adult")
# subset out any combos id1 that are not adult and id2 adult & female
sub2 <- subset(dyads_full_off, dyads_full_off$SampleID_2_Sex=="F" & dyads_full_off$SampleID_2_Age=="adult" & !dyads_full_off$SampleID_1_Age=="adult")
dyads_full_off_only <- rbind(sub1, sub2)
which(duplicated(dyads_full_off_only[,1:2])==T) # none are duplicated pairs

table(dyads_full_off_only$mother_offspring_pair) # 9 m-off pairs, 663 random juvenile adult female pairs

# construct bray curtis dm for this subsample of inds
juv_fem_pairs <- as.character(unique(unlist(dyads_full_off_only[,1:2])))
otu_bact_CSS_juv_fem <- subset_samples(otu_bact_CSS, sample_names(otu_bact_CSS) %in% juv_fem_pairs)
dist_bray_juv_fem <- distance(otu_bact_CSS_juv_fem, method="bray")
dist_bray_juv_fem_melt <- melt(as.matrix(dist_bray_juv_fem)) # n*n 89*89 = 7921
colnames(dist_bray_juv_fem_melt) <- c("SampleID_1", "SampleID_2", "newbray")

# merge dist values with our list
dyads_full_off_only <- merge(dyads_full_off_only, dist_bray_juv_fem_melt, by=c("SampleID_1","SampleID_2"))
setdiff(dyads_full_off_only$value, dyads_full_off_only$newbray) # dissimilarities are the same...

dyads_full_off_chim <- dyads_full_off_only[grepl("Chi", dyads_full_off_only$SampleID_1),]
dyads_full_off_chim <- dyads_full_off_chim[grepl("Chi", dyads_full_off_only$SampleID_2),]
dyads_full_off_chim <- dyads_full_off_chim[!is.na(dyads_full_off_chim$SampleID_1),]

wilcox.test(dyads_full_off_only$value ~ dyads_full_off_only$mother_offspring_pair)
wilcox.test(dyads_full_off_chim$value ~ dyads_full_off_chim$mother_offspring_pair)

dyads_full_off_melt <- melt(dyads_full_off_only[,c("SampleID_1", "SampleID_2", "mother_offspring_pair", "newbray")])
p_diss_off <- ggplot(data=dyads_full_off_melt, aes(x=variable,y=value))
p_diss_off <- p_diss_off + geom_violin(aes(fill=mother_offspring_pair), width=0.5, position=position_dodge(0.6))+ stat_summary(fun.y="median", geom="point", aes(fill=mother_offspring_pair), position=position_dodge(0.6))
p_diss_off  <- p_diss_off + ylab("Bray-Curtis dissimilarity") + xlab("") + theme_bw()
p_diss_off <- p_diss_off + theme(axis.text.x=element_blank(), strip.background=element_rect(size=0), strip.text=element_text(size=10,face="bold"), 
                                     legend.title=element_blank(), legend.text = element_text(size=10, face="bold"), axis.title.y = element_text(size=12, face="bold"),
                                     legend.key.size = unit(1, "cm"))
p_diss_off <- p_diss_off + scale_fill_manual(values=c("pink3","goldenrod2"), labels=c("Adult female \n-juvenile", "Mother-offspring"))
p_diss_off <- p_diss_off + ggtitle("a")

dyads_full_chim_melt <- melt(dyads_full_off_chim[,c("SampleID_1", "SampleID_2", "mother_offspring_pair", "newbray")])
p_diss_chim <- ggplot(data=dyads_full_chim_melt, aes(x=variable,y=value))
p_diss_chim <- p_diss_chim + geom_violin(aes(fill=mother_offspring_pair), width=0.5, position=position_dodge(0.6))+ stat_summary(fun.y="median", geom="point", aes(fill=mother_offspring_pair), position=position_dodge(0.6))
p_diss_chim  <- p_diss_chim + ylab("Bray-Curtis dissimilarity") + xlab("") + theme_bw()
p_diss_chim <- p_diss_chim + theme(axis.text.x=element_blank(), strip.background=element_rect(size=0), strip.text=element_text(size=10,face="bold"), 
                                 legend.title=element_blank(), legend.text = element_text(size=10, face="bold"), axis.title.y = element_text(size=12, face="bold"),
                                 legend.key.size = unit(1, "cm"))
p_diss_chim <- p_diss_chim + scale_fill_manual(values=c("pink3","goldenrod2"), labels=c("Adult female \n-juvenile", "Mother-offspring"))
p_diss_chim <- p_diss_chim + ggtitle("b")
p_diss_chim

# test infants in chimanuka vs others 

## ADONIS model

# use otu_bact_CSS_unique

set.seed(1193)
dist_bray_unique <- distance(otu_bact_CSS_unique, method = "bray") # calc bray curtis dist mat
df_unique <- data.frame(sample_data(otu_bact_CSS_unique)) # make df of sample_data
# adonis test
adonis(dist_bray_unique ~ Social_Group_G + Age_recode + Sex, data = df_unique, permutations = 999)
adonis_mod <- adonis(dist_bray_unique ~ Social_Group_G + Population + Age + Sex, data = df_unique, permutations = 9999)
adonis(dist_bray_unique ~ Age_recode, data = df_unique)
adonis(dist_bray_unique ~ Sex, data = df_unique)
adonis(dist_bray_unique ~  data = df_unique)
adonis(dist_bray_unique ~ v4pcr_series + Social_Group_G + lp_series, data = df_unique)

disper_test <- betadisper(dist_bray_unique, df_unique$Social_Group_G)
permutest(disper_test) 

# compare adonis with all samples (invalid with repeated samples)
dist_bray_test <- distance(otu_bact_CSS, method = "bray") # calc bray curtis dist mat
df_test <- data.frame(sample_data(otu_bact_CSS)) # make df of sample_data
adonis(dist_bray_test ~ Social_Group_G + Population + Age + Sex, data = df_test)

dbrda()
# https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/capscale
# https://stats.stackexchange.com/questions/257042/permanova-for-unbalanced-longitudinal-analyses
# http://onlinetrickpdf.blogspot.se/2015/09/multivariate-repeated-measurements-with.html
# http://onlinelibrary.wiley.com/store/10.1111/j.1466-8238.2009.00490.x/asset/j.1466-8238.2009.00490.x.pdf?v=1&t=j5i95t9f&s=4108767fd7184c1f59f17d685b7c7398a2801f80
# http://www.matthewpintar.net/uploads/2/9/8/5/29857083/pintar_and_resetarits_2017_ecology_beetles.pdf

# try a random variable as permanova sanity check


## MANTEL TESTS

# test for correlation between relatedness matrix and dissimilarity matrix

# we need unique individuals only, for KB
# meta_relate is already subset to KB
# sort by composite individual ID and then sampling date, so the df is ordered by individuals according to date
meta_relate <- arrange(meta_relate, Composite_Indidvidual_ID, date)
# choose first occurrence of unique observations
meta_relate <- meta_relate[!duplicated(meta_relate$Composite_Indidvidual_ID), ]

# update metadata in phyloseq object so it is the correct ids
otu_bact_CSS_relate <- subset_samples(otu_bact_CSS, sample_names(otu_bact_CSS) %in% meta_relate$Sample_ID)

# change names in relatedness matrix to match
meta_relate_2 <-  data.frame(meta_relate[,c("Sample_ID","Composite_Indidvidual_ID","Sample_ID_M")])
merge_relate <- merge(relate, meta_relate_2[,c("Sample_ID","Composite_Indidvidual_ID")], by.x="Ind1", by.y="Composite_Indidvidual_ID")
colnames(merge_relate) <- c("Ind1","Ind2","Relate","ID1")
merge_relate <- merge(merge_relate, meta_relate_2[,c("Sample_ID","Composite_Indidvidual_ID")], by.x="Ind2", by.y="Composite_Indidvidual_ID")
colnames(merge_relate) <- c("Ind2","Ind1","Relate","ID1","ID2")
# need both directions of the pair to build the matrix
merge_relate <- merge_relate[,c(4,5,3)]
# make a copy of this df, with the opposite pair
merge_relate_2 <- merge_relate[,c(2,1,3)]
identical(merge_relate$Relate, merge_relate_2$Relate)
# remove the 'same-individual' pairs from the second df
merge_relate_2 <- merge_relate_2[!merge_relate_2$ID1==merge_relate_2$ID2,]
# make column names the same
colnames(merge_relate_2) <- colnames(merge_relate)
# bind the dfs together
merge_relate_mat <- rbind(merge_relate, merge_relate_2)
# convert to matrix
merge_relate_mat <- acast(merge_relate_mat, ID1~ID2, value.var="Relate")

# compute dissimilarity matrix
dist_bray_relate <- distance(otu_bact_CSS_relate, method = "bray")
dist_bray <- as.matrix(dist_bray_relate)
# order to match relatedness matrix
setdiff(colnames(dist_bray), colnames(merge_relate_mat))
setdiff(colnames(merge_relate_mat), colnames(dist_bray))
merge_relate_mat <- merge_relate_mat[rownames(dist_bray), colnames(dist_bray)]
# test for correlation 
mantel(dist_bray, merge_relate_mat)