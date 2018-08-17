# code for generating plots in Minocher_ESEB_poster.pdf

# run scripts otu_analysis.R and vsearch_beta.R first


      ### Sample map Figure 1 ###

# load ggmap http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
library(ggmap) 

# load map data
mapdata <- read.csv("data/mapdata.csv", sep = ";")
mapdata$Social_Group <- as.character(mapdata$Social_Group)
colnames(mapdata) <- c("group", "lat", "lon", "pop")

mapdata <- mapdata[!mapdata$group == "Langa", ] # remove Langa
mapdata[mapdata$group == "NB", ]$group <- "NRCA"
mapdata[mapdata$group == "Part of Mankoto", ]$group <- "Ma"
mapdata[mapdata$group == "Mpungwe", ]$group <- "Mp"
mapdata[mapdata$group == "Chim", ]$group <- "Ch"

mean_KB <- sapply(mapdata[mapdata$pop == "KB", c(3, 2)], mean)

pdf("poster/plots/poster_plot_map_1.pdf", width = 5, height = 3)
map_KB <- get_map(location = c(lon = 28.70687, lat = -2.333867), 
                  maptype = "terrain", 
                  source = "google",
                  zoom = 13)
map_KB <-
  ggmap(map_KB) + 
  geom_point(data = mapdata,
             color = "indianred2", 
             shape = 17,
             size = 5) +
  geom_text(data = mapdata, 
            aes(label = paste("  ", as.character(group), sep="")),
            size = 4, 
            colour = "black")
map_KB <- map_KB +
  theme_bw() +
  xlab("") + ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
map_KB
dev.off()

pdf("poster/plots/poster_plot_map_2.pdf", width = 8, height = 7)
map_all <- get_map(location = mapdata[8, c(3,2)], 
                   maptype = "terrain", 
                   source = "stamen",
                   zoom = 8)
map_all <-
  ggmap(map_all) + 
  geom_point(data = mapdata[7:8, ],
             shape = 17,
             size = 8,
             mapping = aes(x = lon, y = lat, color = pop))
map_all <- map_all + 
  scale_color_manual(values = c("indianred3", "navy"),
                     labels = c("Kahuzi-Biega \n National Park \n (KBNP)", 
                                "Nkuba Research \n and Conservation Area \n (NRCA)"))
map_all <- map_all + 
  theme_bw() +
  xlab("") + ylab("") +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
map_all
dev.off()



      ### PCoA Figure 2 ###

dist_methods <- unlist(distanceMethodList)
plist <- vector("list", length(dist_methods[8])) # bray
dist_methods <- dist_methods[8]
names(plist) <- dist_methods
for( i in dist_methods ){
  iDist <- distance(otu_bact_CSS, method=i) # build dm
  iPCoA <- ordinate(otu_bact_CSS, "PCoA", distance=iDist) # ordinate
  p <- NULL # don't carry prev plot, plot a temp "p" for each
  p <- plot_ordination(otu_bact_CSS, iPCoA) 
  plist[[i]] <- p # save to file 
}
df <- ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
df[grep("bray", df$distance) ,]$distance <- ""
pcoa_bray <- ordinate(otu_bact_CSS, 
                      "PCoA", 
                      distance = distance(otu_bact_CSS, method = "bray"))
pcoa_bray_plot <- plot_ordination(otu_bact_CSS,
                                  pcoa_bray)

pdf("poster/plots/poster_plot_2.pdf", width = 6.8, height = 5)
df$Population <- revalue(df$Population, c("KB" = "KBNP", 
                                          "NB" = "NRCA"))
pcoa_bray_plot <- ggplot(df, 
                         aes(Axis.1, Axis.2, 
                             color = Population, 
                             shape = Population,
                             fill = Population))
pcoa_bray_plot <- pcoa_bray_plot + 
  geom_point(size = 3.5, alpha = 0.85, 
             aes(color = Population, 
                 shape = Population,
                 fill = Population))
pcoa_bray_plot <- pcoa_bray_plot + 
  facet_wrap(~distance, scales = "free")
pcoa_bray_plot <- pcoa_bray_plot + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 14),
        legend.key.size = unit(1, "line"),
        legend.position = c(1.2, 0.9),
        legend.key.width=unit(1, "lines"), 
        plot.margin = unit(c(0.5, 8, 0.3, 0.5), "lines"))
pcoa_bray_plot <- pcoa_bray_plot + ggtitle("") + xlab("PC.1: 13.2%") + ylab("PC.2: 7.1%")
pcoa_bray_plot <- pcoa_bray_plot + scale_color_manual(name = "Population", values = c("indianred", "navy")) 
pcoa_bray_plot
dev.off()



      ### Lefse bubble plot Figure 3 ###

lefse_mat <- lefse_mat_full[grepl("BS11|S24|Cercis_gigantea_Cercis_gigantea_Cercis_gigantea|Senegalimassilia|Prevotella_7|Lactobacillus|Oribacterium|Faecalibacterium|Ruminiclostridium|UCG_008|UCG_010|Alcaligenaceae_Sutterella|Desulfovibrionaceae_Desulfovibrio|Sphaerochaeta",
                             rownames(lefse_mat_full)),]
rownames(lefse_mat) <- c("Bacteroidetes_Bacteroidales BS11 gut group",
                         "Bacteroidetes_Bacteroidales S24 7 group",
                         "Actinobacteria_Senegalimassilia",
                         "Bacteroidetes_Prevotella 7", 
                         "Cyanobacteria_Chloroplast Cercis gigantea",
                         "Firmicutes_Lactobacillus",
                         "Firmicutes_Oribacterium",
                         "Firmicutes_Faecalibacterium",
                         "Firmicutes_Ruminiclostridium 5", 
                         "Firmicutes_Ruminococcaceae UCG-008", 
                         "Firmicutes_Ruminococcaceae UCG-010", 
                         "Proteobacteria_Sutterella",
                         "Proteobacteria_Desulfovibrio",
                         "Spirochaetae_Sphaerochaeta")


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
melt_lefse_agg$group <- revalue(melt_lefse_agg$group,
                                c("Blank" = "Blank",
                                  "KB" = "KBNP",
                                  "NB" = "NRCA"))

pdf("poster/plots/poster_plot_3.pdf", width = 9, height = 5)
lefse_genus <- ggplot(melt_lefse_agg, aes(group,taxa)) +
  geom_point(aes(fill=factor(phylum), size=average), shape=21) +
  theme_bw() + ylab("") + xlab("") + guides(fill = FALSE) + 
  scale_size_continuous(range = c(0.1,9), breaks = c(0.01,0.1,1,5), labels = c('\u2264 0.01%',"0.1%","1%",'\u2265 5%'), name="Mean\nRelative\nAbundance") +
  theme(legend.position="right", 
        legend.text.align=0, 
        legend.title.align=0.5, 
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill="white", size=.5, linetype="solid", color="black"), 
        axis.text.x=element_text(size=16), axis.text.y = element_text(size=15)) + 
  scale_fill_manual(values=c("forestgreen","seashell3","goldenrod2","pink3","navy","indianred3"))
lefse_genus
dev.off()



      ### Dissimilarity violin plot Figure 4 ###

pdf("poster/plots/poster_plot_4.pdf", width = 6.75, height = 5)
p_diss_group <- ggplot(data = dyads_full_melt, 
                       aes(x = Same_Group, y = value, fill = Same_Group))
p_diss_group <- p_diss_group + 
  geom_violin(width = 0.5, 
              position = position_dodge(0.3),
              trim = FALSE) + 
  stat_summary(fun.y = "median", 
               geom = "point", 
               aes(fill = Same_Group), 
               position = position_dodge(0.3))
p_diss_group <- p_diss_group + 
  scale_fill_manual(name = "Social group",
                    values = alpha(c("indianred", "goldenrod2"), 0.8), 
                    labels = c("same", "different"))
p_diss_group <- p_diss_group + 
  geom_boxplot(aes(x = Same_Group,
                   y = value),
               width = 0.03,
               position = position_dodge(0.3),
               fill = "white")
p_diss_group  <- p_diss_group + 
  ylab("pairwise Bray-Curtis dissimilarity") + 
  xlab("") + 
  theme_classic()
p_diss_group <- p_diss_group + 
  theme(axis.text.x = element_blank(), 
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14), 
        axis.title.y = element_text(size = 14),
        legend.key.size = unit(1, "line"),
        legend.position = c(1.1, 0.8),
        legend.key.width=unit(1, "lines"), 
        plot.margin = unit(c(0.5, 7, 0.3, 0.5), "lines"))
p_diss_group
dev.off()


      ### Relatedness vs Bray Figure 5 ###

melt_dist_bray <- melt(dist_bray)
melt_merge_relate_mat <- melt(merge_relate_mat)
merge_dist_bray_relate <- merge(melt_merge_relate_mat, melt_dist_bray, by = c("Var1", "Var2"))
colnames(merge_dist_bray_relate) <- c("Ind1", "Ind2", "relate", "bray")

# remove same individual pairs (relatedness = 1)
merge_dist_bray_relate <- merge_dist_bray_relate[!merge_dist_bray_relate$relate == 1, ]

pdf("poster/plots/poster_plot_5.pdf", width = 6.75, height = 5)
p_bray_relate <- ggplot(merge_dist_bray_relate,
                        aes(x = relate, y = bray)) +
  geom_point(size = 2, colour = "navy", alpha = 0.3) 
p_bray_relate <- p_bray_relate + 
  geom_smooth(method = lm, colour = "goldenrod2", fill = "blue")
p_bray_relate <- p_bray_relate + 
  theme_classic()
p_bray_relate <- p_bray_relate + 
  ylab("Bray-Curtis dissimilarity") + 
  xlab("Genetic relatedness")
p_bray_relate <- p_bray_relate + 
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))
p_bray_relate  
dev.off()
