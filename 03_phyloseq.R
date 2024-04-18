# phyloseq from dada2 output
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(vegan)
library(DESeq2)
library(reshape2)
source("rarefaction.R")
library("SRS")
library(cowplot)
library(ggtree)
library(grid)
library(RColorBrewer)
library(ggplotify)
library(aplot)
library(ape)
library(phangorn)
library(metR)
library(maps)
library(UpSetR)
library(grImport2)


##############################
# load data
##############################

pcols <- c(`Ascension Island`="orange", `Baffin Island`="dodgerblue2", `Falkland Islands`="darkorchid1")
pcols2 <- c(Brown="darkorange1", Green="green3", Red="red2", Unknown="black")

seqtabs <- list.files(pattern="seqtab.*.RDS", recursive=T, full.names=T)
taxa <- list.files(pattern="taxa.*.RDS", recursive=T, full.names=T)
ncbi <- list.files(pattern="*.taxatable", recursive=T, full.names=T)
trees <- grep("Gachon|rbcL", list.files(pattern="*.phy.treefile", recursive=T, full.names=T), value = T, invert = T)
meta <- read.table("metadata.txt", sep="\t", header=T, na.strings = "", stringsAsFactors = F)
meta$Location <- factor(meta$Location, levels=c("Baffin Island", "Ascension Island", "Falkland Islands"))
meta$MID <- paste0("MID-", meta$MID)
rownames(meta) <- meta$MID

# make phyloseq objects
ps <- lapply(c(1), function(x){
  print(x)
  # sequence table
  seqtab <- readRDS(seqtabs[[x]])
  rownames(seqtab) <- gsub(".fastq.gz", "", rownames(seqtab))
  colnames(seqtab) <- paste0("ASV", 1:ncol(seqtab))
  # SILVA taxonomy
  tax <- readRDS(taxa[[x]])
  rownames(tax) <- paste0("ASV", 1:nrow(tax))
  colnames(tax)[8] <- "Strain"
  # add NCBI taxonomy
  nc <- read.table(ncbi[[x]], sep="\t", stringsAsFactors = F, header=F, row.names=1)
  tax <- cbind(tax, NCBI=nc[match(rownames(tax), rownames(nc)),])
  # tree
  tree <- read_tree(trees[[x]])
  ps <- phyloseq(sample_data(meta), 
                 tax_table(tax), 
                 otu_table(seqtab, taxa_are_rows=F),
                 phy_tree(tree))
  # reroot tree at Ophistokonta
  try(phy_tree(ps) <- ape::root(tree, resolve.root = T, node=ape::getMRCA(tree, taxa_names(subset_taxa(ps, Phylum == "Opisthokonta")))))
  return(ps)
})

# can we identify the species of the Unknown samples
# based on the most abundant ASV?
# un <- otu_table(subset_samples(ps[[1]], Colour == "Unknown"))
# c(apply(un, 1, function(x) head(names(rev(sort(x))))))
# subset_taxa(un, taxa_names %in% c(apply(un, 1, function(x) head(names(rev(sort(x)))))))
# subset_taxa(un, taxa_names(un) %in% c(apply(un, 1, function(x) head(names(rev(sort(x)))))))
# subset_taxa(ps[[1]], taxa_names(ps[[1]]) %in% c(apply(un, 1, function(x) head(names(rev(sort(x)))))))
# tax_table(subset_taxa(ps[[1]], taxa_names(ps[[1]]) %in% c(apply(un, 1, function(x) head(names(rev(sort(x))))))))
# un.tax <- tax_table(subset_taxa(ps[[1]], taxa_names(ps[[1]]) %in% c(apply(un, 1, function(x) head(names(rev(sort(x))))))))

# filter poor samples and taxa
ps.filt <- lapply(ps, function(x){
  tt <- tax_table(x)
  # remove metazoa taxa (contamination)
  # molluscs (Heterochonchia), sponges (Calcaronea), copepods (mostly from epiphytes)
  metazoa <- rownames(tt)[which(tt[,"Order"] == "Metazoa_(Animalia)")]
  # sanity-check where the metazoa ASVs came from
  #print(plot_bar(prune_taxa(taxa_names(x) %in% metazoa, x), fill="Family") + facet_wrap(~Type, scales="free"))
  # remove
  x <- prune_taxa(! taxa_names(x) %in% metazoa, x)
  
  # remove possible host ASVs (check genera)
  host.asvs <- sapply(na.exclude(meta$Host.genus), function(host){
    rownames(tt)[unique(c(grep(host, tt[, "Genus"]), grep(host, tt[, "NCBI"])))]
  })
  host.asvs <- unique(na.exclude(unlist(host.asvs)))
  
  x <- prune_taxa(! taxa_names(x) %in% host.asvs, x)
  # remove low-coverage samples
  x <- prune_samples(sample_sums(x)>=500, x)
  # remove left-over unknown samples
  x <- prune_samples(sample_data(x)$Colour!="Unknown", x)
  # remove rare taxa
  x <- prune_taxa(taxa_sums(x)>=1, x)
  return(x)
})

# manual inspection identified the following as host sequences: "ASV36","ASV13","ASV22","ASV20","ASV24"
ps.filt[[1]] <- prune_taxa(! taxa_names(ps.filt[[1]]) %in% c("ASV36","ASV13","ASV22","ASV20","ASV24"), ps.filt[[1]])

# the metadata table contains the best host descriptions based on morphology and ASV information

###############################
# sample map and host phylogeny
###############################

loc <- data.frame(Site=c("Ascension Island", "Baffin Island", "Falkland Islands"),
                  Lat=c(-7.9467, 65.4215, -51.7963),
                  Lon=c(-14.3559, -70.9654, -59.5236),
                  Label=c("Ascension~Island~(italic(n)==4)", 
                          "Baffin~Island~(italic(n)==17)", 
                          "Falkland~Islands~(italic(n)==17)"))
World <- map_data("world")
gg.map <- ggplot() +
  geom_polygon(data = World, aes(x=long, y = lat, group = group), fill="gray30", alpha=1) +
  geom_segment(data=loc, aes(x=Lon, y=Lat, xend=Lon+20, yend=Lat, colour=Site), size=1) +
  #geom_rect(data=loc, aes(xmin=Lon, ymin=Lat-1, xmax=Lon+50, ymax=Lat+1, fill=Site), colour="white", size=1) +
  geom_point(data=loc, aes(x=Lon, y=Lat, fill=Site, shape=Site), size=8, colour="white") +
  geom_label(data=loc, aes(x=Lon+20, y=Lat, fill=Site, label=Label), 
             colour="white", size=4, parse=T, hjust=0, fontface="bold", label.padding = unit(0.25, "lines")) +
  scale_fill_manual(values=pcols) +
  scale_shape_manual(values=c(21,24,25)) +
  scale_colour_manual(values=pcols) +
  scale_x_longitude(expand=c(0,0)) +
  scale_y_latitude(expand=c(0,0)) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank())

# add barchart
bc <- melt(with(sample_data(ps.filt[[1]]), table(Location, Colour)))
bc <- subset(bc, value>0)
bc$Colour <- factor(bc$Colour, levels=rev(c("Brown", "Red", "Green", "Unknown")))
gg.bc <- ggplot(bc, aes(x=Colour, y=value, fill=Location)) +
  geom_bar(stat="identity", position=position_stack(reverse=T)) +
  geom_text(data=bc, aes(label=value),
            position=position_stack(vjust=0.5, reverse=T), colour="white") +
  scale_fill_manual(values=pcols) +
  scale_colour_manual(values=pcols) +
  coord_flip() +
  labs(x="", y="") +
  #ggtitle("Phylogenetic group") +
  theme_classic() +
  theme(legend.position = "none", 
        plot.background = element_rect(fill="white", colour="white"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title = element_blank())

gg.map <- gg.map + annotation_custom(as.grob(gg.bc), xmin=-178, xmax=-82, ymin=-58, ymax=-16)

# phylogeny
tr.meta <- data.frame(subset(sample_data(ps.filt[[1]]), !is.na(Host.genus)))
tr.meta$MID <- rownames(tr.meta)
tr.meta$Label <- tr.meta$Host.species
tr.meta$Label[is.na(tr.meta$Label)] <- paste(tr.meta$Host.genus[is.na(tr.meta$Label)], "sp.")

# open tree of life  
tol <- read.tree("../../opentree14.9_tree/labelled_supertree/labelled_supertree_ottnames.tre")
# translate species names to MIDs
tree <- tol
for(x in subset(sample_data(ps.filt[[1]]), !is.na(Host.genus))$MID){
  # try and find species first
  sp <- gsub(" ", "_", meta[x, "Host.species"])
  print(x)
  print(sp)
  # else, try genus
  if(is.na(sp) | length(grep(sp, tree$tip.label))==0) {
    sp <- paste0(meta[x, "Host.genus"], "_")
    print(sp)
  }
  print(length(unique(tree$tip.label[grep(sp, tree$tip.label)])))
  sp <- sample(unique(tree$tip.label[grep(sp, tree$tip.label)]), 1)
  print(sp)
  tree$tip.label[grep(sp, tree$tip.label)] <- x
  tr.meta$Label[tr.meta$MID==x] <- paste(tr.meta$Label[tr.meta$MID==x], regmatches(sp, regexpr("ott.*",sp)))
}
tree <- keep.tip(tree, tree$tip.label[grep("MID-", tree$tip.label)])
# root at MID-37 and MID-33 (green/red)
tree <- root(tree, node=getMRCA(tree, c("MID-37", "MID-33")), resolve.root = T)

tree <- groupOTU(tree, split(tr.meta$MID, tr.meta$Colour))
gg.tol <- ggtree(tree, branch.length = "none") %<+% tr.meta +
  geom_tiplab(aes(label=Label), color="black", fontface="italic", align=T, offset=0.5, geom="label", label.size = NA, label.padding=unit(0.01, "lines")) +
  geom_tippoint(aes(fill=Location, shape=Location, colour=Location), size=3) +
  scale_colour_manual(values=c(pcols2, pcols), limits=force, guide = 'none') +
  scale_fill_manual(values=pcols, limits=force, guide = 'none') +
  scale_shape_manual(values=c(24,21,25)) +
  xlim(c(0,17)) +
  #scale_y_continuous(expand=c(0.005,0)) +
  geom_cladelabel(node=getMRCA(tree, subset(tr.meta, Colour=="Green")$MID), label="Green algae", align=TRUE, offset=6.5, offset.text=0.2, extend=0.25,
                  geom='label', color=c(pcols2["Green"], "white"), fill=pcols2["Green"], barsize = 3, label.padding = unit(0.5, "lines")) +
  geom_cladelabel(node=getMRCA(tree, subset(tr.meta, Colour=="Red")$MID), label="Red algae", align=TRUE, offset=6.5, offset.text=0.2, extend=0.25,
                  geom='label', color=c(pcols2["Red"], "white"), fill=pcols2["Red"], barsize = 3, label.padding = unit(0.5, "lines")) +
  geom_cladelabel(node=getMRCA(tree, subset(tr.meta, Colour=="Brown")$MID), label="Brown algae", align=TRUE, offset=6.5, offset.text=0.2, extend=0.25,
                  geom='label', color=c(pcols2["Brown"], "white"), fill=pcols2["Brown"], barsize = 3, label.padding = unit(0.5, "lines")) +
  theme(legend.position = "none")
gg.tol
saveRDS(tree, "tree.TOL.RDS")
grid.arrange(gg.map, gg.tol, ncol=1)

pdf("map.pdf", width=8, height=9.5)
grid.arrange(gg.map, gg.tol, ncol=1, heights=c(0.40, 0.60))
dev.off()

# rarefaction curves
gg.rare <- lapply(ps.filt, function(x){
  gg <- ggrare(x, step = 10, se = FALSE, label="MID", color="Location") +
    facet_grid(~Colour, scales = "free") +
    scale_color_manual(values=pcols, limits=force) +
    scale_x_continuous(expand=expansion(mult=c(0.05,0.3))) +
    scale_y_continuous(expand=expansion(mult=c(0.01,0.05))) +
    ylab("Number of ASVs") + xlab("Number of reads") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
          legend.position = "top",
          panel.grid = element_blank(),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          #axis.title.x = element_blank(),
          strip.text = element_text(size=12),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12))
  return(gg)
})
pdf("rarefaction.pdf", width=12, height=6)
  gg.rare
dev.off()

cbind(Samples=sapply(ps.filt, nsamples),
      ASVs=sapply(ps.filt, ntaxa))
      
# sample sizes
sapply(ps.filt, function(x) table(sample_data(x)$Location))
sapply(ps.filt, function(x){ summary(sample_sums(x)) })
sapply(ps.filt, function(x){ summary(taxa_sums(x)) })
lapply(ps.filt, function(x){ table(sample_data(x)$Location, sample_data(x)$Colour) })



#
# taxonomy
#

# annotation rates
sapply(ps.filt, function(x) colSums(!is.na(tax_table(x)))/ntaxa(x)*100)

# label counts
lapply(ps.filt, function(x) table(tax_table(x)[, "Kingdom"]))
lapply(ps.filt, function(x) table(tax_table(x)[, "Phylum"]))
lapply(ps.filt, function(x) table(tax_table(x)[, "Class"]))
lapply(ps.filt, function(x) table(tax_table(x)[, "Order"]))
lapply(ps.filt, function(x) table(tax_table(x)[, "Family"]))

#
# TAXONOMY plots
#
taxbarplot <- function(p, tax="Phylum"){
  p <- transform_sample_counts(p, function(y) y / sum(y))
  # manually fix some strange taxonomic labels
  pt <- tax_table(p)
  pt[,"Phylum"] <- gsub("Cryptophyceae", "Cryptophyta", pt[,"Phylum"])
  tax_table(p) <- pt
  pal <- brewer.pal(n = length(unique(tax_table(p)[,tax])), name = "Set2")
  pal[2] <- "lightblue1"
  gg.tax1 <- plot_bar(tax_glom(p, tax, NArm=F), fill=tax) +
    facet_wrap(~Colour, scales="free_x", ncol=4) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult=c(0,0))) +
    #scale_fill_viridis(discrete = TRUE, option = "C") +
    #scale_fill_brewer(palette = "Set2") +
    scale_fill_manual(values=pal) +
    #scale_fill_manual(values=terrain.colors(length(unique(tax_table(p)[,tax])))) +
    ylab("Relative abundance") +
    theme_classic() +
    theme(axis.text.x = element_text(size=7, angle=90, hjust=1, vjust=0.5),
          axis.title.x = element_blank(),
          legend.position="left",
          strip.background = element_rect(fill="black"),
          strip.text = element_text(colour="white"))
  gg.tax1 <- ggplot_gtable(ggplot_build(gg.tax1))
  stripr <- which(grepl('strip-t', gg.tax1$layout$name))
  fills <- pcols2
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', gg.tax1$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg.tax1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  gg.tax1
}
grid.draw(taxbarplot(ps.filt[[1]]))

# plot with taxa on Y axis
# boxplot of abundance on x axis
#
taxboxplot <- function(p, tax="Order"){
  p <- tax_glom(p, tax, NArm=F)
  p <- transform_sample_counts(p, function(y) y / sum(y))
  
  df <- melt(otu_table(p))
  colnames(df) <- c("MID", "ASV", "Abundance")
  # merge sample data
  df <- merge(df, data.frame(sample_data(p)), by="MID")
  # merge taxoonmy data
  df <- merge(df, data.frame(as(tax_table(p), "matrix"), ASV=rownames(tax_table(p))), by="ASV")
  df <- df[!is.na(df[, tax]),]
  df <- df[grep("uncultured", df[, tax], invert = T),]
  
  gg.tax3 <- ggplot(subset(df, Abundance>0), aes_string(x=tax, y="Abundance", colour="Colour")) +
    geom_point(alpha=0.5) +
    geom_boxplot(outlier.shape = NA) +
    coord_flip() +
    scale_y_log10(labels=c("0.1%", "1%", "10%", "100%"), breaks=c(0.001, 0.01, 0.1, 1)) +
    scale_x_discrete(limits=rev) +
    facet_wrap(~Colour, nrow=1) +
    scale_colour_manual(values=pcols2, limits=force) +
    ylab("Relative abundance") + xlab("") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(linetype="dashed"),
          legend.position="none",
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          strip.background = element_rect(fill="black"),
          strip.text = element_text(colour="white"))
  gg.tax3 <- ggplot_gtable(ggplot_build(gg.tax3))
  stripr <- which(grepl('strip-t', gg.tax3$layout$name))
  fills <- pcols2
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', gg.tax3$grobs[[i]]$grobs[[1]]$childrenOrder))
    gg.tax3$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  gg.tax3
}
grid.draw(taxboxplot(ps.filt[[1]]))
grid.draw(taxboxplot(ps.filt[[1]], "Class"))

# main figure
library(RColorBrewer)
pl <- plot_grid(taxbarplot(ps.filt[[1]]), taxboxplot(ps.filt[[1]]), ncol=1, 
                labels = c('A','B'), align = 'hv', axis = 'tblr', rel_heights = c(1, 1.2))
ggsave("taxonomy.pdf", pl, width=9, height=8, units="in")




# Peronosporomycetes 
# only keep ASVs with Peronosporomycetes order label
ps.oo <- lapply(ps.filt, function(x){
  x <- subset_taxa(x, Order == "Peronosporomycetes")
  x <- prune_samples(sample_sums(x)>=1, x)
  return(x)
})
sapply(ps.oo, ntaxa)
sapply(ps.oo, nsamples)

# only keep ASVs in the Peronosporomycetes clade
ps.oo.clade <- lapply(ps.filt, function(x){
  tr <- phy_tree(x)
  # most recent common ancestor of all ASVs with oomycete annotations (example: 9 ASVs)
  oo.mrca <- ape::getMRCA(tr, taxa_names(subset_taxa(x, Order == "Peronosporomycetes")))
  # extract full clade from this ancestor
  oo.clade <- ape::extract.clade(tr, oo.mrca)
  # plot trees
  #plot(oo.clade)
  #plot(phy_tree(subset_taxa(x, Order == "Peronosporomycetes")))
  # use tip labels of this clade to make oomycetes phyloseq object
  x <- prune_taxa(taxa_names(x) %in% oo.clade$tip.label, x)
  x <- prune_samples(sample_sums(x)>=1, x)
  return(x)
})
sapply(ps.oo.clade, ntaxa)
sapply(ps.oo.clade, nsamples)
table(tax_table(ps.oo.clade[[1]])[,"NCBI"])

lapply(ps.oo.clade, plot_tree, ladderize="left", label.tips="NCBI")

#
# write oomycetes FASTA file to allow analysis with Gachon et al.
#
fasta <- list.files(pattern="seqtab.*[12].fasta$", recursive=T, full.names=T)
for(a in 1:length(fasta)){
  fas <- readLines(fasta[a])
  fa <- which(gsub(">", "", fas) %in% colnames(otu_table(ps.oo.clade[[a]])))
  fas <- fas[sort(c(fa, fa+1))]
  gachon <- readLines("Gachon_18S_aligned.fasta")
  writeLines(c(fas, gachon), paste0(fasta[a], ".Gachon.fasta"))
}
# then return here for plotting the tree
gachon.trees <- list.files(pattern="seqtab.*Gachon.*treefile$", recursive=T, full.names=T)
gachon.trees <- lapply(gachon.trees, ape::read.tree)
tr <- gachon.trees[[1]]
m <- data.frame(Tip=tr$tip.label, Group="Gachon", Label=tr$tip.label, stringsAsFactors = F)
m$Group[grep("ASV", m$Tip)] <- "ASV"
m$Label[grep("ASV", m$Tip)] <- tax_table(ps.oo.clade[[1]])[grep("ASV", m$Tip, value=T), "NCBI"]

gachon.tr <- ggtree(tr) %<+% m +
  geom_tiplab(aes(colour=Group, label=Label), fontface="italic", align=T, size=3, offset=0.01) +
  geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
             size=3, nudge_x=0.03) +
  #scale_shape_manual(values=c(`TRUE`=22, `FALSE`=NA)) +
  scale_colour_manual(values=c("dodgerblue2", "black")) +
  geom_treescale(width=0.1, x=0, y=20, fontsize=3) +
  xlim(NA, 1.3) +
  theme(legend.position="none")
        #plot.margin = margin(r=-10, l=-10, t=-10, b=-10, unit="pt"))
pdf("gachon.pdf", width=12, height=10)
  gachon.tr
dev.off()

#
# diatoms
# Family Bacillariophyceae in Order Ochrophytes
# 
#sapply(ps.filt, function(x) length(grep("Bacillariophyceae", tax_table(x)[, "Species"])) )
# only keep ASVs with Bacillariophyceae family label
ps.dia <- lapply(ps.filt, function(x){
  x <- subset_taxa(x, Family == "Bacillariophyceae")
  x <- prune_samples(sample_sums(x)>=1, x)
  return(x)
})
sapply(ps.dia, ntaxa)
sapply(ps.dia, nsamples)

# only keep ASVs in the Bacillariophyceae clade
ps.dia.clade <- lapply(ps.filt, function(x){
  tr <- phy_tree(x)
  # most recent common ancestor of all ASVs with diatom annotations (example: 9 ASVs)
  dia.mrca <- ape::getMRCA(tr, taxa_names(subset_taxa(x, Family == "Bacillariophyceae")))
  # extract full clade from this ancestor
  dia.clade <- ape::extract.clade(tr, dia.mrca)
  # plot trees
  #plot(oo.clade)
  #plot(phy_tree(subset_taxa(x, Order == "Peronosporomycetes")))
  # use tip labels of this clade to make oomycetes phyloseq object
  x <- prune_taxa(taxa_names(x) %in% dia.clade$tip.label, x)
  x <- prune_samples(sample_sums(x)>=1, x)
  return(x)
})
sapply(ps.dia.clade, ntaxa)
sapply(ps.dia.clade, nsamples)
table(tax_table(ps.dia.clade[[1]])[,"NCBI"])

lapply(ps.dia.clade[3:4], plot_tree, ladderize="left", label.tips="NCBI")

# what's with the Archaeplastida??
ps.ap <- lapply(ps.filt, function(x){
  x <- subset_taxa(x, Phylum == "Archaeplastida")
  #x <- prune_samples(sample_sums(x)>=1, x)
  return(x)
})
sapply(ps.ap, ntaxa)
tax_table(ps.ap[[1]])[,"Family"]
tax_table(ps.ap[[1]])[,"NCBI"]
tax_table(subset_taxa(ps.ap[[1]], Order=="Chlorophyta"))[,"Genus"]
plot_bar(ps.ap[[1]], fill="Genus", x="Host.genus") + facet_wrap(~Type+Colour, scales="free")
plot_bar(ps.ap[[1]], fill="Family", x="MID") + facet_wrap(~Colour, scales="free")
plot_bar(ps.ap[[1]], fill="Genus", x="MID") + facet_wrap(~Colour, scales="free")
# some endophytes, but also some epiphytic macroalgae... let's leave them in

#
# phylogenetic tree for oomycetes (or diatoms)
#
asv_phylotree <- function(p, off=0.185, nud=0.012, dot=0.225, sc=28, xl=0.4){
  p <- transform_sample_counts(p, function(y) y / sum(y))
  tr <- phy_tree(p)
  m <- data.frame(ASV=taxa_names(p), tax_table(p))
  m$Clade <- "A"
  m$Clade[grep("Eurychasma", m$NCBI)] <- "Eurychasma"
  # for each ASV, check which locations it was observed in
  asvlocs <- t(sapply(rownames(m), function(asv){
    asv <- cbind(otu_table(p)[, asv], sample_data(p))
    return(tapply(asv[,1], asv$Location, sum)>0)
  }))
  
  # do the same for host colour
  asvcols <- t(sapply(rownames(m), function(asv){
    asv <- cbind(otu_table(p)[, asv], sample_data(p))
    return(tapply(asv[,1], asv$Colour, sum)>0)
  }))
  
  library(ggnewscale)
  m <- data.frame(m, asvlocs, asvcols)
  pcols <- c(`Ascension Island`="darkorange", `Baffin Island`="dodgerblue2", `Falkland Islands`="darkorchid1",
             `Dunstaffnage`="black")
  pcols2 <- c(Brown="darkorange", Green="green3", Red="red2", Unknown="black")
  
  ggtree(tr) %<+% m +
    geom_tiplab(aes(label=NCBI, colour=Clade), fontface="italic", align=T, offset=off, hjust=1, geom="label", label.size = NA, label.padding=unit(0.01, "lines")) +
    geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
               size=3, nudge_x=nud) +
    # host colour dots
    geom_tippoint(aes(shape=Brown), x=dot, colour=pcols2[1], fill=pcols2[1], size=3) +
    geom_tippoint(aes(shape=Green), x=dot+0.01, colour=pcols2[2], fill=pcols2[2], size=3) +
    geom_tippoint(aes(shape=Red),   x=dot+0.02, colour=pcols2[3], fill=pcols2[3], size=3) +
    scale_shape_manual(values=c(`TRUE`=21, `FALSE`=NA)) +
    #new_scale("shape") +
    # location dots
    #geom_tippoint(aes(shape=Baffin.Island),    x=dot-0.1, colour=pcols[2], fill=pcols[2], size=2) +
    #geom_tippoint(aes(shape=Ascension.Island), x=dot-0.1+0.01,  colour=pcols[1], fill=pcols[1], size=2) +
    #geom_tippoint(aes(shape=Falkland.Islands), x=dot-0.1+0.02, colour=pcols[3], fill=pcols[3], size=2) +
    #scale_shape_manual(values=c(`TRUE`=22, `FALSE`=NA)) +
    scale_colour_manual(values=c("black", "dodgerblue2")) +
    geom_treescale(width=0.1, x=0, y=sc, fontsize=3) +
    xlim(NA, xl) +
    theme(legend.position="none",
          plot.margin = margin(r=-5, l=-5, t=0, b=0, unit="pt"))
}
gg.tr <- asv_phylotree(ps.oo.clade[[1]])

#
# heat map of ASVs
#
asv_heatmap <- function(p, gg.tr){
  p <- transform_sample_counts(p, function(y) y / sum(y))
  df <- melt(otu_table(p))
  colnames(df) <- c("MID", "ASV", "Abundance")
  # merge sample data
  df <- merge(df, data.frame(sample_data(p)), by="MID")
  # merge taxoonmy data
  df <- merge(df, data.frame(as(tax_table(p), "matrix"), ASV=rownames(tax_table(p))), by="ASV")
  df$ASV <- factor(df$ASV, levels=rev(get_taxa_name(gg.tr)))
  df <- df[order(df$Location),]
  df$MID <- factor(df$MID, levels=unique(df$MID))
  # can we cluster the Eurychasma ASVs by presence/absence across locations?
  #plot(hclust(ade4::dist.binary(t(otu_table(p)), method=1), method="ward.D2"))
  #cl.ord <- get_taxa_name(ggtree(hclust(ade4::dist.binary(t(otu_table(p)), method=1), method="ward.D2")))
  #cl.ord <- get_taxa_name(ggtree(hclust(vegdist(t(otu_table(p))), method="ward.D2")))
  #df$ASV <- factor(df$ASV, levels=cl.ord)
  
  gg.heat <- ggplot(subset(df, Abundance>0), aes(x=MID, y=ASV, fill=Location)) +
    geom_tile(aes(alpha=Abundance)) +
    geom_point(colour="black", size=0.5) +
    facet_wrap(~Location, scales="free_x", ncol=3) +
    scale_fill_manual(values=pcols, limits=force) +
    scale_alpha_continuous(breaks=c(0, 0.01, 0.1, 0.25, 0.5, 0.75), labels=c("0%", "1%", "10%", "25%", "50%", "75%")) +
    scale_y_discrete(labels=as.character(sapply(rev(get_taxa_name(gg.tr)), function(l) subset(df, ASV==l)[1, "NCBI"]))) +
    #scale_y_discrete(labels=as.character(sapply(rev(cl.ord), function(l) subset(df, ASV==l)[1, "NCBI"]))) +
    xlab("Sample") +
    theme_classic() +
    theme(axis.text.y = element_text(face="italic", colour = rep(c("black", "#FF9900"), times=tapply(df$ASV, df$NCBI %in% c("Eurychasma", "Eurychasma dicksonii"), function(x) length(unique(x))))),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position="right",
          strip.text = element_text(colour="white"))
  library("ggh4x")
  bg <- data.frame(Location=names(pcols[1:3]), outline_color=pcols[1:3])
  #df$Abundance <- df$Abundance/sum(df$Abundance)
  #df$Abundance <- log10(df$Abundance)
  ggplot(subset(df, is.finite(Abundance) & Abundance>0), aes(x=MID, y=ASV, fill=Colour)) +
    geom_tile(aes(alpha=Abundance)) +
    geom_point(colour="black", size=0.5) +
    #geom_polygon(aes(x = x,y = y, color = outline_color, fill = NA), data = merge(outline, square)) + 
    facet_nested(~Location+Colour, scales="free_x") +
    scale_fill_manual(values=pcols2, limits=force) +
    scale_alpha_continuous(breaks=c(0.01, 0.1, 0.25, 0.5, 0.75, 1), labels=c("1%", "10%", "25%", "50%", "75%", "100%"), limits=c(0,1),
                           range = c(0.1, 0.9)) +
    scale_y_discrete(labels=as.character(sapply(rev(get_taxa_name(gg.tr)), function(l) subset(df, ASV==l)[1, "NCBI"]))) +
    #scale_y_discrete(labels=as.character(sapply(rev(cl.ord), function(l) subset(df, ASV==l)[1, "NCBI"]))) +
    guides(alpha = guide_legend(), colour =  "none", fill =  "none") +
    xlab("Sample") +
    theme_bw() +
    theme(#axis.text.y = element_text(face="italic", size=12, colour = rep(c("black", "dodgerblue"), times=tapply(df$ASV, df$NCBI %in% c("Eurychasma", "Eurychasma dicksonii"), function(x) length(unique(x))))),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=5, angle=90, hjust=1, vjust=0.5),
      axis.title = element_blank(),
      legend.position="none",
      plot.margin = margin(r=-2, l=-2, t=2, b=2, unit="pt"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linetype="dashed"),
      strip.text = element_text(colour="white"))
}
gg.heat <- asv_heatmap(ps.oo.clade[[1]], gg.tr)

# combined plot for oomycetes
gg <- as.ggplot(gg.heat %>% insert_left(gg.tr, width=0.75))
gg <- ggplot_gtable(ggplot_build(gg))
grid.draw(gg)
fills <- c(pcols[c(2,1,3)], rep("black", times=7))
for(a in 1:10)
{
  gg$grobs[[5]]$children[[4]]$grobs[[25]]$grobs[[a]]$grobs[[1]]$children[[1]]$gp$fill <- fills[a]  
}
pdf("oomycetes.pdf", width=12, height=6)
  grid.draw(gg)
dev.off()

# combined plot for diatoms
gg.tree <- asv_phylotree(ps.dia.clade[[1]], off=0.38, nud=0.012, dot=0.325, sc=80, xl=0.69)
gg.heat <- asv_heatmap(ps.dia.clade[[1]], gg.tree)
gg <- as.ggplot(gg.heat %>% insert_left(gg.tree, width=0.75))
gg <- ggplot_gtable(ggplot_build(gg))
grid.draw(gg)
fills <- c(pcols[c(2,1,3)], rep("black", times=7))
for(a in 1:10)
{
  gg$grobs[[5]]$children[[4]]$grobs[[25]]$grobs[[a]]$grobs[[1]]$children[[1]]$gp$fill <- fills[a]  
}
pdf("diatoms.pdf", width=16, height=12)
  grid.draw(gg)
dev.off()

###################
# alpha diversity
###################

alpha_div <- function(x){
  alpha.stats <- c("Observed", "Chao1", "ACE", "Shannon", "InvSimpson")
  rich <- lapply(x, function(x){ 
    x <- prune_samples(sample_sums(x)>=1, x)
    data.frame(estimate_richness(x, measures=alpha.stats), sample_data(x), stringsAsFactors = F) 
  })
  
  # summary table of all comparisons (kruskal-wallis test)
  kw <- lapply(c("Location", "Colour"), function(f) sapply(rich, function(x) {
    sapply(alpha.stats, function(s) kruskal.test(x[,s], x[,f])$p.value)
  }))
  print(kw)
  wc <- list(
    sapply(rich, function(x) {
      sapply(alpha.stats, function(s) wilcox.test(as.formula(paste0(s, "~Location")), 
                                                  data=subset(x, Location %in% c("Baffin Island", "Falkland Islands")))$p.value)
    }),
    sapply(rich, function(x) {
      sapply(alpha.stats, function(s) { if(any(is.na(subset(x, Colour %in% c("Brown", "Red"))[,s]))) { return(NA) } else { wilcox.test(as.formula(paste0(s, "~Colour")), 
                                                  data=subset(x, Colour %in% c("Brown", "Red")))$p.value }})
    })
  )

  print(wc)  
}
alpha_div(ps.filt)
alpha_div(ps.oo.clade)
alpha_div(ps.dia.clade)

# alpha diversity plots
# three plots horizontally (all taxa, oomycetes, diatoms)
# colour by location (some significant differences)

alpha_plots <- function(x){
  x1 <- ps.filt[[x]]
  x2 <- ps.oo.clade[[x]]
  x3 <- ps.dia.clade[[x]]
  x1 <- prune_samples(sample_sums(x1)>=1, x1)
  x2 <- prune_samples(sample_sums(x2)>=1, x2)
  x3 <- prune_samples(sample_sums(x3)>=1, x3)
  rich1 <- data.frame(estimate_richness(x1, measures=c("Shannon", "Chao1", "InvSimpson")), sample_data(x1), Tax="All")[,-c(2,5)]
  rich2 <- data.frame(estimate_richness(x2, measures=c("Shannon", "Chao1", "InvSimpson")), sample_data(x2), Tax="Oomycetes")[,-c(2,5)]
  rich3 <- data.frame(estimate_richness(x3, measures=c("Shannon", "Chao1", "InvSimpson")), sample_data(x3), Tax="Diatoms")[,-c(2,5)]
  
  # table of significance annotations
  wx <- c(wilcox.test(Chao1~Location, data=subset(rich1, Location %in% c("Baffin Island", "Falkland Islands")))$p.value,
          wilcox.test(Chao1~Location, data=subset(rich2, Location %in% c("Baffin Island", "Falkland Islands")))$p.value,
          wilcox.test(Chao1~Location, data=subset(rich3, Location %in% c("Baffin Island", "Falkland Islands")))$p.value,
          wilcox.test(Shannon~Location, data=subset(rich1, Location %in% c("Baffin Island", "Falkland Islands")))$p.value,
          wilcox.test(Shannon~Location, data=subset(rich2, Location %in% c("Baffin Island", "Falkland Islands")))$p.value,
          wilcox.test(Shannon~Location, data=subset(rich3, Location %in% c("Baffin Island", "Falkland Islands")))$p.value,
          wilcox.test(InvSimpson~Location, data=subset(rich1, Location %in% c("Baffin Island", "Falkland Islands")))$p.value,
          wilcox.test(InvSimpson~Location, data=subset(rich2, Location %in% c("Baffin Island", "Falkland Islands")))$p.value,
          wilcox.test(InvSimpson~Location, data=subset(rich3, Location %in% c("Baffin Island", "Falkland Islands")))$p.value)
  wx[wx<=0.001] <- "***"
  wx[as.numeric(wx)<=0.01] <- "**"
  wx[as.numeric(wx)<=0.05] <- "*"
  wx[as.numeric(wx)>0.05] <- "n.s."
  
  wx <- data.frame(P=wx,
             #P=paste0("italic(P)==", round(wx, 3)), 
             variable=rep(c("Chao1", "Shannon", "InvSimpson"), each=3),
             Tax=rep(c("All", "Oomycetes", "Diatoms"), times=3),
             xmin=rep(c(0.625, 1.625, 2.625), times=3),
             xmax=rep(c(1, 2, 3), times=3),
             y=rep(c(max(c(rich1$Chao1, rich2$Chao1, rich3$Chao1))*1.12,
                     max(c(rich1$Shannon, rich2$Shannon, rich3$Shannon))*1.12,
                     max(c(rich1$InvSimpson, rich2$InvSimpson, rich3$InvSimpson))*1.12), each=3))
  
  df <- melt(rbind(rich1, rich2, rich3))
  df$Location <- factor(df$Location, levels=c("Baffin Island", "Falkland Islands", "Ascension Island"))
  
  # location
  gg1 <- ggplot(df, aes(x=Tax, y=value)) +
    geom_jitter(aes(colour=Location, fill=Location), size=2, colour="white", shape=21, alpha=0.8, position=position_dodge(width=0.75)) +
    geom_boxplot(aes(colour=Location, fill=Location), alpha=0.1, outlier.shape=NA) +
    #geom_text(data=wx, mapping=aes(x=Tax, y=y, label=P), colour="black", fontface="bold") +
    #geom_segment(data=wx, mapping=aes(x=xmin, y=y/1.1*1.05, xend=xmax, yend=y/1.1*1.05), colour="black", size=2) +
    geom_text(data=wx, mapping=aes(x=(xmin+xmax)/2, y=y, label=P), colour="black", fontface="bold") +
    geom_segment(data=wx, mapping=aes(x=xmin, y=y/1.1*1.05, xend=xmax, yend=y/1.1*1.05), colour="black", size=2) +
    facet_wrap(~variable, scales="free_y") +
    scale_colour_manual(values=pcols, limits = force) +
    scale_fill_manual(values=pcols, limits = force) +
    labs(y="Estimate", x="ASVs") +
    theme_bw() +
    theme(legend.position="top",
          panel.grid = element_blank(),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          #axis.title.x = element_blank(),
          strip.text = element_text(size=12),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12))
  
  # host colour
  
  # table of significance annotations
  wx <- c(wilcox.test(Chao1~Colour, data=subset(rich1, Colour %in% c("Brown", "Red")))$p.value,
          wilcox.test(Chao1~Colour, data=subset(rich2, Colour %in% c("Brown", "Red")))$p.value,
          wilcox.test(Chao1~Colour, data=subset(rich3, Colour %in% c("Brown", "Red")))$p.value,
          wilcox.test(Shannon~Colour, data=subset(rich1, Colour %in% c("Brown", "Red")))$p.value,
          wilcox.test(Shannon~Colour, data=subset(rich2, Colour %in% c("Brown", "Red")))$p.value,
          wilcox.test(Shannon~Colour, data=subset(rich3, Colour %in% c("Brown", "Red")))$p.value,
          wilcox.test(InvSimpson~Colour, data=subset(rich1, Colour %in% c("Brown", "Red")))$p.value,
          wilcox.test(InvSimpson~Colour, data=subset(rich2, Colour %in% c("Brown", "Red")))$p.value,
          wilcox.test(InvSimpson~Colour, data=subset(rich3, Colour %in% c("Brown", "Red")))$p.value)
  wx[wx<=0.001] <- "***"
  wx[as.numeric(wx)<=0.01] <- "**"
  wx[as.numeric(wx)<=0.05] <- "*"
  wx[as.numeric(wx)>0.05] <- "n.s."
  
  wx <- data.frame(P=wx,
                   #P=paste0("italic(P)==", round(wx, 3)), 
                   variable=rep(c("Chao1", "Shannon", "InvSimpson"), each=3),
                   Tax=rep(c("All", "Oomycetes", "Diatoms"), times=3),
                   xmin=rep(c(0.6250, 1.6250, 2.6250), times=3),
                   xmax=rep(c(1, 2, 3), times=3),
                   y=rep(c(max(c(rich1$Chao1, rich2$Chao1, rich3$Chao1))*1.12,
                           max(c(rich1$Shannon, rich2$Shannon, rich3$Shannon))*1.12,
                           max(c(rich1$InvSimpson, rich2$InvSimpson, rich3$InvSimpson))*1.12), each=3))
  
  df$Colour <- factor(df$Colour, levels=c("Brown", "Red", "Green", "Unknown"))
  gg2 <- ggplot(df, aes(x=Tax, y=value)) +
    geom_jitter(aes(colour=Colour, fill=Colour), size=2, colour="white", shape=21, alpha=0.8, position=position_dodge(width=0.75)) +
    geom_boxplot(aes(colour=Colour, fill=Colour), alpha=0.1, outlier.shape=NA) +
    geom_text(data=wx, mapping=aes(x=(xmin+xmax)/2, y=y, label=P), colour="black", fontface="bold") +
    geom_segment(data=wx, mapping=aes(x=xmin, y=y/1.1*1.05, xend=xmax, yend=y/1.1*1.05), colour="black", size=2) +
    facet_wrap(~variable, scales="free_y") +
    scale_colour_manual(values=pcols2, limits = force, name = "Host") +
    scale_fill_manual(values=pcols2, limits = force, name = "Host") +
    labs(y="Estimate", x="ASVs") +
    theme_bw() +
    theme(legend.position="top",
          panel.grid = element_blank(),
          axis.title = element_text(size=14),
          axis.text = element_text(size=12),
          strip.text = element_text(size=12),
          legend.text = element_text(size=10),
          legend.title = element_text(size=12))
    plot_grid(gg1, gg2, ncol=1, labels = c('',''), align = 'hv', axis = 'tblr')
}
alpha_plots(1)

pdf("alpha_diversity.pdf", width=10, height=8)
  alpha_plots(1)
dev.off()

##################################
## beta diversity
##################################

# for figure, I want:
# 3 ordination plots horizontally (all, oomycetes, diatoms)
# no dendrograms
# 3 intersection plots of DESeq2 results
# 3 taxonomy plots to summarise taxonomic distribution

# Ordination plots
ord_plot <- function(x=1, m="NMDS", d="bray"){
  #p <- transform_sample_counts(p, function(y) y / sum(y))
  #p <- prune_samples(sample_data(p)$Colour!="Unknown", p)
  # remove ASVs that only occur in 1 sample
  #p <- prune_taxa(colSums(as(otu_table(p), "matrix")>0)>1, p)

  ord <- lapply(list(ps.filt[[x]], ps.oo.clade[[x]], ps.dia.clade[[x]]), function(p) {
    p <- prune_samples(sample_sums(p)>=50, p)
    #p <- transform_sample_counts(p, function(y) y / sum(y))
    
    ord <- ordinate(p, method=m, distance=d, trymax=1000)
    if(m=="NMDS") {
      df <- as.data.frame(ord$points)
    } else {
      df <- as.data.frame(ord$vectors[,1:2])
    }
    df <- cbind(df, sample_data(p))
    colnames(df)[1:2] <- c("Axis1", "Axis2")
    return(df)
  })
  names(ord) <- c("All", "Oomycetes", "Diatoms")
  
  hull.data <- lapply(ord, function(o){
    do.call(rbind, lapply(split(o, o$Location), function(x) x[chull(x[,c("Axis1", "Axis2")]),]))
  })
  hull.data <- cbind(do.call(rbind, hull.data), ASVs=rep(names(hull.data), times=sapply(hull.data, nrow)))
  
  ord <- cbind(do.call(rbind, ord), ASVs=rep(names(ord), times=sapply(ord, nrow)))
  ord$ASVs <- factor(ord$ASVs, levels=c("All", "Oomycetes", "Diatoms"))
  
  #print(df[which.max(df$Axis1^2),])
  ggplot(ord, aes(x=Axis1, y=Axis2, fill=Location)) +
    geom_vline(xintercept=0, linetype="dashed", colour="lightgray") +
    geom_hline(yintercept=0, linetype="dashed", colour="lightgray") +
    geom_polygon(data=hull.data, aes(colour=Location), alpha=0.1) +
    geom_point(aes(colour=Location, shape=Colour), size=4, alpha=1) +
    facet_wrap(~ASVs) +
    scale_shape_manual(name="Host", values=c(Green=15, Brown=17, Red=16, Unknown=13)) +
    scale_fill_manual(values=pcols, limits=force) +
    scale_colour_manual(values=pcols, limits=force) +
    theme_bw() +
    theme(legend.position="top",
          panel.grid = element_blank(),
          #axis.title = element_text(size=14),
          #axis.text = element_text(size=12),
          #strip.text = element_text(size=12),
          #legend.text = element_text(size=10),
          #legend.title = element_text(size=12))
    )
}
ord_plot(1)

# PERMANOVA
ad <- lapply(list(ps.filt[[1]], ps.oo.clade[[1]], ps.dia.clade[[1]]), function(p){
  p <- prune_samples(sample_sums(p)>=50, p)
  p <- prune_samples(sample_data(p)$Colour %in% c("Brown", "Red", "Green"), p)
  #print(table(sample_data(p)$Location, sample_data(p)$Colour))
  return(adonis2(phyloseq::distance(p, method="bray") ~ Location + Colour, by="margin",
                 data = data.frame(sample_data(p)), permutations = 9999))
})
# extract P-values
sapply(ad, function(x) x$`Pr(>F)`[1:2])

###################################
## differential abundance analysis
###################################

diff_abundance <- function(p){
  #p <- prune_samples(sample_sums(p)>=50, p)
  #p <- prune_samples(sample_data(p)$Colour %in% c("Brown", "Red"), p)
  # add a pseudocount to the data, to ensure that DESeq2 can calculate sizeFactors
  otu_table(p) <- otu_table(p)+1
  # run the model (by location, controlling for Colour)
  dd <- phyloseq_to_deseq2(p, ~Colour+Location)
  dd <- DESeq(dd)
  #dd  <- DESeq(dd, test="LRT", reduced=~Colour)
  # extract results (three possible contrasts)
  res1 <- results(dd, contrast=c("Location", "Falkland.Islands", "Ascension.Island"))
  dd$Location <- relevel(dd$Location, ref = "Falkland Islands")
  dd <- nbinomWaldTest(dd)
  res2 <- results(dd, contrast=c("Location", "Baffin.Island", "Ascension.Island"))
  
  dd$Location <- relevel(dd$Location, ref = "Ascension Island")
  dd <- nbinomWaldTest(dd)
  res3 <- results(dd, contrast=c("Location", "Falkland.Islands", "Baffin.Island"))
  
  
  print(summary(res1))
  print(summary(res2))
  print(summary(res3))
  
  # ALTERNATIVE MODEL (Colour, controlling for Location)
  dd <- phyloseq_to_deseq2(p, ~Location+Colour)
  dd  <- DESeq(dd)
  # extract results (all possible contrasts)
  res4 <- results(dd, contrast=c("Colour", "Red", "Brown"))
  res5 <- results(dd, contrast=c("Colour", "Green", "Brown"))
  res6 <- results(dd, contrast=c("Colour", "Green", "Red"))
  #res7 <- results(dd, contrast=c("Colour", "Unknown", "Brown"))
  #res8 <- results(dd, contrast=c("Colour", "Unknown", "Red"))
  #res9 <- results(dd, contrast=c("Colour", "Unknown", "Green"))
  print(summary(res4))
  print(summary(res5))
  print(summary(res6))
  #print(summary(res7))
  #print(summary(res8))
  #print(summary(res9))
  
  r1 <- rownames(subset(res1, padj<=0.05))
  r2 <- rownames(subset(res2, padj<=0.05))
  r3 <- rownames(subset(res3, padj<=0.05))
  r4 <- rownames(subset(res4, padj<=0.05))
  r5 <- rownames(subset(res5, padj<=0.05))
  r6 <- rownames(subset(res6, padj<=0.05))
  
  # for intersection plot
  ups <- data.frame(ASV=rownames(res1), 
                    `Falkland-vs-Baffin`=0,
                    `Falkland-vs-Ascension`=0,
                    `Baffin-vs-Ascension`=0,
                    `Red-vs-Brown`=0,
                    `Green-vs-Brown`=0,
                    `Green-vs-Red`=0)
  rownames(ups) <- ups$ASV
  ups[rownames(subset(res1, padj<=0.05)), 2] <- 1
  ups[rownames(subset(res2, padj<=0.05)), 3] <- 1
  ups[rownames(subset(res3, padj<=0.05)), 4] <- 1
  ups[rownames(subset(res4, padj<=0.05)), 5] <- 1
  ups[rownames(subset(res5, padj<=0.05)), 6] <- 1
  ups[rownames(subset(res6, padj<=0.05)), 7] <- 1
  
  return(list(Location=unique(c(r1, r2, r3)), Colour=unique(c(r4, r5, r6)), ups=ups))
}
dabs <- lapply(list(ps.filt[[1]], ps.oo.clade[[1]], ps.dia.clade[[1]]), diff_abundance)
sapply(dabs, function(x) length(x$Location))
sum(dabs[[2]]$Location %in% dabs[[1]]$Location)
sum(dabs[[3]]$Location %in% dabs[[1]]$Location)
sapply(dabs, function(x) length(x$Colour))
sum(dabs[[2]]$Colour %in% dabs[[1]]$Colour)
sum(dabs[[3]]$Colour %in% dabs[[1]]$Colour)

# simple intersection plot
# (don't split by contrast)
upsplot <- function(d){
  ups <- data.frame(ASV=taxa_names(ps.filt[[1]]), 
                    `Location`=0, `Colour`=0)
  rownames(ups) <- ups$ASV
  ups$Location[ups$ASV %in% d[[1]]] <- 1
  ups$Colour[ups$ASV %in% d[[2]]] <- 1
  
  svg("upsettmp.svg", width=5, height=3)
  print(upset(ups, keep.order=T, sets=c("Colour", "Location"),
        point.size=8, line.size=2, text.scale=2, mb.ratio = c(0.5, 0.5)))
  #c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)))
  
  dev.off()
  # 2) import SVG with grImport2 and capture as grob
  grid::grid.newpage()
  grid.picture(readPicture("upsettmp.svg"), expansion=0.01) 
  gg <- grid::grid.grab(wrap=T)
  gg
}
do.call(arrangeGrob, c(lapply(dabs, upsplot), ncol=3))

upsplot2 <- function(p){
  ups <- data.frame(ASV=taxa_names(ps.filt[[1]]), 
                    `Location`=0, `Host`=0)
  rownames(ups) <- ups$ASV
  ups$Location[ups$ASV %in% taxa_names(p) & ups$ASV %in% dabs[[1]]$Location] <- 1
  ups$Host[ups$ASV %in% taxa_names(p) & ups$ASV %in% dabs[[1]]$Colour] <- 1
  ggplotify::as.ggplot(upset(ups, keep.order=T, sets=c("Host", "Location"), mb.ratio = c(0.6, 0.4),
                             point.size=4, line.size=1, text.scale=c(1, 1.2, 1.2, 1.2, 1.2, 1.2))) + xlim(c(0.25,0.99)) + ylim(c(0.2, 1))
}

# complex intersection plot (contrasts)

gg.ups1 <- ggplotify::as.ggplot(upset(dabs[[1]]$ups[,1:4],
                                text.scale=c(1, 1.3, 1.3, 1.3, 1.3, 1.3), point.size=4, mb.ratio = c(0.5, 0.5), line.size = 1,
                                queries=list(list(query = intersects, params = list("Falkland.vs.Baffin", "Falkland.vs.Ascension"), color = pcols["Falkland Islands"], active = TRUE),
                                             list(query = intersects, params = list("Falkland.vs.Baffin", "Baffin.vs.Ascension"), color = pcols["Baffin Island"], active = TRUE),
                                             list(query = intersects, params = list("Falkland.vs.Ascension", "Baffin.vs.Ascension"), color = pcols["Ascension Island"], active = TRUE),
                                             list(query = intersects, params = list("Falkland.vs.Baffin"), color = pcols["Falkland Islands"], active = TRUE),
                                             #list(query = intersects, params = list("Falkland.vs.Ascension"), color = pcols["Ascension Island"], active = TRUE),
                                             list(query = intersects, params = list("Baffin.vs.Ascension"), color = pcols["Baffin Island"], active = TRUE)))) + xlim(c(0.25,0.99)) + ylim(c(0.2, 1))
gg.ups2 <- ggplotify::as.ggplot(upset(dabs[[1]]$ups[,-c(2:4)],
                                text.scale=c(1, 1.2, 1.2, 1.2, 1.2, 1.2), point.size=4, mb.ratio = c(0.5, 0.5), line.size = 1,
                                queries=list(#list(query = intersects, params = list("Green.vs.Brown", "Green.vs.Red"), color = pcols2["Green"], active = TRUE),
                                             list(query = intersects, params = list("Green.vs.Brown"), color = pcols2["Green"], active = TRUE),
                                             list(query = intersects, params = list("Green.vs.Red"), color = pcols2["Green"], active = TRUE),
                                             list(query = intersects, params = list("Red.vs.Brown"), color = pcols2["Brown"], active = TRUE),
                                             #list(query = intersects, params = list("Red.vs.Brown", "Green.vs.Brown"), color = pcols2["Brown"], active = TRUE),
                                             list(query = intersects, params = list("Red.vs.Brown", "Green.vs.Red"), color = pcols2["Red"], active = TRUE)))) + xlim(c(0.25,0.99)) + ylim(c(0.2, 1))
o <- ord_plot(1)
pdf("betadiv.pdf", width=9, height=8)
  #grid.arrange(o, arrangeGrob(gg.ups1, gg.ups2, ncol=2), ncol=1)
  #grid.arrange(o, arrangeGrob(upsplot2(ps.filt[[4]]),
  #                            upsplot2(ps.oo.clade[[4]]),
  #                            upsplot2(ps.dia.clade[[4]]), ncol=3),
  #                arrangeGrob(gg.ups1, gg.ups2, ncol=2),
  #             ncol=1, heights=c(1,0.47,0.53))
  plot_grid(o, arrangeGrob(upsplot2(ps.filt[[1]]),
                              upsplot2(ps.oo.clade[[1]]),
                              upsplot2(ps.dia.clade[[1]]), ncol=3),
               arrangeGrob(gg.ups1, gg.ups2, ncol=2),
               ncol=1, rel_heights=c(1,0.47,0.53), labels=c("A", NA, "B"))
dev.off()

############################
## phylosymbiosis
############################
#
# distance-based
#
phylosym.mantel <- function(ps, ulv=TRUE){
  # keep samples with host genus
  ps <- prune_samples(!is.na(sample_data(ps)$Host.genus), ps)
  ps <- prune_samples(sample_sums(ps)>=50, ps)
  # remove ASVs that only occur in 1 sample
  #ps <- prune_taxa(colSums(as(otu_table(ps), "matrix")>0)>1, ps)
  #ps <- transform_sample_counts(ps, function(y) y / sum(y))
  
  # microbiome distances
  micro.dist <- phyloseq::distance(ps, method="bray")
  #micro.dist <- phyloseq::distance(ps, method="wunifrac")
  micro.dist <- melt(as.matrix(micro.dist))
  micro.dist[,1] <- as.character(micro.dist[,1])
  micro.dist[,2] <- as.character(micro.dist[,2])
  micro.dist <- micro.dist[micro.dist[,1] != micro.dist[,2],]
  
  # host distances
  # the MIDs involved in the microbiome distances need a corresponding entry for the hosts
  host.dist <- read.table("rbcL.best.mafft.trimal.phy.mldist", fill = T, row.names = 1, stringsAsFactors = F)[-1,]
  colnames(host.dist) <- rownames(host.dist)
  host.dist <- melt(as.matrix(host.dist))
  host.dist[,1] <- as.character(host.dist[,1])
  host.dist[,2] <- as.character(host.dist[,2])
  host.dist <- host.dist[host.dist[,1] != host.dist[,2],]
  
  # for each MID, find best species
  for(x in unique(c(micro.dist[,1], micro.dist[,2]))) {
    # try and find species first
    sp <- gsub(" ", "_", meta[x, "Host.species"])
    # else, try genus
    if(is.na(sp) | sum(host.dist[,1]==sp)==0) {
      ge <- meta[x, "Host.genus"]
      opt <- unique(host.dist[,1][grep(ge, host.dist[,1])])
      if(length(opt)>0) {
        sp <- sample(opt, 1)  
      } else {
        sp <- NA
      }
    }
    if(!is.na(sp)) {
      host.dist[,1][grep(sp, host.dist[,1])] <- x
      host.dist[,2][grep(sp, host.dist[,2])] <- x  
    }
  }
  host.dist <- host.dist[host.dist[,1] != host.dist[,2],]
  host.dist <- host.dist[grep("MID", host.dist$Var1),]
  host.dist <- host.dist[grep("MID", host.dist$Var2),]
  
  # synchronise both distance objects
  mids <- intersect(unique(c(micro.dist[,1], micro.dist[,2])), unique(c(host.dist[,1], host.dist[,2])))
  micro.dist <- subset(micro.dist, Var1 %in% mids & Var2 %in% mids)
  host.dist <- subset(host.dist, Var1 %in% mids & Var2 %in% mids)
  micro.dist <- micro.dist[order(micro.dist[,1], micro.dist[,2]),]
  host.dist <- host.dist[order(host.dist[,1], host.dist[,2]),]
  print(identical(host.dist[,1], micro.dist[,1]))
  print(identical(host.dist[,2], micro.dist[,2]))
  
  host.dist$Comb <- apply(host.dist, 1, function(x) paste(sort(c(x[1], x[2])), collapse="") )
  host.dist <- host.dist[!duplicated(host.dist$Comb),]
  micro.dist$Comb <- apply(micro.dist, 1, function(x) paste(sort(c(x[1], x[2])), collapse="") )
  micro.dist <- micro.dist[!duplicated(micro.dist$Comb),]
  
  print(identical(host.dist[,1], micro.dist[,1]))
  print(identical(host.dist[,2], micro.dist[,2]))
  
  # filter by colour (brown, red)
  lapply(c("*", "Brown", "Red"), function(hcl){
    # always include Ulva (green algae)
    #ulva <- rownames(meta[grep("Ulva", meta$Host.species),])
    ulva <- rownames(meta[grep("Green", meta$Colour),])
    cl <- c(rownames(meta)[grep(hcl, meta$Colour)], ulva) # include Ulva by default
    md <- subset(micro.dist, Var1 %in% cl & Var2 %in% cl)
    hd <- subset(host.dist, Var1 %in% cl & Var2 %in% cl)
    print(identical(hd[,1], md[,1]))
    print(identical(hd[,2], md[,2]))
    
    # remove Ulva?
    if(ulv==FALSE){
      hd <- subset(hd, ! Var1 %in% ulva)
      hd <- subset(hd, ! Var2 %in% ulva)
      md <- subset(md, ! Var1 %in% ulva)
      md <- subset(md, ! Var2 %in% ulva)
    }
    
    # plot dataframe
    df <- data.frame(hd[,1:2], Host=hd[,3], Micro=md[,3], Group="default", stringsAsFactors = F)
    # highlight Ulva
    df[df$Var1 %in% ulva, "Group"] <- "Ulva"
    df[df$Var2 %in% ulva, "Group"] <- "Ulva"
    df <- df[order(df$Group),]
    #print(head(df[order(df$Host, decreasing=T),]))
    
    hd <- with(hd, structure(value, Labels = unique(c(hd[,1], hd[,2])), Size = length(unique(c(hd[,1], hd[,2]))),
                       Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
    
    md <- with(md, structure(value, Labels = unique(c(md[,1], md[,2])), Size = length(unique(c(md[,1], md[,2]))),
                             Diag = FALSE, Upper = FALSE, method = "user", class = "dist"))
    
    mt <- mantel(hd, md, permutations=10000)
    mt.lbl <- ifelse(mt$signif<0.001, paste0("italic(r)==", round(mt$statistic, 3), "*';'~italic(P)<0.001"),
                     paste0("italic(r)==", round(mt$statistic, 3), "*';'~italic(P)==", round(mt$signif, 3)))
    
    #gg <- ggplot(subset(df, Group=="default"), aes(x=Host, y=Micro, colour=Group)) +
    gg <- ggplot(df, aes(x=Host, y=Micro, colour=Group, group=NA)) +
      geom_point(shape=19, size=3, alpha=0.7) +
      geom_smooth(method="gam", colour="dodgerblue2", fill="lightgray") +
      annotate("label", x=Inf, y=-Inf, hjust=1, vjust=0, label=mt.lbl, parse=T, colour="dodgerblue2") +
      scale_colour_manual(values=unname(c(ifelse(is.na(pcols2[hcl]), "black", pcols2[hcl]), "chartreuse2"))) +
      scale_y_continuous(expand=expansion(mult=c(0.2,0.05))) +
      xlab("Host phylogenetic distance") + ylab("Bray-Curtis distance") +
      theme_bw() + 
      theme(legend.position="none",
            panel.grid = element_blank())
    
    return(gg)
  })
}

mantel1 <- phylosym.mantel(ps.filt[[1]])
mantel2 <- phylosym.mantel(ps.oo.clade[[1]])
mantel3 <- phylosym.mantel(ps.dia.clade[[1]])
mantel4 <- phylosym.mantel(ps.filt[[1]], ulv=F)
mantel5 <- phylosym.mantel(ps.oo.clade[[1]], ulv=F)
mantel6 <- phylosym.mantel(ps.dia.clade[[1]], ulv=F)

p1 <- arrangeGrob(grobs=c(mantel1, mantel4), ncol=3, top=textGrob("All ASVs",gp=gpar(fontsize=14,font=2), x=0.01, hjust=0))
p2 <- arrangeGrob(grobs=c(mantel2, mantel5), ncol=3, top=textGrob("Oomycete ASVs",gp=gpar(fontsize=14,font=2), x=0.01, hjust=0))
p3 <- arrangeGrob(grobs=c(mantel3, mantel6), ncol=3, top=textGrob("Diatom ASVs",gp=gpar(fontsize=14,font=2), x=0.01, hjust=0))

pdf("FigS4.pdf", width=12, height=20)
grid.arrange(p1, p2, p3, ncol=1) #, heights=c(1,0.5,1))
dev.off()

#
# tree-based
#

# Robinson-Fould function
RobinsonFould <- function(host.tree, micro.tree, nRandom=10000){
  
  # compute observed RF distance
  obs <- RF.dist(host.tree, micro.tree, rooted=F, check.labels = T, normalize = T)
  
  # generate null trees by shuffling tips of microbiome tree ("null model")
  nullSymbiontDendro=lapply(1:nRandom, function(x) {
    tr <- micro.tree
    tr$tip.label=sample(tr$tip.label)
    return(tr)
  })
  
  # compute RF distance with host tree for each null tree    
  null <- sapply(nullSymbiontDendro, function(x) {
    max(c(RF.dist(host.tree, x, rooted=F, check.labels = T, normalize = T)))
  })
  
  # generate random trees with same tip labels as microbiome tree ("random model")
  tips <- as.phylo(micro.tree)$tip.label
  randomSymbiontDendro <- rmtree(N=nRandom, n=length(tips), tip.label = tips)
  
  # compute RF distance with host tree for each random tree
  rd <- sapply(randomSymbiontDendro, function(x) {
    max(c(RF.dist(host.tree, x, rooted=F, check.labels = T, normalize = T)))
  })
  
  # compute P-values
  # from null trees
  pval <- sum(null<=obs)/(nRandom)
  # from random trees
  pvalRd <- sum(rd<=obs)/(nRandom)
  
  # organise results in small table
  res <- matrix(NA,ncol=2,nrow=2,dimnames = list(c("stat","pval"),c("RF","RFrd")))
  #obs <- obs / (Nnode(host.tree) + Nnode(micro.tree) - 2)
  res["stat",]=rep(obs,2)
  res["pval","RF"]=pval
  res["pval","RFrd"]=pvalRd
  return(res)
}

# could use the Open Tree of Life tree
tol <- tree

# or get host tree from rbcL
tree <- read.tree("rbcL.best.mafft.trimal.phy.treefile")
#tree <- ape::root(tree, outgroup=c("Ulva_clathrata", "Ulva_flexuosa", "Blidingia_minima"), resolve.root=T)
tree <- ape::root(tree, node=ape::getMRCA(tree, c("Ulva_clathrata", "Jania_rubens")), resolve.root=T)
# flip a branch for clarity
tree <- ape::rotate(tree, c("Hypnea_flexicaulis", "Palisada_perforata"))

# microbiome dendrogram
# keep samples with host genus

phylosymbio.RF <- function(ps, tr){
  ps <- prune_samples(!is.na(sample_data(ps)$Host.genus), ps)
  ps <- prune_samples(sample_sums(ps)>=50 | sample_data(ps)$Colour=="Green", ps)
  # remove ASVs that only occur in 1 sample
  #ps <- prune_taxa(colSums(as(otu_table(ps), "matrix")>0)>1, ps)
  ps <- transform_sample_counts(ps, function(y) y / sum(y))
  
  # either only Brown, only red or everything 
  gg.trees <- lapply(list(c("Red", "Brown", "Green"), c("Brown"), c("Red")), function(cl){
    p <- prune_samples(sample_data(ps)$Colour %in% c(cl, "Green"), ps)  
    dg <- ape::as.phylo(hclust(phyloseq::distance(p, method="bray"), method="ward.D2"))
    #dg <- ape::as.phylo(hclust(phyloseq::distance(p, method="wunifrac"), method="ward.D2"))
    # reroot at a green alga (ideally all, but not always possible)
    dg <- ape::root(dg, outgroup=rownames(subset(meta, Host.species=="Ulva clathrata")), resolve.root=T)
    try(dg <- ape::root(dg, outgroup=intersect(rownames(subset(meta, Colour=="Green" & ! is.na(Host.genus))), dg$tip.label)[3], resolve.root=T))
    try(dg <- ape::root(dg, outgroup=intersect(rownames(subset(meta, Colour=="Green" & ! is.na(Host.genus))), dg$tip.label), resolve.root=T))

    # translate species names in host tree to MIDs
    for(x in dg$tip.label) {
      # try and find species first
      sp <- gsub(" ", "_", meta[x, "Host.species"])
      # else, try genus
      if(is.na(sp) | sum(tr$tip.label==sp)==0) {
        ge <- meta[x, "Host.genus"]
        opt <- unique(tr$tip.label[grep(ge, tr$tip.label)])
        if(length(opt)>0) {
          sp <- sample(opt, 1)  
        } else {
          sp <- NA
        }
      }
      if(!is.na(sp)) tr$tip.label[grep(sp, tr$tip.label)] <- x
    }
    
    # remove Green outgroups for brown and red dataset
    if(length(grep("Green", cl))==0) dg <- ape::drop.tip(dg, rownames(subset(meta, Colour=="Green")))
    tr <- ape::keep.tip(tr, intersect(tr$tip.label, dg$tip.label))
    dg <- ape::keep.tip(dg, intersect(tr$tip.label, dg$tip.label))
    
    rf <- RobinsonFould(tr, dg)
    rf.lbl <- ifelse(rf[2,1]<0.001, paste0("italic(RF)==", round(rf[1,1], 3), "*';'~italic(P)<0.001"),
                     paste0("italic(RF)==", round(rf[1,1], 3), "*';'~italic(P)==", round(rf[2,1], 3)))
    
    # tree plot
    tr.meta <- data.frame(sample_data(ps))
    tr.meta$MID <- rownames(tr.meta)
    tr.meta$Label <- tr.meta$Host.species
    tr.meta$Label[is.na(tr.meta$Label)] <- tr.meta$Host.genus[is.na(tr.meta$Label)]
    
    gg.tr <- ggtree(tr, branch.length = "none") %<+% tr.meta +
      geom_tiplab(aes(label=Label), fontface="italic", align=T, offset=5.8, hjust=1, geom="label", label.size = NA, label.padding=unit(0.01, "lines")) +
      geom_tippoint(aes(colour=Colour, fill=Colour, shape=Location), size=3) +
      geom_text2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 80), size=3, nudge_x=0.4) +
      scale_shape_manual(values=c(`Ascension Island`=21,`Baffin Island`=24, `Falkland Islands`=25)) +
      scale_fill_manual(values=pcols2, limits=force, guide = "none") +
      scale_colour_manual(values=pcols2, limits=force, guide = 'none') +
      scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
      theme(legend.position=unlist(ifelse(length(cl)==3, list(c(0.1, 0.8)), "none")),
            legend.title = element_blank())
    
    # dendrogram
    # try and synchronise topologies as far as possible
    try(dg <- ape::root(dg, outgroup=rev(get_taxa_name(gg.tr))[1:3], resolve.root=T))
    try(dg <- ape::root(dg, outgroup=rev(get_taxa_name(gg.tr))[1:2], resolve.root=T))
    try(dg <- ape::root(dg, outgroup=rev(get_taxa_name(gg.tr))[1], resolve.root=T))
    
    gg.ord <- ggtree(dg, branch.length = "none") %<+% tr.meta +
      geom_tiplab(aes(label=Label), fontface="italic", align=T, offset=-5.8, hjust=0, geom="label", label.size = NA, label.padding=unit(0.01, "lines")) +
      geom_tippoint(aes(colour=Colour, fill=Colour, shape=Location), size=3) +
      annotate("label", x=-Inf, y=length(dg$tip.label)-0.5, hjust=1,
               label=rf.lbl, parse=T, colour="dodgerblue2") +
      scale_shape_manual(values=c(`Ascension Island`=21,`Baffin Island`=24, `Falkland Islands`=25)) +
      scale_fill_manual(values=pcols2) +
      scale_colour_manual(values=pcols2) +
      theme(legend.position="none") + scale_x_reverse(expand = expansion(mult = c(0.01, 0)))
    
      arrangeGrob(gg.tr, gg.ord, ncol=2)
  })
  do.call(grid.arrange, gg.trees)
  return(gg.trees)
}

gg.phy1 <- phylosymbio.RF(ps.filt[[1]], tree)
gg.phy1.1 <- phylosymbio.RF(ps.filt[[1]], tol)
# weak evidence of phylosymbiosis for red algae with whole microbiome

gg.phy2 <- phylosymbio.RF(ps.oo.clade[[1]], tree)
gg.phy2.1 <- phylosymbio.RF(ps.oo.clade[[1]], tol)
# no evidence of phylosymbiosis with oomycetes

gg.phy3 <- phylosymbio.RF(ps.dia.clade[[1]], tree)
gg.phy3.2 <- phylosymbio.RF(ps.dia.clade[[1]], tol)
# excellent evidence of phylosymbiosis for red algae with diatoms!

do.call(grid.arrange, c(gg.phy1, gg.phy1.1, ncol=2))
do.call(grid.arrange, c(gg.phy2, gg.phy2.1, ncol=2))
do.call(grid.arrange, c(gg.phy3, gg.phy3.2, ncol=2))


# figures:
# Supplemental: everything
# break up into three pages
#gg.supp <- plot_grid(gg.phy1[[1]], gg.phy1[[2]], gg.phy1[[3]],
#                     gg.phy2[[1]], gg.phy2[[2]], gg.phy2[[3]], 
#                     gg.phy3[[1]], gg.phy3[[2]], gg.phy3[[3]], 
#                     rel_heights = c(0.43, 0.25, 0.12, 0.43, 0.25, 0.12, 0.43, 0.25, 0.12), hjust=0,
#                     ncol=1, labels=c(" A) All ASVs", NA, NA, " B) Oomycete ASVs", NA, NA, " C) Diatom ASVs", NA, NA))
gg.supp <- plot_grid(gg.phy1[[1]], gg.phy1[[2]], gg.phy1[[3]],
                     rel_heights = c(0.43, 0.25, 0.12), hjust=0,
                     ncol=1, labels=c(" A) All ASVs", NA, NA))
ggsave("phylosymbiosis_suppl1.pdf", gg.supp, width=12, height=12)
gg.supp <- plot_grid(gg.phy2[[1]], gg.phy2[[2]], gg.phy2[[3]],
                     rel_heights = c(0.43, 0.25, 0.12), hjust=0,
                     ncol=1, labels=c(" B) Oomycete ASVs", NA, NA))
ggsave("phylosymbiosis_suppl2.pdf", gg.supp, width=12, height=12)
gg.supp <- plot_grid(gg.phy3[[1]], gg.phy3[[2]], gg.phy3[[3]],
                     rel_heights = c(0.43, 0.25, 0.12), hjust=0,
                     ncol=1, labels=c(" C) Diatom ASVs", NA, NA))
ggsave("phylosymbiosis_suppl3.pdf", gg.supp, width=12, height=12)


# Main: only diatoms
# add Mantel plots for diatoms
pdf("phylosymbiosis_main.pdf", width=12, height=13)
#plot_grid(gg.phy3[[1]], gg.phy3[[2]], gg.phy3[[3]], 
#          plot_grid(plotlist=mantel6, ncol=3, align = 'hv', axis = 'tblr'), 
#          ncol=1, rel_heights=c(0.47, 0.33, 0.2, 0.5), labels=c("A", NA, NA, "B"))
plot_grid(plot_grid(plotlist=c(list(mantel3[[1]]), mantel6), ncol=4, align = 'hv', axis = 'tblr'),
          gg.phy3[[1]], gg.phy3[[2]], gg.phy3[[3]], 
          ncol=1, rel_heights=c(0.2, 0.43, 0.25, 0.12), labels=c("A", "B", NA, NA))
dev.off()

