library(ape)
library(phyloseq)
library(ggplot2)
library(ShortRead)
library(data.table)
library(plyr)
library(tidyr)
library(vegan)
library(pairwiseAdonis)
library(grid)
#library(DESeq2)
library(outliers)
library(gridExtra)#for grid.arrange
library(nlme)
library(reshape2)
library(dplyr)
#library(edgeR); packageVersion("edgeR")
#library(arrayQualityMetrics)
#library(microbiome) 
#library(eulerr)
library(dunn.test)
library(ggpubr)
library(purrr)
library(devtools)
#library(ANCOMBC) #- load this only when ready to run ancom as it modifies phyloseq


# Install packages if needed - only do this once for each

#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("ANCOMBC")
#install.packages("pairwiseAdonis")
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

BiocManager::install("phyloseq")


dir <- "~/Dropbox/CarlsLab/ResearchProjects/NOAA_CRCP/PALMATA/"
setwd(dir)


########### Set up for Phylosq ######################

seqtab.nochim=readRDS("Apal_seqtabNochim_9May24.rds")
# Make a copy of this that you can read so you can see how it is set up
#write.csv(seqtab.nochim, file = "C:/Users/elder/OneDrive/Documents/UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/16S/seqtab.nochim.csv")
#Now you can see that you need to split on a character to make your sample names match in the .rds.
rownames(seqtab.nochim) <- rownames(seqtab.nochim) %>% map(function(x) unlist(strsplit(x, "_"))[[1]])

# read in your ASVs
taxa=readRDS("Apal_taxonomy_9May24.rds")
#taxa=replace(taxa, taxa == " ", NA) 
#samdata=read.csv("samdata.csv")

# read in a .csv table with sample names. The names must match the sample names in your sequence file.
samdata=read.csv("sampledata.csv")
head(samdata)
rownames(samdata)=samdata$SampleID


#create phyloseq object

ps1=phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdata), tax_table(taxa))
ps1

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 8102 taxa and 141 samples ]
#sample_data() Sample Data:       [ 141 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 8102 taxa by 6 taxonomic ranks ]

sample_sums(ps1) #definitely some failed samples we will need to prune

############################ Check for mitochondria ###########################

mito=subset_taxa(ps1, Family =="Mitochondria") #8299 ASVs taxa left
mito
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 420 taxa and 141 samples ]
#sample_data() Sample Data:       [ 141 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 420 taxa by 6 taxonomic ranks ]
mito=subset_taxa(mito, Phylum !="Ochrophyta") #protists
mito
mito=subset_taxa(mito, Phylum !="Protalveolata") #protists
mito
mito=subset_taxa(mito, Phylum !="Apicomplexa") #protists
mito

#Merge by type
ps2<-merge_samples(mito,"Genotype")
ps2

pmito=plot_bar(ps2, fill="Family")+  geom_bar(aes( fill=Family), stat="identity", position="stack")+ theme_bw()+
  theme(text = element_text(size = 16))+
  theme(axis.text.x = element_text(face="bold", 
                                   size=10),
        axis.text.y = element_text(face="bold",  
                                   size=14))+ theme_bw() + theme(text = element_text(size = 16))+theme(legend.title = element_text(size = 14),
                                                                                                       legend.text = element_text(size= 14)) + facet_wrap(~Timepoint, scales = "free_x")
pmito

#Visualize and inspect
mitosequences <- Biostrings::DNAStringSet(taxa_names(ps2))
names(mitosequences) <- taxa_names(ps2)
ps2_mito <- merge_phyloseq(ps2, mitosequences)
ps2_mito
#use the ref_seq slots
rs <- refseq(ps2_mito)
rs # inspect
# Get strings with the full taxonomy for each OTU/ASV
tax <- tax_table(ps2_mito)
head(tax) # inspect
tax_strings <- apply(tax, 1, paste, collapse=";")
head(tax_strings) # inspect and BEAUTIFUL

# Create new names for the sequences in the refseq object; these will become
# the sequence headers in the fasta file. Here, I set these to be the OTU/ASV
# name and the tax string separated by one space (determined by the sep
# parameter)

taxa_names(ps2) <- paste0("ASV", seq(ntaxa(ps2)))
new_names <- paste(taxa_names(ps2), tax_strings, sep = " ")
head(new_names) # inspect
# Update the refeq names and save as a fasta file
names(rs) <- new_names
rs # inspect
#create fasta file with Biostrings!!!!
Biostrings::writeXStringSet(rs, "mitochondria.fasta")


######################### Back to our analysis  #######################################

#8102 taxa but what happens when
#filter out chloroplasts, mitochondria, protists and taxa counts 0 and under?

ps1=subset_taxa(ps1, (Class!="Chloroplast") | is.na(Class))
ps1 #still 8102 ASVs
ps1=subset_taxa(ps1, Family !="Mitochondria") 
ps1 #5398 ASVs
sample_sums(ps1) #have lost some samples here we will need to filter

#if interested in chloroplasts or mitochondrial sequences, we can subset and explore later (see above)

ps1=subset_taxa(ps1, Phylum !="Ochrophyta") #protists
ps1=subset_taxa(ps1, Phylum !="Protalveolata") #protists
ps1=subset_taxa(ps1, Phylum !="Apicomplexa") #protists
ps1
sample_sums(ps1)
#5398 taxa left - will record N reads remaining prior to rarefaction
write.csv(sample_sums(ps1), file = "sampsums_10Mar24.csv")

#prune samples with insufficient read data (<1000 reads per sample remaining)
ps1 <- prune_samples(sample_sums(ps1) >= 1000, ps1)
ps1 #down to 118 samples

##################### Rarefy  ####################

#find out how many reads in each sample before rarefying

plot(sort(sample_sums(ps1),TRUE),type="h", ylab="reads")#few high yield dudes, most low 

ps_rare <- phyloseq::rarefy_even_depth(ps1, rngseed = 123, replace = FALSE)
ps_rare #4431 taxa and 118 samples remain; 967 ASVs removed

plot(sort(sample_sums(ps_rare),TRUE),type="h", ylab="reads") 

# subset for later exploration of different time points
ps_rare_T24 = subset_samples(ps_rare, sample_data(ps_rare)$Timepoint != "T0")
ps_rare_T0 = subset_samples(ps_rare, sample_data(ps_rare)$Timepoint != "T24")

ps_rare_T24_site = subset_samples(ps_rare, !is.na(sample_data(ps_rare)$SiteTopBtm))
ps_rare_T24_site #56 samples

##################### T0 analysis  ####################

### What is initial community composition of T0 genets?

write.table(ps_rare_T0 %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              arrange(OTU) %>% rename(ASV = OTU) %>% 
              select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Sample, Abundance) %>%
              spread(Sample, Abundance), 
            file = "ps_rare_T0.relative_abundance.all.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#quick check of dominant taxa
T <- ps_rare_T0 %>% 
  tax_glom(., "Genus") %>% #change taxonomic selection
  transform_sample_counts(function(x)100* x / sum(x)) %>% psmelt() %>% 
  arrange(OTU) %>% rename(OTUsID = OTU) %>% 
  select(OTUsID, Genus, Sample, Abundance) %>% #change taxonomic selection
  spread(Sample, Abundance)

T$Mean <- rowMeans(T[, c(3:ncol(T))])

FAM <- T[, c("Genus", "Mean" ) ] #change taxonomic selection

#order from major to minor 
FAM <- FAM[order(dplyr::desc(FAM$Mean)),]
rownames(FAM) <- NULL
head(FAM,n=10) #this shows most abundant classes from top to bottom
summary(FAM)
sum(FAM$Mean)

tax_glom(ps_rare_T0, taxrank = 'Order') #29 Phyla, 50 classes, 121 Orders

#filter out taxa that don't appear more than 1x in more than 5 samples, so 5/49=0.1020.
# rarefied
filt_ps_rare_T0<-filter_taxa(ps_rare_T0, function(x) sum(x > 1) > (0.102*length(x)), TRUE)
filt_ps_rare_T0

#write out file of relative abundances by taxonomic length 
T <- filt_ps_rare_T0 %>% 
  tax_glom(., "Genus") %>% #change taxonomic selection
  transform_sample_counts(function(x)100* x / sum(x)) %>% psmelt() %>% 
  arrange(OTU) %>% rename(OTUsID = OTU) %>% 
  select(OTUsID, Genus, Sample, Abundance) %>% #change taxonomic selection
  spread(Sample, Abundance)

write.csv(T, "TaxonRelAbunByGenusSample_filt_ps_rare_T0.csv") #for a record or other analysis, like corrplots


#add a tree
random_tree=rtree(ntaxa(filt_ps_rare_T0), rooted=TRUE, tip.label=taxa_names(filt_ps_rare_T0))
#merge the tree
filt_ps_rare_T0=merge_phyloseq(filt_ps_rare_T0, random_tree)
filt_ps_rare_T0

## transform to relative abundance
trans_filt_ps_rare_T0<-transform_sample_counts(filt_ps_rare_T0, function(x) x /sum(x))
sample_sums(trans_filt_ps_rare_T0) #all sum to 1 as expected

##pretty plot color options
library(randomcoloR)
n <- 97
palette <- distinctColorPalette(n)
pie(rep(1, n), col=palette)

colorbfriend= c( "darkgreen", "deepskyblue4", "darkgoldenrod1", "aquamarine3", "coral2", "cadetblue", "plum4", "springgreen4", "aliceblue", "darkcyan", "darksalmon", "darkslateblue", "gold1", "honeydew2", "darkolivegreen4", "indianred3", "bisque1", "palevioletred4", "palegreen4", "maroon", "deepskyblue3", "darkred", "paleturquoise4", "burlywood2", "darkseagreen", "hotpink4", "lightcyan3", "darkslategray", "antiquewhite3", "goldenrod2", "royalblue4", "thistle", "forestgreen", "lightblue", "lightcoral", "lightseagreen", "peachpuff", "midnightblue", "rosybrown2", "mediumseagreen", "tan", "red4")
colorbfriend
####

p1=plot_bar(trans_filt_ps_rare_T0, fill="Class")+ scale_fill_manual(values=colorbfriend)
p1
phylumplot= p1+ geom_bar(aes(fill=Class), stat="identity", position="stack") 
phylumplot + theme_bw()+ 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.text.y = element_text(size=14)) + 
  theme(text = element_text(size = 16))+theme(legend.title = element_text(size = 14))+facet_grid(~Genotype, scale="free_x", drop=TRUE)

#are there any initial differences in diversity among genets (cross year and sym type differs)?

alph0 <- estimate_richness(ps_rare_T0, measures = c("Observed", "Shannon", "Simpson", "Chao1", "Fisher"))
alph0
write.csv(alph0, "alpha_diversity_ps_rare_T0.csv") #for a record or other analysis, like corrplots

#plot richness by measure
plot_richness(ps_rare_T0, x="Genotype")

p1= plot_richness(ps_rare_T0, x="Genotype", measures= c("Observed","Shannon"))+ geom_violin(fill="grey90")+ geom_jitter(width=0.15, alpha=0.8,size=2)+ theme(axis.text.x = element_text( 
  size=20),
  axis.text.y = element_text(  
    size=20))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=24))+
  ylab("Alpha Diversity Measure")+ xlab ("Genet") +theme_bw()

p1

#is the data normal?
hist(alph0$Shannon, main="Shannon diversity", breaks=10) #pretty good 

#Now we want to know if groups are different 

my_factors = data.frame(sample_data(ps_rare_T0))
#summary(ps_rare_T0)
set.seed(123)
aov_alpha_genets=aov(alph0$Shannon ~ my_factors$Genotype)
summary(aov_alpha_genets) #big difference among genets
TukeyHSD(aov_alpha_genets) #driven by X7

######### what about beta diversity among nursery genets?

genotype_colors<-c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon") 
#site_color_specific<-c("Eastern Dry Rocks"="mediumaquamarine", "Marker 32"="mediumseagreen", "W. Sambo"="paleturquoise1", "E. Sambo"="#0099c6", 
       #                "Maryland Shoals"="lightgoldenrod", "Daves Ledge"="lightsalmon", "Looe Key"="rosybrown1", "Big Pine"="light coral", "Raceway"="thistle" )

# plot
set.seed(123)
ordinate<- ordinate(filt_ps_rare_T0, method = "PCoA", distance ="wunifrac")
evals <- ordinate$values$Eigenvalues
filt_ps_rare_T0_betadiv=plot_ordination(filt_ps_rare_T0, ordinate, color = "Genotype") + #, shape = "Timepoint"
  scale_color_manual(values = genotype_colors)+
#  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Genet")+ geom_point(size = 4) + theme_bw() + theme(text = element_text(size = 14))+ 
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
filt_ps_rare_T0_betadiv

write.csv(ordinate$vectors, "beta_diversity_vectors_PCoA_wunifrac_filt_ps_rare_T0.csv") #for a record or other analysis, like corrplots


#use PERMANOVA to test for differences in beta diversity, first confirm data satisfy test assumptions

dfbray=phyloseq::distance(filt_ps_rare_T0, method="bray")
sampledf=data.frame(sample_data(filt_ps_rare_T0))
beta=betadisper(dfbray, sampledf$Genotype)
beta
permutest(beta) #p=0.594, even dispersion

adonis2(formula = dfbray ~ Genotype, data = sampledf)#p=0.001, yes sig differences between genets initially
#but, have 5 levels, therefore many pairwise contrasts to perform

#solution: use pairwiseAdonis
vegan::adonis2(dfbray ~ phyloseq::sample_data(filt_ps_rare_T0)$Genotype) #verifying same result as above

pairwise.adonis(dfbray, phyloseq::sample_data(filt_ps_rare_T0)$Genotype) #now pairwise test

### cool - have 3 groups of genets (1) XK and 14-5; (2) X9 and X5; and (3) X7 
#lets see what taxa differ among these groups using ANCOM-BC 
#see https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html

pseq_genus <- phyloseq::tax_glom(filt_ps_rare_T0, taxrank = "Genus")

output = ancombc2(data = pseq_genus, assay_name = "counts", tax_level = "Genus",
                  fix_formula = "T0_betadivGrp", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                  group = "T0_betadivGrp", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)) #should adjust these if have pre-determined hypotheses to test - see vignette

# Check for Structural zeros
tab_zero = output$zero_ind
head(tab_zero) 

#The primary output of the ANCOM-BC2 methodology identifies taxa with differential abundance 
#based on the chosen covariate. The results include: 1) log fold changes, 
#2) standard errors, 3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators 
#denoting whether the taxon is differentially abundant (TRUE) or not (FALSE), and 
#7) indicators denoting whether the taxon passed the sensitivity analysis (TRUE) or not (FALSE).

res_prim = output$res
head(res_prim)

write.csv(res_prim, "ancombc2_filt_ps_rare_T0_GeneraByBetaDivGrp.csv") #for a record or other analysis, like corrplots

#The ANCOM-BC2 methodology for multiple pairwise comparisons is designed to identify taxa 
#that exhibit differential abundance between any two groups within a set of three or more experimental groups, 
#all while maintaining control over the mdFDR.
#In the subsequent heatmap, each cell represents a log fold-change (in natural log) value. 

res_pair = output$res_pair
head(res_pair)
write.csv(res_pair, "ancombc2_filt_pairwise_ps_rare_T0_GeneraByBetaDivGrp.csv") 

df_fig_pair1 = res_pair %>%
  dplyr::filter(diff_T0_betadivGrpA_X7 == 1 |
                  diff_T0_betadivGrpC_XK145 == 1 | 
                  diff_T0_betadivGrpC_XK145_T0_betadivGrpA_X7 == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_T0_betadivGrpA_X7 == 1, 
                              round(lfc_T0_betadivGrpA_X7, 2), 0),
                lfc2 = ifelse(diff_T0_betadivGrpC_XK145 == 1, 
                              round(lfc_T0_betadivGrpC_XK145, 2), 0),
                lfc3 = ifelse(diff_T0_betadivGrpC_XK145_T0_betadivGrpA_X7 == 1, 
                              round(lfc_T0_betadivGrpC_XK145_T0_betadivGrpA_X7, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)


df_fig_pair1$group = recode(df_fig_pair1$group, 
                           `lfc1` = "X7 vs X5&X9",
                           `lfc2` = "XK&14-5 vs X5&X9",
                           `lfc3` = "XK&14-5 vs X7")

df_fig_pair1$group = factor(df_fig_pair1$group, 
                           levels = c("X7 vs X5&X9",
                                      "XK&14-5 vs X5&X9", 
                                      "XK&14-5 vs X7"))

lo = floor(min(df_fig_pair1$value))
up = ceiling(max(df_fig_pair1$value))
mid = (lo + up)/2

fig_pair = df_fig_pair1 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = 0, limit = c(lo, up), #set midpoint to mid if sensible
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes vs Reference") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
fig_pair

########################################### Okay - lets ask about general diffs between T0 and T24 - 
############################################ what changes when coral move to the field?

### What is community composition of all samples?

write.table(ps_rare %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              arrange(OTU) %>% rename(ASV = OTU) %>% 
              select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Sample, Abundance) %>%
              spread(Sample, Abundance), 
            file = "ps_rare_All.relative_abundance.all.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#quick check of dominant taxa
T <- ps_rare %>% 
  tax_glom(., "Genus") %>% #change taxonomic selection
  transform_sample_counts(function(x)100* x / sum(x)) %>% psmelt() %>% 
  arrange(OTU) %>% rename(OTUsID = OTU) %>% 
  select(OTUsID, Genus, Sample, Abundance) %>% #change taxonomic selection
  spread(Sample, Abundance)

T$Mean <- rowMeans(T[, c(3:ncol(T))])

FAM <- T[, c("Genus", "Mean" ) ] #change taxonomic selection

#order from major to minor 
FAM <- FAM[order(dplyr::desc(FAM$Mean)),]
rownames(FAM) <- NULL
head(FAM,n=10) #this shows most abundant classes from top to bottom
summary(FAM)
sum(FAM$Mean)

tax_glom(ps_rare, taxrank = 'Order') #29 Phyla, 50 classes, 121 Orders

#filter out taxa that don't appear more than 1x in more than 5 samples, so 5/118=0.04237288.
# rarefied
filt_ps_rare<-filter_taxa(ps_rare, function(x) sum(x > 1) > (0.0424*length(x)), TRUE)
filt_ps_rare #445 taxa remain

#write out file of relative abundances by taxonomic length 
T <- filt_ps_rare %>% 
  tax_glom(., "Genus") %>% #change taxonomic selection
  transform_sample_counts(function(x)100* x / sum(x)) %>% psmelt() %>% 
  arrange(OTU) %>% rename(OTUsID = OTU) %>% 
  select(OTUsID, Genus, Sample, Abundance) %>% #change taxonomic selection
  spread(Sample, Abundance)

write.csv(T, "TaxonRelAbunByGenusSample_filt_ps_rare_All.csv") #for a record or other analysis, like corrplots


#add a tree
random_tree=rtree(ntaxa(filt_ps_rare), rooted=TRUE, tip.label=taxa_names(filt_ps_rare))
#merge the tree
filt_ps_rare=merge_phyloseq(filt_ps_rare, random_tree)
filt_ps_rare

## transform to relative abundance
trans_filt_ps_rare<-transform_sample_counts(filt_ps_rare, function(x) x /sum(x))
sample_sums(trans_filt_ps_rare) #all sum to 1 as expected

## now barplot
p2=plot_bar(trans_filt_ps_rare, fill="Class")+ scale_fill_manual(values=colorbfriend) #create palette above
#p2
phylumplot= p2+ geom_bar(aes(fill=Class), stat="identity", position="stack") 
phylumplot + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),axis.title.x = element_blank(), axis.text.y = element_text(size=14)) + 
  theme(text = element_text(size = 16))+theme(legend.title = element_text(size = 14))+facet_grid(~Timepoint, scale="free_x", drop=TRUE)

#are there any initial differences in diversity among genets (cross year and sym type differs)?

alph <- estimate_richness(ps_rare, measures = c("Observed", "Shannon", "Simpson", "Chao1", "Fisher"))
alph
write.csv(alph, "alpha_diversity_ps_rare_All.csv") #for a record or other analysis, like corrplots

#plot richness by measure
plot_richness(ps_rare, x="Timepoint") #doesn't look too convincing

p3= plot_richness(ps_rare, x="Timepoint", measures= c("Observed","Shannon"))+ geom_violin(fill="grey90")+ geom_jitter(width=0.15, alpha=0.8,size=2)+ theme(axis.text.x = element_text( 
  size=20),
  axis.text.y = element_text(  
    size=20))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=24))+
  ylab("Alpha Diversity Measure")+ xlab ("Timepoint") +theme_bw()

p3

#is the data normal?
hist(alph$Shannon, main="Shannon diversity", breaks=10) #a bit left skewed, also repeated measures therefore non-independent 

#Now we want to know if groups are different 

my_factors = data.frame(sample_data(ps_rare))
#summary(ps_rare_T0)
set.seed(123)
kruskal.test(alph$Shannon ~ my_factors$Timepoint) #sig diff


######### what about beta diversity bewteen timepoints?

#genotype_colors<-c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon") 
site_color_specific<-c("Eastern Dry Rocks"="mediumaquamarine", "Marker 32"="mediumseagreen", "W. Sambo"="paleturquoise1", "E. Sambo"="#0099c6", 
               "Maryland Shoals"="lightgoldenrod", "Daves Ledge"="lightsalmon", "Looe Key"="rosybrown1", "Big Pine"="light coral", "Raceway"="thistle" )
site_colors<- c("lightcoral", "lightsalmon", "#0099c6", "darkseagreen1", "rosybrown1", "mediumseagreen", "lightgoldenrod", "thistle", "paleturquoise1")

# plot
set.seed(123)
ordinate<- ordinate(filt_ps_rare, method = "PCoA", distance ="wunifrac")
evals <- ordinate$values$Eigenvalues
filt_ps_rare_betadiv=plot_ordination(filt_ps_rare, ordinate, color = "Sites", shape = "Timepoint") + 
  scale_color_manual(values = site_colors)+
  #  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Site", shape = "Timepoint")+ geom_point(size = 4) + theme_bw() + theme(text = element_text(size = 14))+ 
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
filt_ps_rare_betadiv

write.csv(ordinate$vectors, "beta_diversity_vectors_PCoA_wunifrac_filt_ps_rare_All.csv") #for a record or other analysis, like corrplots


#use PERMANOVA to test for differences in beta diversity, first confirm data satisfy test assumptions

dfbray=phyloseq::distance(filt_ps_rare, method="bray")
sampledf=data.frame(sample_data(filt_ps_rare))
beta=betadisper(dfbray, sampledf$Timepoint)
beta
permutest(beta) #p<0.001, uneven dispersion - need nonparametric test

anosim(dfbray, sampledf$Timepoint, permutations = 1000) #yes, R is 0.8465 and significance=0.000999

### what are predominant taxa driving differences between nursery and field samples? insert this later

#lets see what taxa differ among these groups using ANCOM-BC 
#see https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html

pseq_genus <- phyloseq::tax_glom(filt_ps_rare, taxrank = "Genus")

#must run ancombc bc of only pairwise comparison, not multi-level
output2 = ancombc2(data = pseq_genus, assay_name = "counts", tax_level = "Genus",
                   fix_formula = "Timepoint", rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                   group = "Timepoint", struc_zero = FALSE, neg_lb = FALSE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = NULL) #should adjust these if have pre-determined hypotheses to test - see vignette


#The primary output of the ANCOM-BC2 methodology identifies taxa with differential abundance 
#based on the chosen covariate. The results include: 1) log fold changes, 
#2) standard errors, 3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators 
#denoting whether the taxon is differentially abundant (TRUE) or not (FALSE), and 
#7) indicators denoting whether the taxon passed the sensitivity analysis (TRUE) or not (FALSE).

res_prim = output2$res
head(res_prim)

write.csv(res_prim, "ancombc2_filt_ps_rare_All_GeneraByTimepoint.csv") #for a record or other analysis, like corrplots

df_fig_pair2 = res_prim %>%
  dplyr::filter(diff_TimepointT24 == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_TimepointT24 == 1, 
                              round(lfc_TimepointT24, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)


df_fig_pair2$group = recode(df_fig_pair2$group, 
                            `lfc1` = "T24 vs T0")

df_fig_pair2$group = factor(df_fig_pair2$group, 
                            levels = c("T24 vs T0"))

lo = floor(min(df_fig_pair2$value))
up = ceiling(max(df_fig_pair2$value))
mid = (lo + up)/2

fig_pair = df_fig_pair2 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = 0, limit = c(lo, up), #set midpoint to mid if sensible
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold change T24 vs T0") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
fig_pair



############################## Finally, what taxa are different among sites (and genets) at T24?

### What is initial community composition of T24 genets?

write.table(ps_rare_T24 %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              arrange(OTU) %>% rename(ASV = OTU) %>% 
              select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Sample, Abundance) %>%
              spread(Sample, Abundance), 
            file = "ps_rare_T24.relative_abundance.all.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

#quick check of dominant taxa
T <- ps_rare_T24 %>% 
  tax_glom(., "Class") %>% #change taxonomic selection
  transform_sample_counts(function(x)100* x / sum(x)) %>% psmelt() %>% 
  arrange(OTU) %>% rename(OTUsID = OTU) %>% 
  select(OTUsID, Class, Sample, Abundance) %>% #change taxonomic selection
  spread(Sample, Abundance)

T$Mean <- rowMeans(T[, c(3:ncol(T))])

FAM <- T[, c("Class", "Mean" ) ] #change taxonomic selection

#order from major to minor 
FAM <- FAM[order(dplyr::desc(FAM$Mean)),]
rownames(FAM) <- NULL
head(FAM,n=10) #this shows most abundant classes from top to bottom
summary(FAM)
sum(FAM$Mean)

tax_glom(ps_rare_T24, taxrank = 'Order') #29 Phyla, 50 classes, 121 Orders

#filter out taxa that don't appear more than 1x in more than 5 samples, so 5/69=0.0725.
# rarefied
filt_ps_rare_T24<-filter_taxa(ps_rare_T24, function(x) sum(x > 1) > (0.0725*length(x)), TRUE)
filt_ps_rare_T24 #260 taxa

#write out file of relative abundances by taxonomic length 
T <- filt_ps_rare_T24 %>% 
  tax_glom(., "Genus") %>% #change taxonomic selection
  transform_sample_counts(function(x)100* x / sum(x)) %>% psmelt() %>% 
  arrange(OTU) %>% rename(OTUsID = OTU) %>% 
  select(OTUsID, Genus, Sample, Abundance) %>% #change taxonomic selection
  spread(Sample, Abundance)

write.csv(T, "TaxonRelAbunByGenusSample_filt_ps_rare_T24.csv") #for a record or other analysis, like corrplots


#add a tree
random_tree=rtree(ntaxa(filt_ps_rare_T24), rooted=TRUE, tip.label=taxa_names(filt_ps_rare_T24))
#merge the tree
filt_ps_rare_T24=merge_phyloseq(filt_ps_rare_T24, random_tree)
filt_ps_rare_T24

## transform to relative abundance
trans_filt_ps_rare_T24<-transform_sample_counts(filt_ps_rare_T24, function(x) x /sum(x))
sample_sums(trans_filt_ps_rare_T24) #all sum to 1 as expected

### barplot by genet
p4=plot_bar(trans_filt_ps_rare_T24, fill="Class")+ scale_fill_manual(values=colorbfriend) #generate palette above
#p4
phylumplot= p4+ geom_bar(aes(fill=Class), stat="identity", position="stack") 
phylumplot + theme_bw()+ 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.text.y = element_text(size=14)) + 
  theme(text = element_text(size = 16))+theme(legend.title = element_text(size = 14))+facet_grid(~Genotype, scale="free_x", drop=TRUE)

### barplot by site
p4=plot_bar(trans_filt_ps_rare_T24, fill="Class")+ scale_fill_manual(values=colorbfriend) #generate palette above
#p4
phylumplot= p4+ geom_bar(aes(fill=Class), stat="identity", position="stack") 
phylumplot + theme_bw()+ 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.text.y = element_text(size=14)) + 
  theme(text = element_text(size = 16))+theme(legend.title = element_text(size = 14))+facet_grid(~Sites, scale="free_x", drop=TRUE)

#note that there is one gammaproteobacteria dominated sample per genet (all from EDR)

#are there any differences in diversity among genets or sites?

alph24 <- estimate_richness(ps_rare_T24, measures = c("Observed", "Shannon", "Simpson", "Chao1", "Fisher"))
alph24
write.csv(alph24, "alpha_diversity_ps_rare_T24.csv") #for a record or other analysis, like corrplots

#plot richness by measure
plot_richness(ps_rare_T24, x="Genotype")
plot_richness(ps_rare_T24, x="Sites")

p1= plot_richness(ps_rare_T24, x="Genotype", measures= c("Observed","Shannon"))+ geom_violin(fill="grey90")+ geom_jitter(width=0.15, alpha=0.8,size=2)+ theme(axis.text.x = element_text( 
  size=20),
  axis.text.y = element_text(  
    size=20))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=24))+
  ylab("Alpha Diversity Measure")+ xlab ("Genet") +theme_bw()

p1

p2= plot_richness(ps_rare_T24, x="Sites", measures= c("Observed","Shannon"))+ geom_violin(fill="grey90")+ geom_jitter(width=0.15, alpha=0.8,size=2)+ theme(axis.text.x = element_text( 
  size=20),
  axis.text.y = element_text(  
    size=20))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=24))+
  ylab("Alpha Diversity Measure")+ xlab ("Site") +theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p2

#is the data normal?
hist(alph24$Shannon, main="Shannon diversity", breaks=10) #very left skewed 

#Now we want to know if groups are different 

my_factors = data.frame(sample_data(ps_rare_T24))
#need nonparametric
kruskal.test(alph24$Shannon ~ my_factors$Genotype) #no sig diff among genets
kruskal.test(alph24$Shannon ~ my_factors$Sites) #no sig diff among sites


######### what about beta diversity genets or sites?

genotype_colors<-c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon") 
#site_color_specific<-c("Eastern Dry Rocks"="mediumaquamarine", "Marker 32"="mediumseagreen", "W. Sambo"="paleturquoise1", "E. Sambo"="#0099c6", 
#                "Maryland Shoals"="lightgoldenrod", "Daves Ledge"="lightsalmon", "Looe Key"="rosybrown1", "Big Pine"="light coral", "Raceway"="thistle" )
site_colors<- c("lightcoral", "lightsalmon", "#0099c6", "darkseagreen1", "rosybrown1", "mediumseagreen", "lightgoldenrod", "paleturquoise1")


# plot
set.seed(123)
ordinate<- ordinate(filt_ps_rare_T24, method = "PCoA", distance ="wunifrac")
evals <- ordinate$values$Eigenvalues
filt_ps_rare_T24_betadiv_genet=plot_ordination(filt_ps_rare_T24, ordinate, color = "Genotype") + #, shape = "Timepoint"
  scale_color_manual(values = genotype_colors)+
  #  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Genet")+ geom_point(size = 4) + theme_bw() + theme(text = element_text(size = 14))+ 
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
filt_ps_rare_T24_betadiv_genet

filt_ps_rare_T24_betadiv_site=plot_ordination(filt_ps_rare_T24, ordinate, color = "Sites") + #, shape = "Timepoint"
  scale_color_manual(values = site_colors)+
  #  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Site")+ geom_point(size = 4) + theme_bw() + theme(text = element_text(size = 14))+ 
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
filt_ps_rare_T24_betadiv_site

write.csv(ordinate$vectors, "beta_diversity_vectors_PCoA_wunifrac_filt_ps_rare_T24.csv") #for a record or other analysis, like corrplots


#Test for differences in beta diversity, first confirm data satisfy test assumptions

dfbray=phyloseq::distance(filt_ps_rare_T24, method="bray")
sampledf=data.frame(sample_data(filt_ps_rare_T24))

betaG=betadisper(dfbray, sampledf$Genotype)
betaG
permutest(betaG) #p=0.783, even dispersion among genets

adonis2(formula = dfbray ~ Genotype, data = sampledf)#p=0.282, no sig differences between genets 

betaS=betadisper(dfbray, sampledf$Sites)
betaS
permutest(betaS) #p=0.0.01, very uneven dispersion, must use anosim instead

ano<-anosim(dfbray, sampledf$Sites, permutations = 1000) #is significant difference among sites...

ano

summary(ano)
plot(ano)

str(ano)


################################################### last piece - whats the deal with the site specific diffs in performance?
########## try binning sites into top 3 and bottom 3 by survival and see what pops out

#filter out taxa that don't appear more than 1x in more than 5 samples, so 5/56=0.08929.
# rarefied
filt_ps_rare_T24_site<-filter_taxa(ps_rare_T24_site, function(x) sum(x > 1) > (0.08929*length(x)), TRUE)
filt_ps_rare_T24_site #198 taxa

#add a tree
random_tree=rtree(ntaxa(filt_ps_rare_T24_site), rooted=TRUE, tip.label=taxa_names(filt_ps_rare_T24_site))
#merge the tree
filt_ps_rare_T24_site=merge_phyloseq(filt_ps_rare_T24_site, random_tree)
filt_ps_rare_T24_site

## transform to relative abundance
trans_filt_ps_rare_T24_site<-transform_sample_counts(filt_ps_rare_T24_site, function(x) x /sum(x))
sample_sums(trans_filt_ps_rare_T24_site) #all sum to 1 as expected

#are there any differences in diversity among genets or sites?

alph24s <- estimate_richness(ps_rare_T24_site, measures = c("Observed", "Shannon", "Simpson", "Chao1", "Fisher"))
alph24s
write.csv(alph24s, "alpha_diversity_ps_rare_T24_site.csv") #for a record or other analysis, like corrplots

#plot richness by measure
plot_richness(ps_rare_T24_site, x="SiteTopBtm")

p1= plot_richness(ps_rare_T24_site, x="SiteTopBtm", measures= c("Observed","Shannon"))+ geom_violin(fill="grey90")+ geom_jitter(width=0.15, alpha=0.8,size=2)+ theme(axis.text.x = element_text( 
  size=20),
  axis.text.y = element_text(  
    size=20))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=24))+
  ylab("Alpha Diversity Measure")+ xlab ("Site Bins") +theme_bw()

p1

#is the data normal?
hist(alph24s$Shannon, main="Shannon diversity", breaks=10) #left skewed 

#Now we want to know if groups are different 

my_factors = data.frame(sample_data(ps_rare_T24_site))
#need nonparametric
kruskal.test(alph24s$Shannon ~ my_factors$SiteTopBtm) #no diff among top/btm sites


######### what about beta diversity between site bins

genotype_colors<-c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon") 
#site_color_specific<-c("Eastern Dry Rocks"="mediumaquamarine", "Marker 32"="mediumseagreen", "W. Sambo"="paleturquoise1", "E. Sambo"="#0099c6", 
#                "Maryland Shoals"="lightgoldenrod", "Daves Ledge"="lightsalmon", "Looe Key"="rosybrown1", "Big Pine"="light coral", "Raceway"="thistle" )
site_colors<- c("lightcoral", "lightsalmon", "#0099c6", "darkseagreen1", "rosybrown1", "mediumseagreen", "lightgoldenrod", "paleturquoise1")


# plot
set.seed(123)
ordinate<- ordinate(filt_ps_rare_T24_site, method = "PCoA", distance ="wunifrac")
evals <- ordinate$values$Eigenvalues
filt_ps_rare_T24_betadiv_sitebin=plot_ordination(filt_ps_rare_T24_site, ordinate, color = "SiteTopBtm") + #, shape = "Timepoint"
  scale_color_manual(values = genotype_colors)+
  #  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Site Bins")+ geom_point(size = 4) + theme_bw() + theme(text = element_text(size = 14))+ 
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14))
filt_ps_rare_T24_betadiv_sitebin


write.csv(ordinate$vectors, "beta_diversity_vectors_PCoA_wunifrac_filt_ps_rare_T24.csv") #for a record or other analysis, like corrplots


#Test for differences in beta diversity, first confirm data satisfy test assumptions

dfbray=phyloseq::distance(filt_ps_rare_T24_site, method="bray")
sampledf=data.frame(sample_data(filt_ps_rare_T24_site))

betaG=betadisper(dfbray, sampledf$SiteTopBtm)
betaG
permutest(betaG) #p=0.216, even dispersion among bins

adonis2(formula = dfbray ~ SiteTopBtm, data = sampledf)#p=0.018, !!difs in beta div between top btm sites

### what are predominant taxa driving differences between high/low survival sites?

#lets see what taxa differ among these groups using ANCOM-BC 
#see https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html

pseq_genus <- phyloseq::tax_glom(filt_ps_rare_T24_site, taxrank = "Genus") #none pass FDR

pseq_fam <- phyloseq::tax_glom(filt_ps_rare_T24_site, taxrank = "Family") #none pass FDR

pseq_ord <- phyloseq::tax_glom(filt_ps_rare_T24_site, taxrank = "Order") #none pass FDR

pseq_cls <- phyloseq::tax_glom(filt_ps_rare_T24_site, taxrank = "Class") #none pass FDR

#must run ancombc bc of only pairwise comparison, not multi-level
output3 = ancombc2(data = pseq_genus, assay_name = "counts", tax_level = "Genus",
                   fix_formula = "SiteTopBtm", rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                   group = "SiteTopBtm", struc_zero = FALSE, neg_lb = FALSE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = NULL) #should adjust these if have pre-determined hypotheses to test - see vignette


#The primary output of the ANCOM-BC2 methodology identifies taxa with differential abundance 
#based on the chosen covariate. The results include: 1) log fold changes, 
#2) standard errors, 3) test statistics, 4) p-values, 5) adjusted p-values, 6) indicators 
#denoting whether the taxon is differentially abundant (TRUE) or not (FALSE), and 
#7) indicators denoting whether the taxon passed the sensitivity analysis (TRUE) or not (FALSE).

res_prim = output3$res
head(res_prim)

write.csv(res_prim, "ancombc2_filt_ps_rare_T24_GeneraBySiteBin.csv") #for a record or other analysis, like corrplots
