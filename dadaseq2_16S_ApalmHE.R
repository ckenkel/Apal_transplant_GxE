#####################################

# installing packages: execute just once when first using a script:
#First bioconductor install method is commented out because it is for older than 3.5.0 versions of R. The second version of the install is for recent versions of R
#First I needed to install bioconductor for a new version of R. Now commented out. Only need to do this once.
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.14")

#source("https://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
#OR
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("phyloseq")
#
#INSTALL DADA2 1.4.1
#install.packages("~/Downloads/dada2-1.4.1", #to install from source, just indicate pkg download location
                 #repos = NULL,
                 #type = "source",
                 #dependencies = c("Depends", "Suggests","Imports"))
# Install DESeq2
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
  #BiocManager::install("DESeq2")


#####################################

library(dada2); packageVersion("dada2"); citation("dada2")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
#library(ShortRead)

#Prior to continuing with R analysis, complete the following steps in terminal
#Open Terminal window, navigate to folder with fastq.gz files 

#Trim filenames

#Also - some pre-trimming. Retain only PE reads that match amplicon primer. Remove reads containing Illumina sequencing adapters

#Still in terminal, type the following:
#ls *R1_001.fastq | cut -d '_' -f 1 > samples.list
### bbduk will keep paired reads in two separate files together if you use the in1= in2= out1= out2= function, both will be kept or both will be thrown away.
#for file in $(cat samples.list); do  mv 16S_${file}_*R1.fastq ${file}_R1.fastq; mv 16S_${file}_*R2.fastq ${file}_R2.fastq; done 
#for file in $(cat samples.list); do bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adaptors.fasta k=12 out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq ; done &>bbduk_NoIll.log


############Once above pre-processing steps completed, begin DADA2 anaylsis here

#Set path to trimmed, renamed fastq files
#path <- "C:/Users/User/OneDrive/Documents/UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/16S/SeqData" # CHANGE ME to the directory containing the fastq files after unzipping.
path <- "~/Dropbox/CarlsLab/ResearchProjects/NOAA_CRCP/PALMATA/NoIllAdapReads"
fns <- list.files(path)
fns

################################
##### Trimming/Filtering #######
################################

#forward and reverse fastq filenames
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs <- fastqs[grepl("_R1", fastqs)] # Just the forward read files
fnRs <- fastqs[grepl("_R2", fastqs)] # Just the reverse read files

#Specify the full path to the fnFs and fnRs
fnFs <- sort(list.files(path, pattern="_R1_NoIll.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_NoIll.fastq", full.names = TRUE))

#sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


sample.names
#########Visualize Raw data

#First, lets look at quality profile of R1 reads

plotQualityProfile(fnFs[c(1,2,3,4)]) 
plotQualityProfile(fnFs[c(137,138,139,140)])

#Then look at quality profile of R2 reads

plotQualityProfile(fnRs[c(1,2,3,4)])
plotQualityProfile(fnRs[c(137,138,139,140)])

# The reverse reads for Holland's 16S data are good. Quality is really high till 250bp. I don't see a need to trim...
#Going to argue for taking a few off to improve EE - say 245 bases for reverse

#The distribution of quality scores at each position is shown as a grey-scale heat map, with dark colors corresponding to higher frequency. 
#green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles.
#Recommend trimming where quality profile crashes - in this case, forward reads mostly fine up to 210; for reverse >160 bases it gets below 30; still should leave ~30-50bp overlap

#If using this workflow on your own data: Your reads must still overlap after truncation in order to merge them later! If you are using a less-overlapping primer set, your truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them.

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))

# Filter (NOT USING THIS ANYMORE, SEE FILTER STEP BELOW THAT ACTUALLY CUTS PRIMERS.)
# Used Emily's and the settings given in DaDa2 tutorial for Miseq 2x250bp
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,245),
                     #maxN=0, #DADA does not allow Ns
                     #maxEE=c(2,2), #c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     #truncQ=2,
                     #rm.phix=TRUE, #remove reads matching phiX genome
                     #matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     #compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE



# Filter (This step actually trims the primers, see trimLeft)
# Carly's settings for trimming 16S reads from Miseq 2x250bp
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,220), #seeing some read through - must cut out illumina adapter on far end of read
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(1,1), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     trimLeft=c(20,21), #N nucleotides to remove from the start of each read: 16S primers = F 19bp; R 20bp BUT looking at reads R1 runs into the Illumina adapter at the far end, so need to trim that 25 bases
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
tail(out)

#Hooray! Not many reads lost at all!! :)

#A word on Expected Errors vs a blanket quality threshold
#Take a simple example: a read of length two with quality scores Q3 and Q40, corresponding to error probabilities P=0.5 and P=0.0001. The base with Q3 is much more likely to have an error than the base with Q40 (0.5/0.0001 = 5,000 times more likely), so we can ignore the Q40 base to a good approximation. Consider a large sample of reads with (Q3, Q40), then approximately half of them will have an error (because of the P=0.5 from the Q2 base). We express this by saying that the expected number of errors in a read with quality scores (Q3, Q40) is 0.5.
#As this example shows, low Q scores (high error probabilities) dominate expected errors, but this information is lost by averaging if low Qs appear in a read with mostly high Q scores. This explains why expected errors is a much better indicator of read accuracy than average Q.
 
################################
##### Learn Error Rates #######
################################
#DADA2 learns its error model from the data itself by alternating estimation of the error rates and the composition of the sample until they converge on a jointly consistent solution (this is similar to the E-M algorithm)
#As in many optimization problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).


setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#why do values increase at Q40 in some plots?
#this artefact exists b/c in many sequencing runs there are almost no Q=40 bases. The loess smoothing hits the edge and a lack of observations, causing weird behavior. BUT as there essentially aren't (almost) any Q=40 bases to correct anyway and at the worst, the error rates are overestimated, so it's actually conservative for calling new variants



# DEREPLICATION HAS BEEN BUILT INTO LATER VERSIONS OF DADA2, SO IT IS UNNECESSARY TO DO THIS STEP.
################################
##### Dereplicate reads #######
################################
#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.
#derepFs <- derepFastq(filtFs, verbose=TRUE)
#derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
#names(derepFs) <- sample.names
#names(derepRs) <- sample.names


################################
##### Infer Sequence Variants #######
################################

#apply core sample inference algorithm to dereplicated data
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]
#> dadaFs[[1]]
#dada-class: object describing DADA2 denoising results
#355 sequence variants were inferred from 4223 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
#> dadaRs[[1]]
#dada-class: object describing DADA2 denoising results
#337 sequence variants were inferred from 2939 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[139]]
dadaRs[[139]]



################################
##### Merge paired reads #######
################################

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

###### REDO MAXMISMATCH=0 MORE AGRESSIVELY TRIM READS FIRST, CAN YOU TELL THE CHIMERA SCRIPT TO BE LESS PICKY

#merged paired data
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
summary((mergers[[1]]))
head(mergers[[139]])
summary((mergers[[139]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

################################
##### Construct sequence table #######
################################
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 252-254 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(240,264)] #again, being fairly conservative and still keeping a majority of pairs

table(nchar(getSequences(seqtab2)))
dim(seqtab2)

################################
##### Remove chimeras #######
################################
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)
# 0.9965134 #VERY few chimeras 
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 
#Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though)

################################
##### Track Read Stats #######
################################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab),rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtabRaw","inSilicoCut", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track, file="~/Dropbox/CarlsLab/ResearchProjects/NOAA_CRCP/PALMATA/ReadFilterStats_vTrim.csv", row.names=TRUE, quote = FALSE)


################################
##### Assign Taxonomy #######
################################

#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, to classify sequence variants taxonomically. 
#DADA2 provides a native implementation of the RDP's naive Bayesian classifier. The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences, and outputs the taxonomic assignments with at least minBoot bootstrap confidence.

taxa <- assignTaxonomy(seqtab.nochim, "~/Dropbox/CarlsLab/ResearchProjects/NOAA_CRCP/PALMATA/Silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#inspect taxa
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
tail(taxa.print)

#Now, save outputs so can come back to the analysis stage at a later point if desired
saveRDS(seqtab.nochim, file="~/Dropbox/CarlsLab/ResearchProjects/NOAA_CRCP/PALMATA/Apal_seqtabNochim_9May24.rds")
saveRDS(taxa, file="~/Dropbox/CarlsLab/ResearchProjects/NOAA_CRCP/PALMATA/Apal_taxonomy_9May24.rds")


#tax=readRDS("t1_taxonomy.rds")

#write.csv(tax, "t1_taxonomy.csv")

