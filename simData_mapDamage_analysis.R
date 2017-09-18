
#######################################################

#   Analyse damage patterns in simData files          #

#######################################################

#load packages
library(dplyr)
library(readr)
library(stringr)
library(reshape2)
library(tibble)
library(magrittr)
library(ggplot2)
library(scales)
library(pander)
library(gridExtra)

#set aesthetics
theme_set(theme_bw())
palette12 <- c("#FF3333", "#3333FF", "#009900", "#FF9900", "#990099", 
             "#33CCCC", "#66CC66", "#FFCC66", "#FF99CC", "#3399FF", 
             "#FF6666", "#9966FF")
palette15 <- c("#FF3333", "#3333FF", "#009900", "#FF9900", "#FF99CC", 
               "#3399FF", "#66CC66", "#FFCC66", "#FF6666", "#006699",
               "#336600", "#FFCC99", "#FF0066", "#9966FF", "#33CCCC")

#Create a blank theme to be applied to pie charts
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    plot.title=element_text(size=14, face="plain", hjust = 0.5)
  )

#What genomes are in the simData sets??
simData_setAbundance <- read_csv("SimulatedMetagenome.csv", col_names = FALSE, col_types = "c-n-c-cnc--", skip = 1) %>% 
  set_colnames(c("taxon", "abundance", "genus", "contaminant", "GC", "in_index")) %>% 
  filter(abundance < 1)

#From this can calculate number of reads expected from each taxon 
simData_setAbundance <- simData_setAbundance %>% mutate(No.reads = 1500000*abundance)
#no. expected aligned reads
#simData_setAbundance %>% filter(in_index == TRUE) %>% summarise(sum(.$No.reads)) #answer 5,535,500 wont work as list of what in_index incorrect

#Create a df with taxon and short taxon names for labelling
taxon_labels <- c(`Actinomyces oris strain T14V` = "Actinomyces oris", 
                  `Actinomyces sp. oral taxon 414 strain F0588` = "Actinomyces sp. oral taxon",
                  `Aggregatibacter actinomycetemcomitans strain 624` = "Aggregatibacter actinomycetemcomitans", 
                  `Aggregatibacter aphrophilus strain W10433` = "Aggregatibacter aphrophilus", 
                  `Agrobacterium tumefaciens strain A` = "Agrobacterium tumefaciens", 
                  `Bacillus subtilis BSn5` = "Bacillus subtilis", 
                  `Capnocytophaga haemolytica strain CCUG 32990` = "Capnocytophaga haemolytica", 
                  `Capnocytophaga sp. oral taxon 323 strain F0383` = "Capnocytophaga sp. oral taxon", 
                  `Fusobacterium nucleatum subsp. nucleatum ATCC 25586` = "Fusobacterium nucleatum subsp. nucleatum", 
                  `Fusobacterium nucleatum subsp. polymorphum strain ChDC F306` = "Fusobacterium nucleatum subsp. polymorphum", 
                  `Fusobacterium nucleatum subsp. vincentii 3_1_36A2` = "Fusobacterium nucleatum subsp. vincentii", 
                  `Leptotrichia buccalis DSM 113` = "Leptotrichia buccalis", 
                  `Leptotrichia sp. oral taxon 847` = "Leptotrichia sp. oral taxon", 
                  `Neisseria sicca strain FDAARGOS_2` = "Neisseria sicca", 
                  `Neisseria meningitidis MC58 chromosome` = "Neisseria meningitidis", 
                  `Porphyromonas gingivalis ATCC 33277 DNA` = "Porphyromonas gingivalis", 
                  `Prevotella dentalis DSM 3688` = "Prevotella dentalis", 
                  `Prevotella denticola F0289` = "Prevotella denticola", 
                  `Rothia dentocariosa ATCC 17931` = "Rothia dentocariosa", 
                  `Rothia mucilaginosa DNA complete genome strain: NUM-Rm6536` = "Rothia mucilaginosa", 
                  `Sphingomonas sp. MM-1` = "Sphingomonas sp. MM-1", 
                  `Staphylococcus epidermidis ATCC 12228` = "Staphylococcus epidermidis", 
                  `Streptococcus cristatus AS 1.3089` = "Streptococcus cristatus", 
                  `Streptococcus mitis B6` = "Streptococcus mitis", 
                  `Streptococcus mutans NN202DNA` = "Streptococcus mutans NN202DNA", 
                  `Streptococcus mutans UA159 chromosome` = "Streptococcus mutans UA159", 
                  `Streptococcus oralis Uo5` = "Streptococcus oralis", 
                  `Streptococcus sanguinis SK36` = "Streptococcus sanguinis", 
                  `Veillonella parvula DSM 2008` = "Veillonella parvula")

####plot abundance by taxon and by genus
#plot of abundance by taxon, labelling those that are present in the alignment index
simData_setAbundance %>% ggplot(aes(x=taxon, y=abundance, fill=contaminant)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1)) + 
  scale_x_discrete(labels=taxon_labels) + 
  scale_fill_manual(values = palette12, 
                    name = "", 
                    breaks=c("FALSE", "TRUE"), 
                    labels=c("Oral bacteria", "Common contaminant")) + 
  labs(x="", y="Abundance") + 
  annotate("text", x=c(1,3,6,9,12,14,15,16,18,19,22,24,26,28), 
           y=c(0.03,0.04,0.03,0.10,0.03,0.03,0.03,0.03,0.03,0.04,0.03,0.01,0.05,0.03), label="*", size=8)

#plot of abundance by genus
genusAbundance <- simData_setAbundance %>% select(-taxon, -contaminant, -GC, -in_index) %>% group_by(genus) %>% summarise_each(funs(sum))
genusAbundance <- simData_setAbundance %>% select(genus, contaminant) %>% unique() %>% left_join(genusAbundance)
genusAbundance %>% ggplot(aes(x=genus, y=abundance, fill=contaminant)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1)) + 
  scale_x_discrete(labels=taxon_labels) + 
  scale_fill_manual(values = palette12, 
                    name = "", 
                    breaks=c("FALSE", "TRUE"), 
                    labels=c("Oral bacteria", "Common contaminant")) + 
  labs(x="Genus", y="Abundance") + 
  annotate("text", x=c(1,2,4,6:11,13,14), y=c(0.097,0.08,0.03,0.147,0.06,0.06,0.03,0.06,0.041,0.03,0.195), label="*", size=8)

############# Varied contaminant and damage pattern ###############

#these are the abundances for the endogenous bacterial genomes
simData_endoAbundance <- read_delim("inputGenomes/bact/list", delim = "\t", col_names = FALSE, col_types = "cn") %>% 
  set_colnames(c("species", "abundance")) %>% mutate(source = "endogenous")
simData_endoAbundance$species <- gsub('.fna', '', simData_endoAbundance$species)

#abundances for the env and cont sequences were 0.5
envSpecies <- c("A_tumefaciens", "B_subtilis")
labSpecies <- c("Sphingomonas_sp", "S_epidermidis")
contAbundance <- c(0.5,0.5)
simData_labAbundance <- data_frame(labSpecies, contAbundance) %>% set_colnames(c("species", "abundance")) %>% mutate(source = "lab")
simData_envAbundance <- data_frame(envSpecies, contAbundance) %>% set_colnames(c("species", "abundance")) %>% mutate(source = "environment")

#calcualte abundance endogenous microbial genomes when contaminant low, lowMod, mod and high
simData_variedCont <- simData_endoAbundance %>% 
  mutate(lowCont = abundance*0.85, lowModCont = abundance*0.6, modCont = abundance*0.35, highCont = abundance*0.1)
#do same for envAbundance and labAbundances
simData_variedCont <- simData_envAbundance %>% 
  mutate(lowCont = abundance*0.1, lowModCont = abundance*0.35, modCont = abundance*0.6, highCont = abundance*0.85) %>% 
  bind_rows(simData_variedCont)
simData_variedCont <- simData_labAbundance %>% 
  mutate(lowCont = abundance*0.05, lowModCont = abundance*0.05, modCont = abundance*0.05, highCont = abundance*0.05) %>% 
  bind_rows(simData_variedCont) %>% select(-abundance)

contamination_labels <- c(`lowCont` = "Low Contamination\n(85% endogenous, 5% laboratory, 10% environmental", 
                          `lowModCont` = "Low - Moderate Contamination\n(60% endogenous, 5% laboratory, 35% environmental)", 
                          `modCont` = "Moderate Contamination\n(35% endogenous, 5% laboratory, 35% environmental)", 
                          `highCont` = "High Contamination\n(10% endogenous, 5% laboratory, 35% environmental")
  
#melt this dataframe and plot with facet_wrap to show how proportion of species changes in these simulated dataSets
simData_variedCont %>% melt(id.vars = c("species", "source")) %>% 
  ggplot(aes(x=species, y=value, fill=source)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(~variable, scales = "free_y", labeller = as_labeller(contamination_labels)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1)) + 
  scale_fill_manual(values = c("#FF3333", "#009900", "#3333FF"), 
                    name = "Source", 
                    breaks=c("endogenous", "environment","lab"), 
                    labels=c("Endogenous oral bacteria", "Environmental contaminant", "Laboratory contaminant")) +
  labs(x="", y="Abundance")

################ alignment data ##########################

#Read in txt file with simulated count data and assign to object
#fastqCount <- read_delim("data/fastq_read_count.txt", delim = "\t", 
#                              skip = 1, col_names = FALSE) %>%
#  set_colnames(c("fileName", "fastqCount"))

#import mapped_read_count.txt data
mapCount <- read.csv(file="bwaMapData/bwa_map_count.txt", 
                          sep="", skip = 1, header = FALSE, 
                          col.names = c("count", "MAPQ")) %>% 
  mutate(bam = grepl("bam", count), fileNo = cumsum(bam))

#extract the fileInfo and fileNo information
fileInfo <- mapCount[grep("bam", mapCount$count),] %>% select(count, fileNo)

#rejoin fileInfo as a separate variable
mapCount <- mapCount %>% left_join(fileInfo, by = "fileNo") %>% 
  select(-bam, -fileNo) %>% 
  set_colnames(c("count", "MAPQ", "fileName")) %>% 
  filter(MAPQ != "NA")
#convert factor variables to character
mapCount <- mapCount %>% mutate_if(is.factor, as.character)
#convert splitCount from character to numeric
mapCount$count <- as.numeric(mapCount$count)

#extract data for bwa.bam and rmdup.bam files and assign to separate objects
bwaCount <- mapCount %>% subset(grepl("_bwa.bam", fileName))
rmdupCount <- mapCount %>% subset(grepl("_rmdup.bam", fileName))
#edit fileNames to be identical to fastq
bwaCount$fileName <- gsub('_bwa.bam', '', bwaCount$fileName)
rmdupCount$fileName <- gsub('_rmdup.bam', '', rmdupCount$fileName)

######Create table summarising total fastq, bwa.bam and rmdup.bam for each sample#####
totalRmdupCount <- rmdupCount %>% 
  select(-MAPQ) %>% 
  group_by(fileName) %>% 
  summarise_each(funs(sum)) %>% set_colnames(c("fileName", "nonDupCount"))
totalBwaCount <- bwaCount %>% 
  select(-MAPQ) %>% 
  group_by(fileName) %>% 
  summarise_each(funs(sum)) %>% set_colnames(c("fileName", "alnCount"))
totalCount <- left_join(totalBwaCount, totalRmdupCount, by = "fileName") %>% 
  mutate(fastqCount = 1500000)
  
##From this calculate the total number of duplicate reads identified in each sample and add to totalCount
totalCount <- totalCount %>% mutate(dupCount = alnCount - nonDupCount)
totalCountTable <- totalCount[,c(1,4,2,5,3)]
names(totalCountTable) <- c("fileName", "Processed", "Aligned", "Duplicate", "Remaining")
#totalCountTable %>% pander(caption = "Summary of Count Data for each sample")

#Several of the files are replicate (i.e. almost no difference between 30bp_0.1_undamaged and 30bp_0.5_undamaged and 30bp_Real-profile_undamaged)
#Suggest removing this replicated data, at least for initial count plots

#generate a function edit simDataFileNames to generate more concise and readable names

#USE: df$fileName <- editFileNames(df)

editSimDataFileNames <- function(x){
  x$fileName <- gsub('.bam', '', x$fileName)
  x$fileName <- gsub('\\.b', '_undamaged', x$fileName)
  x$fileName <- gsub('_d', '_damaged', x$fileName)
  x$fileName <- gsub('-damage-no-adapters', '', x$fileName)
  x$fileName <- gsub('Real-damage-profile-', 'Real-profile', x$fileName)
  x$fileName <- gsub('no-adapters', '', x$fileName)
  x$fileName <- gsub('-damage-ACAD-adapters', '', x$fileName)
  x$fileName <- gsub('0-', '0.', x$fileName)
  x$fileName <- gsub('_noDamage_damaged', '_undamaged', x$fileName)
  x$fileName <- gsub('_withDamage', '', x$fileName)
  x$fileName <- gsub('_Real-profile_undamaged', '_undamaged', x$fileName)
}

#apply function to totalCount df's
totalCount$fileName <- editSimDataFileNames(totalCount)

#remove repetitious data
totalCount <- filter(totalCount, !grepl("_..._undamaged", fileName))

#extract length from fileName and add as separate variable
#Extract length from fileName 
extractLength <- function(x){
  x %>% mutate(length = str_extract(fileName,  "(30bp|50bp|70bp|90bp|Empirical)"))
}

#Apply function to totalCount
totalCount <- extractLength(totalCount)

###Plot number of reads aligned and unaligned for each sample with all endogenous content
totalCount %>% select(fileName, fastqCount, alnCount, length) %>% 
  mutate(Unaligned = fastqCount - alnCount) %>% 
  select(-fastqCount) %>% 
  melt(id.vars = c("fileName", "length"), variable.name = "group", value.name = "count") %>% 
  filter(!grepl("endo", fileName)) %>% 
  ggplot(aes(x=fileName, y=count, fill=group)) + 
  geom_bar(width = 1, stat = "identity", position = position_stack(reverse = TRUE), colour = "white") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1)) + 
  scale_y_continuous(labels = scales::comma) + 
  scale_fill_manual(values = c("#FF6666", "#3399FF"), name = "", labels = c("Aligned", "Unaligned")) + 
  facet_grid(.~length, scales = "free_x", space = "free_x")

##Plot number of reads aligned and unaligned for sample with variation in endogenous content
totalCount %>% select(fileName, fastqCount, alnCount, length) %>% 
  mutate(Unaligned = fastqCount - alnCount) %>% 
  select(-fastqCount) %>% 
  melt(id.vars = c("fileName", "length"), variable.name = "group", value.name = "count") %>% 
  filter(grepl("endo", fileName)) %>% 
  ggplot(aes(x=fileName, y=count, fill=group)) + 
  geom_bar(width = 1, stat = "identity", position = position_stack(reverse = TRUE), colour = "white") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 1)) + 
  scale_y_continuous(labels = scales::comma) + 
  scale_fill_manual(values = c("#FF6666", "#3399FF"), name = "", labels = c("Aligned", "Unaligned")) + 
  facet_grid(.~length, scales = "free_x", space = "free_x")


#Calculate % of reads aligned/unaligned, duplicate/nonduplicate and add to totalCountTable
totalCount <- totalCount %>% 
  mutate(propAln = (alnCount/fastqCount)*100, 
         propUnAln = ((fastqCount-alnCount)/fastqCount)*100, 
         propDup = (dupCount/alnCount)*100, 
         propNonDup = (nonDupCount/alnCount)*100)
#round % to 2 decimal places
totalCount[,7:10] <- round(totalCount[,7:10],2)

totalCount %>% 
  select(fileName, propAln, propUnAln) %>% 
  melt(id.vars = "fileName") %>% 
  ggplot(aes(x="", y=value, fill=variable)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0) + 
  blank_theme +
  scale_fill_manual(values = c("#FF6666", "#3399FF"),
                    name = "",
                    breaks=c("propAln", "propUnAln"), 
                    labels=c("Aligned", "Not aligned")) +
  theme(axis.text.x = element_blank(), strip.text = element_text(size = 10), legend.position = "bottom") + 
  geom_text(aes(label = value), position = position_stack(vjust = 0.5)) +
  facet_wrap(~fileName)

################# split count Data ########################

#Extract counts for split.bam files from mapCount
splitCount <- mapCount %>% subset(grepl("_split.bam", fileName))
#edit fileName to remove '_split.bam', and make more concise
splitCount$fileName <- gsub('_split.bam', '', splitCount$fileName)


#split out genome names
extractGenomeName <- function(x){
  x %>% mutate(genome = str_extract(fileName,  "(A.actinomyces|A.oris|A.parvulum|B.subtilis|C.gracilis|C.sporogenes|E.saphenum|E.sulci|F.nucleatum|H.influenza|K.flavida|L.buccalis|M.neoaurum|N.meningitidis|N.sicca|P.denticola|P.florescens|P.gingivalis|P.intermedia|P.propionicum|R.dentocariosa|S.epidermidis|S.gordonii|S.mitis|S.mutans|S.roseum|S.sanguinis|T.denticola|T.forsythia|M.oralis)"), 
               fileName = str_replace(fileName, "_(A.actinomyces|A.oris|A.parvulum|B.subtilis|C.gracilis|C.sporogenes|E.saphenum|E.sulci|F.nucleatum|H.influenza|K.flavida|L.buccalis|M.neoaurum|N.meningitidis|N.sicca|P.denticola|P.florescens|P.gingivalis|P.intermedia|P.propionicum|R.dentocariosa|S.epidermidis|S.gordonii|S.mitis|S.mutans|S.roseum|S.sanguinis|T.denticola|T.forsythia|M.oralis)", ""))
}

#summarise total count for each split.bam file
totalSplitCount <- splitCount %>% 
  select(-MAPQ) %>% 
  group_by(fileName) %>% 
  summarise_each(funs(sum))

totalSplitCount <- extractGenomeName(totalSplitCount)
totalSplitCount$genome <- gsub('A.actinomyces', 'A.actinomycetemcomitans', totalSplitCount$genome)
totalSplitCount$fileName <- editSimDataFileNames(totalSplitCount)

#plot abundance of each species for simulated data
#simData %>% select(abundance, genus) %>% 
#  group_by(genus) %>% 
#  summarise_each(funs(sum)) %>% 
#  ggplot(aes(x="", y=abundance, fill=genus)) + 
#  geom_bar(stat = "identity", colour="white", position = position_dodge()) + 
#  scale_fill_manual(values = palette15) + 
#  labs(x="Simulated data", y="Abundance", fill="Genus")

#plot count for each genome aligned, colouring by genus
#totalSplitCount %>% filter(length == "30bp") %>% 
#  ggplot(aes(x=fileName, y=count, fill=genome)) + 
#  geom_bar(stat = "identity", colour="white", position = position_stack()) + 
  #scale_fill_manual(values = palette15) + 
#  labs(x="", y="Read count", fill="Genome") + 
#  scale_y_continuous(labels = scales::comma)

###Too many files. Can't see anything in plot

##### Plot split count as pie-charts #####

#Add totalCount for each fileName and calculate proportion for each genome
totalSplitCount <- totalCount %>% 
  select(fileName, nonDupCount) %>% 
  left_join(totalSplitCount, by = "fileName") %>% 
  mutate(propAln = count/nonDupCount)
#round % to 4 decimal places
totalSplitCount[,5] <- round(totalSplitCount[,5],4)

palette30 <- c("#E41A1C", "#546C9D", "#459E6F", "#718074", "#B75D70", "#FF9007", "#FFFA31", "#B8782A", "#D87085", "#CC8BAD", 
               "#B43547", "#3983AC", "#4BAB51", "#85658D", "#D46A43", "#FFB315", "#E9D630", "#AB5832", "#EE7CAF", "#B292A3", 
               "#845172", "#3F908E", "#5C9A5C", "#9B4F9D", "#F07816", "#FFD723", "#D0A72D", "#C1645C", "#E685B8", "#999999")

#this is to plot proportion of reads aligning without any filtering.
totalSplitCount %>% filter(grepl("30bp_", fileName)) %>% 
  split(f = .$fileName) %>% 
  lapply(function(x){x %>% 
  ggplot(aes(x="", y=propAln, fill=genome)) + 
  geom_bar(stat = "identity") + 
  coord_polar("y", start = 0) + 
  blank_theme +
  scale_fill_manual(values = palette30, drop = TRUE, name = "") + 
  theme(axis.text.x = element_blank(), strip.text = element_text(size = 10), 
        legend.position = "bottom") + 
      geom_text(aes(label = propAln), position = position_stack(vjust = 0.5)) +
      ggtitle(x$fileName)})

#Re-create splitCount after MAPQ filtering of -q 25
Q25splitCount <- splitCount %>% filter(MAPQ > 24)
#summarise total count for each split.bam file
Q25totalSplitCount <- Q25splitCount %>% 
  select(-MAPQ) %>% 
  group_by(fileName) %>% 
  summarise_each(funs(sum))
#split out genome names
Q25totalSplitCount <- extractGenomeName(Q25totalSplitCount)
Q25totalSplitCount$genome <- gsub('A.actinomyces', 'A.actinomycetemcomitans', Q25totalSplitCount$genome)
#edit fileNames
Q25totalSplitCount$fileName <- editSimDataFileNames(Q25totalSplitCount)
#remove files with fewer than 500 reads
Q25totalSplitCount <- Q25totalSplitCount %>% filter(count >= 500)

#Create df with total remaining count for each file
totalFilteredCount <- Q25totalSplitCount %>% 
  split(f = .$fileName) %>% 
  lapply(function(x){x %>% 
      select(count) %>% 
      summarise_each(funs(sum))}) %>% 
  bind_rows(.id = "fileName") %>% 
  set_colnames(c("fileName", "filteredCount"))

#bind this information to Q25totalSplitCount and calculate proportion aligned
Q25totalSplitCount <- left_join(Q25totalSplitCount, totalFilteredCount, by = "fileName") %>% 
  mutate(propAln = count/filteredCount)
#Generate plotting labels
Q25totalSplitCount <- Q25totalSplitCount %>% mutate(label = round(propAln*100, digits = 1))

#plot new proportion aligned?? will this work since some have been filtered out??
Q25totalSplitCount %>% filter(grepl("30bp_", fileName)) %>% 
  split(f = .$fileName) %>% 
  lapply(function(x){x %>% 
      ggplot(aes(x="", y=propAln, fill=genome)) + 
      geom_bar(stat = "identity") + 
      coord_polar("y", start = 0) + 
      blank_theme +
      scale_fill_manual(values = palette30, drop = TRUE, name = "") + 
      theme(axis.text.x = element_blank(), strip.text = element_text(size = 10), 
            legend.position = "bottom") + 
      geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
      ggtitle(x$fileName)})

### Could also do a plot showing what is not expected to align (e.g. H.influenza)


################### plot fragment lengths #################

#create a list of lgdist files
lgDistFiles <- list.files("mapDamageData", pattern = "lgdistribution.txt", 
                          full.names = TRUE, recursive = TRUE)
#read-in data from each text file and bind into a data frame (modified as some files are empty)
simLengthData <- lapply(lgDistFiles,function(x){
  x <- try(read.table(paste(x,sep=""), head=FALSE, stringsAsFactors=FALSE, sep="\t", 
                      col.names = c("std", "length", "occ"), skip = 4)) %>% mutate(fileName = x)
  if(inherits(x, "try-error"))
    return(NULL)
  else
    return(x)
}) %>% bind_rows

#edit fileName to remove unnecessary information
simLengthData$fileName <- gsub('mapDamageData/', '', simLengthData$fileName)
simLengthData$fileName <- gsub('_split.bam/lgdistribution.txt', '', simLengthData$fileName)

#extract genomes from fileName
simLengthData <- extractGenomeName(simLengthData)
simLengthData$genome <- gsub('A.actinomyces', 'A.actinomycetemcomitans', simLengthData$genome)

#Collate lengths for each sample & plot distributions
simLengthData %>% split(f = .$fileName) %>%
  lapply(function(x){x %>% select(length, occ) %>% 
      group_by(length) %>% 
      summarise_each(funs(sum))}) %>% 
  bind_rows(.id = "fileName") %>% filter(grepl("Empirical", fileName)) %>% 
  ggplot(aes(x=length, y=occ, colour=fileName)) + 
  geom_line() + 
  theme_bw() + 
  labs(x="Read length", y="Number of reads") + 
  facet_wrap(~fileName) + 
  guides(colour=FALSE) + 
  scale_y_continuous(labels = scales::comma) + 
  ggtitle("Distribution of fragment lengths for each simulated file")

################# plot miscoding lesions ###################

#Create a vector with desired column names
subDataColNames <- (c("Chr", "End", "Std", "Pos", "A", "C", "G", "T", "Total", 
                      "GtoA", "CtoT", "AtoG", "TtoC", "AtoC", "AtoT", "CtoG", 
                      "CtoA", "TtoG", "TtoA", "GtoC", "GtoT"))

# create list of all .txt files in folder 
ntSubFiles <- list.files("mapDamageData", pattern = "misincorporation.txt",
                         full.names = TRUE, recursive = TRUE)

# read-in files and bind into a data frame that include the FileName
ntSubData <- ntSubFiles %>% lapply(function(x){
  read_delim(x, delim = "\t", skip = 4, col_names = FALSE, col_type = "cccnnnnnnnnnnnnnnnnnn---------") %>% 
    set_colnames(subDataColNames) %>% filter(Total > 0) %>% filter(Pos < 26) %>%
    mutate(fileName = x) %>% select(-Chr)
}) %>%
  bind_rows

#Edit FileName to include only the sampleID and Genome
ntSubData$fileName <- gsub('mapDamageData/', '', ntSubData$fileName)
ntSubData$fileName <- gsub('_split.bam/misincorporation.txt', '', ntSubData$fileName)


#Extract genome ID from the split fileName, and make fileNames mroe concise
ntSubData <- extractGenomeName(ntSubData)
ntSubData$genome <- gsub('A.actinomyces', 'A.actinomycetemcomitans', ntSubData$genome)
ntSubData$fileName <- editSimDataFileNames(ntSubData)

#Collate counts for pos and neg strand
ntSubData <- ntSubData %>% split(f = .$fileName) %>%
  lapply(function(x){x %>% select(-fileName) %>% split(f = .$End) %>% lapply(function(z){
    z %>% select(-End) %>% split(f = .$genome) %>% 
      lapply(function(a){a %>% select(-genome, -Std) %>% group_by(Pos) %>% summarise_each(funs(sum))}) %>% 
      bind_rows(.id = "genome")}) %>% 
      bind_rows(.id = "End")}) %>% bind_rows(.id = "fileName")


################## Graphing ntSubData ##############

#Add totalFilteredCount to ntSubData
ntSubData <- Q25totalSplitCount %>% select(fileName, genome, count) %>% 
  left_join(ntSubData, by = c("fileName", "genome"))

# create graphing function
ntSub.graph <- function(df, na.rm = TRUE, ...){
  
  # Specify sampleID
  sampleID <- unique(df$fileName)
  
  # create list of genomeID's in data to loop over 
  genomeID_list <- unique(df$genome)
  
  # create list of counts in data to loop over
  genomeCount_list <- unique(df$count)
  
  # create for loop to split data based on sampleID 
  for (i in seq_along(genomeID_list)) {
    
    # create object to store 5p data
    SubFreq_5p <- subset(df, df$genome==genomeID_list[i]) %>% filter(End == "5p") %>% 
      select(-End, -genome, -fileName, -count) %>% 
      #group_by(Pos) %>% summarise_each(funs(sum)) %>% 
      mutate(GtoA = GtoA/G, CtoT = CtoT/C, AtoG = AtoG/A, TtoC = TtoC/T, 
             AtoC = AtoC/A, AtoT = AtoT/A, CtoG = CtoG/C, CtoA = CtoA/C, 
             TtoG = TtoG/T, TtoA = TtoA/T, GtoC = GtoC/G, GtoT = GtoT/G) %>% 
      select(-A, -C, -G, -T, -Total) %>% 
      melt(id.vars = c("Pos"), variable.name = "substitution", value.name = "Frequency")
    
    #create object to store 3p data
    SubFreq_3p <- subset(df, df$genome==genomeID_list[i]) %>% filter(End == "3p") %>% 
      select(-End, -genome, -fileName, -count) %>% 
      #group_by(Pos) %>% summarise_each(funs(sum)) %>% 
      mutate(GtoA = GtoA/G, CtoT = CtoT/C, AtoG = AtoG/A, TtoC = TtoC/T, 
             AtoC = AtoC/A, AtoT = AtoT/A, CtoG = CtoG/C, CtoA = CtoA/C, 
             TtoG = TtoG/T, TtoA = TtoA/T, GtoC = GtoC/G, GtoT = GtoT/G) %>% 
      select(-A, -C, -G, -T, -Total) %>% 
      melt(id.vars = c("Pos"), variable.name = "substitution", value.name = "Frequency")
    
    #plot to object, 5p data
    plot5pData <- SubFreq_5p %>% ggplot(aes(x=Pos, y=Frequency, colour=substitution)) + 
      geom_line() + 
      theme_bw() + 
      scale_y_continuous(limits = c(0, 0.55), position = "left") + 
      ylab("Substitution Frequency") + 
      xlab("Position from the 5' end") + 
      scale_colour_manual(values=c(palette15), name = "Substitution") +
      theme(legend.position = "none")
    
    #plot to object, 3p data 
    plot3pData <- SubFreq_3p %>% ggplot(aes(x=Pos, y=Frequency, colour=substitution)) + geom_line() + theme_bw() +
      theme(axis.title.y=element_blank()) + scale_x_reverse() + 
      scale_y_continuous(limits = c(0, 0.55), position = "right") +
      scale_colour_manual(values=c(palette15), name = "Substitution") +
      ylab("Substitution Frequency") + xlab("Position from the 3' end")
    
    #print plots
    grid.arrange(plot5pData, plot3pData, ncol=2, widths=c(0.8,1), 
                 top = (paste(genomeID_list[i], '-',paste(sampleID,'(No. Reads = ', paste(genomeCount_list[i],')')))))
    
    #End loop
  }
}

#plot
ntSubData %>% filter(fileName == "30bp_0.1_damaged") %>% ntSub.graph()

#Can then plot rest. May be best to split by fileName and run function through lapply

############## plot overall C>T and G>A rate ####################

#Generate a list of sub_freq.txt files
subFreq.Files <- list.files("mapDamageData/", 
                            pattern = "_freq.txt", 
                            full.names = TRUE, 
                            recursive = TRUE)

#turn off scientific notation for numbers 
##(to prevent any values in frequency text that are in sci notation being incorrectly imported)
options(scipen = 999)

#Write a function to load data and edit the fileName to remove unncessesary information
loadMisincorpData <- function(x){
  #load data
  data <- x %>% lapply(function(z){z %>% read_delim(delim = "\t", 
                                                    skip=1, col_names = FALSE, 
                                                    col_types = cols("i", "c"), 
                                                    n_max = 5) %>% 
      set_colnames(c("pos", "freq")) %>% 
      mutate(fileName = z)}) %>% 
    bind_rows()
  #convert character to numeric
  data$freq <- as.numeric(data$freq)
  #extract type of substitution from fileName and insert as new variable
  data <- data %>% mutate(sub = str_extract(fileName,  "(CtoT|GtoA)"))
  #edit fileName to remove top directory
  data$fileName <- gsub('mapDamageData//', '', data$fileName)
  data$fileName <- gsub('_split.bam/5pCtoT_freq.txt', '', data$fileName)
  data$fileName <- gsub('_split.bam/3pGtoA_freq.txt', '', data$fileName)
  return(data)
}

#Apply function to subFreq.Files
subFreqData <- loadMisincorpData(subFreq.Files)

#Extract genome ID from the fileName and place this into data frame as separate variable
subFreqData <- extractGenomeName(subFreqData)
subFreqData$genome <- gsub('A.actinomyces', 'A.actinomycetemcomitans', subFreqData$genome)

#Edit fileNames to be consistent with previous
subFreqData$fileName <- editSimDataFileNames(subFreqData)

#Generate a function for plotting subFreq
subFreqPlot <- function(x, ...){x %>% 
    ggplot(aes(x=genome, y=freq)) + 
    geom_point(size=4) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    facet_wrap(~sub) + 
    scale_y_continuous(limits = c(0,0.5)) +
    labs(x="", y="Misincorporation frequency", title="Substitution frequency at the first and final position of the read", 
         subtitle=paste(x$fileName)) + 
    theme(legend.position = "top", legend.title = NULL)}

#Plot subFreq at position 1 for all genomes on the same axes, using lapply to split by fileName
subFreqData %>% filter(grepl("30bp", fileName)) %>% 
  split(f = .$fileName) %>% lapply(function(x){x %>% 
    filter(pos == "1") %>% subFreqPlot()})

#####Add information about GC content, count and coverage to see if any of these can explain variaiton in deamination freq

##Add information about counts for each file/genome and then cluster (as did with MAPQ and lengths)
subFreqData <- Q25totalSplitCount %>% select(fileName, genome, count) %>% 
  right_join(subFreqData, by = c("fileName", "genome"))
#Collate readCounts into bins
subFreqData$countRange <- cut(subFreqData$count, breaks = c(0,500,1000,5000,10000,50000,100000,500000), 
                            labels = c("0-500", "501-1000", "1001-5000", "5001-10000", "10001-50000", 
                                       "50001-100000", "100001-500000"), right = FALSE)
#re-order columns
subFreqData <- subFreqData[,c(1:3,7,4:6)]


##Replot subFreqData using diff colours for the countRange

#Generate a function for plotting subFreq
subFreqPlot <- function(x, ...){x %>% 
    ggplot(aes(x=genome, y=freq, colour=countRange)) + 
    geom_point(size=4) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.position = "bottom", legend.title = element_blank()) + 
    facet_wrap(~sub) + 
    scale_y_continuous(limits = c(0,0.5)) + 
    scale_colour_manual(values = palette15) + 
    labs(x="", y="Misincorporation frequency", title="Substitution frequency at the first and final position of the read", 
         subtitle=paste(x$fileName))}

#Plot subFreq at position 1 for all genomes on the same axes, using lapply to split by fileName
subFreqData %>% filter(grepl("30bp", fileName)) %>% 
  split(f = .$fileName) %>% lapply(function(x){x %>% 
      filter(pos == "1") %>% 
      subFreqPlot()})

##Add info on GC content. Like before group similar and plot

#Import GC content information
#GC <- read_csv("simFiles_genomeData.csv", col_names = FALSE, col_types = "------c-n-", skip = 1) %>% 
#  set_colnames(c("genome", "GC"))

GC <- read_delim("expandedGenomeList.txt", delim = "\t", col_names = c("genome", "GC"), 
                 col_types = "c-----n-")
GC$genome <- gsub('A.actinomyces', 'A.actinomycetemcomitans', GC$genome)

#bind this info to subFreqData
subFreqData <- left_join(subFreqData, GC, by = "genome")
#Collate GCcontent into bins
subFreqData$GCRange <- cut(subFreqData$GC, breaks = c(20,40,60,80), 
                              labels = c("20-40", "41-60", "61-80"), right = FALSE)
#re-order columns
subFreqData <- subFreqData[,c(1:4,8:9,5:7)]

#Plot subFreq at position 1 for all genomes on the same axes, using lapply to split by fileName
subFreqData %>% filter(grepl("30bp", fileName)) %>% 
  split(f = .$fileName) %>% lapply(function(x){x %>% 
      filter(pos == "1") %>% 
      ggplot(aes(x=genome, y=freq, colour=GCRange)) + 
      geom_point(size=4) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      theme(legend.position = "bottom", legend.title = element_blank()) + 
      facet_wrap(~sub) + 
      scale_y_continuous(limits = c(0,0.5)) + 
      scale_colour_manual(values = palette15) + 
      labs(x="", y="Misincorporation frequency", title="Substitution frequency at the first and final position of the read", 
           subtitle=paste(x$fileName))})

########### Coverage ??? - First need to determine coverage using bash script, import and assign ####
