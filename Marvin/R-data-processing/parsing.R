### Internship
### Laurent Winckers, Marvin Martens - Maastricht University
### Department of Bioinformatics - BiGCaT

##############################
##### SET UP ENVIRONMENT #####
##############################

# clean work space
rm(list=ls())

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# load libraries
## no libraries 

############################
##### PRE-PROCESS DATA #####
############################

# load in files
counts <- read.csv("./GSE147507_RawReadCounts_Human.tsv", sep="\t")

# check colnames 
unique(colnames(counts))

# remove everything from box 3 from Excel sheet -- IMI, ACE, NIC
counts_fil <- counts[,!grepl("A549|Calu3|HealthyLungBiopsy|COVID|Series9", colnames(counts))]
unique(colnames(counts_fil))


########################
##### SPLIT COUNTS #####
########################

# function
split_counts <- function(x) { 
  y <- c("X", paste0("Series1_NHBE_Mock"), x)
  z <- counts_fil[,grepl(paste(y, collapse = "|"), colnames(counts_fil))]
  write.table(z, file = paste(x, ".txt", sep=""), row.names = F, quote = F, sep = "\t")
  }

# use function to split
split_counts(x = "Series1_NHBE_SARS.CoV.2")

# session information
sessionInfo()
