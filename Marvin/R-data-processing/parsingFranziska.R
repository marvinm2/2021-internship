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
counts <- read.csv("./GDS1220_full.csv", sep="\t")

# check colnames 
unique(colnames(counts))

# remove everything from box 3 from Excel sheet -- IMI, ACE, NIC
counts_fil <- counts[,!grepl("GSM49609|GSM49610|GSM49611|GSM49612|GSM49614|GSM49615|GSM49616|GSM49617|GSM49613|Gene.ID", colnames(counts))]
unique(colnames(counts_fil))


########################
##### SPLIT COUNTS #####
########################

# function
split_counts <- function(x) { 
  y <- c("X", paste0("Control"), x)
  z <- counts_fil[,grepl(paste(y, collapse = "|"), colnames(counts_fil))]
  write.table(z, file = paste(x, ".txt", sep=""), row.names = F, quote = F, sep = "\t")
  }


# use function to split
for (item in colnames(counts_fil)) {
split_counts(x = item)
}
# session information
sessionInfo()
