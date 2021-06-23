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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("affyio")


library(affy)
BiocManager::install("affycomp")
BiocManager::install("affydata")
library(tools)


############################
#####   Load CEL files #####
############################

data.raw <- ReadAffy(filenames="GSM352329.CEL", "GSM352330.CEL", "GSM352331.CEL")
                     
eset <- mas5(data.raw)
x <-  exprs(eset)
write.table(data.frame(x,check.names=FALSE),file="testToCSV.csv",sep=",",col.names=NA,quote=FALSE)

############################
##### PRE-PROCESS DATA #####
############################

# load in files
counts <- read.cel("./dataF_unmerged/b/asbes_high_8h/GSM352329.CEL", sep="\t")

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



split_counts(x = "GSM")

# session information
sessionInfo()
