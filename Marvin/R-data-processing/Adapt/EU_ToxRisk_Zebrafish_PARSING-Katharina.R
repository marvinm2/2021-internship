### EU TOxRisk Case study 
### Laurent Winckers, Marvin Martens - Maastricht University
### Department of Bioinformatics - BiGCaT

### Session info:
### R version 4.0.2 (2020-06-22)
### Platform: x86_64-w64-mingw32/x64 (64-bit)
### Running under: Windows 10 x64 (build 18362)
### Packages: -

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
counts <- read.csv("./Counts-per-gene-per-sample_raw_EUT080.csv")

# check colnames 
unique(colnames(counts))

# fix wrong colname
colnames(counts)[which(names(counts) == "S_2.Mhex_EC10_48h_run4_C09")] <- "S_2.Mhex_1.2EX10_48h_run4"

# remove everything from box 3 from Excel sheet -- IMI, ACE, NIC
counts_fil <- counts[,!grepl("IMI|ACE|NIC", colnames(counts))]

# remove S and S_2 from column names
names(counts_fil) <- gsub(pattern = "^S_|S_2.", replacement = "", x = names(counts_fil))

########################
##### SPLIT COUNTS #####
########################

# function
split_counts <- function(x, timepoint) { 
  y <- c("X", paste0("DMSO_",timepoint), paste0("NK_",timepoint), x)
  z <- counts_fil[,grepl(paste(y, collapse = "|"), colnames(counts_fil))]
  write.table(z, file = paste("/katharina/data/", x, ".txt", sep=""), row.names = F, quote = F, sep = "\t")
  }

# use function to split
## Ebut
split_counts(x = "Ebut_1.2EC10_24h", timepoint = "24h")
split_counts(x = "Ebut_EC10_24h", timepoint = "24h")
split_counts(x = "Ebut_EC20_24h", timepoint = "24h")
split_counts(x = "Ebut_1.2EC10_48h", timepoint = "48h")
split_counts(x = "Ebut_EC10_48h", timepoint = "48h")
split_counts(x = "Ebut_EC20_48h", timepoint = "48h")

## VPA
split_counts(x = "VPA_1.2EC10_24h", timepoint = "24h")
split_counts(x = "VPA_EC10_24h", timepoint = "24h")
split_counts(x = "VPA_EC20_24h", timepoint = "24h")
split_counts(x = "VPA_1.2EC10_48h", timepoint = "48h")
split_counts(x = "VPA_EC10_48h", timepoint = "48h")
split_counts(x = "VPA_EC20_48h", timepoint = "48h")

## Mhex
split_counts(x = "Mhex_1.2EC10_24h", timepoint = "24h")
split_counts(x = "Mhex_EC10_24h", timepoint = "24h")
split_counts(x = "Mhex_EC20_24h", timepoint = "24h")
split_counts(x = "Mhex_1.2EC10_48h", timepoint = "48h")
split_counts(x = "Mhex_EC10_48h", timepoint = "48h")
split_counts(x = "Mhex_EC20_48h", timepoint = "48h")


# session information
sessionInfo()
