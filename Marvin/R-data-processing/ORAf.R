### EU TOxRisk Case study - RNAseq 
### Laurent Winckers, Marvin Martens - Maastricht University
### Department of Bioinformatics - BiGCaT

### Session info:
### R version 4.0.2 (2020-06-22)
### Platform: x86_64-w64-mingw32/x64 (64-bit)
### Running under: Windows 10 x64 (build 18362)
### Packages: 

##############################
##### SET UP ENVIRONMENT #####
##############################

# clean work space
rm(list=ls())

# set working directory
getwd()
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# install packages 
#used R version 4.0.5
BiocManager::install("edgeR") 
library(edgeR)
library(tidyr)
library(xlsx)
library(data.table)
library(clusterProfiler)
library(dplyr)
library(ggpubr)


# functions
### filLow
filLow <- function(data){
  
  
  
  myCPM <- cpm(data)
  head(myCPM)
  thresh <- myCPM > 0.5
  table(rowSums(thresh))
  keep <- rowSums(thresh) >= 1
  data <- data[keep,]
  data <- drop_na(data)
}

### exprs
exprs <- function(data,x,y,group){
  y <- DGEList(counts = data[,x:y], group = group)
  y <- calcNormFactors(y)
  design <- model.matrix(~1 + group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- as.data.frame(topTags(qlf,n=Inf))
}

################################
##### LOAD & RUN FUNCTIONS #####
################################

# get filenames
ls <- list.files(path = "./data2")
ls <- gsub(".txt", "", ls)



# loop over files and run functions
for (fileName in ls){

  # load datafile
  data <- read.table(paste0("./data2/", fileName, ".txt"), header = T, sep = "\t")


  # create rownames from probeID column -- column "X"
  row.names(data) <- data$X
  data <- data[-1]

#######################
##### FC ANALYSIS #####
#######################

  # filter lowly expressed genes
  data <- filLow(data = data)

  # run EdgeR 
  dataExpr <- exprs(data = data, x = 1, y =45, 
              group = factor(c(rep("control", 5), rep("case", 40)), 
                             levels = c("control", "case")))
} 


#session information
sessionInfo()





