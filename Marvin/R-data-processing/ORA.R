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
BiocManager::install("edgeR") 
library(edgeR)
library(tidyr)
library(xlsx)
library(data.table)
library(clusterProfiler)
library(dplyr)
library(ggpubr)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# load ID file
IDs <- read.xlsx("./190820 Zebrafish S1500+ Surrogate 1.0 Manifest.xlsx", sheetIndex = 1)

# load WP gmt file and split/clean
### load genesets
wp2gene <- clusterProfiler::read.gmt("wikipathways-20201010-gmt-Danio_rerio.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

# functions
### filLow
filLow <- function(data){
  
  #library(edgeR)
  #library(tidyr)
  
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

### ORA 
ora <- function(wpid2gene,wpid2name, data, fccutoff, file) {
  
  #library(dplyr)
  #library(ggpubr)
  
  buffer <- data %>% dplyr::select(starts_with("entrez") | starts_with("logFC") | starts_with("PValue"))
  
  # ORA
  up.genes <- buffer[buffer[,2] > fccutoff & buffer[,3] < 0.05, 1] 
  dn.genes <- buffer[buffer[,2] < -fccutoff & buffer[,3] < 0.05, 1]
  combined.genes <- c(up.genes,dn.genes)
  
  ewp.up <- clusterProfiler::enricher(
    up.genes,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 1, maxGSSize = 3000,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
  
  write.table(ewp.up, file = paste0("results/up/upORA_",file ,".txt"), sep="\t", quote = F, row.names = F)
  
  ewp.down <- clusterProfiler::enricher(
    dn.genes,
    pvalueCutoff = 1, 
    qvalueCutoff = 1,
    minGSSize = 1, maxGSSize = 3000,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
  
  write.table(ewp.down, file = paste0("results/down/downORA_",file ,".txt"), sep="\t", quote = F, row.names = F)
  
  
  ewp.combined <- clusterProfiler::enricher(
    combined.genes,
    pvalueCutoff = 1, 
    qvalueCutoff = 1,
    minGSSize = 1, maxGSSize = 3000,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
  
  write.table(ewp.combined, file = paste0("results/combined/ORA_",file ,".txt"), sep="\t", quote = F, row.names = F)

}

################################
##### LOAD & RUN FUNCTIONS #####
################################

# get filenames
ls <- list.files(path = "./data")
ls <- gsub(".txt", "", ls)

# loop over files and run functions
for (fileName in ls){

  # load datafile
  data <- read.table(paste0("./data/", fileName, ".txt"), header = T, sep = "\t")
  data <- data %>% select(X, starts_with("DMSO"), everything(), -starts_with("NK"))

  # create rownames from probeID column -- column "X"
  row.names(data) <- data$X
  data <- data[-1]

#######################
##### FC ANALYSIS #####
#######################

  # filter lowly expressed genes
  data <- filLow(data = data)

  # run EdgeR 
  dataExpr <- exprs(data = data, x = 1, y = 10, 
              group = factor(c(rep("control", 5), rep("case", 5)), 
                             levels = c("control", "case")))

############################
##### FINISHING RESULT #####
############################

  # clean result file, rownames (gene IDs) as first column -- match with IDs from ID file and save them
  dataExpr <- setDT(dataExpr, keep.rownames = TRUE)[]
  colnames(dataExpr)[1] <- "IDs"
  dataExpr <- dataExpr %>% tidyr::separate(IDs, c("gene_symbol","ID"), "_")

  vp <- EnhancedVolcano(dataExpr, title = paste0(fileName), lab = dataExpr$ID, 
                      x = 'logFC',y = 'PValue', xlim = c(min(dataExpr$logFC), max(dataExpr$logFC)), FCcutoff = 0.26, pCutoff = 0.05, 
                      col = c("grey30", "orange", "royalblue", "darkorange4"), selectLab = '', ylim = c(0, max(-log10(dataExpr$PValue))+0.5)) 

  svg(file=paste0("./results/volcanoPlots/", fileName, ".svg"), width = 17, height = 18)
  plot(vp)
  dev.off()
  
  dataExpr <- merge.data.frame(dataExpr, IDs, by.x = "ID", by.y = "pid")
  dataExpr <- dataExpr[-2]

  write.table(dataExpr, file = paste0("results/expr/expr_", fileName,".txt"), sep="\t", quote = F, row.names = F)

##################
##### PW ORA #####
##################

  # ORA 
  ora(wpid2gene, wpid2name, data = dataExpr, 0.26, file = paste0(fileName))

}

# session information
sessionInfo()



