#clean workspace
rm(list=ls())
 
#set the working directory
setwd("C:/Users/LARS/Documents/bigcat/Rscripts/")

#for now generate some random data for illustration purposes
dat <- matrix(rnorm(240000),ncol=24)

#create experimental group factor variable
desc <- cbind(group=c(rep("control",12),rep("treated",12)))
desc <- cbind(desc,time=c(rep("t1",6),rep("t2",6),rep("t1",6),rep("t2",6)))
desc <- as.data.frame(cbind(desc,individual=paste("s",rep(c(1:6),4),sep="")))
rownames(desc) <- paste(desc$group,desc$time,desc$individual,sep="_")

#add sample and gene names to dat
colnames(dat) <- rownames(desc)
rownames(dat) <- paste("gene",1:dim(dat)[1],sep="")

#load needed libraries
library(limma)
library(bioDist)
library(gplots)

#include some functions from ArrayAnalysis.org
source ("functions_ArrayAnalysis_v2.R")

#create plots
factors <- c("group","time","individual")
createQCPlots(dat, factors, Table=desc, normMeth="", postfix="")

###########################
##statistical modelling 1##
###########################
	   
#model design: take group_ as a fixed effect, no pairing (so assuming all samples are from different individuals)
design <- model.matrix(~group,data=desc)
colnames(design) <- gsub("[()]","",colnames(design))
#colnames(design) <- gsub(":","_x_",colnames(design))

#fit model
fit <- lmFit(dat,design)
fit <- eBayes(fit)

#extract resulting statistics based on the model, and save those in a table; also save some graphical representations
#topTable(fit, adjust.method="BH", coef=2, number=dim(dat)[1], resort.by="P")
#files <- saveStatOutput(design,fit,postfix="",annotation=NULL)

#create summary table of the model results
#createPvalTab(files,postfix="",html=TRUE)

#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrix <- makeContrasts(
  treated_vs_control = grouptreated,
  levels = colnames(design)
)
#compute the contrast fits
contrast.fit <- contrasts.fit(fit, cont.matrix)
contrast.fit <- eBayes(contrast.fit)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
files.c <- saveStatOutput(cont.matrix,contrast.fit,postfix="",annotation=NULL)

#create summary table of the contrast results
createPvalTab(files.c,postfix="c",html=TRUE)

###########################
##statistical modelling 2##
###########################
	   
#model design: take group_ as a fixed effect, individual_ as the pairing variable
design <- model.matrix(~group,data=desc)
colnames(design) <- gsub("[()]","",colnames(design))
#colnames(design) <- gsub(":","_x_",colnames(design))

#estimate the extra correlation between measurements made on the same subject
corfit <- duplicateCorrelation(dat,design,block=desc$individual)
corfit$consensus
#will be very small value with random data...so no use to do this, but for real data it should be bigger

#fit model
fit <- lmFit(dat,design,block=desc$individual_,correlation=corfit$consensus)
fit <- eBayes(fit)

#extract resulting statistics based on the model, and save those in a table; also save some graphical representations
#topTable(fit, adjust.method="BH", coef=2, number=dim(dat)[1], resort.by="P")
#files <- saveStatOutput(design,fit,postfix="",annotation=NULL)

#create summary table of the model results
#createPvalTab(files,postfix="",html=TRUE)

#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrix <- makeContrasts(
  treated_vs_control = grouptreated,
  levels = colnames(design)
)
#compute the contrast fits
contrast.fit <- contrasts.fit(fit, cont.matrix)
contrast.fit <- eBayes(contrast.fit)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
files.c <- saveStatOutput(cont.matrix,contrast.fit,postfix="",annotation=NULL)

#create summary table of the contrast results
createPvalTab(files.c,postfix="c",html=TRUE)

###########################
##statistical modelling 3##
###########################
	   
#model design: take group_ and time_ as interacting fixed effects, individual_ as the pairing variable
design2 <- model.matrix(~group*time,data=desc)
colnames(design2) <- gsub("[()]","",colnames(design2))
colnames(design2) <- gsub(":","_x_",colnames(design2))

#estimate the extra correlation between measurements made on the same subject
corfit2 <- duplicateCorrelation(dat,design2,block=desc$individual)
corfit2$consensus
#will be very small value with random data...so no use to do this, but for real data it should be bigger

#fit model
fit2 <- lmFit(dat,design2,block=desc$individual_,correlation=corfit2$consensus)
fit2 <- eBayes(fit2)

#extract resulting statistics based on the model, and save those in a table; also save some graphical representations
#topTable(fit2, adjust.method="BH", coef=2, number=dim(dat)[1], resort.by="P")
#files2 <- saveStatOutput(design2,fit2,postfix="",annotation=NULL)

#create summary table of the model results
#createPvalTab(files2,postfix="2",html=TRUE)

#build the (contrast) matrix to compute some group differences of interest based on the model parameters
cont.matrix2 <- makeContrasts(
  treated_vs_control_t1 = grouptreated,
  treated_vs_control_t2 = grouptreated + grouptreated_x_timet2,
  treated_vs_control = grouptreated + 0.5*grouptreated_x_timet2,
  t2_vs_t1_control = timet2,
  t2_vs_t1_treated = timet2 + grouptreated_x_timet2,
  t2_vs_t1 = timet2 + 0.5*grouptreated_x_timet2,
  t2_vs_t1_treated_cor_for_control = grouptreated_x_timet2,
  levels = colnames(design2)
)
#compute the contrast fits
contrast.fit2 <- contrasts.fit(fit2, cont.matrix2)
contrast.fit2 <- eBayes(contrast.fit2)

#extract resulting contrasts based on the model, and save those in a table; also save some graphical representations
files.c2 <- saveStatOutput(cont.matrix2,contrast.fit2,postfix="2",annotation=NULL)

#create summary table of the contrast results
createPvalTab(files.c2,postfix="c2",html=TRUE)
