#########################################################################
# get_DE_genes_july2017.R
#
# D. Lemay 
# R script to obtain differentially regulated genes in Hwang Lab samples
# This version to test time*arm interaction term
##########################################################################

setwd("/Users/danielle.lemay/work/hwang/data")

source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

# Get the raw count data
countsTableAll <- read.table(file="fixed_combinedCounts.txt", header=TRUE, sep="\t", row.names=1)

# Get the phenoData table 
# Each sample (Label) is part of a blinded study arm (Arm): (R, S, or T)
# Each sample has a time point (Time): Fast, 3hr, 6hr
# Each sample has a subject number (Subject): S17, S21, S27, S2, S9
pD.all <- read.table(file="phenoData.txt", header=TRUE, sep="\t", row.names="Label")
pD <- subset(pD.all, select=-Day)

# Keep only rows (genes) with nonzero counts  (25346 reduced to 21028 elements)
countsTableExp <- countsTableAll[ rowSums(countsTableAll) > 1, ]

# make the reference level for Time be the "Fasting" condition
# make the reference level for Arm be "R" (although we don't know what is R)
pD$Time <- relevel(pD$Time, ref="Fast")
pD$Arm <- relevel(pD$Arm, ref="R")

#############################################################################################################
# Statistical analysis strategy
# (1) First independently check 
# (1a) the effect of Subject, 
# (1b) the effect of Time,  
# (1c) the effect of Arm, 
# (1d) the effect of Day,
#
# Then
# (1d) Time:Arm 
# (1e) Time:Subject 
# (1f) Subject:Arm
# (1g) Time:Day
# (1h) Subject:Day
#
# (2) Next, the exploratory analyses (see cluster_analyses_Jan2017.R) showed that Subject had biggest effect
# even though we are most interested in effect of Time and of Arm
# So let's look at 
# (2a) The effect of Time, including the Subject:Time interaction in the model (the subject-specific effects)
# (2b) The effect of Arm, including the Subject:Arm interaction in the model (the subject-specific effects)
# (2c) The effect of Arm, including Time:Arm interaction in the model (using Subjects as replicates)
#############################################################################################################

##########################
# (1a) Effect of Subject (main effects recommended to be tested by Wald test, these direct contrasts)
#########################
ddsSubj <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD)], colData = pD, design= ~Subject)
ddsSubj <- DESeq(ddsSubj)
resultsNames(ddsSubj)
# [1] "Intercept"  "SubjectS17" "SubjectS2"  "SubjectS21" "SubjectS27" "SubjectS9"
# all pairwise contrasts
resSubj.1 <- results(ddsSubj, contrast=list("SubjectS17","SubjectS2"))
resSubj.2 <- results(ddsSubj, contrast=list("SubjectS17","SubjectS21"))
resSubj.3 <- results(ddsSubj, contrast=list("SubjectS17","SubjectS27"))
resSubj.4 <- results(ddsSubj, contrast=list("SubjectS17", "SubjectS9"))
resSubj.5 <- results(ddsSubj, contrast=list("SubjectS2", "SubjectS21"))
resSubj.6 <- results(ddsSubj, contrast=list("SubjectS2", "SubjectS27"))
resSubj.7 <- results(ddsSubj, contrast=list("SubjectS2", "SubjectS9"))
resSubj.8 <- results(ddsSubj, contrast=list("SubjectS21", "SubjectS27"))
resSubj.9 <- results(ddsSubj, contrast=list("SubjectS21", "SubjectS9"))
resSubj.10 <- results(ddsSubj, contrast=list("SubjectS27", "SubjectS9"))

resSig.Subj.1 <- data.frame(subset(resSubj.1, padj < 0.05)) # 8017 rows
resSig.Subj.2 <- data.frame(subset(resSubj.2, padj < 0.05)) # 6500 rows
resSig.Subj.3 <- data.frame(subset(resSubj.3, padj < 0.05)) # 6521 rows
resSig.Subj.4 <- data.frame(subset(resSubj.4, padj < 0.05)) # 5363 rows
resSig.Subj.5 <- data.frame(subset(resSubj.5, padj < 0.05)) # 6734 rows
resSig.Subj.6 <- data.frame(subset(resSubj.6, padj < 0.05)) # 4730 rows
resSig.Subj.7 <- data.frame(subset(resSubj.7, padj < 0.05)) # 8393 rows
resSig.Subj.8 <- data.frame(subset(resSubj.8, padj < 0.05)) # 4073 rows
resSig.Subj.9 <- data.frame(subset(resSubj.9, padj < 0.05)) # 6268 rows
resSig.Subj.10 <- data.frame(subset(resSubj.10, padj < 0.05)) # 6639 rows

# merge a bunch of data frames, keep all rows
BulkMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}

# to get every DE gene, we'd need to take every contrast and then merge them all
all.Subjects <- Reduce(BulkMerge, list(resSig.Subj.1,resSig.Subj.2,resSig.Subj.3,resSig.Subj.4,resSig.Subj.5,resSig.Subj.6,resSig.Subj.7,resSig.Subj.8,resSig.Subj.9,resSig.Subj.10))

#the Bulk Merge results in complaints about column name duplications
#I did verify (by awk and uniquify) that gene names are not duplicated

#13,856 genes are differentially regulated between subjects

write.table(all.Subjects, file="../results/DESeq_allgenes_resSubjects.txt", quote=FALSE, sep="\t")

##########################
# (1b) Effect of Time
#########################
ddsTime <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD)], colData = pD, design= ~Time)
ddsTime <- DESeq(ddsTime)
resultsNames(ddsTime)
##>[1] "Intercept" "TimeFast"  "Time3hr"   "Time6hr" 
resTime.Fastvs3hr <- results(ddsTime, contrast=list("TimeFast", "Time3hr"))
resTime.Fastvs6hr <- results(ddsTime, contrast=list("TimeFast", "Time6hr"))
resTime.3hrvs6hr <- results(ddsTime, contrast=list("Time3hr", "Time6hr"))

resSig.Time.Fv3 <- data.frame(subset(resTime.Fastvs3hr, padj < 0.05)) # 85 genes
resSig.Time.Fv6 <- data.frame(subset(resTime.Fastvs6hr, padj < 0.05)) # 26 genes
resSig.Time.3v6 <- data.frame(subset(resTime.3hrvs6hr, padj < 0.05)) # 3 genes


all.Time <- Reduce(BulkMerge, list(resSig.Time.Fv3, resSig.Time.Fv6, resSig.Time.3v6))
#95 genes

# main effect of time, 95 genes
write.table(all.Time, file="../results/DESeq_allgenes_resTime.txt", quote=FALSE, sep="\t")

##########################
# (1c) Effect of Arm 
#########################
ddsArm <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD)], colData = pD, design= ~Arm)
ddsArm <- DESeq(ddsArm)
resultsNames(ddsArm)
#[1] "Intercept" "ArmR"     "ArmS"     "ArmT"
resArm.RvsS <- results(ddsArm, contrast=list("ArmR", "ArmS"))
resArm.RvsT <- results(ddsArm, contrast=list("ArmR", "ArmT"))
resArm.SvsT <- results(ddsArm, contrast=list("ArmS", "ArmT"))
resSig.Arm.1v2 <- data.frame(subset(resArm.RvsS, padj < 0.05)) # 0 genes
resSig.Arm.1v3 <- data.frame(subset(resArm.RvsT, padj < 0.05)) # 0 genes
resSig.Arm.2v3 <- data.frame(subset(resArm.SvsT, padj < 0.05)) # 0 genes
all.Arm <- Reduce(BulkMerge, list(resSig.Arm.1v2, resSig.Arm.1v3, resSig.Arm.2v3)) #0 genes


################################################################################
# (1d) the main effect of Day
##################################################################################

pD.all$Time <- relevel(pD.all$Time, ref="Fast")
pD.all$Arm <- relevel(pD.all$Arm, ref="R")
pD.all$Day <- relevel(pD.all$Day, ref="D1")


ddsDay <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD.all)], colData = pD.all, design= ~Day)
ddsDay <- DESeq(ddsDay)
resultsNames(ddsDay)
# [1] "Intercept" "DayD1"     "DayD2"     "DayD3"
resDay.D1vsD2 <- results(ddsDay, contrast=list("DayD1", "DayD2"))
resDay.D1vsD3 <- results(ddsDay, contrast=list("DayD1", "DayD3"))
resDay.D2vsD3 <- results(ddsDay, contrast=list("DayD2", "DayD3"))

resSig.Day.D1vD2 <- data.frame(subset(resDay.D1vsD2, padj < 0.05)) # 0 genes
resSig.Day.D1vD3 <- data.frame(subset(resDay.D1vsD3, padj < 0.05)) # 2 genes
resSig.Day.D2vD3 <- data.frame(subset(resDay.D2vsD3, padj < 0.05)) # 0 genes

# merge a bunch of data frames, keep all rows
BulkMerge       <- function(x, y){
  df            <- merge(x, y, by= "row.names", all.x= T, all.y= T)
  rownames(df)  <- df$Row.names
  df$Row.names  <- NULL
  return(df)
}
all.Day <- Reduce(BulkMerge, list(resSig.Day.D1vD2, resSig.Day.D1vD3, resSig.Day.D2vD3))
# 2 genes

# main effect of time, 2 genes
write.table(all.Day, file="../results/DESeq_allgenes_resDay.txt", quote=FALSE, sep="\t")


##################################################################################
# (1d) Time:Arm 
# (1e) Time:Subject 
# (1f) Subject:Arm
# (1g) Time:Day
# (1h) Subject:Day
###################################################################################

################
# (1d) Time:Arm
# effect of Time:Arm interaction
################
ddsAllArm <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD.all)], colData = pD.all, design= ~Arm + Time + Subject + Time:Subject + Arm:Subject + Time:Arm)
ddsAllArm <- DESeq(ddsAllArm, test="LRT", reduced = ~ Arm + Time + Subject + Time:Subject + Arm:Subject)
#116 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT
#Note to self: this is likely do to low-count genes. We could also increase number of iterations (with maxit) if we weren't using the wrapper
# functions (which don't have maxit as a parameter)

# The p-value is based on the LRT (so the same for everything), but the contrasts give the right log2 fold change.
resAllArm.RvS <- results(ddsAllArm, contrast=c("Arm", "R", "S"))
resSig.AllArm.RvS <- data.frame(subset(resAllArm.RvS, padj < 0.05)) 
nrow(resSig.AllArm.RvS) # 0 genes

####################
# (1e) Time:Subject
####################
ddsAllTimeSubj <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD.all)], colData = pD.all, design= ~Arm + Time + Subject + Arm:Subject + Time:Arm + Time:Subject)
ddsAllTimeSubj <- DESeq(ddsAllTimeSubj, test="LRT", reduced = ~ Arm + Time + Subject + Arm:Subject + Time:Arm)
#111 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT

resultsNames(ddsAllTimeSubj)

# The p-value is based on the LRT (so the same for everything), but the contrasts give the right log2 fold change.
resAllTimeSubj.Fv3 <- results(ddsAllTimeSubj, contrast=c("Time", "Fast", "3hr"))
resSig.AllTimeSubj.Fv3 <- data.frame(subset(resAllTimeSubj.Fv3, padj < 0.05)) 
nrow(resSig.AllTimeSubj.Fv3) # 651 genes

write.table(resSig.AllTimeSubj.Fv3, file="../results/DESeq_TimeSubjInteraction.txt", quote=FALSE, sep="\t")

####################
# (1e) Subject:Arm
####################
ddsAllSubjArm <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD.all)], colData = pD.all, design= ~Arm + Time + Subject + Time:Subject + Time:Arm + Arm:Subject)
ddsAllSubjArm <- DESeq(ddsAllSubjArm, test="LRT", reduced = ~ Arm + Time + Subject + Time:Subject + Time:Arm)
#109 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT

# The p-value is based on the LRT (so the same for everything), but the contrasts give the right log2 fold change.
resAllSubjArm.RvS <- results(ddsAllSubjArm, contrast=c("Arm", "R", "S"))
resSig.AllSubjArm.RvS <- data.frame(subset(resAllSubjArm.RvS, padj < 0.05)) 
nrow(resSig.AllSubjArm.RvS) # 2175 genes

write.table(resSig.AllSubjArm.RvS, file="../results/DESeq_SubjArmInteraction.txt", quote=FALSE, sep="\t")

#####################
# (1g) Time:Day
#######################
ddsAllTimeDay <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD.all)], colData = pD.all, design= ~Day + Time + Subject + Time:Subject + Day:Subject + Time:Day)
ddsAllTimeDay <- DESeq(ddsAllTimeDay, test="LRT", reduced = ~ Day + Time + Subject + Time:Subject + Day:Subject)
#111 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT

# The p-value is based on the LRT (so the same for everything), but the contrasts give the right log2 fold change.
resAllTimeDay.D1vD2 <- results(ddsAllTimeDay, contrast=c("Day", "D1", "D2"))
resSig.AllTimeDay.D1vD2 <- data.frame(subset(resAllTimeDay.D1vD2, padj < 0.05)) 
nrow(resSig.AllTimeDay.D1vD2) # 0 genes

#####################
# (1h) Subject:Day
#####################
# effect of Time:Subject interaction
ddsAllSubjDay <- DESeqDataSetFromMatrix(countData = countsTableExp[,rownames(pD.all)], colData = pD.all, design= ~Day + Time + Subject + Time:Subject + Time:Day + Day:Subject)
ddsAllSubjDay <- DESeq(ddsAllSubjDay, test="LRT", reduced = ~ Day + Time + Subject + Time:Subject + Time:Day)
#109 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT

# The p-value is based on the LRT (so the same for everything), but the contrasts give the right log2 fold change.
resAllSubjDay.D1vD2 <- results(ddsAllSubjDay, contrast=c("Day", "D1", "D2"))
resSig.AllSubjDay.D1vD2 <- data.frame(subset(resAllSubjDay.D1vD2, padj < 0.05)) 
nrow(resSig.AllSubjDay.D1vD2) # 2276 genes

write.table(resSig.AllSubjDay.D1vD2, file="../results/DESeq_SubjDayInteraction.txt", quote=FALSE, sep="\t")
