#############################################################################
# partition_variances_feb2019.R
#
# D. Lemay 
# R script to determine drivers of variation in the whole blood transcriptomes
##############################################################################

setwd("/Users/danielle.lemay/work/hwang/data")

source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("variancePartition", version = "3.8")

# Get the count data
countDataAll <- read.table(file="fixed_combinedCounts.txt", header=TRUE, sep="\t", row.names=1)

# Get the phenoData table 
# Each sample is part of a study arm (Arm): (Placebo_SunflowerOil, Placebo_DHA, Blueberry_SunflowerOil)
# Each sample has a time point (Time): Fast, 3hr, 6hr
# Each sample has a subject number (Subject): S17, S21, S27, S2, S9
# Each sample was collected on a different study day: D1, D2, D3 (subjects were randomized to treatment order)
pD.all <- read.table(file="phenoData_unblinded.txt", header=TRUE, sep="\t", row.names="Label")

# make the reference level for Time be "Fast"
# make the reference level for Arm be "Placebo_SunflowerOil"
# make the reference level for Day be "D1"
pD.all$Time <- relevel(pD.all$Time, ref="Fast")
pD.all$Arm <- relevel(pD.all$Arm, ref="Placebo_SunflowerOil")
pD.all$Day <- relevel(pD.all$Day, ref="D1")

####################################################
# Following tutorial for variancePartition
####################################################
library(DESeq2)
### create DESeq2 object from gene-level counts and metadata
dds <- DESeqDataSetFromMatrix(countData = countDataAll[,rownames(pD.all)], colData = pD.all, design= ~1)

# estimate library size correction scaling factors
dds <- estimateSizeFactors(dds)

# identify genes that pass expression cut-off
isexpr <- rowSums(fpm(dds)>1) >= 0.5 * ncol(dds)

# compute log2 Fragments Per Million
quantLog <- log2( fpm( dds )[isexpr,] + 1)

# define formula (using all fixed effects since sample size is small and no variable is continuous)
form <- ~ Time + Day + Arm + Subject 

# load library
source("https://bioconductor.org/biocLite.R")
biocLite("variancePartition")
library('variancePartition')
# Run variancePartition analysis
varPart <- fitExtractVarPartModel( quantLog, form, pD.all)

# sort variables
vp <- sortCols( varPart )

# violin plot of contribution of each variable to total variance
pdf(file="../graphs/partition_variance.pdf")
plotVarPart( vp , col=c("lightblue", "yellow", "green", "orange", "gray"))
dev.off()

summary(vp)
#Subject               Time                Day           
#Min.   :0.0003819   Min.   :0.0000051   Min.   :0.0000036  
#1st Qu.:0.2622109   1st Qu.:0.0141503   1st Qu.:0.0090449  
#Median :0.4399153   Median :0.0354301   Median :0.0245929  
#Mean   :0.4501234   Mean   :0.0529326   Mean   :0.0377759  
#3rd Qu.:0.6285169   3rd Qu.:0.0739707   3rd Qu.:0.0527599  
#Max.   :0.9971730   Max.   :0.7466577   Max.   :0.3909434  
#Arm              Residuals       
#Min.   :1.400e-07   Min.   :0.002355  
#1st Qu.:8.315e-03   1st Qu.:0.264675  
#Median :2.016e-02   Median :0.427327  
#Mean   :3.029e-02   Mean   :0.428878  
#3rd Qu.:4.187e-02   3rd Qu.:0.589499  
#Max.   :2.919e-01   Max.   :0.954807 

# repeat without study day
form2 <- ~ Time + Arm + Subject 
varPart2 <- fitExtractVarPartModel( quantLog, form2, pD.all)
vp2 <- sortCols( varPart2 )
pdf(file="../graphs/partition_variance2.pdf")
plotVarPart( vp2, col=c("lightblue", "yellow", "green", "gray"))
dev.off()

summary(vp2)
#Subject               Time                Arm              Residuals       
#Min.   :0.0003819   Min.   :0.0000051   Min.   :1.400e-07   Min.   :0.002664  
#1st Qu.:0.2622109   1st Qu.:0.0141503   1st Qu.:8.315e-03   1st Qu.:0.296134  
#Median :0.4399153   Median :0.0354301   Median :2.016e-02   Median :0.472477  
#Mean   :0.4501234   Mean   :0.0529326   Mean   :3.029e-02   Mean   :0.466654  
#3rd Qu.:0.6285169   3rd Qu.:0.0739707   3rd Qu.:4.187e-02   3rd Qu.:0.635942  
#Max.   :0.9971730   Max.   :0.7466577   Max.   :2.919e-01   Max.   :0.972279 
