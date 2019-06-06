# whole_blood_RNASeq
R workflow to understand contributors to variance in human whole blood transcriptomes


## Please cite:

Lemay DG, et al. (2019) “Temporal changes in postprandial blood transcriptomes reveal subject-specific pattern of expression of innate immunity genes after a high fat meal” 

## Raw Data:

The raw data set for this project is found in NCBI’s Gene Expression Omnibus database (https://www.ncbi.nlm.nih.gov/geo/) with study accession number GSE127530 and series record GSE127530.

## Processed Data:

fixed_combinedCounts.txt : raw counts file with non-gene features removed (sed '/^N_unmapped/d' | sed '/^N_multimapping/d' | sed '/N_noFeature/d' | sed '/N_ambiguous/d')

phenoData.txt: meta-data file (blinded to treatment arm)

phenoData_unblinded.txt : meta-data file (unblinded to treatment arm)


## R code:

get_DEgenes_july2017_clean.R: R script to obtain differentially regulated genes in whole blood transcriptomes, accounting for multiple drivers of variation

partition_variances_feb2019.R : R script to determine drivers of variation in the whole blood transcriptomes


### License information

This software is a work of the United States Department of Agriculture, Agricultural Research Service and is released under a Creative Commons CC0 public domain attribution.
