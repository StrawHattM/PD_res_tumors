### Martin Gonzalez Fernandez, 2021.06.09


# Initial Dependencies ----------------------------------------------------

library(maftools, quietly = TRUE)
library(mclust, quietly = TRUE)
library(R.utils, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(berryFunctions, quietly = TRUE)
library(purrr, quietly = TRUE)
library(shiny, quietly = TRUE)
library(DT, quietly = TRUE)
library(pheatmap, quietly = TRUE)
library(BiocManager, quietly = TRUE)
library(BSgenome, quietly = TRUE)
library("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE)
library("BSgenome.Mmusculus.UCSC.mm9", quietly = TRUE)
library("BSgenome.Mmusculus.UCSC.mm8", quietly = TRUE)
library(NMF, quietly = TRUE)
library(pheatmap, quietly = TRUE)


# Functions needed for the script to work
source("dependent_functions.R") 




# Setting up the environment ----------------------------------------------

# Set up here your strings with names of MAF groups. 

DRT_mafs <- c("DRT_unfilt", "DRT_10", "DRT_25", "DRT_50", "DRT_75", "DRT_100", "DRT_all")
DST_mafs <- c("DST_unfilt", "DST_10", "DST_25", "DST_50", "DST_75", "DST_100", "DST_all")
all_mafs<- c("D.all_unfilt", "D.all_10", "D.all_100", "D.all_25", "D.all_50", "D.all_75", "D.all_all")

# Assignment and cleanup loops (particular for Genevia data) 
## assumes /MAF_files/mafname_combined.maf structure. Modify if needed.
## Cleanup consists of standardizing Variant Classification names and assignment
## of Protein_Change column according to HGVS.short.
## IMPORTANT: HUGO SYMBOLS ARE IN TITLE CASE FOR MOUSE GENES. solve pipe capHS() 

for(maf in DRT_mafs) {
  assign(x = maf, read.maf(maf = paste0(
    "MAF_files/", paste(maf, "combined.maf", sep = "_")
  )))
  
  tempname <-
    get(maf) %>% subsetMaf(tsb = levels(get(maf)@data$Tumor_Sample_Barcode),
                           mafObj = FALSE) %>%
    mutate(Variant_Classification = as.character(Variant_Classification)) %>%
    mutate(
      Variant_Classification = replace(
        x = Variant_Classification,
        list = Variant_Classification == "frameshift_variant_INS",
        values = "Frame_Shift_Ins"
      )
    ) %>%
    mutate(
      Variant_Classification = replace(
        x = Variant_Classification,
        list = Variant_Classification == "frameshift_variant_DEL",
        values = "Frame_Shift_Del"
      )
    ) %>%
    mutate(Variant_Classification = as.factor(Variant_Classification)) %>%
    read.maf()
  
  tempname <- add.protchange(tempname) %>%
    extract.maf() %>%
    read.maf(clinicalData = "DRT_clinical.tsv")
  
  assign(x = paste0(maf), tempname)
  
  rm(maf, tempname)
}

for(maf in DST_mafs){
  
  assign(x = maf, read.maf(maf = paste0(
    "MAF_files/", paste(maf, "combined.maf", sep = "_")
  )))
  
  tempname <-
    get(maf) %>% subsetMaf(tsb = levels(get(maf)@data$Tumor_Sample_Barcode),
                           mafObj = FALSE) %>%
    mutate(Variant_Classification = as.character(Variant_Classification)) %>%
    mutate(
      Variant_Classification = replace(
        x = Variant_Classification,
        list = Variant_Classification == "frameshift_variant_INS",
        values = "Frame_Shift_Ins"
      )
    ) %>%
    mutate(
      Variant_Classification = replace(
        x = Variant_Classification,
        list = Variant_Classification == "frameshift_variant_DEL",
        values = "Frame_Shift_Del"
      )
    ) %>%
    mutate(Variant_Classification = as.factor(Variant_Classification)) %>%
    read.maf()
  
  tempname <- add.protchange(tempname) %>%
    extract.maf() %>%
    read.maf(clinicalData = "DST_clinical.tsv")
  
  assign(x = paste0(maf), tempname)
  
  rm(maf, tempname)
}

## It's also useful to get names of tumors for each group
## it is however coded into some functions to automatically get them.

DRT1_tumornames <-
  get(DRT_mafs[1])@data$Tumor_Sample_Barcode %>% levels()

DST_tumornames <-
  get(DST_mafs[1])@data$Tumor_Sample_Barcode %>% levels()

# Initial exploration and graphs ------------------------------------------

## With default maftools ranking (top genes = present in most samples)
DRT_unfilt %>% plotmafSummary()
DST_unfilt %>% plotmafSummary()

DRT_unfilt %>% oncoplot()
DST_unfilt %>% oncoplot()

## Summary with basic parameters ordered by Variant_Classification
## Use unfiltered or low filtering data. adjust samples in loop 

for (maf in c("DRT_unfilt", "DRT_10", "DST_unfilt", "DST_10")) {
  
  cat("\n")
  cat(paste("Summary of", maf))
  cat("\n")
  get(maf) %>%
    extract.maf() %>%
    group_by(Variant_Classification) %>%
    summarise(
      Variants = n(),
      Genes = n_distinct(Hugo_Symbol),
      Mean_depth = mean(t_depth),
      Median_depth = median(t_depth),
      Mean_VAF = mean(VAF),
      Median_VAF = median(VAF)
    ) %>%
    print()
  
  rm(maf)
}

##istogram plots, useful to check initial distribution.
## Use mainly unfiltered data, but more can be checked.

## Sequencing Depth

DRT_unfilt %>%
  extract.maf() %>%
  ggplot(aes(t_depth, fill=..x..)) +
  geom_histogram(bins = 100) +
  scale_fill_distiller(palette = "YlOrRd")+
  ggtitle(paste("Sequencing depth in Unfiltered Docetaxel Resistant"))

DST_unfilt %>%
  extract.maf() %>%
  ggplot(aes(t_depth, fill=..x..)) +
  geom_histogram(bins = 100)+
  scale_fill_distiller(palette = "BuPu") +
  ggtitle(paste("Sequencing depth in Unfiltered Docetaxel Sensitive"))

## Variant Allele Fraction (VAF)

DRT_unfilt %>%
  extract.maf() %>%
  ggplot(aes(VAF, fill=..x..)) +
  geom_histogram(bins = 100)+
  scale_fill_distiller(palette = "YlOrBr") +
  ggtitle(paste("Variant Allele Fraction in Unfiltered Docetaxel Resistant"))

DST_unfilt %>%
  extract.maf() %>%
  ggplot(aes(VAF, fill=..x..)) +
  geom_histogram(bins = 100)+
  scale_fill_distiller(palette = "Blues") +
  ggtitle(paste("Variant Allele Fraction in Unfiltered Docetaxel Sensitive"))
## In order to check the actual distribution of the mutations, we should check
## the distribution in completely unfiltered data, including silent mutations 
## and mutations present in the sensitive samples (== before PoN filter)

## It can be useful to check SNP distribution in terms of Transitions vs Transversions with titv: 
## titv objects can be accessed for raw counts and normalized 

DRT.titv <- DRT_all %>% titv(plot = FALSE, useSyn = TRUE)
DST.titv <- DST_all %>% titv(plot = FALSE, useSyn = TRUE)

plotTiTv(DRT.titv, showBarcodes = TRUE)
plotTiTv(DST.titv, showBarcodes = TRUE)

DRT_unfilt.titv <- DRT_unfilt %>% titv(plot = FALSE, useSyn = TRUE)
DST_unfilt.titv <- DST_unfilt %>% titv(plot = FALSE, useSyn = TRUE)

plotTiTv(DRT_unfilt.titv, showBarcodes = TRUE)
plotTiTv(DST_unfilt.titv, showBarcodes = TRUE)


# VAF Thresholding and Heterogeneity Analysis -----------------------------


## After knowing this, we can look for the amount of genes or variants for each
## filtering level with vaf.table(). Indicate if you want genes or variants with
## "unique = " and thresholds with "thresholds = "

vaf.table(DRT_mafs,
          unique = "genes",
          thresholds = seq(0, 1, by = 0.05))
vaf.table(DRT_mafs, unique = "genes")


## Another way to look for a suitable threshold is to use Heterogeneity analyses.
## this will try to cluster mutations by VAF PER SAMPLE:

heterogeneity.10<-inferHeterogeneity(DRT_10, vafCol = 'VAF')

print(heterogeneity.10)
plotClusters(heterogeneity.10)

heterogeneity.10_vaf0.1 <-
  inferHeterogeneity(DRT_10, vafCol = 'VAF', minVaf = 0.1)
print(heterogeneity.10_vaf0.1)
plotClusters(heterogeneity.10_vaf0.1)

### We can manipulate the tumor labels to perform heterogeneity analysis on the
### whole sample:

heterogeneity.10_all <- DRT_10 %>%
  extract.maf() %>%
  mutate(Tumor_Sample_Barcode = "All_Samples") %>%
  read.maf() %>%
  inferHeterogeneity(vafCol = "VAF")

heterogeneity.10_all %>%
  plotClusters(
    genes = DRT_10@data %>%
      filter(!Hugo_Symbol == "UnknownGene") %>%
      arrange(desc(VAF)) %>%
      distinct(Hugo_Symbol, .keep_all = TRUE) %>%
      slice_max(n = 20, order_by = VAF) %>%
      pull(Hugo_Symbol)
  )

## In order to get a more accurate version, we can highlight the top 
## enriched genes in each sample.


for (i in 1:length(DRT_tumornames)) {
  plotClusters(
    heterogeneity.10,
    tsb = DRT_tumornames[i],
    genes = DRT_10 %>%
      extract.maf() %>%
      filter(
        VAF >= 0.2 &
          !Hugo_Symbol == "UnknownGene" &
          Tumor_Sample_Barcode == DRT_tumornames[i]
      ) %>%
      arrange(desc(VAF)) %>%
      distinct(Hugo_Symbol, .keep_all = TRUE) %>%
      slice_max(n = 10, order_by = VAF) %>%
      pull(Hugo_Symbol)
  )
}

for (i in 1:length(DRT_tumornames)) {
  plotClusters(
    heterogeneity.10_vaf0.1,
    tsb = DRT_tumornames[i],
    genes = DRT_10 %>%
      extract.maf() %>%
      filter(
        VAF >= 0.2 &
          !Hugo_Symbol == "UnknownGene" &
          Tumor_Sample_Barcode == DRT_tumornames[i]
      ) %>%
      arrange(desc(VAF)) %>%
      distinct(Hugo_Symbol, .keep_all = TRUE) %>%
      slice_max(n = 10, order_by = VAF) %>%
      pull(Hugo_Symbol)
  )
}


## At this point, we should select a particular threshold of both VAF and depth, 
## where we get a reasonable amount of genes to consider relevant, extract them
## or simply analyze from there on. This is a VAF-threshold-oriented method.
## E.g.: in https://doi.org/10.1371/journal.pone.0206632, average depth was 247; .
## Conchita's group used t_depth>=10, t_alt_depth>=4 and VAF>=0.10.


DRT_enriched_d25_v0.2 <- DRT_25 %>% subsetMaf(query = "VAF>=0.2")


# Clinical Enrichment Analysis -------------------------------------------------

## First we have to merge the MAFs from resistant and sensitive sources and append
## clinical data matching Tumor_Sample_Barcode. This for loop does it automatically for
## all mafs in one of our mafs vectors, but you can do it individually running the code inside.

## You can also provide "clinical_data.tsv" to clinicalEnrichment(annotationDat = ).


for (i in 1:length(DRT_mafs)) {
  temp <- merge_mafs(mafs = c(get(DRT_mafs[i]), get(DST_mafs[i]))) %>%
    extract.maf() %>%
    read.maf(clinicalData = "clinical_data.tsv")
  
  assign(x = str_replace(DRT_mafs, "DRT", "D.all")[i], value = temp)
  
  rm(temp, i)
  
}

all_tumornames<- D.all_all@clinical.data$Tumor_Sample_Barcode

ce.25 <- clinicalEnrichment(
  D.all_25,
  clinicalFeature = "DTX_Response",
  annotationDat = "clinical_data.tsv",
  minMut = 1
)

ce.25$groupwise_comparision %>% summary()
ce.25$groupwise_comparision[p_value <= 0.1]

ce.25$pairwise_comparision %>% summary()
ce.25$pairwise_comparision[fdr <= 0.1]

plotEnrichmentResults(ce.25, pVal = 0.1)

## This analysis suffers from the same issues as maf.compare(), where it only 
## compares presence in each specific sample rather than total # of mutations


# Raw mutation number comparison ------------------------------------------

## In order to compare 

# Signature analysis ------------------------------------------------------

## needs BSgenome and "BSgenome.Mmusculus.UCSC.mm10" installed. Install with
## BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")


## trinucleotideMatrix() gives an error when the start position is lower than 20.
## This is the case in our dataset [row 432487], so we need to filter it out.j
## The call to filter out error causing variants is included in the call.
## Creating the trinucleotide matrix for the samples we want to analyze.


tnm.DRT_all <-
  trinucleotideMatrix(
    DRT_all %>% subsetMaf(query = "Start_Position>20"),
    useSyn = TRUE,
    ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
  )

## Output from previous function will tell whether samples match APOBEC signature.
## If some are, you can plot them with this:

if (sum(tnm.DRT_all$APOBEC_scores$fdr < 0.05) > 0){
  
  tnm.DRT_all %>% plotApobecDiff()
  
}

## Estimation of number of signatures that adequate the most to the samples.
## Indicate nTry as high limit of signatures to try. Max = # Tumor_Sample_Barcodes
## end result is Cophenetic correlation graph: measures goodness of fit to each number of signatures.
## BE CAREFUL when setting parallel: check how many cores you have available.

tnm.DRT_all %>%
  estimateSignatures(nMin = 2,
                     nTry = 3,
                     nrun = 10,
                     parallel = 2) %>%
  plotCophenetic()

## Take the n with the lowest cophenetic metric and pass it as n = arg in extractSignatures
## This object contains the signatures created by NMF (non-negative matrix factorization):
## the TNM is decomposed into n signatures

sig.DRT_all <- tnm.DRT_all %>% extractSignatures(n = 3)

## From here on, we can plot and compare against legacy COSMIC signatures

plotSignatures(sig.DRT_all, sig_db = "legacy", yaxisLim = NA)

sig.DRT_all.orig <- sig.DRT_all %>%
  compareSignatures(sig_db = "legacy")

sig.DRT_all.orig$cosine_similarities %>%
  pheatmap(cluster_rows = FALSE, main = "Cosine similarity against original COSMIC signatures")

# Or the more modern SBS.v3 signatures

plotSignatures(sig.DRT_all, sig_db = "SBS", yaxisLim = NA)

sig.DRT_all.V3 <- sig.DRT_all %>%
  compareSignatures(sig_db = "SBS")

sig.DRT_all.V3$cosine_similarities %>%
  pheatmap(cluster_rows = FALSE, main = "Cosine similarity against SBS V3 signatures")


## We can try to do the same with a MAF including all samples. Ideally, we would fit 2 signatures 
## to match the clinical feature, but it's always better to check cophenetic correlation.

tnm.all <- trinucleotideMatrix(D.all_all %>% subsetMaf(query="Start_Position>20"),
                               useSyn = TRUE,
                               ref_genome = "BSgenome.Mmusculus.UCSC.mm10")

if (sum(tnm.all$APOBEC_scores$fdr < 0.05) > 0) {
  tnm.DRT_all %>% plotApobecDiff()
  
}

tnm.all %>% 
  estimateSignatures(nTry = 5, 
                     nrun = 10,
                     parallel = 2)
  
sig.2all <- tnm.all %>% extractSignatures(n=2)
sig.3all <- tnm.all %>% extractSignatures(n=3)

sig.2all_SBS <- sig.2all %>% compareSignatures(sig_db = "SBS")
sig.3all_SBS <- sig.3all %>% compareSignatures(sig_db = "SBS")

sig.2all %>% plotSignatures(sig_db = "SBS", yaxisLim = NA)
sig.3all %>% plotSignatures(sig_db = "SBS", yaxisLim = NA)

sig.2all_SBS$cosine_similarities %>% pheatmap(cluster_rows=FALSE, main = "Cosine Similarity against SBS V3 signatures")
sig.3all_SBS$cosine_similarities %>% pheatmap(cluster_rows=FALSE, main = "Cosine Similarity against SBS V3 signatures")



## Once this is done, we can get and print each tumor's contribution
## to the signature to see if signatures match clinical features.

sig.2all$contributions

sig2.all.contrib <- sig.2all$contributions %>%
  stack() %>%
  as.data.frame() %>%
  separate(
    col,
    sep = "_",
    remove = FALSE,
    into = c(NA, "DTX_response")
  )

sig2.all.contrib$tumor <-
  factor(sig2.all.contrib$col, levels = all_tumornames)

sig2.all.contrib %>%
  ggplot(aes(tumor, value, fill = DTX_response)) +
  geom_col() +
  facet_wrap( ~ row)

wilcox.test(value ~ DTX_response, data = sig2.all.contrib %>% filter(row == "Signature_1"))
wilcox.test(value ~ DTX_response, data = sig2.all.contrib %>% filter(row == "Signature_2"))



sig.3all$contributions

sig3.all.contrib <- sig.3all$contributions %>%
  stack() %>%
  as.data.frame() %>%
  separate(
    col,
    sep = "_",
    remove = FALSE,
    into = c(NA, "DTX_response")
  )

sig3.all.contrib$tumor <-
  factor(sig3.all.contrib$col, levels = all_tumornames)

sig3.all.contrib %>%
  ggplot(aes(tumor, value, fill = DTX_response)) +
  geom_col() +
  facet_wrap( ~ row)

wilcox.test(value ~ DTX_response, data = sig3.all.contrib %>% filter(row == "Signature_1"))
wilcox.test(value ~ DTX_response, data = sig3.all.contrib %>% filter(row == "Signature_2"))
wilcox.test(value ~ DTX_response, data = sig3.all.contrib %>% filter(row == "Signature_3"))
