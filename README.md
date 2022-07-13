# PD_res_tumors: Platinum Drug Resistant Tumor analysis 

**Author**: Martín González Fernández  
**Group**: Sven Rottenberg, ITPA, University of Bern 

## Data origin

FASTQ and BAM files from Whole Genome Sequencing (WGS) of Carboplatin and Cisplatin resistant tumors, including controls.  
Grown by I. Klebic in ITPA, Bern; DNA isolation by C. Disler and I. Klebic, Sequencing by Novogene, Alignment and Variant Calling by Genevia. 

## QC and alignment 

Raw sequencing reads were initially trimmed using trimGalore (Krueger, 2021) (v. 0.5.0)  in paired-end mode and Cutadapt version 2.4(Martin, 2011), using a PHRED score cut-off of 20. After quality-based trimming, only sequence reads that were still over 20 base pairs in length and with an error rate inferior to 10% were retained for downstream analysis. The quality of the raw and trimmed reads was inspected using FastQC, version 0.11.9 (Andrews, 2010). All reads were then aligned against the FVB/NJ strain reference genome using BWA (Li & Durbin, 2009) (v. 0.7.17) using the default parameters for paired-end reads. Following alignment, PCR and optical duplicates were marked and the base quality scores were recalibrated using PICARD (Broad Institute) (v. 2.18.16) and GATK (van der Auwera et al., 2020) (v. 4.0.11.0), respectively.

## Variant calling

Somatic Mutations were called using GATK’s MuTect2 (Benjamin et al., 2019) using tumor-only parameters. A panel of normals was created using the baseline spleen (D1), naïve (D2) and vehicle-treated (D3) samples; which was used to filter out non-relevant mutations.
MuTect2 called mutations were annotated using Ensembl’s Variant Effect Predictor (VEP) (McLaren et al., 2016) (v. 100.2), and variants annotated as having a Moderate or High impact were retained. 

## Explorative Analysis of Variants

All downstream analysis was performed using R software (R Core Team, 2022) (v. 4.1.2). R package maftools (Mayakonda et al., 2018) was used for explorative analysis of all variant sets at different depth filtering levels (10x, 25x, 50x & 75x). Data wrangling was performed with tidyverse package collection (Wickham et al., 2019), and additional plots were generated with ggplot2 package (Wickham H., 2016). In order to explore the variants detected in the Carboplatin and Cisplatin resistant samples, we generated MAF summaries and oncoplots to explore the most frequent mutations. We also looked at the distribution of transitions and transversions. We also performed positional analysis using rainfall plots to detect kataegis, and detection of mutational clusters using OncodriveCLUSTL (Tamborero et al., 2013). We performed heterogeneity analysis and filtering based on the Variant Allele Fraction (VAF), further selecting based on VAF ____________.  Mutational signature analysis was performed on the samples using mouse genome UCSC release mm10, based on GRCm38.p6 in BSGenome (Bioconductor, 2021). A mutation matrix was to classify nucleotide substitutions into 96 classes based on surrounding bases was generated, and the number of signatures was estimated using non-negative matrix factorization (NMF) (Gaujoux & Seoighe, 2010) . The mutation matrix was decomposed into signatures, which where then compared to two previously characterised COSMIC signature datasets: the legacy 30 signatures and a more recent set of 65 Single Base Substitution (SBS) signatures (Tate et al., 2019). 


## Copy Number Alterations Analysis

Analysis of Copy Number Alterations (CNAs) was performed using CNVkit software (v. 0.9.9) (Talevich et al., 2016). Mouse FVB-NJ genomic sequence, track annotations and target genic regions were downloaded from Ensembl (Cunningham et al., 2022). Following the recommendations for WGS analysis, we limited the analysis to genic regions, we built a reference out of the spleen sample (D1). Then, we calculated bin-level CNAs against the reference, and performed segmentation on the bins using two different algorithms: Circular Binary Segmentation (cbs) (Olshen et al., 2011; Venkatraman & Olshen, 2007) and a Hidden Markov Model (hmm-tumor) (Schreiber, 2017). We compared the lists of segments to identify segments unique to each condition, and generated gene-level reports. Lastly, we generated scatter graphs highlighting the identified segments. 

## Acknowledgements

Calculations were performed on UBELIX (http://www.id.unibe.ch/hpc), the HPC cluster at the University of Bern.
