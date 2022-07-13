# Martín González Fernández,


#main package

# Initial dependencies and stuff ------------------------------------------


if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install("maftools")


#initial load

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


source("dependent_functions.R") #contains functions needed for the script to work. most are also here.


#path to TCGA LAML/BRCA MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
brca.path <-
  system.file("extdata", "brca.maf.gz", package = "maftools")
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
# getting domain list from maftools (useful)
gff = system.file("extdata", "protein_domains.RDs", package = "maftools")
gff = readRDS(file = gff)


# Some functions ----------------------------------------------------------



#function for turning Hugo Symbols inside of a MAF into uppercase. Necessary for some functions
capHS <- function(maf) {
  maf %>%
    subsetMaf(tsb = levels(maf@data$Tumor_Sample_Barcode),
              mafObj = F) %>%
    mutate(Hugo_Symbol = str_to_upper(Hugo_Symbol)) %>%
    read.maf()
}

#function to pipe MAF into tibble extraction

extract.maf <- function(maf) {
  subsetMaf(maf,
            tsb = levels(maf@data$Tumor_Sample_Barcode),
            mafObj = FALSE)
}


# reading the maf with read.maf. Examples are from TCGA
laml <- read.maf(maf = laml.maf, clinicalData = laml.clin)
brca = read.maf(maf = brca.path, verbose = FALSE)


# Assign loops and data cleanup -------------------------------------------


#Set strings containing the names of different MAFs to be analyzed together. Can be multiple.


mafs <-
  c("DRT_unfilt", "DRT_10", "DRT_25", "DRT_50", "DRT_75", "DRT_100")

# The assign call is customized for the MAF_combined.maf naming structure from Genevia.
# REQUIRES A /MAF_files FOLDER IN YOUR CURRENT WORKING DIRECTORY

for (maf in mafs) {
  assign(x = maf, read.maf(maf = paste0(
    "MAF_files/", paste(maf, "combined.maf", sep = "_")
  )))
  
}

#Equally, you can set different strings with the tumor names in each MAF or group of MAFS
tumornames <- levels(get(mafs[1])@data$Tumor_Sample_Barcode)

# This is a loop for cleaning of Genevia Data. Substitutes custom frameshift annotations
# for standardized annotations and adds Protein_Change column for lollipop plots.
# Remove the _clean at the end to just replace the original maf files.
# Requires add.protchange function.

add.protchange <- function(maf) {
  if (maf@maf.silent %>% rownames() %>% is_empty() == FALSE) {
    message("MAF file contains silent mutations that may not be displayed. Careful!!")
    
  }
  
  maf %>%
    subsetMaf(tsb = levels(maf@data$Tumor_Sample_Barcode),
              mafObj = FALSE) %>%
    separate(
      Amino_acids,
      sep = "/",
      into = c("ref.aa", "alt.aa"),
      remove = FALSE
    ) %>%
    mutate(Protein_position = str_replace_all(Protein_position, "-", "_")) %>%
    separate(
      Protein_position,
      sep = "_",
      into = c("pos1", "pos2"),
      remove = FALSE
    ) %>%
    mutate(
      Protein_Change = case_when(
        Variant_Classification == "In_Frame_Ins" &
          !ref.aa == "-" ~ paste0("p.", Protein_position, ref.aa, ">", alt.aa),
        Variant_Classification == "In_Frame_Ins" &
          ref.aa == "-" ~ paste0("p.", Protein_position, "ins", alt.aa),
        Variant_Classification == "In_Frame_Del" &
          !alt.aa == "-" ~ paste0("p.", Protein_position, ref.aa, ">", alt.aa),
        Variant_Classification == "In_Frame_Del" &
          alt.aa == "-" ~ paste0("p.", ref.aa, Protein_position, "del"),
        Variant_Classification == "Frame_Shift_Ins" ~ paste0("p.", ref.aa, pos1, "fs"),
        Variant_Classification == "Frame_Shift_Del" ~ paste0("p.", ref.aa, pos1, "fs"),
        Variant_Classification == "Missense_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Nonsense_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Nonstop_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Splice_Site" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Translation_Start_Site" ~ paste0("p.", ref.aa, Protein_position, alt.aa)
      )
    ) %>%
    select(!ref.aa & !alt.aa & !pos1 & !pos2) %>%
    read.maf()
}

for (maf in mafs) {
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
  
  tempname <- good.protchange(tempname)
  
  assign(x = paste0(maf), tempname)
  
  rm(maf)
}


# Intitial exploration. Summaries, oncoplots etc ----------------------------------------


# Initial exploration of a single MAF. Gonna use DRT_unfilt here just as an example

## Shows sample summary.
getSampleSummary(DRT_unfilt)

## Shows gene summary.
getGeneSummary(DRT_unfilt)

## shows clinical data associated with samples
getClinicalData(DRT_unfilt)

## Writes maf summary to an output file with basename DRT_unfilt. Can add a path to it.
### Also saves a .maf file if using a MAF that only exists in R workspace. Useful to save subsetted MAFs
write.mafSummary(maf = DRT_unfilt, basename = 'DRT_unfilt')


#You can also get information from different fields in @data (see them with getFields()) using classical descriptors:
## This makes you lose information on silent mutations since they are in @maf.silent instead of @data

DRT_unfilt@data$t_depth %>% summary() #summary of sequencing depth

DRT_unfilt@data$VAF %>% summary() #summary of VAF

DRT_unfilt@data$t_depth %>% hist(breaks = 100, main = "Total Depth")    #histogram of sequencing depth

DRT_unfilt@data$VAF %>% hist(breaks = 100, main = "VAF")    #histogram of sequencing depth

DRT_unfilt@data %>%
  filter(VAF > 0.2) %>%
  select(t_depth, VAF) %>%
  summary() #summary of seq. depth of samples with VAF>0.2

DRT_unfilt@data %>%
  filter(VAF > 0.2) %>%
  arrange(desc(VAF)) %>%
  ggplot(aes(t_depth, VAF)) +
  geom_point() +
  ggtitle("VAF ~ total depth") +
  geom_smooth() ## Plots VAF per total reads depth, just an example

#Main Summary plot, "top" mutated genes are the ones present in most number of samples
#e.g. plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

plotmafSummary(
  maf = DRT_unfilt,
  rmOutlier = TRUE,
  addStat = 'median',
  dashboard = TRUE,
  titvRaw = FALSE
)

#Oncoplots produce waterfall plots, better suited for multiple samples. Lots of customization options in doc.
oncoplot(maf = laml, top = 12)
oncoplot(maf = DRT_unfilt, top = 10)


# Transitions and transversion --------------------------------------------

#Plot and observe transitions and transversions

laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
# alternatively, set plot to TRUE or pipe it:
titv(maf = DRT_unfilt, useSyn = FALSE) %>% plotTiTv()



# plotVaf -----------------------------------------------------------------



#We can explore Variant Allele Frequency with plotVaf(), and modify its parameters:

plotVaf(
  DRT_25 %>% subsetMaf(query = "VAF>0.1"),
  vafCol = "VAF",
  genes = DRT_25 %>%
    extract.maf() %>%
    filter(VAF >= 0.2) %>%
    arrange(desc(VAF)) %>%
    distinct(Hugo_Symbol, .keep_all = TRUE) %>%
    filter(!Hugo_Symbol == "UnknownGene") %>%
    slice_max(n = 20, order_by = VAF) %>%
    pull(Hugo_Symbol),
  orderByMedian = TRUE
)

# lollipopPlot raw, fixed and 2 ----------------------------------------------



# lollipopPlot isplays mutation points of a particular gene throughout a MAF file.
# Needs to have the specific (non default) aminoacid change column specified

lollipopPlot(
  maf = laml,
  gene = 'TET2',
  AACol = 'Protein_Change',
  showMutationRate = TRUE
)

## lollipop plots also require a column that specifies the aminoacid mutation.
## I made this function to add that automatically, output is a MAF object with a
## "Protein_Change" (default)

add.protchange <- function(maf) {
  if (maf@maf.silent %>% rownames() %>% is_empty() == FALSE) {
    message("MAF file contains silent mutations that may not be displayed. Careful!!")
    
  }
  
  maf %>%
    subsetMaf(tsb = levels(maf@data$Tumor_Sample_Barcode),
              mafObj = FALSE) %>%
    separate(
      Amino_acids,
      sep = "/",
      into = c("ref.aa", "alt.aa"),
      remove = FALSE
    ) %>%
    mutate(Protein_position = str_replace_all(Protein_position, "-", "_")) %>%
    separate(
      Protein_position,
      sep = "_",
      into = c("pos1", "pos2"),
      remove = FALSE
    ) %>%
    mutate(
      Protein_Change = case_when(
        Variant_Classification == "In_Frame_Ins" &
          !ref.aa == "-" ~ paste0("p.", Protein_position, ref.aa, ">", alt.aa),
        Variant_Classification == "In_Frame_Ins" &
          ref.aa == "-" ~ paste0("p.", Protein_position, "ins", alt.aa),
        Variant_Classification == "In_Frame_Del" &
          !alt.aa == "-" ~ paste0("p.", Protein_position, ref.aa, ">", alt.aa),
        Variant_Classification == "In_Frame_Del" &
          alt.aa == "-" ~ paste0("p.", ref.aa, Protein_position, "del"),
        Variant_Classification == "Frame_Shift_Ins" ~ paste0("p.", ref.aa, pos1, "fs"),
        Variant_Classification == "Frame_Shift_Del" ~ paste0("p.", ref.aa, pos1, "fs"),
        Variant_Classification == "Missense_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Nonsense_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Nonstop_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Splice_Site" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Translation_Start_Site" ~ paste0("p.", ref.aa, Protein_position, alt.aa)
      )
    ) %>%
    select(!ref.aa & !alt.aa & !pos1 & !pos2) %>%
    read.maf()
}

## When using mouse genes (lowercase Hugo Symbol) lollipopPlot doesn't work properly.
## I created lollipopPlot_fixed to fix this and added the dependent functions in
## a file called "dependent_functions.R" that needs to be in work dir and sourced.
## Bear in mind that the domains displayed will be human, not mouse.

source("dependent_functions.R")

# starting from a raw DRT_25, for example
DRT_25 %>%
  add.protchange() %>%
  subsetMaf(query = "VAF>0.2") %>%
  lollipopPlot_fixed(
    gene = "Pgr",
    labelPos = "all",
    labelOnlyUniqueDoamins = FALSE
  )

# Alternatively, you can use capHS() to capitalyze hugo symbols and query in caps:

DRT_25 %>%
  capHS() %>%
  lollipopPlot(gene = "PGR")


#It's also possible to use lollipopPlot2 to plot mutations in 2 different mafs:

lollipopPlot2(
  m1 = DRT_unfilt %>% capHS() ,
  m2 = DST_unfilt %>% capHS(),
  m2_name = "Sensitive",
  m1_name = "Resistant",
  gene = "PGR"
)


# rainfall plots ----------------------------------------------------------



# rainfallPlot computes intermutation distance to detect mutation hotspots (kataegis).
# apart from displaying it in the console, it will make a tsv file in pwd() with kataegis information
rainfallPlot(maf = brca,
             detectChangePoints = TRUE,
             pointSize = 0.4)

## automatically, it will compute the first sample. If you have multiple tumor barcodes, adjust and run the loop

for (tumor in DRT_tumornames) {
  rainfallPlot(maf = DRT_10 %>% subsetMaf(query="VAF>=0.05"),
               detectChangePoints = T,
               tsb = tumor)
}


# tcgaCompare -------------------------------------------------------------




#tcgaCompare allows us to see if the mutational load/burden is similar to that of a tcga cancer type

comparison.tcga.unfilt <-
  tcgaCompare(maf = DRT_unfilt%>% subsetMaf(query="t_alt_count>2"),
              cohortName = "DRT_unfilt",
              logscale = T)

tcgaCompare(maf = DST_unfilt%>% subsetMaf(query="t_depth>=20"),
            cohortName = "DST_unfilt",
            logscale = T)


## We can then extract which cohorts are not detected to be different (~equal) by a t-test.


comparison.tcga.unfilt$pairwise_t_test %>%
  filter(Cohort1 == "DoceRes" | Cohort2 == "DoceRes") %>%
  filter(Pval > 0.05)


# Heterogeneity analysis --------------------------------------------------



# heterogeneity analysis

heterogeneity.25 <-
  inferHeterogeneity(DRT_25, tsb = DRT_tumornames, vafCol = 'VAF')
print(heterogeneity.25)

### plots the clusters highlighting top 10 enriched genes in any sample
plotClusters(
  heterogeneity.25,
  genes = DRT_25@data %>%
    filter(VAF >= 0.2) %>%
    arrange(desc(VAF)) %>%
    distinct(Hugo_Symbol, .keep_all = TRUE) %>%
    slice_max(n = 10, order_by = VAF) %>%
    pull(Hugo_Symbol)
)



# VAF tables --------------------------------------------------------------




# Functions to get total number of genes/variants for each depth and VAF filtering level:
## Both functions require an argument maflist (character string with names of MAF objects)
## Optional argument thresholds requires numerical vector,  defaults from 0 to 0.9

#vaf.Variants and vaf.Genes superseded by vaf.table

vaf.Genes <-
  function(maflist,
           thresholds = seq(0, 0.9, by = 0.1),
           basename = NULL) {
    #Wrapper to create the temp vector with results for each MAF in maflist
    for (maf in maflist) {
      #create empty vector with name for temp list. 1/maf
      assign(paste("genes", maf, sep = "_"), c())
      
      #applies to the maf the same function with different VAF thresholds
      for (vaf in thresholds) {
        # "a" is a temporary variable that serves to store the length output for each VAF threshold in each cycle
        a <- get(maf)@data %>%
          filter(VAF >= vaf) %>%
          pull("Hugo_Symbol") %>%
          unique() %>%
          length()
        
        # glues "a" at the end of the specific maf vector
        assign(x = paste("genes", maf, sep = "_"),
               value = c(get(paste(
                 "genes", maf, sep = "_"
               )), a))
        
      }
    }
    
    #creating an identifiable vector of the temp vectors names
    vafgenestablenames <- paste("genes", maflist, sep = "_")
    
    #skeleton (empty df) of the final table
    vaf.gene.table <- data.frame()
    
    #binds one by one the temp vectors into a data frame, rowwise in the order of maflist
    for (i in 1:length(vafgenestablenames)) {
      vaf.gene.table <- rbind(vaf.gene.table, get(vafgenestablenames[i]))
    }
    
    #setting dimnames for the final table
    rownames(vaf.gene.table) <- maflist
    colnames(vaf.gene.table) <- paste0("VAF>", thresholds)
    
    if (!is.null(basename)) {
      vaf.gene.table %>% write.csv(file = paste0(basename, "_VAFtable.csv"))
    }
    
    #return
    vaf.gene.table
  }
vaf.Variants <-
  function(maflist,
           thresholds = seq(0, 0.9, by = 0.1),
           basename = NULL) {
    for (maf in maflist) {
      assign(paste("variants", maf, sep = "_"), c())
      for (vaf in thresholds) {
        a <- get(maf)@data %>%
          filter(VAF >= vaf) %>%
          pull("Hugo_Symbol") %>%
          length()
        
        assign(x = paste("variants", maf, sep = "_"),
               value = c(get(paste(
                 "variants", maf, sep = "_"
               )), a))
      }
    }
    
    vafvariantstablenames <- paste("variants", maflist, sep = "_")
    vaf.variant.table <- data.frame()
    
    for (i in 1:length(vafvariantstablenames)) {
      vaf.variant.table <-
        rbind(vaf.variant.table, get(vafvariantstablenames[i]))
    }
    
    rownames(vaf.variant.table) <- maflist
    colnames(vaf.variant.table) <- paste0("VAF>", thresholds)
    
    if (!is.null(basename)) {
      vaf.variant.table %>% write.csv(file = paste0(basename, "_VAFtable.csv"))
    }
    
    vaf.variant.table
  }

vaf.table <-
  function(maflist,
           thresholds = seq(0, 0.9, by = 0.1),
           unique = "genes",
           basename = NULL) {
    if (!unique == "genes" & !unique == "variants") {
      stop(
        "Please specify a correct value for \"unique\" argument. default=\"genes\", also admits \"variants\"."
      )
    }
    else {
      if (unique == "genes") {
        #Wrapper to create the temp vector with results for each MAF in maflist
        for (maf in maflist) {
          #create empty vector with name for temp list. 1/maf
          assign(paste("genes", maf, sep = "_"), c())
          
          #applies to the maf the same function with different VAF thresholds
          for (vaf in thresholds) {
            # "a" is a temporary variable that serves to store the length output for each VAF threshold in each cycle
            a <- get(maf)@data %>%
              filter(VAF >= vaf) %>%
              pull("Hugo_Symbol") %>%
              unique() %>%
              length()
            
            # glues "a" at the end of the specific maf vector
            assign(x = paste("genes", maf, sep = "_"),
                   value = c(get(paste(
                     "genes", maf, sep = "_"
                   )), a))
            
          }
        }
        
        #creating an identifiable vector of the temp vectors names
        vafgenestablenames <- paste("genes", maflist, sep = "_")
        
        #skeleton (empty df) of the final table
        vaf.gene.table <- data.frame()
        
        #binds one by one the temp vectors into a data frame, rowwise in the order of maflist
        for (i in 1:length(vafgenestablenames)) {
          vaf.gene.table <- rbind(vaf.gene.table, get(vafgenestablenames[i]))
        }
        
        #setting dimnames for the final table
        rownames(vaf.gene.table) <- maflist
        colnames(vaf.gene.table) <- paste0("VAF>", thresholds)
        
        if (!is.null(basename)) {
          vaf.gene.table %>% write.csv(file = paste0(basename, "_VAFtable.csv"))
        }
        
        #return
        return(vaf.gene.table)
        
      }
      else{
        if (unique == "variants") {
          for (maf in maflist) {
            assign(paste("variants", maf, sep = "_"), c())
            for (vaf in thresholds) {
              a <- get(maf)@data %>%
                filter(VAF >= vaf) %>%
                pull("Hugo_Symbol") %>%
                length()
              
              assign(x = paste("variants", maf, sep = "_"),
                     value = c(get(
                       paste("variants", maf, sep = "_")
                     ), a))
            }
          }
          
          vafvariantstablenames <- paste("variants", maflist, sep = "_")
          vaf.variant.table <- data.frame()
          
          for (i in 1:length(vafvariantstablenames)) {
            vaf.variant.table <-
              rbind(vaf.variant.table, get(vafvariantstablenames[i]))
          }
          
          rownames(vaf.variant.table) <- maflist
          colnames(vaf.variant.table) <- paste0("VAF>", thresholds)
          
          if (!is.null(basename)) {
            vaf.variant.table %>% write.csv(file = paste0(basename, "_VAFtable.csv"))
          }
          
          vaf.variant.table
          
        }
      }
    }
  }


#get the tables

vaf.table(mafs)
vaf.table(mafs, unique = "variants")

if (mean(vaf.table(mafs) <=  vaf.table(mafs, unique = "variants")) != 1) {
  stop("possible error in vaf.table function. More genes altered than total variants")
}


#Subsetting into a new maf to get only variants above a desired VAF.
## You can just not assign them and pipe in front for quick exploration
## mafObj =TRUE coerces into a MAF class object | FALSE makes a df inluding maf.silent (unlike just @data the MAF)
## maftools functions require MAF class, other things require df. Choose accordingly.

DRT_25_0.2 <- subsetMaf(maf = DRT_25,
                        query = "VAF>0.2",
                        mafObj = TRUE)

### Writing into a CSV, mafObj = false to use write.csv after from the df. Use this instead of @data to also include silent mutations.
subsetMaf(maf = DRT_25,
          query = "VAF>0.2",
          mafObj = FALSE) %>%
  write.csv(file = "DRT_25depth_0.2vaf.csv")



# subsetting mafs now for some reason -------------------------------------



### setting mafObj=TRUE allows us to extract specific populations from a MAF file and apply maftools functions to it:
### following code subsets based on different parameters and applies different maftools functions, as examples.

subsetMaf(maf = DRT_25, query = "VAF>0.2") %>%
  plotmafSummary(addStat = 'median',
                 dashboard = TRUE,
                 titvRaw = FALSE)

subsetMaf(maf = DRT_25, query = "VAF>=0.2") %>%
  oncoplot(
    leftBarLims = c(0, 1),
    leftBarData = DRT_25@data %>%
      filter(VAF >= 0.2) %>%
      select(Hugo_Symbol, VAF) %>%
      distinct(Hugo_Symbol, .keep_all = TRUE)
  )

subsetMaf(
  DRT_unfilt,
  mafObj = T,
  genes = DRT_25@data %>%
    filter(VAF >= 0.1) %>%
    select(Hugo_Symbol, VAF) %>%
    slice_max(n = 40, order_by = VAF) %>%
    distinct(Hugo_Symbol, .keep_all = FALSE) %>%
    pull(Hugo_Symbol)
) %>%
  plotVaf(vafCol = "VAF")





for (i in 1:length(DRT_mafs)) {
  temp <- get(DRT_mafs[i]) %>%
    extract.maf() %>%
    ggplot(aes(VAF)) +
    geom_density() +
    ggtitle(paste("Density Plot of Variant Allele Fraction in", DRT_mafs[i]))
  
  assign(x = paste0("a", i), value = temp)
}

multiplot(a1, a2, a3, a4, a5, a6)
rm(a1, a2, a3, a4, a5, a6)



# somatic interactions ----------------------------------------------------



# Somatic interactions require a higher number of samples to work, so don't use.
## can be marginally useful to make coocurrence/exclusive lists and compute event ratio

somaticInteractions(maf = DRT_100,
                    top = 25,
                    pvalue = c(0.05, 0.1))




# OncodriveCLUSTL ---------------------------------------------------------


# Oncodrive can be used to check for clustered mutations in cancer driver genes.
## It assumes that mutations in driver genes will be clustered in specific domains
## Requires Protein_Change column and uppercase Hugo Symbols.
## Also recommended to play around with minMut to establish minimal mutation #

oncodrive(
  maf = DRT_100,
  AACol = "Protein_Change",
  pvalMethod = "zscore",
  minMut = 1
) #gives error

od.unfilt <- DRT_unfilt %>%
  capHS() %>%
  oncodrive(minMut = 3, AACol = "Protein_Change")

od.unfilt %>% plotOncodrive(bubbleSize = 0.5) #not super useful with 3 samples but serves to see distribution

## we can also extract a list of genes significant in oncodrive analysis

clustered <- od.unfilt %>%
  filter(fdr <= 0.05) %>%
  pull(Hugo_Symbol) %>%
  str_to_title()

clustered

# and use it in further analyses

DRT_unfilt %>%
  subsetMaf(genes = clustered,
            mafObj = TRUE) %>%
  plotmafSummary(addStat = "median", showBarcodes = TRUE)

DRT_unfilt %>%
  subsetMaf(genes = clustered,
            mafObj = TRUE) %>%
  oncoplot(
    leftBarLims = c(0, 1),
    leftBarData = DRT_25@data %>%
      filter(Hugo_Symbol %in% clustered) %>%
      select(Hugo_Symbol, VAF) %>%
      arrange(desc(VAF)) %>%
      distinct(Hugo_Symbol, .keep_all = TRUE)
  )


DRT_unfilt %>%
  subsetMaf(genes = clustered) %>%
  plotVaf(
    vafCol = "VAF",
    orderByMedian = TRUE,
    showN = TRUE,
    top = length(clustered)
  )


# pfamDomains (doesn't work) ----------------------------------------------



# pfamDomains allows to add pfam domain information to aminoacid changes.
# Also summarize aa changes into affected domains.
## Doesn't work so far, i have no clue why.

DRT_unfilt %>%
  capHS() %>%
  pfamDomains(DRT_100, top = 10)




# Section to know how to make add.protchange generate standardized --------


varclass <- laml %>%
  subsetMaf(tsb = levels(laml@data$Tumor_Sample_Barcode),
            mafObj = F) %>%
  pull(Variant_Classification) %>%
  levels()

# Getting a sample from laml (canonical)
laml.slice <- data.frame()
for (var in varclass) {
  a <- laml@data %>%
    filter(Variant_Classification == var) %>%
    slice_head(n = 10)
  
  laml.slice <- rbind(laml.slice, a)
  rm(var)
}
laml.slice <-
  laml.slice %>% select(Variant_Classification, Protein_Change)
view(laml.slice)

# Getting a sample from DRT_unfilt (to adjust)
DRT.slice <- data.frame()
for (var in varclass) {
  a <- DRT_unfilt@data %>%
    filter(Variant_Classification == var) %>%
    slice_head(n = 10)
  
  DRT.slice <- rbind(DRT.slice, a)
  rm(var)
}
DRT.slice <-
  DRT.slice %>% select(Variant_Classification, Protein_Change)
view(DRT.slice)


## experiment bench for case_when changes to add.protchang.
good.protchange <- function(maf) {
  if (maf@maf.silent %>% rownames() %>% is_empty() == FALSE) {
    message("MAF file contains silent mutations that may not be displayed. Careful!!")
    
  }
  
  maf %>%
    subsetMaf(tsb = levels(maf@data$Tumor_Sample_Barcode),
              mafObj = FALSE) %>%
    separate(
      Amino_acids,
      sep = "/",
      into = c("ref.aa", "alt.aa"),
      remove = FALSE
    ) %>%
    mutate(Protein_position = str_replace_all(Protein_position, "-", "_")) %>%
    separate(
      Protein_position,
      sep = "_",
      into = c("pos1", "pos2"),
      remove = FALSE
    ) %>%
    mutate(
      Protein_Change = case_when(
        Variant_Classification == "In_Frame_Ins" &
          !ref.aa == "-" ~ paste0("p.", Protein_position, ref.aa, ">", alt.aa),
        Variant_Classification == "In_Frame_Ins" &
          ref.aa == "-" ~ paste0("p.", Protein_position, "ins", alt.aa),
        Variant_Classification == "In_Frame_Del" &
          !alt.aa == "-" ~ paste0("p.", Protein_position, ref.aa, ">", alt.aa),
        Variant_Classification == "In_Frame_Del" &
          alt.aa == "-" ~ paste0("p.", ref.aa, Protein_position, "del"),
        Variant_Classification == "Frame_Shift_Ins" ~ paste0("p.", ref.aa, pos1, "fs"),
        Variant_Classification == "Frame_Shift_Del" ~ paste0("p.", ref.aa, pos1, "fs"),
        Variant_Classification == "Missense_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Nonsense_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Nonstop_Mutation" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Splice_Site" ~ paste0("p.", ref.aa, Protein_position, alt.aa),
        Variant_Classification == "Translation_Start_Site" ~ paste0("p.", ref.aa, Protein_position, alt.aa)
      )
    ) %>%
    select(!ref.aa & !alt.aa & !pos1 & !pos2) %>%
    read.maf()
}

check <- DRT_unfilt %>% good.protchange()
check@data %>% filter(Variant_Classification == "In_Frame_Ins") %>% View()



# drug interactions and druggable genome ----------------------------------



# drugInteractions checks for potentially druggable genes, displaying by categories

DRT_25 %>%
  capHS() %>%
  subsetMaf(query = "VAF>0.1") %>%
  drugInteractions()

## You can also extract specific information about which drugs have been reportedly associated
## with each gene, using drugs=TRUE

res25_0.1.dgi <- DRT_25 %>%
  capHS() %>%
  drugInteractions(
    genes = DRT_25 %>%
      extract.maf() %>%
      filter(VAF > 0.1) %>%
      distinct(Hugo_Symbol) %>%
      pull(Hugo_Symbol) %>%
      str_to_upper(),
    drugs = TRUE
  )
res25_0.1.dgi



# maf/mutCompare: Comparison of 2 MAF files -------------------------------


# Comparison of two different MAF files.

# depth 25
dtxcomp.25 <- mutCompare(
  m1 = DRT_25,
  m2 = DST_25,
  m2Name = "Sensitive",
  m1Name = "Resistant",
  minMut = 3
)
dtxcomp.25
dtxcomp.25 %>% forestPlot(pVal = 0.01, geneFont = 1)


dtxcomp.25_res <-
  list(
    results = dtxcomp.25$results %>% filter(Resistant > Sensitive),
    SampleSummary = dtxcomp.25$SampleSummary
  )
dtxcomp.25_res %>% forestPlot(pVal = 0.05, geneFontSize = 0.6)

#depth 10
dtxcomp.10 <- mutCompare(
  m1 = DRT_10,
  m2 = DST_10,
  m2Name = "Sensitive",
  m1Name = "Resistant",
  minMut = 3
)
dtxcomp.10
dtxcomp.10 %>% forestPlot(pVal = 0.01, geneFont = 1)


dtxcomp.10_res <-
  list(
    results = dtxcomp.10$results %>% filter(Resistant > Sensitive ),
    SampleSummary = dtxcomp.10$SampleSummary
  )
dtxcomp.10_res %>% forestPlot(pVal = 0.05, geneFontSize = 0.5)



# Signature analysis  -----------------------------------------------------




library(BiocManager)
if (!require("BiocManager")) {
  install.packages("BSgenome")
}
library(BSgenome)
library("BSgenome.Mmusculus.UCSC.mm10")

library(NMF)

# needs BSgenome and "BSgenome.Mmusculus.UCSC.mm10" installed. Install with
# install("BSgenome.Mmusculus.UCSC.mm8")


## trinucleotideMatrix() gives an error when the start position is lower than 20.
## This is the case in our dataset [row 432487], so we need to filter it out.

a <- DRT_all %>% extract.maf() %>% filter(Variant_Type == "SNP")
a[432487, ]

## Creating the trinucleotide matrix for the samples we want to analyze.
## The call to filter out error causing variants is included.

tnm.DRT_all.8 <-
  trinucleotideMatrix(
    DRT_all %>% subsetMaf(query = "Start_Position>20"),
    useSyn = TRUE,
    ref_genome = "BSgenome.Mmusculus.UCSC.mm8"
  )
tnm.DRT_all.9 <-
  trinucleotideMatrix(
    DRT_all %>% subsetMaf(query = "Start_Position>20"),
    useSyn = TRUE,
    ref_genome = "BSgenome.Mmusculus.UCSC.mm9"
  )
tnm.DRT_all.10 <-
  trinucleotideMatrix(
    DRT_all %>% subsetMaf(query = "Start_Position>20"),
    useSyn = TRUE,
    ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
  )
## Output from previous function will tell whether samples match APOBEC signature.
## If some are, you can plot them with this:

if (sum(tnm.DRT_all$APOBEC_scores$fdr < 0.05) > 0) {
  tnm.DRT_all %>% plotApobecDiff()
  
}

## Estimation of number of signatures that adequate the most to the samples.
## Indicate nTry as high limit of signatures to try. Max = # Tumor_Sample_Barcodes
## end result is Cophenetic correlation graph: measures goodness of fit to each number of signatures.
## BE CAREFUL when setting parallel: check how many cores you have available.
for(i in 8:10){
  get(paste0("tnm.DRT_all.", i)) %>%
    estimateSignatures(nMin = 2,
                       nTry = 3,
                       parallel = 2)
  mtext(paste0("Signatures against mm",i))
}
## Take the n with the lowest cophenetic metric and pass it as n = arg in extractSignatures
## This object contains the signatures created by NMF (non-negative matrix factorization):
## the TNM is decomposed into n signatures

sig.DRT_all <- tnm.DRT_all %>% extractSignatures(n = 3)


for(i in 8:10){
  
  tempsig<- get(paste0("tnm.DRT_all.", i)) %>%
    extractSignatures(n=2)
  
  assign(x=paste0("sig.DRT_all.", i), value=tempsig)
 
  plotSignatures(get(paste0("sig.DRT_all.", i)), sig_db = "legacy", yaxisLim = NA)
  plotSignatures(get(paste0("sig.DRT_all.", i)), sig_db = "SBS", yaxisLim = NA)
  
  tempcomp.l<- tempsig %>% compareSignatures(sig_db="legacy")
  
  tempcomp.l$cosine_similarities %>% 
    pheatmap(cluster_rows = FALSE, main = paste0("Cosine similarities against original signatures, mm", i))
  
  tempcomp.S<- tempsig %>% compareSignatures(sig_db="SBS")
  
  tempcomp.S$cosine_similarities %>% 
    pheatmap(cluster_rows = FALSE, main = paste0("Cosine similarities against SBS signatures, mm", i))
  
  rm(tempsig, tempcomp.l, tempcomp.S)
}


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
