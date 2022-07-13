library(tidyverse)


# 53BP1 -------------------------------------------------------------------


### just changing the name in the original so that it matches the human gene name in gff 
a<- CARBO_unfilt %>% 
  extract.maf() %>% 
  mutate(Hugo_Symbol=case_when(
    Hugo_Symbol == "Trp53bp1" ~ "Tp53bp1")
  ) %>% 
  read.maf()

b<- CIS_unfilt %>% 
  extract.maf() %>% 
  mutate(Hugo_Symbol=case_when(
    Hugo_Symbol == "Trp53bp1" ~ "Tp53bp1")
  ) %>% 
  read.maf()


lollipopPlot2(m1= capHS(a), m1_name="CARBO_unfilt",
              m2= capHS(b), m2_name="CIS_unfilt",
              gene="TP53BP1")


PD_unfilt %>% extract.maf() %>% filter(Hugo_Symbol=="Trp53bp1") %>% write.csv(file = "53BP1_mutations_in_platdrugs.csv")



# Filtering ---------------------------------------------------------------

nm.PD_unfilt <- PD_unfilt %>% extract.maf()


## This is just to check how many empty elements are in each column, 
## thereby getting to know if they are good for checking duplicity or not.

empt.vector <- c()

for (col in nm.PD_unfilt%>%select(where(is.character))%>%colnames()) {
  
  testcol<-nm.PD_unfilt %>% pull(col) 
  
  tempnum<-sum(!nchar(testcol)) 
  
  empt.vector<-c(empt.vector, tempnum)
  
  rm(testcol, tempnum, col)
  
}

empty_by_columns <- data.frame("columns"=nm.PD_unfilt%>%select(where(is.character))%>%colnames(), "empty_entries"=empt.vector)
empty_by_columns

## I will avoid the ones with empty entries and Protein_Change (since I know there are forced pseudo-NAs in the Splice_Site variatns)



## IMPORTANT: duplicated() only marks the second instance of a duplicate, first one is FALSE. 
## This is intended for the base purpose of the function, filtering out duplicates while leaving the original
## This also means that if we filter based on duplicated(), we can obtain a dataset where one instance will mean
## one duplication, two will mean two duplications, so on and so forth.
## By concatenating duplicated() calls we can pinpoint how many replicates we want to allow


# Present in 2 tumors
dupes <-
  nm.PD_unfilt[duplicated(nm.PD_unfilt[, c("Hugo_Symbol", "Start_Position", "Tumor_Seq_Allele2")])]

# Present in 3 tumors
dupes2 <-
  dupes[duplicated(dupes[, c("Hugo_Symbol", "Start_Position", "Tumor_Seq_Allele2")])]

# Present in 4tumors
dupes3 <-
  dupes2[duplicated(dupes2[, c("Hugo_Symbol", "Start_Position", "Tumor_Seq_Allele2")])]


# Present in 5 tumors
dupes4 <-
  dupes3[duplicated(dupes3[, c("Hugo_Symbol", "Start_Position", "Tumor_Seq_Allele2")])]






table <- nm.PD_unfilt 

i=2

while(is_empty(table)==FALSE){
  
  duptab <- table %>% select("Hugo_Symbol", "Start_Position", "Tumor_Seq_Allele2") %>% duplicated()
  
  table <- table %>% mutate("duped"=duptab) %>% filter(duped==TRUE) 
  
  assign(x=paste0("table", i), value = table)
  
  i = i+1
  
}


duptab <- table %>% select("Hugo_Symbol", "Start_Position", "Tumor_Seq_Allele2") %>% duplicated()

table <- table %>% mutate("duplicated"=duptab) %>% filter(duped==TRUE) 

assign(x=paste0("table", i), value = table)

i = i+1


## Writing CSVs for Carmen. For the


dupes2 %>%
  write.csv(file = "Variants in +3 tumors.csv")

dupes2 %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  write.csv(file = "Genes in +3 tumors.csv")


nm.PD_unfilt %>%
  anti_join(dupes2) %>%
  write.csv(file = "PD_unfilt - Variants in 3 or more tumors filtered out.csv")

nm.PD_unfilt %>%
  anti_join(dupes2) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  write.csv(file = "PD_unfilt - Genes with variants in 3 or more tumors filtered out.csv")
