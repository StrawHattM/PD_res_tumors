###############################################################################
## PD resistant Tumors CNV analysis
## Author: Martín González Fernández
## Date: 11.07.2022
## Sven Rottenberg Lab
###############################################################################

# Dependencies ------------------------------------------------------------


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("karyoploteR", quietly = TRUE))
  BiocManager::install('karyoploteR')

if (!require("CopyNumberPlots", quietly = TRUE))
BiocManager::install('CopyNumberPlots')

library(tidyverse)
library(karyoploteR)
library(CopyNumberPlots)
library(ggfortify)
library(ggpubr)
library(reshape2)


# Import and Wrangling ----------------------------------------------------

cnr_PD <- list.files(path = "copy_number_ratios/", pattern ="PD", full.names = TRUE) 
cns_cbs_PD <- list.files(path = "called_cns/CBS/", pattern ="PD", full.names = TRUE)
cns_hmm_PD <- list.files(path = "called_cns/HMM/", pattern ="PD", full.names = TRUE)

names_PD <- cnr_PD %>% str_sub(20, -1) %>% str_replace("_PD.cnr", "")
names_PD

conds_PD <- c("SPLEEN", "NAIVE", "NAIVE", "CIS", "CIS", "CARBO", "CARBO", "CARBO", "CIS")
sample_PD <- names_PD %>% str_extract_all("D[1234567890]*") %>% unlist()

if (length(cnr_PD) != length(conds_PD)) {warning("length of conditions vector is different than CNR vector")}
if (length(cnr_PD) != length(sample_PD)) {warning("length of samples vector is different than CNR vector")}

for (i in 1:length(cnr_PD)){
  
  temp <- read_table(cnr_PD[i]) 
  cond <- conds_PD[i]
  sample <- sample_PD[i]
  temp <- temp %>% mutate(condition=cond, sample=sample, id=1:dim(temp)[1]) %>% relocate(id, 1)
  
  assign(x = paste0("cnr_",names_PD[i]), value = temp)
  
} 



# PCA plot ----------------------------------------------------------------

PD_log2 <-
  cnr_D2_VEHC %>%
  select(id, D2_VEHC = log2) %>%
  left_join(cnr_D3_NAIVE %>%
              select(id, D3_NAIVE = log2),
            by = "id") %>%
  left_join(cnr_D4_CIS %>%
              select(id, D4_CIS = log2),
            by = "id") %>%
  left_join(cnr_D5_CIS %>%
              select(id, D5_CIS = log2),
            by = "id") %>%
  left_join(cnr_D94_CIS %>%
              select(id, D94_CIS = log2),
            by = "id") %>%
  left_join(cnr_D6_CARBO %>%
              select(id, D6_CARBO = log2),
            by = "id") %>%
  left_join(cnr_D63_CARBO %>%
              select(id, D63_CARBO = log2),
            by = "id") %>%
  left_join(cnr_D7_CARBO %>%
              select(id, D7_CARBO = log2),
            by = "id") %>%
  select(-id) %>% t() %>% as.data.frame()

rownames(PD_log2)

anno.PD_log2 <- PD_log2 %>%
  rownames_to_column(var="id") %>%
  cbind(condition=c("VEHICLE", "NAIVE", rep("CIS", 3), rep("CARBO", 3))) %>% 
  relocate(condition, .after=id)

#"SPLEEN", "VEHICLE", "NAIVE",

prcomp(PD_log2) %>% autoplot(data=anno.PD_log2, colour="condition", label=TRUE, label.repel=TRUE, main="Platinum Resistant Tumors")
prcomp(PD_log2) %>% autoplot(data=anno.PD_log2, colour="condition", label=TRUE, label.repel=TRUE, x=2, y=3, main="Platinum Resistant Tumors, PCs 2&3")
prcomp(PD_log2) %>% autoplot(data=anno.PD_log2, colour="condition", label=TRUE, label.repel=TRUE, x=4, y=3, main="Platinum Resistant Tumors, PCs 4&3")


# Full dataset, log2 density plots and t-tests ----------------------------



## make the full dataset for the comparisons 
## 


PD_allcnr <- bind_rows(
  cnr_D2_VEHC,
  cnr_D3_NAIVE,
  cnr_D4_CIS,
  cnr_D5_CIS,
  cnr_D94_CIS,
  cnr_D6_CARBO,
  cnr_D7_CARBO,
  cnr_D63_CARBO
)

PD_allcnr %>% 
  ggplot(aes(x=factor(condition, levels = c("NAIVE","CIS", "CARBO")), y=log2, fill=condition)) +
  geom_boxplot() + 
  scale_y_continuous(limits = c(-7.5, 5))+
  scale_color_discrete() + 
  ggtitle("Log2 distribution per group", subtitle = "Total ") + xlab("Group")


PD_allcnr %>% 
  ggplot(aes(x=factor(condition, levels = c("NAIVE","CIS", "CARBO")), y=log2, fill=condition)) +
  geom_boxplot() + 
  scale_y_continuous(limits = c(-1, 1))+
  scale_color_discrete() + 
  ggtitle("Log2 distribution per group", subtitle = "Partial, -1 < log2 < 1") + xlab("Group")


PD_allcnr %>% 
  ggplot(aes(x=factor(condition, levels = c("NAIVE","CIS", "CARBO")), y=log2, color=condition)) +
  geom_jitter() + 
  scale_y_continuous(limits = c(-7.5, 5))+
  scale_color_discrete() + 
  ggtitle("Log2 distribution per group", subtitle = "Total ") + xlab("Group")





for (i in 1:dim(cnr_D10SRO212_RTST)[1]) {
  temp_seg <- RT_allcnr %>% filter(id == i)
  
  temp_ttest <- t.test(temp_seg$log2 ~ temp_seg$condition)
  
  
  temp_df <-
    cnr_D10SRO212_RTST %>%
    filter(id == i) %>%
    select(id, chromosome, start, end, gene) %>%
    mutate(
      t_stat = temp_ttest$statistic,
      t_p_val = temp_ttest$p.value,
      mean_DRT = temp_ttest$estimate[1],
      mean_DST = temp_ttest$estimate[2]
    )
  
  if (!exists("RT_allcnr_comparison")) {
    RT_allcnr_comparison <- temp_df
    
  } else {
    RT_allcnr_comparison <- bind_rows(RT_allcnr_comparison, temp_df)
    
  }
  
  print(i)
  
}