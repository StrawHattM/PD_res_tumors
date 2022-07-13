#dependencies

library(tidyverse)

`%nin%` = Negate(`%in%`)


## assigning objects


CIS_files <- list.files("platinum_drugs/annotated_cns/", pattern="CIS", full.names = TRUE)
CARBO_files <- list.files("platinum_drugs/annotated_cns/", pattern="CARBO", full.names = TRUE)


for (i in 1:length(CIS_files)){
  
  temp<- read_delim(file=CIS_files[i])
  
  assign(paste0("CIS_",i),value= temp)
  
  rm(temp)
}


for (i in 1:length(CARBO_files)){
  
  temp<- read_delim(file=CARBO_files[i])
  
  assign(paste0("CARBO_",i),value= temp)

  rm(temp)
}

D2<- read_delim(list.files("platinum_drugs/annotated_cns/", pattern = "D2", full.names = TRUE))
D3<- read_delim(list.files("platinum_drugs/annotated_cns/", pattern = "D3", full.names = TRUE))

PD_control <- union(D2, D3)

# asigning and saving commons

CIScommon<- union(union(CIS_1[1:4], CIS_2[1:4]), CIS_3[1:4]) %>% 
  left_join(CIS_1, 
            by=c("chromosome", "start", "end", "gene")) %>%
  left_join(CIS_2, 
            by=c("chromosome", "start", "end", "gene"),
            suffix=c("", "_2")) %>% 
  left_join(CIS_3, 
            by=c("chromosome", "start", "end", "gene"),
            suffix=c("_1","_3"))  %>% 
  rowwise() %>%
  mutate(log2=mean(c(log2_1, log2_2, log2_3), na.rm=TRUE),
         cn=mean(c(cn_1, cn_2, cn_3), na.rm=TRUE),
         depth=mean(c(depth_1, depth_2, depth_3), na.rm=TRUE),
         p_bintest=mean(c(p_bintest_1, p_bintest_2, p_bintest_3), na.rm=TRUE),
         weight=mean(c(weight_1, weight_2, weight_3), na.rm=TRUE)) %>%
  select(1:4, log2, cn, depth, p_bintest, weight)

CIScommon_noCTRL<- setdiff(CIScommon[1:4], PD_control[1:4]) %>% 
  left_join(CIScommon, by=c("chromosome", "start", "end","gene"))


write_delim(CIScommon_noCTRL, file="custom_cns/PD/CIScommon_noCTRL.cns", delim="\t")

write_delim(CIScommon, file="custom_cns/PD/CIScommon.cns", delim="\t")




CARBOcommon<- union(union(CARBO_1[1:4], CARBO_2[1:4]), CARBO_3[1:4]) %>% 
  left_join(CARBO_1, 
            by=c("chromosome", "start", "end", "gene")) %>%
  left_join(CARBO_2, 
            by=c("chromosome", "start", "end", "gene"),
            suffix=c("", "_2")) %>% 
  left_join(CARBO_3, 
            by=c("chromosome", "start", "end", "gene"),
            suffix=c("_1","_3"))  %>% 
  rowwise() %>%
  mutate(log2=mean(c(log2_1, log2_2, log2_3), na.rm=TRUE),
         cn=mean(c(cn_1, cn_2, cn_3), na.rm=TRUE),
         depth=mean(c(depth_1, depth_2, depth_3), na.rm=TRUE),
         p_bintest=mean(c(p_bintest_1, p_bintest_2, p_bintest_3), na.rm=TRUE),
         weight=mean(c(weight_1, weight_2, weight_3), na.rm=TRUE)) %>%
  select(1:4, log2, cn, depth, p_bintest, weight)

CARBOcommon_noCTRL<- setdiff(CARBOcommon[1:4], PD_control[1:4]) %>% 
  left_join(CARBOcommon, by=c("chromosome", "start", "end","gene"))


write_delim(CARBOcommon_noCTRL, file="custom_cns/PD/CARBOcommon_noCTRL.cns", delim="\t")

write_delim(CARBOcommon, file="custom_cns/PD/CARBOcommon.cns", delim="\t")


## creating and assigning unions 
## 


CISunion<- union(union(CIS_1, CIS_2), CIS_3)


CISunion_noCTRL<- setdiff(CISunion[1:4], PD_control[1:4]) %>% 
  left_join(CISunion, by=c("chromosome", "start", "end","gene"))


write_delim(CISunion_noCTRL, file="custom_cns/PD/CISunion_noCTRL.cns", delim="\t")

write_delim(CISunion, file="custom_cns/PD/CISunion.cns", delim="\t")


CARBOunion<-union(CARBO_1, CARBO_2)

CARBOunion_noCTRL<- setdiff(CARBOunion[1:4], PD_control[1:4])%>% 
  left_join(CARBOunion, by=c("chromosome", "start", "end","gene"))


write_delim(CARBOunion_noCTRL, file="custom_cns/PD/CARBOunion_noCTRL.cns", delim="\t")

write_delim(CARBOunion, file="custom_cns/PD/CARBOunion.cns", delim="\t")



CISunion[1:4] %>% 
  setdiff(CARBOunion[1:4]) %>% 
  left_join(CISunion, by=c("chromosome", "start", "end","gene"))%>% 
  write_delim(file="custom_cns/PD/CIS_Exclusive.cns", delim="\t")

CARBOunion[1:4] %>% 
  setdiff(CISunion[1:4]) %>% 
  left_join(CARBOunion, by=c("chromosome", "start", "end","gene"))%>% 
  write_delim(file="custom_cns/PD/CARBO_Exclusive.cns", delim="\t")


CISunion_noCTRL[1:4] %>% 
  setdiff(CARBOunion_noCTRL[1:4]) %>% 
  left_join(CISunion_noCTRL, by=c("chromosome", "start", "end","gene"))%>% 
  write_delim(file="custom_cns/PD/CIS_noCTRL_Exclusive.cns", delim="\t")

CARBOunion_noCTRL[1:4] %>% 
  setdiff(CISunion_noCTRL[1:4]) %>% 
  left_join(CARBOunion_noCTRL, by=c("chromosome", "start", "end","gene"))%>% 
  write_delim(file="custom_cns/PD/CARBO_noCTRL_Exclusive.cns", delim="\t")