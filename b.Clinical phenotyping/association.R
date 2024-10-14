rm(list=ls())
gc()


# Call libraries needed
library(here)
library(dplyr)
library(tidyverse)
library(foreign)
library(data.table)
library(glmtoolbox)
library(ipw)
library(marginaleffects)

# Set working directory
here::i_am("script/b.Clinical phenotyping/association.R")
setwd(here::here())

# Load datasets
fusion_all<-readRDS("./results/fusion_all.rds")
helix_clinical<-readRDS("./data/helix_clinical.rds")

snf<-fusion_all[[paste0("id_northern_southern")]]
snf<-snf %>% dplyr::mutate(Cluster = ifelse(newgroups==3, "Cluster C",
                                        ifelse(newgroups==2, "Cluster A", 
                                               ifelse(newgroups==1, "Cluster B",NA)))) %>%
  dplyr::mutate(Cluster=as.factor(Cluster)) %>%
  dplyr::select(-c("cohort","newgroups"))

all_clinical<-snf %>%
  dplyr::left_join(helix_clinical, by="HelixID") 
  

# Associations with clinical phenotypes
source("script/b.Clinical phenotyping/!function_association.R")

# definition of outcomes
cont_z_outc<-c("hs_zbmi","hs_zwaist_mets", "hs_z_fatmass",
               "hs_zhdl_mets","hs_ztriglyc_mets","hs_zinsulin_mets",
               "hs_zsysbp_mets","hs_zdiabp_mets", "hs_z_alt",
               "hs_metsscore")

cat_outc<-c("hs_ovob", "metab_health")

all_clinical_nw<- all_clinical %>% dplyr::filter(region=="n/w")
northern<-association("Cluster", all_clinical_nw)
write.csv(northern, "./results/csvs/nothern_association_clinical.csv")

all_clinical_sm<- all_clinical %>% dplyr::filter(region=="s/m")
southern<-association("Cluster",all_clinical_sm)
write.csv(southern, "./results/csvs/southern_association_clinical.csv")

pooled<-association("Cluster", all_clinical)
write.csv(pooled, "./results/csvs/pooled_association_clinical.csv")

# stratification by obesity status
cont_z_outc<-c("hs_zwaist_mets", "hs_z_fatmass",
               "hs_zhdl_mets","hs_ztriglyc_mets","hs_zinsulin_mets",
               "hs_zsysbp_mets","hs_zdiabp_mets", "hs_z_alt",
               "hs_metsscore")
cat_outc<-c("metab_health")

# Overweight/obesity
obesity_northern<-association_strat("Cluster",all_clinical_nw, 1)
write.csv(obesity_northern, "./results/csvs/obesity_strat_northern.csv")
obesity_southern<-association_strat("Cluster",all_clinical_sm, 1)
write.csv(obesity_southern, "./results/csvs/obesity_strat_southern.csv")
obesity_pooled<-association_strat("Cluster",all_clinical,1)
write.csv(obesity_pooled, "./results/csvs/obesity_strat_pooled.csv")

# Normal weight
normalweight_northern<-association_strat("Cluster",all_clinical_nw, 0)
write.csv(normalweight_northern, "./results/csvs/normalweight_strat_northern.csv")
normalweight_southern<-association_strat("Cluster",all_clinical_sm, 0)
write.csv(normalweight_southern, "./results/csvs/normalweight_strat_southern.csv")
normalweight_pooled<-association_strat("Cluster",all_clinical,0)
write.csv(normalweight_pooled, "./results/csvs/normalweight_strat_pooled.csv")


# Association with adolescent anthropometry
adolescent<-readRDS("./data/adolescent.rds") 

adolescent<-adolescent %>%
  dplyr::right_join(all_clinical, by="HelixID") %>%
  dplyr::mutate(attrition=as.factor(ifelse(!is.na(hs2_age),1,
                                           ifelse(is.na(hs2_age),0,NA))))
table(adolescent$attrition)

# Create attrition weights
library(WeightIt)
library(cobalt)

# Northern/Western cohort
adolescent_nw<-adolescent%>% dplyr::filter(region=="n/w") %>% 
  dplyr::select(HelixID,attrition, cohort,h_sex,ethnicity_3cat,h_edumc_None,hs_ovob,hs_age_years,h_age_None, region)
table(adolescent_nw$attrition)

weights_weightit_nw <- weightit(attrition ~ cohort+h_sex+ethnicity_3cat+h_edumc_None+hs_ovob+hs_age_years,
                                data = as.data.frame(adolescent_nw), method="glm")

summary(weights_weightit_nw)
bal.tab(weights_weightit_nw, stats = c("m", "v"), thresholds = c(m = .05)) # all balanced

adolescent_nw$weights_atr<-weights_weightit_nw$weights
adolescent_nw_weight<-adolescent_nw %>%
  dplyr::select(HelixID, weights_atr)

# Southern/Mediterranean cohort
adolescent_sm<-adolescent%>% dplyr::filter(region=="s/m") %>% 
  dplyr::select(HelixID,attrition, cohort,h_sex,ethnicity_3cat,h_edumc_None,hs_ovob,hs_age_years,h_age_None)
table(adolescent_sm$attrition)

weights_weightit_sm <- weightit(attrition ~ cohort+h_sex+h_edumc_None+hs_ovob+hs_age_years,
                                data = as.data.frame(adolescent_sm), method="glm")

summary(weights_weightit_sm)
bal.tab(weights_weightit_sm, stats = c("m", "v"), thresholds = c(m = .05)) # all balanced

adolescent_sm$weights_atr<-weights_weightit_sm$weights
adolescent_sm_weight<-adolescent_sm %>%
  dplyr::select(HelixID, weights_atr)

# combine both weight datasets
weights_all<-rbind(adolescent_nw_weight,adolescent_sm_weight)

# combine initial dataset with weights
adolescent<-adolescent %>%
  dplyr::left_join(weights_all, by="HelixID")

# Run association with adolescent anthropometry

# definition of outcomes
cont_outc_z<-c("hs2_zbmi")
cont_outc_r<- c("hs2_waist", "hs2_FatMass")
cat_outc<-c("hs2_ovob")

adolescent_nw<-adolescent %>% dplyr::filter(region=="n/w")
northern<-association_adolescence("Cluster", adolescent_nw)
write_csv(northern, "./results/csvs/adolescent_northern.csv")

adolescent_sm<-adolescent %>% dplyr::filter(region=="s/m")
southern<-association_adolescence("Cluster", adolescent_sm)
write_csv(southern, "./results/csvs/adolescent_southern.csv")

pooled<-association_adolescence("Cluster", adolescent)
write_csv(pooled, "./results/csvs/adolescent_pooled.csv")
