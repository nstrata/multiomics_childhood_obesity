rm(list=ls())
gc()

# Call libraries needed
library(here)
library(dplyr)
library(tidyverse)
library(foreign)
library(data.table)
library(fastDummies)
library(nnet)
library(ggeffects)

# Set working directory
here::i_am("script/d. Prenatal determinants/probabilities.R")
setwd(here::here())

# load data
exposures_preg<-readRDS("./data/selected_preg.rds")
fusion_all<-readRDS("./results/fusion_all.rds")

snf<-fusion_all[[paste0("id_northern_southern")]]
snf_pooled<-snf %>%
  dplyr::mutate(Cluster = ifelse(newgroups==3, "Cluster C",
                                        ifelse(newgroups==2, "Cluster A", 
                                               ifelse(newgroups==1, "Cluster B",NA)))) %>%
  dplyr::mutate(Cluster=as.factor(Cluster)) %>%
  dplyr::select(-c("cohort","newgroups"))

# Dataset containing pregnancy exposome and cluster
df<-merge(snf_pooled, exposures_preg, by="HelixID")

# Northern/Western cohort
df_nw<-df %>% dplyr::filter(region=="n/w")
df_nw <- dummy_cols(df_nw, select_columns = "h_cohort") %>%
  dplyr::select(-c("h_cohort_BIB"))

# mbmi
mbmi<-multinom(Cluster ~ h_mbmi_None+h_sex+hs_age_years+
                 h_cohort_EDEN+h_cohort_KANC+h_cohort_MOBA, data=df_nw, model = T, Hess=T)
predict(mbmi, df_nw, type="probs")
predictions.mbmi <- ggeffects::ggemmeans(mbmi, terms="h_mbmi_None [all]") %>%
  dplyr::rename(Cluster=response.level)

# pfoa
pfoa<-multinom(Cluster ~ hs_pfoa_m_Log2+h_sex+hs_age_years+
                 h_cohort_EDEN+h_cohort_KANC+h_cohort_MOBA, data=df_nw, model = T, Hess=T)
predict(pfoa, df_nw, type="probs")
predictions.pfoa <- ggeffects::ggemmeans(pfoa, terms="hs_pfoa_m_Log2 [all]") %>%
  dplyr::rename(Cluster=response.level)

# Sourthern/Mediterranean cohort
df_sm<-df%>% dplyr::filter(region=="s/m")
df_sm <- dummy_cols(df_sm, select_columns = "h_cohort") %>%
  dplyr::select(-c("h_cohort_RHEA"))

# hg
hg<-multinom(Cluster ~ hs_hg_m_Log2+h_sex+hs_age_years+
               h_cohort_INMA, data=df_sm, model = T, Hess=T)
predict(hg, df_sm, type="probs")
predictions.hg <- ggeffects::ggemmeans(hg, terms="hs_hg_m_Log2 [all]") %>%
  dplyr::rename(Cluster=response.level)
