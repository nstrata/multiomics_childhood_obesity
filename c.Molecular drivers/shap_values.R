rm(list=ls())
gc()

#################################################
# Call libraries needed
#################################################
library(here)
library(dplyr)
library(tidyverse)
library(foreign)
library(data.table)
library(caret)
library(xgboost)

# Set working directory
here::i_am("script/c.Molecular drivers/shap_values.R")
setwd(here::here())

# Load and merge datasets
fusion_all<-readRDS("./results/fusion_all.rds")
df_omics<-read_rds( "./data/df_omics.rds")
omics<-names(df_omics)[3:1381]

snf<-fusion_all[[paste0("id_northern_southern")]]
snf<-snf %>%
  dplyr::mutate(mocluster = ifelse(newgroups==3, "Cluster_C",
                                   ifelse(newgroups==2, "Cluster_A", 
                                          ifelse(newgroups==1, "Cluster_B",NA)))) %>%
  dplyr::mutate(mocluster=as.factor(mocluster),
                cohort_check = as.factor(cohort)) %>%
  dplyr::select(-c("cohort","newgroups")) %>%
  dplyr::mutate(moclusterC_AB= ifelse(mocluster=="Cluster_C", "Cluster_C",
                                      ifelse(mocluster=="Cluster_A" | mocluster=="Cluster_B", "Cluster_AB",NA)),
                moclusterB_AC= ifelse(mocluster=="Cluster_B", "Cluster_B",
                                      ifelse(mocluster=="Cluster_A" | mocluster=="Cluster_C", "Cluster_AC",NA))) %>%
  dplyr::mutate(moclusterC_AB=as.factor(moclusterC_AB),
                moclusterB_AC=as.factor(moclusterB_AC))


df_all<- df_omics %>% dplyr::right_join(snf, by="HelixID")
rownames(df_all) <-df_all$HelixID

train<-df_all%>%
  dplyr::filter(region=="n/w") %>%
  dplyr::select(mocluster, moclusterC_AB, moclusterB_AC, all_of(omics)) 

test<-df_all%>%
  dplyr::filter(region=="s/m") %>%
  dplyr::select(mocluster, moclusterC_AB, moclusterB_AC, all_of(omics))

#####################
# Cluster C vs AB
####################

# Train
train_CvsAB<- train %>%
  dplyr::select(all_of(omics), moclusterC_AB) %>%
  dplyr::mutate(moclusterC_AB=ifelse(moclusterC_AB=="Cluster_AB",0,
                                     ifelse(moclusterC_AB=="Cluster_C",1,NA)))
train_data_C_AB<-as.matrix(train_CvsAB[,1:1379])
train_label_C_AB=train_CvsAB$moclusterC_AB

# Test
test_CvsAB<- test %>%
  dplyr::select(all_of(omics), moclusterC_AB) %>%
  dplyr::mutate(moclusterC_AB=ifelse(moclusterC_AB=="Cluster_AB",0,
                                     ifelse(moclusterC_AB=="Cluster_C",1,NA)))
test_data_C_AB<-as.matrix(test_CvsAB[,1:1379])
test_label_C_AB=test_CvsAB$moclusterC_AB

# model in Train data
set.seed(2136)
xgbmod_train <- xgboost(data=train_data_C_AB,
                       label=train_label_C_AB,
                       nrounds=260,
                       max_depth = 6,
                       eta = 0.1,
                       gamma = 0.01,
                       colsample_bytree = 0.6,
                       subsample = 0.5,
                       min_child_weight = 0,
                       objective = "binary:logistic")

#obtain shap values
shap_value_train_C_AB<- predict(xgbmod_train, train_data_C_AB, predcontrib=T) %>% as.data.frame() %>%
  dplyr::select(-c(BIAS))
shap_value_train_C_AB<-cbind(train_label_C_AB, shap_value_train_C_AB) 
write_csv(shap_value_train_C_AB, "./results/csvs/shap_value_train_C_AB.csv")


# model in Test data
set.seed(2136)
xgbmod_test <- xgboost(data=test_data_C_AB,
                       label=test_label_C_AB,
                       nrounds=260,
                       max_depth = 6,
                       eta = 0.1,
                       gamma = 0.01,
                       colsample_bytree = 0.6,
                       subsample = 0.5,
                       min_child_weight = 0,
                       objective = "binary:logistic")

#obtain shap values
shap_value_test_C_AB<- predict(xgbmod_test, test_data_C_AB, predcontrib=T) %>% as.data.frame() %>%
  dplyr::select(-c(BIAS))
shap_value_test_C_AB<-cbind(test_label_C_AB, shap_value_test_C_AB) 
write_csv(shap_value_test_C_AB, "./results/csvs/shap_value_test_C_AB.csv")

#####################
# Cluster B vs AC
####################

# Train dataset
train_BvsAC<- train %>%
  dplyr::select(all_of(omics), moclusterB_AC) %>%
  dplyr::mutate(moclusterB_AC=ifelse(moclusterB_AC=="Cluster_AC",0,
                                     ifelse(moclusterB_AC=="Cluster_B",1,NA)))
train_data_B_AC<-as.matrix(train_BvsAC[,1:1379])
train_label_B_AC=train_BvsAC$moclusterB_AC
table(train_label_B_AC)

# Test dataset
test_BvsAC<- test %>%
  dplyr::select(all_of(omics), moclusterB_AC) %>%
  dplyr::mutate(moclusterB_AC=ifelse(moclusterB_AC=="Cluster_AC",0,
                                     ifelse(moclusterB_AC=="Cluster_B",1,NA)))
test_data_B_AC<-as.matrix(test_BvsAC[,1:1379])
test_label_B_AC=test_BvsAC$moclusterB_AC
table(test_label_B_AC)

# model in Train data
set.seed(2136)
xgbmod_train <- xgboost(data=train_data_B_AC,
                       label=train_label_B_AC,
                       nrounds=340,
                       max_depth = 5,
                       eta = 0.05,
                       gamma = 0.5,
                       colsample_bytree = 0.7,
                       subsample = 0.5,
                       min_child_weight = 1,
                       objective = "binary:logistic")

#obtain shap values
shap_value_train_B_AC<- predict(xgbmod_train, train_data_B_AC, predcontrib=T) %>% as.data.frame() %>%
  dplyr::select(-c(BIAS))
shap_value_train_B_AC<-cbind(train_label_B_AC, shap_value_train_B_AC) 
write_csv(shap_value_train_B_AC, "./results/csvs/shap_value_train_B_AC.csv")

# model in Test data
set.seed(2136)
xgbmod_test <- xgboost(data=test_data_B_AC,
                       label=test_label_B_AC,
                       nrounds=340,
                       max_depth = 5,
                       eta = 0.05,
                       gamma = 0.5,
                       colsample_bytree = 0.7,
                       subsample = 0.5,
                       min_child_weight = 1,
                       objective = "binary:logistic")

#obtain shap values
shap_value_test_B_AC<- predict(xgbmod_test, test_data_B_AC, predcontrib=T) %>% as.data.frame() %>%
  dplyr::select(-c(BIAS))
shap_value_test_B_AC<-cbind(test_label_B_AC, shap_value_test_B_AC) 
write_csv(shap_value_test_B_AC, "./results/csvs/shap_value_test_B_AC.csv")


# Calculate mean shap values

# Train data
# CvsAB
shap_train_C_AB<-cbind(colMeans(abs(shap_value_train_C_AB[,2:1380])), colMeans(shap_value_train_C_AB[,2:1380])) %>% as.data.frame %>%
  dplyr::rename(abs_mean_shap_train_C_AB = V1,
                mean_shap_train_C_AB = V2) %>%
  dplyr::mutate(direction_C_AB=ifelse(mean_shap_train_C_AB<0, "negative",
                                      ifelse(mean_shap_train_C_AB>0, "positive", NA)),
                abs_train_C_AB_quar = ntile(abs_mean_shap_train_C_AB,4))
shap_train_C_AB$features=rownames(shap_train_C_AB)

#BvsAC
shap_train_B_AC<-cbind(colMeans(abs(shap_value_train_B_AC[,2:1380])), colMeans(shap_value_train_B_AC[,2:1380])) %>% as.data.frame %>%
  dplyr::rename(abs_mean_shap_train_B_AC = V1,
                mean_shap_train_B_AC = V2) %>%
  dplyr::mutate(direction_B_AC=ifelse(mean_shap_train_B_AC<0, "negative",
                                      ifelse(mean_shap_train_B_AC>0, "positive", NA)),
                abs_train_B_AC_quar = ntile(abs_mean_shap_train_B_AC,4))
shap_train_B_AC$features=rownames(shap_train_B_AC)

# Test data
# CvsAB
shap_test_C_AB<-cbind(colMeans(abs(shap_value_test_C_AB[,2:1380])), colMeans(shap_value_test_C_AB[,2:1380])) %>% as.data.frame %>%
  dplyr::rename(abs_mean_shap_test_C_AB = V1,
                mean_shap_test_C_AB = V2) %>%
  dplyr::mutate(direction_C_AB=ifelse(mean_shap_test_C_AB<0, "negative",
                                      ifelse(mean_shap_test_C_AB>0, "positive", NA)),
                abs_test_C_AB_quar = ntile(abs_mean_shap_test_C_AB,4))
shap_test_C_AB$features=rownames(shap_test_C_AB)

#BvsAC
shap_test_B_AC<-cbind(colMeans(abs(shap_value_test_B_AC[,2:1380])), colMeans(shap_value_test_B_AC[,2:1380])) %>% as.data.frame %>%
  dplyr::rename(abs_mean_shap_test_B_AC = V1,
                mean_shap_test_B_AC = V2) %>%
  dplyr::mutate(direction_B_AC=ifelse(mean_shap_test_B_AC<0, "negative",
                                      ifelse(mean_shap_test_B_AC>0, "positive", NA)),
                abs_test_B_AC_quar = ntile(abs_mean_shap_test_B_AC,4))
shap_test_B_AC$features=rownames(shap_test_B_AC)


# merge and then select the ones with absolute SHAP values at the top quartile and mean shap values at the same direction
# CvsAB
top_shap_all_C_AB<-merge(shap_train_C_AB, shap_test_C_AB, by=c("features","direction_C_AB")) 
names(top_shap_all_C_AB)
top_shap_C_AB<-top_shap_all_C_AB %>% drop_na() %>%
  dplyr::filter(abs_train_C_AB_quar==4 & abs_test_C_AB_quar==4) %>% dplyr::arrange(desc(features))
top_C_AB<-top_shap_C_AB$features
top_C_AB
write_csv(top_shap_C_AB, "./results/csvs/top_shap_C_AB.csv")

# BvsAC
top_shap_all_B_AC<-merge(shap_train_B_AC, shap_test_B_AC, by=c("features","direction_B_AC")) 
names(top_shap_all_B_AC)
top_shap_B_AC<-top_shap_all_B_AC %>% drop_na() %>%
  dplyr::filter(abs_test_B_AC_quar==4 & abs_test_B_AC_quar==4) %>% dplyr::arrange(desc(features))
top_B_AC<-top_shap_B_AC$features
top_B_AC
write_csv(top_shap_B_AC, "./results/csvs/top_shap_B_AC.csv")

