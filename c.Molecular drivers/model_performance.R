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
library(pROC)
library(caret)
library(xgboost)
library(metrica)

# Set working directory
here::i_am("script/c.Molecular drivers/model_performance.R")
setwd(here::here())


# Load datasets
fusion_all<-readRDS("./results/fusion_all.rds")
df_omics<-read_rds( "./data/df_omics.rds")
omics<-names(df_omics)[3:1381]

snf<-fusion_all[[paste0("id_northern_southern")]]
snf<-snf %>% dplyr::mutate(mocluster = ifelse(newgroups==3, "Cluster_C",
                                              ifelse(newgroups==2, "Cluster_A", 
                                                     ifelse(newgroups==1, "Cluster_B",NA)))) %>%
  dplyr::mutate(mocluster=as.factor(mocluster)) %>%
  dplyr::select(-c("cohort","newgroups")) %>%
  dplyr::mutate(moclusterC_AB= ifelse(mocluster=="Cluster_C", "Cluster_C",
                                      ifelse(mocluster=="Cluster_A" | mocluster=="Cluster_B", "Cluster_AB",NA)),
                moclusterB_AC= ifelse(mocluster=="Cluster_B", "Cluster_B",
                                      ifelse(mocluster=="Cluster_A" | mocluster=="Cluster_C", "Cluster_AC",NA))) %>%
  dplyr::mutate(moclusterC_AB=as.factor(moclusterC_AB),
                moclusterB_AC=as.factor(moclusterB_AC))

df_all<- df_omics %>% dplyr::inner_join(snf, by="HelixID")

# train dataset
train<-df_all%>%
  dplyr::filter(region=="n/w") %>%
  dplyr::select(mocluster, moclusterC_AB, moclusterB_AC, all_of(omics))

#test dataset
test<-df_all%>%
  dplyr::filter(region=="s/m") %>%
  dplyr::select(mocluster, moclusterC_AB, moclusterB_AC, all_of(omics))

########################################################################
#### Cluster C vs AB
########################################################################
# Print grdsearch from xgboost model to find the  paramaters
xgbmod <- readRDS("./results/xgbmod_CvsAB.rds")
xgbmod$finalModel
print(xgbmod)
# find best tune based on normalized gini index
xgbmod$bestTune

# Tune grid according to best tune parameters and then run the model in the train dataset to obtain the final model
tuneGridXGB_full <- expand.grid(
  nrounds=260,
  max_depth = 6,
  eta = 0.1,
  gamma = 0.01,
  colsample_bytree = 0.6,
  subsample = 0.5,
  min_child_weight = 0)

#####################
# final datasets
#####################

# Train dataset
train_CvsAB<- train %>%
  dplyr::select(all_of(omics), moclusterC_AB) %>%
  dplyr::mutate(moclusterC_AB=ifelse(moclusterC_AB=="Cluster_AB",0,
                                     ifelse(moclusterC_AB=="Cluster_C",1,NA)))
train_data_C_AB<-as.matrix(train_CvsAB[,1:1379])
train_label_C_AB=train_CvsAB$moclusterC_AB

# Test dataset
test_CvsAB<- test %>%
  dplyr::select(all_of(omics), moclusterC_AB) %>%
  dplyr::mutate(moclusterC_AB=ifelse(moclusterC_AB=="Cluster_AB",0,
                                     ifelse(moclusterC_AB=="Cluster_C",1,NA)))
test_data_C_AB<-as.matrix(test_CvsAB[,1:1379])
test_label_C_AB=test_CvsAB$moclusterC_AB

# Run xgboost model
set.seed(2136)
xgbmod_full <- xgboost(data=train_data_C_AB,
                     label=train_label_C_AB,
                     nrounds=260,
                     max_depth = 6,
                     eta = 0.1,
                     gamma = 0.01,
                     colsample_bytree = 0.6,
                     subsample = 0.5,
                     min_child_weight = 0,
                     objective = "binary:logistic")

xgbclasses <- predict(xgbmod_full, newdata = test_data_C_AB)
xgbclasses <- as.factor(ifelse (xgbclasses > 0.5,1,0))
xgbprobs = predict(xgbmod_full, newdata = test_data_C_AB, type="prob")

test_label_C_AB<-as.factor(test_label_C_AB)
confusionMatrix(data = xgbclasses, test_label_C_AB, positive="1")
metrics_summary(pred=xgbclasses, obs=test_label_C_AB, pos_level = 2,type = "classification")
roc_C <- roc(test_label_C_AB, xgbprobs)
roc_C

#####################################################################################
#### Cluster B vs AC
####################################################################################
# Print grdsearch from xgboost model to find the  paramaters
xgbmod <- readRDS("./results/xgbmod_BvsAC.rds")
xgbmod$finalModel
print(xgbmod)
# find best tune based on normalized gini index
xgbmod$bestTune

# Tune grid according to best tune parameters and then run the model in the train dataset to obtain the final model
tuneGridXGB_full <- expand.grid(
  nrounds=340,
  max_depth = 5,
  eta = 0.05,
  gamma = 0.5,
  colsample_bytree = 0.7,
  subsample = 0.5,
  min_child_weight = 1)

# Train dataset
train_BvsAC<- train %>%
  dplyr::select(all_of(omics), moclusterB_AC) %>%
  dplyr::mutate(moclusterB_AC=ifelse(moclusterB_AC=="Cluster_AC",0,
                                     ifelse(moclusterB_AC=="Cluster_B",1,NA)))
train_data_B_AC<-as.matrix(train_BvsAC[,1:1379])
train_label_B_AC=train_BvsAC$moclusterB_AC

# Test dataset
test_BvsAC<- test %>%
  dplyr::select(all_of(omics), moclusterB_AC) %>%
  dplyr::mutate(moclusterB_AC=ifelse(moclusterB_AC=="Cluster_AC",0,
                                     ifelse(moclusterB_AC=="Cluster_B",1,NA)))
test_data_B_AC<-as.matrix(test_BvsAC[,1:1379])
test_label_B_AC=test_BvsAC$moclusterB_AC

# Run xgboost model
xgbmod_full <- xgboost(data=train_data_B_AC,
                       label=train_label_B_AC,
                       nrounds=340,
                       max_depth = 5,
                       eta = 0.05,
                       gamma = 0.5,
                       colsample_bytree = 0.7,
                       subsample = 0.5,
                       min_child_weight = 1,
                       objective = "binary:logistic")


xgbclasses <- predict(xgbmod_full, newdata = test_data_B_AC)
xgbclasses <- as.factor(ifelse (xgbclasses > 0.5,1,0))
xgbprobs = predict(xgbmod_full, newdata = test_data_B_AC, type="prob")

test_label_B_AC<-as.factor(test_label_B_AC)
confusionMatrix(data = xgbclasses, test_label_B_AC, positive="1")
metrics_summary(pred=xgbclasses, obs=test_label_B_AC, pos_level = 2,type = "classification")
roc_C <- roc(test_label_B_AC, xgbprobs)
roc_C




