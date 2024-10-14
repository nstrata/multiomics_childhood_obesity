rm(list=ls())
gc()

# Call libraries needed
library(here)
library(dplyr)
library(tidyverse)
library(foreign)
library(SNFtool)
library(Spectrum)
library(amap)

# Set working directory
here::i_am("script/a.SNF clusters/propagation.R")
setwd(here::here())

# Load data
df<-read_rds( "./data/df_omics.rds")

# Code to run the label classifier propagation procedure to the Southern European group of children
# FROM SNFtool package https://github.com/cran/SNFtool
# This function implements the label propagation to predict the label(subtype) for new patients.	
# method == 0 indicates to use the local and global consistency method
# method >0 indicates to use label propagation method.
.csPrediction <- function(W,Y0,method){
  alpha=0.9
  P= W/rowSums(W)
  if(method==0){
    Y= (1-alpha)* solve( diag(dim(P)[1])- alpha*P)%*%Y0;
  } else {
    NLabel=which(rowSums(Y0)==0)[1]-1;
    Y=Y0;
    for (i in 1:1000){
      Y=P%*%Y;
      Y[1:NLabel,]=Y0[1:NLabel,];
    }
  }
  return(Y);
}

######## Northern/Western vs Southern Europe
fusion_all=list()
# First run the SNF in the train set
df_train<-df %>% 
  dplyr::filter(region=="n/w")
df_train_id <- df_train %>% 
  dplyr::select(HelixID, cohort)

# It as input a list of matrices with features as rows and IDs in columns
mirna_train<-df_train %>% dplyr::select(starts_with("hsa_")) %>% as.data.frame()
proteins_train<-df_train %>% dplyr::select(starts_with("pro_")) %>% as.data.frame()
metabolites_train<-df_train %>% dplyr::select(starts_with("metser_")) %>% as.data.frame()
cpgs_train<-df_train %>% dplyr::select(starts_with("cg")) %>% as.data.frame()
tcs_train<-df_train %>% dplyr::select(starts_with("TC"))%>% as.data.frame()

# Create affinity matrices
pearson_mirna_train=amap::Dist(mirna_train, "pearson") %>% as.matrix()
eucli_metabolome_train= amap::Dist(metabolites_train, "euclidean") %>% as.matrix()
eucli_proteome_train= amap::Dist(proteins_train, "euclidean") %>% as.matrix()
pearson_cpg_train= amap::Dist(cpgs_train, "pearson") %>% as.matrix()
pearson_tc_train= amap::Dist(tcs_train, "pearson") %>% as.matrix()

set.seed(100)
K=50
sigma=0.8
T=20

W1_train = affinityMatrix(pearson_mirna_train, K, sigma)
W2_train = affinityMatrix(eucli_metabolome_train, K, sigma)
W3_train = affinityMatrix(eucli_proteome_train, K, sigma)
W4_train = affinityMatrix(pearson_cpg_train, K, sigma)
W5_train = affinityMatrix(pearson_tc_train, K, sigma)
W_train= SNF(list(W1_train,W2_train,W3_train,W4_train,W5_train), K, T)
C=3
Cluster_train = spectralClustering(W_train,C)
table(Cluster_train)

fusion_all[[paste0("W_northern")]] = W_train # Similarity network fusion
fusion_all[[paste0("Cluster_northern")]] = as.numeric(Cluster_train)
fusion_all[[paste0("id_northern")]] = cbind(df_train_id, Cluster_train)

# Use the classifier propagation approach

# first subset for the test dataset / the cohort exluded
df_test<-df %>% dplyr::filter(region=="s/m")
df_test_id<-df_test %>% dplyr::select(HelixID, cohort)

mirna_test<-df_test %>% dplyr::select(starts_with("hsa_")) %>% as.data.frame()
proteins_test<-df_test %>% dplyr::select(starts_with("pro_")) %>% as.data.frame()
metabolites_test<-df_test %>% dplyr::select(starts_with("metser_")) %>% as.data.frame()
cpgs_test<-df_test %>% dplyr::select(starts_with("cg")) %>% as.data.frame()
tcs_test<-df_test %>% dplyr::select(starts_with("TC"))%>% as.data.frame()

# Create combined dataset of test and train
mirna_train_test = rbind(mirna_train,mirna_test)
proteins_train_test = rbind(proteins_train,proteins_test)
metabolites_train_test = rbind(metabolites_train,metabolites_test)
cpgs_train_test = rbind(cpgs_train,cpgs_test)
tcs_train_test = rbind(tcs_train,tcs_test)

df_train_test_id<-rbind(df_train_id, df_test_id)

# Calculate distance matrix in combined datasets
pearson_mirna_train_test=amap::Dist(mirna_train_test, "pearson") %>% as.matrix()
eucli_metabolome_train_test= amap::Dist(metabolites_train_test, "euclidean") %>% as.matrix()
eucli_proteome_train_test= amap::Dist(proteins_train_test, "euclidean") %>% as.matrix()
pearson_cpg_train_test= amap::Dist(cpgs_train_test, "pearson") %>% as.matrix()
pearson_tc_train_test= amap::Dist(tcs_train_test, "pearson") %>% as.matrix()

K=50
sigma=0.8
T=20

W1_train_test = affinityMatrix(pearson_mirna_train_test, K, sigma)
W2_train_test = affinityMatrix(eucli_metabolome_train_test, K, sigma)
W3_train_test = affinityMatrix(eucli_proteome_train_test, K, sigma)
W4_train_test = affinityMatrix(pearson_cpg_train_test, K, sigma)
W5_train_test = affinityMatrix(pearson_tc_train_test, K, sigma)
W_train_test= SNF(list(W1_train_test,W2_train_test,W3_train_test,W4_train_test,W5_train_test), K, T)

groups<- as.numeric(as.character(Cluster_train))
table(groups)

Y0=matrix(0,nrow(df_train_test_id), max(groups))
for (i in 1:length(groups)) {
  Y0[i,groups[i]]=1
}
Y=.csPrediction(W_train_test,Y0,1)
newgroups=rep(0,nrow(df_train_test_id))
for (i in 1:nrow(Y)) {
  newgroups[i]=which(Y[i,]==max(Y[i,]))
}

fusion_all[[paste0("W_northern_southern")]]=W_train_test
fusion_all[[paste0("Clusters_northern_southern")]]=newgroups
fusion_all[[paste0("id_northern_southern")]]=cbind(df_train_test_id, newgroups)

saveRDS(fusion_all, "./results/fusion_all.rds")