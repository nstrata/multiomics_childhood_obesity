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


# Set working directory
here::i_am("script/c.Molecular drivers/caret_gridsearch.R")
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

###############################################
# Fine tune the model
###############################################
# normalized gini function
normalizedGini <- function(aa, pp) {
  Gini <- function(a, p) {
    if (length(a) !=  length(p)) stop("Actual and Predicted need to be equal lengths!")
    temp.df <- data.frame(actual = a, pred = p, range=c(1:length(a)))
    temp.df <- temp.df[order(-temp.df$pred, temp.df$range),]
    population.delta <- 1 / length(a)
    total.losses <- sum(a)
    null.losses <- rep(population.delta, length(a)) # Hopefully is similar to accumulatedPopulationPercentageSum
    accum.losses <- temp.df$actual / total.losses # Hopefully is similar to accumulatedLossPercentageSum
    gini.sum <- cumsum(accum.losses - null.losses) # Not sure if this is having the same effect or not
    sum(gini.sum) / length(a)
  }
  Gini(aa,pp) / Gini(aa,aa)
}

# create the normalized gini summary function to pass into caret
giniSummary <- function (data, lev = "Yes", model = NULL) {
  levels(data$obs) <- c('0', '1')
  out <- normalizedGini(as.numeric(levels(data$obs))[data$obs], data[, lev[2]])  
  names(out) <- "NormalizedGini"
  out
}

# BvsAC #####################################################
fit.control = trainControl(
  method = 'cv',
  number = 5,
  summaryFunction = giniSummary,
  classProbs = T,
  savePredictions = T,
  verboseIter = TRUE,
  allowParallel = TRUE)

# create the tuning grid
# https://www.kaggle.com/code/pelkoja/visual-xgboost-tuning-with-caret
tuneGridXGB <- expand.grid(
  nrounds=(1:20)*20,
  max_depth = c(3, 4, 5, 6),
  eta = c(0.025, 0.05, 0.1, 0.5),
  gamma = c(0, 0.01, 0.05, 0.1, 0.5, 0.7, 0.9),
  colsample_bytree = c(0.4,0.6,0.7,0.8),
  subsample = c(0.50, 0.75, 1),
  min_child_weight = c(0,1))

xgbmod <- train(train[,4:1382],train$moclusterB_AC,
                method = 'xgbTree',
                metric = 'NormalizedGini',
                trControl = fit.control,
                tuneGrid = tuneGridXGB,
                verbosity = 0)
write_rds(xgbmod, "./results/xgbmod_BvsAC.rds")

#### C vs AB ##############################################################
fit.control = trainControl(
  method = 'cv',
  number = 5,
  summaryFunction = giniSummary,
  classProbs = T,
  savePredictions = T,
  verboseIter = TRUE,
  allowParallel = TRUE)

# create the tuning grid. Again keeping this small to avoid exceeding kernel memory limits.
# You can expand as your compute resources allow. 
# https://www.kaggle.com/code/pelkoja/visual-xgboost-tuning-with-caret
tuneGridXGB <- expand.grid(
  nrounds=(1:20)*20,
  max_depth = c(3, 4, 5, 6),
  eta = c(0.025, 0.05, 0.1, 0.5),
  gamma = c(0, 0.01, 0.05, 0.1, 0.5, 0.7, 0.9),
  colsample_bytree = c(0.4,0.6,0.7,0.8),
  subsample = c(0.50, 0.75, 1),
  min_child_weight = c(0,1))

# CvsAB #####################################################
# train the xgboost learner
xgbmod <- train(train[,4:1382],train$moclusterC_AB,
                method = 'xgbTree',
                metric = 'NormalizedGini',
                trControl = fit.control,
                tuneGrid = tuneGridXGB,
                verbosity = 0)
write_rds(xgbmod, "./results/xgbmod_CvsAB.rds")

