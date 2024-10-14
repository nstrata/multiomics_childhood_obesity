rm(list=ls())
gc()

# Call libraries needed
library(here)
library(dplyr)
library(tidyverse)
library(foreign)
library(data.table)
library(fastDummies)
library(sharp)

# Set working directory
here::i_am("script/d.Prenatal determinants/lasso_sharp.R")
setwd(here::here())

# load data and merge datasets
exposures_preg<-readRDS("./data/exposures_preg.rds")

fusion_all<-readRDS("./results/fusion_all.rds")
snf<-fusion_all[[paste0("id_northern_southern")]]

snf_pooled<-snf %>%
  dplyr::mutate(Cluster = ifelse(newgroups==3, "Cluster C",
                                 ifelse(newgroups==2, "Cluster A", 
                                        ifelse(newgroups==1, "Cluster B",NA)))) %>%
  dplyr::mutate(Cluster=as.factor(Cluster)) %>%
  dplyr::select(-c("cohort","newgroups"))

df<-merge(snf_pooled, exposures_preg, by="HelixID") 

# Northern/Western cohort
df_nw <- df %>% dplyr::filter(region=="n/w")
df_nw <- dummy_cols(df_nw, select_columns = "h_cohort") %>%
  relocate(c("h_cohort_2","h_cohort_4","h_cohort_5"), .before = "h_sex") %>%
  dplyr::select(-c("h_cohort_1", "h_cohort"))

# Generate outcome dataset
pheno<-data.frame(df_nw[,c(".imp", "Cluster")])
pheno$Cluster<- relevel(pheno$Cluster, ref = "Cluster A")

# Generate a dataframe containing the prenatal variables
exposome<-df_nw %>% dplyr::select(-c("HelixID", "Cluster", ".id", "region")) 

# Run lasso analysis
results_stab <- list()
outcome_columns <- c("Cluster")
for (i in 1:20) {
  exposome_i<-exposome[exposome$.imp==i,] %>%
    dplyr::select(2:44)
  pheno_i<-pheno[pheno$.imp==i,] %>%
    dplyr::select(2)
  # set the penalty
  colcov<-5
  coln<-ncol(exposome_i)-5
  penalty <- c(rep(0, colcov), rep(1, coln))
  penalty
  length(penalty)
  print(paste0("Model: ",names(pheno_i)))
  print(system.time({
    results_stab[[i]] <- VariableSelection(
      xdata = exposome_i, ydata = pheno_i,
      family = "multinomial", penalty.factor = penalty,
      Lambda = LambdaSequence(lmax = 1e-1, lmin = 1e-3, cardinal = 500),
      resampling = "subsampling", tau=0.8, K = 100, n_cat=3 , pi_list=seq(0.51, 0.99, by = 0.01), beep=NULL, verbose=TRUE
    )
  }))
  
}
saveRDS(results_stab, "./results/lasso_northern_western.rds")

# Organizing selected variables according to selection proportion
pi_list=seq(0.51, 0.99, by = 0.01)  
combine_selected <- list()
threshold <- list()
for (i in 1:20){
  m.i<-cbind(names(sort(SelectionProportions(results_stab[[i]]), decreasing = TRUE)),
             sort(SelectionProportions(results_stab[[i]]), decreasing = TRUE))
  combine_selected[[i]]<-m.i[order(m.i[,1],decreasing=TRUE),]
  argmax_id <- ArgmaxId(results_stab[[i]])[2]
  threshold[[i]]<-pi_list[argmax_id]
}
SELECTPROPOR <- as.data.frame(do.call(cbind, combine_selected))
UMBRAL <-as.numeric(do.call(cbind, threshold))
UMBRALMEAN= mean(UMBRAL)
UMBRAL
row.names(SELECTPROPOR)<-SELECTPROPOR$V1
impares <- seq(1,ncol(SELECTPROPOR),by=2)
SELECTPROPOR = SELECTPROPOR[,-impares]
SELECTPROPOR[] <- lapply(SELECTPROPOR, function(x) if(is.character(x)) as.numeric(x) else x)
head(SELECTPROPOR)
SELECTPROPOR$ROW_PROM <- rowMeans(SELECTPROPOR)
SELECTPROPORFILTERED<- SELECTPROPOR[which(SELECTPROPOR$ROW_PROM >= UMBRALMEAN),]
SELECTPROPORFILTERED
selected<-cbind(rownames(SELECTPROPORFILTERED), SELECTPROPORFILTERED$ROW_PROM) %>% as.data.frame()
colnames(selected) <- c("Exposures", "Mean_Proportion")

# Southern/Mediterranean cohort
df_sm <- df %>% dplyr::filter(region=="s/m")
df_sm <- dummy_cols(df_sm, select_columns = "h_cohort") %>%
  relocate(c("h_cohort_6"), .before = "h_sex") %>%
  dplyr::select(-c("h_cohort_3", "h_cohort"))
  
# Generate outcome dataset
pheno<-data.frame(df_sm[,c(".imp", "Cluster")])
pheno$Cluster<- relevel(pheno$Cluster, ref = "Cluster A")

# Generate a dataframe containing the exposome variables
exposome<-df_sm %>% dplyr::select(-c("HelixID", "Cluster", ".id", "region")) 

# Run analysis
library(sharp)
results_stab <- list()
outcome_columns <- c("Cluster")
  for (i in 1:20) {
    exposome_i<-exposome[exposome$.imp==i,] %>%
      dplyr::select(2:42)
    pheno_i<-pheno[pheno$.imp==i,] %>%
      dplyr::select(2)
    # set the penalty
    colcov<-3
    coln<-ncol(exposome_i)-3
    penalty <- c(rep(0, colcov), rep(1, coln)) 
    penalty
    length(penalty)
    print(paste0("Model: ",names(pheno_i)))
    print(system.time({
    results_stab[[i]] <- VariableSelection(
      xdata = exposome_i, ydata = pheno_i,
      family = "multinomial", penalty.factor = penalty,
      Lambda = LambdaSequence(lmax = 1e-1, lmin = 1e-3, cardinal = 500),
      resampling = "subsampling", tau=0.8, K = 100, n_cat=3 , pi_list=seq(0.51, 0.99, by = 0.01), beep=NULL, verbose=TRUE
    )
    }))
  }
saveRDS(results_stab, "./results/lasso_southern_mediterranean.rds")

# Organizing selected variables according to selection proportion
  pi_list=seq(0.51, 0.99, by = 0.01)  
  combine_selected <- list()
  threshold <- list()
  for (i in 1:20){
    m.i<-cbind(names(sort(SelectionProportions(results_stab[[i]]), decreasing = TRUE)),
               sort(SelectionProportions(results_stab[[i]]), decreasing = TRUE))
    combine_selected[[i]]<-m.i[order(m.i[,1],decreasing=TRUE),]
    argmax_id <- ArgmaxId(results_stab[[i]])[2]
    threshold[[i]]<-pi_list[argmax_id]
  }
  SELECTPROPOR <- as.data.frame(do.call(cbind, combine_selected))
  UMBRAL <-as.numeric(do.call(cbind, threshold))
  UMBRALMEAN= mean(UMBRAL)
  UMBRAL
  row.names(SELECTPROPOR)<-SELECTPROPOR$V1
  impares <- seq(1,ncol(SELECTPROPOR),by=2)
  SELECTPROPOR = SELECTPROPOR[,-impares]
  SELECTPROPOR[] <- lapply(SELECTPROPOR, function(x) if(is.character(x)) as.numeric(x) else x)
  head(SELECTPROPOR)
  SELECTPROPOR$ROW_PROM <- rowMeans(SELECTPROPOR)
  SELECTPROPORFILTERED<- SELECTPROPOR[which(SELECTPROPOR$ROW_PROM >= UMBRALMEAN),]
  SELECTPROPORFILTERED
  selected<-cbind(rownames(SELECTPROPORFILTERED), SELECTPROPORFILTERED$ROW_PROM) %>% as.data.frame()
  colnames(selected) <- c("Exposures", "Mean_Proportion")


