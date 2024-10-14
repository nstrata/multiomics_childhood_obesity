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
here::i_am("script/a.SNF clusters/analysis_trainset.R")
setwd(here::here())

# Load data 
df_omics<-read_rds( "./data/df_omics.rds")

#### Northern is the training dataset 
df_omics <- df_omics %>%
  dplyr::filter(region=="n/w") %>%
  dplyr::select(starts_with("hsa_"),starts_with("pro_"),
                starts_with("metser_"), starts_with("TC"),
                starts_with("cg"))
omics<-names(df_omics) # Vector of omics names

# We need ad input a list of matrices with features as rows and IDs in columns
mirna<-df_omics %>%
  dplyr::select(starts_with("hsa_")) %>% as.data.frame()

proteins<-df_omics %>%
  dplyr::select(starts_with("pro_")) %>% as.data.frame()

metabolites<-df_omics %>%
  dplyr::select(starts_with("metser_")) %>% as.data.frame()

cpgs<-df_omics %>%
  dplyr::select(starts_with("cg")) %>% as.data.frame()

tcs<-df_omics %>%
  dplyr::select(starts_with("TC"))%>% as.data.frame()

# Create similarity matrices
pearson_mirna=amap::Dist(mirna, "pearson") %>% as.matrix()
eucli_metabolome= amap::Dist(metabolites, "euclidean") %>% as.matrix()
eucli_proteome= amap::Dist(proteins, "euclidean") %>% as.matrix()
pearson_cpg= amap::Dist(cpgs, "pearson") %>% as.matrix()
pearson_tc= amap::Dist(tcs, "pearson") %>% as.matrix()

####################################################
# Check which hyperparameteres K and sighma to use 
####################################################
# credits to https://github.com/CDNMBioinformatics/LEOCC_metabolomics_SNF/blob/main/2_SNF_parameter_selection.Rmd

filtmetsEuclidDistList <- list(pearson_mirna,eucli_metabolome,eucli_proteome,pearson_cpg,pearson_tc)

for(i in 1:5){
  varResultsDf <- data.frame(KNN=numeric(),alpha=numeric(),variance=numeric())
  for(K in seq(50,90,5)){
    for(alpha in seq(0.3,0.8,0.1)){
      ## Draw similarity graphs
      simGraph <- affinityMatrix(filtmetsEuclidDistList[[i]], K, alpha)
      ## Calculate the variance of the upper triangle
      upperTriVar <- var(simGraph[upper.tri(simGraph)])
      ## Build results data frame
      varResultsDf <- rbind(varResultsDf,
                            data.frame(KNN=K,alpha=alpha,
                                       variance=upperTriVar))
    }
  }
}

###################################################################
# Based on variance use parameters of K=80 and sigma=0.8 for the SNF
###################################################################

#Construct affinity matrices
set.seed(100) 
K=50
sigma=0.8
T=20

W1 = affinityMatrix(pearson_mirna, K, sigma)
W2 = affinityMatrix(eucli_metabolome, K, sigma)
W3 = affinityMatrix(eucli_proteome, K, sigma)
W4 = affinityMatrix(pearson_cpg, K, sigma)
W5 = affinityMatrix(pearson_tc, K, sigma)
W = SNF(list(W1,W2,W3,W4,W5), K, T) # Similarity network fusion


# Calculate optimal number of clusters 
# credits to https://github.com/cran/SNFtool

# Obtain eigenvalues
optk_std<-estimate_k(W, maxk = 7, showplots = TRUE)
optk_std # 3 clusters

# Obtain Rotation cost values
.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}

.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}

rotation_cost_values <- function(W, NUMC) {
  
  W = (W + t(W))/2
  diag(W) = 0
  
  degs = rowSums(W)
  
  # compute unnormalized Laplacian
  degs[degs == 0] = .Machine$double.eps    
  D = diag(degs)    
  L = D - W
  Di = diag(1 / sqrt(degs))
  L = Di %*% L %*% Di
  
  # compute the eigenvectors corresponding to the k smallest
  # eigs$valuess
  eigs = eigen(L)
  eigs_order = sort(eigs$values, index.return=T)$ix
  eigs$values = eigs$values[eigs_order]
  eigs$vectors = eigs$vectors[, eigs_order]
  eigengap = abs(diff(eigs$values))
  eigengap = eigengap * (1 - eigs$values[1:length(eigs$values) - 1] ) / (1 - eigs$values[2:length(eigs$values)])
  
  quality = list()
  for (c_index in 1:length(NUMC)) {
    ck = NUMC[c_index]
    UU = eigs$vectors[, 1:ck]
    EigenvectorsDiscrete <- .discretisation(UU)[[1]]
    EigenVectors = EigenvectorsDiscrete^2
    
    # MATLAB: sort(EigenVectors,2, 'descend');
    temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), function(i) EigenVectors[, i])), ]
    temp1 <- t(apply(temp1, 1, sort, TRUE))  
    
    quality[[c_index]] = (1 - eigs$values[ck + 1]) / (1 - eigs$values[ck]) * 
      sum( sum( diag(1 / (temp1[, 1] + .Machine$double.eps) ) %*% temp1[, 1:max(2, ck-1)] ))
  }
  
  quality= unlist(quality) %>% as.data.frame() %>%
    dplyr::rename("Rotation cost values" = ".")
  
  
  return (quality)
}

rotation<-rotation_cost_values(W, 2:8)
rotation
rotation<-rbind(NA, rotation)

# Combined table of eigenvalues and rotation costs values
stats<-cbind(optk_std[1:8,1:2],rotation)
stats
