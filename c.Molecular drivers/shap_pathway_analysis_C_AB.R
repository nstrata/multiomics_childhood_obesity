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
library(data.table)
library(ggpubr)
library(clusterProfiler)

# Set working directory
here::i_am("script/c.Molecular drivers/shap_pathway_analysis_C_AB.R")
setwd(here::here())

# Load data 
top_shap_C_AB<-read.csv("./results/csvs/top_shap_C_AB.csv")

# First identify the top features  per layer and obtain vectors of features
shap_top<-top_shap_C_AB %>% t() %>% as.data.frame() 
colnames(shap_top) <-shap_top[1,]

# Proteins
shap_pro <- shap_top %>%
  dplyr::select(starts_with("pro")) %>%
  rename_all(~stringr::str_replace_all(.,"pro_",""))
shap_pro<-shap_pro[ , order(names(shap_pro))] # order by column names
proteins<-colnames(shap_pro)

# metabolites
shap_metab <- shap_top %>%
  dplyr::select(starts_with("metser")) %>%
  rename_all(~stringr::str_replace_all(.,"metser_log10.","")) %>% 
  rename_all(~stringr::str_replace_all(.,"_",":"))
shap_metab<-shap_metab[ , order(names(shap_metab))] # order by column names
metabolites<-colnames(shap_metab)

# mirnas
shap_mirnas <- shap_top %>%
  dplyr::select(starts_with("hsa")) %>%
  rename_all(~stringr::str_replace_all(.,"_","-"))
shap_mirnas<-shap_mirnas[ , order(names(shap_mirnas))] # order by column names
mirnas<-colnames(shap_mirnas)

# CpG sites
shap_cpgs <- shap_top %>%
  dplyr::select(starts_with("cg"))
shap_cpgs<-shap_cpgs[ , order(names(shap_cpgs))] # order by column names
cpgs<-colnames(shap_cpgs)

# TCs
shap_tcs <- shap_top %>%
  dplyr::select(starts_with("TC"))
shap_tcs<-shap_tcs[ , order(names(shap_tcs))] # order by column names
tcs<-colnames(shap_tcs)


# Obtain gene annotation for top features and map to entrezid

# Proteins
proteins
load("./data/omics/proteome_subcohort_v5.Rdata")
pro_genes<- proteome_subcohort@featureData@data %>% t() %>% as.data.frame() %>% 
  dplyr::select(any_of(proteins)) %>% t() %>% as.data.frame()
pro_genes<-pro_genes %>%
  dplyr::rename("genes"="Gene_Symbol") %>%
  dplyr::mutate(proteins=rownames(pro_genes)) %>%
  dplyr::select(genes, proteins) %>%
  dplyr::rename("features"="proteins")
# convert gene symbols to entrezid
pro_entrez = bitr(pro_genes$genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
pro_entrez<- pro_entrez %>%
  dplyr::rename("genes"= "SYMBOL")
pro_genes<-merge(pro_genes, pro_entrez, by="genes")

# miRNAs
mirnas
load("./data/omics/mirna_subcohort_notfiltr_v4.Rdata")
mirna_genes<-mirna_subcohort_notfiltr@featureData@data %>% t() %>% as.data.frame() %>%
  dplyr::select(any_of(mirnas)) %>% t() %>% as.data.frame()
mirna_genes[c('first_gene', 'second_gene')] <- str_split_fixed(mirna_genes$Gene_Symbol, ';', 2)
mirna_genes<-mirna_genes %>%
  dplyr::select(first_gene,second_gene) %>%
  dplyr::mutate(mirnas=rownames(mirna_genes))
mirna_genes <- mirna_genes %>% 
  pivot_longer(cols = names(mirna_genes)[1:2],
               values_to = "genes") %>% 
  dplyr::select(mirnas, genes) %>%
  dplyr::rename("features"="mirnas")
mirna_genes<-mirna_genes[!duplicated(mirna_genes$genes),] 
# convert gene symbols to entrezid
mirna_entrez = bitr(mirna_genes$genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
mirna_entrez<- mirna_entrez %>%
  dplyr::rename("genes"= "SYMBOL")
mirna_genes<-merge(mirna_genes, mirna_entrez, by="genes")

# TCs
load("./data/omics/transcriptome_subcohort_f1_v3.Rdata")
tc_genes<-transcriptome_subcohort_f1@featureData@data %>% t() %>% as.data.frame() %>%
  dplyr::select(any_of(tcs)) %>% t() %>% as.data.frame()
tc_genes[c('first_gene', 'second_gene', 'third_gene')] <- str_split_fixed(tc_genes$gene_assignment, ' // ', 3)
tc_genes<-tc_genes %>%
  dplyr::select(second_gene) %>%
  dplyr::mutate(features=rownames(tc_genes)) %>%
  dplyr::rename("genes"="second_gene")
# convert gene symbols to entrezid
tc_entrez = bitr(tc_genes$genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
tc_entrez<- tc_entrez %>%
  dplyr::rename("genes"= "SYMBOL")
tc_genes<-merge(tc_genes, tc_entrez, by="genes")

# CpGs
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")           
library("IlluminaHumanMethylation450kmanifest")
annotation450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# create annotation dataset for methylome
ann450k<-annotation450k %>% 
  as.data.frame() %>% 
  dplyr::rename(chromosome=chr) %>%
  dplyr::filter(chromosome!="chrY") %>%
  dplyr::filter(chromosome!="chrX") %>%
  dplyr::mutate(chromosome_num=chromosome %>%
                  str_remove("chr") %>%
                  as.numeric()) %>% droplevels() %>% rownames_to_column(var="CpG") 
cpg_genes<- ann450k[ann450k$CpG %in% cpgs,] 
cpg_genes[c('first_gene', 'second_gene', 'third_gene', "fourth_gene", "fifth_gene")]<- str_split_fixed(cpg_genes$UCSC_RefGene_Name, ';', 5) 
names(cpg_genes)
cpg_genes<-cpg_genes %>% 
  pivot_longer(cols = names(cpg_genes)[36:40],
               values_to = "genes") %>% 
  dplyr::select(genes,CpG) %>%
  dplyr::rename("features"="CpG")
cpg_genes<-cpg_genes[!duplicated(cpg_genes$genes),]
# convert gene symbols to entrezid
cpg_entrez = bitr(cpg_genes$genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
cpg_entrez<- cpg_entrez %>%
  dplyr::rename("genes"= "SYMBOL")
cpg_genes<-merge(cpg_genes, cpg_entrez, by="genes")

# Run clusterprofiler to conduct pathway over-representation analysis

# Combine all genes in one dataset
features_genes <- rbind(tc_genes,
                        cpg_genes,
                        mirna_genes,
                        pro_genes)
genes<-features_genes$ENTREZID

tc_genes_to_merge<-tc_genes %>% as.data.frame()
write.csv(tc_genes_to_merge, "./results/csvs/tc_genes_C_AB.csv")

# KEGG pathway over-representation analysis
kegg<- enrichKEGG(gene = genes,
                  organism = 'hsa',
                  pvalueCutoff = 0.99,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.1) 
kegg_result<-kegg@result %>%
  dplyr::filter(qvalue<=0.05) %>%
  dplyr::filter(Count>=3) %>%
  dplyr::mutate(GeneRatio= GeneRatio %>% str_replace("/", "_"))
write_csv(kegg_result, "./results/csvs/kegg_result_C_AB.csv")

# GO pathway over-representation analysis
data(geneList, package="DOSE")
go_bp <- enrichGO(gene = genes, 
                  universe = names(geneList),
                  OrgDb = "org.Hs.eg.db",
                  keyType = 'ENTREZID',
                  ont= "BP", 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05,
                  readable = TRUE)

go_bp_result<-go_bp@result %>%
  dplyr::filter(qvalue<0.05) %>%
  dplyr::filter(Count>=3) %>%
  dplyr::mutate(GeneRatio= GeneRatio %>% str_replace("/", "_"))
write.csv(go_bp_result, "./results/csvs/go_result_C_AB.csv")