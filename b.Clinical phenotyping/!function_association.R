##############################################################
##### Function to run association of clusters with outcomes
##############################################################

association<-function(cluster, dataframe) {
  
  # Continuous z outcomes ##############################################################
  cont_z_outc_fit <- lapply(cont_z_outc, 
                            
                            function(y) {
                              model <- as.formula(paste0(y, "~", cluster, "+cohort"))
                              fit <- glm(model, family=gaussian,
                                         data=dataframe) 
                              
                              coefficients<- cbind( y, as.data.frame(summary(fit)$coefficients)[2:3,],
                                                    as.data.frame(summary(fit)$coefficients[2:3,1]-qnorm(0.975)*summary(fit)$coefficients[2:3,2]),
                                                    as.data.frame(summary(fit)$coefficients[2:3,1]+qnorm(0.975)*summary(fit)$coefficients[2:3,2]), 
                                                    as.data.frame(rownames(fit$R)[2:3])
                                                    
                              )
                              return(coefficients) 
                            } )
  
  library(plyr)
  cont_z_outc_cluster<- ldply (cont_z_outc_fit, data.frame)

  results_cont_z_outc_cluster<-cont_z_outc_cluster %>%
    dplyr::rename("outcome"= "y",
                  "estimate" ="Estimate",
                  "se" = "Std..Error",
                  "z.tvalue" = "t.value",
                  "pvalue"= "Pr...t..",
                  "conf.low"= "summary.fit..coefficients.2.3..1....qnorm.0.975....summary.fit..coefficients.2.3..2.",
                  "conf.high"= "summary.fit..coefficients.2.3..1....qnorm.0.975....summary.fit..coefficients.2.3..2..1",
                  "term"= "rownames.fit.R..2.3.") %>% 
    mutate_if(is.numeric, ~round(., 2))
  
  # Categorical outcomes ##############################################################
  cat_outc_fit <- lapply(cat_outc, 
                         
                         function(y) {
                           model <- as.formula(paste0(y, "~", cluster, "+cohort"))
                           fit <- glm(model, family=binomial,
                                      data=dataframe) 
                           
                           coefficients<- cbind( y, as.data.frame(summary(fit)$coefficients)[2:3,],
                                                 as.data.frame(summary(fit)$coefficients[2:3,1]-qnorm(0.975)*summary(fit)$coefficients[2:3,2]),
                                                 as.data.frame(summary(fit)$coefficients[2:3,1]+qnorm(0.975)*summary(fit)$coefficients[2:3,2]),
                                                 as.data.frame(rownames(fit$R)[2:3])
                                                 
                           )
                           return(coefficients) 
                         } )
  
  library(plyr)
  cat_outc_cluster<- ldply (cat_outc_fit, data.frame)
  
  # Rename column names
  results_cat_outc_cluster<-cat_outc_cluster %>%
    dplyr::rename("outcome"= "y",
                  "estimate" ="Estimate",
                  "se" = "Std..Error",
                  "z.tvalue" = "z.value",
                  "pvalue"= "Pr...z..",
                 "conf.low"= "summary.fit..coefficients.2.3..1....qnorm.0.975....summary.fit..coefficients.2.3..2.",
                  "conf.high"= "summary.fit..coefficients.2.3..1....qnorm.0.975....summary.fit..coefficients.2.3..2..1",
                  "term"= "rownames.fit.R..2.3.")
  
  # transform estimates to ORs
  results_cat_outc_cluster <-results_cat_outc_cluster %>%
    dplyr::mutate(estimate=exp(results_cat_outc_cluster$estimate),
                  conf.low=exp(results_cat_outc_cluster$conf.low),
                  conf.high=exp(results_cat_outc_cluster$conf.high))%>% 
    mutate_if(is.numeric, ~round(., 2))
  
  results=rbind(results_cont_z_outc_cluster, results_cat_outc_cluster)
  
  return(results)
  
}

##############################################################
##### Function to run stratified analysis by obesity status ##
##############################################################
association_strat= function(cluster,dataframe, N) {
  
  if (N==0) {
    df<-dataframe %>%
      dplyr::filter(hs_ovob==0)
    normal_weight=association(cluster,df)
    return(normal_weight)
  }else{
    df<-dataframe %>%
      dplyr::filter(hs_ovob==1)
    overweight = association(cluster,df)
    return(overweight)
  }
}

###################################################################
#Function to run association of clusters with adolescent outcomes
###################################################################
association_adolescence<-function(cluster, dataframe) {
  
  # Continuous raw outcomes ###
  cont_outc_fit <- lapply(cont_outc_r, 
                          
                          function(y) {
                            model <- as.formula(paste0(y, "~", cluster, "+cohort+hs2_age+h_sex"))
                            fit <- glm(model, family=gaussian, weights = weights_atr,
                                       data=dataframe) 
                            
                            model_marginal= avg_comparisons(fit,
                                                            variables = cluster,
                                                            vcov = "HC3",
                                                            wts = "weights_atr")
                            
                            
                            coefficients<- cbind( y, as.data.frame(model_marginal[1:2,2:9])
                            )
                            return(coefficients) 
                          } )
  
  library(plyr)
  cont_outc_cluster<- ldply (cont_outc_fit, data.frame) %>% 
    mutate_if(is.numeric, ~round(., 2))
  
  # Continuous z outcomes #########################
  cont_z_outc_fit <- lapply(cont_outc_z, 
                            
                            function(y) {
                              model <- as.formula(paste0(y, "~", cluster, "+cohort"))
                              fit <- glm(model, family=gaussian,weights = weights_atr,
                                         data=dataframe) 
                              
                              
                              model_marginal= avg_comparisons(fit,
                                                              variables = cluster,
                                                              vcov = "HC3",
                                                              wts = "weights_atr")
                              
                              coefficients<- cbind( y, as.data.frame(model_marginal[1:2,2:9])
                              )
                              return(coefficients) 
                            } )
  
  library(plyr)
  cont_z_outc_cluster<- ldply (cont_z_outc_fit, data.frame) %>% 
    mutate_if(is.numeric, ~round(., 2))
  
  count_cont<-rbind(cont_z_outc_cluster, cont_outc_cluster) %>%
    dplyr::select(y, contrast, estimate, p.value, conf.low,	conf.high)
  
  # Categorical outcomes ############################
  cat_outc_fit <- lapply(cat_outc, 
                         
                         function(y) {
                           model <- as.formula(paste0(y, "~", cluster, "+cohort"))
                           fit <- glm(model, family=binomial,weights = weights_atr,
                                      data=dataframe) 
                           
                           model_marginal= avg_comparisons(fit,
                                                           variables = cluster,
                                                           wts = "weights_atr",
                                                           comparison = "lnoravg",
                                                           transform = "exp")
                           
                           
                           coefficients<- cbind( y, as.data.frame(model_marginal[1:2,2:7]))
                           
                           return(coefficients) 
                         } )
  
  library(plyr)
  cat_outc_cluster<- ldply (cat_outc_fit, data.frame) %>% 
    mutate_if(is.numeric, ~round(., 2)) %>%
    dplyr::select(y, contrast, estimate, p.value, conf.low,	conf.high)
  
  
  results=rbind(count_cont, cat_outc_cluster)  %>%
    dplyr::rename("outcome"= "y") %>%
    dplyr::mutate(term=case_when( grepl("Cluster B - Cluster A", contrast) ~ "Cluster B",
                                  grepl("Cluster C - Cluster A", contrast) ~ "Cluster C",
                                  contrast=="ln(odds(Cluster B) / odds(Cluster A))" ~ "Cluster B",
                                  contrast=="ln(odds(Cluster C) / odds(Cluster A))" ~ "Cluster C"))
  
  return(results)
  
}




