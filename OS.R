# Multi-parameter modeling of overall survival in the TCGA and GSE16560 cohorts 
# with log2 expression values of the collagen-related genes

# tools -------

  library(tidyverse)
  library(rlang)
  library(trafo)
  library(stringi)
  
  library(glmnet)
  library(survivalsvm)
  library(randomForestSRC)
  library(gbm)
  library(caret)
  
  library(survival)
  library(survminer)
  library(coxExtensions)
  library(kmOptimizer)
  
  library(furrr)
  library(soucer)
  
  library(ggtext)
  
  explore <- exda::explore
  select <- dplyr::select
  reduce <- purrr::reduce
  set_rownames <- trafo::set_rownames
  
  c('./tools/globals.R', 
    './tools/functions.R', 
    './tools/svm_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)

# analysis globals --------

  insert_msg('Analysis globals')
  
  c('./OS scripts/globals.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# Multi-parameter modeling of overall survival ------
  
  insert_msg('Multi-parameter modeling of overall survival')
  
  ## tuning, training, predictions, and evaluation of performance of 
  ## ML models of BCR-free survival
  
  list(cache_path = c('./cache/elnet_os.RData', 
                      './cache/ridge_os.RData', 
                      './cache/lasso_os.RData', 
                      './cache/svm_os.RData', 
                      './cache/svm_osimp.RData', 
                      './cache/rf_os.RData', 
                      './cache/gbm_os.RData'), 
       script_path = c('./OS scripts/elastic_net.R', 
                       './OS scripts/ridge.R', 
                       './OS scripts/lasso.R', 
                       './OS scripts/svm_survival.R', 
                       './OS scripts/svm_importance.R', 
                       './OS scripts/rf.R', 
                       './OS scripts/gbm.R'), 
       message = c('Loading chached results of Elastic Net Cox modeling', 
                   'Loading chached results of Ridge Cox modeling', 
                   'Loading chached results of LASSO Cox modeling', 
                   'Loading cached results of SVM modeling', 
                   'Loading cached importance testing for the SVM survival score', 
                   'Loading cached results of RF modeling', 
                   'Loading cached results of GBM modeling')) %>% 
    pwalk(access_cache)
  
  ## summary of the modeling results and visualizations
  
  c('./OS scripts/summary.R', 
    './OS scripts/plots.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END -------
  
  insert_tail()