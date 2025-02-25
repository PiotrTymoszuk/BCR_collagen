# Multi-parameter modeling of BCR-free survival with machine learning models 
# trained in the pooled GEO cohort

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
  library(compareC)
  library(survivalROC)

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
  
  c('./BCR scripts/globals.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# multi-parameter modeling of BCR risk -------
  
  insert_msg('Multi-parameter modeling')
  
  ## tuning, training, predictions, and evaluation of performance of 
  ## ML models of BCR-free survival
  
  list(cache_path = c('./cache/elnet_surv.RData', 
                      './cache/ridge_surv.RData', 
                      './cache/lasso_surv.RData', 
                      './cache/svm_surv.RData', 
                      './cache/svm_imp.RData', 
                      './cache/rf_surv.RData', 
                      './cache/gbm_surv.RData', 
                      './cache/surv_combi.RData', 
                      './cache/surv_clin.RData'), 
       script_path = c('./BCR scripts/elastic_net.R', 
                       './BCR scripts/ridge.R', 
                       './BCR scripts/lasso.R', 
                       './BCR scripts/svm_survival.R', 
                       './BCR scripts/svm_importance.R', 
                       './BCR scripts/rf.R', 
                       './BCR scripts/gbm.R', 
                       './BCR scripts/gbm_combi.R', 
                       './BCR scripts/clinic.R'), 
       message = c('Loading chached results of Elastic Net Cox modeling', 
                   'Loading chached results of Ridge Cox modeling', 
                   'Loading chached results of LASSO Cox modeling', 
                   'Loading cached results of SVM modeling', 
                   'Loading cached importance testing for the SVM survival score', 
                   'Loading cached results of RF modeling', 
                   'Loading cached results of GBM modeling', 
                   paste('Loading cached results of GBM modeling,', 
                         'combined collagen/clinic model'), 
                   paste('Loading cached results of GBM modeling', 
                         'with clinical predictors'))) %>% 
    pwalk(access_cache)
  
  ## summary statistics, 
  ## ROC for selected, clinically relevant survival time points
  ## plots of the model evaluation and comparison of model types
  
  c('./BCR scripts/summary.R', 
    './/BCR scripts/roc.R', 
    './BCR scripts/plots.R', 
    './BCR scripts/clinical_plots.R', 
    './BCR scripts/cohort.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END --------
  
  insert_tail()