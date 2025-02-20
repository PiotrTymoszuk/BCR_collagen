# Permutation importance of explanatory variables in a survival SVM score.
# The importance is measured by comparing the CV out-of-fold C-indexes between 
# the full model and models with particular explanatory variables re-shuffled 
# at random

  insert_head()
  
# container ------
  
  svm_imp <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals --------
  
  insert_msg('Analysis globals')
  
  ## tuned parameters of the SVM model
  ## CV folds used also for tuning of the genuine model

  svm_imp$best_tune <- svm_surv$tuning$best_tune
  
  svm_imp$folds <- svm_surv$folds
  
# Importance testing -------
  
  insert_msg('Importance testing')
  
  svm_imp$test <- surv_globals$data$geo_pool %>% 
    column_to_rownames('sample_id') %>% 
    svm_importance(data = ., 
                   time_variable = 'scaled_rfs_months', 
                   event_variable = 'relapse', 
                   folds = svm_imp$folds, 
                   type = svm_imp$best_tune$type[[1]], 
                   diff.meth = svm_imp$best_tune$diff.meth[[1]], 
                   gamma.mu = svm_imp$best_tune$gamma.mu[[1]], 
                   kernel = svm_imp$best_tune$kernel[[1]])
  
# Caching the results --------
  
  insert_msg('Caching the results')
  
  svm_imp <- svm_imp[c("best_tune", "test")]
  
  save(svm_imp, file = './cache/svm_imp.RData')
  
# END ------
  
  plan('sequential')
  
  insert_tail()