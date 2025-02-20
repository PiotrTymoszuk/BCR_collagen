# Permutation importance of explanatory variables in a survival SVM score.
# The importance is measured by comparing the CV out-of-fold C-indexes between 
# the full model and models with particular explanatory variables re-shuffled 
# at random

  insert_head()
  
# container ------
  
  svm_osimp <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals --------
  
  insert_msg('Analysis globals')
  
  ## tuned parameters of the SVM model
  ## CV folds used also for tuning of the genuine model
  
  svm_osimp$data <- os_globals$data$training
  
  svm_osimp$best_tune <- svm_os$tuning$best_tune
  
  svm_osimp$folds <- svm_os$folds
  
# Importance testing -------
  
  insert_msg('Importance testing')
  
  svm_osimp$test <- os_globals$data$training %>% 
    column_to_rownames('sample_id') %>% 
    svm_importance(data = ., 
                   time_variable = 'os_months', 
                   event_variable = 'death', 
                   folds = svm_osimp$folds, 
                   type = svm_osimp$best_tune$type[[1]], 
                   diff.meth = svm_osimp$best_tune$diff.meth[[1]], 
                   gamma.mu = svm_osimp$best_tune$gamma.mu[[1]], 
                   kernel = svm_osimp$best_tune$kernel[[1]])
  
# Caching the results --------
  
  insert_msg('Caching the results')
  
  svm_osimp <- svm_osimp[c("best_tune", "test")]
  
  save(svm_osimp, file = './cache/svm_osimp.RData')
  
# END ------
  
  plan('sequential')
  
  insert_tail()