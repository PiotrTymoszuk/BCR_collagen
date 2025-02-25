# Prediction of BCR at 1, 2, 3, and 5 years after diagnosis by the machine 
# learning models of RFS-free survival. 
# The analysis based on predictor scores is possible only for the Cox-like 
# models, i.e. not possible for the Random Forest model. 
# ROC for survival time points is calculated with the NNE (nearest-neighbor 
# estimate) method implemented in `survivalROC` package.

  insert_head()
  
# container -------
  
  surv_roc <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis data ------
  
  insert_msg('Analysis data')
  
  ## predictor scores of the GBM models, 
  ## RFS-free survival at the identity scale, 
  ## and 0/1-coded BCR event
  
  surv_roc$data <- surv_summary$tertile_data %>% 
    map(map, 
        select, 
        sample_id, predictor_score)
  
  surv_roc$rfs_survival <- surv_globals$data %>% 
    map(select, sample_id, rfs_months, relapse)
  
  for(i in names(surv_roc$data)) {
    
    surv_roc$data[[i]] <- 
      map2(surv_roc$data[[i]], 
           surv_roc$rfs_survival, 
           inner_join, by = 'sample_id')
    
  }
  
  ## transposing: 
  ## performance of different models in one cohort is interesting
  
  surv_roc$data <- transpose(surv_roc$data)
  
  ## time points of interest, months
  
  surv_roc$times <- c(12, 24, 36, 60)
  
  surv_roc$times <- 
    set_names(surv_roc$times, 
              paste0('time_', surv_roc$times))
  
# ROC analysis for survival time points ---------
  
  insert_msg('ROC analysis for survival time points')
  
  for(i in names(surv_roc$times)) {
    
    surv_roc$roc_obj[[i]] <- 
      list(data_lst = surv_roc$data, 
           span = surv_roc$data %>% 
             map(~.x[[1]]) %>% 
             map_dbl(~0.25 * nrow(.x)^(-0.2))) %>% 
      pmap(surv_roc_times, 
           time_variable = 'rfs_months', 
           status_variable = 'relapse', 
           marker_variable = 'predictor_score', 
           predict_time = surv_roc$times[[i]], 
           method = 'NNE')
    
  }

# END -------
  
  surv_roc$data <- NULL
  surv_roc$rfs_survival <- NULL
  
  surv_roc <- compact(surv_roc)
  
  plan('sequential')
  
  insert_tail()