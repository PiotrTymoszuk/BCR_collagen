# Performance statistics, Brier scores and survival in score tertiles 
# for multi-parameter models of overall survival fitted with the Ridge, 
# Elastic Net, LASSO, SVM. Random Forest, and GBM algorithms.

  insert_head()
  
# container -----
  
  os_summary <- list()
  
# summary globals ------
  
  insert_msg('Summary globals')
  
  ## expression for a list with modeling results
  
  os_summary$mod_exprs <- 
    expr(list(ridge = ridge_os, 
              elnet = elnet_os, 
              lasso = lasso_os, 
              svm = svm_os, 
              rf = rf_os, 
              gbm = gbm_os))
  
# Data frames with statistics of performance and variable importance -------
  
  insert_msg('Performance stats and variable importance')
  
  ## performance stats, for Random Forest, CI of the C-index
  ## are not available per default
  
  os_summary$stats <- os_summary$mod_exprs %>% 
    eval %>% 
    map(~.x$stats) %>% 
    map(select, 
        dataset, cohort, 
        c_index, any_of(c('lower_ci', 'upper_ci')), 
        ibs_model, ibs_reference)
  
  os_summary$stats$rf <- os_summary$stats$rf %>% 
    mutate(lower_ci = NA, 
           upper_ci = NA)
  
  ## variable importance
  
  os_summary$importance <- os_summary$mod_exprs %>% 
    eval %>% 
    map(~.x[c('coefs', 'test', 'importance')]) %>% 
    map(compact) 
  
  os_summary$importance$svm <- svm_imp['test']
  
  os_summary$importance$rf <- os_summary$importance$rf$importance["test"]
    
  os_summary$importance <- os_summary$importance %>% 
    map(~.x[[1]]) %>% 
    map(filter, variable != 'full')

# Data frames with Brier scores for unique time points ------
  
  insert_msg('Brier scores for the time points')
  
  os_summary$brier_scores <- os_summary$mod_exprs %>% 
    eval %>% 
    map(~.x$brier_scores)
  
  ## compatible format
  
  os_summary$brier_scores$rf <- 
    os_summary$brier_scores$rf %>% 
    map(mutate, 
        training = test, 
        test = NA)

# Score tertiles and tests for differences between the score tertiles -------
  
  insert_msg('Tertiles and tertile tests')
  
  ## this data is available only for the Cox-like models, SVM, and GBM
  ## but not for the Random Forests

  os_summary$tertile_data <- os_summary$mod_exprs %>% 
    eval %>% 
    map(~.x$score_tbl) %>% 
    map(compact)
  
  os_summary$tertile_data <- 
    os_summary$tertile_data[names(os_summary$tertile_data) != 'rf']
  
  os_summary$tertile_data <- os_summary$tertile_data %>% 
    map2(., 
         c(rep('collagen_score', 3), 
           'svm_score', 
           'gbm_score'), 
         function(data_lst, var) data_lst %>% 
           map(mutate, 
               score_cuts = cut_tertiles(.data[[var]]))) %>% 
    map(map, 
        select, 
        sample_id, os_months, death, score_cuts)
  
  ## tertile N numbers: total and events
  
  os_summary$tertile_n <- os_summary$tertile_data %>% 
    map(map, count, score_cuts, .drop = FALSE) %>% 
    map(map, set_names, c('score_cuts', 'n_total'))
  
  os_summary$tertile_events <- os_summary$tertile_data %>% 
    map(map, filter, death == 1) %>% 
    map(map, count, score_cuts, .drop = FALSE) %>% 
    map(map, set_names, c('score_cuts', 'n_events'))
  
  os_summary$tertile_n <- 
    map2(os_summary$tertile_n, 
         os_summary$tertile_events, 
         function(x, y) map2(x, y, left_join, by = 'score_cuts'))
  
  ## tertile surv_fit objects, median survival times and p values
  
  os_summary$tertile_fits <- os_summary$tertile_data %>% 
    map(map, 
        survminer::surv_fit, 
        formula = Surv(os_months, death) ~ score_cuts)
  
  os_summary$tertile_stats <- os_summary$tertile_fits %>% 
    map(map, surv_median) %>% 
    map(compress, names_to = 'cohort')
  
  os_summary$tertile_test <- os_summary$tertile_fits %>% 
    map(map, surv_pvalue, method = 'S1') %>% 
    map(map, re_adjust, 'pval') %>% 
    map(compress, names_to = 'cohort')
  
# END ------
  
  os_summary$tertile_events <- NULL
  
  os_summary <- compact(os_summary)
  
  insert_tail()