# Performance statistics, Brier scores and survival in score tertiles 
# for multi-parameter
# relapse-free survival models fitted with the Ridge, Elastic Net, LASSO, SVM. 
# Random Forest, and GBM algorithms.

  insert_head()
  
# container -----
  
  surv_summary <- list()
  
# summary globals ------
  
  insert_msg('Summary globals')
  
  surv_summary$mod_expr <-
    expr(list(ridge = ridge_surv, 
              elnet = elnet_surv, 
              lasso = lasso_surv, 
              svm = svm_surv, 
              rf = rf_surv, 
              gbm = gbm_surv, 
              gbm_clinic = surv_clin, 
              gbm_combi = surv_combi))
  
# Data frames with statistics of performance and variable importance -------
  
  insert_msg('Performance stats and variable importance')
  
  ## performance stats, for Random Forest, CI of the C-index
  ## are not available per default
  
  surv_summary$stats <- surv_summary$mod_expr %>% 
    eval %>% 
    map(~.x$stats) %>% 
    map(select, 
        dataset, cohort, 
        c_index, any_of(c('lower_ci', 'upper_ci')), 
        ibs_model, ibs_reference)
  
  surv_summary$stats$rf <- surv_summary$stats$rf %>% 
    mutate(lower_ci = NA, 
           upper_ci = NA)
  
  ## variable importance
  
  surv_summary$importance <- surv_summary$mod_expr %>% 
    eval %>% 
    map(~.x[c('coefs', 'test', 'importance')]) %>% 
    map(compact)
  
  surv_summary$importance$svm <- svm_imp["test"]
  
  surv_summary$importance$rf <- surv_summary$importance$rf$importance["test"]
  
  surv_summary$importance <- surv_summary$importance %>% 
    map(~.x[[1]]) %>% 
    map(filter, variable != 'full')

# Data frames with Brier scores for unique time points ------
  
  insert_msg('Brier scores for the time points')
  
  surv_summary$brier_scores <- surv_summary$mod_expr %>% 
    eval %>% 
    map(~.x$brier_scores)
  
  ## compatible format
  
  surv_summary$brier_scores$rf <- 
    surv_summary$brier_scores$rf %>% 
    map(mutate, 
        training = test, 
        test = NA)

# Normalized predictor scores, score tertiles and tests for differences between the score tertiles -------
  
  insert_msg('Normalized predictor scores, tertiles and tertile tests')
  
  ## this data is available only for the Cox-like models, SVM, and GBM
  ## but not for the Random Forests

  surv_summary$tertile_data <- surv_summary$mod_expr %>% 
    eval %>% 
    map(~.x$score_tbl)
  
  surv_summary$tertile_data <- 
    surv_summary$tertile_data[names(surv_summary$tertile_data) != 'rf']
  
  surv_summary$tertile_data <- surv_summary$tertile_data %>% 
    compact %>% 
    map2(., 
         c(rep('collagen_score', 3), 
           'svm_score', 
           'gbm_score', 
           'clinic_score', 
           'gbm_score'), 
         function(data_lst, var) data_lst %>% 
           map(mutate, 
               score_cuts = cut_tertiles(.data[[var]]))) %>% 
    map(map, 
        select, 
        sample_id, 
        scaled_rfs_months, 
        relapse, 
        score_cuts, 
        any_of(c('collagen_score', 
                 'svm_score', 
                 'clinic_score', 
                 'gbm_score'))) %>% 
    map(map, 
        set_names, 
        c('sample_id', 
          'scaled_rfs_months', 
          'relapse', 
          'score_cuts', 
          'predictor_score'))
  
  ## normalized predictor scores (Z-scores)

  surv_summary$tertile_data <- surv_summary$tertile_data %>% 
    map(map, 
        mutate, 
        z_score = zScores(predictor_score))
    
  ## tertile N numbers: total and events
  
  surv_summary$tertile_n <- surv_summary$tertile_data %>% 
    map(map, count, score_cuts, .drop = FALSE) %>% 
    map(map, set_names, c('score_cuts', 'n_total'))
  
  surv_summary$tertile_events <- surv_summary$tertile_data %>% 
    map(map, filter, relapse == 1) %>% 
    map(map, count, score_cuts, .drop = FALSE) %>% 
    map(map, set_names, c('score_cuts', 'n_events'))
  
  surv_summary$tertile_n <- 
    map2(surv_summary$tertile_n, 
         surv_summary$tertile_events, 
         function(x, y) map2(x, y, left_join, by = 'score_cuts'))
  
  ## tertile surv_fits, median survival times and p values
  
  surv_summary$tertile_fits <- surv_summary$tertile_data %>% 
    map(map, 
        survminer::surv_fit, 
        formula = Surv(scaled_rfs_months, relapse) ~ score_cuts)
  
  surv_summary$tertile_stats <- surv_summary$tertile_fits %>% 
    map(map, surv_median) %>% 
    map(compress, names_to = 'cohort')
  
  surv_summary$tertile_test <- surv_summary$tertile_fits %>% 
    map(map, surv_pvalue, method = 'S1') %>% 
    map(map, re_adjust, 'pval') %>% 
    map(compress, names_to = 'cohort')
  
# comparison of C-indexes between, GBM models (clinic and combined) ------
  
  insert_msg('Comparison of C-indexes, GBM models')
  
  ## with the one-show non-parametric test proposed by Kang 2015
  ## and implemented by compareC package
  
  ## predictor scores 
  
  surv_summary$c_data <- 
    surv_summary$tertile_data[c("gbm_clinic", "gbm_combi")] %>% 
    map(map, select, -score_cuts, -z_score) %>% 
    transpose %>% 
    map(reduce, 
        left_join, 
        by = c('sample_id', 'scaled_rfs_months', 'relapse')) %>% 
    map(set_names, 
        c('sample_id', 
          'scaled_rfs_months', 'relapse', 
          'clinic_score', 'gbm_score'))
  
  ## testing 
  
  surv_summary$c_test <- surv_summary$c_data %>% 
    map(~compareC(timeX = .x$scaled_rfs_months, 
                  statusX = .x$relapse,
                  scoreY = .x$clinic_score, 
                  scoreZ = .x$gbm_score)) %>%
    map(as_tibble) %>% 
    map(mutate, 
        c_index = 1 - est.c)
  
# Inference statistics for the normalized predictor scores --------
  
  insert_msg('Inference stats for normalized predictor scores')

  ## univariable Cox models
  
  for(i in names(surv_summary$tertile_data)) {
    
    surv_summary$z_cox_models[[i]] <- surv_summary$tertile_data[[i]] %>% 
      map(~call2('coxph', 
                 formula = Surv(scaled_rfs_months, relapse) ~ z_score, 
                 data = .x, 
                 x = TRUE, 
                 y = TRUE)) %>% 
      map(eval) %>% 
      map2(surv_summary$tertile_data[[i]], as_coxex)
    
  }
  
  ## the PH assumption is already checked: see the model specific scripts
  ## focusing just on HR and its 95% confidence intervals
  
  surv_summary$z_inference <- surv_summary$z_cox_models %>% 
    map(map, summary, 'inference') %>% 
    map(compress, names_to = 'cohort') %>% 
    compress(names_to = 'algorithm')
  
  ## computation of HRs with 95% confidence intervals
  
  surv_summary$z_inference[, c('hr', 'hr_se', 'hr_lower', 'hr_upper')] <- 
    surv_summary$z_inference[, c("estimate", "se", "lower_ci", "upper_ci")] %>% 
    map_dfc(exp)
  
# END ------
  
  rm(i)
  
  surv_summary$tertile_events <- NULL
  surv_summary$c_data <- NULL
  surv_summary$z_cox_models <- NULL

  surv_summary <- compact(surv_summary)
  
  insert_tail()