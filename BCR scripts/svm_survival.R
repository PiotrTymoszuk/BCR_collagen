# Modeling of biochemical relapse-free survival as a function of expression 
# of the collagen-related genes in the cancer tissue with 
# Support Vector Machines. 
# The explanatory variables are Z-scores of ComBat-corrected log2 expression 
# estimates (first and second order). 
# The SVM procedure generates predictor scores ('svm_score'), whose association 
# with BCR-free survival is evaluated in the training cohort (pooled GEO) 
# and validation collectives is assessed by univariable Cox models.
#
# Pre-processing: normalization, mean centering.
# Tuning and training: in the pooled GEO cohort, repeated cross-validation with 
# C-index as cost function.


  insert_head()
  
# container ------
  
  svm_surv <- list()
  
# parallel backend -----
  
  insert_msg('Parallel backend')

  plan('multisession')
    
# modeling data -------
  
  insert_msg('Modeling data')

  ## analysis tables: from the globals

  svm_surv$data <- surv_globals$data %>% 
    map(column_to_rownames, 'sample_id')
  
# Modeling globals ---------
  
  insert_msg('Modeling globals: CV folds and tuning data frames')
  
  ## folds 
  
  set.seed(1234)
  
  svm_surv$n_rep <- 5
  
  svm_surv$folds <- 1:svm_surv$n_rep %>% 
    map(function(x) createFolds(y = factor(svm_surv$data$geo$relapse), 
                                k = 10, 
                                list = TRUE, 
                                returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:svm_surv$n_rep)) %>% 
    unlist(recursive = FALSE)
  
  ## tune grids: the additive kernel and vanbelle1 method tends to function
  ## the best in terms of C-index as tested per hand. We're tuning
  ## hence just the gamma cost parameter
  
  svm_surv$tune_grid <- surv_globals$svm_grid
  
# Tuning of the SVM models, pooled GEO cohort -------
  
  insert_msg('Tuning of the SVM model in the pooled GEO model')

  svm_surv$tuning <- svm_tune(data = svm_surv$data$geo_pool, 
                              time_variable = 'scaled_rfs_months', 
                              event_variable = 'relapse', 
                              folds = svm_surv$folds, 
                              tune_grid = svm_surv$tune_grid)

# Training of the SVM model in the pooled GEO cohort -------
  
  insert_msg('Training of the GEO model')
  
  svm_surv$svm_model <- 
    survivalsvm(formula = Surv(scaled_rfs_months, relapse) ~ ., 
                data = svm_surv$data$geo, 
                type = svm_surv$tuning$best_tune$type[[1]], 
                diff.meth = svm_surv$tuning$best_tune$diff.meth[[1]], 
                gamma.mu = svm_surv$tuning$best_tune$gamma.mu[[1]], 
                kernel = svm_surv$tuning$best_tune$kernel[[1]])
  
# Predictions: SVM scores --------
  
  insert_msg('Predictions: SVM scores')
  
  ## the scores are actually 'low', 'intermediate', and 'high' risk classes
  
  svm_surv$predictions <- svm_surv$data %>% 
    map(predict, object = svm_surv$svm_model)
  
  svm_surv$score_tbl <- 
    map2(svm_surv$predictions, 
         svm_surv$data, 
         svm_score, 
         time_variable = 'scaled_rfs_months', 
         event_variable = 'relapse')
  
# Univariable Cox models with the tertiles of SVM score as explanatory variable -------
  
  insert_msg('Univariable Cox models')
  
  svm_surv$cox_models <- svm_surv$score_tbl %>% 
    map(~call2(.fn = 'coxph', 
               formula = Surv(scaled_rfs_months, relapse) ~ svm_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., svm_surv$score_tbl, as_coxex)
  
# Assumptions, fit stats and inference -------
  
  insert_msg('Assumptions, fit stats and inference')
  
  ## assumptions: quite severe violation in the training cohort
  
  svm_surv$assumptions <- svm_surv$cox_models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  svm_surv$stats <- svm_surv$cox_models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  svm_surv$inference <- svm_surv$cox_models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  svm_surv[c("assumptions", "stats", "inference")] <- 
    svm_surv[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  ## survival in tertiles of SVM score tertiles
  
  svm_surv$calibration <- svm_surv$cox_models %>% 
    future_map(calibrate.coxex, 
               n = 3, 
               labels = c('low', 'int', 'high'), 
               .options = furrr_options(seed = TRUE))
  
  svm_surv$global_cal <- svm_surv$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
  ## Brier scores for unique time points
  
  svm_surv$brier_scores <- svm_surv$cox_models %>% 
    map(surv_brier)
  
# Differences in survival between the SVM score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  ## median survival in the tertiles and Peto-Peto test
  
  svm_surv$tertile_stats <- svm_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  svm_surv$tertile_test <- svm_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'), 
           p_value = pval) %>% 
    re_adjust
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  svm_surv$data <- NULL
  svm_surv$n_rep <- NULL

  svm_surv <- compact(svm_surv)
  
  save(svm_surv, file = './cache/svm_surv.RData')
  
# END -----
  
  plan('sequential')
  
  insert_tail()