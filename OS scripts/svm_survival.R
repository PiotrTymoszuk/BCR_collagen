# Modeling of overall survival as a function of expression 
# of the collagen-related genes in the cancer tissue with 
# Support Vector Machines. 
# The explanatory variables are Z-scores of ComBat-corrected log2 expression 
# estimates (first and second order). 
# The SVM procedure generates predictor scores ('svm_score'), whose association 
# with BCR-free survival is evaluated in the training ant test portion of 
# the GSE16560 cohort is assessed by univariable Cox models.
#
# Pre-processing: normalization, mean centering.
# Tuning and training: in the training portion of the GSE16560 cohort, 
# repeated cross-validation with C-index as cost function.


  insert_head()
  
# container ------
  
  svm_os <- list()
  
# parallel backend -----
  
  insert_msg('Parallel backend')

  plan('multisession')
    
# modeling data -------
  
  insert_msg('Modeling data')

  ## analysis tables: from the globals

  svm_os$data <- os_globals$data %>% 
    map(column_to_rownames, 'sample_id')
  
# Modeling globals ---------
  
  insert_msg('Modeling globals: CV folds and tuning data frames')
  
  ## folds 
  
  set.seed(1234)
  
  svm_os$n_rep <- 5
  
  svm_os$folds <- 1:svm_os$n_rep %>% 
    map(function(x) createFolds(y = factor(svm_os$data$training$death), 
                                k = 10, 
                                list = TRUE, 
                                returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:svm_os$n_rep)) %>% 
    unlist(recursive = FALSE)
  
  ## tune grids: the additive kernel and vanbelle1 method tends to function
  ## the best in terms of C-index as tested per hand. We're tuning
  ## hence just the gamma cost parameter
  
  svm_os$tune_grid <- os_globals$svm_grid
  
# Tuning of the SVM models -------
  
  insert_msg('Tuning of the SVM model')

  svm_os$tuning <- svm_tune(data = svm_os$data$training, 
                              time_variable = 'os_months', 
                              event_variable = 'death', 
                              folds = svm_os$folds, 
                              tune_grid = svm_os$tune_grid)

# Training of the SVM model -------
  
  insert_msg('Training')
  
  svm_os$svm_model <- 
    survivalsvm(formula = Surv(os_months, death) ~ ., 
                data = svm_os$data$training, 
                type = svm_os$tuning$best_tune$type[[1]], 
                diff.meth = svm_os$tuning$best_tune$diff.meth[[1]], 
                gamma.mu = svm_os$tuning$best_tune$gamma.mu[[1]], 
                kernel = svm_os$tuning$best_tune$kernel[[1]])
  
# Predictions: SVM scores --------
  
  insert_msg('Predictions: SVM scores')
  
  ## the scores are actually 'low', 'intermediate', and 'high' risk classes
  
  svm_os$predictions <- svm_os$data %>% 
    map(predict, object = svm_os$svm_model)
  
  svm_os$score_tbl <- 
    map2(svm_os$predictions, 
         svm_os$data, 
         svm_score, 
         time_variable = 'os_months', 
         event_variable = 'death')
  
# Univariable Cox models with the tertiles of SVM score as explanatory variable -------
  
  insert_msg('Univariable Cox models')
  
  svm_os$cox_models <- svm_os$score_tbl %>% 
    map(~call2(.fn = 'coxph', 
               formula = Surv(os_months, death) ~ svm_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., svm_os$score_tbl, as_coxex)
  
# Assumptions, fit stats and inference -------
  
  insert_msg('Assumptions, fit stats and inference')
  
  ## assumptions: quite severe violation in the training subset
  
  svm_os$assumptions <- svm_os$cox_models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  svm_os$stats <- svm_os$cox_models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  svm_os$inference <- svm_os$cox_models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  svm_os[c("assumptions", "stats", "inference")] <- 
    svm_os[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'training', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  ## survival in tertiles of SVM score tertiles
  
  svm_os$calibration <- svm_os$cox_models %>% 
    future_map(calibrate.coxex, 
               n = 3, 
               labels = c('low', 'int', 'high'), 
               .options = furrr_options(seed = TRUE))
  
  svm_os$global_cal <- svm_os$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'training', 'training', 'test'))
  
  ## Brier scores for unique time points
  
  svm_os$brier_scores <- svm_os$cox_models %>% 
    map(surv_brier)
  
# Differences in survival between the SVM score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  ## median survival in the tertiles and Peto-Peto test
  
  svm_os$tertile_stats <- svm_os$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  svm_os$tertile_test <- svm_os$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'training', 'training', 'test'), 
           p_value = pval) %>% 
    re_adjust
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  svm_os$data <- NULL
  svm_os$n_rep <- NULL

  svm_os <- compact(svm_os)
  
  save(svm_os, file = './cache/svm_os.RData')
  
# END -----
  
  plan('sequential')
  
  insert_tail()