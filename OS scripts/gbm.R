# Survival gradient boosted machines. 
# Explanatory variables: first- and second-order terms of ComBat-adjusted 
# log2 expression of the collagen-related genes. 
# Response: overall survival, with the time variable subjected to min/max 
# scaling. 
# Tuning and training in the TCGA cohort; tuning in a 10-fold CV setting 
# by minimizing deviance of the model.
# The trained model generates 
# predictor scores ('gbm_score') for the training collective and test cohort 
# GSE16560. 
# Evaluation: evaluation of performance of univariable Cox models of BCR-free 
# survival as a function of the GBM predictor score. 
# Variable importance: tree importance; reduction of the model error attributed 
# to the variable included in the learners.

insert_head()

# container -------

  gbm_os <- list()

# Modeling data -------

  insert_msg('Modeling data')
  
  ## analysis tables: from the globals
  
  gbm_os$data <- os_globals$data %>% 
    map(column_to_rownames, 'sample_id')

  ## tuning grid: partially optimized per hand
  
  gbm_os$tune_grid <- os_globals$gbm_grid
  
# CV tuning in the GEO data set -------
  
  insert_msg('CV tuning in the GEO data set')

  gbm_os$tuning <- gbm_tune(data = gbm_os$data$training, 
                            time_variable = 'os_months', 
                            event_variable = 'death', 
                            n_folds = 10, 
                            tune_grid = gbm_os$tune_grid, 
                            distribution = 'coxph')
  
# Training of the GBM model in the GEO cohort ------
  
  insert_msg('Training the GEO model')
  
  set.seed(12345)
  
  gbm_os$gbm_model <- 
    gbm(formula = Surv(os_months, death) ~ ., 
        data = gbm_os$data$training, 
        distribution = 'coxph', 
        n.trees = gbm_os$tuning$best_tune$n.trees[1],
        shrinkage = gbm_os$tuning$best_tune$shrinkage[1],
        interaction.depth = gbm_os$tuning$best_tune$interaction.depth[1],
        n.minobsinnode = gbm_os$tuning$best_tune$n.minobsinnode[1],
        cv.folds = 10) 
  
  ## the optimal iteration number for making predictions
  ## as specified by the minimal CV deviance
  
  gbm_os$best_iter <- gbm.perf(gbm_os$gbm_model, method = 'cv')
  
# Predictions and score tables --------
  
  insert_msg('Predictions and score tables')
  
  gbm_os$score_tbl <- gbm_os$data %>% 
    map(predict, 
        object = gbm_os$gbm_model, 
        n.trees = gbm_os$best_iter) %>% 
    map2(gbm_os$data, ., 
         ~mutate(.x[c('os_months', 'death')], 
                 gbm_score = .y)) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(as_tibble)

# Univariable Cox models -------
  
  insert_msg('Univariable Cox models')
  
  gbm_os$cox_models <- gbm_os$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(os_months, death) ~ gbm_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         gbm_os$score_tbl, 
         as_coxex)
  
# Characteristic of the Cox models: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of GBM scores')
  
  ## assumptions: mild validations in the training and DKFZ cohorts
  
  gbm_os$assumptions <- gbm_os$cox_models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  gbm_os$stats <- gbm_os$cox_models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  gbm_os$inference <- gbm_os$cox_models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  gbm_os[c("assumptions", "stats", "inference")] <- 
    gbm_os[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'training', 'training', 'test'))

# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  gbm_os$calibration <- gbm_os$cox_models %>% 
    map(calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  gbm_os$global_cal <- gbm_os$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'training', 'training', 'test'))
  
  gbm_os$brier_scores <- gbm_os$cox_models %>% 
    map(surv_brier)
  
# Differences in survival between the score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  gbm_os$tertile_stats <- gbm_os$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  gbm_os$tertile_test <- gbm_os$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'training', 'training', 'test')) %>% 
    re_adjust(p_variable = 'pval')
  
# Variable importance ------
  
  insert_msg('Variable importance')
  
  gbm_os$importance <- gbm_os$gbm_model %>% 
    summary(n.trees = gbm_os$best_iter) %>% 
    select(var, rel.inf) %>% 
    set_names(c('variable', 'rel_influence')) %>% 
    as_tibble
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  gbm_os$data <- NULL
  
  gbm_os <- compact(gbm_os)
  
  save(gbm_os, file = './cache/gbm_os.RData')
  
# END ------
  
  insert_tail()