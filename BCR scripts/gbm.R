# Survival gradient boosted machines. 
# Explanatory variables: first- and second-order terms of ComBat-adjusted 
# log2 expression of the collagen-related genes. 
# Response: BCR-free survival, with the time variable subjected to min/max 
# scaling. 
# Tuning and training in the pooled GEO cohort; tuning in a 10-fold CV setting 
# by minimizing deviance of the model.
# The trained model generates 
# predictor scores ('gbm_score') for the training collective and test cohorts 
# (TCGA and DKFZ). 
# Evaluation: evaluation of performance of univariable Cox models of BCR-free 
# survival as a function of the GBM predictor score. 
# Variable importance: tree importance; reduction of the model error attributed 
# to the variable included in the learners.

insert_head()

# container -------

  gbm_surv <- list()

# Modeling data -------

  insert_msg('Modeling data')
  
  ## analysis tables: from the globals
  
  gbm_surv$data <- surv_globals$data %>% 
    map(column_to_rownames, 'sample_id')

  ## tuning grid: partially optimized per hand
  
  gbm_surv$tune_grid <- surv_globals$gbm_grid
  
# CV tuning in the GEO data set -------
  
  insert_msg('CV tuning in the GEO data set')
  
  set.seed(12345)
  
  gbm_surv$tuning <- gbm_tune(data = gbm_surv$data$geo_pool, 
                              time_variable = 'scaled_rfs_months', 
                              event_variable = 'relapse', 
                              n_folds = 10, 
                              tune_grid = gbm_surv$tune_grid, 
                              distribution = 'coxph')
  
# Training of the GBM model in the GEO cohort ------
  
  insert_msg('Training the GEO model')
  
  set.seed(12345)
  
  gbm_surv$gbm_model <- 
    gbm(formula = Surv(scaled_rfs_months, relapse) ~ ., 
        data = gbm_surv$data$geo_pool, 
        distribution = 'coxph', 
        n.trees = gbm_surv$tuning$best_tune$n.trees[1],
        shrinkage = gbm_surv$tuning$best_tune$shrinkage[1],
        interaction.depth = gbm_surv$tuning$best_tune$interaction.depth[1],
        n.minobsinnode = gbm_surv$tuning$best_tune$n.minobsinnode[1],
        cv.folds = 10) 
  
  ## the optimal iteration number for making predictions
  ## as specified by the minimal CV deviance
  
  gbm_surv$best_iter <- gbm.perf(gbm_surv$gbm_model, method = 'cv')
  
# Predictions and score tables --------
  
  insert_msg('Predictions and score tables')
  
  gbm_surv$score_tbl <- gbm_surv$data %>% 
    map(predict, 
        object = gbm_surv$gbm_model, 
        n.trees = gbm_surv$best_iter) %>% 
    map2(gbm_surv$data, ., 
         ~mutate(.x[c('scaled_rfs_months', 'relapse')], 
                 gbm_score = .y)) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(as_tibble)

# Univariable Cox models -------
  
  insert_msg('Univariable Cox models')
  
  gbm_surv$cox_models <- gbm_surv$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(scaled_rfs_months, relapse) ~ gbm_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         gbm_surv$score_tbl, 
         as_coxex)
  
# Characteristic of the Cox models: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of GBM scores')
  
  ## assumptions: mild validations in the training and DKFZ cohorts
  
  gbm_surv$assumptions <- gbm_surv$cox_models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  gbm_surv$stats <- gbm_surv$cox_models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  gbm_surv$inference <- gbm_surv$cox_models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  gbm_surv[c("assumptions", "stats", "inference")] <- 
    gbm_surv[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))

# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  gbm_surv$calibration <- gbm_surv$cox_models %>% 
    map(calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  gbm_surv$global_cal <- gbm_surv$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
  gbm_surv$brier_scores <- gbm_surv$cox_models %>% 
    map(surv_brier)
  
# Differences in survival between the score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  gbm_surv$tertile_stats <- gbm_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  gbm_surv$tertile_test <- gbm_surv$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test')) %>% 
    re_adjust(p_variable = 'pval')
  
# Variable importance ------
  
  insert_msg('Variable importance')
  
  gbm_surv$importance <- gbm_surv$gbm_model %>% 
    summary(n.trees = gbm_surv$best_iter) %>% 
    select(var, rel.inf) %>% 
    set_names(c('variable', 'rel_influence')) %>% 
    as_tibble
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  gbm_surv$data <- NULL
  
  gbm_surv <- compact(gbm_surv)
  
  save(gbm_surv, file = './cache/gbm_surv.RData')
  
# END ------
  
  insert_tail()