# Survival gradient boosted machines: a combined approach with clinical and 
# molecular data. 
# Explanatory variables: first- and second-order terms of ComBat-adjusted 
# log2 expression of the collagen-related genes, Gleason (ISUP), log(PSA) at 
# diagnosis and pT stage (3+ and 2 versus 1).
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
#
# Performance of the GBM model is compared with performance of a 'canonical' 
# Cox model that uses only the clinical predictors (Gleason/ISUP, log PSA, and 
# pT stage): this one is developed with a separate script. 

  insert_head()
  
# container ------
  
  surv_combi <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## analysis data: gene expression and clinical information
  
  surv_combi$data <- 
    map2(surv_globals$data, 
         surv_globals$clinic, 
         inner_join, by = 'sample_id') %>% 
    map(~filter(.x, complete.cases(.x))) %>%
    map(column_to_rownames, 'sample_id')
  
  ## tuning grid for the GBM model: partially optimized per hand
  
  surv_combi$tune_grid <- surv_globals$gbm_grid
  
# CV tuning in the GEO data set -------
  
  insert_msg('CV tuning in the GEO data set')
  
  set.seed(12345)
  
  surv_combi$tuning <- gbm_tune(data = surv_combi$data$geo_pool, 
                              time_variable = 'scaled_rfs_months', 
                              event_variable = 'relapse', 
                              n_folds = 10, 
                              tune_grid = surv_combi$tune_grid, 
                              distribution = 'coxph')
  
# Training of the GBM model in the GEO cohort ------
  
  insert_msg('Training the GEO model')
  
  set.seed(12345)
  
  surv_combi$gbm_model <- 
    gbm(formula = Surv(scaled_rfs_months, relapse) ~ ., 
        data = surv_combi$data$geo_pool, 
        distribution = 'coxph', 
        n.trees = surv_combi$tuning$best_tune$n.trees[1],
        shrinkage = surv_combi$tuning$best_tune$shrinkage[1],
        interaction.depth = surv_combi$tuning$best_tune$interaction.depth[1],
        n.minobsinnode = surv_combi$tuning$best_tune$n.minobsinnode[1],
        cv.folds = 10) 
  
  ## the optimal iteration number for making predictions
  ## as specified by the minimal CV deviance
  
  surv_combi$best_iter <- gbm.perf(surv_combi$gbm_model, method = 'cv')
  
# Predictions and score tables --------
  
  insert_msg('Predictions and score tables')
  
  surv_combi$score_tbl <- surv_combi$data %>% 
    map(predict, 
        object = surv_combi$gbm_model, 
        n.trees = surv_combi$best_iter) %>% 
    map2(surv_combi$data, ., 
         ~mutate(.x[c('scaled_rfs_months', 'relapse')], 
                 gbm_score = .y)) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(as_tibble)
  
# Univariable Cox models -------
  
  insert_msg('Univariable Cox models')
  
  surv_combi$cox_models <- surv_combi$score_tbl %>%
    map(~call2('coxph', 
               formula = Surv(scaled_rfs_months, relapse) ~ gbm_score, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., 
         surv_combi$score_tbl, 
         as_coxex)
  
# Characteristic of the Cox models: assumptions, fit stats and inference -----
  
  insert_msg('Characteristic of the predictor scores')
  
  ## assumptions: mild validations in the training and DKFZ cohorts
  
  surv_combi$assumptions <- surv_combi$cox_models %>% 
    map(summary, type = 'assumptions')
  
  ## fit statistic
  
  surv_combi$stats <- surv_combi$cox_models %>% 
    map(summary, type = 'fit')
  
  ## inference
  
  surv_combi$inference <- surv_combi$cox_models %>%
    map(summary, type = 'inference') %>% 
    map(mutate, 
        estimate = exp(estimate), 
        lower_ci = exp(lower_ci), 
        upper_ci = exp(upper_ci))
  
  ## appending with the cohort information
  
  surv_combi[c("assumptions", "stats", "inference")] <- 
    surv_combi[c("assumptions", "stats", "inference")] %>% 
    map(compress, 
        names_to = 'cohort') %>% 
    map(mutate, 
        dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
# Calibration for the score strata, Nam-D'Agostino method and Brier scores ------
  
  insert_msg('Calibration, D Agostino - Nam and Brier scores')
  
  surv_combi$calibration <- surv_combi$cox_models %>% 
    map(calibrate, n = 3, labels = c('low', 'int', 'high'))
  
  surv_combi$global_cal <- surv_combi$calibration %>% 
    map(summary, type = 'global') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test'))
  
  surv_combi$brier_scores <- surv_combi$cox_models %>% 
    map(surv_brier)
  
# Differences in survival between the score tertiles ------
  
  insert_msg('Differences in survival between the score tertiles')
  
  surv_combi$tertile_stats <- surv_combi$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_median)
  
  surv_combi$tertile_test <- surv_combi$calibration %>% 
    map(~.$surv_fit) %>% 
    map(surv_pvalue, method = 'S1') %>% 
    compress(names_to = 'cohort') %>% 
    mutate(dataset = ifelse(cohort == 'geo_pool', 'training', 'test')) %>% 
    re_adjust(p_variable = 'pval')
  
# Variable importance ------
  
  insert_msg('Variable importance')
  
  surv_combi$importance <- surv_combi$gbm_model %>% 
    summary(n.trees = surv_combi$best_iter) %>% 
    select(var, rel.inf) %>% 
    set_names(c('variable', 'rel_influence')) %>% 
    as_tibble
  
# Caching the results -------
  
  insert_msg('Caching the results')
  
  surv_combi$data <- NULL
  
  surv_combi <- compact(surv_combi)
  
  save(surv_combi, file = './cache/surv_combi.RData')
  
# END ------
  
  insert_tail()