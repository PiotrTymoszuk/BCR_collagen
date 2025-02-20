# Survival random forests. 
# Explanatory variables: first- and second-order terms of ComBat-adjusted 
# log2 expression of the collagen-related genes. 
# Response: overall survival, with the time variable subjected to min/max 
# scaling. 
# Tuning and training in the training portion of the GSE16560 cohort; 
# tuning in an out-of-bag setting by maximizing C-index of the model.
# The trained model generates probabilities of survival at unique time points 
# for the training collective and test cohorts (TCGA and DKFZ). 
# Evaluation: C-index and IBS of the model in the training and test cohorts.  
# Variable importance: permutation importance as proposed by Breiman.

  insert_head()
  
# container -------
  
  rf_os <- list()
  
# Modeling data and analysis globals -------
  
  insert_msg('Modeling data and analysis globals')

  ## analysis tables: from the globals
  
  rf_os$data <- os_globals$data %>% 
    map(column_to_rownames, 'sample_id')

  ## NULL models of BCR-free survival
  
  rf_os$null_models <- rf_os$data %>% 
    map(~call2(.fn = 'coxph', 
               formula = Surv(os_months, death) ~ 1, 
               data = .x, 
               x = TRUE, 
               y = TRUE)) %>% 
    map(eval) %>% 
    map2(., rf_os$data, as_coxex)

  ## tuning grid 
  
  rf_os$tune_grid <- os_globals$rf_grid
  
# Construction of the tuning models --------
  
  insert_msg('construction of the tuning models')
  
  set.seed(1234)
  
  rf_os$tuning$models <- rf_os$tune_grid %>% 
    pmap(rfsrc, 
         formula = Surv(os_months, death) ~ ., 
         data = rf_os$data$training, 
         save.memory = TRUE, 
         ntree = 500)
  
  rf_os$tuning$models <- rf_os$tuning$models %>% 
    set_names(paste0('cond_', 1:length(rf_os$tuning$models)))
  
  rf_os$tune_grid$condition <- 
    paste0('cond_', 1:length(rf_os$tuning$models))
  
# Extraction of the performance stats from the tuning models -------
  
  insert_msg('Extraction of the tuning stats')
  
  ## C-indexes for the OOB predictions: they are stored as 1 - error rate
  
  rf_os$tuning$summary <- rf_os$tuning$models %>% 
    map(~.x$err.rate) %>% 
    map(na.omit) %>% 
    map(as.numeric) %>% 
    map_dbl(function(x) 1 - x) %>% 
    compress(names_to = 'condition', 
             values_to = 'c_index')
  
  rf_os$tuning$summary <- 
    left_join(rf_os$tune_grid, 
              rf_os$tuning$summary, 
              by = 'condition')
  
  ## the best tune
  
  rf_os$tuning$best_tune <- rf_os$tuning$summary %>% 
    filter(c_index == max(c_index))
  
  rf_os$tuning$best_tune <- 
    rf_os$tuning$best_tune[1, ]
  
# Training the RF model --------
  
  insert_msg('Training the RF model')
  
  set.seed(1234)
  
  rf_os$rf_model <- 
    rfsrc(formula = Surv(os_months, death) ~ ., 
          data = rf_os$data$training, 
          save.memory = FALSE, 
          ntree = 500, 
          mtry = rf_os$tuning$best_tune$mtry[1], 
          splitrule = rf_os$tuning$best_tune$splitrule[1],
          nodesize = rf_os$tuning$best_tune$nodesize[1])
  
# Predictions --------
  
  insert_msg('Predictions')
  
  rf_os$predictions <- 
    list(newdata = rf_os$data, 
         outcome = c('train', 'test')) %>% 
    pmap(predict, 
         object = rf_os$rf_model)
  
# Prediction stats ------
  
  insert_msg('Prediction stats')
  
  ## C-indexes, extracted as above for the tuning models
  
  rf_os$stats$c_index <- rf_os$predictions %>% 
    map(~.x$err.rate) %>% 
    map(na.omit) %>% 
    map(as.numeric) %>% 
    map_dbl(function(x) 1 - x)
  
  ## IBS for the models
  
  rf_os$stats$ibs_model <- rf_os$predictions %>% 
    map(get.brier.survival, 
        cens.model = 'rfsrc') %>% 
    map_dbl(~.x$crps.std)
  
  ## IBS for the dummy models: reference
  
  rf_os$stats$ibs_reference <- rf_os$null_models %>% 
    map(summary, 'fit') %>% 
    map_dbl(~.x$ibs_reference[1])
  
  ## a common table
  
  rf_os$stats <- rf_os$stats %>% 
    as_tibble %>%
    mutate(cohort = names(rf_os$stats[[1]]), 
           dataset = ifelse(cohort == 'training', 'training', 'test'))
  
# Brier scores for unique time points -------
  
  insert_msg('Brier scores for the unique time points')
  
  ## Brier scores for the model
  
  rf_os$brier_scores <- rf_os$predictions %>% 
    map(get.brier.survival, 
        cens.model = 'rfsrc') %>% 
    map(~.x$brier.score) %>% 
    map(as_tibble) %>% 
    map(set_names, c('time', 'test'))
  
  ## Brier scores for the reference, i.e. dummy Cox model 
  ## for the training data
  
  rf_os$ref_brier_scores <- rf_os$null_models$training %>% 
    surv_brier %>% 
    select(time, reference)
  
  rf_os$brier_scores <- rf_os$brier_scores %>% 
    map(mutate, 
        training = rf_os$brier_scores$training$test) %>% 
    map(left_join, rf_os$ref_brier_scores, 
        by = 'time') %>% 
    map(as.list) %>% 
    transpose %>% 
    pmap(brier)
  
# variable importance --------
  
  insert_msg('Variable importance')
  
  ## constructing a model with permutation importance measures

  rf_os$importance$model <- 
    rfsrc(formula = Surv(os_months, death) ~ ., 
          data = rf_os$data$training, 
          save.memory = FALSE, 
          ntree = 500, 
          mtry = rf_os$tuning$best_tune$mtry[1], 
          splitrule = rf_os$tuning$best_tune$splitrule[1],
          nodesize = rf_os$tuning$best_tune$nodesize[1], 
          importance = 'permute')
  
  ## extracting the importance stats, i.e. deltas of C-indexes
  
  rf_os$importance$test <- rf_os$importance$model$importance %>% 
    compress(names_to = 'variable', 
             values_to = 'delta_c_index')
  
# Caching the results ------
  
  insert_msg('Caching the results')
  
  rf_os$data <- NULL
  rf_os$ref_brier_scores <- NULL
  rf_os$null_models <- NULL
  
  rf_os$tuning$models <- NULL
  rf_os$tuning <- compact(rf_os$tuning)
  
  rf_os$importance$model <- NULL
  rf_os$importance <- compact(rf_os$importance)

  rf_os <- compact(rf_os)
  
  save(rf_os, file = './cache/rf_os.RData')
  
# END ------
  
  insert_tail()
  
  