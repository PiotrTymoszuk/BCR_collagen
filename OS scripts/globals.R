# COMBAT expression estimates used for multi-parameter modeling of survival
# in the GSE16560 cohort.
#
# Ther cohort is randomly split into a training and test subset in a 2:1 ratio.
#
# Overall survival (+1 months to avoid zero survival times) is min/max scaled.


  insert_head()
  
# container ------
  
  os_globals <- list()
  
# input data: COMBAT estimates and survival times -----
  
  insert_msg('COMBAT expression estimates and survival times')

  ## survival
  
  os_globals$survival <- gse16560$clinic %>% 
    select(sample_id, os_months, death)
  
  ## genes of interest: first and second order terms, 
  ## Z-scores
  
  os_globals$variables <- globals$genes
  
  os_globals$sq_variables <- 
    paste0(os_globals$variables, '_sq')

  os_globals$expression <- gse16560$combat %>% 
    select(sample_id, all_of(os_globals$variables))

  for(i in os_globals$variables) {
    
    os_globals$expression <- os_globals$expression %>% 
      mutate(!!paste0(i, '_sq') := .data[[i]]^2)
    
  }

  os_globals$expression[c(os_globals$variables, os_globals$sq_variables)] <- 
    os_globals$expression[c(os_globals$variables, os_globals$sq_variables)] %>% 
    map_dfc(zScores)

  ## merging with the survival data
  
  os_globals$expression <- 
    inner_join(os_globals$survival, 
               os_globals$expression, 
               by = 'sample_id') %>% 
    filter(complete.cases(.))

  ## data partition: training and test subset
  
  os_globals$train_idx <- gse16560$train_idx
  
  os_globals$data <- 
    list(training = os_globals$train_idx, 
         test = -os_globals$train_idx) %>% 
    map(~os_globals$expression[.x, ])
  
# N numbers -------
  
  insert_msg('N numbers')
  
  ## n numbers: total and events
  
  os_globals$n_numbers <- os_globals$data %>% 
    map(~tibble(n_total = nrow(.x), 
                n_events = sum(.x$death))) %>% 
    compress(names_to = 'cohort')
  
  ## ready-to-use plot captions
  
  os_globals$n_tags <- os_globals$n_numbers %>% 
    mutate(plot_cap = paste0('total: n = ', n_total, 
                             ', events: n = ', n_events)) %>% 
    .$plot_cap %>% 
    set_names(os_globals$n_numbers$cohort)
  
# CV folds used for tuning of the GLMNET models -----
  
  insert_msg('CV folds')
  
  set.seed(1234)
  
  os_globals$n_rep <- 200 
  
  os_globals$folds <- 1:os_globals$n_rep %>% 
    map(function(x) createFolds(y = factor(os_globals$data$training$death), 
                                k = 10, 
                                list = FALSE, 
                                returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:os_globals$n_rep))
  
# Tuning grid for the GBM, RF, and SVM models --------
  
  insert_msg('Tuning grid for the ML models')

  ## SVM 
  
  os_globals$svm_grid <- 
    expand.grid(gamma.mu = c(0.0001, 0.0025, 0.005, 0.01, 
                             0.015, 0.02, 0.025, 0.05, 
                             seq(0.1, 1, by = 0.1)), 
                type = c('vanbelle1'), 
                diff.meth = c('makediff3'), 
                kernel = c('add_kernel'), 
                stringsAsFactors = FALSE)
  
  ## Random Forest
  
  os_globals$rf_grid <- 
    expand.grid(mtry = seq(2, 27, by = 2), 
                splitrule = c('logrank', 'bs.gradient'), 
                nsplit = c(1, 2, 5), 
                nodesize = c(5, 10, 15, 20), 
                stringsAsFactors = FALSE)
  
  ## GBM
  
  os_globals$gbm_grid <- 
    expand.grid(n.trees = c(500, 1000),
                shrinkage = c(0.01, 0.02, 0.03, 0.04, 0.05),
                interaction.depth = c(2, 3, 4),
                n.minobsinnode = c(2, 5))
  
# END -----

  os_globals$survival <- NULL
  os_globals$expression <- NULL
  
  os_globals <- compact(os_globals)
  
  insert_tail()