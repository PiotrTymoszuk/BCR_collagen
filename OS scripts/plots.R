# Plots for the results of multi-parameter survival modeling

  insert_head()
  
# container ------
  
  os_plots <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# Plotting globals ------
  
  insert_msg('Plot globals')
  
  os_plots$plot_titles <- 
    paste0(globals$study_labels["gse16560"], 
           c(', training subset', ', test subset'))

# Tuning process -------
  
  insert_msg('Plots of the tuning process')
  
  ## for the regularized Cox models
  
  os_plots$tuning[c('ridge', 'elnet', 'lasso')] <- 
    list(ridge_os, elnet_os, lasso_os) %>% 
    map(~.x$lambda_tune) %>% 
    map(plot)
    
  os_plots$tuning[c('ridge', 'elnet', 'lasso')] <- 
    list(x = os_plots$tuning[c('ridge', 'elnet', 'lasso')], 
         y = globals$algo_labels[c("ridge", "elnet", "lasso")]) %>% 
    pmap(function(x, y) x + 
           globals$common_theme + 
           labs(title = y, 
                y = 'deviance, repeated CV'))

  ## for SVM: 
  
  os_plots$tuning$svm <- svm_os$tuning$summary %>% 
    mutate(best = ifelse(c_index == max(c_index), 
                         'yes', 'no')) %>% 
    ggplot(aes(x = gamma.mu, 
               y = c_index, 
               fill = best)) + 
    geom_point(shape = 21, 
               size = 2) + 
    scale_x_continuous(trans = 'log', 
                       labels = function(x) signif(x, 2)) + 
    scale_fill_manual(values = c(no = 'steelblue', 
                                 yes = 'coral3'), 
                      name = 'best tune') + 
    globals$common_theme + 
    labs(title = globals$algo_labels["svm"], 
         subtitle = paste('Best tune: \u03B3 =', 
                          signif(svm_os$tuning$best_tune$gamma.mu, 3)), 
         x = expression(gamma),
         y = 'C-index, CV')
  
  ## for RF: 
  
  os_plots$tuning$rf <- rf_os$tuning$summary %>% 
    mutate(best = ifelse(c_index == max(c_index), 
                         'yes', 'no')) %>% 
    ggplot(aes(x = mtry, 
               y = c_index, 
               color = splitrule, 
               fill = best)) + 
    geom_line(aes(group = splitrule)) + 
    geom_point(shape = 21, 
               size = 1, 
               color = 'black') + 
    scale_fill_manual(values = c(no = 'steelblue', 
                                yes = 'coral3'), 
                     name = 'best tune') + 
    facet_grid(nsplit ~ nodesize) + 
    globals$common_theme + 
    labs(title = globals$algo_labels["rf"], 
         subtitle = paste0('Best tune: ', 
                           'mtry = ', signif(rf_surv$tuning$best_tune$mtry[1], 2), 
                           ', splitrule: ', rf_surv$tuning$best_tune$splitrule[1], 
                           ', n splits = ', rf_surv$tuning$best_tune$nsplit[1], 
                           ', node size = ', rf_surv$tuning$best_tune$nodesize[1]), 
         x = 'mtry', 
         y = 'C-index, CV')
  
  ## for the GBM models
  
  os_plots$tuning$gbm <- gbm_os$tuning$summary %>% 
    mutate(best = ifelse(cv_deviance == min(cv_deviance), 
                         'yes', 'no'))
  
  os_plots$tuning$gbm <- os_plots$tuning$gbm %>% 
    ggplot(aes(x = shrinkage, 
               y = cv_deviance, 
               color = factor(n.trees), 
               fill = best)) + 
    geom_line(aes(group = factor(n.trees))) + 
    geom_point(shape = 21, 
               size = 1, 
               color = 'black') + 
    scale_fill_manual(values = c(no = 'steelblue', 
                                 yes = 'coral3'), 
                      name = 'best tune') + 
    facet_grid(interaction.depth ~ n.minobsinnode) + 
    globals$common_theme + 
    labs(title = globals$algo_labels["gbm"], 
         subtitle = paste0('Best tune: ', 
                           'N trees = ', gbm_os$tuning$best_tune$n.trees[1], 
                           ', shrinkage = ', gbm_os$tuning$best_tune$shrinkage[1], 
                           ', interaction depth = ', gbm_os$tuning$best_tune$interaction.depth[1], 
                           ', minimal node size = ', gbm_os$tuning$best_tune$n.minobsinnode[1]), 
         x = 'N trees', 
         y = 'Deviance, CV')

# Performance stats for the algorithms: C-index and IBS, expression models ------
  
  insert_msg('Performance of the algorithms: C-index and IBS, expression models')
  
  ## comparison of the cohorts
  
  os_plots$stat_plots <- 
    list(stats = os_summary$stats, 
         plot_subtitle = globals$algo_labels[names(os_summary$stats)] %>% 
           paste('algorithm')) %>% 
    pmap(plot_surv_stats)
  
  ## comparison of the algorithms within single cohorts, expression-only models
  
  os_plots$algorithm_stat_plots <- 
    os_summary$stats%>% 
    compress(names_to = 'algorithm') %>% 
    mutate(cohort = factor(cohort, os_summary$stats[[1]]$cohort)) %>% 
    blast(cohort) %>% 
    list(stats = ., 
         plot_title = os_plots$plot_titles) %>% 
    pmap(plot_surv_stats, 
         palette = set_names(globals$algo_colors, 
                             unname(globals$algo_labels)), 
         labels = globals$algo_labels, 
         label_variable = 'algorithm', 
         color_variable = 'algorithm', 
         plot_subtitle = NULL) 
  
# Brier scores for the unique time points -------
  
  insert_msg('Brier score for the unique time points')
  
  for(i in names(os_summary$brier_scores)) {
    
    os_plots$bs_plots[[i]] <- os_summary$brier_scores[[i]] %>% 
      map(plot, cust_theme = globals$common_theme)
    
    os_plots$bs_plots[[i]] <- 
      list(x = os_plots$bs_plots[[i]] , 
           y = os_plots$plot_titles) %>% 
      pmap(function(x, y) x + 
             labs(title = y, 
                  subtitle = paste(globals$algo_labels[[i]], 'algorithm'), 
                  x = 'overall survival, months') + 
             scale_color_manual(values = c(reference = 'gray60', 
                                           training = 'coral3'), 
                                labels = c(reference = 'random prediction', 
                                           training = 'model'), 
                                name = ''))
    
  }

# Survival in the score tertiles -------
  
  insert_msg('Kaplan-Meier plots for survival in the score tertiles')
  
  for(i in names(os_summary$tertile_fits)) {
    
    os_plots$km_plots[[i]] <- 
      list(fit = os_summary$tertile_fits[[i]], 
           n_numbers = os_summary$tertile_n[[i]], 
           p_value = os_summary$tertile_test[[i]]$significance, 
           plot_title = os_plots$plot_titles) %>% 
      pmap(plot_tertile_km, 
           x_lab = 'overall survival months')
    
  }
  
# Variable importance -------
  
  insert_msg('Variable importance')
  
  os_plots$importance_plots <- 
    list(data = os_summary$importance, 
         imp_stat = c(rep('coef', 3), 
                      rep('delta_c_index', 2), 
                      'rel_influence'), 
         labeller = c(rep(list(c(positive = 'top unfavorable', 
                                 negative = 'top favorable')), 
                          3), 
                      rep(list(c(positive = 'top improvement', 
                                 negative = 'top worsening')), 
                          3)), 
         plot_title = paste0(globals$algo_labels[names(os_summary$importance)], 
                             ', variable importance'), 
         x_lab = c(list(expression('log HR'[Ridge]), 
                        expression('log HR'[ElasticNet]), 
                        expression('log HR'[LASSO])), 
                   rep(list(expression(Delta * ' C-index')), 2), 
                   list(expression(Delta * ' SSE'))), 
         n_top = c(rep(20, 5), 100), 
         palette = c(rep(list(c(positive = 'firebrick', 
                                negative = 'steelblue', 
                                ns = 'gray70')), 
                         5), 
                     list(c(positive = 'steelblue', 
                            negative = 'steelblue', 
                            ns = 'gray70')))) %>% 
    pmap(plot_surv_importance, 
         form = 'bar') %>% 
    map(~.x + 
          geom_vline(xintercept = 0,
                     linetype = 'dashed') + 
          theme(legend.position = 'none'))
  
# Error structure: plots of square errors of prediction -------
  
  insert_msg('Plots of square errors of prediction')
  
  ## the square errors are available for Cox-like models
  
  os_plots$square_plots <- list(ridge = ridge_os, 
                                elnet = elnet_os, 
                                lasso = lasso_os, 
                                gbm = gbm_os) %>% 
    map(~.x$calibration) %>% 
    future_map(map, 
               plot, 
               type = 'squares', 
               .options = furrr_options(seed = TRUE))
  
  ## plot titles and styling
  ## transposition in a more handy format
  
  for(i in names(os_plots$square_plots)) {
    
    os_plots$square_plots[[i]] <- 
      list(x = os_plots$square_plots[[i]], 
           y = paste(globals$algo_labels[[i]], 
                     os_plots$plot_titles, 
                     sep = ', ')) %>% 
      pmap(function(x, y) x %>% 
             map(~.x + 
                   labs(title = y) + 
                   globals$common_theme + 
                   theme(panel.grid.major.x = element_blank())))
    
  }
  
  os_plots$square_plots <- os_plots$square_plots %>% 
    map(transpose) %>%
    transpose
  
  for(i in names(os_plots$square_plots$time)) {
    
    os_plots$square_plots$time[[i]] <- 
      os_plots$square_plots$time[[i]] %>% 
      map(~.x + 
            labs(x = 'overall survival, months'))
    
  }
  
  os_plots$square_plots$observation <- 
    os_plots$square_plots$observation %>% 
    map(map, 
        ~.x + 
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank()) + 
          geom_hline(yintercept = 0.25, linetype = 'dashed'))
  
# END -------
  
  plan('sequential')
  
  rm(i)
  
  insert_tail()