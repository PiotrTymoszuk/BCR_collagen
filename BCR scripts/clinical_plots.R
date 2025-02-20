# Graphical comparison of performance of two GBM models: 
# 1) with the clinical predictors of BCR-free survival
# 2) with the clinical and transcriptional predictors of BCR-free survival 

  insert_head()
  
# container -----
  
  surv_clplots <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## plotting labels and colors 
  
  surv_clplots$algo_labels <- 
    c(reference = 'random prediction', 
      gbm = paste0(globals$algo_labels["gbm"], ', expression'), 
      globals$algo_labels[c("gbm_clinic", "gbm_combi")])
  
  surv_clplots$algo_colors <- 
    c(reference = 'gray60', 
      globals$algo_colors[c("gbm", "gbm_clinic", "gbm_combi")])
  
  ## global performance stats 
  
  surv_clplots$stats <- 
    surv_summary$stats[c("gbm_clinic", "gbm", "gbm_combi")] %>% 
    compress(names_to = 'algorithm') %>% 
    mutate(cohort_lab = globals$study_labels[cohort], 
           cohort_lab = ifelse(algorithm == 'gbm_clinic', 
                               cohort_lab, NA), 
           algorithm = factor(algorithm, 
                              c('gbm_clinic', 'gbm', 'gbm_combi')))
  
  ## Brier scores for the unique survival time points

  surv_clplots$brier_scores <- 
    surv_summary$brier_scores[c("gbm_clinic", "gbm", "gbm_combi")] %>% 
    map2(., 
         list(c('time', 'reference', 'training'), 
              c('time', 'training'), 
              c('time', 'training')), 
         function(x, y) x %>% 
           map(select, all_of(y))) %>% 
    transpose %>% 
    map(reduce, full_join, by = c('time')) %>% 
    map(set_names, c('time', 'reference', 'gbm_clinic', 'gbm', 'gbm_combi')) %>% 
    map(pivot_longer, 
        cols = c(reference, gbm_clinic, gbm, gbm_combi), 
        names_to = 'algorithm', 
        values_to = 'brier_score') %>% 
    map(mutate, 
        algorithm = factor(algorithm, 
                           c('reference', 
                             'gbm_clinic', 
                             'gbm', 
                             'gbm_combi')))
  
# Plots of the global performance stats -------
  
  insert_msg('Plots of the globals performance stats')
  
  surv_clplots$stat_plot <- surv_clplots$stats %>% 
    ggplot(aes(x = c_index, 
               y = 1 - ibs_model, 
               fill = algorithm)) + 
    geom_vline(xintercept = 0.5, 
               linetype = 'dashed') + 
    geom_hline(yintercept = 0.75, 
               linetype = 'dashed') + 
    geom_line(aes(group = cohort)) + 
    geom_point(shape = 21, 
               size = 2, 
               color = 'gray30') + 
    geom_text_repel(aes(label = cohort_lab), 
                    size = 2.5) + 
    scale_fill_manual(values = surv_clplots$algo_colors, 
                      labels = surv_clplots$algo_labels, 
                      name = '') + 
    globals$common_theme + 
    labs(title = 'GBM model performance', 
         x = 'C-index', 
         y = '1 - IBS')
  
# Plots of Brier scores for the unique time points ------
  
  insert_msg('Plots of Brier scores for the unique time points')
  
  surv_clplots$brier_plots <- 
    list(x = surv_clplots$brier_scores, 
         y = globals$study_labels[names(surv_clplots$brier_scores)]) %>% 
    pmap(function(x, y) x %>% 
           filter(complete.cases(.)) %>% 
           ggplot(aes(x = time, 
                      y = brier_score, 
                      color = algorithm)) + 
           geom_line() + 
           scale_color_manual(values = surv_clplots$algo_colors, 
                              labels = surv_clplots$algo_labels,
                              name = '') + 
           globals$common_theme + 
           labs(title = y, 
                x = 'min/max scaled BCR-free survival', 
                y = 'Brier score'))
  
# END ------
  
  surv_clplots$stats <- NULL
  surv_clplots$brier_scores <- NULL
  
  surv_clplots <- compact(surv_clplots)
  
  insert_tail()
    
    