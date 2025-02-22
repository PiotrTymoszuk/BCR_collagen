# correlation analysis for the pre-operative PSA antigen in blood and 
# expression levels of the PSA-coding gene KLK3 in the cancer tissue 
# (ComBat-adjusted log2 expression). The analysis is performed for the 
# pooled GEO cohort, TCGA, and the DKFZ collective, i.e. the data sets used 
# later for BCR-free survival modeling.
#
# Statistical significance is assessed by Spearman's permutation test.

  insert_head()
  
# container -------
  
  psa <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals --------
  
  insert_msg('analysis globals')
  
  ## ComBat-adjusted KLK3 expression and PSA concentrations
  
  psa$data <- list(geo_pool = geo_pool, 
                   tcga = tcga, 
                   dkfz = dkfz) %>% 
    map(~.x[c('clinic', 'combat')])
  
  for(i in names(psa$data)) {
    
    psa$data[[i]] <- 
      inner_join(psa$data[[i]]$clinic[c('sample_id', 
                                        'rfs_months', 
                                        'relapse', 
                                        'psa_diagnosis')], 
                 psa$data[[i]]$combat[c('sample_id', 
                                        'KLK3')], 
                 by = 'sample_id')
    
  }
  
  ## data sets for the correlation and BCR-free survival 
  ## analyses
  
  psa[c('cor_data', 'surv_data')] <- 
    list(c('KLK3', 'psa_diagnosis'), 
         c('rfs_months', 'relapse', 'KLK3')) %>% 
    map(~map(psa$data, select, all_of(.x))) %>% 
    map(map, ~filter(.x, complete.cases(.x)))
  
  ## minimal strata size for the survival analyis
  
  psa$min_n <-  map_dbl(psa$surv_data, nrow) * 0.25
  
# Correlation analysis -----
  
  insert_msg('Correlation analysis')
  
  psa$cor_test <- psa$cor_data %>% 
    map(~f_cor_test(.x$KLK3, 
                    .x$psa_diagnosis, 
                    type = 'permutation', 
                    method = 'spearman', 
                    as_data_frame = TRUE, 
                    n_iter = 1000)) %>% 
    compress(names_to = 'cohort') %>% 
    as_tibble
  
  ## formatting the testing results: 
  ## plot captions
  
  psa$cor_test <- psa$cor_test %>% 
    re_adjust %>% 
    mutate(plot_cap = paste0('n = ', n, 
                             ', \u03C1 = ', signif(rho, 2),
                             ', p = ', signif(p_adjusted, 2)))

# plots for the correlation analysis  -------
  
  insert_msg('Plots')
  
  psa$cor_plots <- 
    list(x = psa$cor_data, 
         y = globals$study_labels[names(psa$data)], 
         z = psa$cor_test$plot_cap) %>% 
    pmap(function(x, y, z) x %>% 
           plot_correlation(variables = c('KLK3', 'psa_diagnosis'), 
                            type = 'correlation', 
                            cust_theme = globals$common_theme + 
                              theme(axis.title.y = element_markdown()), 
                            plot_title = y, 
                            plot_subtitle = z, 
                            x_lab = 'PSA, ng/mL', 
                            y_lab = paste0(html_italic('KLK3'), 
                                          ', log<sub>2</sub> expression'))) %>% 
    map(~.x + 
          scale_y_continuous(trans = 'log', 
                             labels = function(x) signif(x, 2)) + 
          scale_x_continuous(trans = 'log', 
                             labels = function(x) signif(x, 2))) %>% 
    set_names(names(psa$data))
  
# Association of the KLK3 expression with BCR: optimal cutoffs ----
  
  insert_msg('Association with BCR-free survival')
  
  ## cut objects with the cutoffs corresponding to the largest differences
  ## in BCR-free survival
  
  psa$cut_obj <- list(data = psa$surv_data, 
                      min_n = psa$min_n) %>% 
    future_pmap(find_cutoff, 
                time = 'rfs_months', 
                event = 'relapse', 
                variable = 'KLK3',
                .options = furrr_options(seed = TRUE))
  
  ## testing results
  
  psa$surv_test <- psa$cut_obj %>% 
    map(summary) %>%
    compress(names_to = 'cohort')
  
  ## N numbers of patients and BCR events 
  ## in the expression strata
  
  psa$surv_n_numbers <- psa$cut_obj %>% 
    map(model.frame) %>% 
    map(blast, KLK3_strata) %>% 
    map(map, ~paste0('total: n = ', nrow(.x), 
                     '\nevents: n = ', sum(.x$relapse))) %>% 
    map(~map2(.x, names(.x), 
              function(str, nam) paste(nam, str, sep = '\n')))
  
# Kaplan-Meier plots -------
  
  insert_msg('Kplan-Meier plots')
  
  psa$surv_plots <- psa$cut_obj %>% 
    map(plot) %>% 
    map(~.x$plot)
  
  ## styling of the plots
  
  psa$surv_plots <- 
    list(x = psa$surv_plots, 
         y = psa$surv_n_numbers, 
         z = globals$study_labels[names(psa$surv_plots)]) %>% 
    pmap(function(x, y, z) x + 
           scale_color_manual(values = c('steelblue', 'firebrick'), 
                              labels = y, 
                              name = html_italic('KLK3')) + 
           labs(title = z, 
                x = 'BCR-free survival, months') + 
           globals$common_theme +
           theme(legend.title = element_markdown(), 
                 plot.tag = element_blank()))
  
# END ------
  
  psa$data <- NULL
  psa$cor_data <- NULL
  psa$surv_data <- NULL
  
  psa <- compact(psa)
  
  insert_tail()