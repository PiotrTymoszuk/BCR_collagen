# Viauslization of results of univariable Cox modeling

  insert_head()
  
# container ------
  
  uni_plots <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## the common significant variables
  
  uni_plots$common_significant <- uni_cut$common_significant %>% 
    reduce(union)
  
  ## C-indexes for the common significant variables
  
  uni_plots$data <- uni_cut$modeling %>% 
    map(select, 
        gene_symbol, 
        cutoff, 
        n_total_low, 
        n_total_high, 
        n_events_low, 
        n_events_high, 
        marker, 
        c_index, 
        ibs_model, 
        hr, hr_lower, hr_upper, 
        p_adjusted, 
        significance) %>% 
    compress(names_to = 'cohort') %>% 
    filter(gene_symbol %in% uni_plots$common_significant) %>% 
    mutate(cohort = factor(cohort, rev(names(uni_cut$modeling))), 
           hr_label = paste0(signif(hr, 2), 
                             ' [95%CI: ', signif(hr_lower, 2), 
                             ' to ', signif(hr_upper, 2), ']')) %>% 
    arrange(gene_symbol, cohort)
  
  ## plotting order: by C-index
  
  uni_plots$plot_order <- uni_plots$data %>% 
    mutate(sign = ifelse(marker == 'favorable', 1, -1), 
           c_index = sign * c_index) %>% 
    blast(gene_symbol) %>% 
    map_dbl(~mean(.x$c_index)) %>% 
    sort
  
  uni_plots$data <- uni_plots$data %>% 
    mutate(gene_symbol = factor(gene_symbol, names(uni_plots$plot_order)))
  
  ## subtitles for the KM plots with the optimal cutoff values
  ## and hazard ratios with 95% confidence intervals
  
  uni_plots$plot_subtitles <- uni_plots$data %>% 
    mutate(plot_subtitle = paste0('cutoff = ', signif(cutoff, 2), 
                                  ', HR = ', hr_label)) %>% 
    blast(cohort) %>% 
    map(blast, gene_symbol) %>% 
    map(map, ~.x$plot_subtitle)
  
  ## FDR-corrected p values to be shown in the plots
  
  uni_plots$p_values <- uni_plots$data %>% 
    blast(cohort) %>% 
    map(blast, gene_symbol) %>% 
    map(map, ~.x$significance)
  
  ## legend labels with with N numbers (total and BCR events)

  uni_plots$legend_labels <- uni_plots$data %>% 
    mutate(low_lab = paste0('low\ntotal: n = ', n_total_low, 
                            '\nevents: n = ', n_events_low), 
           high_lab = paste0('high\ntotal: n = ', n_total_high, 
                            '\nevents: n = ', n_events_high), 
           legend_labs = map2(low_lab, high_lab, c)) %>% 
    blast(cohort) %>% 
    map(blast, gene_symbol) %>% 
    map(map, ~.x$legend_labs)
    
# Volcano plots with results of univariable Cox modeling --------
  
  insert_msg('Volcano plots')
  
  uni_plots$volcano_plots <- 
    list(data = uni_cut$modeling %>% 
           map(mutate, beta = log(hr)), 
         plot_title = globals$study_labels[names(uni_cut$modeling)]) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'beta', 
         p_variable = 'p_adjusted', 
         regulation_level = 0, 
         x_lab = 'log HR, high vs low expressors', 
         y_lab = expression('-log'[10] * ' pFDR'), 
         top_significant = 20, 
         top_regulated = 20, 
         label_variable = 'gene_symbol', 
         label_type = 'text', 
         txt_face = 'italic', 
         txt_size = 2, 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          scale_fill_manual(values = c(upregulated = 'firebrick', 
                                       downregulated = 'steelblue', 
                                       ns = 'gray70'), 
                            labels = c(upregulated = 'unfavorable', 
                                       downregulated = 'favorable', 
                                       ns = 'ns'), 
                            name = 'survival marker')) %>% 
    map(tag2subtitle)
  
# Forest plot with HRs with 95% CI for the common significant factors --------
  
  insert_msg('Forest plot')
  
  ## commenting out the text labels: 
  ## HRs will be presented in a table
  
  uni_plots$forest_plot <- uni_plots$data %>% 
    ggplot(aes(x = log(hr), 
               y = cohort, 
               color = marker)) + 
    facet_grid(gene_symbol ~ .) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed') + 
    geom_errorbarh(aes(xmin = log(hr_upper), 
                       xmax = log(hr_lower)), 
                   height = 0, 
                   position = position_dodge(0.9)) +
    geom_point(shape = 16, 
               size = 2, 
               position = position_dodge(0.9)) + 
    #geom_text(aes(label = hr_label), 
     #         size = 2.3, 
      #        hjust = 0.5, 
       #       vjust = -0.8, 
        #      position = position_dodge(0.9), 
         #     show.legend = FALSE) + 
    scale_color_manual(values = c(unfavorable = 'firebrick', 
                                  favorable = 'steelblue', 
                                  ns = 'gray70'), 
                       name = 'association\nwith risk') + 
    scale_y_discrete(labels = globals$study_labels) + 
    globals$common_theme + 
    theme(strip.text.y = element_text(angle = 0, 
                                      hjust = 0, 
                                      face = 'italic')) + 
    labs(title = 'Collagen-related markers, mRNA', 
         subtitle = 'BCR-free survival, shared significant effects', 
         x = 'log HR \u00B1 95% CI, high vs low expressors', 
         y = 'Cohort')
  
# Dot plot of C-indexes ---------
  
  insert_msg('Dot plot of C-indexes')
  
  ## actually: C-index and IBS
  
  uni_plots$c_index_plot <- uni_plots$data %>% 
    ggplot(aes(x = c_index, 
               y = cohort, 
               fill = marker, 
               size = 1 - ibs_model)) + 
    facet_grid(gene_symbol ~ .) + 
    geom_vline(xintercept = 0.5, 
               linetype = 'dashed') + 
    geom_point(shape = 21) + 
    scale_y_discrete(labels = globals$study_labels) + 
    scale_radius(range = c(1, 3),
                 name = '1 - IBS') + 
    scale_fill_manual(values = c(unfavorable = 'firebrick', 
                                 favorable = 'steelblue', 
                                 ns = 'gray70'), 
                      name = 'association with risk') + 
    globals$common_theme + 
    theme(strip.text.y = element_text(angle = 0, 
                                      hjust = 0, 
                                      face = 'italic')) + 
    labs(title = 'Collagen-related markers, mRNA', 
         subtitle = 'BCR-free survival, shared significant effects', 
         x = 'C-index', 
         y = 'Cohort')
  
# Kaplan-Meier plots --------
  
  insert_msg('Kaplan-Meier plots')
  
  ## for the common significant effects 
  
  uni_plots$km_plots <- uni_cut$cutoff_obj %>% 
    map(~.x[uni_plots$common_significant]) %>% 
    map(map, plot) %>% 
    map(map, ~.x$plot)
  
  ## styling
  
  for(i in names(uni_plots$km_plots)) {
    
    ## titles and subtitles, common plot theme, 
    ## displaying p values
    
    uni_plots$km_plots[[i]] <- 
      list(x = uni_plots$km_plots[[i]], 
           y = names(uni_plots$km_plots[[i]]) %>% 
             html_italic %>% 
             paste(globals$study_labels[[i]], sep = ', '), 
           z = uni_plots$plot_subtitles[[i]], 
           v = uni_plots$p_values[[i]]) %>% 
      pmap(function(x, y, z, v) x + 
             labs(title = y, 
                  subtitle = z, 
                  x = 'BCR-free survival, months') + 
             globals$common_theme + 
             theme(plot.title = element_markdown(), 
                   plot.tag = element_blank()) + 
             annotate('text',
                      label = v, 
                      size = 2.75, 
                      x = 0, 
                      y = 0, 
                      hjust = 0, 
                      vjust = 0))
    
    ## color scale with legends containing N number information
    
    uni_plots$km_plots[[i]] <- 
      list(x = uni_plots$km_plots[[i]], 
           y = uni_plots$legend_labels[[i]]) %>% 
      pmap(function(x, y) x +
             scale_color_manual(values = c(low = 'steelblue', 
                                           high = 'firebrick'), 
                                labels = y[[1]], 
                                name = ''))
    
  }
  
  uni_plots$km_plots <- transpose(uni_plots$km_plots)
  
# END ------
  
  uni_plots <- 
    uni_plots[c("volcano_plots", "forest_plot", "c_index_plot", "km_plots")]
  
  insert_tail()