# Assessment of clustering tendency of the collagen data sets

  insert_head()
  
# container -------
  
  clust_tend <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# expression data frames --------
  
  insert_msg('Expression data frames')
  
  clust_tend$expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(filter, tissue_type == 'tumor') %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(globals$genes_interest$gene_symbol)) %>% 
    map(center_data, 'mean')
  
# Numbers of observations -------
  
  insert_msg('Numbers of observations')
  
  ## total observation numbers
  
  clust_tend$n_numbers <- clust_tend$expression %>% 
    map_dbl(nrow) %>% 
    paste('n =', .)
  
  ## numbers of observations used for computation of the Hopkins stats
  ## (25% of total)
  
  clust_tend$n_sample <- clust_tend$expression %>% 
    map_dbl(nrow)
  
  clust_tend$n_sample <- ceiling(0.25 * clust_tend$n_sample)
  
# Hopkins stats and heat maps ------
  
  insert_msg('Hopkins stats and heat maps')
  
  clust_tend$hopkins <- 
    list(data = clust_tend$expression, 
         n = clust_tend$n_sample) %>% 
    future_pmap(get_clust_tendency, 
                .options = furrr_options(seed = TRUE))
  
  ## extracting the numeric stats
  
  clust_tend$stats <- clust_tend$hopkins %>% 
    map(~.x[c('hopkins_stat', 'p_value')]) %>% 
    map(as_tibble) %>% 
    compress(names_to = 'cohort')
  
  ## extracting the visualizations and appending with n numbers
  ## and Hopkins stats
  
  clust_tend$heat_maps <- clust_tend$hopkins %>% 
    map(~.x$plot)
  
  clust_tend$heat_maps <- 
    list(x = clust_tend$heat_maps, 
         y = globals$study_labels[names(clust_tend$heat_maps)], 
         v = signif(clust_tend$stats$hopkins_stat, 2), 
         z = clust_tend$n_numbers) %>% 
    pmap(function(x, y, v, z) x + 
           scale_fill_gradient2(low = 'firebrick', 
                                mid = 'white', 
                                high = 'steelblue', 
                                midpoint = 10, 
                                limits = c(0, 20), 
                                oob = scales::squish) + 
           labs(title = y, 
                subtitle = paste0('H = ', v, ', ', z), 
                x = 'Cancer sample', 
                y = 'Cancer sample',
                fill = 'Euclidean distance') + 
           globals$common_theme + 
           theme(axis.text.x = element_blank(), 
                 axis.text.y = element_blank(), 
                 axis.line = element_blank(), 
                 axis.ticks = element_blank(), 
                 axis.title = globals$common_text))
  
# UMAP -------
  
  insert_msg('UMAP')
  
  ## with cosine distance
  
  clust_tend$umap_objects <- clust_tend$expression %>% 
    map(reduce_data, 
        distance_method = 'cosine', 
        kdim = 2, 
        red_fun = 'umap')
  
# UMAP plots -------
  
  insert_msg('UMAP layout plots')
  
  clust_tend$umap_plots <- 
    list(x = clust_tend$umap_objects, 
         point_color = globals$study_colors[names(clust_tend$umap_objects)]) %>% 
    pmap(plot,
         type = 'scores', 
         cust_theme = globals$common_theme)
  
  ## plot titles and numeric stats
  
  clust_tend$umap_plots <- 
    list(x = clust_tend$umap_plots, 
         y = globals$study_labels[names(clust_tend$heat_maps)], 
         v = signif(clust_tend$stats$hopkins_stat, 2), 
         z = clust_tend$n_numbers) %>% 
    pmap(function(x, y, v, z) x + 
           labs(title = y, 
                subtitle = paste0('H = ', v, ', ', z)))

# END -----
  
  clust_tend$expression <- NULL
  clust_tend$n_numbers <- NULL
  clust_tend$n_sample <- NULL
  clust_tend$hopkins <- NULL
  
  clust_tend <- compact(clust_tend)
  
  plan('sequential')
  
  insert_tail()