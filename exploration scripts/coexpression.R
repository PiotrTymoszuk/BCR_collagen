# Co-regulation of collagen gene expression in the tumor tissue
# Euclidean distances and MDS, as well as by mean-centered PCA. 
# done with Z-scores of ComBat-adjusted log2 expression levels

  insert_head()
  
# container -----
  
  coex <- list()
  
# analysis data frames ------
  
  insert_msg('Analysis data frames')
  
  ## variables of interest
  
  coex$variables <- globals$genes
  
  ## ComBat 
  
  coex$data <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$combat) %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(select, all_of(coex$variables)) %>% 
    map(zScores)
  
# Numbers of observations -------
  
  insert_msg('Numbers of observations')
  
  coex$n_numbers <- coex$data %>% 
    map_dbl(nrow) %>% 
    paste('n =', .)
  
# MDS -----
  
  insert_msg('MDS')
  
  coex$mds_objects <- coex$data %>% 
    map(t) %>% 
    map(reduce_data, 
        kdim = 2, 
        red_fun = 'mds')
  
# Plots ---------
  
  insert_msg('Plots')
  
  coex$mds_plots <- 
    list(x = coex$mds_objects, 
         plot_subtitle =  coex$n_numbers) %>% 
    pmap(plot, 
         cust_theme = globals$common_theme)
  
  ## appending with titles and gene symbols
  
  coex$mds_plots <- 
    list(x = coex$mds_plots, 
         y = globals$study_labels[names(coex$mds_plots)]) %>% 
    pmap(function(x, y) x + 
           geom_text_repel(aes(label = observation), 
                           size = 2.5,
                           fontface = 'italic') + 
           theme(plot.tag = element_blank()) + 
           labs(title = y))
  
# PCA --------
  
  insert_msg('PCA')
  
  ## as inferred from scree plots, the first 4 dimensions
  ## capture 3/4 of the variance; the first two are clearly 
  ## dominant
  
  coex$pca_objects <- coex$data %>% 
    map(reduce_data, 
        kdim = 6, 
        red_fun = 'pca')
  
# PCA plots --------
  
  insert_msg('PCA plots')
  
  ## scree plots of the dimensions' variances
  
  coex$pca_scree_plots <- coex$pca_objects %>% 
    map(plot, 
        type = 'scree', 
        cust_theme = globals$common_theme) %>% 
    map2(., globals$study_labels[names(coex$pca_objects)], 
         ~.x + 
           theme(plot.tag = element_blank()) + 
           labs(title = .y))
  
  ## loadings
  
  coex$pca_loadings_plots <- 
    list(x = coex$pca_objects, 
         plot_subtitle =  coex$n_numbers) %>% 
    pmap(plot, 
         type = 'loadings', 
         label_points = FALSE, 
         cust_theme = globals$common_theme)
  
  ## adding cohort names in the titles and labeling points 
  ## with gene symbols
  
  coex$pca_loadings_plots <- 
    list(x = coex$pca_loadings_plots, 
         y = globals$study_labels[names(coex$pca_loadings_plots)]) %>% 
    pmap(function(x, y, z) x + 
           geom_text_repel(aes(label = variable), 
                           size = 2.5, 
                           fontface = 'italic') + 
           labs(title = y) + 
           theme(plot.tag = element_blank()))
  
# END ------
  
  coex$data <- NULL
  coex$n_numbers <- NULL
  
  coex <- compact(coex)
  
  insert_tail()