# Normalization and batch effect removal with ComBat at the whole transcriptome 
# level: done with the toolset of htGLMNET

  insert_head()

# container ------
  
  combat <- list()
  
# common expression variables and raw expression -------
  
  insert_msg('Common expression variables and raw expression')
  
  ## raw expression
  
  combat$raw_expression <- globals$study_exprs %>% 
    eval %>% 
    map(~.x$expression) %>% 
    map(function(x) if('tissue_type' %in% names(x)) filter(x, tissue_type == 'tumor') else x) %>% 
    compact %>% 
    map(column_to_rownames, 'sample_id') %>% 
    map(~.x[!names(.x) %in% c('patient_id', 'tissue_type')])
  
  ## common genes and collagen-related genes
  
  combat$cmm_genes <- combat$raw_expression %>% 
    map(names) %>% 
    reduce(intersect)
  
  combat$col_genes <- globals$genes
  
  ## matrices in genes in the rows
  
  combat$raw_expression <- combat$raw_expression %>% 
    map(t)

# Batch effect adjustment -------
  
  insert_msg('Batch effect adjustment')
  
  combat$multi_obj <- 
    multi_process(train = combat$raw_expression$gse16560, 
                  test = combat$raw_expression[-1])
  
  ## retrieval of the corrected data sets
  
  combat$expression <- 
    c(list(gse16560 = combat$multi_obj$train), 
      combat$multi_obj$test) %>% 
    map(t) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'sample_id') %>% 
    map(as_tibble)
  
# Pooled data set for the analyses ------
  
  insert_msg('Pooled data set for the analyses')
  
  combat$test_data <- combat$expression %>% 
    map(select, 
        sample_id, all_of(combat$col_genes)) %>% 
    compress(names_to = 'study') %>% 
    mutate(study = globals$study_labels[study], 
           study = unname(study), 
           study = factor(study, 
                          globals$study_labels[names(combat$expression)]))

# Diagnostic stats and tests for the collagen-related genes -------
  
  insert_msg('Diagnostic stats')

  ## descriptive stats
  
  combat$stats <- combat$test_data %>% 
    fast_num_stats(variables = combat$col_genes, 
                   split_fct = 'study')

  ## comparison of the cohorts with one-way ANOVA
  ## and eta-square effect size stat
  
  combat$test <- 
    f_one_anova(combat$test_data[combat$col_genes], 
                f = combat$test_data$study, 
                as_data_frame = TRUE, adj_method = 'BH') %>% 
    re_adjust(method = 'BH') %>% 
    mutate(eff_size = paste('\u03B7\u00B2 =', signif(etasq)), 
           plot_cap = paste(eff_size, significance, sep = ', ')) %>% 
    as_tibble

  ## result table
  
  combat$result_tbl <- 
    left_join(combat$stats, 
              combat$test[c('variable', 'significance', 'eff_size')], 
              by = 'variable') %>% 
    format_summ_tbl
  
# Diagnostic plots -------
  
  insert_msg('Diagnostic plots')
  
  plan('multisession')
  
  combat$plots <- 
    list(variable = combat$col_genes, 
         plot_title = combat$col_genes, 
         plot_subtitle = combat$test$plot_cap) %>% 
    future_pmap(plot_variable, 
                combat$test_data, 
                split_factor = 'study', 
                type = 'violin', 
                cust_theme = globals$common_theme,
                y_lab = expression('normalized log'[2] * ' expression'), 
                x_n_labs = TRUE, 
                .options = furrr_options(seed = TRUE)) %>% 
    map(~.x + 
          theme(plot.title = element_text(face = 'bold.italic'))) %>% 
    set_names(combat$col_genes)
  
  plan('sequential')

# Caching the results ------
  
  insert_msg('Caching the results')
  
  combat <- 
    combat[c("cmm_genes", 
             "expression", 
             "stats", 
             "test", 
             "result_tbl", 
             "plots")]
  
  save(combat, file = './data/combat.RData')
  
# END -----

  insert_tail()