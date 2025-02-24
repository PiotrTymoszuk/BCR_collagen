# Supplementary tables for the revised manuscript 

  insert_head()
  
# container ------
  
  suppl_tab <- list()
  
# Characteristic of the cohorts: pooled GEO, TCGA, DKFZ, and GSE1650 ----
  
  insert_msg('Characteristic of the cohorts')
  
  ## Gleason scores are not shown any more, because they overlap 
  ## with ISUP
  
  suppl_tab$cohorts <- cohorts$stats %>% 
    filter(!stri_detect(Variable, fixed = 'Gleason')) %>% 
    select(Variable, `pooled GEO`, TCGA, DKFZ, GSE16560) %>% 
    as_mdtable(label = 'cohorts', 
               ref_name = 'cohorts', 
               caption = paste('Characteristic of the analyzed cohorts.', 
                               'Numeric variables are presented as medians with', 
                               'interquartile ranges (IQR) and ranges.', 
                               'Qualitative variables are presented as', 
                               'percentages of categories within the', 
                               'complete observation set.'))
  
# Characteristic of the pooled GEO cohort -------
  
  insert_msg('Characteristic of the pooled GEO cohort')
  
  ## Gleason scores are not shown; redundant with ISUP
  
  suppl_tab$pool <- expl_pool$result_tbl %>% 
    filter(!stri_detect(Variable, fixed = 'Gleason')) %>% 
    as_mdtable(label = 'pool', 
               ref_name = 'pool', 
               caption = paste('Characteristic of GEO data sets which constitute', 
                               'the pooled GEO cohort.', 
                               'Numeric variables are presented as medians with', 
                               'interquartile ranges (IQR) and ranges.', 
                               'Qualitative variables are presented as', 
                               'percentages of categories within the', 
                               'complete observation set.'))
  
# Characteristic of the training and test portions of the GSE1650 cohort ------
  
  insert_msg('Characteristic of the training and test portion of the cohort')
  
  ## Gleason scores treated as above
  
  suppl_tab$gse16560 <- expl_os$result_tbl %>% 
    filter(!stri_detect(Variable, fixed = 'Gleason')) %>% 
    as_mdtable(label = 'gse16560', 
               ref_name = 'gse16560', 
               caption = paste('Characteristic of the training and test subsets', 
                               'of the GSE16560 cohort.', 
                               'For modeling of overall survival, the GSE16560', 
                               'cohort was randomly split into a training and', 
                               'a test subset in a  2:1 ratio.', 
                               'Numeric characteristics of the subsets are', 
                               'presented as medians with', 
                               'interquartile ranges (IQR) and ranges.', 
                               'Qualitative variables are presented as', 
                               'percentages of categories within the', 
                               'complete observation set.'))
  
# The collagen-related genes ------
  
  insert_msg('The collagen-related genes')
  
  suppl_tab$genes <- globals$genes_lexicon %>% 
    mutate(entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              column = 'ENTREZID', 
                              keytype = 'SYMBOL'),
           gene_group = factor(gene_group, 
                               c('proline metabolism', 
                                 'collagen modification', 
                                 'ECM component', 
                                 'ECM processing', 
                                 'adhesion'))) %>% 
    arrange(gene_group) %>% 
    select(gene_group, gene_symbol, entrez_id) %>% 
    set_names(c('Functional classification', 'Gene symbol', 'Entrez ID')) %>% 
    mdtable(label = 'genes', 
            ref_name = 'genes', 
            caption = 'Collagen-related genes and their classification.')
  
# Co-expression networks, vertex importance stats --------
  
  insert_msg('Co-expression networks, gene importance')
  
  ## top 5 genes with the largest degree per cohort are shown, 
  ## the complete table as a supplementary file; 
  ## betweenness is min/max scaled
  
  suppl_tab$nets <- net$graph_obj[c("geo_pool", "tcga", "dkfz")] %>% 
    map(get_vertex_attributes) %>% 
    map(select, 
        name, gene_group, degree, hub_score, betweenness, transitivity) %>% 
    map(mutate, 
        betweenness = minMax(betweenness)) %>% 
    compress(names_to = 'cohort') %>% 
    relocate(cohort) %>% 
    mutate(cohort = factor(cohort, c("geo_pool", "tcga", "dkfz"))) %>% 
    arrange(cohort) %>% 
    mutate(cohort = globals$study_labels[as.character(cohort)], 
           betweenness = signif(betweenness, 2), 
           hub_score = signif(hub_score, 2), 
           transitivity = signif(transitivity, 2)) %>% 
    set_names(c('Cohort', 
                'Gene symbol', 'Functional classification', 
                'Degree', 'Hub score', 'Betweenness', 'Transitivity'))
  
  suppl_tab$nets <- suppl_tab$nets %>% 
    as_mdtable(label = 'nets', 
               ref_name = 'nets', 
               caption = paste('Vertex importance statistics for', 
                               'co-expression networks of the collagen-related', 
                               'transcripts.', 
                               'The co-expression networks were built for in the', 
                               'pooled GEO, TCGA, and DKFZ cohorts for', 
                               'transcripts with pairwise correlation of', 
                               'expression levels with', 
                               "Spearman's rho >= 0.5.", 
                               'Metrics of importance of the vertices of the', 
                               'co-expression networks, degree, hub score,', 
                               'betweenness, and transitivity, were computed.', 
                               'Top five vertices with the largest hub score per', 
                               'cohort are presented.', 
                               'The full table is available as a supplementary', 
                               'Excel file.'))
  
# Results of univariable survival modeling --------
  
  insert_msg('Univariable modeling of BCR-free survival')
  
  ## the common markers are shown, the rest is available 
  ## as a supplementary Excel file
  
  suppl_tab$uni_bcr <- uni_cut$modeling %>% 
    compress(names_to = 'cohort') %>% 
    transmute(Cohort = globals$study_labels[cohort], 
              Cohort = factor(Cohort, 
                              globals$study_labels[names(uni_cut$modeling)]), 
              `Gene symbol` = gene_symbol, 
              `Association with BCR risk` = marker, 
              `log2 expression cutoff` = signif(cutoff, 2), 
              `Low expressors, N` = paste0('total: n = ', n_total_low, 
                             '\nBCR: n = ', n_events_low), 
              `High expressors, N` = paste0('total: n = ', n_total_high, 
                              '\nBCR: n = ', n_events_low), 
              `HR, 95%CI` = paste0(signif(hr, 2), '\n[', 
                          signif(hr_upper, 2), ' to ', 
                          signif(hr_lower, 2), ']'), 
              `Significance` = significance, 
              `C-index` = signif(c_index, 2), 
              `IBS` = signif(ibs_model, 2))
  
  suppl_tab$uni_bcr <- suppl_tab$uni_bcr %>% 
    as_mdtable(label = 'uni_bcr', 
               ref_name = 'uni_bcr', 
               caption = paste('Univariable analysis of biochemical relapse', 
                               '(BCR) free survival for the collagen-realted', 
                               'transcripts in the pooled GEO, TCGA, and DKFZ', 
                               'cohort.', 
                               'PCA patients were classified as high and low', 
                               'expressors of the collagen-related transcripts', 
                               'by expression cutoffs corresponding to the', 
                               'largest differences in BCR-free survival', 
                               'between the expression strata.', 
                               'BCR risk of high as compared with low', 
                               'expressors was modeled by univariable Cox', 
                               'proportional hazard regression.', 
                               'P values were corrected for multiple testing', 
                               'with the false discovery rate method.', 
                               'The modeling results are shown for the', 
                               'collagen-related genes associated with', 
                               'significant differences in BCR-free survival', 
                               'in all three investigated cohorts.', 
                               'The full analysis is available as a', 
                               'supplementary Excel file.'))
  
# Tuning parameters, machine learning models ------
  
  insert_msg('Tuning parameters, machine learning models')
  
  ## regularized Cox modeling
  
  suppl_tab$tune[c('ridge', 'elnet', 'lasso')] <- 
    list(ridge_surv, elnet_surv, lasso_surv) %>% 
    map(~.x$lambda_tune$lambda) %>% 
    map(~tibble(criterion = 'minimal deviance,\nrepeated 10-fold cross-validation', 
                best_tune = paste('\u03BB =', signif(.x, 2))))
  
  ## SVM
  
  suppl_tab$tune$svm <- 
    tibble(criterion = 'maximal concordance index,\nrepeated 10-fold cross-validation', 
           best_tune = paste0('SVM model type = ', 
                              svm_surv$tuning$best_tune$type, 
                             '\nkernel = ', 
                             svm_surv$tuning$best_tune$kernel, 
                             '\n\u03B3 = ', 
                             svm_surv$tuning$best_tune$gamma.mu))
  
  ## Random Forest
  
  suppl_tab$tune$rf <- 
    tibble(criterion = 'maximal concordance index,\nout-of-bag predictions', 
           best_tune = paste0('number of variables per tree, mtry = ', 
                              rf_surv$tuning$best_tune$mtry, 
                              '\nsplitting rule = ', 
                              rf_surv$tuning$best_tune$splitrule,
                              '\nminimal node size = ', 
                              rf_surv$tuning$best_tune$nodesize, 
                              '\nnumber of splits = ', 
                              rf_surv$tuning$best_tune$nsplit))
  
  ## GBM models
  
  suppl_tab$tune[c('gbm', 'gbm_clinic', 'gbm_combi')] <- 
    list(gbm_surv, surv_clin, surv_combi) %>% 
    map(~.x$tuning$best_tune) %>% 
    map(~tibble(criterion = 'minimal deviance,\n10-fold cross-validation', 
                best_tune = paste0('number of decision trees = ', 
                                   .x$n.trees, 
                                   '\nshrinkage = ', 
                                   .x$shrinkage, 
                                   '\ninteraction depth = ', 
                                   .x$interaction.depth, 
                                   '\nminimal node size = ', 
                                   .x$n.minobsinnode)))
  
  ## the entire table
  
  suppl_tab$tune <- suppl_tab$tune %>% 
    compress(names_to = 'algorithm') %>% 
    transmute(`BCR model` = globals$algo_labels[algorithm], 
              `BCR model` = stri_replace(`BCR model`, 
                                         regex = '^GBM$', 
                                         replacement = 'GBM, expression'), 
              `Selection criterion` = criterion, 
              `Tuned parameters` = best_tune)
  
  suppl_tab$tune <- suppl_tab$tune %>% 
    as_mdtable(label = 'tune', 
               ref_name = 'tune', 
               caption =  paste('Selection of the optimal parameters of machine', 
                                'learning models of biochemical', 
                                'relapse (BCR) free survival in the pooled GEO', 
                                'training cohort.'))
  
# Evaluation of the BCR-free survival models -------
  
  insert_msg('Evaluation of the BCR-free survival models')
  
  suppl_tab$bcr_eval <- surv_summary$stats %>% 
    compress(names_to = 'algorithm') %>% 
    transmute(`BCR model` = globals$algo_labels[algorithm], 
              `BCR model` = stri_replace(`BCR model`, 
                                         regex = '^GBM$', 
                                         replacement = 'GBM, expression'), 
              `Cohort type` = dataset, 
              Cohort = globals$study_labels[cohort], 
              `C-index` = signif(c_index, 2), 
              IBS = signif(ibs_model, 2))
  
  suppl_tab$bcr_eval <- suppl_tab$bcr_eval %>% 
    as_mdtable(label = 'bcr_eval', 
               ref_name = 'bcr_eval', 
               caption = paste('Performance of machine learning models at', 
                               'prediction of biochemical relapse (BCR) free', 
                               'survival.'))
  
# tuning parameters for the models of overall survival --------
  
  insert_msg('Tuning parameters for the models of overall survival')
  
  ## regularized Cox modeling
  
  suppl_tab$os_tune[c('ridge', 'elnet', 'lasso')] <- 
    list(ridge_os, elnet_os, lasso_os) %>% 
    map(~.x$lambda_tune$lambda) %>% 
    map(~tibble(criterion = 'minimal deviance,\nrepeated 10-fold cross-validation', 
                best_tune = paste('\u03BB =', signif(.x, 2))))
  
  ## SVM
  
  suppl_tab$os_tune$svm <- 
    tibble(criterion = 'maximal concordance index,\nrepeated 10-fold cross-validation', 
           best_tune = paste0('SVM model type = ', 
                              svm_os$tuning$best_tune$type, 
                              '\nkernel = ', 
                              svm_os$tuning$best_tune$kernel, 
                              '\n\u03B3 = ', 
                              svm_os$tuning$best_tune$gamma.mu))
  
  ## Random Forest
  
  suppl_tab$os_tune$rf <- 
    tibble(criterion = 'maximal concordance index,\nout-of-bag predictions', 
           best_tune = paste0('number of variables per tree, mtry = ', 
                              rf_os$tuning$best_tune$mtry, 
                              '\nsplitting rule = ', 
                              rf_os$tuning$best_tune$splitrule,
                              '\nminimal node size = ', 
                              rf_os$tuning$best_tune$nodesize, 
                              '\nnumber of splits = ', 
                              rf_os$tuning$best_tune$nsplit))
  
  ## GBM models
  
  suppl_tab$os_tune$gbm <- 
    tibble(criterion = 'minimal deviance,\n10-fold cross-validation', 
           best_tune = paste0('number of decision trees = ', 
                              gbm_os$tuning$best_tune$n.trees, 
                              '\nshrinkage = ', 
                              gbm_os$tuning$best_tune$shrinkage, 
                              '\ninteraction depth = ', 
                              gbm_os$tuning$best_tune$interaction.depth, 
                              '\nminimal node size = ', 
                              gbm_os$tuning$best_tune$n.minobsinnode))

  ## the entire table
  
  suppl_tab$os_tune <- suppl_tab$os_tune %>% 
    compress(names_to = 'algorithm') %>% 
    transmute(`OS model` = globals$algo_labels[algorithm], 
              `Selection criterion` = criterion, 
              `Tuned parameters` = best_tune)
  
  suppl_tab$os_tune <- suppl_tab$os_tune %>% 
    as_mdtable(label = 'os_tune', 
               ref_name = 'os_tune', 
               caption =  paste('Selection of the optimal parameters of machine', 
                                'learning models of overall survival (OS)', 
                                'in the training subset of the GSE16560 cohort.', 
                                'The models used log2 expression of the', 
                                'collagen-related genes as explanatory', 
                                'variables.'))
  
# Evaluation of the models of overall survival in the GSE16560 cohort -------
  
  insert_msg('Evaluation of performance of the OS models in the GSE16560 cohort')
  
  suppl_tab$os_eval <- os_summary$stats %>% 
    compress(names_to = 'algorithm') %>% 
    transmute(`OS model` = globals$algo_labels[algorithm], 
              `Cohort subset` = cohort, 
              `C-index` = signif(c_index, 2), 
              IBS = signif(ibs_model, 2))
  
  suppl_tab$os_eval <- suppl_tab$os_eval %>% 
    as_mdtable(label = 'os_eval', 
               ref_name = 'os_eval', 
               caption = paste('Performance of machine learning models at', 
                               'prediction of overall survival (OS) in the', 
                               'training and test subsets of the GSE16560', 
                               'cohort.', 
                               'The models used log2 expression of the', 
                               'collagen-related genes as explanatory factors.'))
  
# Saving the tables on the disc -----
  
  insert_msg('Saving the tables on the disc')
  
  suppl_tab %>% 
    save_excel(path = './paper/supplementary_tables.xlsx', 
               prefix = 'Supplementary Table S')
  
# END -------
  
  insert_tail()