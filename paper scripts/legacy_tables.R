# Supplementary tables for the manuscript. 
# As discussed with the study team, we're not showing the GSE16560 cohort, as 
# it does not fulfill the strict selection criteria for the studies (biochemical 
# relapse and relapse-free survival)

  insert_head()
  
# containers ------
  
  suppl_paper_tbl <- list()
  
# globals -------
  
  insert_msg('Globals')
  
  ## metabolic reaction lexicon

  suppl_paper_tbl$reactions <- ana_meta$models$tcga %>% 
    components('regulation') %>% 
    select(react_id) %>% 
    mutate(react_name = annotate_bigg(react_id, annotation_db = Recon2D), 
           react_name = unlist(react_name))
  
# Table S1: cohort characteristic --------
  
  insert_msg('Table S1: cohort characteristic')
  
  suppl_paper_tbl$cohorts <- cohorts$stats %>% 
    select(-GSE16560) %>% 
    set_names(c('Variable', 
                globals$study_labels[names(globals$study_labels) != 'gse16560'])) %>% 
    mdtable(label = 'table_S1_study_cohorts', 
            ref_name = 'cohorts', 
            caption = paste('Characteristic of the analyzed cohorts.', 
                            'Numeric variables are presented as medians with', 
                            'interquartile ranges (IQR) and ranges.', 
                            'Qualitative variables are presented as', 
                            'percentages of categories within the', 
                            'complete observation set.'))
  
# Table S2: collagen genes of interest -----
  
  insert_msg('Table S2: collagen genes of interest')
  
  suppl_paper_tbl$genes <- globals$genes_interest %>% 
    mutate(entrez_id = mapIds(org.Hs.eg.db, 
                              keys = gene_symbol, 
                              column = 'ENTREZID', 
                              keytype = 'SYMBOL'),
           gene_group = car::recode(gene_group, "'Pro' = 'proline turnover'"), 
           gene_group = factor(gene_group, 
                               c('proline turnover', 
                                 'collagen modification', 
                                 'ECM component', 
                                 'ECM processing', 
                                 'adhesion'))) %>% 
    arrange(gene_group) %>% 
    select(gene_group, gene_symbol, entrez_id) %>% 
    set_names(c('Gene group', 'Gene symbol', 'Entrez ID')) %>% 
    mdtable(label = 'table_s2_genes', 
            ref_name = 'genes', 
            caption = 'Collagen-related genes and their classification.')
  
# Table S3: normal-tumor ------
  
  insert_msg('Table S3: normal - tumor')
  
  suppl_paper_tbl$norm_tumor <- norm_tumor$result_tbl %>% 
    mdtable(label = 'table_s3_normal_tumor', 
            ref_name = 'norm_tumor', 
            caption = paste("Expression of the collagen pathway genes", 
                            "in the malignant and benign tissue compared", 
                            "by paired T test with Cohen's d", 
                            "effect size statistic.", 
                            "P values were corrected for multiple testing", 
                            "with the false discovery rate method.", 
                            "log2-transformed expression values are presented", 
                            "as medians with interquartile ranges (IQR)", 
                            "and ranges.", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S4: collagen-related genes and Gleason scores ---------
  
  insert_msg('Table S4: collagen genes and Gleason scores')
  
  suppl_paper_tbl$gleason <- gs_uni$result_tbl %>% 
    filter(Cohort != 'GSE16560') %>% 
    mdtable(label = 'table_s4_gleason', 
            ref_name = 'gleason', 
            caption = paste("Expression of the collagen pathway genes", 
                            "in cancer samples stratified by the ISUP risk system", 
                            "compared by one-way ANOVA with eta-square", 
                            "effect size statistic.", 
                            "P values were corrected for multiple testing", 
                            "with the false discovery rate method.", 
                            "log2-transformed expression values are presented", 
                            "as medians with interquartile ranges (IQR)", 
                            "and ranges.", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S5: expression of the collagen pathway genes in the clusters -------
  
  insert_msg('Table S5: collagen gene expression in the clusters')
  
  suppl_paper_tbl$cluster <- clust_ft$result_tbl %>% 
    compress(names_to = 'cohort') %>% 
    filter(cohort != 'gse16560') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    relocate(cohort) %>% 
    set_names(c('Cohort', 'Variable', 
                'Collagen low', 'Collagen high', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'table_s5_cluster_collagen_genes', 
            ref_name = 'cluster', 
            caption = paste("Expression of the cluster-defining", 
                            "collagen pathway genes in the collagen clusters", 
                            "of prostate cancer.", 
                            "Statistical significance was assessed by two-tailed", 
                            "T test with Cohen's d effect size statistic.", 
                            "P values were corrected for multiple testing", 
                            "with the false discovery rate method.", 
                            "log2-transformed expression values are presented", 
                            "as medians with interquartile ranges (IQR)", 
                            "and ranges.", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S6: clinical characteristic of the collagen clusters -----
  
  insert_msg('Table S6: Clinical characteristic of the clusters')
  
  suppl_paper_tbl$clinic <- ana_clinic$result_tbl %>% 
    filter(cohort != 'GSE16560') %>% 
    set_names(c('Cohort', 'Variable', 
                'Collagen low', 'Collagen high', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'table_s6_cluster_clinical_characteristic', 
            ref_name = 'clinic',
            caption = paste('Clinical characteristic of the collagen clusters.', 
                            'Numeric variables are presented as medians', 
                            'with interquartile ranges (IQR) and ranges.', 
                            'Nominal variables are presented as percentages', 
                            'and counts of categories within the cluster.'))
  
# Table S7: infiltration -------
  
  insert_msg('Table S7: MCP and xcell counter infiltration')
  
  suppl_paper_tbl$infiltration <- list('MCP Counter' = ana_mcp$result_tbl, 
                                       'xCell' = ana_xcell$result_tbl) %>% 
    compress(names_to = 'algorithm') %>% 
    relocate(algorithm) %>% 
    filter(cohort != 'GSE16560') %>% 
    set_names(c('Algorithm', 'Cohort', 
                'Variable', 'Collagen low', 'Collagen hi', 
                'Significance', 'Effect size')) %>% 
    mdtable(label = 'table_s7_mcp_xcell_infiltration', 
            ref_name = 'infiltration', 
            caption = paste('Non-malignant cell numbers predicted for the collagen', 
                            'clusters by the MCP Counter and xCell algorithms.', 
                            'Statistical significance was assessed by Mann-Whitney',
                            'test with r effect size statistic.', 
                            'P values were corrected for multiple testing with the', 
                            'false discovery method.', 
                            'The table is available as a supplementary Excel file.'))

# Table S8: GSVA, Reactome ------
  
  insert_msg('Table S8: GSVA, Reactome')
  
  ## significant signatures shared by at least 5 cohorts
  
  suppl_paper_tbl$reactome <- ana_reactome$test %>% 
    map(filter, response %in% unique(unlist(ana_reactome$common_significant))) %>% 
    map2(., names(.), ~mutate(.x, cohort = globals$study_labels[.y])) %>% 
    do.call('rbind', .) %>% 
    filter(cohort != 'GSE16560')
  
  ## appending with the signature cluster assignment
  
  suppl_paper_tbl$reactome <- 
    left_join(suppl_paper_tbl$reactome, 
              ana_reactome$signature_clusters %>% 
                extract('assignment') %>% 
                set_names(c('response', 'clust_id')), 
              by = 'response') %>% 
    mutate(clust_id = ana_reactome$cluster_labels[as.character(clust_id)])
  
  ## formatting 
  
  suppl_paper_tbl$reactome <- suppl_paper_tbl$reactome %>% 
    mutate(response = exchange(response, 
                               reactome$lexicon),
           clust_id = stri_replace_all(clust_id, 
                                       fixed = '\n', 
                                       replacement = ', ')) %>% 
    select(cohort, clust_id, response, regulation, 
           estimate, lower_ci, upper_ci, 
           p_adjusted, effect_size) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    arrange(cohort, clust_id, response) %>% 
    set_names(c('Cohort', 
                'Signature classification', 
                'Signature', 
                "Regulation, collagen high vs low", 
                'Fold-regulation', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR for fold-regulation', 
                'Effect size')) %>% 
    mdtable(label = 'table_s8_cluster_gsva_reactome', 
            ref_name = 'gsva', 
            caption = paste("Gene set variation", 
                            "analysis with the Reactome pathway gene signatures.", 
                            "Differences in ssGSEA scores between collagen high", 
                            "and collagen low cancers were investigated by", 
                            "two-tailed T test with Cohen's d effect", 
                            "size statistic.", 
                            "Results for signatures significantly regulated",
                            "with at least weak effect size", 
                            "(d at least 0.2) in at least five", 
                            "cohorts are presented.", 
                            "P values were corrected for multiple testing with",
                            "the false discovery rate method (FDR).", 
                            "The table is available as a supplementary", 
                            "Excel file."))
  
# Table S9: differential gene expression -------
  
  insert_msg('Table S9: differential gene expression')
  
  suppl_paper_tbl$dge <- ana_dge$test %>% 
    map(filter, regulation %in% c('upregulated', 'downregulated')) %>% 
    compress(names_to = 'cohort') %>% 
    filter(cohort != 'gse16560') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, 
           gene_symbol, 
           entrez_id, 
           regulation, 
           estimate, 
           lower_ci, 
           upper_ci, 
           p_adjusted, 
           effect_size) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    set_names(c('Cohort', 
                'Gene symbol', 
                'Entrez ID', 
                'Regulation', 
                'log2 fold-regulation, collagen high vs low', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR', 
                'Effect size')) %>% 
    mdtable(label = 'table_s9_clusters_differential_gene_expression', 
            ref_name = 'dge', 
            caption = paste('Genes differentially expressed in the collagen', 
                            'high cluster as compared', 
                            'with collagen low cancers were identified', 
                            'by two-tailed T test', 
                            "with the 1.25-fold regulation cutoff and", 
                            "the Cohen's d effect size statistic of 0.2", 
                            'P values were corrected for multiple testing with', 
                            'the false discovery rate method (FDR).', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S10: GO enrichment -------
  
  insert_msg('Table S10: GO enrichment')
  
  suppl_paper_tbl$go_enrichment <- 
    map2(ana_go$test, ana_go$common_significant, 
         function(x, y) x %>% 
           map(filter, term %in% y)) %>% 
    map(compress, names_to = 'cohort') %>% 
    map(filter, cohort != 'gse16560')
  
  ## appending with the cGO term cluster assignment
  
  suppl_paper_tbl$go_enrichment <- 
    map2(suppl_paper_tbl$go_enrichment, 
         ana_go$go_assignment, 
         left_join, 
         by = c('go_id', 'term')) %>% 
    map2(., ana_go$go_clust_desc, 
         ~mutate(.x, clust_id = .y[as.character(clust_id)])) %>% 
    compress(names_to = 'dge_group')
  
  ## formatting 
  
  suppl_paper_tbl$go_enrichment <- 
    suppl_paper_tbl$go_enrichment %>% 
    mutate(cohort = globals$study_labels[cohort], 
           dge_group = car::recode(dge_group, 
                                   "'downregulated' = 'collagen low'; 
                                   'upregulated' = 'collagen high'")) %>% 
    select(dge_group, cohort, clust_id, go_id, term, or, p_adjusted) %>% 
    set_names(c('Collagen cluster', 'Cohort', 'GO term classification', 
                'GO ID', 'GO name', 'Enrichment OR', 'pFDR')) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    mdtable(label = 'table_s10_go_enrichment', 
            ref_name = 'go_enrichment', 
            caption = paste('Biological process gene ontology (GO) term', 
                            'enrichment within genes differentially regulated', 
                            'in the collagen clusters. The enrichment analysis', 
                            'was performed with goana tool, enrichment p values', 
                            'were corrected for multiple testing with the false', 
                            'discovery rate (FDR) method.', 
                            'Significant enrichment was defined by pFDR < 0.05', 
                            'and odds ratio (OR) for enrichment within', 
                            'differentially regulated genes of at least 1.44.', 
                            'OR for enrichment within genes upregulated', 
                            'in the collagen high and collagen low clusters', 
                            'are presented for significant GO terms shared by', 
                            'at least five cohorts.', 
                            'The table is available as a supplementary Excel', 
                            'file.'))
  
# Table S11: regulons ---------
  
  insert_msg('Table S11: regulons')
  
  suppl_paper_tbl$regulons <- ana_collectri$test %>% 
    compress(names_to = 'cohort') %>% 
    filter(cohort != 'gse16560', 
           source %in% unname(unlist(ana_collectri$common_significant))) %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, regulation, source, score, p_adjusted) %>% 
    set_names(c('Cohort', 
                'Regulation status', 
                'Regulon symbol', 
                'Regulation magnitude, LM score', 
                'pFDR')) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    mdtable(label = 'table_s11_collectri_regulons', 
            ref_name = 'regulons', 
            caption = paste('Activity of transcriptional regulons in the', 
                            'collagen high cluster as compared with the', 
                            'collagen low cluster predicted by the collecTRI', 
                            'model.', 
                            'Regulon activity was estimated with uni-parameter', 
                            'linear modeling with whole-transcriptome', 
                            'effect sizes of differential gene expression,', 
                            'p values were corrected for multiple testing with', 
                            'the false discovery rate (FDR) method.', 
                            'Linear model scores are presented for regulons', 
                            'significantly activated or inhibited in at least', 
                            'five cohorts.', 
                            'The table is available as a supplementary Excel', 
                            'file.'))
  
# Table S12: PROGENy signaling ---------
  
  insert_msg('Table S12: progeny signaling')
  
  suppl_paper_tbl$signaling <- ana_progeny$test %>% 
    compress(names_to = 'cohort') %>% 
    filter(cohort != 'gse16560') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, regulation, source, score, p_adjusted) %>% 
    set_names(c('Cohort', 
                'Regulation status', 
                'Pathway name', 
                'Regulation magnitude, LM score', 
                'pFDR')) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    mdtable(label = 'table_s12_progeny_signaling', 
            ref_name = 'signaling', 
            caption = paste('Activity of signaling pathways in the collagen high', 
                            'cluster as compared with the collagen low cluster', 
                            'predicted by the PROGENy model.', 
                            'Pathway activity was estimated with multi-parameter', 
                            'linear modeling with whole-transcriptome', 
                            'effect sizes of differential gene expression.', 
                            'P values were corrected for multiple testing with', 
                            'the false discovery rate (FDR) method, ', 
                            'linear model scores serve as measures of pathway', 
                            'activity.', 
                            'The table is available as a supplementary Excel', 
                            'file.'))
  
# Table S13: differentially modulated metabolic reactions -------
  
  insert_msg('Table S13: differentially regulated metabolic reactions')
  
  suppl_paper_tbl$metabolism <- ana_meta$models %>% 
    map(components, 'regulation') %>% 
    map(filter, p_adjusted < 0.05) %>% 
    compress(names_to = 'cohort') %>% 
    filter(cohort != 'gse16560') %>% 
    mutate(cohort = globals$study_labels[cohort])
  
  ## annotation of metabolic reactions and formatting
  
  suppl_paper_tbl$metabolism <- 
    left_join(suppl_paper_tbl$metabolism, 
              suppl_paper_tbl$reactions, 
              by = 'react_id') %>% 
    mutate(fold_reg = log2(fold_reg), 
           lower_ci = log2(lower_ci), 
           upper_ci = log2(upper_ci)) %>% 
    select(cohort, subsystem, react_id, react_name, 
           fold_reg, lower_ci, upper_ci, p_adjusted) %>% 
    set_names(c('Cohort', 
                'Recon metbolic subsystem', 
                'Reaction ID', 
                'Reaction name', 
                'log fold-regulation', 
                'lower 95% CI', 
                'upper 95% CI', 
                'pFDR')) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    mdtable(label = 'table_s13_clusters_ciochamical_reactions', 
            ref_name = 'reactions',
            caption = paste('Biochemical reactions predicted to be', 
                            'significantly activated in collagen high', 
                            'as compared with collagen low cancers.', 
                            'Statistical significance was determined by', 
                            'a Monte Carlo simulation.', 
                            'P values were corrected for multiple testing with', 
                            'the false discovery rate (FDR) method.', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S14: biochemical subsystem enrichment -------
  
  insert_msg('Table S14: biochemical subsystems')
  
  suppl_paper_tbl$subsystems <- ana_meta$subsystems %>% 
    compress(names_to = 'cohort') %>% 
    filter(cohort != 'gse16560') %>% 
    mutate(cohort = globals$study_labels[cohort]) %>% 
    select(cohort, subsystem, status, OR, significance) %>% 
    set_names(c('Cohort', 
                'Subsystem', 
                'Reaction activity status', 
                'Enrichment OR', 
                'pFDR')) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>% 
    mdtable(label = 'table_s14_metabolic_subsystem_enrichment', 
            ref_name = 'subsystems', 
            caption = paste('Results of enrichment analysis', 
                            'for significantly activated and inhibited', 
                            'biochemical reactions within the Recon', 
                            'metabolism subsystem.', 
                            'Statistical significance was determined by', 
                            'random sampling from the entire reaction pool', 
                            'p values were corrected for multiple testing with', 
                            'the false discovery rate (FDR) method.', 
                            'Effect size of enrichment of the subsystem in', 
                            'significantly activated or inhibited reactions was', 
                            'measured by odds ratio (OR) statistic.', 
                            'The table is available as a supplementary', 
                            'Excel file.'))
  
# Table S15: tuning of transcriptomic survival models ------
  
  insert_msg('Table S15: tuning of transcriptomic survival models')
  
  suppl_paper_tbl$tuning <- 
    list(ridge = ridge_surv$opt_lambda["lambda"] %>% 
           map_dfc(signif, 3), 
         elnet = elnet_surv$opt_lambda["lambda"] %>% 
           map_dfc(signif, 3), 
         lasso = lasso_surv$opt_lambda["lambda"] %>% 
           map_dfc(signif, 3), 
         svm = svm_surv$tuning$best_tune[c("type", 
                                           "gamma.mu", 
                                           "kernel")], 
         rf = rf_surv$tuning$best_tune[c("mtry", 
                                         "splitrule", 
                                         "nsplit", 
                                         "nodesize")], 
         gbm = gbm_surv$tuning$best_tune[c("n.trees", 
                                           "shrinkage", 
                                           "interaction.depth", 
                                           "n.minobsinnode")]) %>% 
    map(t) %>% 
    map(as.data.frame) %>% 
    map(set_names, 'value') %>% 
    map(rownames_to_column, 'parameter') %>% 
    map(mutate, value = as.character(value)) %>% 
    compress(names_to = 'algorithm') %>% 
    as_tibble
  
  suppl_paper_tbl$tuning  <- suppl_paper_tbl$tuning %>% 
    mutate(parameter = car::recode(parameter, 
                                   "'lambda' = '\u03BB'; 
                                   'type' = 'SVM model type'; 
                                   'gamma.mu' = '\u03B3'; 
                                   'mtry' = 'number of variables per tree, mtry'; 
                                   'splitrule' = 'splitting rule'; 
                                   'nsplit' = 'number of splits'; 
                                   'nodesize' = 'minimal node size'; 
                                   'n.trees' = 'number of decision trees'; 
                                   'n.minobsinnode' = 'minimal node size'; 
                                   'interaction.depth' = 'interaction depth'"), 
           criterion = car::recode(algorithm, 
                                   "'ridge' = 'minimal deviance, repeated 10-fold CV'; 
                                   'elnet' = 'minimal deviance, repeated 10-fold CV'; 
                                   'lasso' = 'minimal deviance, repeated 10-fold CV'; 
                                   'svm' = 'maximal concordance index, repeated 10-fold CV'; 
                                   'rf' = 'maximal concordance index, out-of-bag predictions'; 
                                   'gbm' = 'minimal deviance, 10-fold CV'"), 
           algorithm = surv_globals$algo_labels[algorithm]) %>% 
    select(algorithm, criterion, parameter, value) %>% 
    set_names(c('Algorithm', 'Selection criterion', 'Parameter', 'Value'))
  
  suppl_paper_tbl$tuning <- suppl_paper_tbl$tuning %>% 
    mdtable(label = 'table_s15_survival_model_tuning', 
            ref_name = 'tuning', 
            caption = paste('Selection of the optimal parameters of machine', 
                            'learining survival models of biochemical', 
                            'relapse-free survival with expression of the', 
                            'collagen-related genes.', 
                            'The selection process was accomplished by', 
                            'cross-validation tuning in the pooled GEO cohort.'))
   
# Table S16: performance of the transcriptomic models at prediction of RFS --------
  
  insert_msg('Table S16: RFS prediction by the transcriptomic models')
  
  suppl_paper_tbl$model_survival <- surv_summary$stats %>%
    map(function(x) if('lower_ci' %in% names(x)) {
      
      mutate(x, 
             c_index = paste0(signif(c_index, 2), ' [95% CI: ', 
                              signif(lower_ci, 2), ' to ', 
                              signif(upper_ci, 2), ']'))
      
    } else {
      
      mutate(x, c_index = as.character(signif(c_index, 2)))
      
    }) %>% 
    map(select, dataset, cohort, c_index, ibs_model) %>% 
    compress(names_to = 'algorithm') %>% 
    relocate(algorithm) %>% 
    mutate(cohort = surv_globals$study_labels[cohort], 
           algorithm = surv_globals$algo_labels[algorithm]) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 2) else x) %>% 
    set_names(c('Algorithm', 
                'Data set type', 
                'Cohort', 
                'Concordance index', 
                'Integrated Brier score')) %>% 
    mdtable(label = 'table_s16_rfs_prediction_stats_transcriptomic_models', 
            ref_name = 'model_survival', 
            caption = paste('Performance of machine learning models at', 
                            'prediction of biochemical relapse-free survival', 
                            'with expression of the collagen-related genes.'))
    
# Table S17: univariable survival analysis --------
  
  insert_msg('Table S17: univariable survival analysis')
  
  ## genes found significant in all cohorts are presented
  
  suppl_paper_tbl$rfs_univariable <- rfs_cut$test %>% 
    compress(names_to = 'cohort') %>% 
    select(marker, gene_symbol, cohort, 
           n_total, n_events, 
           cutoff, n_low, n_high, 
           hr, chisq, significance) %>% 
    filter(gene_symbol %in% reduce(rfs_cut$common_significant, union)) %>% 
    mutate(cohort = surv_globals$study_labels[cohort], 
           chisq = signif(chisq, 2), 
           cutoff = signif(cutoff, 3), 
           hr = signif(hr, 2)) %>% 
    arrange(marker, gene_symbol)
  
  suppl_paper_tbl$rfs_univariable <- 
    suppl_paper_tbl$rfs_univariable %>% 
    set_names(c('Association with survival', 
                'Gene symbol', 
                'Cohort', 
                'Total observations, n', 
                'Relapses, n', 
                'Cutoff, normalized expression', 
                'Low expressors, n', 
                'High expressors, n', 
                'Hazard ratio', 
                'Chi-square statistic', 
                'Significance')) %>% 
    mdtable(label = 'table_s17_univariable_rfs_analysis', 
            ref_name = 'rfs_univariable', 
            caption = paste('Results of univariable analysis of biochemical', 
                            'relapse-free survival with expression of the', 
                            'collagen-related genes.', 
                            'Prostate cancer patients were stratified by',
                            'expression cutoffs corresponding to the largest', 
                            'difference in survival assessed by Mentel-Henszel', 
                            'test.', 
                            'Genes found to be significantly associated with', 
                            'the survival in all analyzed cohorts (pooled GEO,', 
                            'TCGA and DKFZ) are presented.', 
                            'The table is available in a supplementary Excel', 
                            'file.'))
  
# Saving the supplementary tables ------
  
  insert_msg('Saving the supplementary tables')
  
  suppl_paper_tbl$reactions <- NULL
  
  suppl_paper_tbl <- compact(suppl_paper_tbl)
  
  suppl_paper_tbl$cover <- 
    tibble(Table = paste0('Supplementary Table S', 1:length(suppl_paper_tbl)), 
           Caption = map_chr(suppl_paper_tbl, attr, 'caption'))
  
  suppl_paper_tbl <- 
    c(suppl_paper_tbl["cover"], 
      suppl_paper_tbl[names(suppl_paper_tbl) != "cover"])
  
  suppl_paper_tbl %>% 
    set_names(c('Cover', paste0('Table S', 1:(length(suppl_paper_tbl) - 1)))) %>% 
    write_xlsx(path = './report/manuscript_supplementary_tables.xlsx')
  
# END ------
  
  insert_tail()