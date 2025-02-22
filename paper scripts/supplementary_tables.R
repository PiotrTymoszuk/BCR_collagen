# Supplementary tables for the revised manuscript 

  insert_head()
  
# container ------
  
  suppl_tab <- list()
  
# Characteristic of the cohorts: pooled GEO, TCGA, DKFZ, and GSE1650 ----
  
  insert_msg('Characteristic of the cohorts')
  
  suppl_tab$cohorts <- cohorts$stats %>% 
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
  
  suppl_tab$pool <- expl_pool$result_tbl %>% 
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
  
  suppl_tab$gse16560 <- expl_os$result_tbl %>% 
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
  
# Tuning parameters, machine learning models ------
  
  insert_msg('Tuning parameters, machine learning models')
  
  
# Saving the tables on the disc -----
  
  insert_msg('Saving the tables on the disc')
  
  suppl_tab %>% 
    save_excel(path = './paper/supplementary_tables.xlsx', 
               prefix = 'Supplementary Table S')
  
# END -------
  
  insert_tail()