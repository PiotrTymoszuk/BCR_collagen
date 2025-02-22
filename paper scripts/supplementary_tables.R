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
  
# Saving the tables on the disc -----
  
  insert_msg('Saving the tables on the disc')
  
  suppl_tab %>% 
    save_excel(path = './paper/supplementary_tables.xlsx', 
               prefix = 'Supplementary Table S')
  
# END -------
  
  insert_tail()