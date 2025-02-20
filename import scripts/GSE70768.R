# Import of the GSE70768 data set

  insert_head()
  
# container -------
  
  gse70768 <- list()
  
# raw data ------
  
  insert_msg('Raw data')
  
  gse70768$raw <- getGEO(GEO = 'GSE70768', destdir = './data/GSE70768')
  
# clinical information -------
  
  insert_msg('Clinical information')
  
  gse70768$clinic <- pData(gse70768$raw[[1]]) %>% 
    as_tibble
  
  gse70768$clinic <- gse70768$clinic %>% 
    transmute(sample_id = geo_accession, 
              patient_id = stri_extract(title, regex = '(AL|AR|NS|TB).*$'), 
              tissue_type = ifelse(`sample type:ch1` == 'Tumour', 
                                   'tumor', 'normal'), 
              tissue_type = factor(tissue_type, c('normal', 'tumor')), 
              age = as.numeric(`age at diag:ch1`), 
              relapse = car::recode(`biochemical relapse (bcr):ch1`, 
                                    "'Y' = 1; 'N' = 0"), 
              relapse = as.numeric(relapse), 
              rfs_months = as.numeric(`time to bcr (months):ch1`), 
              rfs_months = ifelse(is.na(rfs_months), 
                                  as.numeric(`total follow up (months):ch1`), 
                                  rfs_months), 
              ct_stage = stri_extract(`clinical stage:ch1`, regex = 'T\\d{1}'), 
              ct_stage = factor(ct_stage, paste0('T', 1:5)), 
              cn_stage = stri_extract(`clinical stage:ch1`, regex = 'N\\d{1}'), 
              cn_stage = factor(cn_stage, paste0('N', 0:5)), 
              cm_stage = stri_extract(`clinical stage:ch1`, regex = 'M\\d{1}'), 
              cm_stage = factor(cm_stage, paste0('M', 0:5)), 
              pt_stage = stri_extract(`pathology stage:ch1`, regex = 'T\\d{1}'), 
              pt_stage = factor(pt_stage, paste0('T', 1:5)), 
              pn_stage = stri_extract(`pathology stage:ch1`, regex = 'N\\d{1}'), 
              pn_stage = factor(pn_stage, paste0('N', 0:5)), 
              pm_stage = stri_extract(`pathology stage:ch1`, regex = 'M\\d{1}'), 
              pm_stage = factor(pm_stage, paste0('M', 0:5)), 
              surgical_margins = car::recode(`positive surgical margins (psm):ch1`, 
                                             "'Y' = 'positive'; 'N' = 'negative'"), 
              surgical_margins = factor(surgical_margins, c('negative', 'positive')), 
              psa_diagnosis = as.numeric(`psa at diag:ch1`), 
              purity = stri_extract(`tumour %:ch1`, regex = '\\d+'), 
              purity = as.numeric(purity), 
              gleason_sum = stri_extract_all(`tumour gleason:ch1`, 
                                             regex = '\\d{1}'), 
              gleason_major = map_chr(gleason_sum, function(x) if(all(is.na(x))) NA else x[[2]]), 
              gleason_minor = map_chr(gleason_sum, function(x) if(all(is.na(x))) NA else x[[3]]), 
              gleason_sum = map_chr(gleason_sum, function(x) if(all(is.na(x))) NA else x[[1]]), 
              gleason_sum = as.numeric(gleason_sum),
              gleason_simple = cut(gleason_sum, 
                                   c(-Inf, 6, 7, Inf), 
                                   c('ISUP1', 'ISUP2', 'ISUP3+')))
  
  gse70768$clinic[c('gleason_sum', 'gleason_major', 'gleason_minor')] <- 
    gse70768$clinic[c('gleason_sum', 'gleason_major', 'gleason_minor')] %>% 
    map_dfc(factor)
  
  gse70768$clinic <- gse70768$clinic %>% 
    map_dfc(function(x) if(is.factor(x)) droplevels(x) else x)
  
# annotation -------
  
  insert_msg('Annotation')
  
  gse70768$annotation <- fData(gse70768$raw[[1]]) %>%
    as_tibble
  
  gse70768$annotation <- gse70768$annotation %>% 
    transmute(probe_id = ID, 
              entrez_id = as.character(Entrez_Gene_ID)) %>% 
    annotate_raw_symbol
  
# expression --------
  
  insert_msg('Expression')
  
  ## duplicated probes are aggregated by arithmetic means of log2 expression
  ## signals
  
  gse70768$expression <- exprs(gse70768$raw[[1]]) %>% 
    integrate_expression(gse70768$annotation) %>% 
    left_join(gse70768$clinic[c('sample_id', 'patient_id', 'tissue_type')], ., 
              by = 'sample_id')
  
  gse70768$annotation <- gse70768$annotation %>% 
    filter(!duplicated(gene_symbol), 
           !duplicated(entrez_id))
  
# Caching the results ------
  
  insert_msg('Caching the results')
  
  gse70768 <- gse70768[c("clinic", "annotation", "expression")]
  
  save(gse70768, file = './data/gse70768.RData')
  
# END --------
  
  insert_tail()
  