# Import and wrangling of the GSE16560 data set

  insert_head()
  
# container ------
  
  gse16560 <- list()
  
# reading from GEO --------
  
  insert_msg('Reading from GEO')
  
  gse16560$raw <- getGEO(GEO = 'GSE16560', destdir = './data/GSE16560')
  
# clinical data -------
  
  insert_msg('Clinical data')
  
  gse16560$clinic <- pData(gse16560$raw[[1]]) %>% 
    as_tibble
  
  gse16560$clinic <- gse16560$clinic %>% 
    transmute(sample_id = geo_accession, 
              patient_id = title, 
              tissue_type = factor('tumor', c('normal', 'tumor')), 
              age = as.numeric(`age:ch1`), 
              purity = as.numeric(`cancer.percent:ch1`), 
              os_months = as.numeric(`fup.month:ch1`), 
              death = car::recode(`status.all:ch1`, 
                                  "'Alive' = 0; 'Dead' = 1"), 
              death = as.numeric(death), 
              gleason_sum = as.numeric(`gleason:ch1`), 
              gleason_major = as.numeric(`major.gleason:ch1`), 
              gleason_minor = as.numeric(`minor.gleason:ch1`), 
              gleason_simple = cut(gleason_sum, 
                                   c(-Inf, 6, 7, Inf), 
                                   c('ISUP1', 'ISUP2', 'ISUP3+')))
  
  gse16560$clinic[c('gleason_sum', 'gleason_major', 'gleason_minor')] <- 
    gse16560$clinic[c('gleason_sum', 'gleason_major', 'gleason_minor')] %>% 
    map_dfc(factor)
  
# annotation --------
  
  insert_msg('Annotation')
  
  gse16560$annotation <- fData(gse16560$raw[[1]]) %>% 
    as_tibble
  
  gse16560$annotation <- gse16560$annotation %>% 
    transmute(probe_id = ID, 
              entrez_id = as.character(GeneID)) %>% 
    annotate_raw_symbol
  
# expression --------
  
  insert_msg('Expression')
  
  ## duplicated probes are aggregated by arithmetic means of log2 expression
  ## signals
  
  gse16560$expression <- exprs(gse16560$raw[[1]]) %>% 
    integrate_expression(gse16560$annotation) %>% 
    left_join(gse16560$clinic[c('sample_id', 'patient_id', 'tissue_type')], ., 
              by = 'sample_id')
  
  gse16560$annotation <- gse16560$annotation %>% 
    filter(!duplicated(gene_symbol), 
           !duplicated(entrez_id))
  
# Caching the results ------
  
  insert_msg('Caching the results')
  
  gse16560 <- gse16560[c("clinic", "annotation", "expression")]
  
  save(gse16560, file = './data/gse16560.RData')
  
# END --------
  
  insert_tail()