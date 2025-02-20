# Import and wrangling of the TCGA data set (version Cell 2018, cBioportal)

  insert_head()
  
# container ------
  
  tcga <- list()
  
# loading legacy GDC data which contain Gleason scores ------
  
  insert_msg('Legacy GDC data')
  
  ## for source: see our 2020 paper on CAND1
  
  load('./data/TCGA/tcga_gdc.RData')
  
  tcga$clinic$gdc <- tcga_data$clinical
  
  rm(tcga_data)
  
# clinical data -------
  
  insert_msg('Clinical information')
  
  tcga$clinic[c('patient', 'sample', 'normals')] <- 
    c('./data/TCGA/data_clinical_patient.txt', 
      './data/TCGA/data_clinical_sample.txt', 
      './data/TCGA/normals/data_mrna_seq_v2_rsem_normal_samples.txt') %>% 
    list(file = ., 
         skip = c(4, 4, 0)) %>% 
    pmap(read_tsv)
  
  ## patient's data
  
  tcga$clinic$patient <- tcga$clinic$patient %>% 
    transmute(patient_id = PATIENT_ID, 
              age = as.numeric(AGE), 
              neoadjuvant = factor(tolower(HISTORY_NEOADJUVANT_TRTYN), 
                                   c('no', 'yes')), 
              pt_stage = stri_extract(PATH_T_STAGE, regex = 'T\\d{1}'), 
              pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')), 
              pn_stage = stri_extract(PATH_N_STAGE, regex = 'N\\d{1}'), 
              pn_stage = factor(pn_stage, c('N0', 'N1')), 
              race = factor(RACE), 
              radiation_therapy = factor(tolower(RADIATION_THERAPY), 
                                         c('no', 'yes')), 
              death = as.numeric(stri_extract(OS_STATUS, regex = '^\\d{1}')), 
              os_months = as.numeric(OS_MONTHS), 
              tumor_death = as.numeric(stri_extract(DSS_STATUS, regex = '^\\d{1}')), 
              tss_months = as.numeric(DSS_MONTHS), 
              ## progression == biochemical relapse
              relapse = as.numeric(stri_extract(PFS_STATUS, regex = '^\\d{1}')), 
              rfs_months = as.numeric(PFS_MONTHS))
  
  ## patient's data, GDC
  
  tcga$clinic$gdc <- tcga$clinic$gdc %>% 
    transmute(patient_id = bcr_patient_barcode, 
              ## PSA * 10 for the common unit with other data sets 
              psa_diagnosis = as.numeric(psa) * 10, 
              gleason_sum = as.numeric(gleason_score), 
              gleason_major = factor(as.numeric(primary_gleason)), 
              gleason_minor = factor(as.numeric(secondary_gleason)), 
              gleason_simple = cut(gleason_sum, 
                                   c(-Inf, 6, 7, Inf), 
                                   c('ISUP1', 'ISUP2', 'ISUP3+')))
  
  ## sample
  
  tcga$clinic$sample <- tcga$clinic$sample %>% 
    transmute(sample_id = SAMPLE_ID, 
              patient_id = PATIENT_ID,
              tissue_type = factor('tumor', c('normal', 'tumor')), 
              tmb = as.numeric(TMB_NONSYNONYMOUS))
  
  ## normal samples
  
  tcga$clinic$normals <- tcga$clinic$normals %>% 
    select(starts_with('TCGA')) %>%
    names
  
  tcga$clinic$normals <- 
    tibble(sample_id = tcga$clinic$normals, 
           patient_id = stri_extract(tcga$clinic$normals, 
                                     regex = 'TCGA-\\w{2}-\\w{4}'), 
           tissue_type = factor('normal', c('normal', 'tumor')))
  
  ## merging
  
  tcga$clinic <- tcga$clinic[c("patient", "gdc", "sample")] %>% 
    reduce(left_join, by = 'patient_id') %>% 
    full_rbind(., tcga$clinic$normals)
  
# Expression and annotation --------
  
  insert_msg('Expression and annotation')
  
  tcga$expression[c('tumor', 'normal')] <- 
    c('./data/TCGA/data_mrna_seq_v2_rsem.txt', 
      './data/TCGA/normals/data_mrna_seq_v2_rsem_normal_samples.txt') %>% 
    map(read_tsv)
  
  ## annotation
  
  tcga$annotation <- tcga$expression$tumor %>% 
    transmute(probe_id = as.character(1:nrow(.)),
              entrez_id = as.character(Entrez_Gene_Id)) %>% 
    filter(!is.na(entrez_id)) %>% 
    annotate_raw_symbol
  
  ## expression: log2 transformation and aggregatin of duplicated 
  ## genes by arithmetic mean
  
  tcga$expression <- tcga$expression %>% 
    map(select, starts_with('TCGA')) %>% 
    map(as.matrix) %>% 
    map(~set_rownames(.x, as.character(1:nrow(.x)))) %>% 
    map(~log2(.x + 1)) %>% 
    map_dfr(integrate_expression, 
            tcga$annotation)
  
  ## merging with the patient ID and tissue type information
  
  tcga$expression <- 
    inner_join(tcga$clinic[c('sample_id', 'patient_id', 'tissue_type')], 
               tcga$expression, by = 'sample_id')

# Mutations ------
  
  insert_msg('Mutations')
  
  plan('multisession')
  
  ## non-silent, non-synonymous ones, mutations with assigned Entrez ID
  
  tcga$mutations <- 
    read_tsv('./data/TCGA/data_mutations.txt')
  
  tcga$mutations <- tcga$mutations %>% 
    filter(!Consequence %in% c('synonymous_variant'), 
           !is.na(Entrez_Gene_Id), 
           !is.na(Hugo_Symbol)) %>%
    extract_mutations
  
  plan('sequential')
  
# Common samples with clinical and transcriptome -------
  
  insert_msg('Common samples')
  
  tcga$common_samples <- tcga[c("clinic", "expression")] %>% 
    map(~.x$sample_id) %>% 
    reduce(intersect)
  
  tcga[c("clinic", "expression", "mutations")] <- 
    tcga[c("clinic", "expression", "mutations")] %>% 
    map(filter, sample_id %in% tcga$common_samples)
  
# caching the results -------
  
  insert_msg('Caching the results')
  
  tcga <- tcga[c("clinic", "annotation", "expression", "mutations")]
  
  save(tcga, file = './data/tcga.RData')
  
# END -------
  
  insert_tail()