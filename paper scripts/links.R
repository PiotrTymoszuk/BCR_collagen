# Handling of common links and biobiography in the Rmarkdown project

  insert_head()

# container ------

  proj_links <- list()

# Links to development packages --------

  insert_msg('Links to development packages')

  proj_links <-
    list('ExDA' = 'https://github.com/PiotrTymoszuk/ExDA',
         'trafo' = 'https://github.com/PiotrTymoszuk/trafo',
         'figur' = 'https://github.com/PiotrTymoszuk/figur',
         'clustTools' = 'https://github.com/PiotrTymoszuk/clustTools',
         'microViz' = 'https://github.com/PiotrTymoszuk/microViz', 
         'coxExtensions' = 'https://github.com/PiotrTymoszuk/coxExtensions', 
         'kmOptimizer' = 'https://github.com/PiotrTymoszuk/kmOptimizer', 
         'fastTest' = 'https://github.com/PiotrTymoszuk/fastTest', 
         'graphExtra' = 'https://github.com/PiotrTymoszuk/graphExtra', 
         'htGLMNET' = 'https://github.com/PiotrTymoszuk/htGLMNET') %>%
    compress(names_to = 'obj_name',
             values_to = 'x') %>%
    mutate(ref_name = paste0('_', obj_name, '_'))

  proj_links <- proj_links[c('ref_name', 'x')] %>%
    pmap(mdlink) %>%
    set_names(proj_links$obj_name)
  
# links to the data sets --------
  
  data_links <- 
    c(TCGA = 'https://www.cbioportal.org/study/summary?id=prad_tcga_pan_can_atlas_2018', 
      DKFZ = 'https://www.cbioportal.org/study/summary?id=prostate_dkfz_2018', 
      GSE16560 = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16560', 
      GSE54460 = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54460', 
      GSE70768 = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70768', 
      GSE70769 = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70769', 
      GSE220095 = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220095') %>% 
    tibble(ref_name = names(.), 
           x = .)
  
  data_links <- data_links %>%
    pmap(mdlink) %>% 
    set_names(data_links$ref_name)
  
# reading the bibliography -------
  
  insert_msg('Rending the bibliography')
  
  coll_bib <- read_bib('./paper/markdown/coll_biblio.bib')

# END ------

  insert_tail()
