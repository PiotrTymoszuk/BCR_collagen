# Launches the entire analysis pipeline

  library(soucer)
  
  source_all(
    
    c(## import of the transcriptome and clinical information
      'import.R', 
      
      ## exploratory analysis 
      'exploration.R', 
      
      ## modeling of BCR-free survival and overall survival
      'univariable_survival.R', 
      'BCR.R', 
      'OS.R', 
      
      ## manuscript parts
      'paper.R'
      
      ), 
    
    crash = TRUE, message = TRUE)