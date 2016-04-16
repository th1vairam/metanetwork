#!usr/bin/env Rscript

# Submission Script in R
# Clear R console screen output
cat("\014")

# Set library and working directory
.libPaths('/shared/rlibs')
setwd('/shared/Github/metanetwork/R')

# Load libraries
library(plyr)
library(dplyr)
library(synapseClient)

# Login to synapse
key = read.table('/shared/synapseAPIToken')
synapseLogin(username = 'th_vairam',
             apiKey = key$V1)

# Get BIC networks
network.files = synQuery('select * from file where projectId == "syn5584871" and analysisType == "statisticalNetworkReconstruction" and method == "bic" and organism == "HomoSapiens"') %>%
  dplyr::mutate(uniqueName = paste(file.disease, file.tissueTypeAbrv, file.study, file.cogdx, sep = '.'))
  
# Get finished modules
module.files = synQuery('select * from file where projectId == "syn5584871" and analysisType == "statisticalNetworkReconstruction" and method == "bic" and organism == "HomoSapiens" and moduleMethod == "igraph:fast greedy"') 
if (!is.null(module.files)){
  module.files = module.files %>%
    dplyr::mutate(uniqueName = paste(file.disease, file.tissueTypeAbrv, file.study, file.cogdx, sep = '.'))
  
  network.files = network.files %>%
    filter(!(uniqueName %in% module.files$uniqueName))
}
  
# Make directory and write shell scripts for running these files
system('mkdir sgeModuleSubmissions')
fp_all = file(paste('sgeModuleSubmissions/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)

for (i in 1:dim(network.files)[1]){
  id = network.files$file.id[i]
  FNAME = network.files$uniqueName[i]
  fp = file (paste('/shared/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /shared/Github/metanetwork/R/getModules.fastGreedy.R',id,FNAME), 
      file = fp,
      sep = '\n')
  close(fp)
  
  fp_all = file(paste('sgeModuleSubmissions/allSubmissions.sh'),'a+')    
  cat(paste('qsub','-cwd','-V','-l h_vmem=7G', '-l mem_free=7G', paste('/shared/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,sep='.'),
            '-o',paste('/shared/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,'o',sep='.'),
            '-e',paste('/shared/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,'e',sep='.')),
      file = fp_all,
      sep = '\n')
  close(fp_all)
}
