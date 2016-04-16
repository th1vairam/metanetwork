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
network.files = synQuery('select * from file where projectId == "syn5584871" and analysisType == "statisticalNetworkReconstruction" and method == "bic" and organism == "HomoSapiens"')

# Make directory and write shell scripts for running these files
system('mkdir sgeModuleSubmissions')
fp_all = file(paste('sgeModuleSubmissions/allSubmissions.sh'),'w+')    
cat('#!/bin/bash',file=fp_all,sep='\n')
close(fp_all)

for (id in network.files$file.id){
  fp = file (paste('/shared/Github/metanetwork/R/sgeModuleSubmissions/SUB',id,sep='.'), "w+")
  cat('#!/bin/bash', 
      'sleep 30', 
      paste('Rscript /shared/Work/Github/metanetwork/R/getModules.fastGreedy.R',id), 
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
