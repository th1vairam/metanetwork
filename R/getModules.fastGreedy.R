#!/usr/bin/env Rscript

# Function to get modules from network adjacency matrix (from synapse as rda file)
# Get arguments from comman line
args = commandArgs(TRUE)

# Clear R console screen output
cat("\014")

# Set library and working directory
.libPaths('/shared/rlibs')
setwd('/shared/Github/metanetwork/R')

# Load libraries
library(synapseClient)
library(dplyr)
library(WGCNA)
library(tools)
library(stringr)
library(igraph)
library(CovariateAnalysis)

# Needs the dev branch
library(githubr)

# Login to synapse
key = read.table('/shared/synapseAPIToken')
synapseLogin(username = 'th_vairam',
             apiKey = key$V1)

# Get github links for provenance
thisFileName1 = 'makeModuleSubmissionScripts.R'
thisFileName2 = 'getModules.fastGreedy.R'

# Github link
thisRepo <- getRepo(repository = "th1vairam/metanetwork", 
                    ref="branch", 
                    refName='modules')

thisFile1 <- getPermlink(repository = thisRepo,
                         repositoryPath=paste0('R/', thisFileName1))
thisFile2 <- getPermlink(repository = thisRepo,
                         repositoryPath=paste0('R/', thisFileName2))

# Synapse specific parameters
activityName = 'Module Identification'
activityDescription = 'Clustering network genes in to modules using fast greedy algorithm of igraph'

# Get network from synapse (rda format)
NET_OBJ = synapseClient::synGet(args[1])
FNAME = tools::file_path_sans_ext(NET_OBJ$properties$name) %>%
  gsub('Networks','Modules',.)

# Load bic networks
load(NET_OBJ@filePath)

# Convert lsparseNetwork to igraph graph object
g = igraph::graph.adjacency(bicNetworks$rankConsensus$network, mode = 'upper', weighted = T, diag = F)
gc()

# Get modules using fast.greedy method (http://arxiv.org/abs/cond-mat/0408187)
mod = igraph::fastgreedy.community(g)
gc()

# Get individual clusters from the igraph community object
geneModules = igraph::membership(mod) %>% unclass %>%
  CovariateAnalysis::rownameToFirstColumn('EnsembleID') %>%
  dplyr::rename(moduleNumber = DF)

filteredModules = geneModules %>% 
  group_by(moduleNumber) %>%
  summarise(counts = n()) %>%
  filter(counts >= 20)

geneModules$moduleNumber[!(geneModules$moduleNumber %in% filteredModules$moduleNumber)] = 0

# Change cluster number to color labels
geneModules$moduleLabel = WGCNA::labels2colors(geneModules$moduleNumber)

# Write results to synapse
# Create a folder for module results
obj = Folder(name = 'Modules', parentId = NET_OBJ$properties$parentId)
obj = synStore(obj)

write.table(geneModules, paste(FNAME,'tsv',sep='.'), sep='\t', row.names=F, quote=F)
MOD_OBJ = File(paste(FNAME,'tsv',sep='.'), name = FNAME, parentId = obj$properties$id)
annotations(MOD_OBJ) = annotations(NET_OBJ)
MOD_OBJ@annotations$fileType = 'tsv'
MOD_OBJ@annotations$moduleMethod = paste('igraph',algorithm(mod),sep=':')
MOD_OBJ@annotations$modularity = modularity(mod)
MOD_OBJ = synStore(MOD_OBJ, 
                   executed = list(thisFile1, thisFile2),
                   used = NET_OBJ,
                   activityName = activityName,
                   activityDescription = activityDescription)

# Write completed files to synapse
write.table(MOD_OBJ$properties$id, file = 'CompletedModuleIDs.txt', sep='\n', append=T, quote=F, col.names=F, row.names=F)
writeLines(paste('Completed',FNAME,'and stored in',MOD_OBJ$properties$id))
