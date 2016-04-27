#!/usr/bin/env Rscript

# Function to get modules from network adjacency matrix (from synapse as rda file)
# All algorithms implemented here are the part of WGCNA R package
# Get arguments from comman line
args = commandArgs(TRUE) # args[1] = Synapse id of network file, args[2] = filename, args[3] = algorithm name

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
thisFileName2 = 'getModules.dynamicTreeCut.R'

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
activityDescription = 'Clustering network genes in to modules using TOM and dynamic tree cut methodology'

# Get network from synapse (rda format)
NET_OBJ = synapseClient::synGet(args[1])
FNAME = paste(args[2],args[3],sep='.')

# Load bic networks
load(NET_OBJ@filePath)

# Convert lsparseNetwork to adjacency
adjMat = as.matrix(bicNetworks$rankConsensus$network)*1
adjMat = t(adjMat) + adjMat

# Check adjMat is of class matrix
if (!is(adjMat, 'matrix'))
  stop('adjMat should be of class matrix')

# Method specific parameters
TOMType = 'unsigned' #other options are 'signed'
TOMDenom = 'min' #other options are 'mean',
linkageType = 'ward.D' #
distanceType = 'euclidean' #
minClusterSize = 20
deepSplit = T

# Get topological overlap matrix
TOM = WGCNA::TOMdist(adjMat, TOMType, 
                     TOMDenom, verbose = 4)
colnames(TOM) = colnames(adjMat)
rownames(TOM) = rownames(adjMat)
gc()

# Get hierarchical clusters from TOM
d = dist(TOM, method = distanceType)
clust = fastcluster::hclust(d, method = linkageType)
gc()

# Get individual clusters from the hierarchical tree
mod = try(
  dynamicTreeCut::cutreeDynamic(clust,
                                minClusterSize = minClusterSize,
                                method = 'tree',
                                deepSplit = deepSplit))
gc()

if (is(mod,'try-error'))
  mod = rep(0, dim(adjMat)[1])
names(mod) = rownames(adjMat)

# Get individual clusters from the igraph community object
geneModules = mod %>%
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
MOD_OBJ = File(paste(FNAME,'tsv',sep='.'), 
               name = paste('BIC Rank Consensus', 'fast cluster', 
                            TOMType, TOMDenom, linkageType, distanceType), 
               parentId = obj$properties$id)
annotations(MOD_OBJ) = annotations(NET_OBJ)
MOD_OBJ@annotations$fileType = 'tsv'
MOD_OBJ@annotations$analysisType = 'moduleIdentification'
MOD_OBJ@annotations$moduleMethod = paste('fastCluster', TOMType, TOMDenom, linkageType, distanceType,sep=':')

# Convert lsparseNetwork to igraph graph object
g = igraph::graph.adjacency(bicNetworks$rankConsensus$network, mode = 'upper', weighted = T, diag = F)
MOD_OBJ@annotations$modularity = modularity(g, mod)

MOD_OBJ = synStore(MOD_OBJ, 
                   executed = list(thisFile1, thisFile2),
                   used = NET_OBJ,
                   activityName = activityName,
                   activityDescription = activityDescription)

# Write completed files to synapse
write.table(MOD_OBJ$properties$id, file = 'CompletedModuleIDs.txt', sep='\n', append=T, quote=F, col.names=F, row.names=F)
writeLines(paste('Completed',FNAME,'and stored in',MOD_OBJ$properties$id))
