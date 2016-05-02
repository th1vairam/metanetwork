#!usr/bin/env Rscript

#### Code to delploy module preservation analysis ####
# reference: cell type markers
# test: rank consensus based on bic criteria

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

#### Get all the completed networks and modules ####
# Get BIC networks
network.files = synQuery('select * from file where projectId == "syn5584871" and analysisType == "statisticalNetworkReconstruction" and method == "bic" and organism == "HomoSapiens"') %>%
  dplyr::filter(file.cogdx == 1 | file.disease == 'Control') %>%
  dplyr::mutate(uniqueName = paste(file.disease, file.tissueTypeAbrv, file.study, file.cogdx, sep = '.'))

# Get finished modules
module.files = synQuery('select * from file where projectId == "syn5584871" and analysisType == "moduleIdentification" and method == "bic" and organism == "HomoSapiens"')  %>%
  dplyr::filter(file.cogdx == 1 | file.disease == 'Control') %>%
  dplyr::mutate(uniqueName = paste(file.disease, file.tissueTypeAbrv, file.study, file.cogdx, sep = '.'))

#### Get reference network: cell type markers ####
# Download AD related gene sets from synapse
GL_OBJ = synGet('syn4893059');
load(GL_OBJ@filePath) # this will load a list named GeneSets

# Convert hgnc_symbols to ensembl gene ids (also select cell types only)
Hs = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org") # use this one when biomart.org is down
Hs = useDataset("hsapiens_gene_ensembl", Hs)
GeneSets = lapply(GeneSets[grep('Zhang',names(GeneSets),value=T)], function(x, Hs){
  human_ensg2symbol = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                            filters = "hgnc_symbol",                         
                            values = x,
                            mart = Hs)
  return(human_ensg2symbol[,'ensembl_gene_id'])
}, Hs)
ref.g = lapply(GeneSets, function(x){
  adjMat = matrix(rep(1, length(x) * length(x)), length(x), length(x))
  rownames(adjMat) = x; colnames(adjMat) = x
  igraph::graph.adjacency(adjMat, mode = 'undirected', diag = F)
})
ref.gg = do.call(union, ref.g)

ref.mod = lapply(GeneSets, function(x){
  x = data.frame(EnsembleID = x)
}) %>% 
  rbindlist(idcol = 'moduleLabel') %>%
  dplyr::mutate(moduleNumber = as.numeric(factor(moduleLabel)))
  
# Download adjacency matrices from synapse
fileNames = intersect(network.files$uniqueName, module.files$uniqueName)
fileNames = fileNames[3]
for (file.name in fileNames){
  # Get network from synapse (rda format)
  net.obj = network.files$file.id[network.files$uniqueName == file.name] %>%
    synGet
  
  # Load bic networks
  load(net.obj@filePath)
  
  # Convert lsparseNetwork to igraph graph object
  test.g = igraph::graph.adjacency(as(bicNetworks$rankConsensus$network, 'dMatrix'), 
                              mode = 'upper', diag = F)
  gc()
  
  # Get test modules from synapse (tsv format)
  ind = which(module.files$uniqueName == file.name)
  for (j in 1:length(ind)){
    test.mod.obj = module.files$file.id[ind[j]] %>%
      synGet
    
    test.mod = read.table(test.mod.obj@filePath, header=T, sep = '\t')
    
    # Perform module preservation analysis
    results = lapply(unique(ref.mod$moduleLabel), modPreservStats, ref.gg, test.g, ref.mod, test.mod) %>%
      rbindlist
    
    # Create module analysis folders in synapse
    fold.obj = Folder(name = 'Module Analysis', parentId = getParentEntity(test.mod.obj@properties$parentId)@properties$id)
    fold.obj = synStore(fold.obj)
    
    # Write results to synapse
    write.table(results, 
                file = paste(file.name, test.mod.obj@annotations$moduleMethod, 'tsv', sep = '.'), 
                sep = '\t', row.names = F, quote=F)
    file.obj = File(paste(file.name, test.mod.obj@annotations$moduleMethod, 'tsv', sep = '.'),
                    name = test.mod.obj@properties$name,
                    parentId = fold.obj@properties$id)
    annotations(file.obj) =annotations(test.mod.obj)
    file.obj@annotations$analysisType = "modulePreservationAnalysis"
    file.obj@annotations$modularity = NULL
    file.obj = synStore(file.obj, used = c(test.mod.obj$properties$id, 'syn4893059', net.obj@properties$id),
                        activityName = 'Module Preservation Analysis')
    writeLines(paste(file.name, test.mod.obj@annotations$moduleMethod))
  }
}