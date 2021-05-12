# Test whether reduced sampling of mainland species
# and assumption of uniform incomplete taxon sampling erodes signal of 
# mainland effect on diversification rates with hisse

require(ape)
library(hisse)
library(ggplot2)
library(reshape)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Hisse')

# Remove tree and remove outgroup
tree <- read.tree('../FinalData/Anole_TimeCal_NoOut.tree')
tree <- drop.tip(tree, c('agassizi', 'bicaorum', 'lineatus', 'roatanensis',
                         'utilensis', 'pinchoti', 'concolor', 'townsendi'))

# this analysis tries to replicate the finding of Burress and Muñoz 2021:
# they found no effect of an island effect on diversification.
# However, the trees used had at most half as many mainland species as 
# island species and assumed uniform incomplete taxon sampling. 
# In fact, there are more mainland than island species. 
# Thus, rather than simply repeating their analyses using the same trees
# (in effect just replicating their analysis), we want to ask whether 
# 1) reduction in tree size with a bias towards island species
# and 2) assumption of uniform taxon sampling leads to a hisse model that 
# shows no differential diversification

# This will be accomplished by randomly downsampling the tree to equal
# size, and fitting the full hisse model assuming uniform incomplete
# taxon sampling

# Read in data (treating mainland/island as a binary) and remove outgroup. Format for analysis in hisse (vector of states, species as names)
dat <- read.csv('../FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv')
dat <- data.frame('Species' = dat$Species[1:379], 'Island' = rowSums(dat[1:379,33:38]))
dat <- dat[which(dat$Species %in% tree$tip.label),]

m.spp <- as.character(dat$Species[which(dat$Island == 0)])
i.spp <- as.character(dat$Species[which(dat$Island == 1)])

trees <- list()
dats <- list()
for(i in 1:5){
  keep.m <- sample(m.spp, size = 67, replace = F)
  keep.i <- sample(i.spp, size = 124, replace = F)
  keep.spp <- c(keep.i, keep.m)
  trees[[i]] <- 
    drop.tip(tree, 
             tip = tree$tip.label[-which(tree$tip.label 
                                         %in% keep.spp)])
  dats[[i]] <- dat[which(dat$Species %in% trees[[i]]$tip.label),]
}

# Now run the HiSSE model on these trees

hisse.fits <- list()
hisse.recons <- list()
for(i in 1:5){
  ##############################
  # Full HiSSE All Transitions #
  ##############################
  print(i)
  trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
  
  hisse.fits[[i]] <- hisse(phy = trees[[i]], 
                           data = dats[[i]], 
                           f = c(0.45, 0.45), 
                           hidden.states = T, 
                           turnover.anc = c(1,2,3,4),
                           eps.anc = c(1,2,3,4), 
                           trans.rate = trans.rates.hisse, 
                           output.type = "raw")
  hisse.recons[[i]] <- MarginRecon(phy = trees[[i]], 
                                   data = dats[[i]], 
                                   f = c(0.45, 0.45),
                                   pars=hisse.fits[[i]]$solution, 
                                   hidden.states = T, 
                                   aic=hisse.fits[[i]]$AICc, 
                                   n.cores = 1, 
                                   root.p = c(0.5, 0, 0.5, 0))

}





recon_hisse.full_AllTrans <- MarginRecon(tree, dat, f = c(0.82, 0.85), pars=hisse.full_AllTrans$solution, hidden.states = T, aic=hisse.full_AllTrans$AICc, n.cores = 10, root.p = c(0.5, 0, 0.5, 0))

save(hisse.full_AllTrans, file = 'hisse.full_AllTrans_Object.Rsave')
save(recon_hisse.full_AllTrans, file='recon_hisse.full_AllTrans.Rsave')

hisse.full_AllTrans.logL <- hisse.full_AllTrans$loglik
hisse.full_AllTrans.AIC <- hisse.full_AllTrans$AIC
hisse.full_AllTrans.AICc <- hisse.full_AllTrans$AICc
