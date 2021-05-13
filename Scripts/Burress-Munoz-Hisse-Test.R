# this analysis tries to replicate the finding of Burress and Muñoz 2021:
# they found no effect of an island effect on diversification.
# However, the trees used had at most half as many mainland species as 
# island species and assumed uniform incomplete taxon sampling. 
# In fact, there are more mainland than island species. 
# To test whether this species sampling may have influenced their results,
# we fit full hisse models using the largess set of species (from Zheng and Wiens).
# We do so first by pruning Poe's tree down to these species. 

require(ape)
library(hisse)
library(ggplot2)
library(reshape)
library(gghisse)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Hisse')

# Remove tree and remove outgroup
tree <- read.tree('../FinalData/Anole_TimeCal_NoOut.tree')
# we still drop these species because we can't to use our mainland/island 
# designation. In other words, we don't want to change our actual
# data for the binary trait - we just want to change species sampling/tree
tree <- drop.tip(tree, c('agassizi', 'bicaorum', 'lineatus', 'roatanensis',
                         'utilensis', 'pinchoti', 'concolor', 'townsendi'))

# reduce down to the species in zheng & wiens 2016
zheng.spp <- read.table('../zheng-wiens-species.txt')

tree <- drop.tip(tree, tip = tree$tip.label[-which(tree$tip.label %in% zheng.spp$V1)])

# Read in data (treating mainland/island as a binary) and remove outgroup. Format for analysis in hisse (vector of states, species as names)
dat <- read.csv('../FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv')
dat <- data.frame('Species' = dat$Species[1:379], 'Island' = rowSums(dat[1:379,33:38]))
dat <- dat[which(dat$Species %in% tree$tip.label),]

# Now run the HiSSE model on these trees
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)

zheng.poe.tree.fit <- 
  hisse(phy = tree, 
        data = dat, 
        f = c(0.45, 0.45), 
        hidden.states = T, 
        turnover.anc = c(1,2,3,4),
        eps.anc = c(1,2,3,4), 
        trans.rate = trans.rates.hisse, 
        output.type = "raw")
zheng.poe.tree.recon <- 
  MarginRecon(phy = tree, 
              data = dat, 
              f = c(0.45, 0.45),
              pars=zheng.fit$solution, 
              hidden.states = T, 
              aic=zheng.fit$AICc, 
              n.cores = 1, 
              root.p = c(0.5, 0, 0.5, 0))

colfunc<-colorRampPalette(c('#313695', '#4575b4', '#74add1', 
                            '#abd9e9', '#e0f3f8', '#ffffbf', 
                            '#fee090', '#fdae61', '#f46d43', 
                            '#d73027', '#a50026', '#a50026'))
# hisse.contMap <- 
#   plot.hisse.states(zheng.poe.tree.recon, 
#                     rate.param = 'net.div', 
#                     do.observed.only = T, 
#                     rate.colors = colfunc(1000), 
#                     edge.width.state = 1, 
#                     show.tip.label = F, 
#                     legend.kernel.rates = 'gaussian', 
#                     edge.width.rate = 3.5, 
#                     legend = 'all', 
#                     legend.cex = 0.8)

cols <- c("#faae61", "#d72027", "#add8a4", "#2b83bb")

processed_hisse.zheng.poe.tree <- h_process_recon(hisse_recon=zheng.poe.tree.recon)
zheng_rates_poe.tree_plot <- h_scatterplot(
  processed_recon=processed_hisse.zheng.poe.tree,
  parameter="speciation", 
  states_names = c('Mainland', 'Island'),
  colors = c("#d72027", "#2b83bb"))

ggsave(zheng_rates_poe.tree_plot, 
       filename = 'Zheng-Spp_PoeTree-HisseTipRates.pdf')





