##################################################
# Testing for convergence a la Mazel et al. 2017 #
##################################################


# Now for the mainland groups (M1, M2-allopatric, M2-sympatric)

######### Fit models of trait evolution ##########
library(treeplyr)
library(purrr)
library(geiger)
library(mvMORPH)
library(plyr)
library(RRphylo)
library(adephylo)
library(spaa)
library(otuSummary)
library(igraph)
library(scales)
library(plyr)
library(tidyr)
library(stringr)
library(ggplot2)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Convergence/')
source('../ConvergenceDivergenceFunctions.R')

m.data <- read.csv('../FinalData/Male_Transf-Scaled-AllTraitsEcol-FINAL.csv', 
                   header = T)
f.data <- read.csv('../FinalData/Female_Transf-Scaled-AllTraitsEcol-FINAL.csv', 
                   header = T)

# Remove body height as this leads to a substantial reduction in 
# the number of species included. For these convergence analyses,
# we want to ensure we've got as complete a dataset as possible to
# make sure that we've got as many of the 'classic' island ecomorphs
# as possible. We will also remove agasizzi, since this species is a 
# bit obscure for a number of reasons (trees)
m.data <- m.data[,c(1,3:8,10,12:13,23:24)]
m.data <- m.data[-which(m.data$Species %in% c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 'medemi')),]

# We also need to remove lamella counts for females as many species
# are missing data for these traits (for females)
f.data <- f.data[,c(1,3:10,23:24)]
f.data <- f.data[-which(f.data$Species %in% c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 'medemi')),]

# Now remove missing data
m.data <- na.omit(m.data)
f.data <- na.omit(f.data)

# Read in the time-calibrated phylogeny
tree <- read.tree('../FinalData/Anole_TimeCal_NoOut.tree')

# Read in the PSTD results 

m.pstd <- readRDS('Male_mvBM_PSTD.rds')
f.pstd <- readRDS('Female_mvOU_PSTD.rds')

###
###
###
###
###


patry <- read.csv('../FinalData/Allopatric_wRespectTo_Dacty.csv')
m1.patry <- patry[which(patry$Clade =='M1'),]
m2.patry <- patry[-which(patry$Clade =='M1'),]

allopat <- m2.patry[which(m2.patry$AllopatricVersusDacty == 1),]
sympat <- m2.patry[which(m2.patry$AllopatricVersusDacty == 0),]
regions <- c(rep('M1', nrow(m1.patry)), 
             rep('M2-Allopatric', nrow(allopat)), 
             rep('M2-Sympatric', nrow(sympat)))

patry.spp <- c(as.character(m1.patry$Species), 
               as.character(allopat$Species), 
               as.character(sympat$Species))

spp.regions <- data.frame('Species' = patry.spp, 'Region' = regions)
spp.regions.m <- spp.regions[which(spp.regions$Species %in% m.full$tip.label),]
spp.regions.f <- spp.regions[which(spp.regions$Species %in% f.full$tip.label),]


# Now prepare to identify actual convergent/divergent
# species pairs. We will correspond PSTD values and species
# pairs to regions
m.pstd.results <- PreparePSTD(pstd = m.pstd, spp.cats = spp.regions.m)
f.pstd.results <- PreparePSTD(pstd = f.pstd, spp.cats = spp.regions.f)


# Prepare the PSTD values such that they are useful
m.dist.list <- PreparePSTD(pstd = m.pstd, spp.cats = spp.regions.m)
f.dist.list <- PreparePSTD(pstd = f.pstd, spp.cats = spp.regions.f)

# And get the set of convergent/divergent species pairs
m.bott <- GetConvergent(pstd.dists = m.dist.list, 
                        threshold = 0.01, 
                        tree = tree, minDist = 20,
                        dist.type = 'PatristicDist')
m.top <- GetDivergent(pstd.dists = m.dist.list, 
                      threshold = 0.01, 
                      tree = tree)

f.bott <- GetConvergent(pstd.dists = f.dist.list, 
                        threshold = 0.01, 
                        tree = tree, minDist = 20,
                        dist.type = 'PatristicDist')
f.top <- GetDivergent(pstd.dists = f.dist.list, 
                      threshold = 0.01, 
                      tree = tree)

# And summarize these results for plotting as a bar chart. 
m.summs <- SummarizeConvDiv(conv = m.bott, div = m.top, 
                            cats = c('M1', 'M2-Allopatric', 'M2-Sympatric'),
                            dist.list = m.dist.list)
f.summs <- SummarizeConvDiv(conv = f.bott, div = f.top, 
                            cats = c('M1', 'M2-Allopatric', 'M2-Sympatric'),
                            dist.list = f.dist.list)


setwd('Mainland')
capture.output(m.summs, file = 'ConvDiv_mvBM_ByGroup_Scaled_Males-Mainland.txt')
capture.output(f.summs, file = 'ConvDiv_mvOU_ByGroup_Scaled_Females-Mainland.txt')

pdf('ConvDiv_mvBM_Mainland_Scaled_Males.pdf', width = 8, height = 4.66)
m.main.bars <- 
  ggplot(data = m.summs$median, aes(x = Category, y = ObsOverExpPrev)) + 
  geom_bar(stat = "identity", aes(fill = Type), 
           position=position_dodge()) + 
  scale_fill_manual(values = c('#67a9cf', '#ef8a62')) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_bw()
m.main.bars
dev.off()

pdf('ConvDiv_mvOU_Mainland_Scaled_Females.pdf', width=7, height=6)
f.main.bars <- 
  ggplot(data = f.summs$median, aes(x = Category, y = ObsOverExpPrev)) + 
  geom_bar(stat = "identity", aes(fill = Type), 
           position=position_dodge()) + 
  scale_fill_manual(values = c('#67a9cf', '#ef8a62')) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_bw()
f.main.bars
dev.off()


# First we need to define the possible sets of pairs. Because the order can flip, 
# we need to account for this. Likewise, we need to define factor 
# levels for the pairs
pair.types <- list(c('M1--M2-Allopatric', 'M2-Allopatric--M1'), 
                   c('M1--M2-Sympatric', 'M2-Sympatric--M1'),
                   c('M2-Allopatric--M2-Sympatric', 'M2-Sympatric--M2-Allopatric'))

pair.levels <- c('M1--M1', 'M2-Allopatric--M2-Allopatric', 
                 'M2-Sympatric--M2-Sympatric', 'M2-Allopatric--M1', 
                 'M2-Sympatric--M1', 'M2-Sympatric--M2-Allopatric')


Convergent <- m.bott 
Divergent <- m.top
pstdResults <- m.pstd.results
m <- 'median'
threshold <- 0.01


m.heat <-
  PrepareHeatmap(Convergent = m.bott, Divergent = m.top,
                 pstdResults = m.pstd.results, 
                 pair.types = pair.types,
                 pair.levels = pair.levels,
                 threshold = 0.01)
saveRDS(m.heat, file = 'PairwiseConvDiv-Males-Mainland.rds')
capture.output(m.heat, file = 'PairwiseConvDiv-Males-Mainland.txt')


f.heat <-
  PrepareHeatmap(Convergent = f.bott, Divergent = f.top,
                 pstdResults = f.pstd.results, 
                 pair.types = pair.types,
                 pair.levels = pair.levels,
                 threshold = 0.01)
saveRDS(f.heat, file = 'PairwiseConvDiv-Females-Mainland.rds')
capture.output(f.heat, file = 'PairwiseConvDiv-Females-Mainland.txt')


pdf('Male_Convergences-Median_Mainland-Heatmap.pdf', width = 7, height = 7)
m.conv.heat.main <- 
  ggplot(data = m.heat$ObsExp, aes(x = Region1, y = Region2, 
                                   fill = log(ConvMed.ObsExpFreqs))) +
  geom_tile() +
  scale_fill_distiller(type = 'div', palette = 7, limits = c(-1, 1)) +
  guides(fill=guide_legend(title="log(Obs/Exp Convergences)")) +
  theme_bw()  +
  theme(legend.position = c(0.025, 0.975), 
        legend.justification = c(0, 1))
m.conv.heat.main
dev.off()

pdf('Male_Divergences-Median_Mainland-Heatmap.pdf', width = 7, height = 7)
m.div.heat.main <- 
  ggplot(data = m.heat$ObsExp, aes(x = Region1, y = Region2, 
                                   fill = log(DivMed.ObsExpFreqs))) +
  geom_tile() +
  scale_fill_distiller(type = 'div', palette = 7, limits = c(-1, 1)) +
  guides(fill=guide_legend(title="log(Obs/Exp Divergences)")) +
  theme_bw()  +
  theme(legend.position = c(0.025, 0.975), 
        legend.justification = c(0, 1)) 
m.div.heat.main
dev.off()


pdf('Female_Convergences-Median_Mainland-Heatmap.pdf', width = 7, height = 7)
f.conv.heat.main <- 
  ggplot(data = f.heat$ObsExp, aes(x = Region1, y = Region2, 
                                   fill = log(ConvMed.ObsExpFreqs))) +
  geom_tile() +
  scale_fill_distiller(type = 'div', palette = 7, limits = c(-0.75, 0.75)) +
  guides(fill=guide_legend(title="log(Obs/Exp Convergences)")) +
  theme_bw()  +
  theme(legend.position = c(0.025, 0.975), 
        legend.justification = c(0, 1)) 
f.conv.heat.main
dev.off()

pdf('Female_Divergences-Median_Mainland-Heatmap.pdf', width = 7, height = 7)
f.div.heat.main <- 
  ggplot(data = f.heat$ObsExp, aes(x = Region1, y = Region2, 
                                   fill = log(DivMed.ObsExpFreqs))) +
  geom_tile() +
  scale_fill_distiller(type = 'div', palette = 7, limits = c(-2, 2)) +
  guides(fill=guide_legend(title="log(Obs/Exp Divergences)")) +
  theme_bw()  +
  theme(legend.position = c(0.025, 0.975), 
        legend.justification = c(0, 1)) 
f.div.heat.main
dev.off()

library(patchwork)

lay <- "
AABBEE
AABBEE
CCDDFF
CCDDFF
"

pdf('AllMainland-ConvDiv-Plots.pdf', width = 20, height = 14)
m.conv.heat.main + m.div.heat.main + f.conv.heat.main + f.div.heat.main + 
  m.main.bars + f.main.bars + 
  plot_layout(design = lay)
dev.off()







