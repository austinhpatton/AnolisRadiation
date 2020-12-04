require(ape)
setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/Hisse')

# Remove tree and remove outgroup
tree <- read.tree('../FinalData/Anole_TimeCal_NoOut.tree')
tree <- drop.tip(tree, c('agassizi', 'bicaorum', 'lineatus', 'roatanensis',
                         'utilensis', 'pinchoti', 'concolor', 'townsendi'))
# Read in data (treating mainland/island as a binary) and remove outgroup. Format for analysis in hisse (vector of states, species as names)
dat <- read.csv('../FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv')
dat <- data.frame('Species' = dat$Species, 'Island' = rowSums(dat[,33:38]))
dat <- dat[which(dat$Species %in% tree$tip.label),]



## BiSSE Models

library(hisse)
library(ggplot2)
library(reshape)
###############
# BiSSE Model #
###############
# Make Trasition matrix
trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

# Run the model
bisse <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = F, turnover.anc = c(1,2,0,0),
               eps.anc = c(1,2,0,0), trans.rate = trans.rates.bisse, output.type = "raw")

# Get support region by sampling likelihood surface
# recon_bisse <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=bisse$solution, hidden.states = F, aic=bisse$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

# Save results
save(bisse, file = 'bisse_Object-70percMain.Rsave')
# save(recon_bisse, file='recon_bisse-70percMain.Rsave')

# Store likelihood, AIC, and AICc
bisse.logL <- bisse$loglik
bisse.AIC <- bisse$AIC
bisse.AICc <- bisse$AICc

#################################################################
#                   Bisse Null Model - Unequal Trans            #
#################################################################
bisse.null.uneq <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = F, turnover.anc = c(1,1,0,0),
                         eps.anc = c(1,1,0,0), trans.rate = trans.rates.bisse, output.type = "raw")

recon_bisse.null.uneq <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=bisse.null.uneq$solution, hidden.states = F, aic=bisse.null.uneq$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(bisse.null.uneq, file = 'bisse.null.uneq_Object-70percMain.Rsave')
save(recon_bisse.null.uneq, file='recon_bisse.null.uneq-70percMain.Rsave')

bisse.null.uneq.logL <- bisse.null.uneq$loglik
bisse.null.uneq.AIC <- bisse.null.uneq$AIC
bisse.null.uneq.AICc <- bisse.null.uneq$AICc

#################################################################
#                   Bisse Null Model - Equal Trans              #
#################################################################
trans.rate <- TransMatMaker(hidden.states=FALSE)
trans.rate <- ParEqual(trans.rate, c(1,2))

bisse.null.eq <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = F, turnover.anc = c(1,1,0,0),
                       eps.anc = c(1,1,0,0), trans.rate = trans.rate, output.type = "raw")

recon_bisse.null.eq <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=bisse.null.eq$solution, hidden.states = F, aic=bisse.null.eq$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(bisse.null.eq, file = 'bisse.null.eq_Object-70percMain.Rsave')
save(recon_bisse.null.eq, file='recon_bisse.null.eq-70percMain.Rsave')

bisse.null.eq.logL <- bisse.null.eq$loglik
bisse.null.eq.AIC <- bisse.null.eq$AIC
bisse.null.eq.AICc <- bisse.null.eq$AICc

## HiSSE Models
##############################
# Full HiSSE All Transitions #
##############################
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)

hisse.full_AllTrans <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,2,3,4),
                             eps.anc = c(1,2,3,4), trans.rate = trans.rates.hisse, output.type = "raw")

recon_hisse.full_AllTrans <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=hisse.full_AllTrans$solution, hidden.states = T, aic=hisse.full_AllTrans$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(hisse.full_AllTrans, file = 'hisse.full_AllTrans_Object-70percMain.Rsave')
save(recon_hisse.full_AllTrans, file='recon_hisse.full_AllTrans-70percMain.Rsave')

hisse.full_AllTrans.logL <- hisse.full_AllTrans$loglik
hisse.full_AllTrans.AIC <- hisse.full_AllTrans$AIC
hisse.full_AllTrans.AICc <- hisse.full_AllTrans$AICc

####################################
# Full HiSSE No Double Transitions #
####################################
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))

hisse.full_NoDub <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,2,3,4),
                          eps.anc = c(1,2,3,4), trans.rate = trans.rates.hisse, output.type = "raw")

recon_hisse.full_NoDub <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=hisse.full_NoDub$solution, hidden.states = T, aic=hisse.full_NoDub$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(hisse.full_NoDub, file = 'hisse.full_NoDub_Object-70percMain.Rsave')
save(recon_hisse.full_NoDub, file='recon_hisse.full_NoDub-70percMain.Rsave')

hisse.full_NoDub.logL <- hisse.full_NoDub$loglik
hisse.full_NoDub.AIC <- hisse.full_NoDub$AIC
hisse.full_NoDub.AICc <- hisse.full_NoDub$AICc

####################################
# Null 2 All Transitions All Equal #
####################################
trans.rate <- TransMatMaker(hidden.states=TRUE)
trans.rate[!is.na(trans.rate) & !trans.rate == 0] = 1

Null2_AllTrans_Eq <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,1,2,2),
                           eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw")

recon_Null2_AllTrans_Eq <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=Null2_AllTrans_Eq$solution, hidden.states = T, aic=Null2_AllTrans_Eq$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(Null2_AllTrans_Eq, file = 'Null2_AllTrans_Eq_Object-70percMain.Rsave')
save(recon_Null2_AllTrans_Eq, file='recon_Null2_AllTrans_Eq-70percMain.Rsave')

Null2_AllTrans_Eq.logL <- Null2_AllTrans_Eq$loglik
Null2_AllTrans_Eq.AIC <- Null2_AllTrans_Eq$AIC
Null2_AllTrans_Eq.AICc <- Null2_AllTrans_Eq$AICc

##########################################
# Null 2 No Double Transitions All Equal #
##########################################
trans.rate <- TransMatMaker(hidden.states=TRUE)
trans.rate <- ParDrop(trans.rate, c(3,5,8,10))
trans.rate[!is.na(trans.rate) & !trans.rate == 0] = 1

Null2_NoDub_Eq <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,1,2,2),
                        eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw")

recon_Null2_NoDub_Eq <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=Null2_NoDub_Eq$solution, hidden.states = T, aic=Null2_NoDub_Eq$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(Null2_NoDub_Eq, file = 'Null2_NoDub_Eq_Object-70percMain.Rsave')
save(recon_Null2_NoDub_Eq, file='recon_Null2_NoDub_Eq-70percMain.Rsave')

Null2_NoDub_Eq.logL <- Null2_NoDub_Eq$loglik
Null2_NoDub_Eq.AIC <- Null2_NoDub_Eq$AIC
Null2_NoDub_Eq.AICc <- Null2_NoDub_Eq$AICc

#######################################
# Null 2 Three Trans Rates, No Double #
#######################################

# Define transition matrix
trans.rate <- TransMatMaker(hidden.states=TRUE)
trans.rate <- ParDrop(trans.rate, c(3,5,8,10))
to.change <- cbind(c(1,3), c(2,4))
trans.rate[to.change] = 1
# Now set all transitions from 1->0 to be governed by a single rate:
to.change <- cbind(c(2,4), c(1,3))
trans.rate[to.change] = 2
# Finally, set all transitions between the hidden state to be a single rate (essentially giving 
# you an estimate of the rate by which shifts in diversification occur:
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rate[to.change] = 3
trans.rate

Null2_ThreeRate_NoDub <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,1,2,2),
                               eps.anc = c(1,1,2,2), trans.rate = trans.rate, output.type = "raw")

recon_Null2_ThreeRate_NoDub <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=Null2_ThreeRate_NoDub$solution, hidden.states = T, aic=Null2_ThreeRate_NoDub$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(Null2_ThreeRate_NoDub, file = 'Null2_ThreeRate_NoDub_Object-70percMain.Rsave')
save(recon_Null2_ThreeRate_NoDub, file='recon_Null2_ThreeRate_NoDub-70percMain.Rsave')

Null2_ThreeRate_NoDub.logL <- Null2_ThreeRate_NoDub$loglik
Null2_ThreeRate_NoDub.AIC <- Null2_ThreeRate_NoDub$AIC
Null2_ThreeRate_NoDub.AICc <- Null2_ThreeRate_NoDub$AICc

###################################
# HiSSE Null 4 Model, Equal Rates #
###################################
hisse.null4.equal <- hisse.null4(phy = tree, data = dat, f = c(0.70, 0.85),  turnover.anc = rep(c(1,2,3,4),2),
                                 eps.anc = rep(c(1,2,3,4),2), trans.type = "equal", output.type = "raw")

recon_hisse.null4.equal <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=hisse.null4.equal$solution, hidden.states = T, aic=hisse.null4.equal$AICc, n.cores = 2, four.state.null = T, root.p = c(0.5, 0, 0.5, 0))

save(hisse.null4.equal, file = 'hisse.null4.equal_Object-70percMain.Rsave')
save(recon_hisse.null4.equal, file='recon_hisse.null4.equal-70percMain.Rsave')

hisse.null4.equal.logL <- hisse.null4.equal$loglik
hisse.null4.equal.AIC <- hisse.null4.equal$AIC
hisse.null4.equal.AICc <- hisse.null4.equal$AICc

###################################
# HiSSE Null 4 Model, Three Rates #
###################################
hisse.null4.three <- hisse.null4(phy = tree, data = dat, f = c(0.70, 0.85),  turnover.anc = rep(c(1,2,3,4),2),
                                 eps.anc = rep(c(1,2,3,4),2), trans.type = "three.rate", output.type = "raw")

recon_hisse.null4.three <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=hisse.null4.three$solution, hidden.states = T, aic=hisse.null4.three$AICc, n.cores = 2, four.state.null = T, root.p = c(0.5, 0, 0.5, 0))

save(hisse.null4.three, file = 'hisse.null4.three_Object-70percMain.Rsave')
save(recon_hisse.null4.three, file='recon_hisse.null4.three-70percMain.Rsave')

hisse.null4.three.logL <- hisse.null4.three$loglik
hisse.null4.three.AIC <- hisse.null4.three$AIC
hisse.null4.three.AICc <- hisse.null4.three$AICc

##############################
# HiSSE No0B All Transitions #
##############################
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(2,5,12,7,8,9))

no0B_AllTrans <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,2,0,3),
                       eps.anc = c(1,2,0,3), trans.rate = trans.rates.hisse, output.type = "raw")

recon_no0B_AllTrans <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=no0B_AllTrans$solution, hidden.states = T, aic=no0B_AllTrans$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(no0B_AllTrans, file = 'no0B_AllTrans_Object-70percMain.Rsave')
save(recon_no0B_AllTrans, file='recon_no0B_AllTrans-70percMain.Rsave')

no0B_AllTrans.logL <- no0B_AllTrans$loglik
no0B_AllTrans.AIC <- no0B_AllTrans$AIC
no0B_AllTrans.AICc <- no0B_AllTrans$AICc

####################################
# HiSSE No0B No Double Transitions #
####################################
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(2,3,5,7,8,9,10,12))

no0B_NoDub <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,2,0,3),
                    eps.anc = c(1,2,0,3), trans.rate = trans.rates.hisse, output.type = "raw")

recon_no0B_NoDub <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=no0B_NoDub$solution, hidden.states = T, aic=no0B_NoDub$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(no0B_NoDub, file = 'no0B_NoDub_Object-70percMain.Rsave')
save(recon_no0B_NoDub, file='recon_no0B_NoDub-70percMain.Rsave')

no0B_NoDub.logL <- no0B_NoDub$loglik
no0B_NoDub.AIC <- no0B_NoDub$AIC
no0B_NoDub.AICc <- no0B_NoDub$AICc

##############################
# HiSSE No1B All Transitions #
##############################
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(3,6,9,10,11,12))

no1B_AllTrans <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,2,3,0),
                       eps.anc = c(1,2,3,0), trans.rate = trans.rates.hisse, output.type = "raw")

recon_no1B_AllTrans <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=no1B_AllTrans$solution, hidden.states = T, aic=no1B_AllTrans$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(no1B_AllTrans, file = 'no1B_AllTrans_Object-70percMain.Rsave')
save(recon_no1B_AllTrans, file='recon_no1B_AllTrans-70percMain.Rsave')

no1B_AllTrans.logL <- no1B_AllTrans$loglik
no1B_AllTrans.AIC <- no1B_AllTrans$AIC
no1B_AllTrans.AICc <- no1B_AllTrans$AICc

####################################
# HiSSE No1B No Double Transitions #
####################################
trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
trans.rates.hisse <- ParDrop(trans.rates.hisse, drop.par=c(3,5,6,8,9,10,11,12))

no1B_NoDub <- hisse(phy = tree, data = dat, f = c(0.70, 0.85), hidden.states = T, turnover.anc = c(1,2,3,0),
                    eps.anc = c(1,2,3,0), trans.rate = trans.rates.hisse, output.type = "raw")

recon_no1B_NoDub <- MarginRecon(tree, dat, f = c(0.70, 0.85), pars=no1B_NoDub$solution, hidden.states = T, aic=no1B_NoDub$AICc, n.cores = 2, root.p = c(0.5, 0, 0.5, 0))

save(no1B_NoDub, file = 'no1B_NoDub_Object-70percMain.Rsave')
save(recon_no1B_NoDub, file='recon_no1B_NoDub-70percMain.Rsave')

no1B_NoDub.logL <- no1B_NoDub$loglik
no1B_NoDub.AIC <- no1B_NoDub$AIC
no1B_NoDub.AICc <- no1B_NoDub$AICc


# Summarize Model Support
require(qpcR)

logL <- c(Null2_AllTrans_Eq.logL, Null2_NoDub_Eq.logL, Null2_ThreeRate_NoDub.logL, bisse.logL, bisse.null.eq.logL, bisse.null.uneq.logL, hisse.null4.equal.logL, hisse.null4.three.logL, no0B_AllTrans.logL, no0B_NoDub.logL, no1B_AllTrans.logL, no1B_NoDub.logL, hisse.full_AllTrans.logL, hisse.full_NoDub.logL)
AIC <- c(Null2_AllTrans_Eq.AIC, Null2_NoDub_Eq.AIC, Null2_ThreeRate_NoDub.AIC, bisse.AIC, bisse.null.eq.AIC, bisse.null.uneq.AIC, hisse.null4.equal.AIC, hisse.null4.three.AIC, no0B_AllTrans.AIC, no0B_NoDub.AIC, no1B_AllTrans.AIC, no1B_NoDub.AIC, hisse.full_AllTrans.AIC, hisse.full_NoDub.AIC)
AICc <- c(Null2_AllTrans_Eq.AICc, Null2_NoDub_Eq.AICc, Null2_ThreeRate_NoDub.AICc, bisse.AICc, bisse.null.eq.AICc, bisse.null.uneq.AICc, hisse.null4.equal.AICc, hisse.null4.three.AICc, no0B_AllTrans.AICc, no0B_NoDub.AICc, no1B_AllTrans.AICc, no1B_NoDub.AICc, hisse.full_AllTrans.AICc, hisse.full_NoDub.AICc)

results <- as.data.frame(logL, row.names=c('Null2 AllTrans EqRates', 'Null2 NoDub EqRates', 'Null2 ThreeRate NoDub', 'BiSSE', 'BiSSE Null EqTrans', 'BiSSE Null UneqTrans', 'Null4 EqRates', 'Null4 ThreeRate', 'No0B AllTrans', 'No0B NoDub',  'No1B AllTrans', 'No1B NoDub', 'HiSSE Full AllTrans', 'HiSSE Full NoDub'))
results$AIC <- AIC
results$AICc <- AICc

ak.weights <- akaike.weights(AICc)
results$Delta.AIC <- ak.weights$deltaAIC
results$Rel.Log.Lik <- ak.weights$rel.LL
results$Akaike.Weights <- ak.weights$weights
write.csv(results, 'ModelSupport_HiSSE-Assume70percMainSpp.csv')

library(hisse)
#colfunc<-colorRampPalette(c('darkblue', 'royalblue', 'green2', 'yellow2','orange','red'))
colfunc<-colorRampPalette(c('#313695', '#4575b4', '#74add1', 
                            '#abd9e9', '#e0f3f8', '#ffffbf', 
                            '#fee090', '#fdae61', '#f46d43', 
                            '#d73027', '#a50026', '#a50026'))

hisse.list <- list('recon_Null2_AllTrans_Eq' = recon_Null2_AllTrans_Eq, 'recon_Null2_NoDub_Eq' = recon_Null2_NoDub_Eq, 
                   'recon_Null2_ThreeRate_NoDub' = recon_Null2_ThreeRate_NoDub, 'recon_hisse.null4.equal' = recon_hisse.null4.equal, 
                   'recon_hisse.null4.three' = recon_hisse.null4.three, 'recon_no0B_AllTrans' = recon_no0B_AllTrans,
                   'recon_no0B_NoDub' = recon_no0B_NoDub, 'recon_no1B_AllTrans' = recon_no1B_AllTrans, 'recon_no1B_NoDub' =
                     recon_no1B_NoDub, 'recon_hisse.full_AllTrans' = recon_hisse.full_AllTrans, 'recon_hisse.full_NoDub' =
                     recon_hisse.full_NoDub)

pdf('ModAvg_RateTree_NetDiv-Assume70percMainSpp.pdf', height = 7, width=8)
plot.hisse.states(hisse.list, rate.param = 'net.div', do.observed.only = T, rate.colors = colfunc(1000), 
                  edge.width.state = 1, show.tip.label = F, legend.kernel.rates = 'gaussian', 
                  edge.width.rate = 3.5, legend = 'all', legend.cex = 0.8, type = 'phylogram')
dev.off()

hisse.contMap <- plot.hisse.states(hisse.list, rate.param = 'net.div', do.observed.only = T, rate.colors = colfunc(1000), 
                                   edge.width.state = 1, show.tip.label = F, legend.kernel.rates = 'gaussian', 
                                   edge.width.rate = 3.5, legend = 'all', legend.cex = 0.8)

pdf('ModAvg_RateTree_NetDiv-Assume70percMainSpp-ForFig.pdf', height = 10, width=7)
plot(hisse.contMap$rate.tree, outline = F, lwd = 3)
dev.off()

pdf('AncStateRecon-BestFitHisse-Assume70percMainSpp.pdf', width = 8, height = 18)
plot(tree, show.tip.label = T, cex = 0.45, label.offset = 0.5)
tiplabels(pie = recon_no1B_NoDub$tip.mat[,-c(1,5)], piecol = c('#66c2a5', '#fc8d62', '#8da0cb'), cex = 0.2)
nodelabels(pie = recon_no1B_NoDub$node.mat[,-c(1,5)], piecol = c('#66c2a5', '#fc8d62', '#8da0cb'), cex = 0.3)
legend("bottomleft", inset=.02, title="Character State",
       c('Mainland-A', 'Island', 'Mainland-B'), fill=c('#66c2a5', '#fc8d62', '#8da0cb'), horiz=F, cex=1)
dev.off()



