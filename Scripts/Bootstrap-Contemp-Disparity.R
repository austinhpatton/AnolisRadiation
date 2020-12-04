library(dispRity)
library(cowplot)
library(ape)
library(phytools)
library(pairwiseCI)
library(geometry)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')
dat <- read.csv('./FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL-SUBMIT.csv', header = T, row.names = 1)
tree <- read.tree('./FinalData/Anole_TimeCal_NoOut.tree')

dat <- data[,c(2:9,11:12,22)]
eco <- data[,c(14,17,22)]

dat[,1:10] <- apply(X = dat[,1:10], FUN = scale, MARGIN = 2)
eco[,1:2] <- apply(X = eco[,1:2], FUN = scale, MARGIN = 2)

dat <- na.omit(dat)
eco <- na.omit(eco)

dispRity.per.group <-
  function (data, group, metric = c(median, centroids), ...) 
  {
    data_subsetd <- custom.subsets(data, group, ...)
    data_bootstrapped <- boot.matrix(data_subsetd, bootstraps = 1000)
    return(dispRity(data_bootstrapped, metric = metric, ...))
  }

morph.disp.per.group <-
  dispRity.per.group(data = as.matrix(dat[,-11]), 
                     list('I2' = which(dat$Group == 'I2'),
                          'M1' = which(dat$Group == 'M1'),
                          'M2' = which(dat$Group == 'M2')),
                     metric = pairwise.dist)

eco.disp.per.group <-
  dispRity.per.group(data = as.matrix(eco[,-3]), 
                     list('I2' = which(eco$Group == 'I2'),
                          'M1' = which(eco$Group == 'M1'),
                          'M2' = which(eco$Group == 'M2')),
                     metric = pairwise.dist)


morph.disp.broad <-
  dispRity.per.group(data = as.matrix(dat[,-11]), 
                     list('I2' = which(dat$Group == 'I2'),
                          'M' = which(dat$Group %in% c('M1', 'M2'))),
                     metric = pairwise.dist)

eco.disp.broad <-
  dispRity.per.group(data = as.matrix(eco[,-3]), 
                     list('I2' = which(eco$Group == 'I2'),
                          'M' = which(eco$Group %in% c('M1', 'M2'))),
                     metric = pairwise.dist)

avgSq <- function(x){
  dist <- mean(x^2)
  return(dist)
}

i2.disp <- apply(morph.disp.broad$disparity$I2[[2]], FUN = mean, MARGIN = 2)
m.disp <- apply(morph.disp.broad$disparity$M[[2]], FUN = mean, MARGIN = 2)
m1.disp <- apply(morph.disp.per.group$disparity$M1[[2]], FUN = mean, MARGIN = 2)
m2.disp <- apply(morph.disp.per.group$disparity$M2[[2]], FUN = mean, MARGIN = 2)

i2.eco.disp <- apply(eco.disp.broad$disparity$I2[[2]], FUN = mean, MARGIN = 2)
m.eco.disp <- apply(eco.disp.broad$disparity$M[[2]], FUN = mean, MARGIN = 2)
m1.eco.disp <- apply(eco.disp.per.group$disparity$M1[[2]], FUN = mean, MARGIN = 2)
m2.eco.disp <- apply(eco.disp.per.group$disparity$M2[[2]], FUN = mean, MARGIN = 2)

i2.obs.morph <- apply(morph.disp.broad$disparity$I2[[1]], FUN = mean, MARGIN = 2)
m.obs.morph <- apply(morph.disp.broad$disparity$M[[1]], FUN = mean, MARGIN = 2)
m1.obs.morph <- apply(morph.disp.per.group$disparity$M1[[1]], FUN = mean, MARGIN = 2)
m2.obs.morph <- apply(morph.disp.per.group$disparity$M2[[1]], FUN = mean, MARGIN = 2)

i2.obs.eco <- apply(eco.disp.broad$disparity$I2[[1]], FUN = mean, MARGIN = 2)
m.obs.eco <- apply(eco.disp.broad$disparity$M[[1]], FUN = mean, MARGIN = 2)
m1.obs.eco <- apply(eco.disp.per.group$disparity$M1[[1]], FUN = mean, MARGIN = 2)
m2.obs.eco <- apply(eco.disp.per.group$disparity$M2[[1]], FUN = mean, MARGIN = 2)

bs.disp.morph <-
  data.frame('Group' = c(rep('GA', 1000), rep('Mainland', 1000),
                         rep('M1', 1000), rep('M2', 1000)),
             'Traits' = rep('Morphology', 4000),
             'Morphological Disparity' = c(i2.disp, m.disp, 
                             m1.disp, m2.disp))

bs.disp.eco <-
  data.frame('Group' = c(rep('GA', 1000), rep('Mainland', 1000),
                         rep('M1', 1000), rep('M2', 1000)),
             'Traits' = rep('Ecology', 4000),
             'Ecological Disparity' = c(i2.eco.disp, m.eco.disp, 
                             m1.eco.disp, m2.eco.disp))

obs.disp <- 
  data.frame('Group' = c('GA', 'Mainland', 'M1', 'M2'),
             'Traits' = rep('Morphology', 8),
             'Morphological Disparity' = c(i2.obs.morph, m.obs.morph, 
                             m1.obs.morph, m2.obs.morph))

obs.eco <- 
  data.frame('Group' = c('GA', 'Mainland', 'M1', 'M2'),
             'Traits' = rep('Ecology', 8),
             'Ecological Disparity' = c(i2.obs.eco, m.obs.eco, 
                             m1.obs.eco, m2.obs.eco))

morph.disp.broad <-
  data.frame('Group' = c('GA', 'Mainland'), 
             'Morphological Disparity' = c(i2.obs.morph, m.obs.morph),
             'Lower' = c(quantile(i2.disp, 0.025),
                         quantile(m.disp, 0.025)),
             'Upper' = c(quantile(i2.disp, 0.925),
                         quantile(m.disp, 0.925)))

morph.disp <-
  data.frame('Group' = c('GA', 'M1', 'M2'), 
             'Morphological Disparity' = c(i2.obs.morph,
                             m1.obs.morph, m2.obs.morph),
             'Lower' = c(quantile(i2.disp, 0.025),
                         quantile(m1.disp, 0.025),
                         quantile(m2.disp, 0.025)),
             'Upper' = c(quantile(i2.disp, 0.925),
                         quantile(m1.disp, 0.925),
                         quantile(m2.disp, 0.925)))

eco.disp.broad <-
  data.frame('Group' = c('GA', 'Mainland'), 
             'Ecological Disparity' = c(i2.obs.eco, m.obs.eco),
             'Lower' = c(quantile(i2.eco.disp, 0.025),
                         quantile(m.eco.disp, 0.025)),
             'Upper' = c(quantile(i2.eco.disp, 0.925),
                         quantile(m.eco.disp, 0.925)))

eco.disp <-
  data.frame('Group' = c('GA', 'M1', 'M2'), 

             'Ecological Disparity' = c(i2.obs.eco,
                             m1.obs.eco, m2.obs.eco),
             'Lower' = c(quantile(i2.eco.disp, 0.025),
                         quantile(m1.eco.disp, 0.025),
                         quantile(m2.eco.disp, 0.025)),
             'Upper' = c(quantile(i2.eco.disp, 0.925),
                         quantile(m1.eco.disp, 0.925),
                         quantile(m2.eco.disp, 0.925)))

i2.m.pv <- length(which(c(outer(i2.disp, m.disp, '-')) > 0)) / 1000000
i2.m1.pv <- length(which(c(outer(i2.disp, m1.disp, '-')) > 0)) / 1000000
i2.m2.pv <- length(which(c(outer(i2.disp, m2.disp, '-')) > 0)) / 1000000
m1.m2.pv <- length(which(c(outer(m1.disp, m2.disp, '-')) > 0)) / 1000000

i2.m.eco.pv <- length(which(c(outer(i2.eco.disp, m.eco.disp, '-')) > 0)) / 1000000
i2.m1.eco.pv <- length(which(c(outer(i2.eco.disp, m1.eco.disp, '-')) > 0)) / 1000000
i2.m2.eco.pv <- length(which(c(outer(i2.eco.disp, m2.eco.disp, '-')) > 0)) / 1000000
m1.m2.eco.pv <- length(which(c(outer(m1.eco.disp, m2.eco.disp, '-')) > 0)) / 1000000

pvals <- c(i2.m.pv, i2.m1.pv, i2.m2.pv, m1.m2.pv,
           i2.m.eco.pv, i2.m1.eco.pv, i2.m2.eco.pv, m1.m2.eco.pv)
for(i in 1:length(pvals)){
  if(pvals[i] > 0.5){
    pvals[i] <- 1 - pvals[i]
  }
}

pvals <- pvals * 2
  
pvals <- 
  matrix(data = pvals,
         nrow = 2, ncol = 4, byrow = T,
         dimnames = list(c('Morph', 'Eco'),
                         c('I2-M', 'I2-M1', 'I2-M2', 'M1-M2')))

write.table(pvals, file = 'Bootstrapped-Contemp-Disparity-Pvals.txt', 
            sep = '\t', quote = F)


cols <- c("#d72027", "#add8a4", "#2b83bb")

compares <- list('GA', 'Mainland')

broad.morph.disp.p <-
  ggplot(data = morph.disp.broad, aes(x = Group, y = Morphological.Disparity, fill = Group)) +
  geom_point(data = bs.disp.morph[which(bs.disp.morph$Group %in% c('GA', 'Mainland')),], 
             aes(x = Group, y = Morphological.Disparity, fill = Group, color = Group),
             position = position_jitter(width = 0.1), alpha = 0.15) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), size = 1, width = 0.3) +
  geom_point(aes(fill = Group), pch = 21, size = 5, stroke = 1) +
  scale_fill_manual(values = cols[c(1,3)]) +
  scale_color_manual(values = cols[c(1,3)]) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_blank(), 
        legend.position = 'none')

broad.eco.disp.p <-
  ggplot(data = eco.disp.broad, aes(x = Group, y = Ecological.Disparity, fill = Group)) +
  geom_point(data = bs.disp.eco[which(bs.disp.eco$Group %in% c('GA', 'Mainland')),], 
             aes(x = Group, y = Ecological.Disparity, fill = Group, color = Group),
             position = position_jitter(width = 0.1), alpha = 0.15) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), size = 1, width = 0.3) +
  geom_point(aes(fill = Group), pch = 21, size = 5, stroke = 1) +
  scale_fill_manual(values = cols[c(1,3)]) +
  scale_color_manual(values = cols[c(1,3)]) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_blank(), 
        legend.position = 'none')

morph.disp.p <-
  ggplot(data = morph.disp, aes(x = Group, y = Morphological.Disparity, fill = Group)) +
  geom_point(data = bs.disp.morph[-which(bs.disp.morph$Group == 'Mainland'),], 
             aes(x = Group, y = Morphological.Disparity, fill = Group, color = Group),
             position = position_jitter(width = 0.1), alpha = 0.15) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), size = 1, width = 0.3) +
  geom_point(aes(fill = Group), pch = 21, size = 5, stroke = 1) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_blank(), 
        legend.position = 'none')

eco.disp.p <-
  ggplot(data = eco.disp, aes(x = Group, y = Ecological.Disparity, fill = Group)) +
  geom_point(data = bs.disp.eco[-which(bs.disp.eco$Group == 'Mainland'),], 
             aes(x = Group, y = Ecological.Disparity, fill = Group, color = Group),
             position = position_jitter(width = 0.1), alpha = 0.15) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), size = 1, width = 0.3) +
  geom_point(aes(fill = Group), pch = 21, size = 5, stroke = 1) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 20) +
  theme(axis.title.x = element_blank(), 
        legend.position = 'none')

disp.p <- 
  plot_grid(broad.morph.disp.p, morph.disp.p,
            broad.eco.disp.p, eco.disp.p,
            ncol = 2, align = 'V')

pdf('Bootstrapped-Contemporary-Disparity.pdf', width = 6, height = 9.75)
disp.p
dev.off()
