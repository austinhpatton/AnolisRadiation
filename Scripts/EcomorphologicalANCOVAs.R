########### Phylogenetic Ancova #########################################
# Phylogenetic Ancova - ecomorphology

library(geiger)
library(caper)
library(gridExtra)
library(ggpubr)


setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/')

m.dat <- read.csv('./FinalData/Male-ALL-DATA-Trait-Ecol-Regions-FINAL.csv')

# Reduce down to the only groups being used here (greater antilles radiation - here I2, and the two mainland groups)
m.dat <- m.dat[which(m.dat$Group != 'I1'),]

# And remove species that secondarily colonized islands but did not diversify.
m.dat <- 
  m.dat[-which(m.dat$Species %in% c('agassizi', 'bicaorum', 'lineatus', 'roatanensis', 
                                    'concolor', 'pinchoti', 'utilensis', 'townsendi')),]
m.dat$Group <- droplevels(m.dat$Group)


# Set colors for plotting
cols <- c("#faae61", "#d72027", "#add8a4", "#2b83bb")

# Hind Lamellae vs Perch Height
dat <- na.omit(m.dat[,c(1,13,15,23)])

tree <- read.tree('./FinalData/Anole_TimeCal_NoOut.tree')
tree <- drop.tip(tree, which(!tree$tip.label %in% dat$Species))

comDat <- comparative.data(tree, dat, names.col = Species)
inter <- pgls(HindLam ~ Group * PH, data = comDat, 
              lambda = 'ML', delta = 'ML', kappa = 'ML')
inter.anova <- anova(inter)

main <- pgls(HindLam ~ Group + PH, data = comDat, 
             lambda = 'ML', delta = 'ML', kappa = 'ML')
main.anova <- anova(main)

comDat$data$predlm <- main$fitted

p.labs <- c(paste0('Perch Height: P = ', signif(main.anova$`Pr(>F)`[2], 3)),
            paste0('Group: P = ', signif(main.anova$`Pr(>F)`[1], 3)))
lab.y.pos <- max(na.omit(comDat$data$HindLam))

morph.ecol.legend <- ggplot(data = comDat$data, aes(color = Group)) + 
  geom_point(aes(x = HindLam, y = PH)) +
  scale_color_manual(values=cols[2:4])

morph.ecol.legend <- get_legend(morph.ecol.legend)
morph.ecol.legend <- as_ggplot(morph.ecol.legend)

HLam.PHeight <-
  ggplot(data = comDat$data, aes(fill = Group)) + 
  stat_smooth(aes(x = PH, y = predlm, color = Group), size = 2, method = lm, se = F) +
  geom_point(aes(x = PH, y = HindLam), size = 3, pch = 21, stroke = 0.75) +
  theme_bw() +
  theme(legend.position = "none") +
  annotate('text', x=max(comDat$data$PH), 
           y=c(lab.y.pos, lab.y.pos*0.80), 
           hjust=1, label = p.labs) +
  scale_fill_manual(values=cols[2:4]) +
  scale_color_manual(values=cols[2:4]) + 
  ylab('Hindfoot Lamella Count')

################################
# Front Lamellae vs Perch Height
# ~Significant interaction - different lines, different slopes
dat <- na.omit(m.dat[,c(1,12,15,23)])
tree <- read.tree('./FinalData/Anole_TimeCal_NoOut.tree')
tree <- drop.tip(tree, which(!tree$tip.label %in% dat$Species))

comDat <- comparative.data(tree, dat, names.col = Species)
inter <- pgls(ForeLam ~ Group * PH, data = comDat, 
              lambda = 'ML', delta = 'ML', kappa = 'ML')
inter.fit <- anova(inter)
comDat$data$fitted <- inter$fitted

main <- pgls(ForeLam ~ Group + PH, data = comDat, 
             lambda = 'ML', delta = 'ML', kappa = 'ML')
main.fit <- anova(main)
comDat$data$predlm <- main$fitted



p.labs <- c(paste0('Interaction: P = ', signif(inter.fit$`Pr(>F)`[3], 3)),
            paste0('Group: P = ', signif(inter.fit$`Pr(>F)`[1], 3)),
            paste0('Perch Height: P = ', signif(inter.fit$`Pr(>F)`[2], 3)))
lab.x.pos <- max(na.omit(comDat$data$PH))
lab.y.pos <- min(na.omit(comDat$data$ForeLam))

FLam.PHeight <-
  ggplot(data = comDat$data, group = Group) + 
  stat_smooth(aes(x = PH, y = fitted, color = Group), size = 2, method = lm, se = F) +
  geom_point(aes(x = PH, y = ForeLam, fill = Group), pch = 21, size = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  annotate('text', x=lab.x.pos, 
           y=c(lab.y.pos, lab.y.pos*0.90, lab.y.pos*0.80), 
           hjust=1, label = p.labs) +
  scale_fill_manual(values=cols[2:4]) + 
  scale_color_manual(values=cols[2:4]) + 
  ylab('Forefoot Lamella Count')


###################################
# Hindlimb length vs Perch Diameter
dat <- na.omit(m.dat[,c(1,8,18,23)])
tree <- read.tree('./FinalData/Anole_TimeCal_NoOut.tree')
tree <- drop.tip(tree, which(!tree$tip.label %in% dat$Species))

comDat <- comparative.data(tree, dat, names.col = Species)
inter <- pgls(Hind.limb.length ~ Group * PD, data = comDat, 
              lambda = 'ML', delta = 'ML', kappa = 'ML')
inter.fit <- anova(inter)

main <- pgls(Hind.limb.length ~ Group + PD, data = comDat, 
             lambda = 'ML', delta = 'ML', kappa = 'ML')
main.fit <- anova(main)
comDat$data$predlm <- main$fitted


p.labs <- c(paste0('Group: P = ', signif(main.fit$`Pr(>F)`[1], 3)),
            paste0('Perch Diameter: P = ', signif(main.fit$`Pr(>F)`[2], 3)))

lab.y.pos <- min(na.omit(comDat$data$Hind.limb.length))
lab.x.pos <- max(na.omit(comDat$data$PD))

Hlimb.PDiameter <-
  ggplot(data = comDat$data, aes(fill = Group)) +
  stat_smooth(aes(x = PD, y = predlm, color = Group), size = 2, method = lm, se = F) +
  geom_point(aes(x = PD, y = Hind.limb.length), size = 3, pch = 21) +
  theme_bw() +
  theme(legend.position = "none") +
  annotate('text', x=max(comDat$data$PD), 
           y=c(lab.y.pos, lab.y.pos*0.90), 
           hjust=1, label = p.labs) +
  scale_fill_manual(values=cols[2:4]) +
  scale_color_manual(values=cols[2:4]) + 
  ylab('Hindlimb Length')


####################################
# Forelimb length vs Perch Diameter
dat <- na.omit(m.dat[,c(1,10,18,23)])
tree <- read.tree('./FinalData/Anole_TimeCal_NoOut.tree')
tree <- drop.tip(tree, which(!tree$tip.label %in% dat$Species))

comDat <- comparative.data(tree, dat, names.col = Species)
inter <- pgls(Forelimb.length ~ Group * PD, data = comDat, 
              lambda = 'ML', delta = 'ML', kappa = 'ML')
inter.fit <- anova(inter)

p.labs <- c(paste0('Interaction: P = ', signif(inter.fit$`Pr(>F)`[3], 3)),
            paste0('Group: P = ', signif(inter.fit$`Pr(>F)`[1], 3)),
            paste0('Perch Diameter: P = ', signif(inter.fit$`Pr(>F)`[2], 3)))

lab.y.pos <- min(na.omit(comDat$data$Forelimb.length))
lab.x.pos <- max(na.omit(comDat$data$PD))

Flimb.PDiameter <-
  ggplot(data = comDat$data, aes(fill = Group)) +
  stat_smooth(aes(x = PD, y = Forelimb.length, color = Group), size = 2, method = lm, se = F) +
  geom_point(aes(x = PD, y = Forelimb.length), size = 3, pch = 21) +
  theme_bw() +
  theme(legend.position = "none") +
  annotate('text', x=max(comDat$data$PD), 
           y=c(lab.y.pos, lab.y.pos*0.90, lab.y.pos*0.80), 
           hjust=1, label = p.labs) +
  scale_fill_manual(values=cols[2:4]) +
  scale_color_manual(values=cols[2:4]) + 
  ylab('Forelimb Length')

pdf('FourPrimary-Ecomorphological-Relationships.pdf', width = 8, height = 8)
(FLam.PHeight | HLam.PHeight) /
  (Flimb.PDiameter | Hlimb.PDiameter)
dev.off()
