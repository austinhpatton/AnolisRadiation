#####
# Merge and Join BAMM pPC analyses

library(BAMMtools)

setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/BAMM-PCs/')


m.ed.txt<- list()
f.ed.txt <- list()

prefix <- '~/Dropbox/Research/Anole_Diversification/FinalAnalyses/BAMM-PCs/'
traits <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')

# Read in the event data text files, remove the burnin, and rescale to common scale
for(i in 1:5){
  # Males
  path <- paste0(prefix, 'PC', i, '/Males/')
  setwd(path)

  m.ed.txt[[i]] <- 
    read.csv(list.files(path, '*event*'))
  m.ed.txt[[i]] <- m.ed.txt[[i]][which(m.ed.txt[[i]]$generation > 20000000),]
  m.ed.txt[[i]]$generation <- m.ed.txt[[i]]$generation - 20000000

  # Females
  path <- paste0(prefix, 'PC', i, '/Females/')
  setwd(path)
  
 f.ed.txt[[i]] <- 
    read.csv(list.files(path, '*event*'))
 f.ed.txt[[i]] <- f.ed.txt[[i]][which(f.ed.txt[[i]]$generation > 20000000),]
 f.ed.txt[[i]]$generation <- f.ed.txt[[i]]$generation - 20000000
}
  
# Now rescale and then combine them iteratively, treating them as consecutive (won't matter as each will be sampled equally in final analysis)
setwd('~/Dropbox/Research/Anole_Diversification/FinalAnalyses/BAMM-PCs/')

# Males
m.ed.txt[[2]]$generation <- m.ed.txt[[2]]$generation + max(m.ed.txt[[1]]$generation)
m.ed.txt[[3]]$generation <- m.ed.txt[[3]]$generation + max(m.ed.txt[[2]]$generation)
m.ed.txt[[4]]$generation <- m.ed.txt[[4]]$generation + max(m.ed.txt[[3]]$generation)
m.ed.txt[[5]]$generation <- m.ed.txt[[5]]$generation + max(m.ed.txt[[4]]$generation)

m.comb.ed <- rbind(m.ed.txt[[1]], m.ed.txt[[2]], m.ed.txt[[3]], m.ed.txt[[4]], m.ed.txt[[5]])
write.csv(m.comb.ed, 'Male-phylPC-Bamm-Combined-EventData.txt', row.names = F, quote = F)

# Females
f.ed.txt[[2]]$generation <- f.ed.txt[[2]]$generation + max(f.ed.txt[[1]]$generation)
f.ed.txt[[3]]$generation <- f.ed.txt[[3]]$generation + max(f.ed.txt[[2]]$generation)
f.ed.txt[[4]]$generation <- f.ed.txt[[4]]$generation + max(f.ed.txt[[3]]$generation)
f.ed.txt[[5]]$generation <- f.ed.txt[[5]]$generation + max(f.ed.txt[[4]]$generation)

f.comb.ed <- rbind(f.ed.txt[[1]], f.ed.txt[[2]], f.ed.txt[[3]], f.ed.txt[[4]], f.ed.txt[[5]])
write.csv(f.comb.ed, 'Female-phylPC-Bamm-Combined-EventData.txt', row.names = F, quote = F)


########################
# Now read in for proper BAMM analysis 

m.ed <- getEventData('./Male-PC.tree', eventdata = './Male-phylPC-Bamm-Combined-EventData.txt', 
                     nsamples = 25000, type = 'trait', verbose = T)

f.ed <- getEventData('./Female-PC.tree', eventdata = './Female-phylPC-Bamm-Combined-EventData.txt', 
                     nsamples = 25000, type = 'trait', verbose = T)


# Now plot combined PCs for 
pdf('Male-Female-Combined-PCs.pdf', width = 10, height = 8)
par(mfrow = c(1, 2))
plot(m.ed, 2, breaksmethod = 'jenks', labels = T, cex = 0.3, show = T, logcolor = F, legend = T, lwd = 2)
plot(f.ed, 2, breaksmethod = 'jenks', labels = T, cex = 0.3, show = T, logcolor = F, legend = T, lwd = 2)
dev.off()

m.ed.cred <- credibleShiftSet(m.ed, 2)
f.ed.cred <- credibleShiftSet(f.ed, 2)

plot.credibleshiftset(m.ed.cred, 2, breaksmethod = 'jenks', labels = T, cex = 0.3, show = T, logcolor = F, legend = T, lwd = 2)
plot.credibleshiftset(f.ed.cred, 2, breaksmethod = 'jenks', labels = T, cex = 0.3, show = T, logcolor = F, legend = T, lwd = 2)

# Now plot each PC separately
# Males
m.pc1 <- getEventData('./Male-PC.tree', eventdata = './PC1/Males/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)
m.pc2 <- getEventData('./Male-PC.tree', eventdata = './PC2/Males/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)
m.pc3 <- getEventData('./Male-PC.tree', eventdata = './PC3/Males/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)
m.pc4 <- getEventData('./Male-PC.tree', eventdata = './PC4/Males/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)
m.pc5 <- getEventData('./Male-PC.tree', eventdata = './PC5/Males/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)


col.interval <- c(6e-4, 0.012)
pdf('Male-PC-BAMM-Rates-SingleScale.pdf', width = 16, height = 8)
par(mfrow = c(1, 6))
plot(m.ed, 2, breaksmethod = 'jenks', 
     labels = T, cex = 0.5, show = T, logcolor = T, 
     legend = T, lwd = 2, color.interval = col.interval, JenksSubset = 10)
plot(m.pc1, 2, breaksmethod = 'linear', 
     labels = T, cex = 0.5, show = T, logcolor = F, 
     legend = T, lwd = 2, color.interval = col.interval)
plot(m.pc2, 2, breaksmethod = 'linear', labels = T, 
     cex = 0.5, show = T, logcolor = F, 
     legend = T, lwd = 2, color.interval = col.interval)
plot(m.pc3, 2, breaksmethod = 'linear', labels = T, 
     cex = 0.5, show = T, logcolor = F, 
     legend = T, lwd = 2, color.interval = col.interval)
plot(m.pc4, 2, breaksmethod = 'linear', labels = T, 
     cex = 0.5, show = T, logcolor = F, 
     legend = T, lwd = 2, color.interval = col.interval)
plot(m.pc5, 2, breaksmethod = 'quantile', labels = T, 
     cex = 0.5, show = T, logcolor = F, 
     legend = T, lwd = 2, color.interval = col.interval)
dev.off()

m.pc1.best <- getBestShiftConfiguration(m.pc1, 2)
m.pc2.best <- getBestShiftConfiguration(m.pc2, 2)
m.pc3.best <- getBestShiftConfiguration(m.pc3, 2)
m.pc4.best <- getBestShiftConfiguration(m.pc4, 2)
m.pc5.best <- getBestShiftConfiguration(m.pc5, 2)

plot(m.pc1.best, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
plot(m.pc2.best, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
plot(m.pc3.best, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
plot(m.pc4.best, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
plot(m.pc5.best, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)


# Females
f.pc1 <- getEventData('./Female-PC.tree', eventdata = './PC1/Females/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)
f.pc2 <- getEventData('./Female-PC.tree', eventdata = './PC2/Females/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)
f.pc3 <- getEventData('./Female-PC.tree', eventdata = './PC3/Females/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)
f.pc4 <- getEventData('./Female-PC.tree', eventdata = './PC4/Females/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)
f.pc5 <- getEventData('./Female-PC.tree', eventdata = './PC5/Females/Anole_BAMM_P2_6Ch_50m_event_data.txt', 
                      nsamples = 5000, type = 'trait', verbose = T)

pdf('Female-PCs1-5-BAMM.pdf', width = 16, height = 8)
par(mfrow = c(1, 5))
plot(f.pc1, 2, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
plot(f.pc2, 2, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
plot(f.pc3, 2, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
plot(f.pc4, 2, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
plot(f.pc5, 2, breaksmethod = 'jenks', labels = T, cex = 0.5, show = T, logcolor = F, legend = T, lwd = 2)
dev.off()


