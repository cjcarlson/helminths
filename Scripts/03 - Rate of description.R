library(segmented)

npc <- read.csv('~/Github/helminths/Data/timeseries_NPC.csv',
                stringsAsFactor=FALSE)[,-1]
nhm <- read.csv('~/Github/helminths/Data/timeseries_NHM_v2.csv',
                stringsAsFactor=FALSE)[,-1]

npc.df <- data.frame(table(npc$yearone))
names(npc.df) <- c('year','freq')
npc.df$cumsum <- cumsum(npc.df$freq)
npc.df$year <- as.numeric(as.character(npc.df$year))

npc_dates <- seq(min(npc.df$year), max(npc.df$year), 
                  by = 1)
npc_dates <- data.frame(year = npc_dates)
npc_totals <- merge(npc_dates, npc.df, by = "year", 
                     all.x = TRUE)

for (i in 1:nrow(npc_totals)) {
  if(is.na(npc_totals$freq[i])) {
    npc_totals$freq[i] <- 0
  } 
  if(is.na(npc_totals$cumsum[i])) {
    npc_totals$cumsum[i] <- npc_totals$cumsum[i-1] 
  } 
}



nhm.df <- data.frame(table(nhm$year))
names(nhm.df) <- c('year','freq')
nhm.df$cumsum <- cumsum(nhm.df$freq)
nhm.df$year <- as.numeric(as.character(nhm.df$year))

nhm_dates <- seq(min(nhm.df$year), max(nhm.df$year), 
                 by = 1)
nhm_dates <- data.frame(year = nhm_dates)
nhm_totals <- merge(nhm_dates, nhm.df, by = "year", 
                    all.x = TRUE)

for (i in 1:nrow(nhm_totals)) {
  if(is.na(nhm_totals$freq[i])) {
    nhm_totals$freq[i] <- 0
  } 
  if(is.na(nhm_totals$cumsum[i])) {
    nhm_totals$cumsum[i] <- nhm_totals$cumsum[i-1] 
  } 
}



totals <- merge(nhm_totals, npc_totals, by='year', all.x=TRUE)
names(totals) <- c('year','freq.nhm','cumsum.nhm','freq.npc','cumsum.npc')




par(mfrow=c(2,1))
par(mar=c(0,5,4,5))

with(totals,
     barplot(totals$freq.nhm, ylab='Number of species described', axes=FALSE,
             xlab='', ylim=c(0,375)), yaxt='n')
axis(side = 2, cex.axis=1)
par(new = T)
with(totals,
     plot(totals$cumsum.nhm ~ totals$year, ylab='',
          type='l', axes=T,xaxt='n',yaxt='n',xlab='',
          col='blue', lwd=3))
axis(side = 4, cex.axis=1)
#mtext(side = 1, line = 2, 'Year')
mtext(side = 4, line = 3, 'Cumulative total species described')


lfit <- lm(cumsum.nhm ~ 0 + year, data = totals)
sfit <- segmented(lfit, seg.Z = ~ year)
summary(sfit); pscore.test(sfit, ~year)
lines(totals$year,fitted(sfit),col='red',lwd=3.5,lty=2)



par(mar=c(4,5,0,5))
with(totals,
     barplot(totals$freq.npc, ylab='Number of species collected', axes=FALSE,
             xlab='', ylim=c(0,375)), yaxt='n')
axis(side = 2, cex.axis=1)
par(new = T)
with(totals,
     plot(totals$cumsum.npc ~ totals$year, ylab='',
          type='l', axes=T,xaxt='n',yaxt='n',xlab='',
          col='blue', lwd=3))
axis(side = 1, cex.axis=1)
axis(side = 4, cex.axis=1)
mtext(side = 1, line = 2.5, 'Year')
mtext(side = 4, line = 3, 'Cumulative total species collected')

lfit <- lm(cumsum.npc ~ 0 + year, data = totals)
sfit <- segmented(lfit, seg.Z = ~ year)
summary(sfit); pscore.test(sfit, ~year)
lines(totals$year[!is.na(totals$freq.npc)],fitted(sfit),col='red',lwd=3.5,lty=2)



