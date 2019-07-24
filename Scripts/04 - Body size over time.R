library(mgcv)
library(tidyverse)
library(pspline)

associations <- read.csv('~/Github/helminths/Data/associations cleaned v2.csv',
                         stringsAsFactors = FALSE)
associations <- associations[,-1]

associations$year <- ''

for(i in 1:nrow(associations)) {
  x <- strsplit(as.character(associations$ParasiteFull[i]),' ')[[1]]
  associations$year[i] <- gsub(")","",x[length(x)])
  associations$year[i] <- gsub("\\(","",associations$year[i])
  print(i)
}

life <- read.csv('~/Github/helminths/Data/CLC_database_lifehistory.csv')
life <- life[life$Stage=='adult',]

life <- life[life$Parasite.species %in% associations$Parasite, ]
life %>% group_by(Parasite.species) %>% summarize(length = mean(Length)) -> life

life$year <- NA
for(i in 1:nrow(life)){
  life$year[i] <- as.numeric(associations[associations$Parasite==as.character(life$Parasite.species[i]),'year'])[1]
}
life <- life[complete.cases(life),]

MyGAM1 <- gam(log(length) ~ s(year), data=life)
response1 <- predict(MyGAM1, type="response", se.fit=T)

################################################

pantheria <- read.csv('~/Github/helminths/Data/pantheria.csv')
pantheria <- pantheria[,c(6,8)]
names(pantheria) <- c('Host','Bodymass')

associations$bodymass <- NA
pantheria$Host <- as.character(pantheria$Host)

for(i in 1:nrow(associations)) {
  if (associations$Host[i] %in% unique(pantheria$Host)) {
    associations$bodymass[i] <- pantheria[pantheria$Host %in% as.character(associations$Host[i]),]$Bodymass[1]
  }
  print(i)
}

assoc.df <- merge(associations,pantheria,by='Host',all.x=FALSE)

# It's the average host mass because we have no way to know which hosts sampled first

assoc2 <- na.omit(assoc.df)
assoc2 %>% group_by(Parasite,year) %>% summarize(mean.mass = mean(bodymass)) -> assoc2
assoc2 <- assoc2[,c('year','mean.mass')]
assoc2$year <- as.numeric(assoc2$year)

MyGAM2 <- bam(log(mean.mass) ~ s(year), data=assoc2)
response2 <- predict(MyGAM2, type="response", se.fit=T)

par(mfrow=c(1,2))
par(mar=c(4,4,2,2))

plot(0, type="n", bty='n', main="", xlab="Year", ylab="log(Worm length)", lwd=3, ylim=c(0,8), xlim=c(1750,2000))
points(data.frame(life[,3],log(life[,2])), pch=16, col='light blue', cex=0.9)
lines(sm.spline(MyGAM1$model$year , response1$fit) , lwd = 3)
lines(sm.spline(MyGAM1$model$year , response1$fit+1.96*response1$se) , lty = 3 , lwd = 2)
lines(sm.spline(MyGAM1$model$year , response1$fit-1.96*response1$se) , lty = 3 , lwd = 2)

plot(0, type="n", bty='n', main="", xlab="Year", ylab="log(Mammal host mass)", lwd=3, ylim=c(0,20), xlim=c(1750,2000))
points(log(mean.mass) ~ year, data=assoc2, pch=16, col='light blue', cex=0.9)
lines(sm.spline(MyGAM2$model$year , response2$fit) , lwd = 3)
lines(sm.spline(MyGAM2$model$year , response2$fit+1.96*response2$se) , lty = 3 , lwd = 2)
lines(sm.spline(MyGAM2$model$year , response2$fit-1.96*response2$se) , lty = 3 , lwd = 2)
