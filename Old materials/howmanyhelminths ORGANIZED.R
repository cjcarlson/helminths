library(reshape2)
library(vegan)

setwd('C:/Users/cjcar/Documents/Github/helminths')
load('intMats.RData')
View(data.frame(intMats$Iceland))

base <- data.frame(matrix(nrow=0,ncol=3))

for (i in 1:length(intMats)) {
  df <- data.frame(intMats[[i]])
  sp <- colnames(df)
  
  if(length(sp)>0) {
  stuff <- data.frame(matrix(nrow=length(sp),ncol=3))
  stuff[,1] <- names(intMats)[i]
  stuff[,2] <- sp
  
  for (j in 1:length(sp)) {
    stuff[j,3] <- as.numeric(sum(df[,j])>0)
  }
  
  base <- rbind(base,stuff)
  print(i)
  }
}

basemat <- acast(base, X2~X1, value.var="X3")
basemat <- data.frame(basemat)

library(fossil)
sp <- fossil::spp.est(basemat, abund=FALSE)
sp[ncol(basemat),]
sp <- data.frame(sp)
names(sp) <- c('Samples','S.obs','S.obs.u','S.obs.l','Chao2','Chao2.u','Chao2.l','ICE','ICE.u','ICE.l','Jack1','Jack1.u','Jack1.l')

library(ggplot2)
library(cowplot)

g1 <- ggplot(sp, aes(Samples)) + xlim(c(20,400)) + ylim(c(0,50000))+
  geom_line(aes(y=Jack1), colour="red") +
  geom_ribbon(aes(ymin=Jack1.l, ymax=Jack1.u), alpha=0.2) 

g2 <- ggplot(sp, aes(Samples)) + xlim(c(20,400)) + ylim(c(0,50000))+
  geom_line(aes(y=ICE), colour="green") + 
  geom_ribbon(aes(ymin=ICE.l, ymax=ICE.u), alpha=0.2) 

g3 <- ggplot(sp, aes(Samples)) + xlim(c(20,400)) + ylim(c(0,50000))+
  geom_line(aes(y=Chao2), colour="blue") + 
  geom_ribbon(aes(ymin=Chao2.l, ymax=Chao2.u), alpha=0.2) 

g4 <- ggplot(sp, aes(Samples)) + xlim(c(20,400)) + ylim(c(0,50000))+
  geom_line(aes(y=S.obs), colour="black") + 
  geom_ribbon(aes(ymin=S.obs.l, ymax=S.obs.u), alpha=0.2) 

plot_grid(g4,g1,g2,g3,labels=c('A','B','C','D'))


#

load('/Users/Colin/Dropbox/helminths/NHMwithCites.RData')

basemat2 <- acast(retDF,Parasite~country,value.var='CitationNumber',sum)
basemat2 <- data.frame(basemat2)

sp2 <- fossil::spp.est(basemat2, abund=FALSE, counter=TRUE)
sp2[ncol(basemat2),]


sp2 <- data.frame(sp2)
names(sp2) <- c('Samples','S.obs','S.obs.u','S.obs.l','Chao2','Chao2.u','Chao2.l','ICE','ICE.u','ICE.l','Jack1','Jack1.u','Jack1.l')

library(ggplot2)
library(cowplot)

g5 <- ggplot(sp2, aes(Samples)) + xlim(c(20,250)) + ylim(c(0,50000))+
  geom_line(aes(y=Jack1), colour="red") +
  geom_ribbon(aes(ymin=Jack1.l, ymax=Jack1.u), alpha=0.2) 

g6 <- ggplot(sp2, aes(Samples)) + xlim(c(20,250)) + ylim(c(0,50000))+
  geom_line(aes(y=ICE), colour="green") + 
  geom_ribbon(aes(ymin=ICE.l, ymax=ICE.u), alpha=0.2) 

g7 <- ggplot(sp2, aes(Samples)) + xlim(c(20,250)) + ylim(c(0,50000))+
  geom_line(aes(y=Chao2), colour="blue") + 
  geom_ribbon(aes(ymin=Chao2.l, ymax=Chao2.u), alpha=0.2) 

g8 <- ggplot(sp2, aes(Samples)) + xlim(c(20,250)) + ylim(c(0,50000))+
  geom_line(aes(y=S.obs), colour="black") + 
  geom_ribbon(aes(ymin=S.obs.l, ymax=S.obs.u), alpha=0.2) 

plot_grid(g4,g1,g2,g3,g8,g5,g6,g7,nrow = 2,ncol=4,labels=c('A','B','C','D','E','F','G','H'))


####

basemat2 <- acast(retDF,Parasite~Host,value.var='CitationNumber',sum)
basemat2 <- data.frame(basemat2)

sp2 <- fossil::spp.est(basemat2, abund=FALSE, counter=TRUE)
sp2[ncol(basemat2),]


sp2 <- data.frame(sp2)
names(sp2) <- c('Samples','S.obs','S.obs.u','S.obs.l','Chao2','Chao2.u','Chao2.l','ICE','ICE.u','ICE.l','Jack1','Jack1.u','Jack1.l')

library(ggplot2)
library(cowplot)

g9 <- ggplot(sp2, aes(Samples)) + xlim(c(20,12500)) + ylim(c(0,50000))+
  geom_line(aes(y=Jack1), colour="red") +
  geom_ribbon(aes(ymin=Jack1.l, ymax=Jack1.u), alpha=0.2) 

g10 <- ggplot(sp2, aes(Samples)) + xlim(c(20,12500)) + ylim(c(0,50000))+
  geom_line(aes(y=ICE), colour="green") + 
  geom_ribbon(aes(ymin=ICE.l, ymax=ICE.u), alpha=0.2) 

g11 <- ggplot(sp2, aes(Samples)) + xlim(c(20,12500)) + ylim(c(0,50000))+
  geom_line(aes(y=Chao2), colour="blue") + 
  geom_ribbon(aes(ymin=Chao2.l, ymax=Chao2.u), alpha=0.2) 

g12 <- ggplot(sp2, aes(Samples)) + xlim(c(20,12500)) + ylim(c(0,50000))+
  geom_line(aes(y=S.obs), colour="black") + 
  geom_ribbon(aes(ymin=S.obs.l, ymax=S.obs.u), alpha=0.2) 

plot_grid(g9,g10,g11,g12,labels=c('A','B','C','D'))


### 
country.all <- unique(retDF$country)
country.sub <- c()
for (i in 1:length(country.all)) {
  ret.sub <- retDF[retDF$country==country.all[i],]
  if (nrow(ret.sub)>19) {country.sub <- c(country.sub,country.all[i])}
}
ret.s <- retDF[retDF$country %in% country.sub,]

d <- data.frame(unique(ret.s$country))
d[,2:13] <-0
names(d)<-c('country','S.obs','S.obs.u','S.obs.l','Chao2','Chao2.u','Chao2.l','ICE','ICE.u',
               'ICE.l','Jack1','Jack1.u','Jack1.l')

for (i in c(1:length(country.sub))){
  ret.sub <- retDF[retDF$country==country.sub[i],]
  basemat2 <- acast(ret.sub,Parasite~Host,value.var='CitationNumber',sum)
  basemat2 <- data.frame(basemat2)
  basemat2[is.na(basemat2)] <- 0
  
  sp2 <- fossil::spp.est(basemat2, abund=FALSE, counter=TRUE)
  d[i,2:13]<-sp2[ncol(basemat2),2:13]
  print(i)
  print(country.sub[i])
}


#
#
#
#
# Functional vegan scripts

basemat <- acast(base, X1~X2, value.var="X3",fun=sum)
basemat <- data.frame(basemat)
sp1 <- specaccum(basemat)
summary(sp1)
plot(sp1)

mod1 <- fitspecaccum(sp1,'lomolino')
mod2 <- fitspecaccum(sp1,'gleason')
mod3 <- fitspecaccum(sp1,'asymp')
mod4 <- fitspecaccum(sp1,'gompertz')
mod5 <- fitspecaccum(sp1,'michaelis-menten')
mod6 <- fitspecaccum(sp1,'logis')
mod7 <- fitspecaccum(sp1,'weibull')

# Many of these crash

