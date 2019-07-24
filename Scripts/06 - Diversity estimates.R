
# RUN 03 - RATE OF DESCRIPTION.R FIRST

library(codependent)
library(tidyverse)
library(reshape2)

set.seed(8675309)

setwd('~/Github/helminths')
all.data <- read_delim('./data/associations cleaned v2.csv',delim=',')[,-1]

all.data$hostgroup[all.data$hostgroup=='Chondrostei'] <- 'Osteichthyes'
all.data$hostgroup[all.data$hostgroup=='Holostei'] <- 'Osteichthyes'
all.data$hostgroup[all.data$hostgroup=='Cladistei'] <- 'Osteichthyes'
all.data$hostgroup[all.data$hostgroup=='Teleostei'] <- 'Osteichthyes'





#################
acan.df <- all.data[all.data$group=='Acanthocephalans',c('Parasite','hostgroup')]
acan.df$value =1 
acan.df <- unique(na.omit(acan.df))
acan.df <- acast(acan.df , Parasite ~ hostgroup)
acan.df[is.na(acan.df)] <- 0
acan.df <- acan.df[,c('Amphibia', 'Aves','Mammalia','Reptilia','Osteichthyes','Chondrichthyes')] # ADD FISH BACK - BUT > 2 FISH GROUPS CURRENTLY
acan.df <- acan.df[rowSums(acan.df)>0,]

all.data %>% filter(group=='Acanthocephalans') %>% filter(hostgroup=='Amphibia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> acan.amph
all.data %>% filter(group=='Acanthocephalans') %>% filter(hostgroup=='Aves') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> acan.aves
all.data %>% filter(group=='Acanthocephalans') %>% filter(hostgroup=='Mammalia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> acan.mamm
all.data %>% filter(group=='Acanthocephalans') %>% filter(hostgroup=='Reptilia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> acan.rept
all.data %>% filter(group=='Acanthocephalans') %>% filter(hostgroup=='Osteichthyes') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> acan.oste
all.data %>% filter(group=='Acanthocephalans') %>% filter(hostgroup=='Chondrichthyes') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> acan.chon

est.a.amp <- copredict(acan.amph, n=7302, iter=100)[[1]][1]
est.a.ave <- copredict(acan.aves, n=10425, iter=100)[[1]][1]
est.a.mam <- copredict(acan.mamm, n=5513, iter=100)[[1]][1]
est.a.rep <- copredict(acan.rept, n=10038, iter=100)[[1]][1]
est.a.ost <- copredict(acan.oste, n=28000, iter=100)[[1]][1]
est.a.con <- copredict(acan.chon, n=1111, iter=100)[[1]][1]
est.a.fake <- c(est.a.amp,est.a.ave,est.a.mam,est.a.rep,est.a.ost,est.a.con)

acan.df <- cbind(data.frame(parasite=c(1:nrow(acan.df))), acan.df)
sum(est.a.fake)
acan.corr <- multigroup(acan.df, orders=c("A","B","C","D","E","F"), 
     est.a.fake)


#################
cest.df <- all.data[all.data$group=='Cestodes',c('Parasite','hostgroup')]
cest.df$value =1 
cest.df <- unique(na.omit(cest.df))
cest.df <- acast(cest.df , Parasite ~ hostgroup)
cest.df[is.na(cest.df)] <- 0
cest.df <- cest.df[,c('Amphibia', 'Aves','Mammalia','Reptilia','Osteichthyes','Chondrichthyes')] # ADD FISH BACK - BUT > 2 FISH GROUPS CURRENTLY
cest.df <- cest.df[rowSums(cest.df)>0,]

all.data %>% filter(group=='Cestodes') %>% filter(hostgroup=='Amphibia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> cest.amph
all.data %>% filter(group=='Cestodes') %>% filter(hostgroup=='Aves') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> cest.aves
all.data %>% filter(group=='Cestodes') %>% filter(hostgroup=='Mammalia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> cest.mamm
all.data %>% filter(group=='Cestodes') %>% filter(hostgroup=='Reptilia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> cest.rept
all.data %>% filter(group=='Cestodes') %>% filter(hostgroup=='Osteichthyes') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> cest.oste
all.data %>% filter(group=='Cestodes') %>% filter(hostgroup=='Chondrichthyes') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> cest.chon

est.c.amp <- copredict(cest.amph, n=7302, iter=100)[[1]][1]
est.c.ave <- copredict(cest.aves, n=10425, iter=100)[[1]][1]
est.c.mam <- copredict(cest.mamm, n=5513, iter=100)[[1]][1]
est.c.rep <- copredict(cest.rept, n=10038, iter=100)[[1]][1]
est.c.ost <- copredict(cest.oste, n=28000, iter=100)[[1]][1]
est.c.con <- copredict(cest.chon, n=1111, iter=100)[[1]][1]
est.c.fake <- c(est.c.amp,est.c.ave,est.c.mam,est.c.rep,est.c.ost,est.c.con)

cest.df <- cbind(data.frame(parasite=c(1:nrow(cest.df))), cest.df)
sum(est.c.fake)
cest.corr <- multigroup(cest.df, orders=c("A","B","C","D","E","F"), 
                   est.c.fake)


#################
nema.df <- all.data[all.data$group=='Nematodes',c('Parasite','hostgroup')]
nema.df$value =1 
nema.df <- unique(na.omit(nema.df))
nema.df <- acast(nema.df , Parasite ~ hostgroup)
nema.df[is.na(nema.df)] <- 0
nema.df <- nema.df[,c('Amphibia', 'Aves','Mammalia','Reptilia','Osteichthyes','Chondrichthyes')] # ADD FISH BACK - BUT > 2 FISH GROUPS CURRENTLY
nema.df <- nema.df[rowSums(nema.df)>0,]

all.data %>% filter(group=='Nematodes') %>% filter(hostgroup=='Amphibia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> nema.amph
all.data %>% filter(group=='Nematodes') %>% filter(hostgroup=='Aves') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> nema.aves
all.data %>% filter(group=='Nematodes') %>% filter(hostgroup=='Mammalia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> nema.mamm
all.data %>% filter(group=='Nematodes') %>% filter(hostgroup=='Reptilia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> nema.rept
all.data %>% filter(group=='Nematodes') %>% filter(hostgroup=='Osteichthyes') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> nema.oste
all.data %>% filter(group=='Nematodes') %>% filter(hostgroup=='Chondrichthyes') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> nema.chon

est.n.amp <- copredict(nema.amph, n=7302, iter=100)[[1]][1]
est.n.ave <- copredict(nema.aves, n=10425, iter=100)[[1]][1]
est.n.mam <- copredict(nema.mamm, n=5513, iter=100)[[1]][1]
est.n.rep <- copredict(nema.rept, n=10038, iter=100)[[1]][1]
est.n.ost <- copredict(nema.oste, n=28000, iter=100)[[1]][1]
est.n.con <- copredict(nema.chon, n=1111, iter=100)[[1]][1]
est.n.fake <- c(est.n.amp,est.n.ave,est.n.mam,est.n.rep,est.n.ost,est.n.con)

nema.df <- cbind(data.frame(parasite=c(1:nrow(nema.df))), nema.df)
sum(est.n.fake)
nema.corr <- multigroup(nema.df, orders=c("A","B","C","D","E","F"), 
                   est.n.fake)



#################
trem.df <- all.data[all.data$group=='Trematodes',c('Parasite','hostgroup')]
trem.df$value =1 
trem.df <- unique(na.omit(trem.df))
trem.df <- acast(trem.df , Parasite ~ hostgroup)
trem.df[is.na(trem.df)] <- 0
trem.df <- trem.df[,c('Amphibia', 'Aves','Mammalia','Reptilia','Osteichthyes','Chondrichthyes')] # ADD FISH BACK - BUT > 2 FISH GROUPS CURRENTLY
trem.df <- trem.df[rowSums(trem.df)>0,]

all.data %>% filter(group=='Trematodes') %>% filter(hostgroup=='Amphibia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> trem.amph
all.data %>% filter(group=='Trematodes') %>% filter(hostgroup=='Aves') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> trem.aves
all.data %>% filter(group=='Trematodes') %>% filter(hostgroup=='Mammalia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> trem.mamm
all.data %>% filter(group=='Trematodes') %>% filter(hostgroup=='Reptilia') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> trem.rept
all.data %>% filter(group=='Trematodes') %>% filter(hostgroup=='Osteichthyes') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> trem.oste
all.data %>% filter(group=='Trematodes') %>% filter(hostgroup=='Chondrichthyes') %>% 
  select(Host, Parasite) %>% data.frame() %>% unique() -> trem.chon

est.t.amp <- copredict(trem.amph, n=7302, iter=100)[[1]][1]
est.t.ave <- copredict(trem.aves, n=10425, iter=100)[[1]][1]
est.t.mam <- copredict(trem.mamm, n=5513, iter=100)[[1]][1]
est.t.rep <- copredict(trem.rept, n=10038, iter=100)[[1]][1]
est.t.ost <- copredict(trem.oste, n=28000, iter=100)[[1]][1]
est.t.con <- copredict(trem.chon, n=1111, iter=100)[[1]][1]
est.t.fake <- c(est.t.amp,est.t.ave,est.t.mam,est.t.rep,est.t.ost,est.t.con)

trem.df <- cbind(data.frame(parasite=c(1:nrow(trem.df))), trem.df)
sum(est.t.fake)
trem.corr <- multigroup(trem.df, orders=c("A","B","C","D","E","F"), 
                   est.t.fake)



################################################




# INTERMISSION




################################################


true.total <- acan.corr+cest.corr+nema.corr+trem.corr

sum(est.a.fake,
    est.c.fake,
    est.n.fake,
    est.t.fake)



#################### 

round(100*length(unique(rownames(acan.df[acan.df$Chondrichthyes==1,])))/est.a.con)
round(100*length(unique(rownames(acan.df[acan.df$Osteichthyes==1,])))/est.a.ost)
round(100*length(unique(rownames(acan.df[acan.df$Amphibia==1,])))/est.a.amp)
round(100*length(unique(rownames(acan.df[acan.df$Reptilia==1,])))/est.a.rep)
round(100*length(unique(rownames(acan.df[acan.df$Aves==1,])))/est.a.ave)
round(100*length(unique(rownames(acan.df[acan.df$Mammalia==1,])))/est.a.mam)


round(100*length(unique(rownames(cest.df[cest.df$Chondrichthyes==1,])))/est.c.con)
round(100*length(unique(rownames(cest.df[cest.df$Osteichthyes==1,])))/est.c.ost)
round(100*length(unique(rownames(cest.df[cest.df$Amphibia==1,])))/est.c.amp)
round(100*length(unique(rownames(cest.df[cest.df$Reptilia==1,])))/est.c.rep)
round(100*length(unique(rownames(cest.df[cest.df$Aves==1,])))/est.c.ave)
round(100*length(unique(rownames(cest.df[cest.df$Mammalia==1,])))/est.c.mam)


round(100*length(unique(rownames(nema.df[nema.df$Chondrichthyes==1,])))/est.n.con)
round(100*length(unique(rownames(nema.df[nema.df$Osteichthyes==1,])))/est.n.ost)
round(100*length(unique(rownames(nema.df[nema.df$Amphibia==1,])))/est.n.amp)
round(100*length(unique(rownames(nema.df[nema.df$Reptilia==1,])))/est.n.rep)
round(100*length(unique(rownames(nema.df[nema.df$Aves==1,])))/est.n.ave)
round(100*length(unique(rownames(nema.df[nema.df$Mammalia==1,])))/est.n.mam)


round(100*length(unique(rownames(trem.df[trem.df$Chondrichthyes==1,])))/est.t.con)
round(100*length(unique(rownames(trem.df[trem.df$Osteichthyes==1,])))/est.t.ost)
round(100*length(unique(rownames(trem.df[trem.df$Amphibia==1,])))/est.t.amp)
round(100*length(unique(rownames(trem.df[trem.df$Reptilia==1,])))/est.t.rep)
round(100*length(unique(rownames(trem.df[trem.df$Aves==1,])))/est.t.ave)
round(100*length(unique(rownames(trem.df[trem.df$Mammalia==1,])))/est.t.mam)


round(100*(length(unique(rownames(acan.df[acan.df$Chondrichthyes==1,])))+
           length(unique(rownames(cest.df[cest.df$Chondrichthyes==1,])))+
           length(unique(rownames(nema.df[nema.df$Chondrichthyes==1,])))+
           length(unique(rownames(trem.df[trem.df$Chondrichthyes==1,]))))/
        (est.a.con + est.c.con + est.n.con + est.t.con))

round(100*(length(unique(rownames(acan.df[acan.df$Osteichthyes==1,])))+
             length(unique(rownames(cest.df[cest.df$Osteichthyes==1,])))+
             length(unique(rownames(nema.df[nema.df$Osteichthyes==1,])))+
             length(unique(rownames(trem.df[trem.df$Osteichthyes==1,]))))/
        (est.a.ost + est.c.ost + est.n.ost + est.t.ost))

round(100*(length(unique(rownames(acan.df[acan.df$Amphibia==1,])))+
             length(unique(rownames(cest.df[cest.df$Amphibia==1,])))+
             length(unique(rownames(nema.df[nema.df$Amphibia==1,])))+
             length(unique(rownames(trem.df[trem.df$Amphibia==1,]))))/
        (est.a.amp + est.c.amp + est.n.amp + est.t.amp))

round(100*(length(unique(rownames(acan.df[acan.df$Reptilia==1,])))+
             length(unique(rownames(cest.df[cest.df$Reptilia==1,])))+
             length(unique(rownames(nema.df[nema.df$Reptilia==1,])))+
             length(unique(rownames(trem.df[trem.df$Reptilia==1,]))))/
        (est.a.rep + est.c.rep + est.n.rep + est.t.rep))

round(100*(length(unique(rownames(acan.df[acan.df$Aves==1,])))+
             length(unique(rownames(cest.df[cest.df$Aves==1,])))+
             length(unique(rownames(nema.df[nema.df$Aves==1,])))+
             length(unique(rownames(trem.df[trem.df$Aves==1,]))))/
        (est.a.ave + est.c.ave + est.n.ave + est.t.ave))

round(100*(length(unique(rownames(acan.df[acan.df$Mammalia==1,])))+
             length(unique(rownames(cest.df[cest.df$Mammalia==1,])))+
             length(unique(rownames(nema.df[nema.df$Mammalia==1,])))+
             length(unique(rownames(trem.df[trem.df$Mammalia==1,]))))/
        (est.a.mam + est.c.mam + est.n.mam + est.t.mam))


# row percentages

round(100*length(unique(rownames(acan.df)))/acan.corr)
round(100*length(unique(rownames(cest.df)))/cest.corr)
round(100*length(unique(rownames(nema.df)))/nema.corr)
round(100*length(unique(rownames(trem.df)))/trem.corr)

round(100*(length(unique(rownames(acan.df)))+
           length(unique(rownames(cest.df)))+
           length(unique(rownames(nema.df)))+
           length(unique(rownames(trem.df))))/(acan.corr+cest.corr+nema.corr+trem.corr))

known.nhm <- length(unique(rownames(acan.df)))+
  length(unique(rownames(cest.df)))+
  length(unique(rownames(nema.df)))+
  length(unique(rownames(trem.df)))

#####

## ALL THESE NUMBERS NEED TO BE REDONE AND MAYBE ANNOTATED 

13426/true.total
known.nhm/true.total

(true.total-13426) * (1/npc.rate)
(true.total-known.nhm) * (1/nhm.rate)

(1 + 2.6)*acan.corr
(1 + 2.4)*cest.corr
(1 + 1.2)*nema.corr
(1 + 3.1)*trem.corr

(1 + 2.6)*acan.corr + (1 + 2.4)*cest.corr + (1 + 1.2)*nema.corr + (1 + 3.1)*trem.corr -> total.cryptic

(total.cryptic-13426) * (1/npc.rate)
(total.cryptic-known.nhm) * (1/nhm.rate)
