
hostest <- function(codobj, rate) {(rate*codobj[[1]][[1]]/codobj[[2]][[1]])^(1/codobj[[2]][[2]])}

hostest(est.a.amp, 0.5)

##########


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

est.c.amp <- copredict(cest.amph, n=7302, iter=100)
est.c.ave <- copredict(cest.aves, n=10425, iter=100)
est.c.mam <- copredict(cest.mamm, n=5513, iter=100)
est.c.rep <- copredict(cest.rept, n=10038, iter=100)
est.c.ost <- copredict(cest.oste, n=28000, iter=100)
est.c.con <- copredict(cest.chon, n=1111, iter=100)


hostest(est.c.amp, 0.5) - length(unique(cest.amph$Host))
hostest(est.c.ave, 0.5) - length(unique(cest.aves$Host))
hostest(est.c.mam, 0.5) - length(unique(cest.mamm$Host))
hostest(est.c.rep, 0.5) - length(unique(cest.rept$Host))
hostest(est.c.ost, 0.5) - length(unique(cest.oste$Host))
hostest(est.c.con, 0.5) - length(unique(cest.chon$Host))

