
library(codependent)
library(tidyverse)
library(reshape2)

set.seed(8675309)

setwd('~/Github/helminths')
all.data <- read_delim('./data/associations cleaned v2.csv',delim=',', col_types='cccccc')[,-1]

all.data$hostgroup[all.data$hostgroup=='Chondrostei'] <- 'Osteichthyes'
all.data$hostgroup[all.data$hostgroup=='Holostei'] <- 'Osteichthyes'
all.data$hostgroup[all.data$hostgroup=='Cladistei'] <- 'Osteichthyes'
all.data$hostgroup[all.data$hostgroup=='Teleostei'] <- 'Osteichthyes'

for(i in unique(all.data$group)) {
  print(i)
  k <- all.data[all.data$group==i,]
  for(j in c('Chondrichthyes','Osteichthyes','Amphibia','Reptilia','Aves','Mammalia')) {
    l <- k[k$hostgroup==j,]
    l <- na.omit(l[,1:2])
    print(round(mean(table(unique(l)[,'Parasite'])),2))
  }
}

unique(all.data$group)
