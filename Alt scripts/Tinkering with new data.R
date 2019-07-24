load("C:/Users/cjcar/Downloads/nhmByParasite.RData")
helmDF <- unique(helmDF[,1:2])

library(plyr)

associations <- read.csv('~/Github/helminths/Data/associations cleaned.csv')
associations <- unique(associations[,2:3])

h1 <- data.frame(table(helmDF$Parasite))
h2 <- data.frame(table(associations$Parasite))

h3 <- dplyr::left_join(h1, h2, by="Var1")

plot(h3[,3:2], xlab='Old scrape', ylab='New scrape')
hist(h3[,2]-h3[,3])

length(unique(helmDF$Parasite[!(helmDF$Parasite %in% associations$Parasite)]))
