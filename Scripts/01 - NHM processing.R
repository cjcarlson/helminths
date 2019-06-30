library(plyr)

associations <- read.csv('~/Github/helminths/Data/associations cleaned.csv')
associations <- associations[,-1]

associations$genus <- ''
associations$year <- ''

for(i in 1:nrow(associations)) {

x <- strsplit(as.character(associations$ParasiteFull[i]),' ')[[1]]
associations$year[i] <- gsub(")","",x[length(x)])
associations$year[i] <- gsub("\\(","",associations$year[i])

x <- strsplit(as.character(associations$Parasite[i]),' ')[[1]]
associations$genus[i] <- gsub(" ","",x[1])

print(i)
}

genera <- unique(associations[,c(2,7,8)])



generalism <- data.frame(table(unique(associations[,c(1,2,4,6)])[,c(2,3,4)]))
generalism <- generalism[!generalism$Freq==0,]
generalism$hostgroup <- as.character(generalism$hostgroup)
generalism[generalism$hostgroup %in% c("Chondrostei", "Holostei", "Cladistei", "Teleostei"),]$hostgroup <- "Osteichthys" 


hostpar <- associations[,c('Host','Parasite')]
hostpar <- unique(hostpar)
thostpar <- table(hostpar$Parasite)
genera$hs <- 0


#  genera$hs[i] <- sum(generalism[generalism$Parasite==genera$Parasite[i],]$Freq)

for (i in 1:nrow(genera)) {
  genera$hs[i] <- thostpar[genera$Parasite[i]][[1]]
  print(i)
}

genera$year <- as.numeric(genera$year)
genera$year[genera$year>2018] <- NA
genera$year[genera$year<1700] <- NA



genlist <- unique(genera$genus)

genera$first <- ''

for (i in 1:length(genlist)) {
  working <- genera[genera$genus==genlist[i],]
  
  working <- na.omit(working)
  
  if(nrow(working)>0) {
  
  if(length(na.omit(working$year))>1) {
    genera[genera$genus==genlist[i],]$first <- 'not the first'
    genera[genera$Parasite %in% working[working$year==min(working$year),]$Parasite,]$first <- 'first'
    
    
  } else {
    genera[genera$genus==genlist[i],]$first <- 'exclude'
    
  }
    
  } else {
    genera[genera$genus==genlist[i],]$first <- 'exclude'
  }
}

write.csv(genera,'~/Github/helminths/Data/timeseries_NHM.csv')
