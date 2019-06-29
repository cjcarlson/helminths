library(plyr)
library(helminthR)
library(mgcv)

setwd('C:/Users/cjcar/Dropbox/helminths/')

npc <- read.csv('USNPC September 2017.csv')

yearfinder <- function(str) {
  str <- gsub('/',' ',str)
  str <- gsub('\\.',' ',str)
  str <- gsub('-',' ',str)
  candidates <- strsplit(str,' ')
  candidates <- candidates[[1]][nchar(candidates[[1]])==4]
  list <- c(as.numeric(as.character(candidates)))
  list <- list[!is.na(list)]
  if(length(list)==0) return(NA) else {
    list <- list[list>1500]
    list <- list[list<2020]
    return(year <- min(list))
  }
}

for (i in 1:nrow(npc)) {
  npc$colyear[i] <- as.numeric(as.character(yearfinder(npc$ColDateVisitedFrom[i])))
  npc$dbyear[i] <- as.numeric(as.character(yearfinder(npc$CatDateCataloged[i])))
}

# So, we want to use colyear

tax.filter <- c('Platyhelminthes','Acanthocephala','Nematoda')

npc.sub <- npc[npc$ClaPhylum %in% tax.filter, ]

nameparse <- function(str) {
  str <- as.character(str)
  if(!is.na(str)) {
    candidates <- strsplit(str,split=' ')
    name <- paste(candidates[[1]][1],candidates[[1]][2], sep=' ')
    return(name)
  }
}

for(i in 1:nrow(npc.sub)){ 
  npc.sub$nameparse[i] <- nameparse(npc.sub$SummaryData[i])
}

npc.sub <- npc.sub[,c('nameparse','colyear')]
npc.sub <- npc.sub[!is.na(npc.sub$colyear),]

# Year of first appearance

yofa <- plyr::ddply(npc.sub, c('nameparse'), summarize,
                    yearone = min(colyear))
yofa <- yofa[-1,]

# Cumulative

totals <- data.frame(table(yofa$yearone))
totals$cumsum <- cumsum(totals$Freq)
totals$Var1 <- as.numeric(as.character(totals$Var1))
plot(totals$cumsum ~ totals$Var1)

# How many hosts?

howmanyhosts <- function(name) {
  g <- strsplit(name,' ')[[1]][1]
  s <- strsplit(name,' ')[[1]][2]
  hosts <- tryCatch(findParasite(genus=g,species=s,speciesOnly=TRUE,removeDuplicates=TRUE),error=function(e) data.frame(NA))
  if(nrow(hosts)==0) {
    return(NA)
  } else {
    return(nrow(hosts))
  }
}
howmanyhosts('Taenia solium')
howmanyhosts('Taenia smurfium')

yofa$hosts <- 0

for(i in 1:nrow(yofa)) {
  yofa$hosts[i] <- howmanyhosts(yofa$nameparse[i])
  if(i %in% (c(1:130)*100)) {print(i)}
}

write.csv(yofa, '~/Github/helminths/Data/timeseries_NPC.csv')