
load('C:/Users/cjcar/Documents/Github/helminthLDR/analysis/figshare/ldr.Rdata')


sub.acan <- data.frame(group='',Host='',Parasite='',ParasiteFull='',locality='')
sub.acan <- sub.acan[-1,]
for (i in 1:length(acans)) {
  if(!is.null(acans[[i]])) {
    sub.sub.acan <- acans[[i]]
    sub.sub.acan$group='Acanthocephala'
    sub.sub.acan$locality=names(acans)[[i]]
    sub.acan <- rbind(sub.acan,sub.sub.acan)
  }
}


sub.cest <- data.frame(group='',Host='',Parasite='',ParasiteFull='',locality='')
sub.cest <- sub.cest[-1,]
for (i in 1:length(cests)) {
  if(!is.null(cests[[i]])) {
    sub.sub.cest <- cests[[i]]
    sub.sub.cest$group='Cestoda'
    sub.sub.cest$locality=names(cests)[[i]]
    sub.cest <- rbind(sub.cest,sub.sub.cest)
  }
}


sub.nema <- data.frame(group='',Host='',Parasite='',ParasiteFull='',locality='')
sub.nema <- sub.nema[-1,]
for (i in 1:length(nemas)) {
  if(!is.null(nemas[[i]])) {
    sub.sub.nema <- nemas[[i]]
    sub.sub.nema$group='Nematoda'
    sub.sub.nema$locality=names(nemas)[[i]]
    sub.nema <- rbind(sub.nema,sub.sub.nema)
  }
}


sub.trem <- data.frame(group='',Host='',Parasite='',ParasiteFull='',locality='')
sub.trem <- sub.trem[-1,]
for (i in 1:length(trems)) {
  if(!is.null(trems[[i]])) {
    sub.sub.trem <- trems[[i]]
    sub.sub.trem$group='Trematoda'
    sub.sub.trem$locality=names(trems)[[i]]
    sub.trem <- rbind(sub.trem,sub.sub.trem)
  }
}

associations <- rbind(sub.acan,sub.cest,sub.nema,sub.trem)

## CLEANING TIME ##

associations <- associations[!is.na(associations$Host),]
associations <- associations[!is.na(associations$Parasite),]

associations <- associations[!associations$Host=='',]
associations <- associations[!associations$Parasite=='',]

spex <- function(x) {return(!(length(grep("sp.",x))==1))}

sps <- apply(associations,c(1,2),spex)
associations <- associations[sps[,1]==TRUE,]
sps <- apply(associations,c(1,2),spex)
associations <- associations[sps[,2]==TRUE,]

associations <- associations[!(substr(associations$ParasiteFull,1,1)=="("),]
associations <- associations[!(substr(associations$Host,1,1)=="("),]

twoname <- function(x) {
  y <- strsplit(x," ")
  return(length(y[[1]][!y[[1]]==''])>1)
}

sps <- apply(associations,c(1,2),twoname)
associations <- associations[sps[,1]==TRUE,]
sps <- apply(associations,c(1,2),twoname)
associations <- associations[sps[,2]==TRUE,]

revis <- function(name) {
  
  list <- strsplit(name,' ')[[1]]
  
  if(substr(list[2],1,1)=='(') {
    return(paste(list[1],list[3]))
  } else {
    if(substr(list[2],1,1)=='['){
      return(paste(list[1],list[3]))
    } else {
      return(paste(list[1],list[2]))
    }
  }
}

for (i in 1:nrow(associations)) {
  associations$Parasite[i] <- revis(associations$ParasiteFull[i])
}

for (i in 1:nrow(associations)) {
  associations$Host[i] <- revis(associations$Host[i])
}


sps <- apply(associations,c(1,2),spex)
associations <- associations[sps[,2]==TRUE,]

spq <- function(x) {return(!(length(grep("\\?",x))==1))}
sps <- apply(associations,c(1,2),spq)
associations <- associations[sps[,1]==TRUE,]
sps <- apply(associations,c(1,2),spq)
associations <- associations[sps[,2]==TRUE,]





classifier <- function(name) { 
  name2 <- tryCatch(gnr_resolve(name, best_match_only = TRUE)$matched_name,error=function(e) {NULL})
  if(is.null(name2)) {
    return(c(name,NA))
  } else {
    ids <- tryCatch(suppressMessages(get_tsn_(name2,accepted=FALSE)[[1]]$tsn),error=function(e) {NA})
    if (length(ids)==0) { return (c(name,NA)) } else {
      ids <- ids[[1]]
      ids <- as.numeric(as.character(ids))
      x <- classification(ids, db = 'itis')[[1]]
      if(nrow(x)==1) {
        ids <- tryCatch(synonyms(ids,db='itis')[[1]]$acc_tsn[1], error=function(e) {NA})
        if(is.na(ids)) {
          return(c(name2,NA))
        } else {
          x <- classification(ids, db = 'itis')[[1]]
        }
      }
      if(nrow(x[x$rank=='class',])>0) {
        species <- x[x$rank=='species',]$name
        class <- x[x$rank=='class',]$name
        return(c(species,class))
      } else {
        return(c(name2,NA))
      }
    }
  }
}

dictionary <- data.frame(hosts = unique(associations$Host), hostclean='', group='')
dictionary$hosts <- as.character(dictionary$hosts)
dictionary$hostclean <- as.character(dictionary$hostclean)
dictionary$group <- as.character(dictionary$group)

for(i in 1:nrow(dictionary)){
  pb <- txtProgressBar(min = 0, max = nrow(dictionary), style = 3)
  dictionary[i,2:3] <- classifier(dictionary$hosts[i])
  setTxtProgressBar(pb, i)
}

write.csv(dictionary,'dictionary.csv')
write.csv(associations,'associations uncleaned.csv')

## dictionary <- read.csv('~/Github/helminths/Data/dictionary.csv')

# associations <- read.csv('~/Github/helminths/Data/associations uncleaned.csv',
#                         stringsAsFactors = FALSE)[,-1]


associations$hostgroup <- NA

for (i in 1:nrow(associations)){
  associations$Host[i] <- dictionary[dictionary$hosts==associations$Host[i],]$hostclean
  
  if(!length(dictionary$group[dictionary$hosts==associations$Host[i]])==0) {
    
    associations$hostgroup[i] <- as.character(dictionary[dictionary$hosts==associations$Host[i],]$group)
    
  }
  print(i)
}

write.csv(associations,'associations cleaned.csv')