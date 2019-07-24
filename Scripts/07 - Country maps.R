library(ggplot2)
library(rgdal)
library(sf)

# RUN "DIVERSITY ESTIMATES.R" FIRST

library(codependent)
library(tidyverse)
library(reshape2)

setwd('~/Github/helminths')
all.data <- read_delim('./data/associations cleaned.csv',delim=',')[,-1]

## CLEAN UP THE GEOGRAPHY ONCE AND FOR ALL 


america.names <- c('Alabama', 'Alaska', 'Arizona', 'Arkansas', 'Colorado',
                   "District of Columbia/Delaware/New Jersey/Maryland",
                   'Florida', 'Georgia (USA)', 'Hawaiian Islands', 'Idaho',
                   'Illinois','Indiana', 'Kansas', 'Kentucky', 'Louisiana',
                   "Maine",'Michigan','Minnesota','Missouri','Montana',"Mississippi",
                   "Nebraska", "New Hampshire/Vermont/Massachusetts/Connecticut/Rhode Is.",
                   'Ohio', 'Oklahoma', 'Oregon', 'Pacific coast of USA', "Pennsylvania",
                   "Puerto Rico", "Tennessee", 'Texas', 'Virginia', 'Washington',
                   "West Virginia", 'Wisconsin', 'Wyoming', "Nevada", "New Mexico", "North Dakota",
                   'New York State', "Utah", "California", "Iowa","North Carolina",
                   "South Carolina", "South Dakota")
all.data$locality[all.data$locality %in% america.names] <- 'United States'

canada.names <- c('Alberta', 'Labrador', 'Labrador + Newfoundland', 'Manitoba',
                  "New Brunswick","New Brunswick + Nova Scotia", "Newfoundland", 
                  "Northwest Territories", "Nova Scotia", "Ontario", "Quebec",
                  "Saskatchewan", "Yukon", "Yukon + Northwest Territories",
                  "British Columbia")
all.data$locality[all.data$locality %in% canada.names] <- 'Canada'

aus.names <- c("New South Wales", "Western Australia", 'Queensland',
               "South Australia", 'Victoria')
all.data$locality[all.data$locality %in% aus.names] <- 'Australia'

all.data <- separate_rows(all.data, locality, sep='\\+', convert=TRUE)
unique(all.data$locality)

all.data$locality <- recode(all.data$locality,
                            Azerbaidzan = 'Azerbaijan',
                            Azores = 'Portugal',
                            ' Qatar ' = 'Qatar',
                            ' United Arab Emirates' = 'United Arab Emirates',
                            "Balearic Islands" = "Spain",
                            " Luxembourg" = "Luxembourg",
                            'Belgium ' = "Belgium",
                            "Borneo" = "Malaysia+Brunei+Indonesia",
                            "Byelorussia" = "Belarus",
                            "Central Africa"='Central African Republic',
                            "Corsica" = "France",
                            "Eire " = "Northern Ireland",
                            " N. Ireland" = "Northern Ireland",
                            "England, Wales " = "United Kingdom",
                            ' Isle of Man'="United Kingdom",
                            " Latvia"="Latvia",
                            "Fijian Islands"='Fiji',
                            "Georgia (USSR)"='Georgia',
                            "Iraq "="Iraq",
                            " Kuwait"="Kuwait",
                            "Israel "="Israel",
                            " Jordan " = "Jordan",
                            " Lebanon"="Lebanon",
                            'Kirgizia'='Kyrgyzstan',
                            "Mali "='Mali',
                            ' Burkina Faso'='Burkina Faso',
                            "Natal"='South Africa',
                            "Sardinia"="Italy",
                            "Senegal "='Senegal',
                            ' The Gambia '='The Gambia',
                            ' Guinea-Bassau'='Guinea-Bassau',
                            "Spain "='Spain',
                            ' Andalusia'="Spain",
                            'Switzerland '='Switzerland',
                            " Lichtenstein"="Lichtenstein",
                            'Tadzhikistan'='Tajikistan',
                            'Uganda '='Uganda',
                            ' Burundi '='Burundi',
                            ' Rwanda'="Rwanda",
                            "Ukraine, incl. Moldavia"='Ukraine',
                            "Benin "='Benin',
                            ' Togo '='Togo',
                            ' Ghana'='Ghana',
                            'Burma [Myanmar]'='Myanmar',
                            ' Haiti'="Haiti",
                            'Dominican Republic '="Dominican Republic",
                            'Congo'='Republic of Congo',
                            "Ethiopia (incl. Dhibouti)"='Ethiopia+Djibouti',
                            'Gabon '='Gabon',
                            ' Congo '='Republic of Congo',
                            ' Equatorial Guinea'='Equatorial Guinea',
                            'Guinea '='Guinea',
                            'Sierra Leone '='Sierra Leone',
                            ' Liberia'='Liberia',
                            'Sabah'='Malaysia',
                            'Sicily '="Italy",
                            ' Malta'="Italy",
                            'Society Islands (incl. Tahiti)'='France',
                            "Southern Yemen"='Yemen',
                            "Southern Yemen "="Yemen",
                            " Yemen"="Yemen",
                            "Spitzbergen (Svalbard)"="Norway",
                            'Sumatra'='Indonesia',
                            "Transvaal"='South Africa',
                            "Galapagos Islands"='Ecuador',
                            'Gibraltar'='Spain',
                            'Orange Free State'='South Africa',
                            'Kampuchea' = "Cambodia",
                            "Java" = "Indonesia",
                            "Finno-Karelian ASSR" ='Russia',
                            "Russia (Asian)" = "Russia",
                            "Russia (European)" = "Russia",
                            'Yakut ASSR' = "Russia")

all.data <- unique(separate_rows(all.data, locality, sep='\\+', convert=TRUE))


#########

load("C:/Users/cjcar/Dropbox/HowManyHelminths2019/mammalRichnessReduced.RData")

iso$Name[iso$Name=="Congo, the Democratic Republic of the"] <- "Democratic Republic of the Congo"
iso$Name[iso$Name=="Iran, Islamic Republic of"] <- "Iran"
iso$Name[iso$Name=="Korea, Republic of"] <- "South Korea"
iso$Name[iso$Name=="Lao People's Democratic Republic"] <- "Laos"
iso$Name[iso$Name=="Macedonia, the former Yugoslav Republic of"] <- 'Macedonia'
iso$Name[iso$Name=="Korea, Democratic People's Republic of" ] <- "North Korea"
iso$Name[iso$Name=="Russian Federation"] <- 'Russia' 
iso$Name[iso$Name=="Venezuela, Bolivarian Republic of"] <- 'Venezuela' 
iso$Name[iso$Name=="Viet Nam"] <- 'Vietnam'



### Pred 
pred.2 <- function(model, n.indep) {
  q <- stats::coef(model)
  est <- q["b"] * (n.indep)^(q["z"])
  return(est)
}
b.a <- binera(acan.mamm, iter=100)
b.c <- binera(cest.mamm, iter=100)
b.n <- binera(nema.mamm, iter=100)
b.t <- binera(trem.mamm, iter=100)

iso$ac <- sapply(iso$mammalRichness, pred.2, model=b.a)
iso$ce <- sapply(iso$mammalRichness, pred.2, model=b.c)
iso$ne <- sapply(iso$mammalRichness, pred.2, model=b.n)
iso$tr <- sapply(iso$mammalRichness, pred.2, model=b.t)
iso$helminths <- iso$ac + iso$ce + iso$ne + iso$tr

####

iso$ac.known <- 0
all.data %>% filter(group=='Acanthocephala') %>% filter(hostgroup=='Mammalia') %>% 
  data.frame() %>% unique() -> sub.df

for (i in 1:length(iso$Name)) {
  name <- iso$Name[i]
  if(name %in% sub.df$locality) {
    iso$ac.known[i] <- length(unique(sub.df[sub.df$locality==name,'Parasite']))
  }
}

iso$ce.known <- 0
all.data %>% filter(group=='Cestoda') %>% filter(hostgroup=='Mammalia') %>% 
  data.frame() %>% unique() -> sub.df

for (i in 1:length(iso$Name)) {
  name <- iso$Name[i]
  if(name %in% sub.df$locality) {
    iso$ce.known[i] <- length(unique(sub.df[sub.df$locality==name,'Parasite']))
  }
}

iso$ne.known <- 0
all.data %>% filter(group=='Nematoda') %>% filter(hostgroup=='Mammalia') %>% 
  data.frame() %>% unique() -> sub.df

for (i in 1:length(iso$Name)) {
  name <- iso$Name[i]
  if(name %in% sub.df$locality) {
    iso$ne.known[i] <- length(unique(sub.df[sub.df$locality==name,'Parasite']))
  }
}

iso$tr.known <- 0
all.data %>% filter(group=='Trematoda') %>% filter(hostgroup=='Mammalia') %>% 
  data.frame() %>% unique() -> sub.df

for (i in 1:length(iso$Name)) {
  name <- iso$Name[i]
  if(name %in% sub.df$locality) {
    iso$tr.known[i] <- length(unique(sub.df[sub.df$locality==name,'Parasite']))
  }
}

iso$helm.known <- iso$ac.known + iso$ce.known + iso$ne.known + iso$tr.known
iso$ac.desc <- round((iso$ac-iso$ac.known)/(iso$ac)*100)
iso$ce.undesc <- round((iso$ce-iso$ce.known)/(iso$ce)*100)
iso$ne.undesc <- round((iso$ne-iso$ne.known)/(iso$ne)*100)
iso$tr.undesc <- round((iso$tr-iso$tr.known)/(iso$tr)*100)
iso$helm.undesc <- round((iso$helminths-iso$helm.known)/(iso$helminths)*100)

#sp <- as_Spatial(st_cast(iso))
#writeOGR(sp, dsn='.', layer='helminths2.shp', driver='ESRI Shapefile')
#save.image("C:/Users/cjcar/Dropbox/HowManyHelminths2019/Maps workspace backup.RData")

g1 <- ggplot(data = iso) +
  geom_sf(aes(fill = helm.undesc)) + coord_sf(datum = NA) +
  scale_fill_viridis_c(option = "plasma") + 
  labs(fill = "Undocumented") + theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.background = element_blank())
ggsave("Figure4c.pdf", plot=g1, width = 11, height = 5)

g2 <- ggplot(data = iso) +
  geom_sf(aes(fill = helminths)) + coord_sf(datum = NA) +
  scale_fill_viridis_c(option = "plasma", trans='sqrt') + theme_void() + 
  labs(fill = "Possible    ") + theme(panel.background = element_blank())
ggsave("Figure4a.pdf", plot=g2, width = 11, height = 5)

g3 <- ggplot(data = iso) +
  geom_sf(aes(fill = helm.known)) + coord_sf(datum = NA) +
  scale_fill_viridis_c(option = "plasma") + theme_void() + 
  labs(fill = "Known      ") + theme(panel.background = element_blank())
ggsave("Figure4b.pdf", plot=g3, width = 11, height = 5)
