
# Prepare MCC-data for Epi analysis
# author: Samuel Lüthi
# date: 18 October 2022

##### 1. LOAD MCC DATASET
# DEFINE DIRECTORY DATASET
dirpre <- "~/OneDrive - ETH Zurich/WCR/Projects/Heat/Data/MCC/"

# LOAD DATA (LAST VERSION UPDATED FRANCESCO)
load(paste(dirpre,"MCCdata_20211007.RData",sep="/"))

# DEFINE DIRECTORY TO OUTPUT
dirout <- "~/OneDrive - ETH Zurich/WCR/Projects/Heat/Data/MCC/clean_data/"

##### 2. SELECT COUNTRIES
# updated country list (13.12.21)
subcountry <- sort(c("arg0515","aus8809","bra9718","can8615","chl0414","chi9615",
                     "col9813","crc0017","cze9415", "ecu1418","est9718","fnl9414","fra0015", "fracar0015","fragui0015", "frareu0015", "ger9315", 
                     "grc0110","gua0916","irl8407","irn0215", "isr8520", "ita8710", "jap7215","kor9718","kuw0016","mex9814","mld0110","nor6918","net9516r","pan1316",
                     "par0419","per0814","phi0619","por8018", "pue0916","rom9416", "sa9713","spa9014", "sui9513","swe9016","tha9908","twn9414","uk9016",
                     "uru1216","usa7306","vie0913"))

dlist <- dlist[cities$country%in%subcountry]
cities <- cities[cities$country%in%subcountry,]
countries <- countries[countries$country%in%subcountry,]

names(dlist) <- cities$city

##### 3. EXCLUDE STRANGE LOCATIONS (BASED ON PRELIMINARY ANALYSIS)

# EXCLUDE TANGSHAN IN CHINA, WITH WEIRD DEATH COUNTS DISTRIBUTION
dlist[cities$city=="tngs.chi9615"] <- NULL
cities <- cities[-which(cities$city=="tngs.chi9615"),]

# EXCLUDE NANJING IN CHINA, WITH NO COUNTS FOR ALL/NONEXTERNAL MORTALITY
dlist[cities$city=="nnjn.chi9615"] <- NULL
cities <- cities[-which(cities$city=="nnjn.chi9615"),]

# EXCLUDE 
dlist[cities$city=="amjb.sa9713"] <- NULL
cities <- cities[-which(cities$city=="amjb.sa9713"),]

# EXCLUDE 
dlist[cities$city=="jhtg.sa9713"] <- NULL
cities <- cities[-which(cities$city=="jhtg.sa9713"),]

# EXCLUDE 
dlist[cities$city=="sdbn.sa9713"] <- NULL
cities <- cities[-which(cities$city=="sdbn.sa9713"),]

# EXCLUDE 
dlist[cities$city=="xhrp.sa9713"] <- NULL
cities <- cities[-which(cities$city=="xhrp.sa9713"),]

# EXCLUDE 
dlist[cities$city=="zlln.sa9713"] <- NULL
cities <- cities[-which(cities$city=="zlln.sa9713"),]

# EXCLUDE 
dlist[cities$city=="wstr.sa9713"] <- NULL
cities <- cities[-which(cities$city=="wstr.sa9713"),]

# EXCLUDE 
dlist[cities$city=="ilmb.sa9713"] <- NULL
cities <- cities[-which(cities$city=="ilmb.sa9713"),]

# EXCLUDE 
dlist[cities$city=="riet.ita8710"] <- NULL
cities <- cities[-which(cities$city=="riet.ita8710"),]

# EXCLUDE 
dlist[cities$city=="namp.usa7306"] <- NULL
cities <- cities[-which(cities$city=="namp.usa7306"),]

# EXCLUDE 
dlist[cities$city=="phbn.tha9908"] <- NULL
cities <- cities[-which(cities$city=="phbn.tha9908"),]

# EXCLUDE 
dlist[cities$city=="gngz.chi9615"] <- NULL
cities <- cities[-which(cities$city=="gngz.chi9615"),]

# INCLUDE (fix region for France Reunion)
levels(cities$region) <- c(levels(cities$region), "Indian Ocean")
cities$region[cities$countryname == "FranceReunion"] <- "Indian Ocean"

# SET OUTLIERS TO MISSING
for(i in seq(nrow(cities))) {
  if(is.null(dlist[[i]]$all)) dlist[[i]]$nonext[dlist[[i]]$outlierm==1] <- NA else
    dlist[[i]]$all[dlist[[i]]$outlierm==1] <- NA
  dlist[[i]]$tmean[dlist[[i]]$outliert==1] <- NA
}

rownames(cities) <- seq(nrow(cities))
cities[1:4] <- as.data.frame(as.matrix(cities[1:4]))
rownames(countries) <- seq(nrow(countries))

# SAVE DATA
save.image(paste(dirout,"MCCdata_clean_220103.RData",sep="/"))
