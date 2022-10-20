# Prepare MCC-data for Epi analysis
# author: Samuel Lüthi & Ana Vicedo (original script)
# date: 18 October 2022

### MAIN ANALYSIS
# In this code, we estimate the exposure-response function for each
# location using the two-stage approach.
# The specifications are the ones of the temp var paper of Antonio.

#### 1. LOAD LIBRARIES 
library(dlnm) ; library(mixmeta) ; library(splines) ; library(tsModel)
library(mgcv) ; library(lubridate)


#### 2. LOAD MCC DATASET (PREPARED) AND SET WORKING DIRECTORIES
# DEFINE DIRECTORY DATASET
dirpre <- "~/OneDrive - ETH Zurich/WCR/Projects/Heat/Data/MCC/"

# DEFINE DIRECTORY TO OUTPUT
dirout <- "~/OneDrive - ETH Zurich/WCR/Projects/Heat/Data/MCC/clean_data/"

# LOAD DATA (LAST VERSION UPDATED JANUARY 2022)
load(paste(dirout,"MCCdata_clean_220103.RData",sep="/"))


#### 3. DEFINITION OF THE PARAMETERS
# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "ns"
vardegree <- NULL
varper <- c(10,75,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 21
lagnk <- 3

# DEGREE OF FREEDOM FOR TREND
dftrend <- 8

# COMPUTE PERCENTILES
per <- t(sapply(dlist,function(x) 
  quantile(x$tmean,c(2.5,10,25,50,75,90,97.5)/100,na.rm=T)))

# CREATE INDICATOR FOR ONLY NON-EXTERNAL MORTALITY
indnonext <- sapply(dlist,function(x) !"all"%in%names(x))

# DEFINE THE OUTCOME
out <- "all"

# MODEL FORMULA
formula <- y ~ cb + dow + ns(date,df=round(dftrend*length(date)/365.25))



#### 4. RUN THE 1ST STAGE MODEL
# CREATE THE OBJECTS TO STORE THE RESULTS
ncoef <- length(varper) + ifelse(varfun=="bs",vardegree,1)
coefall <- matrix(NA,nrow(cities),ncoef,dimnames=list(cities$city))
vcovall <- vector("list",nrow(cities))
names(vcovall) <- cities$city

tmeanrange <- matrix(NA,nrow(cities),2)

for(i in seq(length(dlist))) {

  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  
  # DEFINE THE OUTCOME
  data$y <- if(indnonext[i]) as.integer(data$nonext) else data[[out]]
  
  # STORE RANGE TEMPERATURE
  tmeanrange[i,] <- range(data$tmean, na.rm=T)

  # DEFINE THE CROSSBASIS
  argvar <- list(fun=varfun,knots=quantile(data$tmean,varper/100,na.rm=T), 
                 Bound=range(data$tmean,na.rm=T))
  arglag <- list(knots=logknots(lag,lagnk))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  
  cb <- crossbasis(data$tmean,lag=lag,argvar=argvar,arglag=arglag) 
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")

  # REDUCTION TO OVERALL CUMULATIVE
  redall <- crossreduce(cb,model,cen=mean(data$tmean,na.rm=T))
  coefall[i,] <- coef(redall)
  vcovall[[i]] <- vcov(redall)

}

rownames(tmeanrange) <- names(dlist)

save.image(paste(dirout,"Objects/1ststage_run_220328.RData",sep="/"))

#### 5. SECOND STAGE
load(paste(dirout,"Objects/1ststage_run_220328.RData",sep="/"))

# CREATE AVERAGE TEMPERATURE AND RANGE AS META-PREDICTORS
avgtmean <- sapply(dlist,function(x) mean(x$tmean,na.rm=T))
rangetmean <- sapply(dlist,function(x) IQR(x$tmean,na.rm=T)) ## IQR!!
cntkgcl <- factor(paste(cities$country, cities$kgclzone,sep="-"))

# META-REGRESSION INCL. RANDOM TERM
mvmlall <- mixmeta(coefall ~ kgclzone1 + region + I(gdp/10000) + avgtmean + rangetmean,
random=~1|cntkgcl/city,vcovall,cities,control=list(showiter=T, igls.inititer=10),method="ml")
summary(mvmlall)

# OBTAIN BLUPS
blupall <- blup(mvmlall,vcov=T)

#### 6. ESTIMATE THE MMT
# GENERATE THE MATRIX FOR STORING THE RESULTS
minperccity <- mintempcity <- rep(NA,length(dlist))
names(mintempcity) <- names(minperccity) <- cities$city

# DEFINE MINIMUM MORTALITY VALUES: EXCLUDE LOW AND VERY HOT TEMPERATURE
for(i in seq(length(dlist))) {
  # EXTRACT THE DATA
  print(i)
  data <- dlist[[i]]
  predvar <- quantile(data$tmean,2:98/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
    knots=quantile(data$tmean,varper/100,na.rm=T),
    Bound=range(data$tmean,na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis,argvar)
  minperccity[i] <- (2:98)[which.min((bvar%*%blupall[[i]]$blup))]
  mintempcity[i] <- quantile(data$tmean,minperccity[i]/100,na.rm=T)
}

# OVERALL AND COUNTRY-SPECIFIC POINTS OF MINIMUM MORTALITY
(minperccountry <- round(tapply(minperccity,cities$country,median)))
(minperctot <- median(minperccity))
save.image(paste(dirout,"Objects/2ndstage_run_220328.RData",sep="/"))

#### 7. PLOT EXPOSURE-RESPONSE FROM BLUPS + 1ST STAGE (BLUE)
pdf(paste(dirout,"Plots/expresphist_selblups_ml_run_220328.pdf",sep="/"),width=9,height=13)
layout(matrix(seq(7*4),nrow=7,byrow=T))
par(mar=c(4,3.8,3,2.4),mgp=c(2.5,1,0),las=1)

for(i in seq(length(dlist))) {

  data <- dlist[[i]]
  argvar <- list(x=dlist[[i]]$tmean,fun=varfun,
    knots=quantile(data$tmean,varper/100,na.rm=T),
    Bound=range(data$tmean,na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis,argvar)
  pred1 <- crosspred(bvar,coef=coefall[i,],vcov=vcovall[[i]],
    model.link="log",cen=mintempcity[i],by=0.1, from=range(data$tmean, na.rm=T)[1],
    to=range(data$tmean, na.rm=T)[2])
  pred2 <- crosspred(bvar,coef=blupall[[i]]$blup,vcov=blupall[[i]]$vcov,
    model.link="log",cen=mintempcity[i],by=0.1, from=range(data$tmean, na.rm=T)[1],
    to=range(data$tmean, na.rm=T)[2])
  
  plot(pred2,ylim=c(0.5,2.5),yaxt="n",lab=c(6,5,7),xlab="Tmean",ylab="RR",
    main=cities$cityname[i],frame.plot=F)
  lines(pred1$predvar,pred1$allRRfit,col=4,lwd=1.5, lty="dashed")

  mtext(cities$countryname[i],cex=0.7,line=0)
  axis(1)
  axis(2,at=1:5*0.5)
}
dev.off()


pdf(paste(dirout,"Plots/expresphist_blups_ml_run_220328.pdf",sep="/"),width=9,height=13)
layout(matrix(seq(7*4),nrow=7,byrow=T))
par(mar=c(4,3.8,3,2.4),mgp=c(2.5,1,0),las=1)

for(i in seq(length(dlist))) {

  data <- dlist[[i]]
  argvar <- list(x=dlist[[i]]$tmean,fun=varfun,
    knots=quantile(data$tmean,varper/100,na.rm=T),
    Bound=range(data$tmean,na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis,argvar)
  pred2 <- crosspred(bvar,coef=blupall[[i]]$blup,vcov=blupall[[i]]$vcov,
    model.link="log",cen=mintempcity[i],by=0.1, from=range(data$tmean, na.rm=T)[1],
    to=range(data$tmean, na.rm=T)[2])
  plot(pred2,ylim=c(0.5,2.5),yaxt="n",lab=c(6,5,7),xlab="Tmean",ylab="RR",
    main=cities$cityname[i],frame.plot=F)
  mtext(cities$countryname[i],cex=0.7,line=0)
  axis(1)
  axis(2,at=1:5*0.5)
}
dev.off()



#### 8. SAVE RR (INCL EXTRAPOLATION) AND TEMPERATURE-MORTALITY DATA TO XLSX
# (simple way to transfer to Python for probabilistic analysis)

library(openxlsx)
# RR data
wb_rr <- createWorkbook()
wb_cdat <- createWorkbook()

for(i in seq(length(dlist))) {

  # PRINT
  cat(i,"")
  data <- dlist[[i]]
  argvar <- list(x=dlist[[i]]$tmean,fun=varfun,
    knots=quantile(data$tmean,varper/100,na.rm=T),
    Bound=range(data$tmean,na.rm=T))
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(onebasis,argvar)
  pred2 <- crosspred(bvar,coef=blupall[[i]]$blup,vcov=blupall[[i]]$vcov,
    model.link="log",cen=mintempcity[i],by=0.1, from=range(data$tmean, na.rm=T)[1]-3.0,
    to=range(data$tmean, na.rm=T)[2]+6.0) # we extrapolate up to 6 degrees

  # get RR only
  RR <- data.frame(pred2$predvar, unname(pred2$allRRfit),
    unname(pred2$allRRlow), unname(pred2$allRRhigh))
  names(RR) <- c("temp","RRfit","RRlow","RRhigh")
  addWorksheet(wb=wb_rr, sheetName = cities$city[i])
  writeData(wb_rr, sheet = i, RR)
  
  # daily mort and temp data
  if("all" %in% colnames(data)){
    cdat <- data[c("date", "all","tmean")]
  } else if ("nonext" %in% colnames(data)) {
    cdat <- data[c("date", "nonext","tmean")]
  }
  names(cdat) <- c("date","deaths","temp")
  addWorksheet(wb=wb_cdat, sheetName = cities$city[i])
  writeData(wb_cdat, sheet = i, cdat)
}

saveWorkbook(wb_rr, paste(dirout,"all_RR_run_220328.xlsx", sep=""), overwrite = FALSE)
saveWorkbook(wb_cdat, paste(dirout,"all_MORT_run_220328.xlsx", sep=""), overwrite = FALSE)

# write temperature of minimum mortality

city <- names(mintempcity)
TMM <- unname(mintempcity)
data_tmm <- data.frame(city, TMM)
print(data_tmm)
write.csv(data_tmm, paste(dirout,"all_TMM_220328.csv", sep=""))

#### 9. SUMMARY TABLE FOR APPENDIX
loc_summary <- cities
t_range <- matrix(NA,nrow(cities),3)
mort_range <- matrix(NA,nrow(cities),3)
year_range <- matrix(NA,nrow(cities),2)
tot_mort  <- matrix(NA,nrow(cities),1)

for(i in seq(nrow(cities))) {
  d <- dlist[[i]]
  t_range[i, 1] <- median(d$tmean,na.rm=T)
  t_range[i, 2] <- quantile(d$tmean, 0.25, na.rm=T)
  t_range[i, 3] <- quantile(d$tmean, 0.75, na.rm=T)
  
  m <- if(indnonext[i]) as.integer(d$nonext) else d[[out]]
  mort_range[i, 1] <- median(m,na.rm=T)
  mort_range[i, 2] <- quantile(m, 0.25, na.rm=T)
  mort_range[i, 3] <- quantile(m, 0.75, na.rm=T)
  
  tot_mort[i] <- sum(m, na.rm=T)
  
  year_range[i,1] <- min(d$year, na.rm=T)
  year_range[i,2] <- max(d$year, na.rm=T)
}

loc_summary['T_median'] <- t_range[,1]
loc_summary['T_25'] <- t_range[,2]
loc_summary['T_75'] <- t_range[,3]

loc_summary['Mort_median'] <- mort_range[,1]
loc_summary['Mort_25'] <- mort_range[,2]
loc_summary['Mort_75'] <- mort_range[,3]

loc_summary['year_start'] <- year_range[,1]
loc_summary['year_end'] <- year_range[,2]

loc_summary['tot_mort'] <- tot_mort
