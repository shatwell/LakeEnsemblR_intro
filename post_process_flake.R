# post process flake results


## set working directory to location the R script is stored - only in RStudio
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Rappbode reservoir ------------------------------------------------------

modname <- "rappbode_flake.rslt"
obsname <- "rappbode_temp.csv"

## read results file. Note we have to skip the first line of the text file
dat <- read.table(paste0("flake/", modname), skip=1, header=TRUE)


## Clean up the model data
spinup <- 365 # The model was run with a 365 day spinup period
dat <- dat[spinup+1:(nrow(dat)-spinup),] # remove the spinup period
start <- as.POSIXct("1979-01-01", tz="UTC") # start of the simulation period
dat$time <- 0:(nrow(dat)-1) # reset the time
dat$date <- dat$time*3600*24+start # add date column


## plot the results
plot(Ts~date,dat,type="l")
lines(Tb~date,dat)

## Compare with observed data
### read observed data
obs <- read.csv(paste0("flake/", obsname))
obs$date <- as.POSIXct(obs$date,tz="UTC")

### plot the results for 2015 and 2016
plot(m0~date,obs, col="red")
lines(Ts~date,dat, col="red")
points(m25~date,obs,col="blue")
lines(Tb~date,dat,col="blue")

### calculate mixed layer depth
library(rLakeAnalyzer)
# Loop through all the profiles and calculate MLD for each one using the rLakeAnalyzer function thermo.depth()
hmix <- NULL
for(ii in 1:nrow(obs)) {
  hmix <- c(hmix, 
            thermo.depth(wtr = as.numeric(obs[ii,-1]), 
                         depths = 0:50)  
  )
}

### plot mixed layer depth
plot(hmix~obs$date, ylim=c(30,0))
lines(h_ML~date, dat)

# Lake Erken --------------------------------------------------------------

### Hint: you can also run flake from within R instead of a shell:
# setwd("flake")
# system2(command = "flake", args = "erken.nml") # only Windows
# system2(command = "./flake", args = "erken.nml") # only Linux or Mac
# setwd("../")

modname <- "erken_flake.rslt"
obsname <- "Erken_temp_daily.csv"

## read results file. Note we have to skip the first line of the text file
dat <- read.table(paste0("flake/", modname), skip=1, header=TRUE)


## Clean up the model data
spinup <- 365 # The model was run with a 365 day spinup period
dat <- dat[spinup+1:(nrow(dat)-spinup),] # remove the spinup period
start <- as.POSIXct("1979-01-01", tz="UTC") # start of the simulation period
dat$time <- 0:(nrow(dat)-1) # reset the time
dat$date <- dat$time*3600*24+start # add date column


## plot the results
plot(Ts~date,dat,type="l", col="red")
lines(Tb~date,dat, col="blue")

## Compare with observed data
### observed data
obs <- read.csv(paste0("flake/", obsname)) # read file
obs$date <- strptime(obs$TIMESTAMP, # parse timestamp into date POSIX format
                     format="%Y%m%d",
                     tz="UTC")
obs$year <- obs$date$year+1900 # add year column to make plotting easier
obs$date <- as.POSIXct(obs$date) # change to plottable date format

### plot the results for 2015 and 2016
plot(WTEMP~date,obs, subset=DEPTH<=1.0 & year %in% 2015:2016, 
     col="red", cex=0.5)
lines(Ts~date,dat, col="red3", lwd=1.5)
points(WTEMP~date,obs,subset=DEPTH==15, col="lightblue", cex=0.5)
lines(Tb~date,dat,col="blue3", lwd=1.5)
legend("topleft", c("Ts","Tb"), col=c("red3","blue3"), 
       lwd=1.5, bty="n")

