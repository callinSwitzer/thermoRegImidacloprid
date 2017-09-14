#Manually navigate to working directory
setwd("/Users/james/Dropbox/Work/Neonicotinoids/thermoregulationExpts/dataForCallin")

kk <- 1
# morning <- c(4,8)
# evening <- c(18,22)
#Create empty table for separate temperature outputs
cnames <- c('colony', 'treatment', 'dayTime', 'day', 'cohort','airTPerf', 'broodTPerf')
thermPerfData <- as.data.frame(matrix(nrow = 1,ncol = length(cnames)))
colnames(thermPerfData) <- cnames

vis = 0
#Calculate therm performance 
	sp <- 32.5 # Define temperature set point
	
	
data <- read.csv('summaryDataC1.csv')
cohort <- 1
#Adjust time to count from 6 in the morning
treatList <- c('c', 't', 'c', 't', 't', 'c')
head(data)
#Remove transition times
# data$hour <- data$time - floor(data$time)
# morningInd <- data$hour > morning[1]/24 & data$hour < morning[2]/24
# eveningInd <- data$hour > evening[1]/24 & data$hour < evening[2]/24
# keeperInd <- !(morningInd | eveningInd)
# data <- data[keeperInd,]

data$time <- data$time + 6/24
data$day <- floor(data$time)
days <- unique(data$day)
data$hour <- data$time - data$day #Rewrite data$hour
#Calculate rounded brood temperature
int <- 0.5
data$ambientRnd <- round((data$ambient/int))*int
daytimeRange <- c(0, 16/24) #Slightly tricky here - this is now isolating the hours from 6 am to 8 pm, but adjusted relatively to new timing (with 0 being 6 am on day one)
for(i in 1:17){
	
	curDat <- subset(data, day == i)
		 dayInd <- curDat$hour > daytimeRange[1] & curDat$hour < daytimeRange[2]
	 
	 #Separate night and day chunks
	dayDat <- curDat[dayInd,]
	nightDat <- curDat[!dayInd,]
	
	for(j in 1:6){
	 bi <- paste('c', j, 'brood', sep = "")
	 ai <- paste('c', j, 'air', sep = "")

	
	#Daytime data
	{
	tmp <- aggregate(dayDat[,c(ai,bi, 'ambient')], list(dayDat$ambientRnd), mean)
	colony <- j
	treatment <- treatList[j]
	dayTime <- 1
	day <- i
	
	x <- tmp$ambient - sp
	y <- tmp[,ai] - sp
	airModel <- lm(y~x + 0)
	
	x <- tmp$ambient - sp
	y <- tmp[,bi] - sp
	broodModel <- lm(y~x + 0)
	
	tmp <- thermPerfData[1,]
	tmp$colony <- colony
	tmp$treatment <- treatment
	tmp[,1:5] <- c(colony, treatment, dayTime, day, cohort)
	tmp$airTPerf <- summary(airModel)$coefficients[1]
	tmp$broodTPerf <- summary(broodModel)$coefficients[1]

	 thermPerfData <- rbind(thermPerfData, tmp)
	}
	
	#Nightime data
	{
	tmp <- aggregate(nightDat[,c(ai,bi, 'ambient')], list(nightDat$ambientRnd), mean)
	colony <- j
	treatment <- treatList[j]
	dayTime <- 0
	day <- i
	
x <- tmp$ambient - sp
	y <- tmp[,ai] - sp
	airModel <- lm(y~x + 0)
	
	x <- tmp$ambient - sp
	y <- tmp[,bi] - sp
	broodModel <- lm(y~x + 0)
	
	tmp <- thermPerfData[1,]
	tmp$colony <- colony
	tmp$treatment <- treatment
	tmp[,1:5] <- c(colony, treatment, dayTime, day, cohort)
	tmp$airTPerf <- summary(airModel)$coefficients[1]
	tmp$broodTPerf <- summary(broodModel)$coefficients[1]

	 thermPerfData <- rbind(thermPerfData, tmp)
	}
	
	}
	
# tmp <- aggregate(dayDat, list(dayDat$ambientRnd), mean)
# thermPerfData <- rbind(thermPerfData, tmp)
	
}


data <- read.csv('summaryDataC2.csv')
cohort <- 2
treatList <- c('c', 't', 't', 't', 'c', 'c')
#Remove transition times
# data$hour <- data$time - floor(data$time)
# morningInd <- data$hour > morning[1]/24 & data$hour < morning[2]/24
# eveningInd <- data$hour > evening[1]/24 & data$hour < evening[2]/24
# keeperInd <- !(morningInd | eveningInd)
# data <- data[keeperInd,]
#Adjust time to count from 6 in the morning
data$time <- data$time + 6/24
data$day <- floor(data$time)
days <- unique(data$day)
data$hour <- data$time - data$day
#Calculate rounded brood temperature
data$ambientRnd <- round((data$ambient/int))*int
daytimeRange <- c(0, 16/24) #Slightly tricky here - this is now isolating the hours from 6 am to 8 pm, but adjusted relatively to new timing (with 0 being 6 am on day one)
for(i in 1:9){
	
	curDat <- subset(data, day == i)
		 dayInd <- curDat$hour > daytimeRange[1] & curDat$hour < daytimeRange[2]
	 
	 #Separate night and day chunks
	dayDat <- curDat[dayInd,]
	nightDat <- curDat[!dayInd,]
	
	for(j in 1:6){
	 bi <- paste('c', j, 'brood', sep = "")
	 ai <- paste('c', j, 'air', sep = "")

	
	#Daytime data
	{
	tmp <- aggregate(dayDat[,c(ai,bi, 'ambient')], list(dayDat$ambientRnd), mean)
	colony <- j
	treatment <- treatList[j]
	dayTime <- 1
	day <- i
	
	x <- tmp$ambient - sp
	y <- tmp[,ai] - sp
	airModel <- lm(y~x + 0)
	
	x <- tmp$ambient - sp
	y <- tmp[,bi] - sp
	broodModel <- lm(y~x + 0)
	
	tmp <- thermPerfData[1,]
	tmp$colony <- colony
	tmp$treatment <- treatment
	tmp[,1:5] <- c(colony, treatment, dayTime, day, cohort)
	tmp$airTPerf <- summary(airModel)$coefficients[1]
	tmp$broodTPerf <- summary(broodModel)$coefficients[1]

	 thermPerfData <- rbind(thermPerfData, tmp)
	}
	
	#Nightime data
	{
	tmp <- aggregate(nightDat[,c(ai,bi, 'ambient')], list(nightDat$ambientRnd), mean)
	colony <- j
	treatment <- treatList[j]
	dayTime <- 0
	day <- i
	
x <- tmp$ambient - sp
	y <- tmp[,ai] - sp
	airModel <- lm(y~x + 0)
	
	x <- tmp$ambient - sp
	y <- tmp[,bi] - sp
	broodModel <- lm(y~x + 0)
	
	tmp <- thermPerfData[1,]
	tmp$colony <- colony
	tmp$treatment <- treatment
	tmp[,1:5] <- c(colony, treatment, dayTime, day, cohort)
	tmp$airTPerf <- summary(airModel)$coefficients[1]
	tmp$broodTPerf <- summary(broodModel)$coefficients[1]

	 thermPerfData <- rbind(thermPerfData, tmp)
	}
	
	}
# tmp <- aggregate(dayDat, list(dayDat$ambientRnd), mean)
# thermPerfData <- rbind(thermPerfData, tmp)
	
}



data <- read.csv('summaryDataC3.csv')
cohort <- 3
treatList <- c('t', 'c', 'c', 'c', 't', 't')
# #Remove transition times
# data$hour <- data$time - floor(data$time)
# morningInd <- data$hour > morning[1]/24 & data$hour < morning[2]/24
# eveningInd <- data$hour > evening[1]/24 & data$hour < evening[2]/24
# keeperInd <- !(morningInd | eveningInd)
# data <- data[keeperInd,]
#Adjust time to count from 6 in the morning
data$time <- data$time + 6/24
data$day <- floor(data$time)
days <- unique(data$day)
data$hour <- data$time - data$day
#Calculate rounded brood temperature
data$ambientRnd <- round((data$ambient/int))*int
daytimeRange <- c(0, 16/24) #Slightly tricky here - this is now isolating the hours from 6 am to 8 pm, but adjusted relatively to new timing (with 0 being 6 am on day one)
for(i in 1:9){
	
	curDat <- subset(data, day == i)
		 dayInd <- curDat$hour > daytimeRange[1] & curDat$hour < daytimeRange[2]
	 
	 #Separate night and day chunks
	dayDat <- curDat[dayInd,]
	nightDat <- curDat[!dayInd,]
	
	for(j in 1:6){
	 bi <- paste('c', j, 'brood', sep = "")
	 ai <- paste('c', j, 'air', sep = "")

	
	#Daytime data
	{
	tmp <- aggregate(dayDat[,c(ai,bi, 'ambient')], list(dayDat$ambientRnd), mean)
	colony <- j
	treatment <- treatList[j]
	dayTime <- 1
	day <- i
	
	x <- tmp$ambient - sp
	y <- tmp[,ai] - sp
	airModel <- lm(y~x + 0)
	
	x <- tmp$ambient - sp
	y <- tmp[,bi] - sp
	broodModel <- lm(y~x + 0)
	
	tmp <- thermPerfData[1,]
	tmp$colony <- colony
	tmp$treatment <- treatment
	tmp[,1:5] <- c(colony, treatment, dayTime, day, cohort)
	tmp$airTPerf <- summary(airModel)$coefficients[1]
	tmp$broodTPerf <- summary(broodModel)$coefficients[1]

	 thermPerfData <- rbind(thermPerfData, tmp)
	}
	
	#Nightime data
	{
	tmp <- aggregate(nightDat[,c(ai,bi, 'ambient')], list(nightDat$ambientRnd), mean)
	colony <- j
	treatment <- treatList[j]
	dayTime <- 0
	day <- i
	
x <- tmp$ambient - sp
	y <- tmp[,ai] - sp
	airModel <- lm(y~x + 0)
	
	x <- tmp$ambient - sp
	y <- tmp[,bi] - sp
	broodModel <- lm(y~x + 0)
	
	tmp <- thermPerfData[1,]
	tmp$colony <- colony
	tmp$treatment <- treatment
	tmp[,1:5] <- c(colony, treatment, dayTime, day, cohort)
	tmp$airTPerf <- summary(airModel)$coefficients[1]
	tmp$broodTPerf <- summary(broodModel)$coefficients[1]

	 thermPerfData <- rbind(thermPerfData, tmp)
	}
	
	}
# tmp <- aggregate(dayDat, list(dayDat$ambientRnd), mean)
# thermPerfData <- rbind(thermPerfData, tmp)
	
}



thermPerfData <- thermPerfData[complete.cases(thermPerfData),]










### Plotting and run models
cls <- c('darkgreen', 'darkred', 'green', 'red')
#Invert thermPerfs (to 1-slope)
thermPerfData[,c('airTPerf', 'broodTPerf')] <-  1-thermPerfData[,c('airTPerf', 'broodTPerf')]
thermPerfData$uniqueCol <- as.factor(paste(thermPerfData$colony, thermPerfData$cohort))

thermPerfData$day <- as.numeric(thermPerfData$day)
subdata <- subset(thermPerfData, day > 0 & day < 5)

#pdf('/Users/james/Dropbox/Work/Neonicotinoids/thermoregulationExpts/figures/thermPerf.pdf')
boxplot(broodTPerf~treatment*dayTime, data = subdata, col = c('green', 'red'), xlab = "Treatment", ylab = "Brood thermoregulation performance", boxwex = 0.7, at = c(1,2,4,5), outline = FALSE, whisklty = 1, axes = FALSE, boxlwd = 0.1, whisklwd = 0.2)
axis(2)
#dev.off()

library(lme4)
library(lmerTest)
model <- lmer(broodTPerf~treatment*dayTime+(1|day) +(1|uniqueCol) + (1|cohort), data = subdata)
plot(model)

#Significant interaction between treatment and nightime, so build separate models for each time period
modelDay <- lmer(broodTPerf~treatment+(1|day) +(1|uniqueCol) + (1|cohort), data = subset(subdata, dayTime == 1))
modelNight <- lmer(broodTPerf~treatment+(1|day) +(1|uniqueCol) + (1|cohort), data = subset(subdata, dayTime == 0))






# #### Plotting stuff
# data <- read.csv('/Users/james/Dropbox/Work/Neonicotinoids/thermoregulationExpts/data/cohort2/summaryData.csv')

# subdata = subset(data, time > 1)
# times <- subdata$time
# xl <- c(10,35)
# for(i in seq(1,length(times),10)){
	
	# plot(control[i]~ambient[i], data = subdata, xlim = xl, ylim = xl, pch = 19, cex = 0.4, axes = FALSE, ann = FALSE)
	# if(i == 1){
		# axis(1)
		# axis(2)
	# }
# par(new = TRUE)
	
	
# }


# #png('/Users/james/Dropbox/Work/Neonicotinoids/thermoregulationExpts/figures/ControlTemps.png')
# par(mfcol = c(1,2))
# data <- read.csv('/Users/james/Dropbox/Work/Neonicotinoids/thermoregulationExpts/data/cohort2/summaryData.csv')
# subdata = subset(data, time > 1)
	# plot(control~ambient, data = subdata, xlim = xl, ylim = xl, pch = 19, cex = 0.3, axes = FALSE, col = rgb(1, 0.5, 0.02, 0.2),xlab = "Ambient Temperature (C)", ylab = "Sham Brood Temperature (C)")
	# axis(1)
	# axis(2)
	# abline (a = 0, b = 1, lty = 2, lwd = 2, col = 'grey50')
		# #lines(c1brood~ambient, data = subdata, col = rgb(0.5, 1, 0.5, 0.3))

	
# data <- read.csv('/Users/james/Dropbox/Work/Neonicotinoids/thermoregulationExpts/data/cohort1/summaryData.csv')
# subdata = subset(data, time > 1)
	# plot(control~ambient, data = subdata, xlim = xl, ylim = xl, pch = 19, cex = 0.3, axes = FALSE, col = rgb(1, 0.5, 0.02, 0.2),xlab = "Ambient Temperature (C)", ylab = "Sham Air Temperature (C)")
	# axis(1)
	# axis(2)
	# abline (a = 0, b = 1, lty = 2, lwd = 2, col = 'grey50')
	
	# #lines(c1air~ambient, data = subdata, col = rgb(0.5, 1, 0.5, 0.3))
# #dev.off()
