# Callin Switzer
# 14 Sept 2017
# Checking over James' code for analysis of thermoregulation data


library(ggplot2)
library(plyr)
library(viridis)

# set ggplot theme
theme_set(theme_bw())



#Manually navigate to working directory
setwd("data")

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
treatList1 <- c('c', 't', 'c', 't', 't', 'c')



data$time <- data$time + 6/24
# now, does time == 0.0 mean midnight, or 6am?
data$day <- floor(data$time)
days <- unique(data$day)
data$hour <- data$time - data$day #Rewrite data$hour

#Calculate rounded brood temperature
int <- 0.5
data$ambientRnd <- round((data$ambient/int))*int
daytimeRange <- c(0, 16/24) #Slightly tricky here - this is now isolating the hours from 6 am to 8 pm, but adjusted relatively to new timing (with 0 being 6 am on day one)
data$cohort = 1


table(data$day)
# why are you excluding day 0 and 18? (in the loop for 1:17)


data2 <- read.csv('summaryDataC2.csv')
cohort <- 2
treatList2 <- c('c', 't', 't', 't', 'c', 'c')
#Adjust time to count from 6 in the morning
data2$time <- data2$time + 6/24
data2$day <- floor(data2$time)
days <- unique(data2$day)
data2$hour <- data2$time - data2$day
#Calculate rounded brood temperature
data2$ambientRnd <- round((data2$ambient/int))*int
daytimeRange <- c(0, 16/24) #Slightly tricky here - this is now isolating the hours from 6 am to 8 pm, but adjusted relatively to new timing (with 0 being 6 am on day one)
data2$cohort = 2
head(data2)


data3 <- read.csv('summaryDataC3.csv')
cohort <- 3
treatList3 <- c('t', 'c', 'c', 'c', 't', 't')
#Adjust time to count from 6 in the morning
data3$time <- data3$time + 6/24
data3$day <- floor(data3$time)
days <- unique(data3$day)
data3$hour <- data3$time - data3$day
#Calculate rounded brood temperature
data3$ambientRnd <- round((data3$ambient/int))*int
daytimeRange <- c(0, 16/24) #Slightly tricky here - this is now isolating the hours from 6 am to 8 pm, but adjusted relatively to new timing (with 0 being 6 am on day one)
data3$cohort = 3
head(data3)


combData <- rbind(data, data2, data3)

head(combData)

# make new time with random dates
combData$time2 <- as.POSIXct(Sys.Date() + combData$time + combData$cohort * 20)
indxs <- seq(from = 1, to = nrow(combData), length.out = 1000)
plot(combData$c1air[indxs], x = combData$time2[indxs], type = 'l')


# now convert to long format
head(combData)
library(tidyr)


data_long <- gather(combData, condition, temp, c1brood:c6air, factor_key=TRUE)
head(data_long)


data_long$colony = paste(data_long$cohort,  substr(data_long$condition,start = 2,2 ), sep = "_")
data_long$location = substr(data_long$condition,start = 3,1000 )

#paste(sort(unique(data_long$colony)), collapse = "\' , \'")

data_long$treatment <- mapvalues(data_long$colony, from = c('1_1' , '1_2' , '1_3' , '1_4' , '1_5' , '1_6' , '2_1' , '2_2' , '2_3' , '2_4' , '2_5' , '2_6' , '3_1' , '3_2' , '3_3' , '3_4' , '3_5' , '3_6'), to = c(treatList1, treatList2, treatList3))

indxs2 <- seq(from = 1, to = nrow(data_long), length.out = 5000)

ggplot(data_long[indxs2, ], aes(x = time2, y = temp, color = colony)) + 
  geom_line() + 
  facet_wrap(~location)

ggplot(data_long[indxs2, ], aes(x = ambient, y = temp, color = treatment)) + 
  geom_point() + 
  facet_wrap(~location) + 
  geom_hline(aes(yintercept = sp)) + 
  scale_color_viridis(discrete = TRUE) + 
  geom_smooth()






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
	
	
	# plot(dayDat$c1brood, ylim = c(10, 40))
	# points(dayDat$control, col = 'red')
	# points(dayDat$c1air)
	# points(dayDat$ambient)
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
	
	
	x <- tmp$ambient 
	y <- tmp[,ai]
	airModel <- lm(y~x + 0)
	
	plot(x, y); abline(airModel)
	
	x <- tmp$ambient - sp
	y <- tmp[,bi] - sp
	broodModel <- lm(y~x + 0)
	
	plot(x, y); abline(broodModel)
	
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



head(thermPerfData)


ggplot(thermPerfData, aes(x = airTPerf, y = broodTPerf, color = treatment)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~dayTime)


plot(thermPerfData$airTPerf, thermPerfData$broodTPerf, col = as.factor(thermPerfData$treatment))






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
summary(model)
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
