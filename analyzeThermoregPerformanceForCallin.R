# Callin Switzer
# 5 Oct 2017
# Looking over James' code for analysis of thermoregulation data


ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("ggplot2", "plyr", "viridis", "lme4", "tidyr", 
              "dplyr", "lubridate", "signal", "zoo", "effects", 
              "gamm4", "MuMIn")
ipak(packages)


# set ggplot theme
theme_set(theme_classic())

# set wd
setwd("/Users/cswitzer/Documents/GitRepos/thermoRegImidacloprid")

# Define temperature set point
sp <- 32.5 

#___________________________________________
# Load data from each cohort
#___________________________________________
# cohort 1
data <- read.csv('data/summaryDataC1.csv')
treatList1 <- c('c', 't', 'c', 't', 't', 'c')
data$cohort = 1


# cohort #2
data2 <- read.csv('data/summaryDataC2.csv')
treatList2 <- c('c', 't', 't', 't', 'c', 'c')
data2$cohort = 2

# cohort #3
data3 <- read.csv('data/summaryDataC3.csv')
treatList3 <- c('t', 'c', 'c', 'c', 't', 't')
data3$cohort = 3

#___________________________________________
# combine cohorts
#___________________________________________
combData <- rbind(data, data2, data3)

# converts time to a proportion of a day (0 is midnight, 0.5 is noon)
combData$time1 <- (combData$time) %% 1

# add day
combData$day = floor(combData$time)

# Make the "time" variable into a datetime format -- 
# dates are all set 1970-01-01
# 86400 = 24 * 60 * 60, which converts days to seconds
combData$time2 <- format(as.POSIXct((combData$time1) * 86400, 
                                    origin = "1970-01-01", tz = "UTC"), 
                         "%H:%M:%S")

# time3 is the hour as a number, for example 14.5 is equal to 14:30
combData$time3 <- as.numeric(substr(combData$time2, 1,2)) +
  as.numeric(substr(combData$time2, 4,5))/60 + 
  as.numeric(substr(combData$time2, 7,8)) / 60 / 60

#___________________________________________________________
# smooth ambient temperature and calculate increasing vs.
# decreasing
#___________________________________________________________

dataIncDec <- lapply(unique(combData$cohort), function(ii){
  foo <- combData[combData$cohort == ii, ]
  
  plot(foo$ambient, x=foo$time, pch = ".")
  foo$rollMeanAmbient <- rmn <- rollapply(foo$ambient, width = 482, FUN = mean, align = "center", partial = TRUE)
  lines(y = rmn, x = foo$time, col = 'orange', lwd = 1)
  
  # what is 482 points?  -- answer: 3 hours
  print (mean(sapply(1:100, function(jj){
    foo$time3[jj + 482] - foo$time3[jj]})))
  
  foo$tempDer <- c(NA, diff(rmn)) # add NA to the beginning, since diff reduces size by 1
  foo$tempIncrease <- ifelse(foo$tempDer >= 0, yes = "Increasing", no = "Decreasing")
  
  foo$tempIncrease_noSM <- c(NA, diff(foo$ambient))
  foo$tempIncrease_noSmooth <- ifelse(foo$tempIncrease_noSM >= 0, yes = "Increasing", no = "Decreasing")
  
  return(foo)
})

# merge back into dataframe
df3 <- rbind.fill(dataIncDec)
combData <- merge(combData, df3, all.x = TRUE)

# visualize increasing vs. decreasing
p1 <- ggplot(combData, aes(x = time, y = ambient)) + 
  geom_line(size = 0.1) + 
  geom_point(aes(y = rollMeanAmbient, color = tempIncrease), size = 0.01, alpha = 0.1) + 
  facet_grid(cohort~., labeller = labeller(.rows = label_both, .cols = label_both)) + 
  scale_color_viridis(discrete = TRUE,begin = 0.4, end = 1) + 
  theme(legend.position = "none") + 
  labs(x = "Time (days)", y = "Ambient Temperature (C)")
p1

png("figures/ambientIncreasingDecreasing.png", width = 10, height= 7, units = "in", res = 250)
p1
dev.off()

pdf("figures/ambientIncreasingDecreasing.pdf", width = 10, height= 7)
p1
dev.off()
  
#___________________________________________
# convert combined dataset to long format
# and add some extra variables
#___________________________________________
data_long <- gather(combData, condition, temp, c1brood:c6air, factor_key=TRUE)


# the colony variable includes cohort_colony information
data_long$colony = paste(data_long$cohort,  
                         substr(data_long$condition,start = 2,2 ), 
                         sep = "_")

# location is either air or brood
data_long$location = substr(data_long$condition,start = 3,1000 )

# dayTime is either day or night
data_long$dayTime = ifelse(data_long$time3> 6 & data_long$time3 < 20, "day", "night")

# insert treatment
data_long$treatment <- mapvalues(data_long$colony, 
                                 from = c('1_1' , '1_2' , '1_3' , '1_4' , '1_5' , 
                                          '1_6' , '2_1' , '2_2' , '2_3' , '2_4' , 
                                          '2_5' , '2_6' , '3_1' , '3_2' , '3_3' , 
                                          '3_4' , '3_5' , '3_6'), 
                                 to = c(treatList1, treatList2, treatList3))

data_long$treatment <- mapvalues(data_long$treatment, from = c("c", "t"), 
                                 to = c("control_grp", "treatment_grp"))


#drop NA's -- 36 is from 18 colonies and air/brood measurements in each colony
sum(is.na(data_long$tempIncrease)) 
data_long <- data_long[!is.na(data_long$tempIncrease), ]

#___________________________________________
# visualize brood temperature vs. 
# ambient temperature
#___________________________________________
# new dataset that includes only brood temp (not air)
brooddta <- data_long[data_long$location == "brood", ]
airdta <- data_long[data_long$location == "air", ]


#___________________________________________
# conduct very simple analysis
# is avg brood temp lower in treated colonies than non-treated colonies?
#___________________________________________

# calculate colony means for ambient temp - brood temp for all days
colonyMeans <- data.frame(tapply((brooddta$temp - brooddta$ambient), INDEX = brooddta$colony, mean))
colonyMeans$colony = rownames(colonyMeans)
colonyMeans$treatment <- mapvalues(colonyMeans$colony, 
                                   from = c('1_1' , '1_2' , '1_3' , '1_4' , '1_5' , 
                                            '1_6' , '2_1' , '2_2' , '2_3' , '2_4' , 
                                            '2_5' , '2_6' , '3_1' , '3_2' , '3_3' , 
                                            '3_4' , '3_5' , '3_6'), 
                                   to = c(treatList1, treatList2, treatList3))
colonyMeans$treatment2 <- mapvalues(colonyMeans$treatment, from = c("c", "t"), 
                                 to = c("No Imidacloprid", "Nectar + Imidacloprid"))
colonyMeans$treatment2 <- relevel(as.factor(colonyMeans$treatment2), ref = "No Imidacloprid")
colnames(colonyMeans)[1] <- "broodTemp"
cm1 <- colonyMeans[order(colonyMeans$broodTemp, decreasing = FALSE), ]

# do the same for air temperature
colonyMeans <- data.frame(tapply((airdta$temp - airdta$ambient), INDEX = airdta$colony, mean))
colonyMeans$colony = rownames(colonyMeans)
colnames(colonyMeans)[1] = "airTemp"

mdf <- merge(cm1, colonyMeans, by = c("colony"))

# permuation test function
permFunc <- function(treatment, temp){
  # get a different permutation of treatment
  t2 <- sample(treatment, replace = FALSE)
  diff(tapply(temp, t2, mean))
}

# conduct permutation test for air
perms <- sapply(1:100000, FUN = function(i) permFunc(mdf$treatment, mdf$airTemp))
actualDiff <- diff(tapply(mdf$airTemp, mdf$treatment, mean))

# p-value for randomization test for airTemp
mean(perms < actualDiff | perms > - actualDiff) # ~0.00822


# conduct permutation test for brood
perms <- sapply(1:100000, FUN = function(i) permFunc(mdf$treatment, mdf$broodTemp))
actualDiff <- diff(tapply(mdf$broodTemp, mdf$treatment, mean))

# p-value for randomization test for Brood Temp
mean(perms < actualDiff | perms > - actualDiff) # ~0.00491

#___________________________________________
# Model relationship between ambient temp and
# brood temp, without transforming the response variable
#___________________________________________
# make a new variable that represents "day" as an integer
# the dayInt variable is unique among cohorts
# cohort 1 starts with day 100, cohort 2 with 200, and 3 with 300
brooddta$dayInt <- as.numeric(brooddta$day) + as.numeric(brooddta$cohort) * 100

# 20K points is about the max I can do in a reasonable time.
set.seed(17444)
brooddta_sm <- sample_n(brooddta, 5000, replace = FALSE)
brooddta_sm$treatment <- relevel(as.factor(brooddta_sm$treatment), ref = "control_grp")

# make this into a factor
brooddta_sm$tempIncrease <- as.factor(brooddta_sm$tempIncrease)

# make interaction into a single factor
brooddta_sm$treatTempIncrInt <- interaction(brooddta_sm$treatment, brooddta_sm$tempIncrease)

g01 <- gamm4(temp ~ s(ambient, by = treatment, k = 5) + s(time1, by = treatment, k = 5) + treatTempIncrInt, random = ~ (1|colony) + (1|day), data = brooddta_sm, REML = FALSE)

# use AIC to compare models
# note: REML must be FALSE to get an accurate AIC for comparison
# note: difference between AIC and AICc is negligible when sample size is large
AICc(g01$mer) 


g02 <- gamm4(temp ~ s(ambient, by = treatTempIncrInt, k = 5) + s(time1, by = treatTempIncrInt, k = 5) + treatTempIncrInt, random = ~ (1|colony) + (1|day), data = brooddta_sm, REML = FALSE)
AICc(g02$mer) # this AICc is lower than above

# refit with REML=TRUE
g02 <- gamm4(temp ~ s(ambient, by = treatTempIncrInt, k = 5) + s(time1, by = treatTempIncrInt, k = 5) + treatTempIncrInt, random = ~ (1|colony) + (1|day), data = brooddta_sm, REML = TRUE)


pdf("figures/smoothplot.pdf", width = 15, height = 15)
par(mfrow = c(3,4))
aab <- plot(g02$gam, all.terms = TRUE, rug = FALSE, shade = TRUE)
summary(g02$gam) # Summary for paper -- brood
summary(g02$mer)
dev.off()


plot(g02$mer)

# residuals are spread a little wide.
mean(residuals(g02$mer) > 2 | residuals(g02$mer) < -2)


# plot actual vs. predicted
brooddta_sm$preds1 <-  predict(g02$mer, type = 'response')
brooddta_sm$preds <-  predict(g02$mer, type = 'response', re.form = NA)
ggplot(brooddta_sm, aes(x = time, y= preds1)) + 
  facet_wrap(~colony + treatment) + 
 geom_line(aes(y = temp), color = 'black') +  # black is actual
 geom_line(alpha = 0.5, color = 'red')  + # red is predicted
  geom_line(alpha = 0.5,aes(y = preds),  color = 'green') # green is predictions without random effects
  

## Overall, model is not obtimal, but it is more interpretable because
# we didn't have to transform the y-variable

# visualize differences
# plot raw data, and prediction, while holding time constant
nd <- brooddta_sm[, c("ambient", "treatment", "time1", "tempIncrease", "treatTempIncrInt")]
nd$time1 = 0.5
nd$ambient <- mean(nd$ambient)

nd <- nd[!(duplicated(nd)), ]
nd

nd2 <- expand.grid(ambient = c( 10, 20 , 30), 
                   treatTempIncrInt = levels(interaction(levels(nd$treatment), levels(nd$tempIncrease))), 
                   time1 = seq(0,1, length.out = 30))

nd2$preds1 <-  predict(g02$gam, newdata = nd2, type = 'response', re.form = NA)
nd2$se1 <- predict(g02$gam, newdata = nd2, type = 'response', re.form = NA, se = TRUE)$se

pdf('figures/facetedPrededSmooths.pdf', width = 10, height = 6)
ggplot(nd2, aes(x = time1,  y= preds1, color = as.factor(ambient))) + 
  geom_ribbon(aes(ymin = preds1 - 1.96*se1, ymax = preds1 + 1.96*se1), alpha = 0.4, color = NA) +
  facet_grid(ambient~treatTempIncrInt, labeller = labeller(.rows = label_both, .cols = label_value)) +
  geom_line(aes(y = preds1)) + 
  #geom_errorbar(aes(ymin = preds1 - 1.96*se1, ymax = preds1 + 1.96*se1), width = 0.05, position = position_dodge(width = 0.02))+ 
  scale_color_viridis(discrete = TRUE, option = "C", end = 0.7, name = "Ambient\nTemp (C)") +
  labs(x = "Time of day (0=midnight, 0.5 = noon)", y="Predicted brood temperature (C)")
dev.off()


# save increasing vs. decreasing plots
pdf('figures/gamSmoothOnIncreasingDecreasing.pdf', width = 10, height = 7)
ggplot(sample_n(brooddta, 50000, replace = FALSE), aes(x = ambient,  y= temp, color = interaction( treatment))) + 
  facet_grid( ~ tempIncrease, labeller = labeller(.rows = label_both, .cols = label_both)) +
  geom_point(aes(y = temp), alpha = 0.05) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE) + 
  scale_color_viridis(discrete = TRUE, option = "C", end = 0.7, name = "Treatment group") +
  labs(x = "Ambient temp (C)", y="Brood temperature (C)")
dev.off()



#___________________________________________
# Model relationship between ambient temp and
# AIR temp, 
#___________________________________________
# make a new variable that represents "day" as an integer
# the dayInt variable is unique among cohorts
# cohort 1 starts with day 100, cohort 2 with 200, and 3 with 300
airdta$dayInt <- as.numeric(airdta$day) + as.numeric(airdta$cohort) * 100

# 20K points is about the max I can do in a reasonable time.
set.seed(17444)
airdta_sm <- sample_n(airdta, 5000, replace = FALSE)
airdta_sm$treatment <- relevel(as.factor(airdta_sm$treatment), ref = "control_grp")

# make this into a factor
airdta_sm$tempIncrease <- as.factor(airdta_sm$tempIncrease)

# make interaction into a single factor
airdta_sm$treatTempIncrInt <- interaction(airdta_sm$treatment, airdta_sm$tempIncrease)

g01 <- gamm4(temp ~ s(ambient, by = treatment, k = 5) + s(time1, by = treatment, k = 5) + treatTempIncrInt, random = ~ (1|colony) + (1|day), data = airdta_sm, REML = FALSE)

# use AIC to compare models
# note: REML must be FALSE to get an accurate AIC for comparison
# note: difference between AIC and AICc is negligible when sample size is large
AICc(g01$mer) 


g02 <- gamm4(temp ~ s(ambient, by = treatTempIncrInt, k = 5) + s(time1, by = treatTempIncrInt, k = 5) + treatTempIncrInt, random = ~ (1|colony) + (1|day), data = airdta_sm, REML = FALSE)
AICc(g02$mer) # this AICc is lower than above

# refit with REML=TRUE
g02 <- gamm4(temp ~ s(ambient, by = treatTempIncrInt, k = 5) + s(time1, by = treatTempIncrInt, k = 5) + treatTempIncrInt, random = ~ (1|colony) + (1|day), data = airdta_sm, REML = TRUE)


pdf("figures/smoothplot_air.pdf", width = 15, height = 15)
par(mfrow = c(3,4))
aab <- plot(g02$gam, all.terms = TRUE, rug = FALSE, shade = TRUE)
summary(g02$gam) # Summary for paper -- air temp
summary(g02$mer)
dev.off()


plot(g02$mer)

# residuals are fine and dandy here!
mean(residuals(g02$mer) > 2 | residuals(g02$mer) < -2)


# plot actual vs. predicted
airdta_sm$preds1 <-  predict(g02$mer, type = 'response')
airdta_sm$preds <-  predict(g02$mer, type = 'response', re.form = NA)
ggplot(airdta_sm, aes(x = time, y= preds1)) + 
  facet_wrap(~colony + treatment) + 
  geom_line(aes(y = temp), color = 'black') +  # black is actual
  geom_line(alpha = 0.5, color = 'red')  + # red is predicted
  geom_line(alpha = 0.5,aes(y = preds),  color = 'green') # green is predictions without random effects


## Overall, this model is easier to fit than brood temp

# visualize differences
# plot raw data, and prediction, while holding time constant
nd <- airdta_sm[, c("ambient", "treatment", "time1", "tempIncrease", "treatTempIncrInt")]
nd$time1 = 0.5
nd$ambient <- mean(nd$ambient)

nd <- nd[!(duplicated(nd)), ]
nd

nd2 <- expand.grid(ambient = c( 10, 20 , 30), 
                   treatTempIncrInt = levels(interaction(levels(nd$treatment), levels(nd$tempIncrease))), 
                   time1 = seq(0,1, length.out = 30))

nd2$preds1 <-  predict(g02$gam, newdata = nd2, type = 'response', re.form = NA)
nd2$se1 <- predict(g02$gam, newdata = nd2, type = 'response', re.form = NA, se = TRUE)$se

pdf('figures/facetedPrededSmooths_Air.pdf', width = 10, height = 6)
ggplot(nd2, aes(x = time1,  y= preds1, color = as.factor(ambient))) + 
  geom_ribbon(aes(ymin = preds1 - 1.96*se1, ymax = preds1 + 1.96*se1), alpha = 0.4, color = NA) +
  facet_grid(ambient~treatTempIncrInt, labeller = labeller(.rows = label_both, .cols = label_value)) +
  geom_line(aes(y = preds1)) + 
  #geom_errorbar(aes(ymin = preds1 - 1.96*se1, ymax = preds1 + 1.96*se1), width = 0.05, position = position_dodge(width = 0.02))+ 
  scale_color_viridis(discrete = TRUE, option = "C", end = 0.7, name = "Ambient\nTemp (C)") +
  labs(x = "Time of day (0=midnight, 0.5 = noon)", y="Predicted brood temperature (C)")
dev.off()


# save increasing vs. decreasing plots
pdf('figures/gamSmoothOnIncreasingDecreasing_Air.pdf', width = 10, height = 7)
ggplot(sample_n(airdta, 50000, replace = FALSE), aes(x = ambient,  y= temp, color = interaction( treatment))) + 
  facet_grid( ~ tempIncrease, labeller = labeller(.rows = label_both, .cols = label_both)) +
  geom_point(aes(y = temp), alpha = 0.05) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE) + 
  scale_color_viridis(discrete = TRUE, option = "C", end = 0.7, name = "Treatment group") +
  labs(x = "Ambient temp (C)", y="Brood temperature (C)")
dev.off()


