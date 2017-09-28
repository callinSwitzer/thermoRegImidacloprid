# Callin Switzer
# 14 Sept 2017
# Looking over James' code for analysis of thermoregulation data



ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("ggplot2", "plyr", "viridis", "lme4", "tidyr", 
              "dplyr", "lubridate", "signal", "zoo")
ipak(packages)


# set ggplot theme
theme_set(theme_classic())

#set working directory
setwd("data")

# Define temperature set point
sp <- 32.5 

#___________________________________________
# Load data from each cohort
#___________________________________________
# cohort 1
data <- read.csv('summaryDataC1.csv')
treatList1 <- c('c', 't', 'c', 't', 't', 'c')
data$cohort = 1


# cohort #2
data2 <- read.csv('summaryDataC2.csv')
treatList2 <- c('c', 't', 't', 't', 'c', 'c')
data2$cohort = 2

# cohort #3
data3 <- read.csv('summaryDataC3.csv')
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

# create a variable for temp rising or falling
# smooth ambient temperature
dataIncDec <- lapply(unique(combData$cohort), function(ii){
  foo <- combData[combData$cohort == ii, ]
  plot(foo$ambient, x=foo$time, pch = ".")
  
  # 482 points is about 3 hours
  foo$rollMeanAmbient <- rmn <- rollapply(foo$ambient, width = 482, FUN = mean, align = "center", partial = TRUE)
  lines(y = rmn, x = foo$time, col = 'green', lwd = 1)
  
  foo$tempDer <- c(NA, diff(rmn)) # add NA to the beginning, since diff reduces size by 1
  foo$tempIncrease <- ifelse(foo$tempDer >= 0, yes = "Increasing", no = "Decreasing")
  
  foo$tempIncrease_noSM <- c(NA, diff(foo$ambient))
  foo$tempIncrease_noSmooth <- ifelse(foo$tempIncrease_noSM >= 0, yes = "Increasing", no = "Decreasing")
  
  return(foo)
})




# merge back into dataframe
df3 <- rbind.fill(dataIncDec)

combData <- merge(combData, df3, all.x = TRUE)

head(combData)

ggplot(combData, aes(x = time, y = ambient)) + 
  geom_line(size = 0.1) + 
  geom_point(aes(y = rollMeanAmbient, color = tempIncrease), size = 0.01, alpha = 0.1) + 
  facet_grid(cohort~., labeller = labeller(.rows = label_both, .cols = label_both)) + 
  scale_color_viridis(discrete = TRUE,begin = 0.4, end = 1) + 
  theme(legend.position = "none")
  
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

# new dataset that includes only brood temp (not air)
brooddta <- data_long[data_long$location == "brood", ]

#___________________________________________
# Analysis to visualize the distance less than 32.5
# How good are bees at maintaining temperature at 32.5?
#___________________________________________


# make a new variable that represents "day" as an integer
# the dayInt variable is unique among cohorts
# cohort 1 starts with day 100, cohort 2 with 200, and 3 with 300
brooddta$dayInt <- as.numeric(brooddta$day) + as.numeric(brooddta$cohort) * 100


# create new variables that are distance from 32.5
brooddta$brood_dist_from_32<- (sp - brooddta$temp)
brooddta$amb_dist_from_32 <- (sp - brooddta$ambient)



# subset data
set.seed(28383)
brooddta_sm <- sample_n(brooddta, 5000, replace = FALSE)
brooddta_sm$treatment <- relevel(as.factor(brooddta_sm$treatment), ref = "control_grp")

# make this into a factor
brooddta_sm$tempIncrease <- as.factor(brooddta_sm$tempIncrease)


bb <- ggplot(brooddta_sm, aes(x = ambient, y = temp)) + 
  geom_point(aes(color = treatment))



# visualize where the origin is going to go
bb + geom_hline(aes(yintercept = sp)) + geom_vline(aes(xintercept = sp)) + 
  labs(y = "Brood temp (C)", x = "ambient temp (C)") + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7)


# visualize reorientation of data
cc <- ggplot(brooddta_sm, aes(x = amb_dist_from_32, y = brood_dist_from_32, color = treatment)) + 
  geom_point(alpha = 0.5) + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7) + 
  xlab("Number of degrees below from 32.5 C (Ambient)") + 
  ylab("Number of degrees below 32.5 C (Brood)") + 
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) 
cc


dd <- ggplot(brooddta_sm, aes(x = amb_dist_from_32, y = abs(brood_dist_from_32), color = treatment)) + 
  geom_point(alpha = 0.5) + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7) + 
  xlab("Number of degrees below 32.5 C (Ambient)") + 
  ylab("Absolute number of degrees from 32.5 C (Brood)") + 
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0))
dd


#___________________________________________
# Fit a Generalized Linear Mixed Effects Model
# How good are bees at maintaining temperature at 32.5?
# Holding the temperature constant, does is night time associated
# with bees not being able to keep brood temperature near 32.5
#___________________________________________

brooddta$brood_abs_dist_from_32 <- abs(brooddta$brood_dist_from_32) + 0.1 
# the 0.1 makes the dataset have no zeros (requirement for gamma glm)
brooddta$amb_abs_dist_from_32 <- abs(brooddta$amb_dist_from_32)

# get scaled amb_abs_dist_from_32, which will be used for modeling
brooddta$scaled_amb_abs_dist_from_32 <- scale(brooddta$amb_abs_dist_from_32, center = TRUE, scale = TRUE)

# small dataset for temps below setpoint
set.seed(193984)
bsmall <- sample_n(brooddta[brooddta$ambient < sp, ], 10000, replace = FALSE)
head(bsmall)

bsmall$dayInt <- as.factor(bsmall$dayInt)

gm1 <- glmer(brood_abs_dist_from_32 ~ scaled_amb_abs_dist_from_32  * treatment + scaled_amb_abs_dist_from_32*  
               tempIncrease + scale(time3) +  
               I((scale(time3))^2) + I(scale(time3)^3) + (1|dayInt) + 
                (1|colony), data = bsmall, family = Gamma("sqrt"), 
             control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(gm1)

plot(residuals(gm1, type = "pearson"), x = bsmall$time3)

ndf <- data.frame(resds = residuals(gm1, type = 'pearson'), fitd = predict(gm1, type = 'response'))
ggplot(ndf, aes(y = resds, x = fitd)) + 
  geom_point() + 
  geom_smooth()

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

vif.mer(gm1) # doesn't seem too high, when there are no interactions

bsmall$preds <- predict(gm1, type = 'response', re.form=NA)
predDF<- data.frame(x = bsmall$scaled_amb_abs_dist_from_32, y = bsmall$preds, treatment = bsmall$treatment, colony = bsmall$colony, dayTime = bsmall$dayTime)
predDF <- predDF[order(predDF$x, predDF$colony), ]

centt <- attr(brooddta$scaled_amb_abs_dist_from_32, which = "scaled:center")
scle <- attr(brooddta$scaled_amb_abs_dist_from_32, which = "scaled:scale")


ee <- ggplot(bsmall, aes(x = abs(amb_dist_from_32), y = brood_abs_dist_from_32, color = treatment)) + 
  geom_point(alpha = 0.2) + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7) + 
  xlab("Absolute number of degrees from 32.5 C (Ambient)") + 
  ylab("Absolute number of degrees from 32.5 C (Brood)") + 
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0))
ee

ee + geom_line(aes(x = amb_abs_dist_from_32, y= preds, linetype = interaction(treatment, tempIncrease)))


ee +  geom_smooth(aes(linetype =  tempIncrease, color = treatment), size = 1)

# plot actual vs. predicted
bsmall$preds1 <- predict(gm1, type = 'response')
ggplot(bsmall, aes(x = time, y= preds1)) + 
  facet_wrap(~colony + treatment) + 
  geom_line(aes(y = brood_abs_dist_from_32), color = 'black') +  # black is actual
  geom_line(alpha = 0.5, color = 'red') # red is predicted


dev.off()




library(gamm4)
head(bsmall)
bsmall$treatment <- as.factor(bsmall$treatment)
g1 <- gamm4(brood_abs_dist_from_32 ~ s(scaled_amb_abs_dist_from_32, by = treatment, k = 5) + 
              treatment + s(time1, k = 10) + treatment*tempIncrease, 
            random =  ~(1|dayInt) +  (1|colony), 
            data = bsmall,family = Gamma('log'), REML = TRUE)
summary(g1$gam)
summary(g1$mer)





par(mfrow = c(2,3))
plot(g1$gam, all.terms = TRUE, rug = FALSE)
plot(g1$mer)

bsmall$residsHigh <- residuals(g1$mer, type= "pearson") > 2

# residuals are spread a little wide.
mean(residuals(g1$mer) > 2 | residuals(g1$mer) < -2)

dev.off()
plot(y = residuals(g1$mer, type= "pearson"), x = predict(g1$mer, type = 'link'))

bsmall[bsmall$residsHigh, ]




# plot actual vs. predicted
bsmall$preds1 <-  predict(g1$mer, type = 'response')
ggplot(bsmall, aes(x = time, y= preds1)) + 
  facet_wrap(~colony + treatment) + 
  geom_line(aes(y = brood_abs_dist_from_32), color = 'black') +  # black is actual
  geom_line(alpha = 0.5, color = 'red') # red is predicted



gg <- ggplot(bsmall, aes(x = abs(amb_dist_from_32), y = brood_abs_dist_from_32, color = residsHigh)) + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~residsHigh) + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7) + 
  xlab("Absolute number of degrees from 32.5 C (Ambient)") + 
  ylab("Absolute number of degrees from 32.5 C (Brood)") + 
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0))
gg






