# Callin Switzer
# 14 Sept 2017
# Looking over James' code for analysis of thermoregulation data


library(ggplot2)
library(plyr)
library(viridis)
library(lme4)
library(tidyr)
library(dplyr)
library(lubridate)
library(timeSeries)

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
  
#___________________________________________
# convert combined dataset to long format
# and add some extra variables
#___________________________________________
data_long <- gather(combData, condition, temp, c1brood:c6air, factor_key=TRUE)
head(data_long)

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

#___________________________________________
# double check to make sure the hours seem right
#___________________________________________
ggplot(sample_n(data_long[data_long$location == 'air', ], size = 1000), 
       aes(x = dayTime, y = temp)) +
  geom_boxplot()

# visualize times vs. temps, to see if the times are aligned
ggplot(sample_n(data_long[data_long$location == 'air', ], size = 1000), 
       aes(x = time3, y = temp, color = as.factor(cohort))) +
  geom_point() + 
  geom_smooth(aes(group = cohort)) + 
  labs(x = "Hour of the day", y = 'air temp')


#___________________________________________
# Visualize a few different smoothing options
#___________________________________________
# generate a subsample of data to speed up visualization
indxs2 <- sample(1:nrow(data_long), size  = 10000)

aa <- ggplot(data_long[indxs2, ], aes(x = ambient, y = temp, color = treatment)) + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~location) + 
  geom_hline(aes(yintercept = sp)) + 
  scale_color_viridis(discrete = TRUE) 

# polynomial smooth
aa + 
  stat_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1)

# log
aa + 
  stat_smooth(method = "lm", formula = y ~ log(x), size = 1)

# generalized additive model smooth
aa + 
  stat_smooth(method = "gam", formula = y ~ s(x, k = 5), size = 1)

# loess
aa + 
  stat_smooth(method = "loess", size = 2)

#___________________________________________
# visualize brood temperature vs. 
# ambient temperature
#___________________________________________
# new dataset that includes only brood temp (not air)
brooddta <- data_long[data_long$location == "brood", ]
nrow(brooddta)


# create a variable for temp rising or falling
# smooth ambient temperature
set.seed(123)
sm <- sample_n(brooddta, 100000)
sm <- sm[order(sm$cohort, sm$time), ]

dataIncDec <- lapply(unique(brooddta$cohort), function(ii){
  foo <- sm[sm$cohort == ii, ]
  plot(foo$ambient, x=foo$time)
  
  ambTime <- timeSeries(foo$ambient)
  
  foo$ss <- as.data.frame(smoothSpline(ambTime,spar = 0.3))$spline
  
  
  foo$tempDer <- c(NA, diff(foo$ss)) # add NA to the beginning, since diff reduces size by 1
  foo$tempIncrease <- ifelse(foo$tempDer >= 0, yes = "Increasing", no = "Decreasing")
  lines(x = foo$time, y = foo$ss, col = 'green')
  
  return(foo)
})


# ggplot(foo, aes(x = time, y = ambient, color = tempIncrease)) + 
#   geom_point(size = 0.1) + 
#   geom_line(aes(y = ss), color = 'black', alpha = 0.5) + 
#   scale_color_viridis(discrete = TRUE)


# merge back into dataframe
df3 <- rbind.fill(dataIncDec)

broodsm <- merge(sm, df3)

head(broodsm)

ggplot(broodsm, aes(x = time, y = ambient, color = tempIncrease)) + 
  geom_point(size = 0.1) + 
  facet_grid(cohort~., labeller = labeller(.rows = label_both, .cols = label_both)) + 
  scale_color_viridis(discrete = TRUE)



# # visualize brood temp, faceted by treatment, with hexbin plot
# ggplot(brooddta, aes(x = ambient, y = temp)) +
#   geom_hex() +
#   facet_grid(~treatment) +
#   geom_hline(aes(yintercept = sp), linetype = 2) +
#   scale_fill_viridis(option = "A", direction = -1) + 
#   stat_smooth(data = brooddta_sm, method = "loess", se= FALSE, color = 'black') + 
#   labs(x =  expression("Ambient Temperature " ( degree*C)), 
#        y = expression("Brood Temperature " ( degree*C))) + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.border = element_rect(colour = "black"))

# # visualize density on the same plot
# set.seed(123)
# brooddta_sm <- sample_n(brooddta, 10000, replace = FALSE)
# brooddta_sm$treatment <- relevel(as.factor(brooddta_sm$treatment), ref = "control_grp")
# 
# # don't use ggsave(), because those pdf's aren't editable!
# pdf("~/Desktop/gamSmooth3.pdf", width = 7, height = 5)
# ggplot(brooddta_sm, aes(x = ambient, y = temp)) +
#   geom_point(alpha = 0.2, aes(color = treatment), shape = 20, stroke = 0, size = 1) + 
#   geom_hline(aes(yintercept = sp), linetype = 2) +
#   scale_fill_viridis(option = "C", discrete = TRUE, end = 0.7)  +
#   scale_color_viridis(option = "C", discrete = TRUE, end = 0.7)  +
#   stat_smooth(method = "gam", formula = y ~ s(x, k = 4), data = brooddta, 
#               se= FALSE, aes(color = treatment), size = 1) + 
#   labs(x =  expression("Ambient Temperature " ( degree*C)), 
#        y = expression("Brood Temperature " ( degree*C))) + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         panel.border = element_rect(colour = "black"))
# dev.off()


#___________________________________________
# Model relationship between ambient temp and
# brood temp, without transforming the response variable
#___________________________________________


library(gamm4)
# make a new variable that represents "day" as an integer
# the dayInt variable is unique among cohorts
# cohort 1 starts with day 100, cohort 2 with 200, and 3 with 300
broodsm$dayInt <- as.numeric(broodsm$day) + as.numeric(broodsm$cohort) * 100

# drop NA's -- there will be three from the numerical derivative above
broodsm <- broodsm[!is.na(broodsm$tempIncrease), ]



# 20K points is about the max I can do in a reasonable time.
set.seed(456)
brooddta_sm <- sample_n(broodsm, 5000, replace = FALSE)
brooddta_sm$treatment <- relevel(as.factor(brooddta_sm$treatment), ref = "control_grp")

# make this into a factor
brooddta_sm$tempIncrease <- as.factor(brooddta_sm$tempIncrease)
#g1 <- gamm4(temp ~ s(ambient, k = 10) + treatment +  + s(time1, k = 5), random = ~ (1|colony) + (1|dayInt), data = brooddta_sm)

g1 <- gamm4(temp ~ s(ambient, by = treatment, k = 5) + s(time1, k = 5) + treatment * tempIncrease, random = ~ (1|colony) + (1|dayInt), data = brooddta_sm, REML = FALSE)

deviance(g1$mer)

g2 <- gamm4(temp ~ s(ambient, by = treatment, k = 5) + s(time1, k = 5) + treatment + tempIncrease, random = ~ (1|colony) + (1|dayInt), data = brooddta_sm, REML = FALSE)

deviance(g2$mer)


anova(g2$mer, g1$mer)


g3 <- gamm4(temp ~ s(ambient, by = treatment, k = 5) + s(time1, k = 5) + treatment * tempIncrease + cohort, random = ~ (1|colony) + (1|dayInt), data = brooddta_sm, REML = FALSE)
summary(g3$mer)

anova(g1$mer, g3$mer) # drop cohort


g1 <- gamm4(temp ~ s(ambient, by = treatment, k = 5) + s(time1, k = 5) + treatment * tempIncrease, random = ~ (1|colony) + (1|dayInt), data = brooddta_sm, REML = TRUE)
summary(g1$gam)
summary(g1$mer)


par(mfrow = c(2,3))
plot(g1$gam, all.terms = TRUE, rug = FALSE)
plot(g1$mer)

# residuals are spread a little wide.
mean(residuals(g1$mer) > 2 | residuals(g1$mer) < -2)


# plot actual vs. predicted
brooddta_sm$preds1 <-  predict(g1$mer, type = 'response')
ggplot(brooddta_sm, aes(x = time, y= preds1)) + 
  facet_wrap(~colony + treatment) + 
  geom_line(aes(y = temp), color = 'black') +  # black is actual
  geom_line(alpha = 0.5, color = 'red') # red is predicted


## Overall, model is not obtimal,but we didn't have to transform the y-variable


# plot raw data, and prediction, while holding time constant
nd <- brooddta_sm
#nd$time1 = 0
#nd$ambient <- mean(nd$ambient)

tapply(brooddta_sm$temp, INDEX = brooddta_sm$treatment, mean)

brooddta_sm$preds1 <-  predict(g1$gam, newdata = nd, type = 'response', re.form = NA)
#brooddta_sm$preds1 <-  predict(g1$mer, type = 'response', re.form = NA)

ggplot(brooddta_sm, 
       aes(x = ambient, y= temp, color = treatment)) + 
  geom_point(alpha = 0.1) + 
  geom_line(aes(y = preds1)) + facet_wrap(~treatment)

unique(brooddta_sm$preds1)


# could calculate out of sample accuracy to show significance....just an idea.

#________________________________________________________________
### end of gamm4
#________________________________________________________________


# confidence bands, on a daily level

ggplot(sample_n(brooddta, 50000), aes(x = time3, y = temp)) + 
  geom_point(aes(color = treatment), alpha = 0.01) + 
  facet_grid(cohort~.) + 
  stat_smooth(aes(color = treatment), method = "gam", formula = y ~ s(x ,k = 5))




# log-transformed model
expMod <- lmer(temp ~ log(ambient) * treatment + scale(time3) + I(scale(time3)^2) + I(scale(time3)^3) + dayTime +  (1|colony) + (1|dayInt), data = brooddta)

expMod <- lmer( temp ~ log(ambient, base = 40) + treatment +  (1|colony) + (1|dayInt), data = brooddta)

expMod <- lmer( temp ~ log(ambient) + I(log(ambient)^2) + I(log(ambient)^3) + treatment +  (1|colony) + (1|dayInt), data = brooddta)


car::vif(lm(temp ~ ambient) + treatment + scale(time3) + I(scale(time3)^2) + I(scale(time3)^3),data = brooddta))

summary(expMod)
plot(expMod) # somewhat troubling -- residuals are very high, and fan shaped

bb <- ggplot(brooddta, aes(x = ambient, y = temp, color = treatment)) +
  geom_point(alpha = 0.5) + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7)

bb



# predict and visualize from log-transormed model
brooddta$expModPreds <- predict(expMod, re.form = NA)
# note: re.form means ignore random effects, and predict an overall average

bbb <- ggplot(brooddta,  aes(x = ambient, y = temp, color = treatment)) + 
  geom_point(alpha = 0.1) + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7) + 
  xlab("Temp (Ambient)") + 
  ylab("Temp (Brood)") + 
  geom_hline(aes(yintercept = sp)) + geom_vline(aes(xintercept = sp))
bbb +  geom_line(aes( y=expModPreds, linetype = dayTime), size = 1)

bbb +  stat_smooth(method = "loess", se = F,  span = 0.99, aes(linetype = dayTime), size = 2)



bb + facet_wrap(~treatment + dayTime) + geom_line(data = brooddta, aes(y = expModPreds), col = 'black')

bbb + geom_line(data = brooddta, aes(y = expModPreds, linetype =  dayTime), size = 1)

bb + geom_line(data = brooddta, aes(y = expModPreds, color = interaction(treatment, dayTime)), size = 3) + 
  facet_wrap(~dayTime)


#___________________________________________
# Analysis to visualize the distance less than 32.5
# How good are bees at maintaining temperature at 32.5?
#___________________________________________
brooddta$brood_dist_from_32<- (sp - brooddta$temp)

brooddta$amb_dist_from_32 <- (sp - brooddta$ambient)


# visualize where the origin is going to go
bb + geom_hline(aes(yintercept = sp)) + geom_vline(aes(xintercept = sp)) + 
  labs(y = "Brood temp (C)", x = "ambient temp (C)")


# visualize reorientation of data
cc <- ggplot(brooddta, aes(x = amb_dist_from_32, y = brood_dist_from_32, color = treatment)) + 
  geom_point(alpha = 0.5) + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7) + 
  xlab("Number of degrees below from 32.5 C (Ambient)") + 
  ylab("Number of degrees below 32.5 C (Brood)") + 
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0)) 
cc


dd <- ggplot(brooddta, aes(x = amb_dist_from_32, y = abs(brood_dist_from_32), color = treatment)) + 
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

bsmall <- brooddta #sample_n(brooddta, 10000)
head(bsmall)

bsmall$dayInt <- as.factor(bsmall$dayInt)

gm1 <- glmer(brood_abs_dist_from_32 ~ scaled_amb_abs_dist_from_32  + treatment + dayTime + scale(time3) +  
               I((scale(time3))^2) + I(scale(time3)^3) + 
               (1|dayInt) + (1|colony), data = bsmall, family = Gamma("sqrt"))


# gm11 <- lm(brood_abs_dist_from_32 ~ time3 +  I((time3)^2) +I(time3^3),  data = bsmall)
# summary(gm11)
# plot(brood_abs_dist_from_32 ~  time3,  data = bsmall)
# points(x = bsmall$time3, y = predict(gm11), col = 'red' )



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
  geom_point(alpha = 0.02) + 
  scale_color_viridis(option = "C", discrete = TRUE, end = 0.7) + 
  xlab("Absolute number of degrees from 32.5 C (Ambient)") + 
  ylab("Absolute number of degrees from 32.5 C (Brood)") + 
  geom_hline(aes(yintercept = 0)) + geom_vline(aes(xintercept = 0))
ee

ee + geom_line(aes(x = amb_abs_dist_from_32, y= preds, linetype = interaction(treatment, dayTime)))


ee +  geom_smooth(aes(linetype =  dayTime, color = treatment), size = 1)

# plot actual vs. predicted
bsmall$preds1 <- predict(gm1, type = 'response')
ggplot(bsmall, aes(x = time, y= preds1)) + 
  facet_wrap(~colony + treatment) + 
  geom_line(aes(y = brood_abs_dist_from_32), color = 'black') +  # black is actual
  geom_line(alpha = 0.5, color = 'red') # red is predicted



## convert back to origin at (0,0)




dev.off()
