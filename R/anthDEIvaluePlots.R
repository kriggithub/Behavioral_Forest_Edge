# 10/20/25 
# Lowest AIC Models Selection (Anth)
library(tidyverse)
library(ggpubr)
library(segmented)
library(strucchange)
library(chngpt)
library(minpack.lm)
library(investr)


anthBinData <- read.csv("data/anthBinData.csv")
# DEI proportions
deiVals <- c(0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9)


########################################################################################
# DistNN
anthBinDataDistNNSub <- anthBinData %>% 
  filter(!is.na(wtAvgDistNN))
########################################################################################

# Power Model AIC 65.95
powerabcDistNN <- nlsLM(wtAvgDistNN ~ a * ((wtAvgAnthDist/400)^b) + c, data = anthBinDataDistNNSub, 
                        start = list(a = 0, b = 1, c = 3), weights = nMonkeys, control = nls.lm.control(maxiter = 1000))
powerabcDistNN$call$control$maxiter <- 1000
powerabcDistNN$call$start <- as.list(coef(powerabcDistNN))



# INCREASING

# max observed x
DistNNxmax <- max(anthBinDataDistNNSub$wtAvgAnthDist, na.rm = TRUE)
DistNNxmin <- min(anthBinDataDistNNSub$wtAvgAnthDist, na.rm = TRUE)



# compute min/max predicted y
DistNNymax <- predict(powerabcDistNN, newdata = data.frame(wtAvgAnthDist = DistNNxmax))
DistNNymin <- predict(powerabcDistNN, newdata = data.frame(wtAvgAnthDist = DistNNxmin))




# Make a results dataframe
DistNNresults <- map_dfr(deiVals, function(d) {
  # define target y0 at that DEI level
  y0 <- DistNNymin + d * (DistNNymax - DistNNymin)
  
  inv <- invest(powerabcDistNN,
                seed = 123,
                y0 = y0,
                interval = "percentile",
                nsim = 10000,
                boot.type = "nonparametric",
                progress = TRUE)
  
  tibble(
    DEI = d,
    estimate = as.numeric(inv$estimate),
    lower = as.numeric(inv$lower),
    upper = as.numeric(inv$upper)
  )
})




DistNNplot <- ggplot(DistNNresults, aes(x = DEI, y = estimate)) +
  geom_point(color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, color = "red") +
  scale_x_continuous(breaks = DistNNresults$DEI,
                     labels = DistNNresults$DEI) +
  labs(x = "Selected DEI threshold",
       y = "Distance from Anthropogenic Edge (m)",
       title = "Distance from Nearest Neighbors DEI Thresholds") +
  theme_bw()









########################################################################################
# NumNN
anthBinDataNumNNSub <- anthBinData %>% 
  filter(!is.na(wtAvgNumNN))
########################################################################################


# linear
linearNumNN <- lm(data = anthBinDataNumNNSub, formula = wtAvgNumNN ~ wtAvgAnthDist, weights = nMonkeys)


# DECREASING

# max and min observed x
NumNNxmin <- min(anthBinDataNumNNSub$wtAvgAnthDist, na.rm = TRUE)
NumNNxmax <- max(anthBinDataNumNNSub$wtAvgAnthDist, na.rm = TRUE)


# model predicted y at that x
NumNNymax <- predict(linearNumNN,
                     newdata = data.frame(wtAvgAnthDist = NumNNxmin))
NumNNymin <- predict(linearNumNN,
                     newdata = data.frame(wtAvgAnthDist = NumNNxmax))



# Make a results dataframe
NumNNresults <- map_dfr(deiVals, function(n) {
  # define target y0 at that DEI level
  y0 <- NumNNymax - n * (NumNNymax - NumNNymin)
  
  inv <- invest(linearNumNN,
                seed = 123,
                y0 = y0,
                interval = "percentile",
                nsim = 10000,  
                boot.type = "nonparametric", 
                progress = TRUE)
  
  tibble(
    DEI = n,
    estimate = as.numeric(inv$estimate),
    lower = as.numeric(inv$lower),
    upper = as.numeric(inv$upper)
  )
})




NumNNplot <- ggplot(NumNNresults, aes(x = DEI, y = estimate)) +
  geom_point(color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, color = "red") +
  scale_x_continuous(breaks = NumNNresults$DEI,
                     labels = NumNNresults$DEI) +
  labs(x = "Selected DEI threshold",
       y = "Distance from Anthropogenic Edge (m)",
       title = "Number of Nearest Neighbors DEI Thresholds") +
  theme_bw()






########################################################################################
# % Time Feeding
anthBinDataFeedingSub <- anthBinData %>% 
  filter(!is.na(wtAvgFeedingPct))
########################################################################################

# linear
linearFeedingPct <- lm(data = anthBinDataFeedingSub, formula = wtAvgFeedingPct ~ wtAvgAnthDist, weights = nMonkeys)



# DECREASING

# max and min observed x
FeedingPctxmin <- min(anthBinDataFeedingSub$wtAvgAnthDist, na.rm = TRUE)
FeedingPctxmax <- max(anthBinDataFeedingSub$wtAvgAnthDist, na.rm = TRUE)


# model predicted y at that x
FeedingPctymax <- predict(linearFeedingPct,
                          newdata = data.frame(wtAvgAnthDist = FeedingPctxmin))
FeedingPctymin <- predict(linearFeedingPct,
                          newdata = data.frame(wtAvgAnthDist = FeedingPctxmax))




# Make a results dataframe
feedingResults <- map_dfr(deiVals, function(f) {
  # define target y0 at that DEI level
  y0 <- FeedingPctymax - f * (FeedingPctymax - FeedingPctymin)
  
  inv <- invest(linearFeedingPct,
                seed = 123,
                y0 = y0,
                interval = "percentile",
                nsim = 10000,  
                boot.type = "nonparametric", 
                progress = TRUE)
  
  tibble(
    DEI = f,
    estimate = as.numeric(inv$estimate),
    lower = as.numeric(inv$lower),
    upper = as.numeric(inv$upper)
  )
})




FeedingPlot <- ggplot(feedingResults, aes(x = DEI, y = estimate)) +
  geom_point(color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, color = "red") +
  scale_x_continuous(breaks = feedingResults$DEI,
                     labels = feedingResults$DEI) +
  labs(x = "Selected DEI threshold",
       y = "Distance from Anthropogenic Edge (m)",
       title = "Feeding % DEI Thresholds") +
  theme_bw()





# plot title
plotTitle <- ggplot() +
  theme_void() +
  annotate("text",
           x = 0, y = 0,
           label = "Anthropogenic Edge DEI Thresholds", 
           hjust = 0.5, vjust = 0, size = 5)




anthDEIthresholdPlots <- ggarrange(DistNNplot,
                                   NumNNplot,
                                   FeedingPlot,
                             plotTitle, ncol = 2, nrow = 2)


#ggexport(anthDEIthresholdPlots, filename = "output/DEIplots/anthDEIthresholdPlots.pdf", height = 10, width = 10)


