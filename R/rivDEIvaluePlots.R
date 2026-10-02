# 10/20/25 
# Lowest AIC Models Selection (Riv)
library(tidyverse)
library(ggpubr)
library(segmented)
library(strucchange)
library(chngpt)
library(minpack.lm)
library(investr)
library(msm)

rivBinData <- read.csv("data/rivBinData.csv")
# DEI proportions
deiVals <- c(0.1, 0.25, 0.33, 0.5, 0.66, 0.75, 0.9)

########################################################################################
# DistNN
rivBinDataDistNNSub <- rivBinData %>% 
  filter(!is.na(wtAvgDistNN))
########################################################################################

# linear
linearDistNN <- lm(data = rivBinDataDistNNSub, formula = wtAvgDistNN ~ wtAvgRivDist, weights = nMonkeys)


# INCREASING  

# max and min observed x
DistNNxmin <- min(rivBinDataDistNNSub$wtAvgRivDist, na.rm = TRUE)
DistNNxmax <- max(rivBinDataDistNNSub$wtAvgRivDist, na.rm = TRUE)

# model predicted y at that x
DistNNymax <- predict(linearDistNN,
                      newdata = data.frame(wtAvgRivDist = DistNNxmax))
DistNNymin <- predict(linearDistNN,
                      newdata = data.frame(wtAvgRivDist = DistNNxmin))





# Make a results dataframe
DistNNresults <- map_dfr(deiVals, function(d) {
  # define target y0 at that DEI level
  y0 <- DistNNymin + d * (DistNNymax - DistNNymin)
  
  inv <- invest(linearDistNN,
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
       y = "Distance from River Edge (m)",
       title = "Distance from Nearest Neighbors DEI Thresholds") +
  theme_bw()







########################################################################################
# % Time Resting
rivBinDataRestSub <- rivBinData %>% 
  filter(!is.na(wtAvgRestPct))
########################################################################################


# linear
linearRestPct <- lm(data = rivBinDataRestSub, formula = wtAvgRestPct ~ wtAvgRivDist, weights = nMonkeys)



# Decreasing  

# max and min observed x
RestPctxmin <- min(rivBinDataRestSub$wtAvgRivDist, na.rm = TRUE)
RestPctxmax <- max(rivBinDataRestSub$wtAvgRivDist, na.rm = TRUE)

# model predicted y at that x
RestPctymax <- predict(linearRestPct,
                       newdata = data.frame(wtAvgRivDist = RestPctxmin))
RestPctymin <- predict(linearRestPct,
                       newdata = data.frame(wtAvgRivDist = RestPctxmax))



# Make a results dataframe
restingResults <- map_dfr(deiVals, function(r) {
  # define target y0 at that DEI level
  y0 <- RestPctymax - r * (RestPctymax - RestPctymin)
  
  inv <- invest(linearRestPct,
                seed = 123,
                y0 = y0,
                interval = "percentile",
                nsim = 10000,  
                boot.type = "nonparametric", 
                progress = TRUE)
  
  tibble(
    DEI = r,
    estimate = as.numeric(inv$estimate),
    lower = as.numeric(inv$lower),
    upper = as.numeric(inv$upper)
  )
})




RestingPlot <- ggplot(restingResults, aes(x = DEI, y = estimate)) +
  geom_point(color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, color = "red") +
  scale_x_continuous(breaks = restingResults$DEI,
                     labels = restingResults$DEI) +
  labs(x = "Selected DEI threshold",
       y = "Distance from Anthropogenic Edge (m)",
       title = "Feeding % DEI Thresholds") +
  theme_bw()




# plot title
plotTitle <- ggplot() +
  theme_void() +
  annotate("text",
           x = 0, y = 0,
           label = "River Edge DEI Thresholds", 
           hjust = 0.5, vjust = 0, size = 5)




rivDEIthresholdPlots <- ggarrange(DistNNplot,
                                   RestingPlot,
                                   plotTitle, ncol = 2, nrow = 2)


#ggexport(rivDEIthresholdPlots, filename = "output/DEIplots/rivDEIthresholdPlots.pdf", height = 10, width = 10)

