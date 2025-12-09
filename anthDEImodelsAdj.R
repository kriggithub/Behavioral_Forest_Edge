# 10/20/25 
# Lowest AIC Models Selection (Anth)
library(tidyverse)
library(ggpubr)
library(segmented)
library(strucchange)
library(chngpt)
library(minpack.lm)
library(investr)
library(msm)



anthBinData <- read.csv("anthBinData.csv")

# create prediction dataframe

predData <- data.frame(
  wtAvgAnthDist = seq(min(anthBinData$wtAvgAnthDist, na.rm = T),
                      max(anthBinData$wtAvgAnthDist, na.rm = T),
                      length.out = 200)
)



########################################################################################
# DistNN
anthBinDataDistNNSub <- anthBinData %>% 
  filter(!is.na(wtAvgDistNN))
########################################################################################

# Power Model AIC 65.95
powerabcDistNN <- nlsLM(wtAvgDistNN ~ a * ((wtAvgAnthDist/400)^b) + c, data = anthBinDataDistNNSub, 
                        start = list(a = 0, b = 1, c = 3), weights = nMonkeys, control = nls.lm.control(maxiter = 1000))
powerabcDistNNAIC <- AIC(powerabcDistNN)


predData$powerabcDistNN <- predict(powerabcDistNN, newdata = predData)


powerabcDistNN$call$control$maxiter <- 1000
powerabcDistNN$call$start <- as.list(coef(powerabcDistNN))

# INCREASING

# max observed x
DistNNxmax <- max(anthBinDataDistNNSub$wtAvgAnthDist, na.rm = TRUE)
DistNNxmin <- min(anthBinDataDistNNSub$wtAvgAnthDist, na.rm = TRUE)


# model predicted y at that x
DistNNymax <- predict(powerabcDistNN,
                         newdata = data.frame(wtAvgAnthDist = DistNNxmax))
DistNNymin <- predict(powerabcDistNN,
                      newdata = data.frame(wtAvgAnthDist = DistNNxmin))
# 
DistNNydei <- DistNNymin + (2/3)*(DistNNymax - DistNNymin)

DistNNinvest <- invest(powerabcDistNN,
                  seed = 123,
                  y0 = DistNNydei,
                  interval = "percentile",
                  nsim = 10000,  
                  boot.type = "nonparametric", 
                  progress = TRUE)


DistNNPE <- as.numeric(DistNNinvest$estimate)
DistNNlower <- as.numeric(DistNNinvest$lower)
DistNNupper <- as.numeric(DistNNinvest$upper)

DistNNciLabel <- paste0("95% CI = (", round(DistNNlower, 1), ", ", round(DistNNupper, 1), ")")
DistNNpeLabel <- paste0("Point estimate = ", round(DistNNPE, 1))

DistNNciBand <- data.frame(
  xmin = DistNNlower,
  xmax = DistNNupper,
  ymin = -Inf,
  ymax = Inf,
  type = DistNNciLabel
)

DistNNpointLine <- data.frame(
  xintercept = DistNNPE,
  type = DistNNpeLabel
)


axis_theme <- theme(
  axis.title = element_text(size = 16),  # axis labels
  axis.text  = element_text(size = 12)   # tick labels
)

spacing_theme <- theme(
  plot.margin = margin(t = 20, r = 20, b = 20, l = 20)  # units are points
)




powerabcDistNNplot <- ggplot(anthBinDataDistNNSub, aes(x = wtAvgAnthDist, y = wtAvgDistNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgDistNN - wtSeDistNN, ymax = wtAvgDistNN + wtSeDistNN)) +
  labs(x = "Distance to Anthropogenic Edge (m)", 
       y = "Mean Distance to Nearest Neighbors" 
       # title = paste0("Distance from Nearest Neighbors (Power AIC = ", round(powerabcDistNNAIC, 2),")")
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = powerabcDistNN)) +
  geom_rect(data = DistNNciBand,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_vline(data = DistNNpointLine,
             aes(xintercept = xintercept, color = type),
             linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = setNames("red", DistNNciLabel)) +
  scale_color_manual(values = setNames("red", DistNNpeLabel)) +
  labs(fill = NULL, color = NULL) +
  theme(
    # legend.position = "bottom",               
    # legend.justification = "center",          
    # legend.direction = "horizontal",          
    # legend.box = "horizontal",                
    # legend.background = element_blank(),
    # legend.key = element_blank()
    legend.position = "none"
  )

powerabcDistNNplot <- powerabcDistNNplot + axis_theme + spacing_theme


########################################################################################
# NumNN
anthBinDataNumNNSub <- anthBinData %>% 
  filter(!is.na(wtAvgNumNN))
########################################################################################


# linear
linearNumNN <- lm(data = anthBinDataNumNNSub, formula = wtAvgNumNN ~ wtAvgAnthDist, weights = nMonkeys)
linearNumNNAIC <- AIC(linearNumNN)
predData$linearNumNN <- predict(linearNumNN, newdata = predData)



# DECREASING

# max and min observed x
NumNNxmin <- min(anthBinDataNumNNSub$wtAvgAnthDist, na.rm = TRUE)
NumNNxmax <- max(anthBinDataNumNNSub$wtAvgAnthDist, na.rm = TRUE)


# model predicted y at that x
NumNNymax <- predict(linearNumNN,
                     newdata = data.frame(wtAvgAnthDist = NumNNxmin))
NumNNymin <- predict(linearNumNN,
                      newdata = data.frame(wtAvgAnthDist = NumNNxmax))
# 
NumNNydei <- NumNNymax - (2/3)*(NumNNymax - NumNNymin)

NumNNinvest <- invest(linearNumNN,
                       seed = 123,
                       y0 = NumNNydei,
                       interval = "percentile",
                       nsim = 10000,  
                       boot.type = "nonparametric", 
                       progress = TRUE)


NumNNPE <- as.numeric(NumNNinvest$estimate)
NumNNlower <- as.numeric(NumNNinvest$lower)
NumNNupper <- as.numeric(NumNNinvest$upper)

NumNNciLabel <- paste0("95% CI = (", round(NumNNlower, 1), ", ", round(NumNNupper, 1), ")")
NumNNpeLabel <- paste0("Point estimate = ", round(NumNNPE, 1))

NumNNciBand <- data.frame(
  xmin = NumNNlower,
  xmax = NumNNupper,
  ymin = -Inf,
  ymax = Inf,
  type = NumNNciLabel
)

NumNNpointLine <- data.frame(
  xintercept = NumNNPE,
  type = NumNNpeLabel
)




linearNumNNplot <-ggplot(anthBinDataNumNNSub, aes(x = wtAvgAnthDist, y = wtAvgNumNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgNumNN - wtSeNumNN, ymax = wtAvgNumNN + wtSeNumNN)) +
  labs(x = "Distance to Anthropogenic Edge (m)", 
       y = "Mean Number of Nearest Neighbors" 
       # title = paste0("Number of Nearest Neighbors (Linear AIC = ", round(linearNumNNAIC, 2),")")
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = linearNumNN)) +
  geom_rect(data = NumNNciBand,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_vline(data = NumNNpointLine,
             aes(xintercept = xintercept, color = type),
             linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = setNames("red", NumNNciLabel)) +
  scale_color_manual(values = setNames("red", NumNNpeLabel)) +
  labs(fill = NULL, color = NULL) + 
  theme(
    # legend.position = "bottom",               
    # legend.justification = "center",          
    # legend.direction = "horizontal",          
    # legend.box = "horizontal",                
    # legend.background = element_blank(),
    # legend.key = element_blank()
    legend.position = "none"
  )

linearNumNNplot <- linearNumNNplot + axis_theme + spacing_theme


########################################################################################
# % Time Feeding
anthBinDataFeedingSub <- anthBinData %>% 
  filter(!is.na(wtAvgFeedingPct))
########################################################################################

# linear
linearFeedingPct <- lm(data = anthBinDataFeedingSub, formula = wtAvgFeedingPct ~ wtAvgAnthDist, weights = nMonkeys)
linearFeedingPctAIC <- AIC(linearFeedingPct)
predData$linearFeedingPct <- predict(linearFeedingPct, newdata = predData)




# DECREASING

# max and min observed x
FeedingPctxmin <- min(anthBinDataFeedingSub$wtAvgAnthDist, na.rm = TRUE)
FeedingPctxmax <- max(anthBinDataFeedingSub$wtAvgAnthDist, na.rm = TRUE)


# model predicted y at that x
FeedingPctymax <- predict(linearFeedingPct,
                     newdata = data.frame(wtAvgAnthDist = FeedingPctxmin))
FeedingPctymin <- predict(linearFeedingPct,
                     newdata = data.frame(wtAvgAnthDist = FeedingPctxmax))
# 
FeedingPctydei <- FeedingPctymax - (2/3)*(FeedingPctymax - FeedingPctymin)

FeedingPctinvest <- invest(linearFeedingPct,
                      seed = 123,
                      y0 = FeedingPctydei,
                      interval = "percentile",
                      nsim = 10000,  
                      boot.type = "nonparametric", 
                      progress = TRUE)


FeedingPctPE <- as.numeric(FeedingPctinvest$estimate)
FeedingPctlower <- as.numeric(FeedingPctinvest$lower)
FeedingPctupper <- as.numeric(FeedingPctinvest$upper)

FeedingPctciLabel <- paste0("95% CI = (", round(FeedingPctlower, 1), ", ", round(FeedingPctupper, 1), ")")
FeedingPctpeLabel <- paste0("Point estimate = ", round(FeedingPctPE, 1))

FeedingPctciBand <- data.frame(
  xmin = FeedingPctlower,
  xmax = FeedingPctupper,
  ymin = -Inf,
  ymax = Inf,
  type = FeedingPctciLabel
)

FeedingPctpointLine <- data.frame(
  xintercept = FeedingPctPE,
  type = FeedingPctpeLabel
)




linearFeedingPctplot <-ggplot(anthBinDataFeedingSub, aes(x = wtAvgAnthDist, y = wtAvgFeedingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgFeedingPct - wtSeFeedingPct, ymax = wtAvgFeedingPct + wtSeFeedingPct)) +
  labs(x = "Distance to Anthropogenic Edge (m)", 
       y = "Mean Percent Time Spent Feeding" 
       # title = paste0("Feeding % (Linear AIC = ", round(linearFeedingPctAIC, 2),")")
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = linearFeedingPct)) +
  geom_rect(data = FeedingPctciBand,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_vline(data = FeedingPctpointLine,
             aes(xintercept = xintercept, color = type),
             linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = setNames("red", FeedingPctciLabel)) +
  scale_color_manual(values = setNames("red", FeedingPctpeLabel)) +
  labs(fill = NULL, color = NULL) + 
  theme(
    # legend.position = "bottom",               
    # legend.justification = "center",          
    # legend.direction = "horizontal",          
    # legend.box = "horizontal",                
    # legend.background = element_blank(),
    # legend.key = element_blank()
    legend.position = "none"
  )

linearFeedingPctplot <- linearFeedingPctplot + axis_theme + spacing_theme


########################################################################################
# % Time Moving
anthBinDataMovingSub <- anthBinData %>% 
  filter(!is.na(wtAvgMovingPct))
########################################################################################

# null
nullMovingPct <- lm(data = anthBinDataMovingSub, formula = wtAvgMovingPct ~ 1, weights = nMonkeys)
nullMovingPctAIC <- AIC(nullMovingPct)
predData$nullMovingPct <- predict(nullMovingPct, newdata = predData)




nullMovingPctplot <- ggplot(anthBinDataMovingSub, aes(x = wtAvgAnthDist, y = wtAvgMovingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgMovingPct - wtSeMovingPct, ymax = wtAvgMovingPct + wtSeMovingPct)) +
  labs(x = "Distance to Anthropogenic Edge (m)", 
       y = "Mean Percent Time Spent Moving" 
       # title = paste0("Moving % (Null AIC = ", round(nullMovingPctAIC, 2), ")"),
       # caption = "No DEI Effects"
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = nullMovingPct)) +
  theme(
    # legend.position = "bottom",
    # legend.box = "horizontal",
    # legend.key = element_blank(),
    # legend.background = element_blank(),                           # remove legend
    # plot.caption = element_text(hjust = 0.5, size = 10)  
    legend.position = "none"
  )

nullMovingPctplot <- nullMovingPctplot + axis_theme + spacing_theme


########################################################################################
# % Time Resting
anthBinDataRestSub <- anthBinData %>% 
  filter(!is.na(wtAvgRestPct))
########################################################################################



# logistic
logisticRestPct <- nlsLM(wtAvgRestPct ~ a/(1+(b * exp(-c*(wtAvgAnthDist-100)/400))) + d, data = anthBinDataRestSub, 
                         start = list(a = 4, b = 240, c = 105, d = 80), weights = nMonkeys, control = nls.lm.control(maxiter = 1000))
logisticRestPctAIC <- AIC(logisticRestPct)
predData$logisticRestPct <- predict(logisticRestPct, newdata = predData)




# model coefficients to find point of inflection for logistic model
RestCoefs <- coef(logisticRestPct)
restb <- RestCoefs[["b"]]
restc <- RestCoefs[["c"]]
resta <- RestCoefs[["a"]]
restd <- RestCoefs[["d"]]

restd

restd + resta


vcovMat <- vcov(logisticRestPct)



# lower asymptote d
se_lower_dm <- deltamethod(~ x4,
                           mean = RestCoefs,
                           cov  = vcovMat)
lower_est <- restd
lower_CI  <- lower_est + c(-1.96, 1.96) * se_lower_dm

# upper asymptote a + d
se_upper_dm <- deltamethod(~ x1 + x4,
                           mean = RestCoefs,
                           cov  = vcovMat)
upper_est <- restd + resta
upper_CI  <- upper_est + c(-1.96, 1.96) * se_upper_dm






# Delta method
RestPCTse <- deltamethod(~ 100 + (400/x3) * log(x2), 
                          mean = RestCoefs, 
                          cov = vcovMat)




RestPctPE <- ((log(restb)*400)/restc) + 100 # Point of inflection calculation from model coefficients
RestPctlower <- RestPctPE - 1.96*RestPCTse
RestPctupper <- RestPctPE + 1.96*RestPCTse

RestPctciLabel <- paste0("95% CI = (", round(RestPctlower, 1), ", ", round(RestPctupper, 1), ")")
RestPctpeLabel <- paste0("Point estimate = ", round(RestPctPE, 1))

RestPctciBand <- data.frame(
  xmin = RestPctlower,
  xmax = RestPctupper,
  ymin = -Inf,
  ymax = Inf,
  type = RestPctciLabel
)

RestPctpointLine <- data.frame(
  xintercept = RestPctPE,
  type = RestPctpeLabel
)




logisticRestPctplot <-ggplot(anthBinDataRestSub, aes(x = wtAvgAnthDist, y = wtAvgRestPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgRestPct - wtSeRestPct, ymax = wtAvgRestPct + wtSeRestPct)) +
  labs(x = "Distance to Anthropogenic Edge (m)", 
       y = "Mean Percent Time Spent Resting"
       # title = paste0("Resting % (Logistic AIC = ", round(logisticRestPctAIC, 2), ")")
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = logisticRestPct)) +
  geom_rect(data = RestPctciBand,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_vline(data = RestPctpointLine,
             aes(xintercept = xintercept, color = type),
             linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = setNames("red", RestPctciLabel)) +
  scale_color_manual(values = setNames("red", RestPctpeLabel)) +
  labs(fill = NULL, color = NULL) +
  theme(
    # legend.position = "bottom",               
    # legend.justification = "center",          
    # legend.direction = "horizontal",          
    # legend.box = "horizontal",                
    # legend.background = element_blank(),
    # legend.key = element_blank()
    legend.position = "none"
  )


logisticRestPctplot <- logisticRestPctplot + axis_theme + spacing_theme







allDEIplotsAnth <- ggarrange(powerabcDistNNplot,
                             linearNumNNplot,
                             linearFeedingPctplot,
                             nullMovingPctplot,
                             logisticRestPctplot,
                             ncol = 2, nrow = 3,
                             labels   = c("a", "b", "c", "d", "e"),        # panel letters
                             label.x  = 0.02,                               # a little inset from left
                             label.y  = 0.98,                               # near the top
                             hjust    = 0,                                  # left aligned
                             vjust    = 1,                                  # top aligned
                             font.label = list(size = 20, face = "bold"))


allDEIplotsAnth

# save.image(file = "anthDEImodels.RData")


ggexport(allDEIplotsAnth, filename = "allDEIplotsAnthAdj.pdf", height = 15, width = 11)





