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

# create prediction dataframe

predData <- data.frame(
  wtAvgRivDist = seq(min(rivBinData$wtAvgRivDist, na.rm = T),
                     max(rivBinData$wtAvgRivDist, na.rm = T),
                     length.out = 200)
)


########################################################################################
# DistNN
rivBinDataDistNNSub <- rivBinData %>% 
  filter(!is.na(wtAvgDistNN))
########################################################################################

# linear
linearDistNN <- lm(data = rivBinDataDistNNSub, formula = wtAvgDistNN ~ wtAvgRivDist, weights = nMonkeys)
linearDistNNAIC <- AIC(linearDistNN)
predData$linearDistNN <- predict(linearDistNN, newdata = predData)


# INCREASING  

# max and min observed x
DistNNxmin <- min(rivBinDataDistNNSub$wtAvgRivDist, na.rm = TRUE)
DistNNxmax <- max(rivBinDataDistNNSub$wtAvgRivDist, na.rm = TRUE)

# model predicted y at that x
DistNNymax <- predict(linearDistNN,
                     newdata = data.frame(wtAvgRivDist = DistNNxmax))
DistNNymin <- predict(linearDistNN,
                     newdata = data.frame(wtAvgRivDist = DistNNxmin))
# 
DistNNydei <- DistNNymin + (2/3)*(DistNNymax - DistNNymin)


DistNNinvest <- invest(linearDistNN,
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



linearDistNNplot <-ggplot(rivBinDataDistNNSub, aes(x = wtAvgRivDist, y = wtAvgDistNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgDistNN - wtSeDistNN, ymax = wtAvgDistNN + wtSeDistNN)) +
  labs(x = "Distance to River Edge (m)", 
       y = "Mean Distance to Nearest Neighbors" 
       # title = paste0("Distance from Nearest Neighbors (Linear AIC = ", round(linearDistNNAIC, 2), ")")
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = linearDistNN)) +  
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


linearDistNNplot <- linearDistNNplot + axis_theme + spacing_theme



########################################################################################
# NumNN
rivBinDataNumNNSub <- rivBinData %>% 
  filter(!is.na(wtAvgNumNN))
########################################################################################



# logistic
logisticNumNN <- nlsLM(wtAvgNumNN ~ a/(1+(b * exp(-c*(wtAvgRivDist)/400))) + d, data = rivBinDataNumNNSub, 
                       start = list(a = 10, b = 10, c = 10, d = 1), weights = nMonkeys, control = nls.lm.control(maxiter = 1000))
logisticNumNNAIC <- AIC(logisticNumNN)
predData$logisticNumNN <- predict(logisticNumNN, newdata = predData)



# model coefficients to find point of inflection for logistic model
NumNNCoefs <- coef(logisticNumNN)
NumNNb <- NumNNCoefs[["b"]]
NumNNc <- NumNNCoefs[["c"]]
# for asymptotes
NumNNa <- NumNNCoefs[["a"]]
NumNNd <- NumNNCoefs[["d"]]

NumNNd # upper
NumNNa + NumNNd # lower

NumNNvcovMat <- vcov(logisticNumNN)



# lower asymptote d
se_lower_dm <- deltamethod(~ x4,
                           mean = NumNNCoefs,
                           cov  = NumNNvcovMat)
lower_est <- NumNNd
lower_CI  <- lower_est + c(-1.96, 1.96) * se_lower_dm

# upper asymptote a + d
se_upper_dm <- deltamethod(~ x1 + x4,
                           mean = NumNNCoefs,
                           cov  = NumNNvcovMat)
upper_est <- NumNNa + NumNNd
upper_CI  <- upper_est + c(-1.96, 1.96) * se_upper_dm








# Delta method
NumNNse <- deltamethod(~ (400/x3) * log(x2), 
                         mean = NumNNCoefs, 
                         cov = NumNNvcovMat)




NumNNPE <- ((log(NumNNb)*400)/NumNNc) # Point of inflection calculation from model coefficients
NumNNlower <- NumNNPE - 1.96*NumNNse
NumNNupper <- NumNNPE + 1.96*NumNNse

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




logisticNumNNplot <-ggplot(rivBinDataNumNNSub, aes(x = wtAvgRivDist, y = wtAvgNumNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgNumNN - wtSeNumNN, ymax = wtAvgNumNN + wtSeNumNN)) +
  labs(x = "Distance to River Edge (m)", 
       y = "Mean Number of Nearest Neighbors" 
       # title = paste0("Number of Nearest Neighbors (Logistic AIC = ", round(logisticNumNNAIC, 2), ")")
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = logisticNumNN)) +
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


logisticNumNNplot <- logisticNumNNplot + axis_theme + spacing_theme



########################################################################################
# % Time Feeding
rivBinDataFeedingSub <- rivBinData %>% 
  filter(!is.na(wtAvgFeedingPct))
########################################################################################


# null
nullFeedingPct <- lm(data = rivBinDataFeedingSub, formula = wtAvgFeedingPct ~ 1, weights = nMonkeys)
nullFeedingPctAIC <- AIC(nullFeedingPct)
predData$nullFeedingPct <- predict(nullFeedingPct, newdata = predData)

nullFeedingPctplot <- ggplot(rivBinDataFeedingSub, aes(x = wtAvgRivDist, y = wtAvgFeedingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgFeedingPct - wtSeFeedingPct, ymax = wtAvgFeedingPct + wtSeFeedingPct)) +
  labs(x = "Distance to River Edge (m)", 
       y = "Mean Percent Time Spent Feeding" 
       # title = paste0("Feeding % (Null AIC = ", round(nullFeedingPctAIC, 2), ")"),
       # caption = "No DEI Effects"
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = nullFeedingPct)) +
  theme(
    legend.position = "none"                            # remove legend
    # plot.caption = element_text(hjust = 0.5, size = 10)  
  )

nullFeedingPctplot <- nullFeedingPctplot + axis_theme + spacing_theme


########################################################################################
# % Time Moving
rivBinDataMovingSub <- rivBinData %>% 
  filter(!is.na(wtAvgMovingPct))
########################################################################################



# segmented
linearMovingPct <- lm(data = rivBinDataMovingSub, formula = wtAvgMovingPct ~ wtAvgRivDist, weights = nMonkeys)
segmentedMovingPct <- segmented(linearMovingPct, seg.Z = ~ wtAvgRivDist, psi = 250)
segmentedMovingPctAIC<- AIC(segmentedMovingPct)
predData$segmentedMovingPct <- predict(segmentedMovingPct, newdata = predData)




MovingPctPE <- confint(segmentedMovingPct)[1]
MovingPctlower <- confint(segmentedMovingPct)[2]
MovingPctupper <- confint(segmentedMovingPct)[3]



MovingPctciLabel <- paste0("95% CI = (", round(MovingPctlower, 1), ", ", round(MovingPctupper, 1), ")")
MovingPctpeLabel <- paste0("Point estimate = ", round(MovingPctPE, 1))

MovingPctciBand <- data.frame(
  xmin = MovingPctlower,
  xmax = MovingPctupper,
  ymin = -Inf,
  ymax = Inf,
  type = MovingPctciLabel
)

MovingPctpointLine <- data.frame(
  xintercept = MovingPctPE,
  type = MovingPctpeLabel
)




segmentedMovingPctplot <-ggplot(rivBinDataMovingSub, aes(x = wtAvgRivDist, y = wtAvgMovingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgMovingPct - wtSeMovingPct, ymax = wtAvgMovingPct + wtSeMovingPct)) +
  labs(x = "Distance to River Edge (m)", 
       y = "Mean Percent Time Spent Moving" 
       # title = paste0("Moving % (Segmented AIC = ", round(segmentedMovingPctAIC, 2), ")")
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = segmentedMovingPct)) +
  geom_rect(data = MovingPctciBand,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = type),
            alpha = 0.2, inherit.aes = FALSE) +
  geom_vline(data = MovingPctpointLine,
             aes(xintercept = xintercept, color = type),
             linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = setNames("red", MovingPctciLabel)) +
  scale_color_manual(values = setNames("red", MovingPctpeLabel)) +
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


segmentedMovingPctplot <- segmentedMovingPctplot + axis_theme + spacing_theme


########################################################################################
# % Time Resting
rivBinDataRestSub <- rivBinData %>% 
  filter(!is.na(wtAvgRestPct))
########################################################################################


# linear
linearRestPct <- lm(data = rivBinDataRestSub, formula = wtAvgRestPct ~ wtAvgRivDist, weights = nMonkeys)
linearRestPctAIC <- AIC(linearRestPct)
predData$linearRestPct <- predict(linearRestPct, newdata = predData)



# Decreasing  

# max and min observed x
RestPctxmin <- min(rivBinDataRestSub$wtAvgRivDist, na.rm = TRUE)
RestPctxmax <- max(rivBinDataRestSub$wtAvgRivDist, na.rm = TRUE)

# model predicted y at that x
RestPctymax <- predict(linearRestPct,
                      newdata = data.frame(wtAvgRivDist = RestPctxmin))
RestPctymin <- predict(linearRestPct,
                      newdata = data.frame(wtAvgRivDist = RestPctxmax))
# 
RestPctydei <- RestPctymax - (2/3)*(RestPctymax - RestPctymin)

RestPctinvest <- invest(linearRestPct,
                       seed = 123,
                       y0 = RestPctydei,
                       interval = "percentile",
                       nsim = 10000,  
                       boot.type = "nonparametric", 
                       progress = TRUE)


RestPctPE <- as.numeric(RestPctinvest$estimate)
RestPctlower <- as.numeric(RestPctinvest$lower)
RestPctupper <- as.numeric(RestPctinvest$upper)

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




linearRestPctplot <-ggplot(rivBinDataRestSub, aes(x = wtAvgRivDist, y = wtAvgRestPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgRestPct - wtSeRestPct, ymax = wtAvgRestPct + wtSeRestPct)) +
  labs(x = "Distance to River Edge (m)", 
       y = "Mean Percent Time Spent Resting" 
       # title = paste0("Resting % (Linear AIC = ", round(linearRestPctAIC, 2), ")")
  ) +
  theme_bw() +
  geom_line(data = predData, aes(y = linearRestPct)) +
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


linearRestPctplot <- linearRestPctplot + axis_theme + spacing_theme






allDEIplotsRiv <- ggarrange(linearDistNNplot,
                            logisticNumNNplot,
                            nullFeedingPctplot,
                            segmentedMovingPctplot,
                            linearRestPctplot,
                            ncol = 2, nrow = 3,
                            labels   = c("a", "b", "c", "d", "e"),        # panel letters
                            label.x  = 0.02,                               # a little inset from left
                            label.y  = 0.98,                               # near the top
                            hjust    = 0,                                  # left aligned
                            vjust    = 1,                                  # top aligned
                            font.label = list(size = 20, face = "bold"))


# save.image(file = "output/rivDEImodels.RData")


ggexport(allDEIplotsRiv, filename = "output/allDEIplotsRivAdj.pdf", height = 15, width = 11)



