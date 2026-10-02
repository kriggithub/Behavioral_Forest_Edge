library(ggplot2)
library(ggpubr)


# monkeyIdData, AnthDist
ma1 <- ggplot(monkeyIdData, mapping = aes(x = AnthDist, y = RestPct,)) +
  geom_point() +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "% Time Resting", 
       title = "% Time Resting") +
  theme_bw()

ma2 <- ggplot(monkeyIdData, mapping = aes(x = AnthDist, y = MovingPct,)) +
  geom_point() +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "% Time Moving", 
       title = "% Time Moving") +
  theme_bw()


ma3 <- ggplot(monkeyIdData, mapping = aes(x = AnthDist, y = FeedingPct,)) +
  geom_point() +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "% Time Feeding", 
       title = "% Time Feeding") +
  theme_bw()




ma4 <- ggplot(monkeyIdData, mapping = aes(x = AnthDist, y = AvgNumNN,)) +
  geom_point() +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "Mean # of Nearest Neighbors", 
       title = "Average Number of Nearest Neighbors") +
  theme_bw()

ma5 <- ggplot(monkeyIdData, mapping = aes(x = AnthDist, y = AvgDistNN,)) +
  geom_point() +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "Mean Distance from Nearest Neighbor", 
       title = "Average Distance from Nearest Neighbors") +
  theme_bw()


matitle <- ggplot() +
  theme_void() +
  annotate("text",
           x = 0, y = 0,
           label = "Anthropogenic Edge Effects per Monkey", 
           hjust = 0.5, vjust = 0, size = 5)

monkeyAnthPlots <- ggarrange(ma1, ma2, ma3, ma4, ma5, matitle, ncol = 3, nrow = 2)
ggexport(monkeyAnthPlots, filename = "output/initialPlots/monkeyAnthPlots.pdf", height = 10, width = 15)





# monkeyIdData, RivDist
mr1 <- ggplot(monkeyIdData, mapping = aes(x = RivDist, y = RestPct,)) +
  geom_point() +
  labs(x = "Distance from River Edge (m)", 
       y = "% Time Resting", 
       title = "% Time Resting") +
  theme_bw()

mr2 <- ggplot(monkeyIdData, mapping = aes(x = RivDist, y = MovingPct,)) +
  geom_point() +
  labs(x = "Distance from River Edge (m)", 
       y = "% Time Moving", 
       title = "% Time Moving") +
  theme_bw()


mr3 <- ggplot(monkeyIdData, mapping = aes(x = RivDist, y = FeedingPct,)) +
  geom_point() +
  labs(x = "Distance from River Edge (m)", 
       y = "% Time Feeding", 
       title = "% Time Feeding") +
  theme_bw()


mr4 <- ggplot(monkeyIdData, mapping = aes(x = RivDist, y = AvgNumNN,)) +
  geom_point() +
  labs(x = "Distance from River Edge (m)", 
       y = "Mean # of Nearest Neighbors", 
       title = "Average Number of Nearest Neighbors") +
  theme_bw()

mr5 <- ggplot(monkeyIdData, mapping = aes(x = RivDist, y = AvgDistNN,)) +
  geom_point() +
  labs(x = "Distance from River Edge (m)", 
       y = "Mean Distance from Nearest Neighbor", 
       title = "Average Distance from Nearest Neighbors") +
  theme_bw()

mrtitle <- ggplot() +
  theme_void() +
  annotate("text",
           x = 0, y = 0,
           label = "River Edge Effects per Monkey", 
           hjust = 0.5, vjust = 0, size = 5)


monkeyRivPlots <- ggarrange(mr1, mr2, mr3, mr4, mr5, mrtitle, ncol = 3, nrow = 2)
ggexport(monkeyRivPlots, filename = "output/initialPlots/monkeyRivPlots.pdf", height = 10, width = 15)






# anthBinData

ba1 <- ggplot(anthBinData, aes(x = wtAvgAnthDist, y = wtAvgRestPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgRestPct - wtSeRestPct, ymax = wtAvgRestPct + wtSeRestPct)) +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "(Weighted) Mean % Time Spent Resting", 
       title = "% Time Spent Resting") +
  theme_bw()

ba2 <- ggplot(anthBinData, aes(x = wtAvgAnthDist, y = wtAvgMovingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgMovingPct - wtSeMovingPct, ymax = wtAvgMovingPct + wtSeMovingPct)) +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "(Weighted) Mean % Time Spent Moving", 
       title = "% Time Spent Moving") +
  theme_bw()



ba3 <- ggplot(anthBinData, aes(x = wtAvgAnthDist, y = wtAvgFeedingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgFeedingPct - wtSeFeedingPct, ymax = wtAvgFeedingPct + wtSeFeedingPct)) +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "(Weighted) Mean % Time Spent Feeding", 
       title = "% Time Spent Feeding") +
  theme_bw()



ba4 <- ggplot(anthBinData, aes(x = wtAvgAnthDist, y = wtAvgNumNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgNumNN - wtSeNumNN, ymax = wtAvgNumNN + wtSeNumNN)) +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "(Weighted) Mean # of Nearest Neighbors", 
       title = "Average Number of Nearest Neighbors") +
  theme_bw()

ba5 <- ggplot(anthBinData, aes(x = wtAvgAnthDist, y = wtAvgDistNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgDistNN - wtSeDistNN, ymax = wtAvgDistNN + wtSeDistNN)) +
  labs(x = "Distance from Anthropogenic Edge (m)", 
       y = "(Weighted) Mean Distance from Nearest Neighbor", 
       title = "Average Distance from Nearest Neighbors") +
  theme_bw()

batitle <- ggplot() +
  theme_void() +
  annotate("text",
           x = 0, y = 0,
           label = "Anthropogenic Edge Effects per 15 meters from Edge", 
           hjust = 0.5, vjust = 0, size = 5)


binAnthPlots <- ggarrange(ba1, ba2, ba3, ba4, ba5, batitle, ncol = 3, nrow = 2)
ggexport(binAnthPlots, filename = "output/initialPlots/binAnthPlots.pdf", height = 10, width = 15)



# rivBinData


ggplot(rivBinData, aes(x = bin, y = wtAvgRestPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgRestPct - wtSeRestPct, ymax = wtAvgRestPct + wtSeRestPct))

ggplot(rivBinData, aes(x = bin, y = wtAvgMovingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgMovingPct - wtSeMovingPct, ymax = wtAvgMovingPct + wtSeMovingPct))

ggplot(rivBinData, aes(x = bin, y = wtAvgFeedingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgFeedingPct - wtSeFeedingPct, ymax = wtAvgFeedingPct + wtSeFeedingPct))

ggplot(rivBinData, aes(x = bin, y = wtAvgNumNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgNumNN - wtSeNumNN, ymax = wtAvgNumNN + wtSeNumNN))

ggplot(rivBinData, aes(x = bin, y = wtAvgDistNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgDistNN - wtSeDistNN, ymax = wtAvgDistNN + wtSeDistNN))




br1 <- ggplot(rivBinData, aes(x = wtAvgRivDist, y = wtAvgRestPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgRestPct - wtSeRestPct, ymax = wtAvgRestPct + wtSeRestPct)) +
  labs(x = "Distance from River Edge (m)", 
       y = "(Weighted) Mean % Time Spent Resting", 
       title = "% Time Spent Resting") +
  theme_bw()

br2 <- ggplot(rivBinData, aes(x = wtAvgRivDist, y = wtAvgMovingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgMovingPct - wtSeMovingPct, ymax = wtAvgMovingPct + wtSeMovingPct)) +
  labs(x = "Distance from River Edge (m)", 
       y = "(Weighted) Mean % Time Spent Moving", 
       title = "% Time Spent Moving") +
  theme_bw()



br3 <- ggplot(rivBinData, aes(x = wtAvgRivDist, y = wtAvgFeedingPct)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgFeedingPct - wtSeFeedingPct, ymax = wtAvgFeedingPct + wtSeFeedingPct)) +
  labs(x = "Distance from River Edge (m)", 
       y = "(Weighted) Mean % Time Spent Feeding", 
       title = "% Time Spent Feeding") +
  theme_bw()



br4 <- ggplot(rivBinData, aes(x = wtAvgRivDist, y = wtAvgNumNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgNumNN - wtSeNumNN, ymax = wtAvgNumNN + wtSeNumNN)) +
  labs(x = "Distance from River Edge (m)", 
       y = "(Weighted) Mean # of Nearest Neighbors", 
       title = "Average Number of Nearest Neighbors") +
  theme_bw()

br5 <- ggplot(rivBinData, aes(x = wtAvgRivDist, y = wtAvgDistNN)) +
  geom_point() + 
  geom_errorbar(aes(ymin = wtAvgDistNN - wtSeDistNN, ymax = wtAvgDistNN + wtSeDistNN)) +
  labs(x = "Distance from River Edge (m)", 
       y = "(Weighted) Mean Distance from Nearest Neighbor", 
       title = "Average Distance from Nearest Neighbors") +
  theme_bw()

brtitle <- ggplot() +
  theme_void() +
  annotate("text",
           x = 0, y = 0,
           label = "River Edge Effects per 15 meters from Edge", 
           hjust = 0.5, vjust = 0, size = 5)


binRivPlots <- ggarrange(br1, br2, br3, br4, br5, brtitle, ncol = 3, nrow = 2)
ggexport(binRivPlots, filename = "output/initialPlots/binRivPlots.pdf", height = 10, width = 15)





