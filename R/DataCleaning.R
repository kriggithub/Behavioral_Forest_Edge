#### 9/5/2025
#### Kurt Riggin
#### Data Subsetting Scripts

# read in packages
library(dplyr)
library(Hmisc)
library(stringr)



# read in working data
primateData <- read.csv("data/BehavioralData.csv", header = T)

# capitalize all IDs to standardize
primateData$ID <- toupper(primateData$ID)



### Clean data

primateData$Focal_Desc <- str_to_upper(str_trim(primateData$Focal_Desc))
primateData$ID <- str_to_upper(primateData$ID)
primateData$Activity <- str_to_upper(str_trim(primateData$Activity))
primateData$Activity[primateData$Activity == "T"] <- "L"
primateData$Activity[primateData$Activity == "OV" | primateData$Activity == ""] <- NA
primateData$NN_Desc <- str_trim(primateData$NN_Desc)
primateData$EdgeType <- str_trim(primateData$EdgeType)
primateData$Collector <- str_trim(primateData$Collector)
primateData$Rest <- as.numeric(primateData$Activity =="R")
primateData$Feed <- as.numeric(primateData$Activity =="F")
primateData$Move <- as.numeric(primateData$Activity =="L")
primateData$Social <- as.numeric(primateData$Activity =="S")


### Data Frame 1: Summarized by ID (monkey)

monkeyIdData <- primateData %>% 
  group_by(ID) %>% 
  dplyr::summarize(RestPct = mean(Activity == "R")*100,
            MovingPct = mean(Activity == "L")*100,
            FeedingPct = mean(Activity == "F")*100,
            OtherPct = mean(!(Activity %in% c("R", "L", "F")))*100,
            AvgNumNN = mean(NumNN, na.rm = T),
            AvgDistNN = mean(DistNN, na.rm = T),
            RivDist = mean(RivDist),
            AnthDist = mean(AnthDist),
            nObs = n())

# write.csv(monkeyIdData, file = "data/monkeyIdData.csv")



### Data Frame 2: Summarized within band (Anth)

anthBinData <- primateData %>% 
  mutate(bin = as.factor(ceiling(AnthDist / 15))) %>% 
  group_by(ID) %>% 
  mutate(RestPct = mean(Activity == "R")*100,
         MovingPct = mean(Activity == "L")*100,
         FeedingPct = mean(Activity == "F")*100,
         OtherPct = mean(!(Activity %in% c("R", "L", "F")))*100,
         AvgNumNN = mean(NumNN, na.rm = T),
         AvgDistNN = mean(DistNN, na.rm = T),
         RivDist = mean(RivDist),
         AnthDist = mean(AnthDist),
         nObs = n()) %>% 
  filter(row_number()==1) %>%
  ungroup() %>% 
  group_by(bin) %>% 
  dplyr::summarize(nMonkeys = n(),
                   wtAvgRestPct = weighted.mean(RestPct, nObs, na.rm = T),
                   wtSdRestPct = sqrt(wtd.var(RestPct, nObs, na.rm = T)),
                   wtSeRestPct = wtSdRestPct/sqrt(nMonkeys),
                   wtAvgMovingPct = weighted.mean(MovingPct, nObs, na.rm = T),
                   wtSdMovingPct = sqrt(wtd.var(MovingPct, nObs, na.rm = T)),
                   wtSeMovingPct = wtSdMovingPct/sqrt(nMonkeys),
                   wtAvgFeedingPct = weighted.mean(FeedingPct, nObs, na.rm = T),
                   wtSdFeedingPct = sqrt(wtd.var(FeedingPct, nObs, na.rm = T)),
                   wtSeFeedingPct = wtSdFeedingPct/sqrt(nMonkeys),
                   wtAvgOtherPct = weighted.mean(OtherPct, nObs, na.rm = T),
                   wtSdOtherPct = sqrt(wtd.var(OtherPct, nObs, na.rm = T)),
                   wtSeOtherPct = wtSdOtherPct/sqrt(nMonkeys),
                   wtAvgNumNN = weighted.mean(AvgNumNN, nObs, na.rm = T),
                   wtSdNumNN = sqrt(wtd.var(AvgNumNN, nObs, na.rm = T)),
                   wtSeNumNN = wtSdNumNN/sqrt(nMonkeys),
                   wtAvgDistNN = weighted.mean(AvgDistNN, nObs, na.rm = T),
                   wtSdDistNN = sqrt(wtd.var(AvgDistNN, nObs, na.rm = T)),
                   wtSeDistNN = wtSdDistNN/sqrt(nMonkeys),
                   wtAvgAnthDist = weighted.mean(AnthDist, nObs, na.rm = T),
                   wtSdAnthDist = sqrt(wtd.var(AnthDist, nObs, na.rm = T)),
                   wtSeAnthDist = wtSdAnthDist/sqrt(nMonkeys))


# SE = Sd / sqrt (n Monkeys in bin)



# write.csv(anthBinData, file = "data/anthBinData.csv")



### Data Frame 3: Summarized within band (Riv)

rivBinData <- primateData %>% 
  mutate(bin = as.factor(ceiling(RivDist / 15))) %>% 
  group_by(ID) %>% 
  mutate(RestPct = mean(Activity == "R")*100,
         MovingPct = mean(Activity == "L")*100,
         FeedingPct = mean(Activity == "F")*100,
         OtherPct = mean(!(Activity %in% c("R", "L", "F")))*100,
         AvgNumNN = mean(NumNN, na.rm = T),
         AvgDistNN = mean(DistNN, na.rm = T),
         RivDist = mean(RivDist),
         AnthDist = mean(AnthDist),
         nObs = n()) %>% 
  filter(row_number()==1) %>%
  ungroup() %>% 
  group_by(bin) %>% 
  dplyr::summarize(nMonkeys = n(),
                   wtAvgRestPct = weighted.mean(RestPct, nObs, na.rm = T),
                   wtSdRestPct = sqrt(wtd.var(RestPct, nObs, na.rm = T)),
                   wtSeRestPct = wtSdRestPct/sqrt(nMonkeys),
                   wtAvgMovingPct = weighted.mean(MovingPct, nObs, na.rm = T),
                   wtSdMovingPct = sqrt(wtd.var(MovingPct, nObs, na.rm = T)),
                   wtSeMovingPct = wtSdMovingPct/sqrt(nMonkeys),
                   wtAvgFeedingPct = weighted.mean(FeedingPct, nObs, na.rm = T),
                   wtSdFeedingPct = sqrt(wtd.var(FeedingPct, nObs, na.rm = T)),
                   wtSeFeedingPct = wtSdFeedingPct/sqrt(nMonkeys),
                   wtAvgOtherPct = weighted.mean(OtherPct, nObs, na.rm = T),
                   wtSdOtherPct = sqrt(wtd.var(OtherPct, nObs, na.rm = T)),
                   wtSeOtherPct = wtSdOtherPct/sqrt(nMonkeys),
                   wtAvgNumNN = weighted.mean(AvgNumNN, nObs, na.rm = T),
                   wtSdNumNN = sqrt(wtd.var(AvgNumNN, nObs, na.rm = T)),
                   wtSeNumNN = wtSdNumNN/sqrt(nMonkeys),
                   wtAvgDistNN = weighted.mean(AvgDistNN, nObs, na.rm = T),
                   wtSdDistNN = sqrt(wtd.var(AvgDistNN, nObs, na.rm = T)),
                   wtSeDistNN = wtSdDistNN/sqrt(nMonkeys),
                   wtAvgRivDist = weighted.mean(RivDist, nObs, na.rm = T),
                   wtSdRivDist = sqrt(wtd.var(RivDist, nObs, na.rm = T)),
                   wtSeRivDist = wtSdRivDist/sqrt(nMonkeys))


# write.csv(rivBinData, file = "data/rivBinData.csv")


# save.image("output/DataCleaning.RData")
