library(gtools)
library(Rmisc) 
library(dplyr) 
library(tidyr)
library(data.table)
library(ggplot2)
library(stats)
library(lubridate) 
library(RColorBrewer)
library(magrittr)
library(lme4)
library(emmeans)
library(stringr)
library(pracma)
require(purrr)
library(factoextra)
library(ggcorrplot)
library(ggbeeswarm)
library(plotly)
library(zoo)
library(ggbeeswarm)
library(plotrix)
library(lmerTest)
library(broom)
library(broom.mixed)
library(sjPlot)

#Working directory for mac
#setwd("/Volumes/huddlingvidmicro/neurophotometrics/fp_Data/processedFiles_Homecage_2024")
setwd("/Volumes/cluster/alcova/huddlingvidmicro/neurophotometrics/fp_Data/processedFiles_Homecage_2024")
getwd()

#PLOT DIRECTORY 
#"/plots_6sds_2pd/"; "/plots/"; "/plots_7sds_2pd/"
plotdir = paste(getwd(),"/plots_6sds_2pd/", sep = "") 
plotdir

######################################################################## 
########################################################################
# 
# FUNCTIONS
# 
########################################################################
########################################################################
#extract events 
extract.with.context <- function(x, colname, rows, after = 0, before = 0) {
  match.idx  <- which(x[[colname]] %in% rows)
  span       <- seq(from = -before, to = after)
  extend.idx <- c(outer(match.idx, span, `+`))
  extend.idx <- Filter(function(i) i > 0 & i <= nrow(x), extend.idx)
  extend.idx <- sort(unique(extend.idx))
  return(x[extend.idx, , drop = FALSE]) 
}

#create data list 
create_data_list <- function(extracted, rows, start_idx, span) {
  dataList <- list()
  
  #create a list of vectors for each plotgroup. This should avoid collisions due to unique plotGroups being too close to eachother. 
  for (j in 1:length(start_idx)) {
    df <- extracted
    match.idx <- start_idx[j]
    extend.idx <- c(outer(match.idx, span, `+`))
    
    # Handle special cases for the first and last start index
    if (j == 1) {
      extend.idx_2 <- extend.idx[extend.idx > 0 & extend.idx < start_idx[j + 1]]
      span_2 <- span[which(extend.idx %in% extend.idx_2)]
    } else if (j == length(start_idx)) {
      extend.idx_2 <- extend.idx[extend.idx > start_idx[j - 1] & extend.idx <= nrow(df)]
      span_2 <- span[which(extend.idx %in% extend.idx_2)]
    } else {
      extend.idx_2 <- extend.idx[extend.idx > start_idx[j - 1] & extend.idx < start_idx[j + 1]]
      span_2 <- span[which(extend.idx %in% extend.idx_2)]
    }
    
    df2 <- df[extend.idx_2, , drop = FALSE]
    df2$plotGroup <- j
    
    dataList[[j]] <- cbind(df2, span_2)
  }
  
  return(dataList)
}
#dataList <- create_data_list(extracted=extracted, rows=rows, span=span, start_idx)
#every_data = do.call(rbind, dataList)

######################################################################## 
########################################################################
# READ IN LOOKUP TABLE (lut) FOR METADATA.
# Assign colors according to factor levels 
# Assign variables 
########################################################################
########################################################################
lut <- read.csv(file = 'lookup_photom.csv',
                colClasses=c("Date"="character", "Sex" = "character",
                             "CageID" = "character", "mouseID" = "character",
                             "Ambient" = "character"))
str(lut)
#lut2 <- read.csv(file='lookup_injectionTimes.csv', colClasses=c("Injection.time.stamp"="POSIXct"))

head(lut)

thermoRate = "minute"

######################################################################## 
########################################################################
# READ IN DATA AND LOOK AT BASIC PROPERTIES 
########################################################################
########################################################################
# read in processedHuddl files that have been manually pasted into a common folder for an experiment
processed.files <- list.files(pattern = "processed")
#read in data
processed_content <-
  processed.files %>%
  lapply(read.table,
         sep=',', 
         header = TRUE,
         colClasses=c("Time.stamp"="POSIXct",  #note that Time.stamp = fp; dateTime=noldus. 
                      "noldus"="POSIXct",
                      "Time.stamp.round"="POSIXct",
                      "Activity"="numeric",
                      "lmQuotient"="numeric",
                      "deltaFslide"="numeric",
                      "timeS"="numeric",
                      "Temp"="numeric",
                      "hoboC"="numeric"
                      #"date"="character",
                      #"mouseID" = "character","cageID" = "character", "Temp" = "numeric"
                      ),
         stringsAsFactors = F,
         na.strings = c("-"))


# read file names
processed_filenames <- processed.files %>% basename() %>% as.list()
# combine file content list and file name list
all_man <- mapply(c, processed_content, processed_filenames, SIMPLIFY = FALSE)
# unlist all lists and change column name
all_man <- rbindlist(all_man, fill = T)

#rename V1
all_man <- all_man %>% dplyr::rename(Noldus_file = V1)

# convert to tibble 
combo <- as_tibble(all_man)

#Convert date to Date Format
#foo$date <- as.Date(as.character(foo$date), format="%Y-%m-%d")
combo$date <- as.Date(as.character(combo$dateTime), format="%Y-%m-%d")
lut$Date <- as.Date(as.character(lut$Date), format="%d/%m/%Y")

#Extract time 
combo$time <- strftime(combo$dateTime, format="%H:%M:%S") 

#merge combo and lut. Cols that match are Pi and Date & mouseID.
str(lut)
str(combo)
combo1 = combo %>% left_join(lut, by = c("Noldus_file")) #note that Date is startDate and date is actual date 
str(combo1)
combo1$Ambient <- as.factor(combo1$Ambient)
combo1$dateTime <- as.POSIXct(combo1$dateTime)
combo1 <- combo1 %>%
  unite("expID", c("date","Time"), remove=FALSE)

# round the data to the nearest X 
# alignRate = ".1s" # 1 minute, 2 minutes, 30 seconds, 
# combo2 = combo1 %>%
#   mutate(Time.round.1Sec = lubridate::round_date(dateTime,alignRate)) #

#get frame rate 
fps = combo1 %>% dplyr::group_by(date = format(dateTime, "%Y-%m-%d"),
                               hour = format(dateTime, "%H"),
                               minute = format(dateTime, "%M"),
                               second = format(dateTime, "%S"),
                               Noldus_file) %>%
  dplyr::summarize(count=n())
fps = mean(fps$count)
fps = round(fps,0)
fps

#Clean up combo1 data
str(combo1)
combo1 <- combo1 %>%
  mutate(mouseID = factor(mouseID),
         expID = factor(expID),
         pairedSolo = factor(pairedSolo)) %>%
  select(-any_of(c("Result.1","X","Region0G","Region0G.415","peak",
                   "deltaFslide",
                   "All.Mice.In.Dome",
                   "Injection.Human.contact.","No.contact.Human.contact.",
                   "ActivityTrim","grp","robustfit","Region0G.415.bc","Track_file","NP_file","comments",
                   "Eating.or.Drinking.Social.states.","Grooming.Other.Social.states.",
                   "Grooming.Social.states.",
                   "Time",
                   "Stationary.Social.states.","Nesting.or.Building.Nest.Social.states.","Nesting.or.Nest.Building.Social.states.")))

#plot lmQ.bc
p1 = ggplot(data=combo1, aes(x=Ambient, y=lmQ.bc, group=Ambient)) +
  geom_boxplot()
ggsave(filename = paste(plotdir,"lmQ.bc.box", ".png", sep=""), p1, width = 10, height = 10, units = 'cm')
p2 = qplot(Ambient, lmQ.bc, data = combo1, geom = "violin")
ggsave(filename = paste(plotdir,"lmQ.bc.violin", ".png", sep=""), p2, width = 10, height = 10, units = 'cm')

#plot nF.bc 
p1 = ggplot(data=combo1, aes(x=Ambient, y=nF.bc , group=Ambient)) +
  geom_boxplot()
ggsave(filename = paste(plotdir,"nF.bc.box", ".png", sep=""), p1, width = 10, height = 10, units = 'cm')
p2 = qplot(Ambient, nF.bc , data = combo1, geom = "violin")
ggsave(filename = paste(plotdir,"nF.bc .violin", ".png", sep=""), p2, width = 10, height = 10, units = 'cm')

#########
#data massaging for next steps
#########
# get rid of rows with NA in, eg, lmQuotient 
combo2 = combo1 %>%
  drop_na("lmQuotient") %>%
  drop_na("lmQ.bc")

combo2$Ambient <- as.factor(combo2$Ambient)

# make a counter per group. 
combo2 = combo2 %>% group_by(expID) %>% 
  mutate(counter = as.numeric(row_number(expID))) %>%
  ungroup()

#round off Temp, hoboC, and activity for the sake of plotting 
combo2$TempRound <- round(combo2$Temp, digits=1)
#combo2$hoboCRound <- round(combo2$hoboC, digits=1)
combo2$ActivityRound <- round(combo2$Activity, digits=1)

#########
# Make minute and second columns for downstream processing
#########
#timeS is the cumulative minute
#make a column called timeM: the cumulative minute
combo2 = combo2 %>%
  dplyr::arrange(expID, timeS) %>% 
  dplyr::mutate(timeM = timeS/60) %>%
  dplyr::mutate(timeMround = round(timeM, digits=0)) %>%
  ungroup()

#########
#str formatting
#########
# cagemate needs to be numeric 
combo2$cagemate <- as.numeric(combo2$cagemate)

########################################
# 
# Find and differentiate Activehuddle and Nesting around sleeHud
# 
########################################

#Create combo3 for testing and define parameter
combo3 <- combo2
periQuiescentRows = fps*60 #amount of rows pulled from before/after sleeHud. Originally 100
unique(combo3$huddleState)
Max <- 60 #maximum run to be kept when organizing epochs. Originally 60

#Find difference between rows and then group - Nest
#Cuts short Epochs using Max for actiHud and Nest
Nest <- subset(combo3, huddleState == "Nest")

Trial.time <- as.numeric(Nest$Time.stamp.round)
Nest_Dif <- Nest %>%
  dplyr::mutate(Time.Elapsed = c(0,diff(Trial.time)))
Nest_Dif$Time.Elapsed <- replace(dplyr::lag(Nest_Dif$Time.Elapsed, n=1), Nest_Dif$Time.Elapsed > Max, NA)
Nest_Dif$Time.Elapsed <- dplyr::lead(Nest_Dif$Time.Elapsed, n=1)

Nest_Grouped <- Nest_Dif %>% 
  dplyr::mutate(group = cumsum(lag(is.na(Time.Elapsed), default = TRUE)) )
Nest_Grouped = Nest_Grouped %>%
  tidyr::unite("groupState", c(huddleState,group), remove=FALSE) 

#Find difference between rows and then group - actiHud
#Cuts short Epochs using Max for actiHud and Nest
actiHud <- subset(combo3, huddleState == "actiHud")

Trial.time <- as.numeric(actiHud$Time.stamp.round)
actiHud_Dif <- actiHud %>%
  dplyr::mutate(Time.Elapsed = c(0,diff(Trial.time)))
actiHud_Dif$Time.Elapsed <- replace(dplyr::lag(actiHud_Dif$Time.Elapsed, n=1), actiHud_Dif$Time.Elapsed > Max, NA)
actiHud_Dif$Time.Elapsed <- dplyr::lead(actiHud_Dif$Time.Elapsed, n=1)

actiHud_Grouped <- actiHud_Dif %>% 
  dplyr::mutate(group = cumsum(lag(is.na(Time.Elapsed), default = TRUE)) )
actiHud_Grouped = actiHud_Grouped %>%
  tidyr::unite("groupState", c(huddleState,group), remove=FALSE)

#cbind subsets back to combo3
combo3 = combo3 %>%
  dplyr::left_join(Nest_Grouped, by='Time.stamp.round', suffix = c("", ".y"), multiple = "first") %>%
  dplyr::select(-ends_with(".y")) %>%
  
  dplyr::left_join(actiHud_Grouped, by='Time.stamp.round', suffix = c("", ".y"), multiple = "first") %>%
  dplyr::mutate(groupState = coalesce(groupState, groupState.y)) %>% 
  dplyr::mutate(Time.Elapsed = coalesce(Time.Elapsed, Time.Elapsed.y)) %>% 
  dplyr::select(-ends_with(".y"))

#Generate rlid to use for incorporating full epochs later
combo3 <- combo3 %>% 
  dplyr::mutate(rlid = data.table::rleid(groupState))
combo3 <- within(combo3, rlid[is.na(groupState)] <- NA)

#Create prequiescent
combo3$preQuiescent <- 1

#Find the indices of "start" occurrences for sleeHud_start 
start_indices <- which(combo3$bStartStop == "sleeHud_start" | combo3$bStartStop == "sleeSol_start")

#Iterate over the start indices and fill in sleeHud 
for (i in start_indices) {
  range_start <- max(i - periQuiescentRows, 1)
  range_end <- i - 1
  combo3$preQuiescent[range_start:range_end] <- combo3$huddleState[range_start:range_end]
}

#convert 1s to NAs
combo3 = combo3 %>%
  dplyr::mutate(preQuiescent = replace(preQuiescent, preQuiescent == 1, NA))

#Create postquiescent
combo3$postQuiescent <- 1

#Find the indices of "stop" occurrences for sleeHud_stop
stop_indices <- which(combo3$bStartStop == "sleeHud_stop"| combo3$bStartStop == "sleeSol_stop")

#Iterate over the stop indices and fill in pre sleeHud
for (i in stop_indices) {
  range_start <- min(i + 1, nrow(combo3))
  range_end <- min(i + periQuiescentRows, nrow(combo3))
  if (range_end > length(combo3$postQuiescent)) {
    range_end <- length(combo3$postQuiescent)
  }
  combo3$postQuiescent[range_start:range_end] <- combo3$huddleState[range_start:range_end]
}

#convert 1s to NAs
combo3 = combo3 %>%
  dplyr::mutate(postQuiescent = replace(postQuiescent, postQuiescent == 1, NA))

#Use rlid to complete epochs
combo3.1 = combo3 %>%
  drop_na(rlid) %>%
  dplyr::group_by(rlid) %>%
  #fill(preQuiescent) %>%
  fill(preQuiescent, .direction = "up") %>% 
  fill(postQuiescent, .direction = "down") %>%
  ungroup()

combo4 = combo3 %>%
  dplyr::left_join(combo3.1, by='Time.stamp.round', suffix = c("", ".y"), multiple = "first") %>%
  dplyr::select(-ends_with(".y"))

#Define both and neither
combo4 = combo4 %>%
  dplyr::mutate(Both = case_when(
    !is.na(preQuiescent) & !is.na(postQuiescent) & (grepl("Nest", huddleState) | grepl("actiHud", huddleState)) ~ paste("both_",preQuiescent),
    is.na(preQuiescent) & is.na(postQuiescent) & (grepl("Nest", huddleState) | grepl("actiHud", huddleState)) ~ paste("neither_",huddleState),
    TRUE ~ NA_character_))

#Define pre and post sleeHud
combo4 = combo4 %>%
  dplyr::mutate(huddleState2 = case_when(
    !is.na(combo4$preQuiescent) & is.na(Both) & (grepl("Nest", huddleState) | grepl("actiHud", huddleState)) ~ paste("pre_", preQuiescent),
    !is.na(combo4$postQuiescent) & is.na(Both) & (grepl("Nest", huddleState) | grepl("actiHud", huddleState)) ~ paste("post_", postQuiescent),
    !is.na(Both) ~ Both,
    TRUE ~ NA_character_))

#Limit number of columns
combo4 <- combo4 %>%
  dplyr::select(-group, -preQuiescent, -postQuiescent, -Both, -rlid)

#Reassign to combo2 for now to test, will go through and change later
combo2 <- combo4

unique(combo2$huddleState2)
############################
#End workspace here
############################

#order solo before paired
combo2$pairedSolo <- factor(combo2$pairedSolo, levels = c("solo","paired"))

########################################
# Plot calcium per behavioral state using PEAKS
# PART 1: lmQ.bc
# this may be the best metric for comparing between solo and paired
########################################
##########
# define parameters
##########
sds = 6
peakdistance = fps*2

###################### 
#Z-score lmQ.bc 
# Note to self: change the name to lmQ.bc.Z. Systematically go through and replace as needed
# This is important because if you rerun this line, it will re-Zscore the the Zscored data
######################
combo2 = combo2 %>%
  dplyr::group_by(expID) %>%
  dplyr::mutate(lmQ.bc.Z = (lmQ.bc - mean(lmQ.bc))/sd(lmQ.bc)) %>%
  ungroup()

###################### 
#"Maximum normalization" : lmQ.bc.MM 
######################
combo2 = combo2 %>%
  dplyr::group_by(mouseID, pairedSolo) %>%
  dplyr::mutate(lmQ.bc.MM = (lmQ.bc / (max(lmQ.bc)))) %>%
  ungroup()

###################### 
# check on a subset first: can do this for either lmQ.bc.Z or lmQ.bc.MM
###################### 
#plot a subset of the data 
head = combo2[15000:60000, ]
g = ggplot(head, aes(x=timeS, y = lmQ.bc.Z)) + geom_line() #Time.round.1Sec #Time.round1sec
g2 = g + 
  scale_x_continuous(breaks = round(seq(min(head$timeS), max(head$timeS), by = 10),digits = 0)) + 
  theme_classic()
g2
#ggplotly(g2)

h = ggplot(head, aes(x=timeS, y = lmQuotient)) + geom_line() #Time.round.1Sec #Time.round1sec
h2 = h + 
  scale_x_continuous(breaks = round(seq(min(head$timeS), max(head$timeS), by = 10),digits = 0)) + 
  theme_classic()
h2
#ggplotly(h2)

#examine standard deviation in solo and paired animals 
hsem <- summarySE(combo2, measurevar="lmQ.bc", groupvars=c("mouseID","pairedSolo","Ambient"), na.rm=TRUE) 
ggplot(hsem, aes(x=pairedSolo, y=sd)) + 
  facet_grid(cols = vars(Ambient))+
  geom_point()
#Zscore 
hsem <- summarySE(combo2, measurevar="lmQ.bc.Z", groupvars=c("mouseID","pairedSolo","Ambient"), na.rm=TRUE) 
ggplot(hsem, aes(x=pairedSolo, y=sd)) + 
  facet_grid(cols = vars(Ambient))+
  geom_point()

############################################  
# FIND PEAKS  
############################################ 
# first a subset of data (aka head)
peaks <- findpeaks(head$lmQ.bc.Z,
                   minpeakdistance = peakdistance,
                   #minpeakheight = quantile(head$deltaFslide, na.rm=TRUE)[[4]]*1.20,
                   minpeakheight = sd(head$lmQ.bc.Z)*sds)
head$new_col <- rep(0, nrow(head))
head$new_col <- replace(head$new_col, peaks[,2], values=1) #the second the position/index is the maximum is reached

pHead = ggplot(head, aes(x=timeS, y = lmQ.bc.Z)) + 
  geom_line() + 
  geom_point(data= head[head$new_col ==1, ], color = "red")
pHead
#ggplotly(pHead)

# On whole dataset 
# First separate by mouse ID, get SD for that data, 
# then find peaks for each experiment based on mouse's SD.
# This conserves expID and mouseID in final output
res <- data.frame()
for(m in unique(combo2$mouseID)) {
  sub = combo2 %>% dplyr::filter(mouseID == m)
  sub_SD = sd(sub$lmQ.bc.Z)
  for(e in unique(sub$expID)) {
    sub2 = sub %>% dplyr::filter(expID == e)
    peaks <- findpeaks(sub2$lmQ.bc.Z,
                       minpeakdistance = peakdistance,
                       minpeakheight = sub_SD * sds) 
    peaks <- data.frame(peak_height = peaks[, 1], peak_index = peaks[, 2])
    peaks1 <- peaks %>% 
      mutate(mouseID = m,
             expID = e)
    res <- rbind(res, peaks1)
  }
}

# #this calculates SD based on mouseID and solo vs paired
# res <- data.frame()
# for(m in unique(combo2$mouseID)) {
#   for(sp in unique(combo2$pairedSolo)) {
#   sub = combo2 %>% dplyr::filter(mouseID == m & pairedSolo == sp)
#   sub_SD = sd(sub$lmQ.bc.Z)
#   for(e in unique(sub$expID)) {
#     sub2 = sub %>% dplyr::filter(expID == e)
#     peaks <- findpeaks(sub2$lmQ.bc.Z,
#                        minpeakdistance = peakdistance,
#                        minpeakheight = sub_SD * sds)
#     peaks <- data.frame(peak_height = peaks[, 1], peak_index = peaks[, 2])
#     peaks1 <- peaks %>%
#       mutate(mouseID = m,
#              expID = e)
#     res <- rbind(res, peaks1)
#   }
#   }
# }
# #this calculates SD based experimentID
# res <- data.frame()
# for(e in unique(combo2$expID)) {
#   sub = combo2 %>% dplyr::filter(expID == e)
#   sub_SD = sd(sub$lmQ.bc.Z)
#   peaks <- findpeaks(sub$lmQ.bc.Z,
#                        minpeakdistance = peakdistance,
#                        minpeakheight = sub_SD * sds)
#   peaks <- data.frame(peak_height = peaks[, 1], peak_index = peaks[, 2])
#   peaks1 <- peaks %>%
#     mutate(expID = e)
#   res <- rbind(res, peaks1)
#   }

# join res back to combo2, then get rid of NAs
res$peak = res$peak_index
combo2peaks = combo2 %>% 
  left_join(res, by=c("mouseID"="mouseID", "expID" = "expID", "counter"="peak")) #"mouseID"="mouseID", 

#need a binary for rows containing a peak
combo2peaks = combo2peaks %>% 
  dplyr::mutate(peakBinary = case_when(is.na(peak_index) ~ 0,
                                       !is.na(peak_index) ~ 1))

#Double check that lmQ.bc.Z matches peakHeight
peaks_sub <- combo2peaks %>% dplyr::filter(peak_index != "NA")
peakSE <- summarySE(combo2, measurevar = "lmQ.bc.Z", groupvars = c("expID", "mouseID", "Ambient", "pairedSolo"))

#Create plotting parameters
ymin = min(combo2peaks$lmQ.bc.Z)
ymax = max(combo2peaks$lmQ.bc.Z)

#plot traces for each ambient temperature, SUBSETTED PER EXPERIMENT
for(i in unique(combo2peaks$pairedSolo)){
  sub <- combo2peaks %>% dplyr::filter(pairedSolo == i)
  p3 = ggplot(sub, aes(x=timeS, y = lmQ.bc.Z)) + 
    facet_grid(cols=vars(expID), rows = vars(mouseID)) +
    geom_line() + 
    geom_point(aes(x = timeS, y = peak_height), color = "red") + 
    ylim(ymin, ymax)  +
    theme_minimal()
  #p3
  ggsave(filename = paste(plotdir,"peaks/", "peaks","_lmQ.bc.Z", "_peakDistance",peakdistance,"_",i,"_SDs", sds, ".pdf", sep=""), 
         p3,width = 25, height = 9, units = 'cm')
}

#drop NAs from combo2peaks and plot just the peaks 
combo2peaks2 = combo2peaks %>% drop_na("peak_index")

###############
# Paired and Solo data plots  
###############
#summarize and plot peak values: huddleState
for(i in levels(combo2peaks2[["pairedSolo"]])) {
  df0 = combo2peaks2[combo2peaks2$pairedSolo == i, ]
  df0.1 <- summarySE(df0, measurevar="lmQ.bc.Z", groupvars=c("huddleState", "mouseID"), na.rm=TRUE) 
  
  #huddlestate
  pp = ggplot(df0, aes(x=huddleState, y=lmQ.bc.Z)) + #color = mouseID
    geom_beeswarm(color = "grey60", size=0.5,  alpha=0.5) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1)) +
    theme(legend.position="none") + 
    theme(axis.title.x = element_blank())
  print(pp)
  ggsave(filename = paste(plotdir,"peakValues/", "peakValues_by_huddleState_",i,"_","lmQ.bc.Z",
                          "_peakDistance",peakdistance,"_SDs", sds, ".pdf", sep=""),pp,
         width = 5, height = 4, units = 'cm')
  
  # Ambient 
  hsem <- summarySE(df0, measurevar="lmQ.bc.Z", groupvars=c("Ambient"), na.rm=TRUE) 
  ggplot(df0, aes(x=Ambient, y=lmQ.bc.Z)) +
    geom_beeswarm(color = "grey60", alpha = 0.7, size = 1.5, corral.width = 1.5) +
    geom_errorbar(data = hsem, aes(ymin=lmQ.bc.Z-se, ymax=lmQ.bc.Z+se), colour="black", width=.4) +
    theme_classic() +
    theme(axis.title.x = element_blank())
  ggsave(filename = paste(plotdir,"peakValues/","peakValues_by_Ambient_",i,"_","lmQ.bc.Z",
                          "_peakDistance",peakdistance,"_SDs", sds, ".pdf", sep=""),
         width = 4, height = 4, units = 'cm')
  
  #cumulative distribution function
  ggplot(df0, aes(x = lmQ.bc.Z, col = huddleState))+
    stat_ecdf(alpha=0.6) + 
    theme_classic()
  ggsave(filename = paste(plotdir,"peakValues/","CDF","_lmQ.bc.Z_",i,"_",
                          "_peakDistance",peakdistance,"_SDs", sds,"_", ".pdf", sep=""),
         width = 8, height = 5, units = 'cm')
  
}

###############
# Solo vs paired data 
###############
#summarize and plot peak values: solo vs paired 
hsem <- summarySE(combo2peaks2, 
                  measurevar="lmQ.bc.Z", groupvars=c("pairedSolo"), na.rm=TRUE) 
ggplot(combo2peaks2, aes(x=pairedSolo, y=lmQ.bc.Z)) + 
  geom_beeswarm(color = "grey60", alpha = 0.7, size = 1.5, corral.width = 1.5) +
  #geom_errorbar(data = hsem, aes(ymin=lmQ.bc.Z-se, ymax=lmQ.bc.Z+se), colour="black", width=.2, size = 0.5,position=pd) +
  stat_summary(fun.data = "mean_se", 
               geom = "errorbar", #"pointrange"
               width = 0.4,
               size = 0.3,
               lwd = 0.5,
               #position = position_nudge(x= c(0.2,0.2))
  ) +
  theme_classic() +
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"peakValues/","peakValues_by_soloPaired_","lmQ.bc.Z","_peakDistance",peakdistance,"_SDs", sds, "allData", ".pdf", sep=""),
       width = 4, height = 4, units = 'cm')

######################################################################## 
########################################################################
# Spaghetti plot it!! Peak counts and peak amp (combo2peaks2)
########################################################################
########################################################################

################## 
# spag plot for peak amplitude (combo2peaks2 data)
################## 
peakcount <- combo2peaks2 %>%
  group_by(mouseID, pairedSolo, Ambient)%>%
  summarise(counts = n())

#reference file to include all mouseIDs, file names and temps
peakcount2 <- combo2peaks %>%
  group_by(mouseID, pairedSolo, Ambient, File.Name)%>%
  summarise(ignore = n())

#Join the two and convert NA to 0
peakcount <- peakcount2 %>%
  left_join(peakcount, by=c("mouseID"="mouseID", "pairedSolo" = "pairedSolo", "Ambient"="Ambient"))
peakcount$counts[is.na(peakcount$counts)] <- 0

peakamp <- combo2peaks2 %>%
  group_by(mouseID, pairedSolo, Ambient, File.Name) %>%
  summarise(meanpeakamp = mean(lmQ.bc.Z)) %>%
  ungroup()

#Join to peakcount2 for NA-->0 
peakamp <- peakcount2 %>%
  left_join(peakamp, by=c("mouseID"="mouseID", "pairedSolo" = "pairedSolo", "Ambient"="Ambient", "File.Name"))
peakamp$meanpeakamp[is.na(peakamp$meanpeakamp)] <- 0

#plot
pos <- position_nudge(x=c(0.2,0.2,0.2))
peakampSE <- summarySE(peakamp, measurevar = "meanpeakamp", groupvars = c("pairedSolo","Ambient"), na.rm = TRUE)
p = ggplot(data=peakampSE, aes(x = pairedSolo, y = meanpeakamp)) + 
  facet_grid(cols = vars(Ambient)) + 
  geom_line(data = peakamp, aes(x=pairedSolo, y= meanpeakamp, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
  geom_point(data= peakamp, aes(x=pairedSolo, y = meanpeakamp, color= mouseID), width = 0.1, alpha=0.7) +
  geom_errorbar(aes(ymin=meanpeakamp-se, ymax=meanpeakamp +se), width=.3, color="black", position = pos)+
  theme_classic()
p
ggsave(filename = paste(plotdir,"summary/", "summary_peak_amp_mouseID_spag_combo2peaks2", peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""),
       p, width = 13, height = 7, units = 'cm') 

#plot peak amp by mouseIDs: SoloPaired
hsem4.pa <- summarySE(peakamp, measurevar = "meanpeakamp", groupvars = c("mouseID", "pairedSolo"), na.rm = TRUE)
ggplot(data= peakamp, aes(x=pairedSolo, y=meanpeakamp)) + 
  geom_point(data= hsem4.pa, aes(x=pairedSolo, y = meanpeakamp, color= mouseID), width = 0.1, alpha=0.7) +
  geom_line(data = hsem4.pa, aes(x=pairedSolo, y= meanpeakamp, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
  stat_summary(fun.data = "mean_se", 
               geom = "pointrange",
               width = 0.01,
               size = 0.05,
               lwd = 0.5,
               position = position_nudge(x= c(-0.2,0.2))) +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + 
  theme(legend.position="none")
ggsave(filename = paste(plotdir,"summary/", "summary_peak_amp_social_spag_", peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""),
       width = 4, height = 4, units = 'cm')

#plot peak amp by ambient for each pairedSolo
for (i in levels(peakamp[["pairedSolo"]])) {
  df = peakamp[peakamp$pairedSolo == i, ]
  dfSmry =  summarySE(df[df$pairedSolo == i, ],
                      measurevar = "meanpeakamp", groupvars = c("mouseID", "Ambient"), na.rm = TRUE)
  
  p2 = ggplot(data= df, aes(x=Ambient, y=meanpeakamp)) + 
    geom_point(data= dfSmry, aes(x=Ambient, y = meanpeakamp, color= mouseID), width = 0.1, alpha=0.7) +
    geom_line(data = dfSmry, aes(x=Ambient, y= meanpeakamp, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 #width = 0.01,
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(0.2,0.2,0.2))) +
    theme_classic() + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none")
  print(p2)
  ggsave(filename = paste(plotdir,"summary/", "summary_peak_amp_ambient_spag_", i, "_",
                          peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""), p2,
         width = 4, height = 4, units = 'cm')
}

#plot peak amp by pairedSolo for each ambient
for (i in levels(peakamp[["Ambient"]])) {
  df = peakamp[peakamp$Ambient == i, ]
  dfSmry =  summarySE(df[df$Ambient == i, ],
                      measurevar = "meanpeakamp", groupvars = c("mouseID", "pairedSolo"), na.rm = TRUE)
  
  p2.2 = ggplot(data= df, aes(x=pairedSolo, y=meanpeakamp)) + 
    geom_point(data= dfSmry, aes(x=pairedSolo, y = meanpeakamp, color= mouseID), width = 0.1, alpha=0.7) +
    geom_line(data = dfSmry, aes(x=pairedSolo, y= meanpeakamp, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 #width = 0.01,
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(0.2))) +
    theme_classic() + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none")
  print(p2.2)
  ggsave(filename = paste(plotdir,"summary/", "summary_peak_amp_ambient_x_pairedSolo_spag_", i, "_",
                          peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""), p2.2,
         width = 4, height = 4, units = 'cm')
}

# plot peak amp by ambient: ALL
hsem4.pa <- summarySE(peakamp, measurevar = "meanpeakamp", groupvars = c("mouseID", "Ambient"), na.rm = TRUE)
ggplot(data= peakamp, aes(x=Ambient, y=meanpeakamp)) + 
  geom_point(data= hsem4.pa, aes(x=Ambient, y = meanpeakamp, color= mouseID), width = 0.1, alpha=0.7) +
  geom_line(data = hsem4.pa, aes(x=Ambient, y= meanpeakamp, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
  stat_summary(fun.data = "mean_se", 
               geom = "pointrange",
               width = 0.01,
               size = 0.05,
               lwd = 0.5,
               position = position_nudge(x= c(0.2,0.2,0.2))) +
  theme_classic() + 
  theme(axis.title.x = element_blank()) + 
  theme(legend.position="none")
ggsave(filename = paste(plotdir,"summary/", "summary_peak_amp_ambient_spag_all_", peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""),
       width = 4, height = 4, units = 'cm')

################## 
# peak amplitude per huddle state: pairedSolo
################## 
for(i in levels(combo2peaks2[["pairedSolo"]])) {
  peakamp.paired = combo2peaks2[combo2peaks2$pairedSolo == i,]
  
  peakamp.paired <- peakamp.paired %>%
    group_by(mouseID, pairedSolo, Ambient, huddleState, File.Name)%>%
    dplyr::summarise(
      meanAMP = mean(lmQ.bc.Z),
      counts = n()) %>%
    dplyr::mutate(weight = counts / (sum(.$counts))) %>%
    dplyr::mutate(meanAMPw = meanAMP*weight) %>%
    ungroup()
  
  #reference file to include all mouseIDs, file names and temps
  peakcount.paired2 = combo2peaks[combo2peaks$pairedSolo == i,]
  
  peakcount.paired3 <- peakcount.paired2 %>%
    group_by(mouseID, Ambient, File.Name, huddleState, pairedSolo)%>%
    summarise(ignore = n())
  
  #Join the two
  peakamp.paired <- peakcount.paired3 %>%
    left_join(peakamp.paired, by=c("mouseID"="mouseID", "Ambient"="Ambient", 
                                     "huddleState" = "huddleState", "File.Name" = "File.Name", "pairedSolo" = "pairedSolo"))
  peakamp.paired$meanAMP[is.na(peakamp.paired$counts)] <- 0
  peakamp.paired$meanAMPw[is.na(peakamp.paired$counts)] <- 0
  
  # plot per-mouse meanAMP
  peakcount.pairedSE <- summarySE(peakamp.paired, 
                                  measurevar="meanAMPw", groupvars=c("huddleState", "mouseID"), na.rm=TRUE)
  p3 = ggplot(data = peakcount.pairedSE, aes(x=huddleState, y=meanAMPw)) + 
    geom_point(data= peakcount.pairedSE, aes(x=huddleState, y = meanAMPw, color= mouseID), size = 1, alpha=0.5) +
    geom_line(data = peakcount.pairedSE, aes(x=huddleState, y= meanAMPw, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(0))) +
    theme_classic() + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none") + 
    labs(y ="Mean amp") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  print(p3)
  ggsave(filename = paste(plotdir,"summary/", "peakamp_by_huddleState_", i, "_","lmQ.bc.Z",
                          peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""), p3,
         width = 6, height = 4, units = 'cm')
  
  mod = lmer(meanAMPw ~ 0 + huddleState + (1|mouseID), data = peakcount.pairedSE)
  modPs = summary(mod)$coefficients
  write.csv(modPs,  file = paste(plotdir,"summary/", "peakamp_by_huddleState_", i, "_","lmQ.bc.Z",
                                 peakdistance, "_peakdistance", "_sds", sds, ".csv", sep = ""))
  
  # plot per-mouse, per-ambient meanAMP
  peakcount.pairedSE <- summarySE(peakamp.paired, 
                                  measurevar="meanAMPw", groupvars=c("huddleState", "mouseID", "Ambient"), na.rm=TRUE)
  p4 = ggplot(data = peakcount.pairedSE, aes(x=huddleState, y=meanAMPw)) + 
    geom_point(data= peakcount.pairedSE, aes(x=huddleState, y = meanAMPw, color= mouseID), size = 1, alpha=0.5) +
    geom_line(data = peakcount.pairedSE, aes(x=huddleState, y= meanAMPw, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(0))) +
    theme_classic() + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none") + 
    labs(y ="Mean amp") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_grid(cols = vars(Ambient))
  print(p4)
  ggsave(filename = paste(plotdir,"summary/", "peakamp_by_huddleState_by_AMBIENT_", i, "_","lmQ.bc.Z",
                          peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""), p4,
         width = 12, height = 7, units = 'cm')
  
  mod = lmer(meanAMPw ~ 0 + huddleState + Ambient + (1|mouseID), data = peakcount.pairedSE)
  summary(mod)
  modPs = summary(mod)$coefficients
  write.csv(modPs,  file = paste(plotdir,"summary/", "peakamp_by_huddleState_by_AMBIENT_", i, "_","lmQ.bc.Z",
                                 peakdistance, "_peakdistance", "_sds", sds, ".csv", sep = ""))
  
}

################## 
# peak counts per huddle state
################## 
for(i in levels(combo2peaks2[["pairedSolo"]])) {
  peakcount.paired = combo2peaks2[combo2peaks2$pairedSolo == i,]
  
  peakcount.paired <- peakcount.paired %>%
    group_by(mouseID, pairedSolo, Ambient, huddleState, File.Name)%>%
    summarise(counts = n())
  
  #reference file to include all mouseIDs, file names and temps
  peakcount.paired2 = combo2peaks[combo2peaks$pairedSolo == i,]
  
  peakcount.paired3 <- peakcount.paired2 %>%
    group_by(mouseID, Ambient, File.Name, huddleState, pairedSolo)%>%
    summarise(ignore = n())
  
  #Join the two
  peakcount.paired <- peakcount.paired3 %>%
    left_join(peakcount.paired, by=c("mouseID"="mouseID", "Ambient"="Ambient", 
                                     "huddleState" = "huddleState", "File.Name" = "File.Name", "pairedSolo" = "pairedSolo"))
  peakcount.paired$counts[is.na(peakcount.paired$counts)] <- 0
  
  # plot per-mouse counts
  peakcount.pairedSE <- summarySE(peakcount.paired, 
                                  measurevar="counts", groupvars=c("huddleState", "mouseID"), na.rm=TRUE)
  p3 = ggplot(data = peakcount.pairedSE, aes(x=huddleState, y=counts)) + 
    geom_point(data= peakcount.pairedSE, aes(x=huddleState, y = counts, color= mouseID), size = 1, alpha=0.5) +
    geom_line(data = peakcount.pairedSE, aes(x=huddleState, y= counts, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(0))) +
    theme_classic() + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none") + 
    labs(y ="Mean peak count") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  print(p3)
  ggsave(filename = paste(plotdir,"summary/", "peakcounts_by_huddleState_", i, "_","lmQ.bc.Z",
                          peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""), p3,
         width = 6, height = 4, units = 'cm')
  #get estimates for each effect: https://stackoverflow.com/questions/68424747/how-to-test-if-linear-mixed-effects-model-lmer-is-greater-than-1-in-r 
  mod = lmer(counts ~ 0 + huddleState + (1|mouseID), data = peakcount.pairedSE)
  modPs = summary(mod)$coefficients
  write.csv(modPs,  file = paste(plotdir,"summary/", "peakcounts_by_huddleState_", i, "_","lmQ.bc.Z",
                                 peakdistance, "_peakdistance", "_sds", sds, ".csv", sep = ""))
  
  # plot per-mouse counts
  peakcount.pairedSE <- summarySE(peakcount.paired, 
                                  measurevar="counts", groupvars=c("huddleState", "mouseID", "Ambient"), na.rm=TRUE)
  p4.2 = ggplot(data = peakcount.pairedSE, aes(x=huddleState, y=counts)) + 
    geom_point(data= peakcount.pairedSE, aes(x=huddleState, y = counts, color= mouseID), size = 1, alpha=0.5) +
    geom_line(data = peakcount.pairedSE, aes(x=huddleState, y= counts, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(0))) +
    theme_classic() + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none") + 
    labs(y ="Mean peak count") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    facet_grid(cols = vars(Ambient))
  print(p4.2)
  ggsave(filename = paste(plotdir,"summary/", "peakcounts_by_huddleState_by_AMBIENT_", i, "_","lmQ.bc.Z",
                          peakdistance, "_peakdistance", "_sds", sds, ".pdf", sep = ""), p4.2,
         width = 12, height = 7, units = 'cm')
  #get estimates for each effect: https://stackoverflow.com/questions/68424747/how-to-test-if-linear-mixed-effects-model-lmer-is-greater-than-1-in-r 
  mod = lmer(counts ~ 0 + huddleState + (1|mouseID), data = peakcount.pairedSE)
  modPs = summary(mod)$coefficients
  write.csv(modPs,  file = paste(plotdir,"summary/", "peakcounts_by_huddleState_by_AMBIENT_", i, "_","lmQ.bc.Z",
                                 peakdistance, "_peakdistance", "_sds", sds, ".csv", sep = ""))
}

################## 
# Calculate and plot peak frequency
################## 
expdurations <- combo2peaks %>%
  group_by(mouseID, File.Name, pairedSolo, Ambient, time)%>%
  summarise(duration = n()) %>%
  summarise(duration=sum(duration/10)) %>%
  ungroup()

expduration.pc <- expdurations %>%
  left_join(peakcount, by=c("File.Name" = "File.Name", "pairedSolo" = "pairedSolo", "Ambient"="Ambient", "mouseID"="mouseID"))

expduration.pc2 <- expduration.pc %>%
  group_by(mouseID, pairedSolo, Ambient)%>%
  summarise(duration = mean(duration), counts = mean(counts)) %>%
  ungroup()

expduration.pc1 <- expduration.pc %>%
  mutate(peakfreq = (counts/duration)*1000)

#plot peak frequency per mouseID and ambient
expduration.se <- summarySE(expduration.pc1, measurevar = "peakfreq", groupvars = c("pairedSolo","Ambient"), na.rm = TRUE)
p = ggplot(data=expduration.pc1, aes(x = pairedSolo, y = peakfreq)) + 
  facet_grid(cols = vars(Ambient)) + 
  geom_line(data = expduration.pc1, aes(x=pairedSolo, y= peakfreq, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
  geom_jitter(data= expduration.pc1, aes(x=pairedSolo, y = peakfreq, color= mouseID), width = 0.1, alpha=0.7) +
  geom_errorbar(data=expduration.se, aes(ymin=peakfreq -se, ymax=peakfreq +se), width=.3, color="black", position = pos)+
  theme_classic() +
  labs(x = "Assay", y= "Peak Frequency (arb)")
p
scale_x_discrete(guide = guide_axis(angle = 45)) 
ggsave(filename = paste(plotdir,"summary/", "summary_peakfreq_mouseID_combo2peaks2", "_peakdistance", peakdistance, "_sds", sds, ".pdf", sep = ""),
       p, width = 13, height = 7, units = 'cm')

#plot peak freq by mouseIDs: ambient
expduration.se <- summarySE(expduration.pc1, measurevar = "peakfreq", groupvars = c("mouseID", "Ambient"), na.rm = TRUE)
p = ggplot(data=,expduration.pc1,aes(x = Ambient, y = peakfreq)) + 
  geom_line(data = expduration.se, aes(x=Ambient, y= peakfreq, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
  geom_point(data= expduration.se, aes(x=Ambient, y = peakfreq, color= mouseID), width = 0.1, alpha=0.7) +
  stat_summary(fun.data = "mean_se", 
               geom = "pointrange",
               #width = 0.01,
               size = 0.05,
               lwd = 0.5,
               position = position_nudge(x= c(0.2,0.2,0.2))) +
  theme_classic() +
  theme(axis.title.x = element_blank()) + 
  theme(legend.position="none") +
  labs(x = "Assay", y= "Peak frequency (au)")
p
ggsave(filename = paste(plotdir,"summary/", "summary_peakfreq_ambient_all", "_peakdistance", peakdistance, "_sds", sds, ".pdf", sep = ""),
       p, width = 4, height = 4, units = 'cm')

# plot peak frequency per ambient for each pairedSolo
for (i in levels(expduration.pc1[["pairedSolo"]])) {
  df2 = expduration.pc1[expduration.pc1$pairedSolo == i, ]
  dfSmry2 = summarySE(df2, measurevar = "peakfreq", groupvars = c("mouseID", "Ambient"), na.rm = TRUE)
  
  p4 = ggplot(data=df2, aes(x = Ambient, y = peakfreq)) + 
    geom_line(data = dfSmry2, aes(x=Ambient, y= peakfreq, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
    geom_point(data= dfSmry2, aes(x=Ambient, y = peakfreq, color= mouseID), alpha=0.7) +
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 #width = 0.01,
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(0.2,0.2,0.2))) +
    theme_classic() + 
    labs(y= "Peak Frequency (arb)") + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none")
  print(p4)
  ggsave(filename = paste(plotdir,"summary/", "summary_peakfreq_ambient_",i, "_peakdistance", peakdistance, "_sds", sds, ".pdf", sep = ""),p4,
         width = 4, height = 4, units = 'cm')
  
  meanPF = mean(df2$peakfreq)
  mod = lmer(peakfreq ~ Ambient + (1|mouseID), data = df2)
  summary(mod)
  modPs = summary(mod)$coefficients
  #write.csv(modPs,  file = paste(plotdir,"summary/", "summary_peakfreq_ambient_",i, "_peakdistance", peakdistance, "_sds", sds,".csv",sep = ""))
  
  modTK = emmeans(mod, list(pairwise ~ Ambient), adjust = "tukey")
  write.csv(tidy(modTK$`pairwise differences of Ambient`),
            paste(plotdir,"summary/", "summary_peakfreq_ambient_",i, "_peakdistance", peakdistance, "_sds", sds,
                  ".csv",sep = ""))
}

# plot peak frequency per pariedSolo for each ambient
for (i in levels(expduration.pc1[["Ambient"]])) {
  df2 = expduration.pc1[expduration.pc1$Ambient == i, ]
  dfSmry2 = summarySE(df2, measurevar = "peakfreq", groupvars = c("mouseID", "pairedSolo"), na.rm = TRUE)
  
  p4.2 = ggplot(data=df2, aes(x = pairedSolo, y = peakfreq)) + 
    geom_line(data = dfSmry2, aes(x=pairedSolo, y= peakfreq, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
    geom_point(data= dfSmry2, aes(x=pairedSolo, y = peakfreq, color= mouseID), alpha=0.7) +
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 #width = 0.01,
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(-0.2, 0.2))) +
    theme_classic() + 
    labs(y= "Peak Frequency (arb)") + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none")
  print(p4.2)
  ggsave(filename = paste(plotdir,"summary/", "summary_peakfreq_pairedSolo_x_Ambient_",i, "_peakdistance", peakdistance, "_sds", sds, ".pdf", sep = ""),p4.2,
         width = 4, height = 4, units = 'cm')
}

#plot peak freq by mouseIDs: pairedSolo
expduration.se <- summarySE(expduration.pc1, measurevar = "peakfreq", groupvars = c("mouseID", "pairedSolo"), na.rm = TRUE)
p = ggplot(data=,expduration.pc1,aes(x = pairedSolo, y = peakfreq)) + 
  geom_line(data = expduration.se, aes(x=pairedSolo, y= peakfreq, group=mouseID, color=mouseID), alpha=0.25, color="gray") +
  geom_point(data= expduration.se, aes(x=pairedSolo, y = peakfreq, color= mouseID), width = 0.1, alpha=0.7) +
  stat_summary(fun.data = "mean_se", 
               geom = "pointrange",
               #width = 0.01,
               size = 0.05,
               lwd = 0.5,
               position = position_nudge(x= c(-0.2,0.2))) +
  theme_classic() +
  theme(axis.title.x = element_blank()) + 
  labs(x = "Assay", y= "Peak frequency (au)") + 
  theme(legend.position="none")
p
ggsave(filename = paste(plotdir,"summary/","summary_peakfreq_pairedSolo", "_peakdistance", peakdistance, "_sds", sds, ".pdf", sep = ""),
       p, width = 4, height = 4, units = 'cm')
#ggplotly(p)

mod = lmer(peakfreq ~ pairedSolo + (1|mouseID), data = expduration.se)
summary(mod)
anova(mod)

######################################################
# PERI BEHAVIOR TIME HISTROGRAMS AND SNAPSHOTS
# UPDATE: see solution from sean harrington "time_series_center_SH_ACN_2022_11_23.R"
# Goals: 
# 1) avoid collisions by creating an independent result for each run of e.g. actiHud_start, and then bind all the results together
# 2) Put plotGroup inside the result as well. 
######################################################
# remove columns 
#remove unwanted columns 
combo2peaks = combo2peaks %>%
  select(-any_of(c("lmDifference","deltaFslide",
                   "File.Name",
                   "Timestamp",
                   "Grooming.Social.states.","Time",
                   "Strain","CageID",
                   "dateTime", "date", "time", "TempRound", "hoboCRound","ActivityRound")))

# indicate whether any given second of data contains a peak (to be used for peak modeling and logistic regression below)
peakValue <- 1
combo2peaks = combo2peaks %>% 
  group_by(Time.stamp.round) %>% 
  mutate(PeakType = case_when(
    any(peakBinary %in% peakValue) ~ 1,
    T ~ 0
  )) %>%
  ungroup()

rm(every_data)
rm(extracted)

#################
# run the extracted and create_data_list functions to get a list of bouts 
#################
#define after & before. 
# before <- fps*15 and after <- fps*180
# before <- fps*15 and after <- fps*60
before <- fps*15
after <- fps*60

#select behavior state start/stop 
# "sleeHud_start","sleeHud_stop" 
# "sleeSol_start", "sleeSol_stop"
# "actiHud_start" "actiHud_stop"
# "Nest_start"     "Nest_stop"
unique(combo2peaks$huddleState)
unique(combo2peaks$bStartStop)
colname = "bStartStop"    #bStartStop will be used to get the indexes of the beginning and end
rows <- "actiHud_start"

#for (k in levels(combo2peaks[["pairedSolo"]])) {
#c2p = combo2peaks[combo2peaks$pairedSolo == k , ]
#exracted function
extracted = extract.with.context(x=combo2peaks, colname=colname, rows=rows, after = after, before = before)
extracted = extracted %>% dplyr::mutate(huddleState = replace_na(as.character(huddleState), "null"))

# plotGroup counter 
start_idx <- which(extracted$bStartStop == rows)
#dataList = list()
span <- seq(from = -before, to = after) # this is the same through the loop, can be outside it

#create dataList and every_data using functions defined above 
dataList <- create_data_list(extracted=extracted, rows=rows, span=span, start_idx)
every_data = do.call(rbind, dataList)
str(every_data)

#Plot all the data columns at once
colNames = c("lmQ.bc.Z")
#colNames = names(df)[1:3]
for(i in colNames){
  plt <- ggplot(data=every_data, aes_string(x="span_2", y = i)) + 
    geom_line(data=every_data, aes_string(x="span_2", y = i, group = "plotGroup"), color = "gray50") +
    geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=1) + 
    geom_smooth(data=every_data, aes_string(x="span_2", y = i), span = 0.1) +
    ggtitle(rows) +
    theme_classic()
  print(plt)
  ggsave(filename = paste(plotdir,"logistic/","periTime_",rows,"_", "_value_",i, "_before",before,"_after",after,"_", ".pdf", sep=""), plt, width = 7, height = 7, units = 'cm')
  Sys.sleep(2)
}

##############################################
# LOGISTIC REGRESSION
# https://bookdown.org/egarpor/PM-UC3M/glm-prediction.html 
# http://www.sthda.com/english/articles/36-classification-methods-essentials/151-logistic-regression-essentials-in-r/ 
##############################################

#logistic regression: https://stats.oarc.ucla.edu/r/dae/logit-regression/
mylogit <- glm(PeakType ~ span_2, data = every_data, family = "binomial")

# Summarize the model
summary(mylogit)

newdata3 <- cbind(every_data, predict(mylogit, newdata = every_data, type = "link",se = TRUE))
newdata3 <- within(newdata3, {
  PredictedProb <- plogis(fit)
  LL <- plogis(fit - (1.96 * se.fit))
  UL <- plogis(fit + (1.96 * se.fit))
})

#plot
# multiplty the y axis by 10 (10 fps) and then by 100 (percent), for a total of 1000
formatter1000 <- function(){
  function(x) x*1000}
formatter0.1 <- function(){
  function(x) x*0.1}
plr = ggplot(newdata3, aes(x = span_2, y = PredictedProb)) + 
  geom_ribbon(aes(ymin = LL,ymax = UL), alpha = 0.2) + 
  geom_line(linewidth = 0.5) + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=1) +
  scale_y_continuous(labels=formatter1000(),limits = c(0, 0.025),breaks = seq(0, 0.025, by = 0.010)) + 
  scale_x_continuous(labels=formatter0.1()) + 
  labs(x = "Time (s)", y= "Peak prob./sec") +
  theme_classic()
plr
ggsave(filename = paste(plotdir,"logistic/", "logisticRegression_", rows,"_", "_before",before,"_after",after,"_",
                        "_peakdistance", peakdistance, "_sds", sds,".pdf", sep=""), plr,
       width = 5, height = 5, units = 'cm')

#save output for downstream model comparison
#assign(paste0("state_",rows,"_",k, sep=""), newdata3)
assign(paste0("state_",rows, sep=""), newdata3)

#remove every_data for next round
rm(every_data); rm(extracted); rm(newdata3)
#}

#compare models

#condense predicted values to the second-level 
state_sleeHud_stop_smry = state_sleeHud_stop %>% #state_sleeHud_stop
  dplyr::group_by(span_2) %>%
  dplyr::summarize(peakBnry = max(peakBinary)) %>%
  dplyr::mutate(state = "sleeHud_stop") %>%
  ungroup()
sum(state_sleeHud_stop_smry$peakBnry)

state_sleeSol_stop_smry = state_sleeSol_stop %>%
  dplyr::group_by(span_2) %>%
  dplyr::summarize(peakBnry = max(peakBinary)) %>%
  dplyr::mutate(state = "sleeSol_stop") %>%
  ungroup()
sum(state_sleeSol_stop_smry$peakBnry)

join = state_sleeHud_stop_smry %>%
  full_join(state_sleeSol_stop_smry) %>%
  dplyr::mutate_at(vars(state), factor)
str(join)
levels(join$state)

ggplot(data=join, aes(x=span_2, y=peakBnry, group = state, color = state)) + 
  geom_jitter() + 
  geom_smooth()

logit <- glm(peakBnry ~ state*span_2, data = join, family = "binomial")
summary(logit)

#for active huddle and nesting #state_actiHud_start and state_Nest_start
state_actiHud_start_smry = state_actiHud_start %>%
  dplyr::group_by(span_2,pairedSolo) %>%
  dplyr::summarize(peakBnry = sum(PeakType)) %>%
  dplyr::mutate(state = "actiHud_start") %>%
  ungroup()
sum(state_actiHud_start_smry$peakBnry)

state_Nest_start_smry = state_Nest_start %>%
  dplyr::group_by(span_2, pairedSolo) %>%
  dplyr::summarize(peakBnry = sum(PeakType)) %>%
  dplyr::mutate(state = "Nest_start") %>%
  ungroup()
sum(state_Nest_start_smry$peakBnry)

join = state_actiHud_start_smry %>%
  full_join(state_Nest_start_smry) %>%
  dplyr::mutate_at(vars(state), factor)
str(join)
levels(join$state)

joinPaired = join[join$pairedSolo == "paired",]
ggplot(data=joinPaired,
       aes(x=span_2, y=peakBnry, group = state, color = state)) + 
  geom_jitter() + 
  geom_smooth()

logit <- glm(peakBnry ~ state*span_2, family = poisson(), data = joinPaired)
summary(logit)

# reduce data to the level of one second
# take the mean per second
# every_data2 = every_data %>%
#   group_by(Time.stamp.round) %>%
#   arrange(Timestamp) %>%
#   #filter(row_number()==1)
#   filter(row_number()==1 | row_number()==n())

# every_data2$rankP <- predict(mylogit, newdata = every_data2, type = "response")
# ggplot(every_data2, aes(x = span_2, y = rankP)) +
#   geom_smooth()
# ggplot(every_data2, aes(span_2, PeakType)) +
#   geom_point(alpha = 0.2) +
#   geom_smooth(method = "glm", method.args = list(family = "binomial"))

##############################################
#
# ACTIVITY BEFORE AND AFTER CALCIUM PEAKS
#
##############################################
rm(every_data); rm(extracted)

#remove outlier activity points
sd(combo2peaks$Activity, na.rm = TRUE)
mean(combo2peaks$Activity, na.rm = TRUE)
combo2peaks <-  combo2peaks %>% mutate(Activity = replace(Activity, Activity > 13, NA))

#zscore activity 
combo2peaks = combo2peaks %>%
  dplyr::group_by(expID) %>% 
  dplyr::mutate(ActivityZ = c(scale(Activity))) %>%
  ungroup()
hist(combo2peaks$ActivityZ)

#remove outlier activity points
#combo2peaks <-  combo2peaks %>% mutate(ActivityZ = replace(ActivityZ, ActivityZ > 13, NA))

#check mean activity in each experiment 
ap = ggplot(combo2peaks, aes(expID, ActivityZ)) + 
  stat_summary(fun.y = mean, geom = "bar") + 
  stat_summary(fun.data = mean_se, geom = "errorbar") + 
  scale_x_discrete(guide = guide_axis(angle = 45)) 
ap
ggsave(filename = paste(plotdir,"ActivityZ_eachExpMean",".pdf", sep=""),width = 15, height = 15, units = 'cm')

#smooth activity 
#make some rolling averages. 20fps * 0.2 = 4
testRoll = head %>%
  mutate(roll20 = rollapply(Activity, width = 20, FUN = mean,fill = "extend"), 
         roll10 = rollapply(Activity, width = 10, FUN = mean,fill = "extend"), 
         roll9 = rollapply(Activity, width = 9, FUN = mean,fill = "extend"), 
         roll5 = rollapply(Activity, width = 5, FUN = mean,fill = "extend"),
         roll4 = rollapply(Activity, width = 4, FUN = mean,fill = "extend"), 
         roll3 = rollapply(Activity, width = 3, FUN = mean,fill = "extend"))
# and plot one or some of them to pick the best "width" 
testRoll %>% 
  pivot_longer(cols=c(Activity, roll9),
               names_to="dataType",
               values_to = "value") %>% 
  ggplot(aes(x=timeS, y=value, color=dataType)) + 
  geom_line()

#smooth the data 
smooth = 8
combo2peaks = combo2peaks %>%
  dplyr::group_by(expID) %>% 
  dplyr::mutate(ActivitySmooth = rollapply(ActivityZ, width = smooth, FUN = mean, fill = "extend")) %>%
  dplyr::ungroup()

#set the span
before <- fps*400
after <- fps*400
#grab peaks
colname = "peakBinary"    #bStartStop will be used to get the indexes of the beginning and end
rows <- 1

#use extracted function to get peaks and spans
extracted = extract.with.context(x=combo2peaks, colname=colname, rows=rows, after = after, before = before)

#put a counter on each span
start_idx <- which(extracted$peakBinary == rows)
#dataList = list()
span <- seq(from = -before, to = after) # this is the same through the loop, can be outside it

#create dataList and every_data using functions defined above 
dataList <- create_data_list(extracted=extracted, rows=rows, span=span, start_idx)
every_data = do.call(rbind, dataList)
str(every_data)

#make a bee swarm with "before" and "after," i.e., negative and positive 
every_data = every_data %>%
  mutate(NegPos = case_when(span_2 <= 0 ~ "before",
                            span_2 > 0 ~ "after"))

#order factor levels 
every_data$NegPos <- factor(every_data$NegPos, levels = c("before","after"))

hist(every_data$ActivitySmooth)
hist(log10(every_data$ActivitySmooth))
hist(log10(every_data$Activity))

#swarm plot all data 
# for (j in levels(every_data[["pairedSolo"]])) {
#   newdf = every_data[every_data$pairedSolo == j, ]
#   pp = ggplot(data = newdf, aes(x=NegPos, y=log10(ActivitySmooth))) +
#     geom_quasirandom(method = "smiley", alpha = 0.05, size=.1) + 
#     theme_classic()
#   ggsave(filename = paste(plotdir,"ActivitySmooth_beeswarm_",sds,"sds_",j, "_",rows,"_value_", "_before",before,"_after",after,"_",".pdf", sep=""), 
#          pp, width = 9, height = 7, units = 'cm')
# }

#normalize the activity value by activity prior to the calcium event 
every_data = every_data %>%
  dplyr::group_by(expID, plotGroup) %>% 
  dplyr::mutate(ActivitySmPreNorm = ActivitySmooth/mean(ActivitySmooth[NegPos == 'before'])) %>%
  ungroup()

#plot activity according to calcium events (span_2)
# for (j in levels(every_data[["pairedSolo"]])) {
#     df = every_data[every_data$pairedSolo == j, ]
#     plt <- ggplot(data=df, aes(x=span_2, y = ActivitySmooth)) + 
#     geom_line(data=df, aes(x=span_2, y = ActivitySmooth, group = "plotGroup"), 
#               color = "gray70", size = 0.001, alpha = 0.5) +
#       scale_size_manual(0.0001) +
#     geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=1) + 
#     geom_smooth(data=df, aes(x=span_2, y = ActivitySmooth), span = 0.01, linewidth =0.5) +
#     theme_classic()
#   print(plt)
#   ggsave(filename = paste(plotdir,"activity/", "Activity_periTime_",sds,"sds_",j, "_","ActivitySmooth", "_before",before,"_after",after,
#                           "_",".pdf",sep=""), 
#          plt, width = 9, height = 7, units = 'cm')
#   Sys.sleep(2)
# }

# get activity before and after peak, per mouseID
#every_data$ActivitySmoothLog <- log10(every_data$ActivitySmooth)
for (i in levels(every_data[["pairedSolo"]])) {
  hsemAct <- summarySE(every_data[every_data$pairedSolo == i, ],
                       measurevar = "ActivitySmooth", groupvars = c("mouseID","NegPos"), na.rm = TRUE)
  px = ggplot(data=hsemAct, aes(x=NegPos, y=ActivitySmooth)) + 
    geom_point(data=hsemAct, aes(x=NegPos, y=ActivitySmooth, group = mouseID, color = mouseID)) + 
    geom_line(data=hsemAct, aes(x=NegPos, y=ActivitySmooth, group = mouseID, color = mouseID)) + 
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 #width = 0.01,
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(-0.2,0.2))
    ) + 
    theme_classic() +
    labs(y="Activity z score") + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none")
  print(px)
  ggsave(filename = paste(plotdir,"activity/","Summary_Activity_periTime_",i,"_",sds,"sds_","_","ActivitySmooth",
                          "_before",before,"_after",after,"_",".pdf", sep=""), px, width = 4, height = 4, units = 'cm')
  
  #ADD STATS SAVE HERE 
}

############
#activity: get per seconds means 
############
# PER SECOND
every_data = every_data %>%
  dplyr::mutate(spanSec = round(span_2/10, digits = 0))

for (i in levels(every_data[["pairedSolo"]])) {
  hsemAct <- summarySE(every_data[every_data$pairedSolo == i, ],
                       measurevar = "ActivitySmooth", groupvars = c("mouseID","NegPos","spanSec"), na.rm = TRUE)
  hsemAct2 <- summarySE(every_data[every_data$pairedSolo == i, ],
                       measurevar = "ActivitySmooth", groupvars = c("NegPos","spanSec"), na.rm = TRUE)
  py = ggplot() +
    #geom_line(data=hsemAct,aes(x=spanSec, y=ActivitySmooth, group = mouseID), size = 0.1, color = "grey90") + 
    geom_vline(data=hsemAct,aes(x=spanSec, y=ActivitySmooth), xintercept = 0, linetype="solid", color = "gray30", linewidth=0.4) + 
    geom_ribbon(data = hsemAct2, aes(x = spanSec, ymin = ActivitySmooth - ci, ymax = ActivitySmooth + ci), fill = "grey60") +
    geom_line(data=hsemAct2,aes(x=spanSec, y=ActivitySmooth), size = 0.09, color = "black") + 
    #geom_smooth(data=hsemAct2,aes(x=spanSec, y=ActivitySmooth), span = 0.4, linewidth = 0.2) +
    theme_classic() +
    labs(y="Activity z score") + 
    theme(axis.title.x = element_blank())
  print(py)
  ggsave(filename = paste(plotdir,"activity/", "Activity_periTime_perSecAvg_",i,"_",sds,"sds_","_","ActivitySmooth",
                          "_before",before,"_after",after,"_",".pdf", sep=""),py,
         width = 5, height = 4, units = 'cm')
}

#collapse the per-second means to pairedSolo means
for (i in levels(every_data[["pairedSolo"]])) {
  hsemAct <- summarySE(every_data[every_data$pairedSolo == i, ],
                       measurevar = "ActivitySmooth", groupvars = c("mouseID","NegPos","spanSec"), na.rm = TRUE)
  hsemAct2 <- summarySE(hsemAct, measurevar = "ActivitySmooth", groupvars = c("mouseID","NegPos"), na.rm = TRUE)
  pz = ggplot(data=hsemAct2, aes(x=NegPos, y=ActivitySmooth)) + 
    geom_point(data=hsemAct2, aes(x=NegPos, y=ActivitySmooth, group = mouseID, color = mouseID)) + 
    geom_line(data=hsemAct2, aes(x=NegPos, y=ActivitySmooth, group = mouseID, color = mouseID)) + 
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 width = 0.01,
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(-0.2,0.2))
    ) + 
    theme_classic() +
    scale_y_continuous(limits = c(-0.45, 0.2),breaks = seq(-0.4, 0.2, by = 0.2)) +
    labs(y="Activity z score") + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position="none")
  print(pz)
  ggsave(filename = paste(plotdir,"activity/","Summary_Activity_periTime_perSecond",i,"_",sds,"sds_","_","ActivitySmooth",
                          "_before",before,"_after",after,"_",".pdf", sep=""), pz,
         width = 4, height = 4, units = 'cm')
  
  mod = lmer(ActivitySmooth ~ NegPos + (1|mouseID), data = hsemAct2)
  summary(mod)
  modPs = summary(mod)$coefficients
  write.csv(modPs,  file = paste(plotdir,"activity/","Summary_Activity_periTime_perSecond",i,"_",sds,"sds_","_","ActivitySmooth",
                                 "_before",before,"_after",after,"_",
                                 ".csv",sep = ""))
}

mod = lmer(ActivitySmooth ~ NegPos + (1|mouseID), data = hsemAct)
summary(mod)
#remove every_data to clear next run 
#rm(every_data)
#rm(extracted)

######################
######################
# Tb before and after calcium peaks--continues from previous section
######################
######################
str(every_data)

#need cols: two Tb datasets; span_2 and spanSec; calcium; mouseID; expID; NegPos; pairedSolo
#https://www.datanovia.com/en/blog/dplyr-how-to-compute-summary-statistics-across-multiple-columns/

every_dataTb = every_data %>%
  #tidyr::unite("expID_mouseID_Ambient_pairedSolo", expID,mouseID,Ambient,pairedSolo) %>%
  dplyr::group_by(spanSec,pairedSolo,expID,mouseID) %>%
  summarise(across(
    .cols = c("Temp", "cagemate"), 
    .fns = list(Mean = mean, SD = sd, SEM = std.error), na.rm = TRUE, 
    .names = "{col}_{fn}"
  )) %>%
  ungroup()
every_dataTb = every_dataTb %>% drop_na(Temp_Mean)

#hsemAct <- summarySE(every_dataTb, measurevar = "Temp", groupvars = c("mouseID","NegPos","pairedSolo", "spanSec"), na.rm = TRUE)
#plot Tb according to calcium events (span_2)
for (j in levels(every_dataTb[["pairedSolo"]])) {
  df = every_dataTb[every_dataTb$pairedSolo == j, ]
  dfSUMM <- summarySE(df,measurevar = "Temp_Mean", groupvars = c("spanSec"), na.rm = TRUE)
  
  plt <- ggplot(data=df, aes(x=spanSec, y = Temp_Mean)) + 
    #geom_line(data=df, aes(x=spanSec, y = Temp_Mean, group = "plotGroup"), color = "gray70", size = 0.1, alpha = 0.5) +
    scale_size_manual(0.01) +
    geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=1) + 
    geom_smooth(data=df, aes(x=spanSec, y = Temp_Mean), span = 0.8, linewidth =0.5, color="magenta") +
    #geom_smooth(data=df, aes(x=spanSec, y = cagemate2), span = 0.01, linewidth =0.5, color="gray60") +
    theme_classic()
  print(plt)
  ggsave(filename = paste(plotdir,"bodyTemp/", "BodyTemp_periTime_",sds,"sds_",j, "_","Tb", "_before",before,"_after",after,"_",".pdf", sep=""), 
         plt, width = 5, height = 5, units = 'cm')
  Sys.sleep(2)
  
  plt2 <- ggplot() + 
    geom_ribbon(data = dfSUMM, aes(x=spanSec, ymin = Temp_Mean - se, ymax = Temp_Mean + se), color = "gray90", alpha = 0.5) + 
    geom_line(data = dfSUMM, aes(x=spanSec, y=Temp_Mean), color = "black") + 
    theme_classic()
  print(plt2)
}

#extract per-individual means
for (j in levels(every_dataTb[["pairedSolo"]])) {
  df = every_dataTb[every_dataTb$pairedSolo == j, ]
  dfSUMM <- summarySE(df,measurevar = "Temp_Mean", groupvars = c("spanSec", "mouseID"), na.rm = TRUE)
  
  numSeconds = 50
  numExps = length(unique(dfSUMM$mouseID))
  dfSUMM$block = rep(1:ceiling(nrow(dfSUMM)/(numSeconds*numExps)), each = numSeconds*numExps)[1:nrow(dfSUMM)]
  
  dfSUMM2 <- summarySE(dfSUMM,measurevar = "Temp_Mean", groupvars = c("block", "mouseID"), na.rm = TRUE)
  
  plt <- ggplot() + 
    geom_line(data=dfSUMM2, aes(x=block, y = Temp_Mean, group = mouseID, color = mouseID), size = 0.5, alpha = 0.5) +
    stat_summary(data=dfSUMM2, aes(x=block, y = Temp_Mean),
                 fun.data = "mean_se", #mean_cl_normal
                 geom = "pointrange",
                 size = 0.1,
                 lwd = 0.6, 
                 color = "black", 
                 position = position_dodge(width = 1)) + 
    theme_classic()
  plt
  
  dfSUMM2_subset = dfSUMM2 %>%
    dplyr::mutate(test = case_when(block == 7 ~ "before",
                                   block == 8 ~ "before", 
                                   block == 16 ~ "post",
                                   block == 17 ~ "post")) %>%
    na.omit()
  
  dfSUMM2_subset_mean = summarySE(dfSUMM2_subset,measurevar = "Temp_Mean", groupvars = c("test", "mouseID"), na.rm = TRUE)
  
  b = ggplot(data = dfSUMM2_subset_mean, aes(x=test, y=Temp_Mean)) + 
    geom_point(data = dfSUMM2_subset_mean, aes(x=test, y=Temp_Mean, group=mouseID)) +
    geom_line(data = dfSUMM2_subset_mean, aes(x=test, y=Temp_Mean, group=mouseID)) + 
    stat_summary(fun.data = "mean_se", 
                 geom = "pointrange",
                 size = 0.05,
                 lwd = 0.5,
                 position = position_nudge(x= c(-0.2,0.2))
    ) 
  b 
}

#plot Tb according to calcium events (span_2) CAGEMATE
for (j in levels(every_dataTb[["pairedSolo"]])) {
  df = every_dataTb[every_dataTb$pairedSolo == j, ]
  plt <- ggplot(data=df, aes(x=spanSec, y = cagemate_Mean)) + 
    #geom_line(data=df, aes(x=spanSec, y = Temp2, group = "plotGroup"), color = "gray70", size = 0.001, alpha = 0.5) +
    #scale_size_manual(0.0001) +
    geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=1) + 
    geom_smooth(data=df, aes(x=spanSec, y = cagemate_Mean), span = 0.8, linewidth =0.5, color="magenta") +
    #geom_smooth(data=df, aes(x=spanSec, y = cagemate2), span = 0.01, linewidth =0.5, color="gray60") +
    theme_classic()
  print(plt)
  ggsave(filename = paste(plotdir,"bodyTemp/","BodyTemp_periTime_cagemate",sds,"sds_",j, "_","Tb", "_before",before,"_after",after,"_",".pdf", sep=""), 
         plt, width = 9, height = 7, units = 'cm')
  Sys.sleep(2)
}

#plot Tb of fp and cagemage: all data 
str(combo2peaks)
for (j in levels(combo2peaks[["pairedSolo"]])) {
  df2 = combo2peaks[combo2peaks$pairedSolo == j, ]
  #df2 = df2 %>%dplyr::filter(lmQ.bc.Z > -3 )
  df2 = df2 %>% dplyr::filter(PeakType == 1)
  df.hsem = summarySE(df2, measurevar="lmQ.bc.Z", groupvars=c("Temp"), na.rm=TRUE) 
  plt2 = ggplot(data= df2, aes(x= Temp, y = lmQ.bc.Z)) +
    #geom_point(data= df2, aes(x= cagemate, y = lmQ.bc.Z), color = "gray60", alpha = 0.4, size = 0.5) + 
    geom_point(data= df2, aes(x= Temp, y = lmQ.bc.Z), color="magenta", alpha = 0.4, size = 0.5) + #color="magenta" # color=Ambient
    #geom_ribbon(data = df.hsem, aes(ymin = lmQ.bc.Z - ci, ymax = lmQ.bc.Z + ci), fill = "grey70") + 
    #geom_line(data = df.hsem, aes(x = Temp, y = lmQ.bc.Z), linewidth = 0.1) + 
    geom_smooth(span=0.05, linewidth = 0.5) + 
    theme_classic() + 
    scale_x_continuous(breaks = seq(35.0, 38.5, by = 0.5), limits = c(35.0, 38.5)) + 
    scale_y_continuous(breaks = seq(0, 40, by = 10), limits = c(-3, 40)) + 
    labs(x = "Peripeak Temp", y = "dF/F")
  plt2
  ggsave(filename = paste(plotdir,"bodyTemp/","BodyTemp_lmQ.bc.Z_fp_",sds,"sds_",j, "_","Tb",
                          "_before",before,"_after",after,"_",".pdf", sep=""), 
         plt2, width = 9, height = 7, units = 'cm')
  
}

df3 = combo2peaks %>% dplyr::filter(PeakType == 1)
max(combo2peaks$Temp, na.rm = T)
min(combo2peaks$Temp, na.rm = T)
mu <- ddply(df3, "pairedSolo", summarise, grp.mean=mean(Temp, na.rm = TRUE))
head(mu)
for (i in unique(df3$pairedSolo)) {
  df.full = combo2peaks[combo2peaks$pairedSolo == i, ]
  df.peaks = df.full[df.full$PeakType == 1, ]
  mu <- ddply(df.peaks,"pairedSolo", summarise, grp.mean=mean(Temp, na.rm = TRUE))
  mu.full <- ddply(df.full, "pairedSolo", summarise, grp.mean=mean(Temp, na.rm = TRUE))
  p3.2 = ggplot() + #color=pairedSolo, #fill=pairedSolo
    geom_histogram(data = df.full, aes(x=Temp,after_stat(ndensity)), 
                   position = "identity", alpha=0.4, bins = 20, color = "white", fill="gray40", linewidth =0.1) +
    geom_histogram(data = df.peaks, aes(x=Temp,after_stat(ndensity)), 
                   position="identity", alpha=0.4, bins = 20, color = "white", fill="red", linewidth =0.1) + 
    #scale_x_continuous(breaks = seq(35.0, 38.0, by = 1), limits = c(35.0, 38.5)) + 
    geom_vline(data=mu, aes(xintercept=grp.mean),linetype="dashed", color = "black") + 
    geom_vline(data=mu.full, aes(xintercept=grp.mean),linetype="dashed", color = "black") + 
    theme_classic() + 
    labs(x = "Body temperature", y = "Density (scaled)") + 
    theme(legend.position="none")
  print(p3.2)
  ggsave(filename = paste(plotdir,"bodyTemp/","BodyTemp___periPeak_Tb_",sds,"sds_",i, ".pdf", sep=""),p3.2,
         width = 5, height = 5, units = 'cm')
  
  mod = lmer(Temp ~ PeakType + (1|mouseID), data =  df.full)
  summary(mod)
  modPs = summary(mod)$coefficients
  write.csv(modPs,  file = paste(plotdir,"bodyTemp/","BodyTemp___periPeak_Tb_",sds,"sds_",i, ".csv", sep=""))
}

#plot Tb of fp and cagemage: means at times of peaks 
#take Temp and cagemate from wide to long 
pairedPeaks = combo2peaks2[combo2peaks2$pairedSolo == "paired" , ]
pairedPeaks = pivot_longer(pairedPeaks, c(Temp, cagemate), values_to = "Temp", names_to = "logger")
hsem = summarySE(pairedPeaks, measurevar="Temp", groupvars=c("mouseID", "logger"), na.rm=TRUE) 
ggplot(data = pairedPeaks, aes(x = logger, y = Temp)) + 
  geom_point(color = "gray40", size =0.5, alpha = 0.5) + 
  stat_summary(fun.data = "mean_se", 
               geom = "pointrange",
               width = 0.01,
               size = 0.05,
               lwd = 0.5,
               position = position_nudge(x= c(0.2,0.2))
  ) + 
  theme_classic()
ggsave(filename = paste(plotdir,"bodyTemp/","BodyTemp_lmQ.bc.Z_fp_cm_peaksOnly_paired",sds,"sds_", "_","Tb",
                        "_before",before,"_after",after,"_",".pdf", sep=""),
       width = 3.5, height = 5, units = 'cm')


#plot relationship between activity and Tb 
# smooth out activity in combo3 (or combo4?) 
# by taking mean activity per minute (could be based on timeMround)
# to eliminate outliers 
# then regress this against Tb, with a line/grouping per mouseID

ggplot(data = combo2peaks, aes(x=ActivitySmooth, y = Temp)) + 
  geom_point()
hist(combo2peaks$ActivitySmooth)
#remove every_data to clear next run 
rm(every_data)
rm(extracted)
rm(dataList) 

##############################################
############################################## 
#
# RUN LENGTH ENCODING OF BEHAVIOR BOUTS AND PEAKS
#
##############################################
##############################################
#optional: eliminate very short bouts. (A value of 1 is a placeholder and should not eliminate any bouts--none are that short)
trimLength <- 1
combo2peaks$huddleState <- inverse.rle(within.list(rle(as.character(combo2peaks$huddleState)), values[lengths<trimLength] <- NA))

#make a column that gives the start frame for each event 
combo3peaks = combo2peaks %>% 
  dplyr::arrange(Time.stamp) %>%
  #group_by(huddleState) %>% 
  dplyr::group_by(grp = rleid(huddleState)) %>%
  dplyr::mutate(startStop = case_when(
    row_number() == n() ~ 'stop',
    row_number() == 1 ~ 'start')) %>%
  as_tibble() %>%
  tidyr::unite("bStartStop", c(huddleState,startStop), remove=FALSE) %>%
  ungroup()

#put a counter (in minutes) going (1) forward from the bout-start, and (2) going backwards from bout-end
combo3peaks = combo3peaks %>% 
  dplyr::group_by(grp) %>%
  #600 is 10fps * 60s
  dplyr::mutate(cum_Rows = row_number() / 600) %>% 
  dplyr::mutate(reverseCount = (row_number() - n()) / 600) %>% 
  #round to the nearest minute
  dplyr::mutate(cum_RowsRnd = round(cum_Rows, digits = 0)) %>%
  dplyr::mutate(reverseCountRnd = round(reverseCount, digits = 0)) %>%
  ungroup()

#####################
#plot bout length for each huddleState2 behavior (i.e. pre-post quiesence) using cum_Rows
#####################
boutLengths_prePost = combo3peaks %>%
  dplyr::group_by(mouseID, Ambient, pairedSolo, huddleState2, grp) %>%
  dplyr::summarise(
    maxLength = max(cum_Rows),
    peakCnt = sum(peakBinary)
  )

#get rid of  NA, and remove spaces 
boutLengths_prePost = boutLengths_prePost[!grepl("NA", boutLengths_prePost$huddleState2), ]
boutLengths_prePost = boutLengths_prePost %>% drop_na(huddleState2)
boutLengths_prePost$huddleState2 = str_replace_all(boutLengths_prePost$huddleState2, pattern=" ", repl="")

#split up huddleState2
boutLengths_prePost = boutLengths_prePost %>% separate(huddleState2, into = c("periQui", "behavior"), sep="_(?=[^_]+$)")

#reorder factor levels 
unique(boutLengths_prePost$periQui)
boutLengths_prePost$periQui = factor(boutLengths_prePost$periQui, levels = c("pre", "post", "both", "neither"))

nestingSOLO = boutLengths_prePost[boutLengths_prePost$pairedSolo == "solo" & boutLengths_prePost$behavior == "Nest", ]
for (i in unique(boutLengths_prePost$behavior)) { 
  newdf = boutLengths_prePost[boutLengths_prePost$behavior == i, ]
  pp = ggplot(data=newdf, aes(x=periQui, y=peakCnt)) + #group = pairedSolo, color=pairedSolo
    #geom_violin(color = "gray60",fill = "gray50", trim = TRUE,alpha = 0.5) +
    #geom_jitter(position=position_jitter(width=0.1, height = 0.01), size = 0.5,fill = "gray70", color = "gray40",alpha=0.4, ) + 
    stat_summary(fun.data = "mean_se", #mean_cl_normal
                 geom = "pointrange",
                 size = 0.01,
                 #alpha=0.5,
                 lwd = 0.6, 
                 color = "black", 
                 position = position_dodge(width = 0.8)
    ) + 
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    theme(legend.position="none", 
          axis.title.x=element_blank())
  print(pp)
  ggsave(filename = paste(plotdir,"bouts/","boutLength_peakCount_periQuiesence",i,"_", sds,"sds_", ".pdf", sep=""),pp,
         width = 3, height = 4, units = 'cm')
}

ggplot(data=boutLengths_prePost,aes(x=peakCnt, fill=periQui)) + 
  geom_histogram() + facet_grid(cols = vars(behavior))

#model
AH_paired = boutLengths_prePost[boutLengths_prePost$behavior == "actiHud" ,]
mod = lmer(peakCnt ~ 0 + periQui + (1|mouseID), data =  AH_paired)
summary(mod)
modPs = summary(mod)$coefficients
write.csv(modPs,  file = paste(plotdir,"bouts/","boutLength_peakCount_periQuiesence_AH_paired",sds,"sds_", 
                               ".csv", sep = ""))

nesting = boutLengths_prePost[boutLengths_prePost$behavior == "Nest" ,]
mod = lmer(peakCnt ~ 0 + periQui + (1|mouseID), data =  nesting)
summary(mod)
modPs = summary(mod)$coefficients
write.csv(modPs,  file = paste(plotdir,"bouts/","boutLength_peakCount_periQuiesence_Nesting",sds,"sds_", 
                               ".csv", sep = ""))


#####################
#determine bout durations of huddleState2
#####################
boutLengths = combo3peaks %>%
  dplyr::group_by(mouseID, Ambient, pairedSolo, huddleState2, grp) %>%
  dplyr::summarise(
    maxLength = max(cum_Rows),
    peakCnt = sum(peakBinary)
  )

#get rid of 'instantaneous' columns and NA 
boutLengths = boutLengths[!grepl("NA", boutLengths$huddleState2), ]
boutLengths = boutLengths[!grepl("cont", boutLengths$huddleState2), ]
boutLengths = boutLengths[!grepl("Sta", boutLengths$huddleState2), ]
boutLengths = boutLengths %>% dplyr::filter(grepl('actiHud|Nest', huddleState2))

#clean up, same as above. 
boutLengths$huddleState2 = str_replace_all(boutLengths$huddleState2, pattern=" ", repl="")
boutLengths = boutLengths %>% separate(huddleState2, into = c("periQui", "behavior"), sep="_(?=[^_]+$)", remove = FALSE)
boutLengths$periQui = factor(boutLengths$periQui, levels = c("pre", "post", "both", "neither"))

boutSums = boutLengths %>%
  dplyr::group_by(mouseID, Ambient, pairedSolo, huddleState2, periQui, behavior) %>%
  dplyr::summarise(
    boutSum = sum(maxLength),
    peakSum = sum(peakCnt)
  )

ggplot(data=boutSums, aes(x=periQui, y=boutSum)) + 
  facet_grid(cols = vars(Ambient)) +
  stat_summary(
               fun.data = mean_se,
               geom = "col",
               width = 0.8,
               color="gray70",
               fill="gray90",
               position = position_dodge(1),
               alpha = 0.3) +
  geom_jitter(width=0.05, size = 0.8, color="gray70", fill="white") + 
  #geom_violin()+
  #geom_point(position=position_dodge(width=1),alpha = 0.9, size = 0.7) + 
  stat_summary(fun.data = "mean_se", 
               geom = "pointrange",
               size = 0.05,
               lwd = 0.8, 
               color = "black", 
               position = position_dodge(width = 1)
  ) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  facet_grid(cols = vars(behavior)) + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","boutCumulative_huddlestate2_",sds,"sds_", ".pdf", sep=""),
       width = 7, height = 6, units = 'cm')

mod = lmer(boutSum ~ periQui + (1|mouseID), data = boutSums[boutSums$behavior == "actiHud" , ])
modTK = emmeans(mod, list(pairwise ~ periQui), adjust = "tukey")
write.csv(tidy(modTK$`pairwise differences of periQui`),
          paste(plotdir,"bouts/","boutCumulative_huddlestate2_activeHuddle",sds,"sds_", ".csv", sep=""))

mod = lmer(boutSum ~ periQui + (1|mouseID), data = boutSums[boutSums$behavior == "Nest" , ])
modTK = emmeans(mod, list(pairwise ~ periQui), adjust = "tukey")
write.csv(tidy(modTK$`pairwise differences of periQui`),
          paste(plotdir,"bouts/","boutCumulative_huddlestate2_activeHuddle",sds,"sds_", ".csv", sep=""))

#####################
#plot bout length for each behavior using cum_Rows i.e. huddleState
#####################
boutLengths = combo3peaks %>%
  dplyr::group_by(mouseID, Ambient, pairedSolo, huddleState, grp) %>%
  dplyr::summarise(
    maxLength = max(cum_Rows),
    peakCnt = sum(peakBinary)
  )

#get rid of 'instantaneous' columns and NA 
boutLengths = boutLengths[!grepl("NA", boutLengths$huddleState), ]
boutLengths = boutLengths[!grepl("cont", boutLengths$huddleState), ]
boutLengths = boutLengths[!grepl("Sta", boutLengths$huddleState), ]

ggplot(data=boutLengths, aes(x=huddleState, y=maxLength, group = pairedSolo, color=pairedSolo)) + 
  geom_jitter(position = position_jitterdodge(), alpha = 0.4) + 
  #geom_bar(stat = "identity") + 
  stat_summary(fun.data = "mean_se", 
               geom = "pointrange",
               size = 0.1,
               lwd = 0.5, 
               color = "black", 
               position = position_dodge(width = 0.8)
  ) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(cols = vars(Ambient))
ggsave(filename = paste(plotdir,"bouts/","boutLength",sds,"sds_", ".pdf", sep=""),
       width = 15, height = 9, units = 'cm')

boutSums = boutLengths %>%
  dplyr::group_by(mouseID, Ambient, pairedSolo, huddleState) %>%
  dplyr::summarise(
    boutSum = sum(maxLength),
    peakSum = sum(peakCnt)
  )

ggplot(data=boutSums, aes(x=huddleState, y=boutSum, group = pairedSolo, color=pairedSolo, fill = pairedSolo)) + 
  facet_grid(cols = vars(Ambient)) +
  stat_summary(aes(fill = pairedSolo),
               fun.data = mean_se,
               geom = "col",
               width = 0.8,
               position = position_dodge(1),
               alpha = 0.3) +
  #geom_jitter(position = position_jitterdodge(), alpha = 0.9) + 
  geom_point(position=position_dodge(width=1),alpha = 0.9, size = 0.7) + 
  stat_summary(fun.data = "mean_se", 
               geom = "pointrange",
               size = 0.1,
               lwd = 0.1, 
               color = "black", 
               position = position_dodge(width = 1)
  ) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(filename = paste(plotdir,"bouts/","boutCumulative",sds,"sds_", ".pdf", sep=""),
       width = 15, height = 9, units = 'cm')

unique(boutSums$huddleState)
boutSumsReduced = boutSums %>%
  dplyr::filter(grepl('actiHud|Nest|sleeHud|sleeSol', huddleState))

#plot relationship between boutSum and peakSum.
#note you cannot subset on ambient: there are too few points to do a regression
xMax = max(boutSumsReduced$boutSum*1.2)
yMax = max(boutSumsReduced$peakSum*1.2)
for (i in unique(boutSumsReduced$huddleState)) { 
  dat = boutSumsReduced[boutSumsReduced$huddleState == i , ]
  #dat = boutSumsReduced[boutSumsReduced$huddleState == i & boutSumsReduced$pairedSolo == "solo", ]
  ptrio = ggplot(data=dat, aes(x=boutSum, y = peakSum)) +
    #facet_grid(cols = vars(huddleState)) +
    geom_point(aes(color = Ambient), size = 2, alpha=0.7) + 
    scale_color_manual(breaks = c("15", "23", "29"),values=c("#1A85FF", "#1AFF1A", "#D41159")) +
    geom_smooth(method="lm", alpha = 0.2, linewidth = 0.4, color = "gray50") +
    #xlim(min = -3, max = xMax) +
    #ylim(min = -3, max = yMax) + 
    theme_classic() +
    theme(legend.position="none")
  print(ptrio)
  ggsave(filename = paste(plotdir,"bouts/","boutSum_peakSum_",i,"_",
                          sds,"sds_", ".pdf", sep=""),ptrio,
         width = 5, height = 5, units = 'cm')
  
  dat$Ambient = as.numeric(levels(dat$Ambient))[dat$Ambient]
  mod = lmer(peakSum ~ boutSum * Ambient + (1|mouseID), data = dat)
  mod = lm(peakSum ~ boutSum * Ambient, data = dat)
  summary(mod)
  modPs = summary(mod)$coefficients
  write.csv(modPs,  file = paste(plotdir,"bouts/","boutSum_peakSum_",i,"_",
                                 sds,"sds_", 
                                 ".csv", sep = ""))
  
  ptrioEff = plot_model(mod, type = "std",
                        title = "",
                        vline.color = "black") + 
    theme_sjplot() +
    coord_cartesian()  
    #label_angle(angle.x = 90)
  print(ptrioEff)
  ggsave(filename = paste(plotdir,"bouts/","boutSum_peakSum_",i,"_",
                          sds,"sds_","forest_", ".pdf", sep=""),ptrioEff,
         width = 4, height = 4, units = 'cm')
  }

#####################
#set parameters
#####################
alpha = 0.4
span = 0.7
#####################
#Sleep huddling data. Summarize peak counts
#####################
sleeHud = combo3peaks[combo3peaks$huddleState == "sleeHud", ]

#loop through each ambient
for(i in levels(sleeHud[["Ambient"]])) { 
  df3 = sleeHud[sleeHud$Ambient == i, ]
  
  #reverseCountRnd df3
  sleeHudRev = df3 %>%
    group_by(reverseCountRnd) %>%
    dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
    ungroup()
  p5 = ggplot(data = sleeHudRev, aes(x=reverseCountRnd, y=numPeaks)) + 
    geom_col(colour="gray40", fill = "white") + 
    geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
    geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
    scale_y_continuous(limits = c(-2, 10),breaks = seq(0, 10, by = 2)) + 
    theme_classic() + 
    theme(axis.title.x = element_blank())
  print(p5)
  ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_sleeHud_end_","ambient", i,"_", sds,"sds_",".pdf", sep=""),p5,
         width = 5, height = 5, units = 'cm')
  
  # cum_RowsRnd df3
  sleeHudFor = df3 %>%
    group_by(cum_RowsRnd) %>%
    dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
    ungroup()
  p6 = ggplot(data = sleeHudFor, aes(x=cum_RowsRnd, y=numPeaks)) + 
    geom_col(colour="gray40", fill = "white") + 
    geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
    geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
    scale_y_continuous(limits = c(-2, 10),breaks = seq(0, 10, by = 2)) + 
    theme_classic() +
    theme(axis.title.x = element_blank())
  print(p6)
  ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_sleeHud_start_","ambient", i,"_", sds,"sds_",".pdf", sep=""),p6,
         width = 5, height = 5, units = 'cm')
}

#reverseCountRnd sleeHud
sleeHudRev = sleeHud %>%
  group_by(reverseCountRnd) %>%
  dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
  ungroup()
ggplot(data = sleeHudRev, aes(x=reverseCountRnd, y=numPeaks)) + 
  geom_col(colour="gray40", fill = "white") + 
  geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
  scale_y_continuous(limits = c(-2, 18),breaks = seq(0, 15, by = 5)) + 
  theme_classic() + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_sleeHud_end_",sds,"sds_",".pdf", sep=""),
       width = 5, height = 5, units = 'cm')
# cum_RowsRnd sleeHud
sleeHudFor = sleeHud %>%
  group_by(cum_RowsRnd) %>%
  dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
  ungroup()
ggplot(data = sleeHudFor, aes(x=cum_RowsRnd, y=numPeaks)) + 
  geom_col(colour="gray40", fill = "white") + 
  geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
  scale_y_continuous(limits = c(-2, 18),breaks = seq(0, 15, by = 5)) +
  theme_classic() +
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_sleeHud_start_",sds,"sds_",".pdf", sep=""),
       width = 5, height = 5, units = 'cm')

#####################
#sleeSol data. #summarize peak counts
#####################
sleeSol = combo3peaks[combo3peaks$pairedSolo == "solo" & combo3peaks$huddleState == "sleeSol", ]

#loop through each ambient
for(i in levels(sleeSol[["Ambient"]])) { 
  df3 = sleeSol[sleeSol$Ambient == i, ]
  
  #reverseCountRnd df3
  sleeSolRev = df3 %>%
    group_by(reverseCountRnd) %>%
    dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
    ungroup()
  p5 = ggplot(data = sleeSolRev, aes(x=reverseCountRnd, y=numPeaks)) + 
    geom_col(colour="gray40", fill = "white") + 
    geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
    geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
    scale_y_continuous(limits = c(-2, 13),breaks = seq(0, 12, by = 4)) + 
    theme_classic() + 
    theme(axis.title.x = element_blank())
  print(p5)
  ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_sleeSol_end_","ambient", i,"_", sds,"sds_",".pdf", sep=""),p5,
         width = 5, height = 5, units = 'cm')
  
  # cum_RowsRnd df3
  sleeSolFor = df3 %>%
    group_by(cum_RowsRnd) %>%
    dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
    ungroup()
  p6 = ggplot(data = sleeSolFor, aes(x=cum_RowsRnd, y=numPeaks)) + 
    geom_col(colour="gray40", fill = "white") + 
    geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
    geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
    scale_y_continuous(limits = c(-2, 13),breaks = seq(0, 12, by = 4)) + 
    theme_classic() +
    theme(axis.title.x = element_blank())
  print(p6)
  ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_sleeSol_start_","ambient", i,"_", sds,"sds_",".pdf", sep=""),p6,
         width = 5, height = 5, units = 'cm')
}

sleeSolRev = sleeSol %>%
  group_by(reverseCountRnd) %>%
  dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
  ungroup()
ggplot(data = sleeSolRev, aes(x=reverseCountRnd, y=numPeaks)) + 
  geom_col(colour="gray40", fill = "white") + 
  geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
  scale_y_continuous(limits = c(-2, 18),breaks = seq(0, 15, by = 5)) +
  theme_classic() + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_sleeSol_stop_",sds,"sds_",".pdf", sep=""),
       width = 5, height = 5, units = 'cm')

sleeSolFor = sleeSol %>%
  group_by(cum_RowsRnd) %>%
  dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
  ungroup()
ggplot(data = sleeSolFor, aes(x=cum_RowsRnd, y=numPeaks)) + 
  geom_col(colour="gray40", fill = "white") + 
  geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
  scale_y_continuous(limits = c(-2, 18),breaks = seq(0, 15, by = 5)) +
  theme_classic() +
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_sleeSol_start_",sds,"sds_",".pdf", sep=""),
       width = 5, height = 5, units = 'cm')

#####################
#Active huddling data. #summarize peak counts
#####################
actiHud = combo3peaks[combo3peaks$huddleState == "actiHud", ]

actiHudRev = actiHud %>%
  group_by(reverseCountRnd) %>%
  dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
  ungroup()
ggplot(data = actiHudRev, aes(x=reverseCountRnd, y=numPeaks)) + 
  geom_col(colour="gray40", fill = "white") + 
  geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
  scale_y_continuous(limits = c(-4, 20),breaks = seq(0, 15, by = 5)) + 
  theme_classic() + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_actiHud_stop_",sds,"sds_",".pdf", sep=""),
       width = 5, height = 5, units = 'cm')

actiHudFor = actiHud %>%
  group_by(cum_RowsRnd) %>%
  dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
  ungroup()
ggplot(data = actiHudFor, aes(x=cum_RowsRnd, y=numPeaks)) + 
  geom_col(colour="gray40", fill = "white") + 
  geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
  scale_y_continuous(limits = c(-4, 20),breaks = seq(0, 15, by = 5)) + 
  theme_classic() +
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_actiHud_start_",sds,"sds_",".pdf", sep=""),
       width = 5, height = 5, units = 'cm')

#####################
#Nest data. #summarize peak counts
#####################
nestn = combo3peaks[combo3peaks$huddleState == "Nest", ]

nestnRev = nestn %>%
  group_by(reverseCountRnd) %>%
  dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
  ungroup()
ggplot(data = nestnRev, aes(x=reverseCountRnd, y=numPeaks)) + 
  geom_col(colour="gray40", fill = "white") + 
  geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
  scale_y_continuous(limits = c(-10, 45),breaks = seq(0, 40, by = 10)) + 
  theme_classic() + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_nest_stop_",sds,"sds_",".pdf", sep=""),
       width = 5, height = 5, units = 'cm')

nestnFor = nestn %>%
  group_by(cum_RowsRnd) %>%
  dplyr::summarize(numPeaks  = sum(peakBinary == 1)) %>%
  ungroup()
ggplot(data = nestnFor, aes(x=cum_RowsRnd, y=numPeaks)) + 
  geom_col(colour="gray40", fill = "white") + 
  geom_smooth(span = span, alpha = alpha, color="magenta4", fill = "magenta") + 
  geom_vline(xintercept = 0, linetype="dotted", color = "gray10", linewidth=0.5) +
  scale_y_continuous(limits = c(-10, 45),breaks = seq(0, 40, by = 10)) + 
  theme_classic() +
  theme(axis.title.x = element_blank())
ggsave(filename = paste(plotdir,"bouts/","summary_periCounts_nest_start_",sds,"sds_",".pdf", sep=""),
       width = 5, height = 5, units = 'cm')

#########################
# determine cumulative time for each bout. 
# determine relationship between bout length and calcium peaks
#########################
combo3peaks = combo3peaks %>%
  tidyr::unite(behaviorBout, c("huddleState", "grp"))

#step 2 determine elapsed time 
str(combo3peaks)
combo3peaks <- combo3peaks %>% 
  arrange(Time.stamp) %>%
  group_by(Noldus_file,behaviorBout) %>% 
  mutate(diff = Time.stamp - lag(Time.stamp),
         diff_secs = as.numeric(diff, units = 'secs')) %>%
  mutate(elapSec = cumsum(replace_na(diff_secs, 0))) %>%
  mutate(elapSecRound = ceiling(elapSec)) %>%
  ungroup()

#summarize bout length and number of peaks
combo3peaksSmry= combo3peaks %>% 
  group_by(mouseID, pairedSolo, Noldus_file,behaviorBout) %>%
  summarise(peakCount= sum(peakBinary>0),
            boutLength= max(elapSec, na.rm = T)) %>%
  separate(behaviorBout, into = c("behavior", "bout"), sep="_(?=[^_]+$)") %>%
  ungroup()
#loop through pairedSocial and states to plot peakcount and boutlength
for (i in levels(combo3peaksSmry[["pairedSolo"]])) {
  for (j in levels(as.factor(combo3peaksSmry[["behavior"]]))) {
    dat = combo3peaksSmry[combo3peaksSmry$pairedSolo == i & combo3peaksSmry$behavior == j, ]
    
    gg = ggplot(data = dat, aes(x=boutLength, y = peakCount)) + 
      geom_point(colour="gray40",)+
      ggtitle(paste(i,j,sep="_")) + 
      geom_smooth(span = 0.8, alpha = 0.3, color="magenta4", fill = "magenta") +
      theme_classic()
    print(gg) 
    ggsave(filename = paste(plotdir,"/bouts/", "bouts_",i,"_", j, "_", sds,"sds_",".pdf", sep=""),gg, 
           width = 5, height = 5, units = 'cm')
  }
  
}

##############################################
##############################################
# PLOT "STACKED" PEAKS FOR CALCIUM DYNAMICS
# Determine 'full width at half maximum' (FWHM)
# peak shape in solo/paired, and in ambient 
##############################################
##############################################
#set the span
before <- fps*10
after <- fps*10
#grab peaks
colname = "peakBinary"    #bStartStop will be used to get the indexes of the beginning and end
rows <- 1
#use extracted function to get peaks and spans
extracted = extract.with.context(x=combo2peaks, colname=colname, rows=rows, after = after, before = before)

#put a counter on each span
start_idx <- which(extracted$peakBinary == rows)
#dataList = list()
span <- seq(from = -before, to = after) # this is the same through the loop, can be outside it

#create dataList and every_data using function defined above 
dataList <- create_data_list(extracted=extracted, rows=rows, span=span, start_idx)
every_data = do.call(rbind, dataList)
str(every_data)

#plot all events 
ggplot(data = every_data, aes(x=span_2, y=lmQ.bc.Z, 
                              #group=pairedSolo, color=pairedSolo
                              )) + 
  #geom_line(color="gray60", size = 0.01) + 
  #geom_smooth(span=0.00001, method = "gam") + 
  geom_smooth(span=0.3, method = "loess") + 
  scale_x_continuous(labels=formatter0.1()) + 
  labs(x = "Time (s)") + 
  theme_classic()

# group by span_2 to calculate mean, sem, sd
# and fwhm (full width half maximum)
# https://stackoverflow.com/questions/75458148/calculating-fwhm-from-spectra-in-r
fwhm <- function(x, y) {
  peak <- which.max(y)
  pre  <- seq(peak)
  post <- peak:length(x)
  half <- y[peak]/2
  xvals <- suppressWarnings(c(approx(y[pre], x[pre], xout = half)$y, 
                              approx(y[post], x[post], xout = half)$y))
  
  ss <- which(x > min(xvals) & x < max(xvals))
  data.frame(x = c(min(xvals), x[ss], max(xvals)),
             y = c(half, y[ss], half))
}
half_width <- fwhm(hsem$span_2, hsem$lmQ.bc.Z)
fwhmValue = (abs(half_width[1,1]) + tail(half_width[1], n= 1))  /  10

hsem <- summarySE(every_data, measurevar="lmQ.bc.Z", groupvars=c("span_2"), na.rm=TRUE) 
ggplot(data = hsem, aes(x=span_2, y=lmQ.bc.Z)) + 
  geom_ribbon(aes(ymin = lmQ.bc.Z - se, ymax = lmQ.bc.Z + se), fill = "grey70") + 
  geom_line() + 
  geom_area(data = half_width, aes(x, y), fill = "red", alpha = 0.2) +
  scale_x_continuous(labels=formatter0.1()) + 
  labs(x = "Time (s)") +
  theme_classic()
ggsave(filename = paste(plotdir,"_overview_avgTrace_allTraces",sds,"sds_",".pdf", sep=""),
       width = 9, height = 7, units = 'cm')

# pairedSolo 
hsem <- summarySE(every_data, measurevar="lmQ.bc.Z", groupvars=c("span_2", "pairedSolo"), na.rm=TRUE) 
ggplot(data = hsem, aes(x=span_2, y=lmQ.bc.Z, group=pairedSolo, fill = pairedSolo)) + 
  geom_ribbon(aes(ymin = lmQ.bc.Z - se, ymax = lmQ.bc.Z + se), alpha = 0.3) + 
  geom_line(aes(color=pairedSolo)) + 
  theme_classic()
ggsave(filename = paste(plotdir,"_overview_avgTrace_pairedSolo",sds,"sds_",".pdf", sep=""),
       width = 9, height = 7, units = 'cm')

# Ambient 
hsem <- summarySE(every_data, measurevar="lmQ.bc.Z", groupvars=c("span_2", "Ambient"), na.rm=TRUE) 
ggplot(data = hsem, aes(x=span_2, y=lmQ.bc.Z, group=Ambient, fill = Ambient)) + 
  geom_ribbon(aes(ymin = lmQ.bc.Z - se, ymax = lmQ.bc.Z + se), alpha = 0.3) + 
  geom_line(aes(color=Ambient)) + 
  theme_classic()
ggsave(filename = paste(plotdir,"_overview_avgTrace_Ambient",sds,"sds_",".pdf", sep=""),
       width = 9, height = 7, units = 'cm')

  
  
  
  
  
  































#############################
# Subset data before and after peaks, to measure Tb, behavior, etc.
#############################
#Define time before and after peak to be extracted
before = 1200 # 120 sec
after = 1800 # 180 sec

#Remove unused columns 
combo2peaks.1 <- combo2peaks %>%
  select(-any_of(c("LD","Time.round.1Sec","Time.stamp.round","Region0G","File.Path",
                   "Region0G.415","Noldus_file","robustfit"
                   )))


combo2peaks.1$expID <- as.factor(combo2peaks.1$expID) 
res$expID <- as.factor(res$expID)

############################# If getting injection, remove injection times
peaks_no_inject <- combo2peaks.1 %>%
  mutate(lmQ.bc.sub = case_when(
    injectionState == "injection" ~ NA_real_,  # Replace injection with NA
    TRUE ~ lmQ.bc  # Use TRUE as the default case
  )) %>%
  mutate(peakMid = case_when(
    injectionState == "injection" ~ NA_real_,  # Replace injection with NA
    TRUE ~ peakMid  # Use TRUE as the default case
  ))
#################################################

peak_subsets <- data.frame()
for (exp in unique(res$expID)) {
  #sub <- peaks_no_inject %>%
  sub <- combo2peaks.1 %>%
    filter(expID == exp)
  res_sub <- res %>%
    filter(expID == exp)
  #It is grabbing the same peaks for each exp 
  for (peak in unique(res_sub$peak)) {.   #JASON, I  CHANGED THIS FROM res_sub$peakMid
    min = peak - before
    max = peak + after 
    sub1 <- sub %>%
      filter(counter >= min & counter <= max,) %>%
      mutate(peakID = peak)
    sub1$time.norm <- sub1$counter - peak
    peak_subsets <- rbind(peak_subsets, sub1)
}}

count.peaks <- peak_subsets %>%
  group_by(expID, mouseID,  Ambient) %>%
  mutate(count_0 = sum(time.norm == 0, na.rm = TRUE)) %>%
  ungroup()
count.peaks <- summarySE(count.peaks, measurevar="count_0", 
          groupvars=c("Ambient", "expID"), na.rm=TRUE)

peak_subsets$Activity = as.numeric(peak_subsets$Activity)

#############################################
# JASON, NOTE THAT I THINK MIN-MAX NORMALIZATION SHOULD BE PER MOUSE 
# NOW DEFINED AS lmQ.bc.MM
# also note that the stacked histogram should be subsetted on paired or solo
#############################################
peak_subsets1 <- peak_subsets %>%
  #Min max normalize calcium traces
  #filter(Temp != "NA") %>%
  #filter(lmQ.bc.sub != "NA") %>%
  mutate(lmQ.bc.norm = (lmQ.bc - min(lmQ.bc)) / (max(lmQ.bc) - min(lmQ.bc))) %>%
  #mutate(Temp.norm = (Temp - min(Temp)) / (max(Temp) - min(Temp))) %>%
  #filter(Temp.norm != "NA") %>%
  mutate(Activity.norm = (Activity - min(Activity)) / (max(Activity) - min(Activity)))
  

peaks_summ <- summarySE(peak_subsets1, measurevar="lmQ.bc", 
                        groupvars=c("mouseID", "time.norm", "expID", "Ambient"), na.rm=TRUE) 
#Temp_summ <- summarySE(peak_subsets1, measurevar="Temp.norm", groupvars=c("mouseID.y", "expID", "time.norm"), na.rm=TRUE) 
Activity_summ <- summarySE(peak_subsets1, measurevar="Activity.norm", 
                           groupvars=c("mouseID", "expID", "time.norm",  "Ambient"), na.rm=TRUE) 

#Plot summarized data 
#Plot is working properly, and peaks are being subsetted. Prolem is coming from peaks being close together,
# and causing multiple to be plotted on same plot. 
############################# Need bigger threshold? ########################################## 
ggplot(peak_subsets1, aes(x=time.norm, y=lmQ.bc)) + 
  #facet_grid(cols = vars(Ambient, expID), rows = vars(peakID)) + 
  facet_grid(cols = vars(Ambient), rows = vars(mouseID)) +
  geom_line(aes(group = mouseID), color = "black") + 
  #geom_line(data = Activity_summ, aes(x = time.norm, y = Activity.norm / 8, group = expID ), color = "green") +
  #stat_smooth(data = Temp_summ, aes(x=time.norm, y=Temp.norm /2.5), 
  #            method="loess", span = 0.3, se=T, color = "red", fill="red", alpha=0.2) +
  #stat_smooth(data = Activity_summ, aes(x=time.norm, y=Activity.norm * 4), 
  #            method="loess", span = 0.2, se=T, color = "darkgreen", fill="darkgreen", alpha=0.2) +
  scale_y_continuous("Calcium",  
    sec.axis = sec_axis(~ . *3, name = "Activity.norm")) + 
  theme_minimal()
ggsave(filename = paste(plotdir,"Peaks_lmQ.bc_raw_withActivity_SD.",sds,"_peakdist.",peakdistance,".pdf", sep=""),
       width = 15, height = 7, units = 'cm')


#### Plot behavior as a histogram over the same time course
peak_subsets1 <- peak_subsets1 %>%
  filter(huddleState != "NA" & pairedSolo == "paired") 
ggplot(peak_subsets1, aes(x=time.norm, fill = huddleState)) +
  geom_histogram(position = 'fill', bins = 150) + 
  geom_vline(xintercept = 0, color = "gray50") + 
  theme_minimal()
#ggsave(filename = paste(plotdir,"Hist_behavior_by_Peaks",".pdf", sep=""),width = 10, height = 6, units = 'cm')
 
###Control for previous data
combo2_NArm <- combo2 %>%
  filter(huddleState != "NA")
ggplot(combo2_NArm, aes(x=timeS, fill = huddleState)) +
  geom_histogram(position = 'fill', bins = 1) + 
  #geom_vline(xintercept = 0, color = "black") + 
  theme_minimal()
ggsave(filename = paste(plotdir,"Hist_behavior_control",".pdf", sep=""),
       width = 8, height = 7, units = 'cm')

##############################
#WHEN WE HAVE INJECTION DATA
#Subset in the same way as peaks, but for injection time 
##############################
#Define time before and after peak to be extracted
before.min = 10
before = 10 * 60 * before.min  #10fps * 60 sec * x min
after.min = 60
after = 10 * 60 * after.min #10 fps * 60 sec * 30 min

combo2$expID <- as.factor(combo2$expID)
combo2$injectionState <- as.factor(combo2$injectionState)
combo2.1 <- combo2 %>% 
  filter(expID != "NA")

inj_subsets <- data.frame()

# Loop over each group
for (grp in unique(combo2.1$expID)) {
  # Subset data for the current group
  exp_sub <- combo2.1[combo2.1$expID == grp, ]
  
  # Find the indices where the event starts and ends for the current group
  event_start <- which(exp_sub$injectionState == "injection")[1]
  event_end <- tail(which(exp_sub$injectionState == "injection"), 1)
  
  # Define the window around the event for the current group
  window_before <- event_start - before
  window_after <- event_end + after
  
  # Subset the data based on the window for the current group
  group_subset <- exp_sub[window_before:window_after, ]
  
  # Store the subset in the list
  inj_subsets <- rbind(inj_subsets, group_subset)
}

inj_subsets <- inj_subsets %>%
  #Min max normalize calcium traces
  #mutate(lmQ.bc.norm = (lmQ.bc - min(lmQ.bc)) / (max(lmQ.bc) - min(lmQ.bc))) %>%
  filter(lmQ.bc != "NA",
         lmQuotient != "NA") %>%
  #mutate(lmQ.norm = (lmQuotient - min(lmQuotient)) / (max(lmQuotient) - min(lmQuotient))) %>%
  mutate(lmQ.sub = lmQuotient) %>%
  mutate(lmQ.sub = case_when(
    injectionState == "injection" ~ NA_real_,  # Replace injection with NA
    TRUE ~ lmQ.sub  # Use TRUE as the default case
  )) %>%
  #mutate(Temp.norm = (Temp - min(Temp)) / (max(Temp) - min(Temp))) #%>%
  mutate(Activity.norm = (Activity - min(Activity)) / (max(Activity) - min(Activity)))

#Normalize time for each exp 
inj_subsets = inj_subsets %>% 
  group_by(expID) %>% 
  mutate(counter = as.numeric(row_number(expID))) %>%
  ungroup()

#bc.summ <- summarySE(inj_subsets, measurevar="lmQ.bc.norm", 
#                        groupvars=c("expID", "counter", "Injection", "Ambient", "injectionState"), na.rm=TRUE) 
sub.summ <- summarySE(inj_subsets, measurevar="lmQ.sub", 
                       groupvars=c("expID", "counter", "Injection", "Ambient", "injectionState"), na.rm=TRUE)
full.summ <- summarySE(inj_subsets %>% filter(injectionState == "injection"),
                       measurevar="lmQuotient", 
                       groupvars=c("expID", "counter", "Injection", "Ambient", "injectionState"), na.rm=TRUE)
#Temp_summ <- summarySE(inj_subsets, measurevar="Temp.norm", groupvars=c("mouseID.y", "time.norm"), na.rm=TRUE) 
Activity_summ <- summarySE(inj_subsets, measurevar="Activity.norm", 
                           groupvars=c("expID", "counter", "Injection",  "Ambient", "injectionState"), na.rm=TRUE) 

#Plot summarized data 
ggplot(sub.summ, aes(x=counter, y=lmQ.sub)) + 
  facet_grid(rows = vars(Injection), cols = vars(Ambient)) +
  geom_vline(data = subset(full.summ, injectionState %in% c("injection")),
             aes(xintercept = counter), color = "darkgray", size = 0.5, alpha = 0.3) +
  geom_line(aes(group = Injection), color = "black") +
  geom_line(data = full.summ, aes(x=counter, y=lmQuotient), color = "darkred") +
  #geom_line(data = Activity_summ, aes(x = counter, y = Activity.norm *1.5 ), color = "darkgreen") +
  stat_smooth(method="lm",  se=T, alpha=0.4, color = "blue") +
  #stat_smooth(data = Activity_summ, aes(x=counter, y=Activity.norm * 15),  color = "darkgreen", fill = "darkgreen",
  #            method="lm",  se=T, alpha=0.4) +
  #scale_y_continuous("Calcium",  
  #                   sec.axis = sec_axis(~ . /1.5 , name = "Activity.norm")) + 
  theme_minimal()
ggsave(filename = paste(plotdir,"Injection_raw_injRemoved_before.",before.min,"_after.",after.min,".pdf", sep=""),
       width = 15, height = 9, units = 'cm')

#Plot linear models of data 
sub.summ <- sub.summ %>% unite("amb_inj", Ambient:Injection, remove=FALSE)
ggplot(sub.summ, aes(x=counter, y=lmQ.sub, group = amb_inj, color = amb_inj)) + 
  geom_vline(data = subset(full.summ, injectionState %in% c("injection")),
             aes(xintercept = counter), color = "darkgray", size = 0.5, alpha = 0.3) +
  #stat_smooth(method="lm",  se=T, alpha=0.4) +
  stat_smooth(method="loess", span = 0.5,  se=T, alpha=0.4) +
  #stat_smooth(data = Activity_summ, aes(x=counter, y=Activity.norm * 15),  color = "darkgreen", fill = "darkgreen",
  #            method="lm",  se=T, alpha=0.4) +
  #scale_y_continuous("Calcium",  
  #                   sec.axis = sec_axis(~ . /1.5 , name = "Activity.norm")) + 
  theme_minimal()
ggsave(filename = paste(plotdir,"Injection_amb.inj_loess.5_injRemoved_before.",before.min,"_after.",after.min,".pdf", sep=""),
       width = 15, height = 9, units = 'cm')

###############################
#Set end of sleep huddle as 0 x value and plot calcium accordingly
###############################



##############################
# Create regression model between calcium and Tb, activity, behavior, etc.
##############################
# Remove na's from Tb
# and normalize to min-max to get rid of negative values for glmer
mod_data <- combo2peaks %>%
  filter(Temp != "NA") %>%
  mutate(lmQ.bc.norm = (lmQ.bc - min(lmQ.bc)) / (max(lmQ.bc) - min(lmQ.bc))) %>%
  mutate(Activity.log = log(Activity))

y = mod_data$Activity.log  # Activity, Temp, 

ggplot(mod_data, aes(x=lmQ.bc.norm, y=y)) + 
  geom_point() + 
  geom_smooth(method="lm", se=TRUE, color = "red", fill="blue")
ggsave(filename = paste(plotdir,"Scatter_lmQ_by_","Activity.log",".pdf", sep=""),
              width = 15, height = 7, units = 'cm')

m1 <- lmer(lmQ.bc.norm ~ Activity.log * Temp + (1|mouseID.y),
            data = mod_data)
summary(m1)

#Cannot get glmer to run properly, takes > 2 hours
m2 <- glmer(lmQ.bc.norm ~ Activity.log + (1|mouseID.y),
            data = mod_data, family = poisson)

AIC(m1)

##############################################
# Try a PCA to see if behavior clusters the data 
##############################################
#First get rid of NA's if necessary 
colSums(is.na(mod_data))
pca_data0 <- mod_data %>%
  filter(huddleState != "NA")

#Remove non-numeric columns
pca_data <- subset(pca_data0, select = c(Activity.log, Temp, lmQ.bc.norm))

#Correlation matrix
pca_corr <- cor(pca_data)
ggcorrplot(pca_corr)
ggsave(filename = paste("Correlation_matrix_", var, ".pdf"), width = 10, height = 10, units = 'cm')

#Applying the PCa
data.pca <- prcomp(pca_data, scale = T)
summary(data.pca)

#Loading values
data.pca$rotation[, 1:2]

#Scree plot
fviz_eig(data.pca, addlabels = TRUE)
ggsave(filename = paste("Scree_plot_", var, ".pdf"), width = 10, height = 10, units = 'cm')

#Biplot
fviz_pca_var(data.pca, col.var = "black")
ggsave(filename = paste("Biplot_", var, ".pdf"), width = 10, height = 10, units = 'cm')

#Contribution of each variable
fviz_cos2(data.pca, choice = "var", axes = 1:2)
ggsave(filename = paste("Contributions_", var, ".pdf"), width = 10, height = 10, units = 'cm')

#Combined biplot and contribution
fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("red", "orange", "green"),
             repel = TRUE)
ggsave(filename = paste("Combined_biplot_contribution_", var, ".pdf"), width = 10, height = 10, units = 'cm')


#Scatterplot for visualization
pca_result <- prcomp(pca_data, 
                 scale=F)
pca_data <- as.data.frame(pca_result$x)
pca_data$huddleState <- pca_data0$huddleState

ggplot(pca_data, aes(x = PC1, y = PC2, color = huddleState)) +
  geom_point() +
  labs(x = "Principal Component 1", y = "Principal Component 2", color = "huddleState")
ggsave(filename = paste("Scatterplot_coloredBY", "Strain_gs", "_", var, ".pdf"), width = 10, height = 10, units = 'cm')






