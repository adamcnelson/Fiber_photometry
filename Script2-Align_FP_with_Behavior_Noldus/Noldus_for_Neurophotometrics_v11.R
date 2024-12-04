#libraries 
library(tidyverse)
library(Rmisc)
library(stringi)
library(gtools)
library(data.table)
library(reshape2)
library(zoo)
library(caTools) 
library(lubridate)
library(grid)
library(cowplot)
library(egg)
library(data.table)
library(pracma)
library(plotly)

#DIRECTORIES 
setwd("/Volumes/cluster/alcova/huddlingvidmicro/neurophotometrics/fp_Data/oxycre_3645_july2024/07-26-24_paired_29")

#PLOT DIRECTORY 
plotdir = paste(getwd(),"/behavior_plots/", sep = "")
plotdir 

########################################################################
########################################################################
# PART 1: Read in fiber photometry fp data  
# 1. Make a rounded time stamp
# 2. Get time0--the time stamp of when bonsai/video/fp was started
# 3. 
########################################################################
########################################################################
temp = list.files()
temp = temp[grep("NP_processed", temp)]
fp <- read.csv(file= temp, colClasses=c("dateTime"="POSIXct"))
str(fp)

# round the data to the nearest X 
options(digits.secs = 3)
alignRate = ".1s" # 1 minute, 2 minutes, 30 seconds, 
fp = fp %>%
  #WARNING! round_date doesn't round *exactly* right.
  mutate(Time.round.1Sec = lubridate::round_date(dateTime,alignRate))
  # this is the best alternative I could find 
  #dplyr::mutate(Time.round.1Sec = format(fp$dateTime, "%Y-%m-%d %H:%M:%OS1")) %>% 
  #dplyr::mutate(Time.round.1Sec = strptime(Time.round.1Sec, format = "%Y-%m-%d %H:%M:%OS"))
attr(fp$Time.round.1Sec, "tzone") <- "America/Denver"
print(fp$Time.round.1Sec)
str(fp)

#NEXT: 'COMPRESS' fp data to a value to every 0.1s for all the data columns
fp2 = fp[!duplicated(fp$Time.round.1Sec),]

#get Time0: the first dateTime of the Bonsai run, including fp data and the video. 
time0 = fp$dateTime[1]
options(digits.secs=3) 
#options(op) #this resets options! 
time0 <- strptime(time0, "%Y-%m-%d %H:%M:%OS", tz = "America/Denver")
time0

########################################################################
########################################################################
# PART 2: DATA MASSAGING
# APPEND AUTO/MANUAL DATA TOGETHER, WITH FILENAME INCLUDED AS A COLUMN. 

########################################################################
########################################################################
# list file names and sort them in order, so they're in the right sequence
files.floor <- list.files(pattern = "Arena-Subject") #testNoldus_Excel #Arena-Subject #adjusted

#body temp 
therm <- list.files(pattern = c("DAT")) #thermo loggers, if relevant

# read floor file content
floor_content <-
  files.floor %>%
  lapply(read.table,
         sep=',',      #warning, make sure is comma ',' sep and not tab '\t' sep
         skip = 37,     #have to get rid of the metadata
         header = TRUE,
         fill = TRUE,
         skipNul = TRUE,
         #encoding = "UTF-16", # UTF-16LE #UCS-2LE #UTF-16
         #as.is = FALSE,
         #stringsAsFactors = F,
         na.strings = c("-"))
# read thermologger content
therm_content <-
  therm %>%
  lapply(read.table,
         sep=';',            #warning, thermo data is ; separated
         header = FALSE,     #warning, celcius character is alien
         #colClasses= c("Date-Time" = "POSIXct"),
         stringsAsFactors = F,
         skip = 1,     #skip the headers bc they have weird characters 
         na.strings = c("-"))
head(therm_content)

# get rid of first line in each element of the list (which is just units)
floor_content <- lapply(floor_content, function(x) x[-1,])

# read file names
floor_filenames <- files.floor %>% basename() %>% as.list()
therm_filenames <- therm %>% basename() %>% as.list() 

# combine file content list and file name list
all_floor <- mapply(c, floor_content, floor_filenames, SIMPLIFY = FALSE)
all_therm <- mapply(c, therm_content, therm_filenames, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_floor <- rbindlist(all_floor, fill = T)
all_therm <- rbindlist(all_therm, fill = T)

# update thermo colnames  
names(all_therm) <- c("number","Time.stamp","Temp","mouseID")

# convert to tibble 
all_floor <- as_tibble(all_floor)
all_therm <- as_tibble(all_therm)

# change column name (of the last column)
#names(all_floor)[length(names(all_floor))]<-"File.Path" 
all_floor <- all_floor %>%
  dplyr::rename(File.Name = V1)

#make sure you have numeric values (they can sometimes be read in as character)
cols.num <- c("Trial.time","Recording.time", "Activity")
all_floor[cols.num] <- sapply(all_floor[cols.num],as.numeric)

########################################
########################################
# NOLDUS data
########################################
########################################
# add index line: this will be what you use to sort the flow of time (in arbitrary units)
# add frameTime: a cumulative sequence of time according to the framerate
all_floor <- all_floor %>% dplyr::mutate(ArbTime = row_number())
#frameTime
frameLength = as.numeric(nth(all_floor$Trial.time, 2))
all_floor$frameTime <- seq(from = 0,  by = frameLength, length.out = nrow(all_floor))

#convert multiple columns to numeric 
str(all_floor)
cols.num = c("Trial.time","Recording.time","X.center","Y.center","Area","Areachange","Elongation", "Activity")
#make sure 5C pedestal exps are numeric
all_floor[cols.num] <- sapply(all_floor[cols.num],as.numeric)
str(all_floor)

#activity on floor (dont include zero's here, bc that's basically when no mice are present)
ggplot(data = all_floor, aes(x=Activity)) +
  geom_histogram(binwidth=0.2) + 
  theme_classic() +
  scale_x_log10()
ggsave(filename = paste(plotdir,"activity", ".tiff", sep=""),width = 10, height = 7, units = 'cm')

########################################
########################################
# TIME STAMP MANAGEMENT in Noldus, hobo loggers, and Tb loggers
########################################
########################################
# get time stamp and day/night 
# 1. add time0 to Trial.time to get Time.stamp (note there's a problem for recordings that are broken into "clips". For those cases make new column that additively increments by frame.rate for the full set of clips) )
# 2. add an hour column
# 3. add light dark (LD) 
all_floor = all_floor %>%
  mutate(Time.stamp = all_floor$frameTime + time0) %>% 
  mutate(date = date(Time.stamp),
         hour = hour(Time.stamp),
         LD = case_when(hour >= 6 & hour <18 ~ "light", 
                        TRUE ~ "dark")) 
attr(all_floor$Time.stamp, "tzone") <- "America/Denver"

########################################
########################################
# CAUTION, MAKE SURE thermoRate MATCHES LOGGER READ FREQUENCY
thermoRate = "5 seconds" # 1 minute, 2 minutes, 30 seconds,
# caution! make sure Time.stamp.round in the huddling data is consistent with thermo frequency 
# make a rounded time stamp to align with Thermo temp data using thermoRate; 
all_floor = all_floor %>% mutate(Time.stamp.round = lubridate::round_date(Time.stamp,thermoRate)) #30 seconds 
########################################
########################################

# # TEMPERATURE DATA: Tb body temp loggers. Merge with all data 

# Tb body temp thermo logger 
# convert thermo time to POSIXct
#warning, timestamps exported variously as "%d-%m-%Y %H:%M:%S", "%d.%m.%y %H:%M:%S", "%Y.%m.%d %H:%M:%S", "%Y.%d.%m %H:%M:%S"
all_therm$Time.stamp = as.POSIXct(all_therm$Time.stamp, tz = "America/Denver","%d.%m.%Y %H:%M:%S")
#change comma to period, if applicable (eg. 35,93 should be 35.93)
all_therm = all_therm %>% dplyr::mutate(Temp = as.numeric(gsub(",", ".", Temp)))
# convert thermo data long to wide, with IDs as cols and temp as rows 
all_therm <- all_therm %>%
  dplyr::select(-number) %>%
  tidyr::pivot_wider(names_from = c(mouseID), values_from = c(Temp))
#remove all text after underscore
names(all_therm)[2:ncol(all_therm)] <- sub("\\_.*", "", names(all_therm)[2:ncol(all_therm)])
#rename fpmouse to Temp so it has the same col name in the output file 
all_therm = rename(all_therm, Temp = fpmouse)

###########
# create the all_floor2 df, depending on what elements you already have
###########
#if all_floor2 does not exist yet, but HOBO does, join hobo2 with all_floor
if (!exists('all_floor2') & exists('hobo2') & exists('all_therm')) {
  all_floor2 = all_floor %>%
    left_join(hobo2, by = c("Time.stamp.round" = "hoboTime"))
  all_floor2 = all_floor2 %>%
    left_join(all_therm, by = c("Time.stamp.round" = "Time.stamp"))
} else if (!exists('all_floor2') & exists('all_therm')) {
  all_floor2 = all_floor %>%
    left_join(all_therm, by = c("Time.stamp.round" = "Time.stamp"))
} else if (!exists('all_floor2') & !exists('hobo2') & !exists('all_therm')) {
  all_floor2 = all_floor
} 

########################################
########################################
# ASSIGN BEHAVIOR STATES
# for manually scored, HC, social behavior data 
# this section needs to account for the fact that different states have slightly different names 
# and for the fact that there are presence/absence issues with the different states
########################################
########################################
#make a column for the ethogram 
str(all_floor2)

#first, create huddleState by converting behaviors that have variable names in noldus output files 
str(all_floor2)
activeHuddleCols = c("Active.huddle.Social.states.", "Active.huddle.Social.States.")
sleepHuddleCols = c("Sleep.Huddle.Social.states.", "Sleep.huddle.Social.states.")
groomSocialCols = c("Groom.Social.states.","Grooming.Other.Social.states.")
nestingCols = c("Nesting or Nest Building.Social states.",
                "Nesting.or.Building.Nest.Social.states.",
                "Nesting.or.Nest.Building.Social.states.")

# Filter the dataframe to include only existing columns
activeHuddleCols <- intersect(activeHuddleCols, names(all_floor2))
sleepHuddleCols <- intersect(sleepHuddleCols, names(all_floor2))
groomSocialCols <- intersect(groomSocialCols, names(all_floor2))
nestingCols <- intersect(nestingCols, names(all_floor2))

# Function to filter out non-existent columns
filter_existing_cols <- function(df, cols) {
  intersect(cols, names(df))
}

#
all_floor2 = all_floor2 %>%
  mutate(huddleState = 
           ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, activeHuddleCols)] == 1, na.rm = TRUE) > 0, "actiHud",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, sleepHuddleCols)] == 1, na.rm = TRUE) > 0, "sleeHud",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, groomSocialCols)] == 1, na.rm = TRUE) > 0, "groomSoc",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, nestingCols)] == 1, na.rm = TRUE) > 0, "Nest",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, "Contact.initiated.Social.states.")] == 1, na.rm = TRUE) > 0, "contactI", 
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, "Contact.received.Social.states.")] == 1, na.rm = TRUE) > 0, "contactR", 
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, "Solo.Activity.Social.states.")] == 1, na.rm = TRUE) > 0, "soloAct",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, "Grooming.Self.Social.states.")] == 1, na.rm = TRUE) > 0, "groomSel",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, "Solo.Sleep.Social.states.")] == 1, na.rm = TRUE) > 0, "sleeSol",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, "Locomotion.Social.states.")] == 1, na.rm = TRUE) > 0, "LMA",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, "Eating.or.Drinking.Social.states.")] == 1, na.rm = TRUE) > 0, "EatDri",
                  ifelse(rowSums(all_floor2[, filter_existing_cols(all_floor2, "Stationary.Social.states.")] == 1, na.rm = TRUE) > 0, "Sta",
                  NA_character_)))))))))))))
unique(all_floor2$huddleState)

########################################
########################################
# DETERMINE START AND END INDICES OF CONSECUTIVE RUNS (RUN LENGTH ENCODING)
# Goal will be to get n rows of fp before an after an event  
# for (a) manual score HC exps, or for (b) pedestal exps
########################################
########################################

##############
#optional: trim down short runs
##############
trimLength <- 1
all_floor2$huddleState <- inverse.rle(within.list(rle(as.character(all_floor2$huddleState)), values[lengths<trimLength] <- NA))

#make a column that gives the start frame for each event 
all_floor2 = all_floor2 %>% 
  arrange(Trial.time) %>%
  #group_by(huddleState) %>% 
  dplyr::group_by(grp = rleid(huddleState)) %>%
  #dplyr::mutate(Indicator = ifelse(row_number() == 1, 'start', '')) %>% 
  dplyr::mutate(startStop = case_when(
                row_number() == n() ~ 'stop',
                row_number() == 1 ~ 'start')) %>%
  as_tibble() %>%
  tidyr::unite("bStartStop", c(huddleState,startStop), remove=FALSE)
str(all_floor2)

########################################
########################################
# merge all_floor2 with fp data
#  
########################################
########################################
#merge fp and all_floor2 join by "Time.round.1Sec" with Time.stamp
# round the data to the nearest X 
alignRate
all_floor3 = all_floor2 %>%
  mutate(Time.round.1Sec = lubridate::round_date(Time.stamp,alignRate)) #
attr(fp$Time.round.1Sec, "tzone") <- "America/Denver"

#drop File.Path from all_floor3 if it exists because it will interfere with File.Path from fp2 (which is what we want)
all_floor3 = all_floor3 %>% select(-one_of("File.Path"))

#join all_floor3 and fp2
FP.noldus = all_floor3 %>%
  left_join(fp2, by = c("Time.round.1Sec")) #added File.Path on a different machine 

#get rid of rows with NA in Region0G
FP.noldus2 = FP.noldus[!is.na(FP.noldus$Region0G), ]

#get frames per second fps frame rate
fps = FP.noldus2 %>% dplyr::group_by(date = format(dateTime, "%Y-%m-%d"),
                               hour = format(dateTime, "%H"),
                               minute = format(dateTime, "%M"),
                               second = format(dateTime, "%S"),
                               #File.Path
                               ) %>%
  dplyr::summarize(count=n())
fps = mean(fps$count)
fps = round(fps,0)
fps
#or hard code it: 
#fps = 10

#but first get rid of unwanted columns. This df is too big. 
str(FP.noldus2)
FP.noldus2 = FP.noldus2 %>%
  dplyr::select(-any_of(c("Trial.time", "Recording.time", 
                "X.center", "Y.center", "Area", "Areachange", "Elongation", 
                "Activity.state.Highly.active.","Activity.state.Moderately.active.", "Activity.state.Inactive.",
                "Contact.initiated.Social.states.", "Contact.received.Social.states.",
                "Active.huddle.Social.states.","Sleep.Huddle.Social.states.",
                "Grooming.Self.Social.states.","Solo.Sleep.Social.states.",
                "Injection.Human.contact.","No.contact.Human.contact.",
                "Sleep.huddle.Social.states.",
                "Solo.Activity.Social.states.", 
                "Groom.Social.states.",
                "Locomotion.Social.states.",
                "Result.1",
                "ArbTime","frameTime","date", "hour", "grp","FrameCounter",
                "File.Path",
                "FrameCounter.y", ".resid",
                "idkTime", "robustfit")))

#histograms
hist(FP.noldus2$normalizedF)
hist(FP.noldus2$lmQuotient)
hist(FP.noldus2$Temp)

########################
#order levels from least active to most active
########################
unique(FP.noldus2$huddleState)

# PAIRED EXPERIMENTS 
FP.noldus2 = FP.noldus2 %>%
  dplyr::mutate(huddleState= 
                  factor(huddleState,
                         levels = c("sleeHud", "actiHud", 
                                    "sleeSol", "Nest","groomSel", "groomSoc",
                                    "contactI","contactR", "EatDri", "LMA")))


##########
# Z score lmQ.bc
##########
FP.noldus2 = FP.noldus2 %>%
  dplyr::group_by(File.Name) %>%
  dplyr::mutate(lmQ.bc.Z = (lmQ.bc - mean(lmQ.bc))/sd(lmQ.bc)) %>%
  ungroup()
unique(FP.noldus2$File.Name)
str(FP.noldus2)

########################################
######################################## 
#  temperature plots 
########################################
########################################
#since temp is only recorded every second, need to have one sample per second 

FP.noldus2.1 = FP.noldus2 %>%
  group_by(Time.stamp.round) %>%
  dplyr::filter(row_number()==1) %>%
  ungroup()

p.Tb.DeltaF = ggplot(FP.noldus2.1, aes(x=Temp, y=lmQ.bc )) +
  geom_point(color="grey80") + 
  geom_smooth()+
  theme_classic()
p.Tb.DeltaF
ggsave(filename = paste(plotdir,"bodyTemp by lmQ.bc", ".pdf", sep=""),p.Tb.DeltaF, width = 12, height = 10, units = 'cm')

p.Tb.cm.DeltaF = ggplot(FP.noldus2.1, aes(x=cagemate, y=lmQ.bc )) +
  geom_point(color="grey80") + 
  geom_smooth()+
  theme_classic()
p.Tb.cm.DeltaF
ggsave(filename = paste(plotdir,"bodyTempCageMate by lmQ.bc", ".pdf", sep=""),p.Tb.cm.DeltaF, width = 12, height = 10, units = 'cm')


########################################
########################################
# make plots of behavior and temperature data                             #WARNING, NEED TO GO THROUGH AND CHANGE huddlestate TO behaviorstate
########################################
########################################
#PLOT BEHAVIOR
p.ethog = ggplot(FP.noldus2, aes(x = Time.round.1Sec, y = huddleState, fill = huddleState, color = huddleState), ) + 
  #geom_bar(position="stack", stat="identity")
  #geom_line(size = 8) +
  geom_tile()+
  #labs(x = '', y = '', title = '') +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_x_continuous(expand = c(0, 0)) +
  guides(legend.position ="none") + 
  theme(axis.title.x = element_blank())
p.ethog
ggsave(filename = paste(plotdir,"p.ethog", ".pdf", sep=""),width = 15, height = 4, units = 'cm')

p.ethog.blank = p.ethog + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#plot temp for all mice with loggers
p.Tb.pair = ggplot(FP.noldus2, aes(x=Time.round.1Sec, y=Temp)) +
  geom_line(aes(x=Time.round.1Sec, y=cagemate), color = "black") + 
  geom_line(color="magenta") + 
  labs(x = '', y = 'Body temp', title = '') +
  theme(legend.position = "none") +
  theme(legend.position="none") +  
  theme_classic()
p.Tb.pair
ggsave(filename = paste(plotdir,"bodyTemp_pair", ".pdf", sep=""),p.Tb.pair, width = 12, height = 8, units = 'cm')

#just the fp mouse
p.Tb = ggplot(FP.noldus2, aes(x=Time.round.1Sec, y=Temp)) +
  geom_line() + 
  labs(x = '', y = 'Body temp', title = '') +
  theme(legend.position = "none") +
  theme(legend.position="none") +  
  theme_classic()
p.Tb
ggsave(filename = paste(plotdir,"bodyTemp", ".pdf", sep=""),p.Tb, width = 12, height = 8, units = 'cm')
#cagemate
p.Tb.cm = ggplot(FP.noldus2, aes(x=Time.round.1Sec, y=cagemate)) +
  geom_line() + 
  labs(x = '', y = 'Body temp', title = '') +
  theme(legend.position = "none") +
  theme(legend.position="none") +  
  theme_classic()
ggsave(filename = paste(plotdir,"bodyTempCageMate", ".pdf", sep=""),p.Tb.cm, width = 12, height = 8, units = 'cm')

########################################
########################################
# Plot calcium per behavioral state using PEAKS
# 
########################################
########################################
dataType = "lmQuotient"
dataType = "lmQ.bc"
dataType = "lmQ.bc.Z"

#plot a subset of the data 
#head = FP.noldus2 %>% slice(10000:12000)
head = FP.noldus2[5000:30000, ]
#head[[dataType]]
ggplot(head, aes_string(x="Time.round.1Sec", y = dataType)) +  #Time.round.1Sec #Time.round1sec
  geom_line()  

##########
# define parameters
##########
peakdistance = fps*2
SDs = 6

#Find peaks in one column on a subset of data (aka head)
#TESTING: FOR LMQUOTIENT, DEFINE PEAKS ON THE BASELINECORRECTED DATA, THEN REVERT. THIS GETS AROUND SHIFTING BASELINE
peaks <- findpeaks(head$lmQ.bc.Z,                       
                   minpeakdistance = peakdistance,
                   minpeakheight = sd(head$lmQ.bc.Z)*SDs 
                   ) 
head$new_col <- rep(0, nrow(head))
head$new_col <- replace(head$new_col, peaks[,2], values=1) #the second the position/index is the maximum is reached
str(head)
head[head$new_col ==1, ]
ggplot(head, aes_string(x="Time.stamp", y = dataType)) + 
  geom_line() + 
  geom_point(data= head[head$new_col ==1, ], color = "red")

#full dataset for one column (lmQuotient)
peaks <- findpeaks(FP.noldus2$lmQ.bc.Z, 
                   minpeakdistance = peakdistance,
                   minpeakheight = sd(FP.noldus2$lmQ.bc.Z)*SDs
                   ) 
FP.noldus2$peak <- rep(0, nrow(FP.noldus2))
FP.noldus2$peak <- replace(FP.noldus2$peak, peaks[,2], values=1)

FP.noldus2[FP.noldus2$peak ==1, ]
p.peaks = ggplot(FP.noldus2, aes_string(x="timeS", y = dataType)) + 
  geom_line() + 
  geom_point(data= FP.noldus2[FP.noldus2$peak ==1, ], color = "red",size = 0.6) + 
  theme_classic()
p.peaks
ggsave(filename = paste(
  plotdir,"peaks",peakdistance,"_trimlength",trimLength,"_",SDs,"SDs_", dataType, ".png", sep=""), p.peaks, width = 12, height = 9, units = 'cm')

#summarize and plot peak values
hsem <- summarySE(FP.noldus2[FP.noldus2$peak ==1, ], 
                  measurevar=dataType, groupvars=c("huddleState"), na.rm=TRUE) 

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.1) # move them .05 to the left and right
ggplot(FP.noldus2[FP.noldus2$peak ==1, ], aes_string(x="huddleState", y=dataType)) + 
  geom_jitter(color = "grey90") + 
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  #geom_errorbar(data = hsem, aes_string(ymin=dataType-ci, ymax=dataType+ci), colour="black", width=.5, position=pd) +
  theme_classic()
ggsave(filename = paste(
  plotdir,"peak_values_by_state_",dataType,"_peakDistance",peakdistance,"_trimlength",trimLength,"_",SDs,"SDs_",
  ".pdf", sep=""), width = 9, height = 9, units = 'cm')

#plot Tb x DeltaF for specific states
unique(FP.noldus2$huddleState)
str(FP.noldus2)
subset = FP.noldus2[which (FP.noldus2$huddleState == "sleeHud" & FP.noldus2$peak ==1), ]
subset = FP.noldus2[which (FP.noldus2$huddleState == "sleeHud"), ]
#subset = FP.noldus2[which (FP.noldus2$huddleState == "actiHud" & FP.noldus2$peak ==1), ]

p.Tb.DeltaF.state = ggplot(data = subset, aes(x=Temp, y=lmQ.bc )) +
  geom_point(size = 0.5) + 
  #geom_smooth()+
  theme_classic()
p.Tb.DeltaF.state
ggsave(filename = paste(plotdir,"bodyTemp_lmQ_state","_sleeHud", ".png", sep=""),p.Tb.DeltaF.state, width = 8, height = 8, units = 'cm')

p.Tb.cm.DeltaF.state = ggplot(data = subset, aes(x=cagemate, y=lmQ.bc )) +
  geom_point(size = 0.5) + 
  #geom_smooth()+
  theme_classic()
p.Tb.cm.DeltaF.state
ggsave(filename = paste(plotdir,"bodyTemp_lmQ_state","_sleeHud_cagemate", ".png", sep=""),p.Tb.cm.DeltaF.state, width = 8, height = 8, units = 'cm')

p.Tb.DeltaF = ggplot(FP.noldus2[FP.noldus2$peak ==1, ],
                     aes(x=Temp, y=lmQ.bc )) +
  geom_point(color="grey80") + 
  geom_smooth()+
  theme_classic()
p.Tb.DeltaF
ggsave(filename = paste(plotdir,"bodyTemp by lmQ.bc peaks", ".pdf", sep=""),p.Tb.DeltaF, width = 12, height = 10, units = 'cm')

########################################
########################################
# Extract N rows before and after an event, aka peri-stimulus time histogram 
# for (a) manual score HC exps, or for (b) pedestal exps
########################################
########################################

#define after & before. 
# 10 fps
# For before and after the start/end of a quiescence bout: 30 and 180s or 240s
before <- fps*30
after <- fps*60

# fp$huddleState cannot have NA for the lines that follow 
#FP.noldus2.trim = FP.noldus2.trim %>% dplyr::mutate(huddleState = replace_na(as.character(huddleState), "null"))

#################
# select behaviors from (1) social or (2) hot/cold
#################
# (1) social behavior start
#     social: "actiHud_start","actiHud_stop","LMA_start"","LMA_stop","groom_start","groom_stop",
#     "soloAct_start","soloAct_stop","sleeHud_start","sleeHud_stop" 
unique(FP.noldus2$bStartStop)
unique(FP.noldus2$huddleState)
colname = "bStartStop"    #bStartStop will be used to get the indexes of the beginning and end
rows <- "sleeHud_start"
#select behavior: GROUPED: "sleeHud", "actiHud", "groom", "soloAct", "LMA". SOLO: "sleep"
#behavior <- "sleeSol"

extract.with.context <- function(x, colname, rows, after = 0, before = 0) {
  match.idx  <- which(x[[colname]] %in% rows)
  span       <- seq(from = -before, to = after)
  extend.idx <- c(outer(match.idx, span, `+`))
  extend.idx <- Filter(function(i) i > 0 & i <= nrow(x), extend.idx)
  extend.idx <- sort(unique(extend.idx))
  return(x[extend.idx, , drop = FALSE]) 
}
extracted = extract.with.context(x=FP.noldus2, colname=colname, rows=rows, after = after, before = before)
extracted = extracted %>% dplyr::mutate(huddleState = replace_na(as.character(huddleState), "null"))

#################
# 1) avoid collisions by creating an independent result for each run of e.g. actiHud_start, and then bind all the results together
# 2) Put plotGroup inside the result as well. 
#################
# plotGroup counter 
start_idx <- which(extracted$bStartStop == rows)
#dataList = list()
span <- seq(from = -before, to = after) # this is the same through the loop, can be outside it

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

#create dataList and every_data using functions defined above 
dataList <- create_data_list(extracted=extracted, rows=rows, span=span, start_idx)
every_data = do.call(rbind, dataList)
str(every_data)

#Plot all the data columns at once
colNames = c("lmQ.bc","nF.bc","lmQuotient","normalizedF")
#colNames = names(df)[1:3]
for(i in colNames){
  plt <- ggplot(data=every_data, aes_string(x="span_2", y = i)) + 
    geom_line(data=every_data, aes_string(x="span_2", y = i, group = "plotGroup"), color = "gray90") +
    geom_vline(xintercept = 0, linetype="dotted", color = "gray80", linewidth=1) + 
    geom_smooth(data=every_data, aes_string(x="span_2", y = i), span = 0.1) +
    theme_classic() +
    labs(x = "time since event (deciseconds)", y = "dF/F")
  print(plt)
  ggsave(filename = paste(plotdir,rows,"_value_",i, "_before",before,"_after",after,"_", "trimLength_", trimLength, ".png", sep=""), plt, width = 9, height = 7, units = 'cm')
  Sys.sleep(2)
}

rm(every_data)
rm(extracted)

########################################
########################################
# PLot: behavior state, fp, activity, hobo loggers; Tb loggers;
########################################
########################################
p.hudd = ggplot(FP.noldus2, aes(Time.round.1Sec, huddleState)) + 
  #geom_bar(position="stack", stat="identity")
  #geom_line(size = 8) +
  geom_tile(aes(fill = huddleState),show.legend = FALSE)+
  labs(x = '', y = 'huddleState', title = '') +
  theme(legend.position = "none") +
  theme(legend.position="none") +  
  theme_classic() 
#p.hudd

p.lmQ = ggplot(FP.noldus2, aes(Time.round.1Sec, lmQuotient)) +
  geom_line(size = 0.4) +
  labs(x = 'Time', y = 'lmQuotient', title = '') + 
  #scale_x_continuous(expand = c(0, 0)) +
  theme_classic()
#p.lmQ
ggsave(filename = paste(plotdir,"p.lmQuotient", ".pdf", sep=""),p.lmQ, width = 15, height = 7, units = 'cm')

p.norm = ggplot(FP.noldus2, aes(Time.round.1Sec, normalizedF)) +
  geom_line(size = 0.4) +
  labs(x = 'Time', y = 'lmQuotient', title = '') + 
  theme_classic()
p.norm
ggsave(filename = paste(plotdir,"p.norm", ".pdf", sep=""),p.norm, width = 15, height = 7, units = 'cm')

p.lmQ.bc = ggplot(FP.noldus2, aes(Time.round.1Sec, lmQ.bc)) +
  geom_line(size = 0.4) +
  labs(x = '', y = 'lmQuotient', title = '') + 
  #scale_x_continuous(expand = c(0, 0)) +
  theme_classic()
p.lmQ.bc
#ggplotly(p.lmQ.bc)
ggsave(filename = paste(plotdir,"p.lmQ.bc", ".pdf", sep=""),p.lmQ.bc, width = 15, height = 7, units = 'cm')

p.lmQ.bc.Z = ggplot(FP.noldus2, aes(Time.round.1Sec, lmQ.bc.Z)) +
  geom_line(size = 0.4) +
  labs(x = '', y = 'lmQuotient', title = '') + 
  #scale_x_continuous(expand = c(0, 0)) +
  theme_classic()
p.lmQ.bc.Z
#ggplotly(lmQ.bc.Z)
ggsave(filename = paste(plotdir,"p.lmQ.bc.Z", ".pdf", sep=""),p.lmQ.bc.Z, width = 15, height = 7, units = 'cm')

p.deltaFslide = ggplot(FP.noldus2, aes(Time.round.1Sec, deltaFslide)) +
  geom_line(size = 0.4) +
  labs(x = 'Time', y = 'lmQuotient', title = '') + 
  #scale_x_continuous(expand = c(0, 0)) +
  theme_classic()
#p.deltaFslide
ggsave(filename = paste(plotdir,"p.deltaFslide", ".pdf", sep=""),p.deltaFslide, width = 15, height = 7, units = 'cm')

p.lmQ.bc.Temp = ggplot(FP.noldus2, aes(Temp, lmQ.bc)) +
  geom_jitter(size = 0.4) +
  geom_smooth() +
  theme_classic() 
p.lmQ.bc.Temp
ggsave(filename = paste(plotdir,"p.lmQ.bc.Temp", ".png", sep=""),p.lmQ.bc.Temp, width = 15, height = 7, units = 'cm')

p.lmQ.bc.Temp.cm = ggplot(FP.noldus2, aes(cagemate, lmQ.bc)) +
  geom_jitter(size = 0.4) +
  #labs(x = 'Time', y = '', title = '') + 
  geom_smooth() +
  theme_classic() 
p.lmQ.bc.Temp.cm
ggsave(filename = paste(plotdir,"p.lmQ.bc.Temp.cm", ".pdf", sep=""),p.lmQ.bc.Temp.cm, width = 15, height = 7, units = 'cm')

p.nF.bc = ggplot(FP.noldus2, aes(Time.round.1Sec, nF.bc)) +
  geom_line(size = 0.4) +
  labs(x = '', y = 'normalizedF', title = '') + 
  theme_classic()
#p.nF.bc
ggsave(filename = paste(plotdir,"p.nF.bc", ".pdf", sep=""),p.nF.bc, width = 15, height = 7, units = 'cm')

p.delta = ggplot(FP.noldus2, aes(Time.round.1Sec, deltaFq)) +
  geom_line(size = 0.4) +
  labs(x = 'Time', y = 'deltaF/F', title = '') + 
  geom_smooth() +
  theme_classic()
#p.delta
#ggsave(filename = paste(plotdir,"p.deltaF", ".pdf", sep=""),p.delta, width = 15, height = 7, units = 'cm')


#Activity: get rid of artifacts (i.e. above 7)
#FP.noldus2 = mutate(FP.noldus2, ActivityTrim = case_when(Activity > 9 ~ NA_real_, TRUE ~ Activity))
p.activ = ggplot(FP.noldus2, aes(Time.round.1Sec, Activity)) +
  geom_line(size = 0.4) +
  labs(x = '', y = 'Activity (pixels)', title = '') + 
  theme_classic()
p.activ
ggsave(filename = paste(plotdir,"p.activ", ".pdf", sep=""),p.activ, width = 15, height = 7, units = 'cm')

p.activ.lmQ.bc = ggplot(FP.noldus2, aes((Activity), lmQ.bc)) +
  geom_jitter(size = 0.4) +
  geom_smooth() +
  theme_classic()
p.activ.lmQ.bc
ggsave(filename = paste(plotdir,"p.activ.lmQ.bc", ".pdf", sep=""),p.activ.lmQ.bc, width = 10, height = 10, units = 'cm')

# all_floor2.2
#put these together on the same plot 
# use ggplotGrob, grid.draw, and rbind to equalize the x-axes 
gA <- ggplotGrob(p.lmQ)
gB <- ggplotGrob(p.hudd)
png(paste(plotdir,"p.lmQ.hudd", ".png", sep=""), width = 20, height = 15, units = 'cm', res = 500)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))
dev.off()

########################################
########################################
# PLot: multi-panel plots
########################################
########################################

###################
# down sample df 
# Downsampling the data
###################
df_sampled <- FP.noldus2 %>% sample_frac(0.25)  # Keep X% of the points

p.lmQ.bc.Z_sampled = ggplot(df_sampled, aes(Time.round.1Sec, lmQ.bc.Z)) +
  geom_line(size = 0.4) +
  labs(x = '', y = 'lmQuotient', title = '') + 
  theme_classic()
p.lmQ.bc.Z_sampled
#ggplotly(lmQ.bc.Z)
ggsave(filename = paste(plotdir,"p.lmQ.bc.Z_sampled", ".pdf", sep=""),p.lmQ.bc.Z_sampled, width = 15, height = 7, units = 'cm')

p.activ.sampled = ggplot(df_sampled, aes(Time.round.1Sec, Activity)) +
  geom_line(size = 0.4) +
  labs(x = '', y = 'Activity (pixels)', title = '') + 
  theme_classic()
p.activ.sampled

p.Tb.sampled = ggplot(df_sampled, aes(x=Time.round.1Sec, y=Temp)) +
  geom_line() + 
  labs(x = '', y = 'Body temp', title = '') +
  theme(legend.position = "none") +
  theme(legend.position="none") +  
  theme_classic()
p.Tb.sampled

####################
# DeltaF & behavior 
####################
#DOWNSAMPLED
# p.ethog, p.lmQ.bc
x4 = egg::ggarrange(p.ethog+
                      theme(axis.text.x=element_blank()) + 
                      theme(axis.title.y = element_blank()),
                    ggplot() + theme_void(),
                    p.lmQ.bc.Z_sampled+
                      theme(axis.title.y = element_blank()),
                    heights = c(0.5*(5/8),  #multiply by *(5/8) to make solo compatible with paired
                                -0.09, #-0.09
                                0.5))
x4
ggsave(filename = paste(plotdir, "p.ethog,p.lmQ.bc.Z_sampled", ".pdf", sep=""), useDingbats = TRUE,
       x4, width = 8, height = 6, 
       units = 'cm')
ggsave(filename = paste(plotdir, "p.ethog,p.lmQ.bc.Z_sampled", ".png", sep=""), x4, width = 10, height = 8, units = 'cm')

#p.ethog, p.lmQ.bc, p.Tb, p.activ.sampled
x4 = egg::ggarrange(p.ethog+
                      theme(axis.text.x=element_blank()) + 
                      theme(axis.title.y = element_blank()),
                    ggplot() + theme_void(),
                    p.Tb.sampled + 
                      #theme(axis.title.y = element_blank()) +
                      theme(axis.text.x=element_blank()) +
                    labs(y = "Tb(C)"),
                    ggplot() + theme_void(),
                    p.activ.sampled +
                      #theme(axis.title.y = element_blank()) + 
                      theme(axis.text.x=element_blank()) +
                      labs(y = "Locomotion"),
                    ggplot() + theme_void(),
                    p.lmQ.bc.Z_sampled+
                      #theme(axis.title.y = element_blank()) +
                      labs(y = "dF/F"),
                    heights = c(0.4,  #multiply by *(5/8) to make solo compatible with paired
                                -0.25, #-0.09
                                0.4,
                                -0.25,
                                0.4,
                                -0.25,
                                0.4
                                ))
x4
ggsave(filename = paste(plotdir, "p.ethog,p.lmQ.bc.Z_sampled,p.Tb, p.activ.sampled", ".pdf", sep=""), useDingbats = TRUE,
       x4, width = 9, height = 9, 
       units = 'cm')
ggsave(filename = paste(plotdir, "p.ethog,p.lmQ.bc.Z_sampled,p.Tb, p.activ.sampled", ".png", sep=""), x4, width = 10, height = 8, units = 'cm')

#
x4 = egg::ggarrange(p.ethog.blank, p.lmQ.bc,  heights = c(0.5, 0.5))
x4
ggsave(filename = paste(plotdir, "p.ethog,p.lmQ.bc", ".pdf", sep=""), x4, width = 15, height = 10, units = 'cm')

#alternative approach: p.ethog, p.lmQ, Activity
#cowplot::plot_grid(p.ethog, p.lmQ,  p.activ, align = "v", ncol = 1, rel_heights = c(0.5, 0.2))
x4 = egg::ggarrange(p.ethog, p.lmQ, p.activ,  heights = c(0.5, 0.5,0.5))
ggsave(filename = paste(plotdir, "p.ethog,p.lmQ,p.activ", ".pdf", sep=""), x4, width = 15, height = 10, units = 'cm')

#alternative approach: p.ethog, p.norm, Activity
#cowplot::plot_grid(p.ethog, p.norm,  p.activ, align = "v", ncol = 1, rel_heights = c(0.5, 0.2))
x4 = egg::ggarrange(p.ethog, p.norm, p.activ,  heights = c(0.5, 0.5,0.5))
ggsave(filename = paste(plotdir, "p.ethog,p.norm,p.activ", ".pdf", sep=""), x4, width = 15, height = 10, units = 'cm')

#alternative approach: p.ethog, p.nF.bc
#cowplot::plot_grid(p.ethog, p.nF.bc,  align = "v", ncol = 1, rel_heights = c(0.5, 0.2))
x4 = egg::ggarrange(p.ethog, p.nF.bc,  heights = c(0.5, 0.5))
ggsave(filename = paste(plotdir, "p.ethog,p.nF.bc", ".pdf", sep=""), x4, width = 15, height = 10, units = 'cm')

#alternative approach: p.ethog, p.deltaFslide
#cowplot::plot_grid(p.ethog, p.deltaFslide,  align = "v", ncol = 1, rel_heights = c(0.5, 0.2))
x4 = egg::ggarrange(p.ethog, p.deltaFslide,  heights = c(0.5, 0.5))
ggsave(filename = paste(plotdir, "p.ethog, p.deltaFslide", ".pdf", sep=""), x4, width = 15, height = 10, units = 'cm')

####################
# DeltaF,behavior,therm
####################
#alternative approach: p.hudd, p.Tb, p.lmQ.bc
x4 = egg::ggarrange(p.ethog+ theme(axis.text.x=element_blank()),
                    ggplot() + theme_void(),
                    p.Tb.pair+ theme(axis.text.x=element_blank()),
                    ggplot() + theme_void(),
                    p.lmQ.bc,  
                    heights = c(0.5, -0.1, 0.5,-0.2, 0.5))
x4
ggsave(filename = paste(plotdir, "p.ethog,p.Tb.pair,p.lmQ.bc", ".pdf", sep=""), x4, width = 15, height = 10, units = 'cm')


x4 = egg::ggarrange(p.Tb.pair+ theme(axis.text.x=element_blank()),
                    ggplot() + theme_void(),
                    p.lmQ.bc,  
                    heights = c(0.5, -0.1, 0.5))
ggsave(filename = paste(plotdir, "p.Tb.pair,p.lmQ.bc", ".pdf", sep=""), x4, width = 15, height = 10, units = 'cm')

#solo mice
x4 = egg::ggarrange(p.ethog+ theme(axis.text.x=element_blank()),
                    ggplot() + theme_void(),
                    p.Tb+ theme(axis.text.x=element_blank()),
                    ggplot() + theme_void(),
                    p.lmQ.bc, 
                    heights = c(0.5, -0.1, 0.5,-0.2, 0.5))
ggsave(filename = paste(plotdir, "p.ethog,p.Tb,p.lmQ.bc", ".pdf", sep=""), x4, width = 15, height = 10, units = 'cm')

########################################
# save FP.noldus2 as a csv
########################################
#get rid of extraneous columns
str(FP.noldus2)

#make an export file 
export = FP.noldus

#extract a unique identifier from File.Path [note there may/maynot be a .y added to col ]
identifier = str_extract(export$File.Path[1], "(?<=alldata)(.+)(?=\\.)")

export = FP.noldus2 %>%
  dplyr::select(-any_of(c(contains("states"), 
                contains("Dome"),
                contains("X"),
                "In.zone", 
                "IN.ARENA",
                #- Contact.initiated.Social.states.,
                #- Contact.received.Social.states.,
                "Active.huddle.Social.states.",
                #-Sleep.huddle.Social.states.,
                #-Sleep.Huddle.Social.states.,
                #-Solo.Activity.Social.states., 
                #-Groom.Social.states.,
                "Locomotion.Social.states.",
                #-Result.1,
                contains("File.Path.x"),
                "index",
                "LedState",
                ".fitted",
                "FP.lin_fit",
                "scaled415.470",
                "slideMean",
                "deltaFq",
                "deltaFfit"
  )))

# save gcamp.fitted2 or gcamp.fitted3 as a csv 
write.table(export, file=paste(plotdir, "NP_noldus_processed_",identifier, ".csv", sep=""), sep=",", row.names = FALSE)


