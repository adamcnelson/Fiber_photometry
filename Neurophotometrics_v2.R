library(tidyverse)
library(gtools) 
library(Rmisc) 
library(data.table)
library(stats)
library(lubridate) 
library(RColorBrewer)
library(magrittr)
library(lme4)
library(emmeans)
library(broom)
require(MASS)
require(caTools)
library(zoo)

library(plotly)
 
#10042022_nonHomecage and tempsweep (TS)
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10042022_nonHomecage")
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10042022_tempsweep")
#10_12_2022
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10122022_homecageprac")
#10142022_temp_redo/5C
setwd('/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10142022_temp_redo/5C')
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10142022_temp_redo/40C")
#10_28_2022 5, TS
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/10_28_22/10_28_22_TempSweep")
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/10_28_22/10_28_22_5C")
#11_01_2022 5, 40, TS
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/10282022_11222022_JR/11_1_22/11_1_22_TS")
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_1_22/11_1_22_5C")
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_1_22/11_1_22_40C")
#11_03_2022 HC
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_3_22_HC")
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_3_22_HC_CHOPPED")
#11_8_22
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_8_22/11_8_22_5C")
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_8_22/11_8_22_40C")
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_8_22/11_8_22_TS")
#11_10_22
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_10_22_HC")
#11_17_22
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_17_22_HC") 
#11_22_22
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_22_22/TS")
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/10282022_11222022_JR/11_22_22/5C")
#Oxycre_Jan_Feb_JR
 #1_31_23_HP 
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/1_31_23_HP/5C") # -
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/1_31_23_HP/40C") # -
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/1_31_23_HP/TS") # -
 #02_02_23_HC
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_02_23_HC") #  
 #02_07_23_HP
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_07_23_HP/5C") #-
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_07_23_HP/40C") #-
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_07_23_HP/TS") #-
 #02_09_23_HC
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_09_23_HC")
 #02_14_23_HP
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/fp_Data/Oxycre_Jan_Feb_JR/02_14_23_HP/5C")
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_14_23_HP/40C")
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_14_23_HP/TS")
 #02_16_23_HC
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_16_23_HC")
 #02_21_23_HC
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_21_23_HC")
 #02_23_23_Cold_HC
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre_Jan_Feb_JR/02_23_23_Cold_HC")
 #03_02_23_HC_Solo
 setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/03_02_23_HC_Solo")
 #5_17_23_HC_Cold
 setwd("//alcova.arcc.uwyo.edu/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/Oxycre1208_May_June/5_17_23_HC_Cold")
 

#TEST_02_13_23
setwd("/Volumes/project/huddlingVideosMicroscopy/neurophotometrics/Bonsai_RawData/TEST_02_13_23")

#PLOT DIRECTORY 
plotdir = paste(getwd(),"/plots/", sep = "")

##########################################S##############################
########################################################################
# PART 1: DATA MASSAGING
# 
########################################################################
########################################################################
# 1. Read in NP data, timestamps, and keydow 
# list file names and sort them in order, so they're in the right sequence
# NP raw data 
alldata <- list.files(pattern = "alldata")
#Aligned time stamp from system time 
computerClock <- list.files(pattern = "data_computerClock")
# keydown
keydown <- list.files(pattern = "keydown")
# in the event of multiple files, keep everything ordered
alldata <- mixedsort(alldata)
computerClock <- mixedsort(computerClock)

# read NP content
NPraw <-
  alldata %>%
  lapply(read.table,
         sep=',',       #warning, make sure is comma ',' sep and not tab '\t' sep
         #skip = 0,      
         header = TRUE,
         #encoding = "UTF-16", 
         #as.is = FALSE,
         stringsAsFactors = F,
         na.strings = c("-"))
# read timestamp content
clockTime <-
  computerClock %>%
  lapply(read.table,
         sep=',',       #warning, make sure is comma ',' sep and not tab '\t' sep
         #skip = 0,      
         header = TRUE,
         #encoding = "UTF-16", 
         #as.is = FALSE,
         stringsAsFactors = F,
         na.strings = c("-"))
# read keydown content
keyDown <-
  keydown %>%
  lapply(read.table,
         sep=',',       #warning, make sure is comma ',' sep and not tab '\t' sep
         #skip = 0,      
         header = FALSE, #keydown does not have a header! 
         #encoding = "UTF-16", 
         #as.is = FALSE,
         stringsAsFactors = F,
         na.strings = c("-"))

# read alldata file names. Filename will allow us to combine separate experiments later on. 
all_filenames <- alldata %>% basename() %>% as.list()
# combine alldata and file name lists
all_dat <- mapply(c, NPraw, all_filenames, SIMPLIFY = FALSE)
# unlist all lists and change column name
all_dat <- rbindlist(all_dat, fill = T)
clockTime <- rbindlist(clockTime, fill = T)
keyDown <- rbindlist(keyDown, fill = T)

# change column name (of the last column) 
names(all_dat)[length(names(all_dat))]<-"File.Path" 
# convert to tibble 
all_dat <- as_tibble(all_dat)

#drop unused columns 
drop.cols <- c('Stimulation', 'Output0', 'Output1', 'Input0', 'Input1')
all_dat = all_dat %>% dplyr::select(-one_of(drop.cols))

######################
# time management
######################
# clockTime: first column is frame index and equivalent to FrameCounter
colnames(clockTime)[1] = "FrameCounter"
colnames(clockTime)[2] = "idkTime"
colnames(clockTime)[3] = "dateTime"

#join alldata and clockTime, and convert time to POSIXct
all_dat2 = all_dat %>% left_join(clockTime, by = "FrameCounter")

# convert dateTime from chr to POSIXct 
#op <- options(digits.secs=3) 
options(digits.secs=6) 
#options(op)

all_dat3 = all_dat2 %>%
  mutate(dateTime = str_replace(dateTime, "T", " ")) %>%
  mutate(dateTime = strptime(dateTime, "%Y-%m-%d %H:%M:%OS")) %>% #this gets up POSIXlt
  mutate(dateTime = as.POSIXct(dateTime, tz = "MST","%Y-%m-%d %H:%M:%S"))

######################
# NP management
######################
# Save mean F value when LEDs are off. Can be used as an underestimation of baseline for dF/F calcs later on
FP.no_led = all_dat3[1,"Region0G"] 
# (1) get rid of first row in all_dat3 (since it's the LED-off baseline)
# (2) make sure that there are an equal number of frames for each channel!!! Imbalance can cause trouble downstream
# Note! The number of frames for channels 1 and 2 will either be equal, OR, they will differ by 1. 

#first observe the number of rows per channel (channel 7 is the baseline)
rowCounts = all_dat3 %>% dplyr::group_by(LedState) %>% dplyr::summarise(n = n())
rowCounts
#conditionally remove rows from all_dat3 
all_dat3 = 
  if(rowCounts[1,2] == rowCounts[2,2]){
    all_dat3[-c(1), ]
  } else {
    all_dat3[-c(1,2), ]
  }
#double check that the number of rows for channels 1 and 2 are equal, and confirm there is no channel 7 anymore 
rowCounts = all_dat3 %>% dplyr::group_by(LedState) %>% dplyr::summarise(n = n())
rowCounts

# number of channels 
FP.n_channels = length(unique(all_dat3$LedState))

#?? get rid of null frames? (there should be none)
all_dat3 = all_dat3 %>%
  dplyr::filter(LedState != 0)

######################
# Eliminate approximately first 1 to 5000 rows of data for each channel. There is typically only noise and artifacts here.
# Or, since converting to NA causes all kinds of problems with the nls model, replace with a mean value
# update: solution found to retaining NAs in nls model: specify "newdata" in augment. E.g. fitted = augment(fit, newdata = sensor1.2)
######################
trim = 2000
# at 20 Hz, 5000 rows = 250 seconds, or 4.2 minutes 
# at 20 Hz, 4000 rows = 200 seconds, or 3.33 min
all_dat3[all_dat3$LedState == 1, ] <- all_dat3[all_dat3$LedState == 1, ] %>% 
  dplyr::mutate(Region0G = ifelse(row_number() <= trim, NA, Region0G))  #mean(head(all_dat3$Region0G, 2000L)) or #NA
all_dat3[all_dat3$LedState == 2, ] <- all_dat3[all_dat3$LedState == 2, ] %>% 
  dplyr::mutate(Region0G = ifelse(row_number() <= trim, NA , Region0G)) #mean(head(all_dat3$Region0G, 2000L)) of # NA

#all_dat3FOO <- all_dat3 %>% dplyr::mutate(Region0G = ifelse(row_number(c(1:100)), mean(all_dat3$Region0G), Region0G))
#all_dat3 %>% dplyr::mutate(Region0G = na_if(Region0G, row_number() <= 10000))
#all_dat3[all_dat3$Region0G [1:5],]

######################
# Smoothing
######################
# Sabatini paper (Suk Joon Lee et al) does a 200ms window across 1000frames per second (fps), for a factor of 0.2*fps
#take a fraction of the data to aid in visualization 
head = all_dat3[all_dat3$LedState == 1, ] %>% slice(20000:22000) #10000:12000
ggplot(head, aes(x=FrameCounter, y = Region0G)) + 
  geom_line() 

#make some rolling averages. 20fps * 0.2 = 4
testRoll = head %>%
  mutate(roll10 = rollapply(Region0G, width = 10, FUN = mean,fill = "extend"), 
         roll5 = rollapply(Region0G, width = 5, FUN = mean,fill = "extend"),
         roll4 = rollapply(Region0G, width = 4, FUN = mean,fill = "extend"), 
         roll3 = rollapply(Region0G, width = 3, FUN = mean,fill = "extend"))
# and plot one or some of them to pick the best "width" 
testRoll %>% 
  pivot_longer(cols=c(Region0G, roll3),
    names_to="dataType",
    values_to = "value"
  ) %>% 
  ggplot(aes(x=FrameCounter, y=value, color=dataType)) + 
  geom_line()

#smooth the data. 
all_dat3 = all_dat3 %>%
  dplyr::group_by(LedState) %>% 
  dplyr::mutate(Region0G = rollapply(Region0G, width = 3, FUN = mean, fill = "extend")) %>%
  dplyr::ungroup()

#max values 
all_dat3 %>% 
  dplyr::group_by(LedState) %>%
  dplyr::filter(Region0G == max(Region0G, na.rm=TRUE))

######################
# Plot Raw data 
######################
#%plot deinterleaved data to evaluate the raw signal
#%note: tilelayout allows us to link the x axes of the two graphs so you can zoom in on putative motion artifacts

# calcium independent 
ci = ggplot(data = all_dat3[all_dat3$LedState == 1, ], aes(FrameCounter, Region0G)) + 
  geom_line(linewidth = 0.4) +
  labs(x = 'ArbTime', y = '415-signal raw isosbestic', title = '') + 
  theme_classic()
ci

# calcium dependent 
cd = ggplot(data = all_dat3[all_dat3$LedState == 2, ], aes(FrameCounter, Region0G)) + 
  geom_line(linewidth = 0.4) +
  labs(x = 'ArbTime', y = '470-signal raw calcium dependent', title = '') + 
  theme_classic()
cd
ggplotly(cd)

#put these together on the same plot 
# use ggplotGrob, grid.draw, and rbind to equalize the x-axes 
gA <- ggplotGrob(ci)
gB <- ggplotGrob(cd)
png(paste(plotdir,"fp_raw_Smoothing", ".png", sep=""), width = 20, height = 15, units = 'cm', res = 500)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))
dev.off()

######################
# Normalize data 
# Correcting for photobleaching and outliers
######################
# From Sage 
# order of operations:
# 1) deinterleave data by flag (LED) and save into a matrix
# 2) fit isosbestic signal with a biexponential decay -- the shape of this
# decay is a good approximation of the CONCENTRATION of GCaMP molecules
# underneath your fiber. it decreases as the photobleach. the amplitude,
# however, is tiny. to adjust for this, we:
# 3) linearly scale the fitted decay to the 470 data using robust fit
# 4) divide the raw 470 data by this scale fit to get a corrected signal
# 5) note: this isn't dF/F but it is INTERNALLY reliable -- that is, you can
# compare the beginning of the recording to the end of the recording. dF/F
# requires a good approximation of baseline. you can use the FP.no_led as an
# underestimation of this -- or determine it empirically. it often isn't
# critical to a sound analysis.

# %fit raw isosbestic with biexponential
# %FP.fit1 = fit(temp_x,FP.data(:,3),'exp2'); %fit raw isosbestic with biexponential
# first define isosbestic vs gcamp (i.e. deinterleave)
isos = all_dat3[all_dat3$LedState == 1, ]
gcamp = all_dat3[all_dat3$LedState == 2, ]

#get frame rate 
fps = isos %>% dplyr::group_by(date = format(dateTime, "%Y-%m-%d"),
                               hour = format(dateTime, "%H"),
                               minute = format(dateTime, "%M"),
                               second = format(dateTime, "%S"),
                               File.Path) %>%
  dplyr::summarize(count=n())
fps = mean(fps$count)
fps = round(fps,0)

#parameter fitting: https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates 
c.0 <- min(all_dat3$Region0G) * 0.5
model.0 <- lm(log(Region0G - c.0) ~ FrameCounter, 
              data = isos, 
              na.action = na.exclude)
all(is.na(isos$Region0G))
summary(model.0); coef(model.0)

# Fit exponential decay model 
# Approach 1: from this blog post: https://martinlab.chem.umass.edu/r-fitting-data/
# this post also does biexponential decay (aka double exp decay), but it's very hard to fit perameters--model throws errors!
# eDecay <- function(t, ampl, tau) {
#   ampl*exp(-t/tau) }
# mode1 <- nls(Region0G ~ eDecay(FrameCounter,myA,myT),
#              data = all_dat3[all_dat3$LedState == 1, ],
#              start=list(myA=0.016,myT=0.013))

# Approach 2: approach: https://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/
# use nls to fit isosbestic with exponential decay. 
# Here's the model we're fitting: y(t)∼yf+(y0−yf)e^−αt
# The measured value y starts at y0 and decays towards yf at a rate α.
# SSasymp is a shortcut that guesses it's own parameters, and instead of fitting the rate constant α directly (as above), 
# it searches for the logarithm of α: y(t)∼yf+(y0−yf)e^−exp(logα)t
# "The only benefit of na.exclude over na.omit is that the former will retain the original number of rows in the data. "
fit <- nls(Region0G ~ SSasymp(FrameCounter, yf, y0, log_alpha), 
           data = isos, 
           na.action=na.exclude) #if it doesnt run with na.exclude, try na.omit

fitted = augment(fit, newdata=isos)
paraters = tidy(fit)

#testing below testing below testing below testing below 
# this page may be helpful when nls SSasymp isnt working: 
# https://stackoverflow.com/questions/18364402/r-nls-singular-gradient
#option 
# starting values
# f <- function(x,a,b) {a * exp(b * x)}
# fm0 <- nls(log(Region0G) ~ log(f(FrameCounter, a, b)), isos, start = c(a = 0.001, b = 0.01))
# nls(y ~ f(x, a, b), dat2, start = coef(fm0))
# #option
# nls(Region0G ~ exp(b * FrameCounter), isos, start = c(b = 0.009), alg = "plinear")
# #option
# nls(Region0G ~ exp(loga + b * FrameCounter), data=isos, start = list(loga = log(1), b = 0.009))
#testing above testing above testing above testing above

#plot isosbesticExponentialDecayFit
qplot(FrameCounter, Region0G, data = augment(fit)) + 
  geom_line(aes(y = .fitted), color = "red") + 
  theme_classic()
ggsave(filename = paste(plotdir,"isosbesticExponentialDecayFit", ".pdf", sep=""), width = 15, height = 4, units = 'cm')

#%linearly scale fit to 470 data using robustfit
#FP.fit2 = robustfit(FP.fit1(temp_x),FP.data(:,2),'bisquare'); %scale using robust fit -- note you can "tune" this function -- but the default is often groovy enough
#FP.lin_fit = FP.fit1(temp_x)*FP.fit2(2)+FP.fit2(1); %save resulting scaled fit

#robust regression (rr) in R: https://stats.oarc.ucla.edu/r/dae/robust-regression/ 
fitted <- tibble::rowid_to_column(fitted, "index")
gcamp <- tibble::rowid_to_column(gcamp, "index")
nrow(fitted); nrow(gcamp)

#rename Region0G in fitted to join with gcamp, and pair it down
fitted <- fitted %>% dplyr::rename("Region0G.415" = "Region0G") %>% 
  dplyr::select(c(index, .fitted, Region0G.415))

#join gcamp and fitted 
gcamp.fitted = gcamp %>% 
  dplyr::left_join(fitted, by="index")

#robust regression with bisquare weighting function
rr.bisquare <- rlm(Region0G ~ .fitted, data=gcamp.fitted, psi = psi.bisquare, na.action = na.exclude)
summary(rr.bisquare)
rr.bisquare$fitted.values
length(rr.bisquare$fitted.values) #this will be too few values if there are NAs in the data@ 
#use augment + newdata to retain NAs (ie full dataset)
rlm.fitted = augment(rr.bisquare, newdata=gcamp.fitted)
rlm.fitted$.fitted

rr.bisquare$coefficients[[1]] #first coefficient is the y intercept 
rr.bisquare$coefficients[[2]] #second coefficient is the first regression coefficient 
#put robustfit fitted values back into dataframe 
gcamp.fitted2 = cbind(gcamp.fitted, robustfit = rlm.fitted$.fitted)

#FP.lin_fit = scale the isos exponential decay with coefficients from robust fit
gcamp.fitted2$FP.lin_fit = gcamp.fitted2$.fitted * rr.bisquare$coefficients[[2]] + rr.bisquare$coefficients[[1]]

# These coefficients below are from Sage's matlab robust fit
#gcamp.fitted2$FP.lin_fit = gcamp.fitted2$.fitted * 0.578916855375855 + 0.00587958385858820

#normalizedF: %divide 470 data by scaled fit
gcamp.fitted2$normalizedF = gcamp.fitted2$Region0G/gcamp.fitted2$FP.lin_fit

# plot 415 data fit with biexponential
p1 = ggplot() + 
  geom_line(data=gcamp.fitted2, aes(x= index, y = Region0G.415)) + 
  geom_line(data=gcamp.fitted2, aes(x=index, y = .fitted), color = "red") + 
  ggtitle("415 data fit with biexponential")
p1  
ggsave(filename = paste(plotdir,"415 data fit with biexponential", ".pdf", sep=""), width = 15, height = 4, units = 'cm')

# plot linearly scaled biexponential fit over 470 data
p2 = ggplot() + 
  geom_line(data=gcamp.fitted2, aes(x=index, y = Region0G)) + 
  geom_line(data=gcamp.fitted2, aes(x=index, y = FP.lin_fit), color = "red") + 
  ggtitle("linearly scaled biexponential fit over 470 data")
p2
ggsave(filename = paste(plotdir,"linearly scaled biexponential fit over 470 data", ".pdf", sep=""), width = 15, height = 4, units = 'cm')

# plot normalzedF 
p3 = ggplot() + 
  #geom_line(data=gcamp.fitted2, aes(x=index, y = Region0G)) + 
  geom_line(data=gcamp.fitted2, aes(x=index, y = normalizedF), color = "black", size = 0.1) + 
  ggtitle("normalizedF")
p3
ggplotly(p3)
ggsave(filename = paste(plotdir,"normalizedF", ".pdf", sep=""), width = 15, height = 4, units = 'cm')

# %alternative normalization method -- subtracting or dividing scaled 415 from 470
# From sage: 
# FP.sfit = fitlm(FP.data(:,3),FP.data(:,2));
#FP.data(:,5) = FP.data(:,2)./FP.sfit.Fitted;
#FP.data(:,6) = FP.data(:,2)-FP.sfit.Fitted;

# Warning! 
#If the bi-exponential model doesnt work, then you wont be able to perform this operation on gcamp.fitted2
#In that case this alternative method would have to be on all_dat3
#therefore make all_dat3 and gcamp.fitted2 interachangeable 

#If biexponential doesnt work, first take all_dat3 long to wide 
all_dat3W = all_dat3 %>%
  pivot_wider(names_from = LedState,
              values_from = Region0G) %>%
  dplyr::rename(Region0G = `2`) %>%
  dplyr::rename(Region0G.415 = `1`) %>%
  dplyr::mutate(Region0G.415 = lag(Region0G.415, n=1, order_by=FrameCounter)) #this is lead or lag 
str(all_dat3W)

#get rid of rows with NA in the data channels
all_dat3W = all_dat3W[!with(all_dat3W,is.na(Region0G)| is.na(Region0G.415)),]

#give an index line so it's compatible with gcamp.fitted2
all_dat3W = tibble::rowid_to_column(all_dat3W, "index")

#now make all_dat3W equivalent to gcamp.fitted3 so it feeds into the next sequence 
gcamp.fitted3 <- all_dat3W

#or change gcamp.fitted2 to gcamp.fitted3 !!!!
gcamp.fitted3 <- gcamp.fitted2

#linear fit of 415 and 470 
scaled415.470 = lm(Region0G ~ Region0G.415 , data = gcamp.fitted3)
scaled415.470$fitted.values
length(scaled415.470$fitted.values) #this will be too few rows if there are NAs in the data 
#use augment to retain NAs (i.e. full dataset) in fitted values
lm.fitted = augment(scaled415.470, newdata=gcamp.fitted3)
lm.fitted$.fitted

#putscaled415.470 fitted values back into dataframe
gcamp.fitted3 = cbind(gcamp.fitted3, scaled415.470 = lm.fitted$.fitted)
# divide 470 by scaled415.470
gcamp.fitted3 = gcamp.fitted3 %>% 
  dplyr::mutate(lmQuotient = Region0G/scaled415.470) %>%
  dplyr::mutate(lmDifference = Region0G-scaled415.470)

# deltaf/f based on 465/405
# Zhe Zhang. F = 465 nm LED / 405 nm LED. And then calculated as (F-F0)/F0 in which F0 was the mean of all data points during the trial. 
meanRegion0G = mean(gcamp.fitted3$Region0G, na.rm=TRUE)
gcamp.fitted3 = gcamp.fitted3 %>%
  dplyr::mutate(deltaFq = (Region0G/Region0G.415)-meanRegion0G/Region0G)

# deltaFfitted. See https://www.nature.com/articles/s41598-021-03626-9.pdf?origin=ppub 
# (Signal - Fitted Control)/Fitted Control
gcamp.fitted3 = gcamp.fitted3 %>%
  dplyr::mutate(deltaFfit = (Region0G - scaled415.470)/scaled415.470)

# deltaFslide
# Based on TDT: dF/F is a relative change metric, uses sliding average window as the baseline signal Fo.
#'Window Duration' can be 3 to 120 seconds. The dF/F calculation, which is (F - Fo)/ Fo is performed on each demodulated stream before differencing occurs. 
# 10fps * 30s = 300
windowRoll = fps*60
gcamp.fitted3 = gcamp.fitted3 %>% 
  dplyr::mutate(slideMean = caTools::runmean(x = gcamp.fitted3$Region0G, k = windowRoll,endrule="keep")) %>%
  dplyr::mutate(deltaFslide = (Region0G - slideMean)/slideMean)

#reproducing sages final plot 
x0 = ggplot() + 
  geom_line(data=gcamp.fitted3, aes(x=index, y = Region0G))  +
  geom_line(data=gcamp.fitted3, aes(x=index, y = scaled415.470), color = "red") + 
  ggtitle("470Raw and scaled415.470")
x0
ggsave(filename = paste(plotdir,"470Raw and scaled415.470", ".pdf", sep=""), x0, width = 15, height = 4, units = 'cm')

x1 = ggplot() + 
  geom_line(data=gcamp.fitted3, aes(x=index, y = Region0G-scaled415.470))  +
  ggtitle("470Raw-scaled415.470")
x1
ggsave(filename = paste(plotdir,"470Raw-scaled415.470", ".pdf", sep=""), x1, width = 15, height = 4, units = 'cm')

x2 = ggplot() + 
  #geom_line(data=gcamp.fitted3, aes(x=index, y = Region0G)) + 
  geom_line(data=gcamp.fitted3, aes(x=index, y = lmQuotient), color = "red", size = 0.1) + 
  ggtitle("470 lmQuotient")
x2
ggplotly(x2)
ggsave(filename = paste(plotdir,"470 lmQuotient", ".pdf", sep=""), x2, width = 15, height = 4, units = 'cm')

#lmquotient by timestamp
p = ggplot() + 
  geom_line(data=gcamp.fitted3, aes(x=dateTime, y = lmQuotient), color = "red", size = 0.1) + 
  ggtitle("470 lmQuotient")
ggplotly(p)

#remember, can't plot raw 470 against normalized because values are divergent! 
x3 = ggplot() + 
  #geom_line(data=gcamp.fitted3, aes(x=index, y = Region0G))  +
  geom_line(data=gcamp.fitted3, aes(x=index, y = normalizedF), color = "red") + 
  ggtitle("470Raw and normalizedF")
x3
ggsave(filename = paste(plotdir,"normalizedF", ".pdf", sep=""), x3, width = 15, height = 4, units = 'cm')

x4 = ggplot() + 
  #geom_line(data=gcamp.fitted3, aes(x=index, y = Region0G))  +
  geom_line(data=gcamp.fitted3, aes(x=index, y = deltaFq), color = "red") + 
  ggtitle("deltaFq")
x4
ggsave(filename = paste(plotdir,"deltaFq", ".pdf", sep=""), x4, width = 15, height = 4, units = 'cm')

x5 = ggplot() + 
  geom_line(data=gcamp.fitted3, aes(x=index, y = deltaFfit), color = "red") + 
  ggtitle("deltaFfit")
x5
ggsave(filename = paste(plotdir,"deltaFfit", ".pdf", sep=""), x5, width = 15, height = 4, units = 'cm')

x6 = ggplot() + 
  #geom_line(data=head(gcamp.fitted3, n = 5000), aes(x=index, y = deltaFslide), color = "red") + 
  geom_line(data=gcamp.fitted3, aes(x=index, y = deltaFslide), color = "red") + 
  ggtitle("deltaFslide")
x6
ggsave(filename = paste(plotdir,"deltaFslide", ".pdf", sep=""), x6, width = 15, height = 4, units = 'cm')

# save gcamp.fitted2 or gcamp.fitted3 as a csv
write.table(gcamp.fitted3, file=paste( "NP_processed", ".csv", sep=""), sep=",", row.names = FALSE)



