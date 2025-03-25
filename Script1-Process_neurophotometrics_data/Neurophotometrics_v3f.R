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
library(ggprism)
require(baseline)
library(signal) 


##########
# Automatically set the main directory to one level above the script's current location 
maindir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data directory (directly inside Fiber_photometry-main)
datadir <- file.path(maindir, "datafiles")

# Create a directory to store plots (if it doesn't already exist)
plotdir <- file.path(datadir, "plots/")
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}

# Print directories to verify
print(paste("Main directory:", maindir))
print(paste("Data directory:", datadir))
print(paste("Plot directory:", plotdir))

########################################################################
########################################################################
# PART 1: DATA MASSAGING
# 
########################################################################
########################################################################
# 1. Read in NP data, timestamps 
setwd(datadir)
# list file names and sort them in order, so they're in the right sequence
# NP raw data 
alldata <- list.files(pattern = "alldata")
#Aligned time stamp from system time 
computerClock <- list.files(pattern = "data_computerClock")

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

# read alldata file names. Filename will allow us to combine separate experiments later on. 
all_filenames <- alldata %>% basename() %>% as.list()
# combine alldata and file name lists
all_dat <- mapply(c, NPraw, all_filenames, SIMPLIFY = FALSE)
# unlist all lists and change column name
all_dat <- rbindlist(all_dat, fill = T)
clockTime <- rbindlist(clockTime, fill = T)

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
#butterworth filter
######################
# a very simple function that pads the series with a reversed values before/after and then 
# subsets the result to the original range of interest can solve this problem.
# https://stackoverflow.com/questions/53950001/why-is-this-butterworth-filter-presenting-different-results-in-r-and-matlab 
ButterEndEffect <- function(filt,x) {
  signal::filtfilt(filt,c(rev(x),x,rev(x)))[(length(x) + 1):(2 * length(x))]
}

#define filter and parameters
# see https://www.mathworks.com/help/signal/ref/buttord.html#d126e10450
order <- 4  # Filter order
Hz <- 30 
cutoff <- 2/(Hz/2) # Cutoff frequency
b <- butter(order, cutoff, type="low")  # Design the filter coefficients

#examine on subset 
head = all_dat3[all_dat3$LedState == 2, ] %>% slice(30000:70000) #10000:12000
ggplot(head, aes(x=FrameCounter, y = Region0G)) + 
  geom_line() +
  theme_classic()
ggsave(filename = paste(plotdir,"butterworth_before", ".png", sep=""), width = 8, height = 5, units = 'cm')

#spectrum(head$Region0G)
#apply the filter 
head = head %>%
  dplyr::group_by(LedState) %>% 
  dplyr::mutate(Region0G = ButterEndEffect(b,Region0G))  %>%
  dplyr::ungroup()
ggplot(head, aes(x=FrameCounter, y = Region0G)) + 
  geom_line() +
  theme_classic()
ggsave(filename = paste(plotdir,"butterworth_after", ".png", sep=""), width = 8, height = 5, units = 'cm')
#spectrum(head$Region0G)

#incorporate butterworth filter
all_dat3 = all_dat3 %>%
  dplyr::group_by(LedState) %>% 
  dplyr::mutate(Region0G = ButterEndEffect(b,Region0G))  %>%
  dplyr::ungroup() 

######################
# Eliminate approximately first 1 to 5000 rows of data for each channel. There is typically only noise and artifacts here.
# Or, since converting to NA causes all kinds of problems with the nls model, replace with a mean value
# update: solution found to retaining NAs in nls model: specify "newdata" in augment. E.g. fitted = augment(fit, newdata = sensor1.2)
######################
#30 Hz
# at 30 Hz per channel, 3000 rows = 100 seconds (3000/30 = 100)
# at 30 Hz, 2000 rows = 66 seconds 
# at 30 Hz, 20,000 rows = 666 seconds, or 11 minutes

# Define the trimming function
trim_region0G <- function(data, led_state, trim, trim2) {
  data %>%
    dplyr::filter(LedState == led_state) %>%
    mutate(Region0G = ifelse(row_number() <= trim | row_number() >= trim2, NA, Region0G))
}

# Parameters
trim <- 78000
expLength = nrow(all_dat3[all_dat3$LedState == 1, ])
trim2 <- expLength - 5400

# Apply the function for both LedState values and combine the results
all_dat3 <- bind_rows(
  trim_region0G(all_dat3, led_state = 1, trim = trim, trim2 = trim2),
  trim_region0G(all_dat3, led_state = 2, trim = trim, trim2 = trim2),
  dplyr::filter(all_dat3, !LedState %in% c(1, 2)) # Keep rows with other LedState values
)

######################
# Add relative time in seconds 
######################
str(all_dat3)
all_dat3$FrameCounter

all_dat3 = all_dat3 %>%
  dplyr::mutate(diff_row = Timestamp - lag(Timestamp)) %>%
  dplyr::mutate(timeS = cumsum(ifelse(is.na(diff_row), 0, diff_row)) + diff_row*0) %>%
  dplyr::select(-diff_row)

######################
# Smoothing
######################
# Sabatini paper (Suk Joon Lee et al) does a 200ms window across 1000frames per second (fps), for a factor of 0.2*fps
#take a fraction of the data to aid in visualization 
head = all_dat3[all_dat3$LedState == 2, ] %>% slice(90000:110000) #10000:12000
ggplot(head, aes(x=FrameCounter, y = Region0G)) + 
  geom_line() 

######
# spectral density analysis
######
# spectrum(head$Region0G)
# perio <- spec.pgram(fast = TRUE, head$Region0G, spans = c(3, 5, 3), taper = 0.2,
#                     plot = TRUE)
# acf(head$FrameCounter, lag.max = 200)
# mspect <- spectrum(head$Region0G, log="yes", spans=c(2,2), plot=TRUE)
# delta <- 10/60
# specx <- mspect$freq/delta
# specy <- 2*mspect$spec
# plot(specx, specy, xlab="Period", ylab="Spectral Density", type="l")
# max(specy)
# #butterworth filter
# bf <- butter(1, 2.376219e-06, type="high")
# b <- filter(bf, head$Region0G)
# plot(head$FrameCounter, b, col="black", pch=20)

######
# rolling averages
######
#make some rolling averages. 20fps * 0.2 = 4
testRoll = head %>%
  mutate(roll20 = rollapply(Region0G, width = 20, FUN = mean,fill = "extend"), 
         roll10 = rollapply(Region0G, width = 10, FUN = mean,fill = "extend"), 
         roll9 = rollapply(Region0G, width = 9, FUN = mean,fill = "extend"), 
         roll5 = rollapply(Region0G, width = 5, FUN = mean,fill = "extend"),
         roll4 = rollapply(Region0G, width = 4, FUN = mean,fill = "extend"), 
         roll3 = rollapply(Region0G, width = 3, FUN = mean,fill = "extend"))
# and plot one or some of them to pick the best "width" 
testRoll %>% 
  pivot_longer(cols=c(Region0G, roll9),
    names_to="dataType",
    values_to = "value"
  ) %>% 
  ggplot(aes(x=timeS, y=value, color=dataType)) + 
  geom_line()

#smooth the data. 
#Traditionally set to 3 (As of June 2023)
#At 20 fps, 1 frame = 0.05s. A width of 5 gives 0.05*5 = 0.25 seconds, or 250 ms. 
#At 30 fps, 1 frame = 0.033s. A width of 8 gives 0.033*8 = 0.266s, or 266ms
smooth = 9
all_dat3 = all_dat3 %>%
  dplyr::group_by(LedState) %>% 
  dplyr::mutate(Region0G = rollapply(Region0G, width = smooth, FUN = mean, fill = "extend")) %>%
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

# downsample so you can embed pdfs in adobe illustrator 
df_sampled <- all_dat3 %>% sample_frac(0.30)  # Keep X% of the points

# calcium independent 
ci = ggplot(data = df_sampled[df_sampled$LedState == 1, ], aes(timeS, Region0G)) + 
  geom_line(linewidth = 0.4, color = "magenta3") +
  labs(x = 'timeS', y = '415-isosbestic', title = '') + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 25)) +
  guides(x = guide_prism_minor()) +
  theme_classic()
ci

# calcium dependent 
cd = ggplot(data = df_sampled[df_sampled$LedState == 2, ], aes(timeS, Region0G)) + 
  geom_line(linewidth = 0.4, color="green3") +
  labs(x = 'timeS', y = '470-GCaMP', title = '') + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 25)) +
  guides(x = guide_prism_minor()) +
  theme_classic()
cd
#ggplotly(cd)

#put these together on the same plot 
# use ggplotGrob, grid.draw, and rbind to equalize the x-axes 
gA <- ggplotGrob(ci)
gB <- ggplotGrob(cd)
png(paste(plotdir,"fp_raw_Smoothing","_lowPass",cutoff, "_.png", sep=""), width = 12, height = 10, units = 'cm', res = 300)
#pdf(paste(plotdir,"fp_raw_Smoothing","_lowPass",cutoff, "_.pdf", sep=""), width = 5, height = 4)
grid::grid.newpage()
grid::grid.draw(rbind(gA, gB))
dev.off()

######################
# Normalize data 
# Correcting for photobleaching and outliers
######################
# From Sage Aronson
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
fps

# Fit exponential decay model 

# Approach 1: 
# https://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/
# use nls to fit isosbestic with exponential decay. 
# Here's the model we're fitting: y(t)∼yf+(y0−yf)e^−αt
# The measured value y starts at y0 and decays towards yf at a rate α.
# SSasymp is a shortcut that guesses it's own parameters, and instead of fitting the rate constant α directly (as above), 
# it searches for the logarithm of α: y(t)∼yf+(y0−yf)e^−exp(logα)t
# "The only benefit of na.exclude over na.omit is that the former will retain the original number of rows in the data. "
fit <- nls(Region0G ~ SSasymp(timeS, yf, y0, log_alpha), 
           data = isos, 
           na.action=na.exclude) 

fitted = augment(fit, newdata=isos)
parameters = tidy(fit)

#plot isosbesticExponentialDecayFit
qplot(timeS, Region0G, data = augment(fit)) + 
  geom_line(aes(y = .fitted), color = "red") + 
  theme_classic()
ggsave(filename = paste(plotdir,"isosbesticExponentialDecayFit", ".pdf", sep=""), width = 15, height = 4, units = 'cm')

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
#rr.bisquare$fitted.values
length(rr.bisquare$fitted.values) #this will be too few values if there are NAs in the data@ 
#use augment + newdata to retain NAs (ie full dataset)
rlm.fitted = augment(rr.bisquare, newdata=gcamp.fitted)
#rlm.fitted$.fitted

rr.bisquare$coefficients[[1]] #first coefficient is the y intercept 
rr.bisquare$coefficients[[2]] #second coefficient is the first regression coefficient 
#put robustfit fitted values back into dataframe 
gcamp.fitted2 = cbind(gcamp.fitted, robustfit = rlm.fitted$.fitted)

#FP.lin_fit = scale the isos exponential decay with coefficients from robust fit
gcamp.fitted2$FP.lin_fit = gcamp.fitted2$.fitted * rr.bisquare$coefficients[[2]] + rr.bisquare$coefficients[[1]]

#normalizedF: %divide 470 data by scaled fit
gcamp.fitted2$normalizedF = gcamp.fitted2$Region0G/gcamp.fitted2$FP.lin_fit

# plot 415 data fit with biexponential
p1 = ggplot() + 
  geom_line(data=gcamp.fitted2, aes(x= timeS, y = Region0G.415)) + 
  geom_line(data=gcamp.fitted2, aes(x=timeS, y = .fitted), color = "red") + 
  ggtitle("415 data fit with biexponential") + 
  theme_classic()
#p1  
ggsave(filename = paste(plotdir,"415 data fit with biexponential", ".png", sep=""),p1, width = 15, height = 4, units = 'cm')

# plot linearly scaled biexponential fit over 470 data
p2 = ggplot() + 
  geom_line(data=gcamp.fitted2, aes(x=timeS, y = Region0G)) + 
  geom_line(data=gcamp.fitted2, aes(x=timeS, y = FP.lin_fit), color = "red") + 
  ggtitle("linearly scaled biexponential fit over 470 data") + 
  theme_classic()
#p2
ggsave(filename = paste(plotdir,"linearly scaled biexponential fit over 470 data", ".png", sep=""),p2, width = 15, height = 4, units = 'cm')

# plot normalzedF 
p3 = ggplot() + 
  #geom_line(data=gcamp.fitted2, aes(x=timeS, y = Region0G)) + 
  geom_line(data=gcamp.fitted2, aes(x=timeS, y = normalizedF), color = "black", size = 0.1) + 
  ggtitle("normalizedF") + 
  theme_classic()
#p3
ggsave(filename = paste(plotdir,"normalizedF", ".png", sep=""), p3, width = 15, height = 4, units = 'cm')

# Warning! 
#If the bi-exponential model doesnt work, then you wont be able to perform this operation on gcamp.fitted2
#In that case this alternative method would have to be on all_dat3.
#therefore make all_dat3 and gcamp.fitted2 interachangeable.
# %alternative normalization method -- dividing scaled 415 from 470

#If biexponential doesnt work, first take all_dat3 long to wide 
all_dat3W = all_dat3 %>%
  pivot_wider(names_from = LedState,
              values_from = Region0G) %>%
  dplyr::rename(Region0G = `2`) %>%
  dplyr::rename(Region0G.415 = `1`) %>%
  dplyr::mutate(Region0G.415 = lag(Region0G.415, n=1, order_by=timeS)) #this is lead or lag 
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
#lm.fitted$.fitted

#putscaled415.470 fitted values back into dataframe
gcamp.fitted3 = cbind(gcamp.fitted3, scaled415.470 = lm.fitted$.fitted)
# divide 470 by scaled415.470
gcamp.fitted3 = gcamp.fitted3 %>% 
  dplyr::mutate(lmQuotient = Region0G/scaled415.470) %>%
  dplyr::mutate(lmDifference = Region0G-scaled415.470)

#plot scaled415.470 
#y=scaled415.470$coefficients[[1]] #first coefficient is the y intercept 
#fc=scaled415.470$coefficients[[2]] #second coefficient is the first regression coefficient 
ggplot() + 
  #geom_line(data=gcamp.fitted3, aes(x=timeS, y = scaled415.470), color = "red", size = 0.1) + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = lmQuotient), color = "black", size = 0.1) + 
  #ggtitle("scaled415.470") + 
  ggtitle("GCaMP signal divided by the scaled isosbestic") +
  theme_classic()
ggsave(filename = paste(plotdir,"lmQuotient", ".png", sep=""), width = 15, height = 4, units = 'cm')

#OTHER METHODS OF NORMALIZATION
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
windowRoll = fps*20
gcamp.fitted3 = gcamp.fitted3 %>% 
  dplyr::mutate(slideMean = caTools::runmean(x = gcamp.fitted3$Region0G, k = windowRoll,endrule="keep")) %>%
  dplyr::mutate(deltaFslide = (Region0G - slideMean)/slideMean)

######################
# Baseline correction 
#
######################
# following Luping Yin...Dayu Lin Et al 2022
# 'Matlab funciton msbackadj with a moving window of 25% of the total recording duration was applied to the raw signal
# to obtain the instantaneous baseline signal 
# try this: https://cran.r-project.org/web/packages/baseline/baseline.pdf 

# use a pause to ensure the base R plots get saved appropriately
testit <- function(x){
  p1 <- proc.time()
  Sys.sleep(x)
  proc.time() - p1 # The cpu usage should be negligible
}

# Isolate lmQuotient, then do a baseline correction 
metric = gcamp.fitted3$lmQuotient
metric[is.na(metric)]<-mean(metric[1:trim*1.01],na.rm=TRUE)
metric.bc = baseline(matrix(metric, nrow=1), method='irls')
# Plot with base R by opening a jpeg file, then close the file.
plot(metric.bc)
testit(7)
png(paste(plotdir, "baselineCorrected_compare_lmQ.png", sep=""), width = 500, height = 300)
plot(metric.bc)
dev.off()
#get corrected data back into df
corrected = getCorrected(metric.bc)
#make sure the corrected vector length is same length as the df
length(corrected); length(gcamp.fitted3$lmQuotient)
gcamp.fitted3["lmQ.bc"] = corrected[1,]
#plot 
x8 = ggplot() + 
  #geom_line(data=head(gcamp.fitted3, n = 5000), aes(x=timeS, y = deltaFslide), color = "red") + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = lmQ.bc), color = "red") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 200)) +
  #guides(x = guide_prism_minor()) +
  ggtitle("lmQ baselineCorrected")
#x8
ggsave(filename = paste(plotdir,"baselineCorrected.lmQ", ".png", sep=""), x8, width = 15, height = 4, units = 'cm')

#down-sample so you can embed in adobe illustrator 
# downsample so you can embed pdfs in adobe illustrator 
df_sampled <- gcamp.fitted3 %>% sample_frac(0.05)  # Keep X% of the points
#z score 
df_sampled = df_sampled %>%
  dplyr::mutate(lmQ.bc.Z = (lmQ.bc - mean(lmQ.bc))/sd(lmQ.bc))

x8.1 = ggplot() + 
  geom_line(data=df_sampled, aes(x=timeS, y = lmQ.bc.Z), color = "black") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 200)) +
  theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(filename = paste(plotdir,"baselineCorrected.lmQ.Z_sampled", ".pdf", sep=""), x8.1, width = 10, height = 8, units = 'cm')

# Alternatively, isolate normalizedF, and then do a baseline correction 
metric = gcamp.fitted3$normalizedF
metric[is.na(metric)]<-mean(metric[1:trim*1.01],na.rm=TRUE)
metric.bc = baseline(matrix(metric, nrow=1), method='irls')
# Plot with base R by opening a jpeg file, then close the file.
plot(metric.bc)
testit(7)
png(paste(plotdir, "baselineCorrected_compare_normalizedF.png", sep=""), width = 500, height = 300)
plot(metric.bc)
dev.off()

#get corrected data back into df
corrected = getCorrected(metric.bc)
#make sure the corrected vector length is same length as the df
length(corrected); length(gcamp.fitted3$normalizedF)
gcamp.fitted3["nF.bc"] = corrected[1,]
#plot 
x9 = ggplot() + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = nF.bc), color = "red") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 200)) +
  #guides(x = guide_prism_minor()) +
  ggtitle("lmQ baselineCorrected")
#x9
ggsave(filename = paste(plotdir,"baselineCorrected.nF", ".pdf", sep=""), x9, width = 15, height = 4, units = 'cm')


######################
# Plot everything 
#
######################
xRaw = ggplot() + 
  #geom_line(data=gcamp.fitted3, aes(x=timeS, y = Region0G)) + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = Region0G), color = "red", size = 0.1) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 200)) +
  ggtitle("Region0G")
#xRaw 

#reproducing sages final plot 
x0 = ggplot() + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = Region0G))  +
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = scaled415.470), color = "red") + 
  ggtitle("470Raw and scaled415.470")
#x0
ggsave(filename = paste(plotdir,"470Raw and scaled415.470", ".pdf", sep=""), x0, width = 15, height = 4, units = 'cm')
x1 = ggplot() + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = Region0G-scaled415.470))  +
  ggtitle("470Raw-scaled415.470")
#x1
ggsave(filename = paste(plotdir,"470Raw-scaled415.470", ".pdf", sep=""), x1, width = 15, height = 4, units = 'cm')

x2 = ggplot() + 
  #geom_line(data=gcamp.fitted3, aes(x=timeS, y = Region0G)) + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = lmQuotient), color = "red", size = 0.1) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 200)) +
  ggtitle("470 lmQuotient")
#x2
#ggplotly(x2)
ggsave(filename = paste(plotdir,"470 lmQuotient", ".pdf", sep=""), x2, width = 15, height = 4, units = 'cm')


#remember, can't plot raw 470 against normalized because values are divergent! 
x3 = ggplot() + 
  #geom_line(data=gcamp.fitted3, aes(x=timeS, y = Region0G))  +
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = normalizedF), color = "red") + 
  ggtitle("470Raw and normalizedF")
#x3
ggsave(filename = paste(plotdir,"normalizedF", ".pdf", sep=""), x3, width = 15, height = 4, units = 'cm')

x4 = ggplot() + 
  #geom_line(data=gcamp.fitted3, aes(x=timeS, y = Region0G))  +
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = deltaFq), color = "red") + 
  ggtitle("deltaFq")
#x4
ggsave(filename = paste(plotdir,"deltaFq", ".pdf", sep=""), x4, width = 15, height = 4, units = 'cm')

x5 = ggplot() + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = deltaFfit), color = "red") + 
  ggtitle("deltaFfit")
#x5
ggsave(filename = paste(plotdir,"deltaFfit", ".pdf", sep=""), x5, width = 15, height = 4, units = 'cm')

x6 = ggplot() + 
  #geom_line(data=head(gcamp.fitted3, n = 5000), aes(x=timeS, y = deltaFslide), color = "red") + 
  geom_line(data=gcamp.fitted3, aes(x=timeS, y = deltaFslide), color = "red") + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 200)) +
  #guides(x = guide_prism_minor()) +
  ggtitle("deltaFslide")
#x6
ggsave(filename = paste(plotdir,"deltaFslide", ".pdf", sep=""), x6, width = 15, height = 4, units = 'cm')

#extract a unique identifier from File.Path
identifier = str_extract(gcamp.fitted3$File.Path[1], "(?<=alldata)(.+)(?=\\.)")

# save gcamp.fitted2 or gcamp.fitted3 as a csv 
write.table(gcamp.fitted3, file=paste(plotdir, "NP_processed_",identifier, "_trim",trim,"_trim2",trim2, "_smooth",smooth, ".csv", sep=""), sep=",", row.names = FALSE)






