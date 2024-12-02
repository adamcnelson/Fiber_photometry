R script for analysis of fiber photometry, behavior, and body temperature
================


## Introduction

Fiber photometry is a way to record neuronal calcium events in animals during naturalistic behavioral and physiological responses. These three scripts serve three function. (1) Process fiber photometric data to yield different options of dF/F. (2) Align the dF/F data with behavior and body temperature data steams. (3) Aggregate multiple experiments for perform statistics.

## Installation
* One option is to install these scripts by cloning them with Git :

    For this, run the following commands :

    ```
    cd folder/to/clone/into
    git clone https://github.com/adamcnelson/Fiber_photometry
    ```

## Script 1: processing FP data from the Neurophotometrics system.
 
 Features:

-   Bleach detrending 
-   Data smoothing
-   Remove high-frequency noise 
-   Trim off fluorescence values from the first last/minutes of recording (user defined)
-   Baseline correction
-   Date-time timestamp (POSIXct)
-   Find peaks
-   Calculate dF/F using a number of different methods

### Example dataset
The following data was read in from the [Neurophotometrics system using a Bonsai workflow](https://neurophotometrics.com/bonsai-manual). Note that in this workflow the fluorescence data and associated timestamp values are read in as separate files.
Isosbestic channel: LedState == 1
GCaMP channel:  LedState == 2

``` r
all_dat
```
<img src="README_images/script1/Screenshot 2024-12-02 at 11.04.21 AM.png" width="550" />

``` r
computerClock
```
<img src="README_images/script1/Screenshot 2024-12-02 at 12.06.08 PM.png" width="300" />

### Smooth data 
The ButterEndEffect function removes high-frequency noise
``` r
ButterEndEffect
```
Before

<img src="README_images/script1/butterworth_before.png" width="400" />

After

<img src="README_images/script1/butterworth_after.png" width="400" />

Alternatively, use the `rollapply` function to smooth the data with a rolling average

### Convert fluorescence values at beginning / end of recording to `NA`
The beginning and end of a fiber photometry recording is very often noisy due to artifacts. The `trim_region0G` function converts calcium fluorescence values to NA for a user-defined number of minutes. 

``` r
trim_region0G <- function(data, led_state, trim, trim2) {
  data %>%
    dplyr::filter(LedState == led_state) %>%
    mutate(Region0G = ifelse(row_number() <= trim | row_number() >= trim2, NA, Region0G))
}
```

We often remove fluorescence from at least the first 10 minutes (~20,000 rows at 30 Hz). 

We can then align and visualize the calcium-dependent (GCaMP) and calcium-independent (isosbestic) signals using `ggplotGrob`
<img src="README_images/script1/fp_raw_Smoothing_lowPass0.133333333333333_.png" width="400" />

### Normalize data 
Here we present a few options for detrending the data and dealing with outlier values. 
#### Option 1: Correct the GCaMP data with an exponential decay model fit to the isosbestic data 
 -  Helpful tutorial from [Douglas Watson](https://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/)
 -  Use `nls`, the R base function to fit non-linear equations, and `SSasymp`, self-starting function that guesses its own start parameters.
 -  Use `na.exclude` to retain the original number of rows in the data 
 -  Use Broom's `augment` funciton to extract the predicted values. 
 
``` r
fit <- nls(Region0G ~ SSasymp(timeS, yf, y0, log_alpha), 
           data = isos, 
           na.action=na.exclude) #if it doesnt run with na.exclude, try na.omit

fitted = augment(fit, newdata=isos)
paraters = tidy(fit)
```
-   Linearly scale the fitted decay to the GCaMP (470) data using robust fit with a bisquare weighting function (rr.bisquare)
-   Then scale the fitted decay with coefficients from the robust fit (FP_lit_fit).
-   Finally, derive a normalized fluoresence value by dividing the GCaMP (470) values by the scaled fit (FP_lit_fit). 
``` r
rr.bisquare <- rlm(Region0G ~ .fitted, data=gcamp.fitted, psi = psi.bisquare, na.action = na.exclude)
```
``` r
gcamp.fitted2$FP.lin_fit = gcamp.fitted2$.fitted * rr.bisquare$coefficients[[2]] + rr.bisquare$coefficients[[1]]
```
``` r
gcamp.fitted2$normalizedF = gcamp.fitted2$Region0G/gcamp.fitted2$FP.lin_fit
```
- Plot the linearly scaled exponential decay fit over the GCaMP (470) data
<img src="README_images/script1/linearly scaled biexponential fit over 470 data.png" width="400" />

-   Plot the normalized data 
<img src="README_images/script1/normalizedF.png" width="400" />

------------------------------------------------------------------------

Authors: Adam Nelson
