R script for analysis of fiber photometry, behavior, and body temperature
================


## Introduction

Fiber photometry is a way to record neuronal calcium events in animals during naturalistic behavioral and physiological responses. These three scripts (1) process fiber photometric data to yield different options of DF/F, (2) align the DF/F data with behavior and body tempeature data steams, and (3) aggregate multiple experiments for perform statistics.

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
-   Remove high-frequency noise using Butterworth filter 
-   Trim off flouresence values from the beginning and end of the recording (user defined)
-   Baseline correction
-   Find peaks
-   Calculate dF/F using a number of different methods

### Downloading the test dataset
The following data was read in from the [Neurophotometrics system using a Bonsai workflow](https://neurophotometrics.com/bonsai-manual).

<img src="README_images/script1/Screenshot 2024-12-02 at 11.04.21 AM.jpeg" width="576" />
------------------------------------------------------------------------

Authors: Adam Nelson
