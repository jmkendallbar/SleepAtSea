# Sleep at Sea
---

This directory contains:

- ## [scripts](/scripts)
  Here, we provide code required to reproduce analyses specific to this paper: 
- ## [data](/data)
  Here, we provide the data required to reproduce figures & analyses specific to this paper.

## Sharing/Access information

This repository contains the data and code for our paper:

> Kendall-Bar, JM; Williams, TM; Mukherji, R; Lozano, DA; Pitman, JK; Holser, RR; Keates, T; Beltran, RS; Robinson, PW; Crocker, DE; Adachi, T; Lyamin, OI; Vyssotski, AL; Costa, DP (2023). *Brain activity of diving seals reveals short sleep cycles at depth*. Science.
> <https://doi.org/xxx/xxx>

### How to cite

Please cite this compendium as:

> Kendall-Bar, JM; Williams, TM; Mukherji, R; Lozano, DA; Pitman, JK; Holser, RR; Keates, T; Beltran, RS; Robinson, PW; Crocker, DE; Adachi, T; Lyamin, OI; Vyssotski, AL; Costa, DP (2023). *Compendium of R code and data for Brain activity of diving seals reveals short sleep cycles at depth*. Accessed 24 Oct 2022. Online at
> <https://doi.org/xxx/xxx>

## Contents

- ## [scripts](/scripts)
  Here, we provide code required to reproduce analyses specific to this paper: 
  
  <details>
  <summary> Read more. </summary> 

    Raw electrophysiological data were processed according to methods outlined in [Kendall-Bar et al., 2022](https://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-022-00287-x). Most notably, we use the [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) package to apply [independent component analysis](https://sccn.ucsd.edu/~arno/eeglab/auto/runica.html) to aid in the identification and removal of heart rate artifacts in the EEG signal. Intertial motion sensing data were processed according to methods outlined in [Cade et al. 2021](https://animalbiotelemetry.biomedcentral.com/articles/10.1186/s40317-021-00256-w) and tools freely available through the associated [CATS Toolbox Wiki](https://github.com/wgough/CATS-Methods-Materials). 

  </details>

  - ### **Step 00 - Speed Estimation:** `00_Speed-Estimation.m` - MATLAB code to estimate speed based on vertical speed + pitch, stroke rate, and body rotation. To be used in conjunction with the [CATS Toolbox](https://github.com/wgough/CATS-Methods-Materials) to create a dead-reckoned 3D track. 
    - **Input data:** 1Hz dataset with Depth (in meters), Pitch (in degrees), Roll (in degrees), and Stroke_Rate (in strokes per minute). `01_Hypnotrack_1Hz_ALL_ANIMALS.csv` (described below) is an example of an acceptable input data format.
    - **Output data:** 1Hz or interpolated 5Hz vector with speed estimate & plot to preview speed estimate output.
  - ### **Step 01 - Sleep Spiral Analysis:** `01_SleepSpiralAnalysis.m` - Code to identify sleep spirals and generate summary statistics regarding the prevalence of Slow Wave Sleep (SWS) and Rapid-Eye-Movement (REM) sleep within sleep spirals.
    - ### Part A: Identify Sleep Spirals
        - **Input data:** 1Hz dataset with sleep state and motion data (x, y, z [Depth], pitch, roll, heading). `01_Hypnotrack_1Hz_ALL_ANIMALS.csv` (described below) contains an acceptable input format, but the column *is_SleepSpiral* is generated through this step and is not required as an input.
        - **Output data:** 1Hz dataset with binary designation *is_SleepSpiral* and other summary tables including a list of *Naps*, *Loops*, and *Spirals* and associated data.
    - ### Part B: Summarize Sleep Spiral Data across Animals
        - **Input data:** `01_Hypnotrack_1Hz_ALL_ANIMALS.csv` - Sleep, motion, and sleep spiral data for all animals
        - **Output data:** Statistics on sleep spirals for the paper.
  - ### **Step 02 - Sleep Estimation & Summary** 
    - ### Part A: Sleep Estimation in MATLAB: `02A_Sleep-Estimates.m` 
      Code to estimate sleep using Time-Depth Records. Code generates plots and CSVs to summarize daily activity and sleep identification results.
      - **Input data:** Time-depth data with intervals of 8 or 10s (see example datasets below: `02_00_SLEEP_TOPPID_testNN_10_NewRaw.csv`)
      - **Output data:** Summary plots (for sleep data only) and CSVs to summarize diving/sleep activity. Examples show intermediate daily activity outputs for two adult females (correspond to Seals 1 & 2 from Manuscript's Figure 4)
        - `2011034_2036_Daily_Activity.csv` & `2011020_1015_Daily_Activity.csv`
    - ### Part B: Compile Results in R: `02B_CompileResults.Rmd`
      Compiles results from individual seals' CSV files into single spreadsheets for all seals.
      - **Input data:** References Adult Female metadata file (`02_00_AdultFemaleData_Metadata.csv`) & individual output files from previous step.
      - **Output data:** CSVs summarizing seal metadata, total dives, and daily activity.
        - `02_01_SleepEstimates_ALL_SealsUsed.csv`
        - `02_01_SleepEstimates_ALL_DailyActivity.csv`
    - ### Part C: Summary in R: `02C_Sleep-Estimate-Summary.Rmd`
      Code to summarize sleep estimation results 
      - **Input data:** Output data from previous step.
      - **Output data:** Two CSVs used to create maps in the paper:
        - `02_02_AdultFemaleData_DailySleepEstimates_long.csv`
        - `02_02_AdultFemaleData_DailySleepEstimates_wide_goodtracks.csv`
  - ### **Step 04 - Visualization**
    - ### Part A: Connect swim & glide attributes: `04A_swim_and_glide.py`
      Follow the steps outlined in [this tutorial](https://github.com/jmkendallbar/VisualizingLifeintheDeep/blob/main/00_3D_00_StrokeData_to_Swim.md) to animate a swim cycle and connect those animated keyframes to these swim and glide attributes.
    - ### Part B: Connect data to animation: `04B_InputData-ForAnimation.py`
      Run this code using PyMel and Maya commands in the python Script Editor in Autodesk Maya to control the animation of a 3D rigged animal with this data.
  - ### **MATLAB Functions**
        This repository uses a couple MATLAB functions, download the repo and add these to your path to allow MATLAB to reference these in the scripts above.


- ## [data](/data)
  Data used in the analysis. Large data files can be found in our Dryad repository (LINK).
    - ## *00 Metadata and Summary Data files:*
      - ### ***Results Summary*** - Wide-format table with metadata for EEG animals and summarizing all sleep data by location and sleep stage. 
        This table shows calculations for Total Sleep Time with and without putative REM sleep. 
        `00_00_Sleep-Results-Summary-Table.xlsx`

      - ### ***EEG Metadata*** - Wide-format table with EEG recording metadata.
        `00_00_Sleep-Recording-Metadata.xlsx`
        <details>
        <summary> Column descriptions. </summary>
          - **Deployment** - Chronological EEG deployment sequence (1-13)
          - **TestNumber** - Recording ID number
          - **SealID** - Unique identifier for each seal
          - **Recording.ID** - identifier combining the location (in the lab [CAPTIVE], in the wild [WILD], or translocated [XLOC]), age (in years [yr] or months [mo]), and age class (juvenile or weanling) of the seal
          - **TOPPID** - Unique ID to match to `00_Sleep-Results-Summary-Table.xlsx` and TOPP database ('20' stands for *Mirounga angustirostris*, next two digits represent year, next three digits represent deployment number per year).
          - **StartLogger_DateTime** - start date & time (format: 'YYYY-MM-DD HH:MM:SS') for the recording
          - **OnAnimal_DateTime** - date & time logger was attached to animal (as detected by ECG)
          - **Duration_OnAnimal** - Duration of recording in hours (after OnAnimal_DateTime)
          - **ChannelConfiguration** - Vector of Channel #s for Raw EDF files that correspond to the vector of channel names: [LEOG REOG LEMG REMG LEEG1 REEG2 LEEG3 REEG4]
        </details>

      - ### ***Adult Female Metadata*** 
        `02_00_AdultFemaleData_Metadata.csv` - Metadata table for adult female deployments
        <details>
        <summary> Column descriptions. </summary>
          - **TOPPID** - Unique ID per deployment matching those in the TOPP database
          - **Year** - Year of start of deployment
          - **Season** - Season (Post-Breeding or Post-Molt)
          - **TDR_QC** - Binary to designate whether the time-depth record was of sufficient quality to run the sleep identification model.
          - **Track_QC** - Binary to designate whether the track was of sufficient quality (and length) to visualize spatial sleep results in summary figure
        </details>

    - ## *01 Processed data files:*
      - ### ***Hypnograms*** - Processed Sleep Scoring (lab, wild, & at sea)
        `00_Hypnogram_30s_ALL_ANIMALS.csv` - Processed sleep scoring data for 30s epochs for all animals.
        <details>
        <summary> Column descriptions. </summary>

          - **timebins** - Time in R format for the beginning of the 30s epoch
          - **SealID** - unique identifier for each seal
          - **Recording.ID** - identifier combining the location (in the lab [CAPTIVE], in the wild [WILD], or translocated [XLOC]), age (in years [yr] or months [mo]), and age class (juvenile or weanling) of the seal
          - **ID** - in the lab [CAPTIVE], in the wild [WILD], or translocated [XLOC]
          - **Sleep.Code** - Specific sleep state designation: 
            - ***Active Waking***
            - ***Quiet Waking*** 
            - ***Drowsiness*** - Intermittent slow waves
            - ***LV Slow Wave SLeep*** - Low-voltage slow wave sleep
            - ***HV Slow Wave Sleep*** - High-voltage slow wave sleep
            - ***Certain REM Sleep*** - Rapid-Eye-Movement (REM) Sleep scored with high confidence (high degree of Heart Rate Variability [HRV])
            - ***Putative REM Sleep*** - REM Sleep scored with low confidence (low HRV)
            - ***Unscorable*** - Data not scorable due to interference, motion artifacts, or signal quality
          - **Simple.Sleep.Code** - Simplified sleep state designation: 
            - ***Active Waking***
            - ***Quiet Waking*** 
            - ***Drowsiness*** - Intermittent slow waves
            - ***SWS*** - Slow wave sleep (LV & HV combined)
            - ***REM*** - REM Sleep (certain and putative combined)
            - ***Unscorable*** - Data not scorable due to interference, motion artifacts, or signal quality
          - **Resp.Code** - Respiratory state designation:
            - ***Eupnea*** - between first breath and last breath
            - ***transition to Eupnea*** - transition to tachycardia
            - ***Apnea*** - between last breath and first breath
            - ***transition to Apnea*** - transition to bradycardia
            - ***Unscorable*** - not scorable due to noise obscuring HR detection
          - **Water.Code** - Location of animal
            - ***LAND*** - on land (in pen in the lab or on beach in the wild)
            - ***SHALLOW WATER*** - in water < 2m deep (in pool in the lab or in the lagoon at Ano Nuevo)
            - ***DEEP WATER*** - animal traversing the continental shelf (< 200 m / in water shallow enough that the animal can rest / travel along bottom)
            - ***OPEN OCEAN*** - animal in water deeper than 200 m / in water deep enough that the animal cannot rest / travel along bottom
          - **Time_s_per_day** - Time of day in seconds (out of 86400)
          - **Day** - Day of the recording
        </details>

      - ### ***Hypnotracks*** - Processed Sleep Scoring & Motion Data (3D tracks & Sleep State for seals at sea)
        `01_Hypnotrack_1Hz_ALL_ANIMALS.csv` - Timeseries data at 1Hz showing sleep state and processed motion data.
        <details>
        <summary> Column descriptions. </summary>
          
          *Note:* **Sleep_Num**, **Simple_Sleep_Num**, **Water_Num**, and **Resp_Num** redundantly code categorical/string data into numerical values for ease of analysis and plotting.

          - **Seconds** - Seconds elapsed for each recording
          - **R_Time** - Local time [PST] in R Format (YYYY-MM-DD HH:MM:SS)
          - **SealID** - unique identifier for each seal
          - **Recording_ID** - identifier combining the location (in the lab [CAPTIVE], in the wild [WILD], or translocated [XLOC]), age (in years [yr] or months [mo]), and age class (juvenile or weanling) of the seal
          - **ID** - in the lab [CAPTIVE], in the wild [WILD], or translocated [XLOC]
          - **Sleep_Code** - Specific sleep state designation: 
            - ***Active Waking***
            - ***Quiet Waking*** 
            - ***Drowsiness*** - Intermittent slow waves
            - ***LV Slow Wave SLeep*** - Low-voltage slow wave sleep
            - ***HV Slow Wave Sleep*** - High-voltage slow wave sleep
            - ***Certain REM Sleep*** - Rapid-Eye-Movement (REM) Sleep scored with high confidence (high degree of Heart Rate Variability [HRV])
            - ***Putative REM Sleep*** - REM Sleep scored with low confidence (low HRV)
            - ***Unscorable*** - Data not scorable due to interference, motion artifacts, or signal quality
          - **Simple_Sleep_Code** - Simplified sleep state designation: 
            - ***Active Waking***
            - ***Quiet Waking*** 
            - ***Drowsiness*** - Intermittent slow waves
            - ***SWS*** - Slow wave sleep (LV & HV combined)
            - ***REM*** - REM Sleep (certain and putative combined)
            - ***Unscorable*** - Data not scorable due to interference, motion artifacts, or signal quality
          - **Resp_Code** - Respiratory state designation:
            - ***Eupnea*** - between first breath and last breath
            - ***transition to Eupnea*** - transition to tachycardia
            - ***Apnea*** - between last breath and first breath
            - ***transition to Apnea*** - transition to bradycardia
            - ***Unscorable*** - not scorable due to noise obscuring HR detection
          - **Water_Code** - Location of animal
            - ***LAND*** - on land (in pen in the lab or on beach in the wild)
            - ***SHALLOW WATER*** - in water < 2m deep (in pool in the lab or in the lagoon at Ano Nuevo)
            - ***DEEP WATER*** - animal traversing the continental shelf (< 200 m / in water shallow enough that the animal can rest / travel along bottom)
            - ***OPEN OCEAN*** - animal in water deeper than 200 m / in water deep enough that the animal cannot rest / travel along bottom
          - **DN** - Matlab date number for local time
          - **pitch** - angle downward (-) or upward (+) in radians
          - **roll** - angle of roll to the right (+) or left (-) in radians
          - **heading** - angle of rotation to the left (counterclockwise: +) or right (clockwise: -) in radians
          - **x** - pseudotrack's x-position translation from the origin/start of deployment in meters
          - **y** - pseudotrack's y-position translation from the origin/start of deployment in meters
          - **z** - pseudotrack's z-position translation from the origin/start of deployment in meters
          - **geoX** - x position translation from the origin/start of deployment in meters using geo-referenced pseudotrack
          - **geoY** - y position translation from the origin/start of deployment in meters using geo-referenced pseudotrack
          - **Depth** - depth in meters (same as z except *(-1))
          - **speed** - estimated speed in meters per second
          - **Lat** - Latitude in Decimal Degrees
          - **Long** - Longitude in Decimal Degrees
          - **Stroke_Rate** - automated peak detection result for stroke frequency in strokes per minute
          - **Heart_Rate** - automated peak detection result for heart rate in beats per minute
          - **L_EEG_Delta** - Delta power (0.5-4Hz) for Left Hemisphere electroencephalogram (EEG)
          - **R_EEG_Delta** - Delta power (0.5-4Hz) for Right Hemisphere electroencephalogram (EEG)
          - **HR_VLF_Power** - Very low frequency (0-0.005 Hz) power for Heart Rate (quantification of HRV)
          - **is_SleepSpiral** - binary (0 or 1) showing which segments qualify as a sleep spiral (> 2 consecutive loops in the same direction during sleep)

        </details>
      
      - ### ***Hypnotrack Excerpt*** - 1Hz Excerpt from a sleep dive at sea (used in data visualization)
        `01_02_AnimationExcerpt_Hypnotrack_1Hz.csv` - Timeseries data at 1Hz with x y z positions and sleep data.
        <details>
        <summary> Column descriptions. </summary>

          - **Seconds** - Seconds elapsed for each recording
          - **R_Time** - Local time [PST] in R Format (YYYY-MM-DD HH:MM:SS)
          - **SealID** - unique identifier for each seal
          - **Recording_ID** - identifier combining the location (in the lab [CAPTIVE], in the wild [WILD], or translocated [XLOC]), age (in years [yr] or months [mo]), and age class (juvenile or weanling) of the seal
          - **ID** - in the lab [CAPTIVE], in the wild [WILD], or translocated [XLOC]
          - **Sleep_Code** - Specific sleep state designation: 
            - ***Active Waking***
            - ***Quiet Waking*** 
            - ***Drowsiness*** - Intermittent slow waves
            - ***LV Slow Wave SLeep*** - Low-voltage slow wave sleep
            - ***HV Slow Wave Sleep*** - High-voltage slow wave sleep
            - ***Certain REM Sleep*** - Rapid-Eye-Movement (REM) Sleep scored with high confidence (high degree of Heart Rate Variability [HRV])
            - ***Putative REM Sleep*** - REM Sleep scored with low confidence (low HRV)
            - ***Unscorable*** - Data not scorable due to interference, motion artifacts, or signal quality
          - **Simple_Sleep_Code** - Simplified sleep state designation: 
            - ***Active Waking***
            - ***Quiet Waking*** 
            - ***Drowsiness*** - Intermittent slow waves
            - ***SWS*** - Slow wave sleep (LV & HV combined)
            - ***REM*** - REM Sleep (certain and putative combined)
            - ***Unscorable*** - Data not scorable due to interference, motion artifacts, or signal quality
          - **Resp_Code** - Respiratory state designation:
            - ***Eupnea*** - between first breath and last breath
            - ***transition to Eupnea*** - transition to tachycardia
            - ***Apnea*** - between last breath and first breath
            - ***transition to Apnea*** - transition to bradycardia
            - ***Unscorable*** - not scorable due to noise obscuring HR detection
          - **Water_Code** - Location of animal
            - ***LAND*** - on land (in pen in the lab or on beach in the wild)
            - ***SHALLOW WATER*** - in water < 2m deep (in pool in the lab or in the lagoon at Ano Nuevo)
            - ***DEEP WATER*** - animal traversing the continental shelf (< 200 m / in water shallow enough that the animal can rest / travel along bottom)
            - ***OPEN OCEAN*** - animal in water deeper than 200 m / in water deep enough that the animal cannot rest / travel along bottom
          - **DN** - Matlab date number for local time
          - **pitch** - angle downward (-) or upward (+) in radians
          - **roll** - angle of roll to the right (+) or left (-) in radians
          - **heading** - angle of rotation to the left (counterclockwise: +) or right (clockwise: -) in radians
          - **x** - pseudotrack's x-position translation from the origin/start of deployment in meters
          - **y** - pseudotrack's y-position translation from the origin/start of deployment in meters
          - **z** - pseudotrack's z-position translation from the origin/start of deployment in meters
          - **geoX** - x position translation from the origin/start of deployment in meters using geo-referenced pseudotrack
          - **geoY** - y position translation from the origin/start of deployment in meters using geo-referenced pseudotrack
          - **Depth** - depth in meters (same as z except *(-1))
          - **speed** - estimated speed in meters per second
          - **Lat** - Latitude in Decimal Degrees
          - **Long** - Longitude in Decimal Degrees
          - **Stroke_Rate** - automated peak detection result for stroke frequency in strokes per minute
          - **Heart_Rate** - automated peak detection result for heart rate in beats per minute
          - **L_EEG_Delta** - Delta power (0.5-4Hz) for Left Hemisphere electroencephalogram (EEG)
          - **R_EEG_Delta** - Delta power (0.5-4Hz) for Right Hemisphere electroencephalogram (EEG)
          - **HR_VLF_Power** - Very low frequency (0-0.005 Hz) power for Heart Rate (quantification of HRV)

        </details>

      - ### ***Higher-resolution Rotation and Swimming Data*** - 10Hz Excerpt from a sleep dive at sea (rotation and heart rate data for visualization)
        `01_03_AnimationExcerpt_RotationSwim_10Hz.csv` - Timeseries data at 10Hz showing rotation and swimming behavior. Stroke rate and glide controller data were processed using the methods demonstrated in [Kendall-Bar et al. 2021](https://ieeexplore.ieee.org/document/9622956) through this [Github Wiki](https://github.com/jmkendallbar/VisualizingLifeintheDeep/blob/main/00_3D_00_StrokeData_to_Swim.md).  
        <details>
        <summary> Column descriptions. </summary>

          - **Seconds** - Seconds elapsed for each recording
          - **ECG** - Electrocardiogram (ECG) data in microvolts
          - **pitch** - angle downward (-) or upward (+) in degrees
          - **roll** - angle of roll to the right (+) or left (-) in degrees
          - **heading** - angle of rotation to the left (+) or right (-) in degrees
          - **GyrZ** - gyroscope data (angular acceleration) showing stroking behavior
          - **Glide_Controller** - glide controller for animation (1 for gliding 0 for stroking) and smoothed over a 5-second window
          - **Depth** - depth in meters
          - **Heart_Rate** - heart rate in beats per minute
          - **Stroke_Rate** - stroke rate in stroked per minute
          - **Heart_Detected** - binary (0 or 1- heartbeat detected)
          - **Stroke_Detected** - binary (0 or 1- stroke detected)

        </details>

      - ### ***EEG/ECG Excerpt for a Sleep Dive*** - Processed Sleep Scoring & Motion Data (3D tracks & Sleep State for seals at sea)
        `01_04_AnimationExcerpt_ECGEEGs_50Hz.csv` - Timeseries data at 50Hz with ECG, LEEG, and REEG during a sleeping dive.
        <details>
        <summary> Column descriptions. </summary>

          - **Seconds** - Time in R format for the beginning of the 30s epoch
          - **ECG** - Electrocardiogram (ECG) data (in microVolts)
          - **LEEG** - Electroencephalogram (EEG) data (in microVolts)
          - **REEG** - Electroencephalogram (EEG) data (in microVolts)

        </details>

      - ### ***Time-depth Data for Sleep Estimation*** - Processed Data for Sleep Estimation Script
        `02_00_SLEEP_TOPPID_testNN_10_NewRaw.csv` - Timeseries data at 8 second intervals with sleep scoring information, geographic locations, and time-depth data.
        <details>
        <summary> Column descriptions. </summary>

          - **Columns same as Hypnotrack**
          Plus the required columns for this script:
          - **time** - Matlab date number for time 
          - **CorrectedDepth** - Zero-offset corrected depth

        </details>
      
      - ### ***Sleep Estimation Script Ouput***
        `02_01_2011020_1015_Daily_Activity.csv` / `02_01_2011034_2036_Daily_Activity.csv` - Wide format observations of each day for two specific deployments. Activity labels are in individual columns (compared with long format above).
        <details>
        <summary> Column descriptions. </summary>

        - Columns:
          - **TOPPID** - unique identifier for each instrument deployment
          - **SEALID** - unique identifier for each seal
          - **unique_Days** - Matlab date number for the day of the observation
          - **Days_Elapsed** - Number of days into the trip
          - **Percent_of_Trip** - The percent of the trip that has elapsed 
          - ***Daily_recording*** - number of recording hours per day (should be 24 or close to it) 
          - ***Daily_diving*** - number of diving hours per day (time spent below 2 m)
          - ***Daily_long_SI*** - number of hours in an extended surface interval (at surface for > 10 min)
          - ***Daily_filtered_long_drift_long_SI*** - number of hours of estimated sleep (includes potential sleep while drifting, on the ocean floor, and in extended surface intervals)
          - ***Dailydive_long_glide*** - number of hours spent in an extended glide (to be roughly compared to estimated sleep (Daily_filtered_long_drift_long_SI)).
          - **Lat** - Latitude in decimal degrees
          - **Long** - Longitude in decimal degrees
          - **Lon360** - Longitude in decimal degrees from 0 to 360 (no negative values)

        </details>

      - ### ***Compiled Sleep Estimate Data*** 
        Intermediate outputs with summarized data across seals
        `02_01_SleepEstimates_ALL_SealsUsed.csv`
        <details>
        <summary> Column descriptions. </summary>

          - **TOPPID** - Unique deployment ID matching TOPP Database	
          - **haveStrokes** - Binary code designating the presence of stroke rate data
          -	**haveSleep**	- Binary code designating the presence of sleep data
          - **haveLatLong** - Binary code designating the presence of LatLong data
          	SIs_long	Drifts_long	Flats_long	Filtered_Drifts_long	Season_Code**
          Plus the required columns for this script:
          - **Dives** - Number of total dives (for each seal)
          - **SIs_long** - Number of extended surface intervals (>10 min)
          - **Flats_long** - Number of estimated naps on the sea floor
          - **Filtered_Drifts_long** - Number of estimated naps 
          - **Season_Code** - Season designation (PB- post-breeding or PM- post molt)

        </details>


      - ### ***Adult Female Sleep Estimate Data***
        `02_01_AdultFemaleData_DailySleepEstimates_long.csv` - Long format observations of hours per day performing different activities for Sleep Estimate analysis for all animals. Observations for the same seal (across multiple deployments) were later grouped and averaged for each trip percentile to compare individuals. Activity labels are in column `DailyActivity_label` and values are in column `h_per_day`.
        <details>
        <summary> Column descriptions. </summary>

        - Columns:
          - **TOPPID** - unique identifier for each instrument deployment
          - **SEALID** - unique identifier for each seal
          - **Season_Code** - PB [Post-breeding (short trip ~ 2 months)] or PM [Post-molt (long trip ~ 7 months)]
          - **triprecord_days** - total number of days for the deployment 
          - **deploys_per_seal** - number of deployments for the seal 
          - **unique_Days** - Matlab date number for the day of the observation
          - **Days_Elapsed** - Number of days into the trip
          - **Percent_of_Trip** - The percent of the trip that has elapsed 
          - **DailyActivity_label** - Labels for the type of observation in each row:
            - ***daily_recording*** - number of recording hours per day (should be 24 or close to it) 
            - ***daily_diving*** - number of diving hours per day (time spent below 2 m)
            - ***daily_long_SI*** - number of hours in an extended surface interval (at surface for > 10 min)
            - ***daily_filtered_long_drift_long_SI*** - number of hours of estimated sleep (includes potential sleep while drifting, on the ocean floor, and in extended surface intervals)
            - ***daily_long_glide*** - where available, number of hours of gliding as measured with a kami kami stroke rate logger
          - **h_per_day** - Hours per day (out of 24) for each of the activity labels listed above

        </details>

      - ### ***Adult Female Sleep Estimate Data for Geospatial Analysis***
        `02_02_AdultFemaleData_DailySleepEstimates_wide_goodtracks.csv` - Wide format observations of each day for deployments with tracks of sufficient quality to be included in geospatial sleep analysis. Activity labels are now in individual columns (compared with long format above).
        <details>
        <summary> Column descriptions. </summary>

        - Columns:
          - **TOPPID** - unique identifier for each instrument deployment
          - **SEALID** - unique identifier for each seal
          - **Season_Code** - PB [Post-breeding (short trip ~ 2 months)] or PM [Post-molt (long trip ~ 7 months)]
          - **unique_Days** - Matlab date number for the day of the observation
          - **Days_Elapsed** - Number of days into the trip
          - **Percent_of_Trip** - The percent of the trip that has elapsed 
          - ***Daily_recording*** - number of recording hours per day (should be 24 or close to it) 
          - ***Daily_diving*** - number of diving hours per day (time spent below 2 m)
          - ***Daily_long_SI*** - number of hours in an extended surface interval (at surface for > 10 min)
          - ***Daily_filtered_long_drift_long_SI*** - number of hours of estimated sleep (includes potential sleep while drifting, on the ocean floor, and in extended surface intervals)
          - **Lat** - Latitude in decimal degrees
          - **Long** - Longitude in decimal degrees
          - **Lon360** - Longitude in decimal degrees from 0 to 360 (no negative values)

        </details>

    - ## *Raw data files:*
      - ### ***EEG Raw Data***
        - Raw EEG, EMG, EOG, and motion sensor data for all deployments, labeled by *TestNumber* (see `00_Sleep-Recording-Metadata.xlsx` for related metadata including start times and channel configuration details).

----

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/) Please provide attribution (Jessica Kendall-Bar, 2023) when using our figures.

**Code :** Please cite the Zenodo DOI provided above when using our code.

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse. Please cite the Dryad DOI provided above when using our code.
