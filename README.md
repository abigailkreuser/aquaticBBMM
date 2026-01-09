# aquaticBBMM

This is the code used to create occurrence distributions from intermittent sightings of individuals in a manuscript in review. Right whale data can be requested from the North Atlantic Right Whale Consortium https://www.narwc.org/. I will try to create a dummy dataset or publicly available data in the future for examples. I used a lot more code for visualizations and summarizing the data, but wanted this repository to focus on the model and outputs itself. I'm interested in making this a package and more user friendly one day, so updates are likely to come! Feel free to reach out with any questions. 


Code for this model was written in R Statistical Software (v4.2.2; R Core Team 2022) was adapted from source code in the R packages move (v4.1.6; Kranstauber et al. 2021) and BBMM (v3.0; Nielsen et al. 2013).
The network which is essential for routing paths around land was inspired by inspired by the pathroutr package (v0.1.1-beta, London 2020)


## Data preparation

1. `20220427_data_pull_n_shape.R`
Links pairs of sightings to create track segments which are then used to create the csv files below. 
  - creates `mvmnt` which holds most meta data for the tracks. 
    - saved as: `20250203mvmnt_withAreasLats_jDates.csv`
  - creates .csv files in `*date*_indvTrack_csv` folder that needs to be copied over to the super computer and read into the bbmm scripts. 
    -saved as: `2025.02.04_indvTrack_csv`
  - Single location variant of this script is: `2025.01.22_singleSightingsNOTinTracks.R` 
    - also includes some exploratory figures right now

2. `20240331_MoveOrStay.R`
  - assigns area and latitudinal zones to the sight_1 and sight_2 locations
  - complete in the middle of the `pull_n_shape` script. 
  - ordinal dates and seasonal ordinal dates are also assigned here. (seasonal dates Nov1 starts the year)
  
3. Copy `date_indv_csv` to super computer's `data/processedData/`. 

## Network 
Again data is not sharable here to recreate, and not be possible for data-poor study systems. The network could be supplemented with opportunistic detections of the species of interest or similar species (e.g., inaturalist, gbif, movebank) or simulated using habitat constraints. Unstructured mesh grids or networks are often used to model physical oceanography dynamics, and nodes from pre-existing grids of study areas could be integrated with known locations to create a biologically relevant network.

`2024.10.16_NetworkAroundLand_plusFigures.R`

## Model scripts and calculations

1. Motion variance script  
`2024.06.18_FINALmotionVar_with_Network_4Hyperion.R`
- need to update track path to the indv csv's of tracks 
- update to the mvmnt csv file with > 2 locations. 
`20250203_mvmnt_with_n_greaterthan2.csv`

2. Run Model Scripts 
  - bbmm for all tracks 
`2025.09.26_bbmmDaily_steps4short8long_res4km_gridDist1c_part1.R`
- bbmm for single locations 
`2025.02.05_bbmm_singleLocs_Daily_steps4short8long_res4km_gridDist1c.R`

3. Run Whos missing scripts 
`2024.08.11_missingTrackIDs.R`

4. copy files for the single locations into the larger folder 

5. calculate the daily rasters for every single year in model  
`2024.06.22_yearBYyearJUSTdailyRasters_1c_part1.R`
- a list is created for the files corresponding to a given day, it only needs to be made once and then can be used for other years to reference. 

6. calculate the monthly rasters for every time period  
`2025.03.06_decadalRastersMonthly.R`
- a list is created for the files corresponding to a given month, it only needs to be made once and then can be used for other years to reference. 
- these outputs are used to create seasonal rasters

7. Calculating daily sum whale days or occurence probabilities of the modeled indv
`2025.02.08_HabUse_yearlyPop_NOindv_singlesAdded.R`






