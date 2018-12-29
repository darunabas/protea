
library(raster)
library(tidyverse)


# load functions
source("code/functions.R")


# load protea occurrences and climate rasters
d <- read.csv("data/protea_records_cleaned.csv", stringsAsFactors = F)
temp <- stack("data/temperature.tif")
precip <- stack("data/precipitation.tif")


# calculate flowering year variables
d <- split(d, d$species) %>%
      lapply(function(x){
            fy <- flowering_year(x$julian, peak_mos=3, low_mos=6)
            x$peak <- fy$peak
            x$low <- fy$low
            x$aseasonality <- fy$aseasonality
            x$dofy <- fy$dofy
            return(x)
            }) %>%
      do.call("rbind", .) %>%
      mutate(anchor_month = peak)
      
# calculate spatial and temporal anomalies
coordinates(d) <- c("lon", "lat")
d <- split(d, d$species) %>%
      lapply(anomalies, climate=temp_rst, varname="temperature", lag=-8:3) %>%
      lapply(anomalies, climate=ppt_rst, varname="precipitation", lag=-8:3) %>%
      lapply(as.data.frame) %>%
      do.call("rbind", .) %>%
      na.omit() %>%
      group_by(species) %>%
      mutate(dofy_anomaly = dofy - mean(dofy))


## mixed effects model predicting flowering date

# full model. REML=F to allow model testing.
mfull <- lmer(dofy_anomaly ~ temperature_spatial_anomaly + precipitation_spatial_anomaly + temperature_temporal_anomaly + precipitation_temporal_anomaly + 
                    (0 + temperature_spatial_anomaly + precipitation_spatial_anomaly + temperature_temporal_anomaly + precipitation_temporal_anomaly|species), 
              data=d, REML=F)

# four reduced models
mtn <- lmer(dofy_anomaly ~ precipitation_spatial_anomaly + temperature_temporal_anomaly + precipitation_temporal_anomaly + 
                  (0 + precipitation_spatial_anomaly + temperature_temporal_anomaly + precipitation_temporal_anomaly|species), 
            data=d, REML=F)
mpn <- lmer(dofy_anomaly ~ temperature_spatial_anomaly + temperature_temporal_anomaly + precipitation_temporal_anomaly + 
                  (0 + temperature_spatial_anomaly + temperature_temporal_anomaly + precipitation_temporal_anomaly|species), 
            data=d, REML=F)
mta <- lmer(dofy_anomaly ~ temperature_spatial_anomaly + precipitation_spatial_anomaly + precipitation_temporal_anomaly + 
                  (0 + temperature_spatial_anomaly + precipitation_spatial_anomaly + precipitation_temporal_anomaly|species), 
            data=d, REML=F)
mpa <- lmer(dofy_anomaly ~ temperature_spatial_anomaly + precipitation_spatial_anomaly + temperature_temporal_anomaly + 
                  (0 + temperature_spatial_anomaly + precipitation_spatial_anomaly + temperature_temporal_anomaly|species), 
            data=d, REML=F)

# log likelihood significance testing
anova(mtn, mfull)
anova(mta, mfull)
anova(mpn, mfull)
anova(mpa, mfull)



