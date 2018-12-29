
# sliding window calculations of peak and low flowering seasons, seasonality, and day of flowering year
flowering_year <- function(j, # vectory of julian days in range 1:365
                           peak_mos=3, # length of peak flowering season, in months
                           low_mos=6 # length of low flowering season, in months
){
      
      # wrapped sliding window observation counts for every day of year
      counts_peak <- c()
      counts_low <- c()
      ji <- j
      for(i in 1:365){
            counts_peak <- rbind(counts_peak, sum(ji %in% 1:(30*peak_mos)))
            counts_low <- rbind(counts_low, sum(ji %in% 1:(30*low_mos)))
            ji <- ji - 1
            ji[ji==0] <- 365
      }
      
      # peak -- center of window with max observations
      peak <- which(counts_peak==max(counts_peak))
      if(length(peak)>1) warning("Multiple windows tied for peak flowering -- returning the earliest one")
      peak <- peak[1] + peak_mos * 30 / 2
      peak[peak>365] <- peak[peak>365] - 365
      
      # low -- center of window with min observations
      low <- which(min(counts_low)==counts_low)
      if(length(low)>1) warning("Multiple windows tied for low flowering -- returning the earliest one")
      low <- low[1] + low_mos * 30 / 2
      low[low>365] <- low[low>365] - 365
      
      # normalized version of julian day based on low season
      dofy <- j - low
      dofy[dofy < 1] <-dofy[dofy < 1] + 365
      
      # index of how restricted flowering is
      aseasonality <- min(counts_low) / max(counts_peak)
      
      list(peak=peak, # julian day of the center of peak flowering season
           low=low, # julian day of the center of low flowering season
           aseasonality=aseasonality, # aseasonality index, lower values are higher seasonality
           dofy=dofy) # day of flowrering year, i.e. julian days after center of low flowering season
}


# This function summarizes climate for a portion of the focal year trailing a focal date
# Its odd API is designed to be used in conjunction with raster::calc
window_anomaly <- function(x, # a vector of the focal year, focal julian day, and then followed by climate time series 
                           lag=0, # vector of month offsets. 0 for the focal month, 1:3 for the three mos preceeding, etc
                           fun=mean,
                           start_year=1950, # these must match the climate data date range
                           end_year=2010
){
      # sort components
      if(is.na(x[3])) return(NA)
      y <- x[1]
      j <- x[2]
      x <- x[3:length(x)]
      if(y < start_year+1 | y > end_year-1) return(NA) # buffer needed so offsets can spill into adjoining years
      
      # monthly means
      x <- matrix(x, ncol=length(x)/12)
      m <- apply(x, 1, mean)
      
      # anomalies
      x <- x[,y-start_year+1 + c(-1,0,1)] # isolate data from focal year and neighbors
      x <- x - cbind(m, m, m) # compute anomaly
      x <- as.vector(x) # 36-month window, needed so offsets can spill into adjoining years
      
      # isolate values for focal months
      j <- floor(j/365*12)
      x <- x[j-lag+12]
      fun(x)
}


# derive spatial and temporal climate anomalies for each occurrence
anomalies <- function(points, # SpatialPointsDataFrame of occurrences for one species
                      climate, # year-by-month climate time series raster stack
                      varname, # label to identify climate variable in output
                      ... # arguments passed to window_anomaly
){
      require(raster)
      
      # extract climate time series for each occurrence point
      clim <- extract(climate, points)
      
      # calculate raw temporal anomaly using window_anomaly function
      points$temporal_anomaly <- apply(cbind(points$year, points$anchor_month, clim), 1, FUN=window_anomaly, ...)
      
      # center spatial and temporal anomalies on zero
      mean_clim <- apply(clim, 1, FUN=mean, na.rm=T)
      if(varname=="precipitation") mean_clim <- log10(mean_clim)
      demean <- function(x) x - mean(x, na.rm=T)
      points$spatial_anomaly <- demean(mean_clim)
      points$temporal_anomaly <- demean(points$temporal_anomaly)
      
      # append variable name to new columns
      names(points)[names(points)=="spatial_anomaly"] <- paste0(varname, "_spatial_anomaly")
      names(points)[names(points)=="temporal_anomaly"] <- paste0(varname, "_temporal_anomaly")
      
      # return the input SpatialPointsDataFrame, with two new columns
      return(points) 
}
