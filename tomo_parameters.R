# ----------------------------------------------------------------------------------------
# General parameters
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# Initialise variables
# ----------------------------------------------------------------------------------------
tomo  <- list()
# Create domain and prior "objects" with default parameters
# tomo$domain <- get_domain_param(print_out = FALSE)
# tomo$prior  <- get_prior_param()

tomo$minstep <- 5
tomo$utils$ncores <- ncores
tomo$utils$results_directory <- results_directory
tomo$utils$data_directory    <- data_directory
tomo$utils$SAVE_RESULTS <- SAVE_RESULTS

# Additional_label is added to all filenames
tomo$utils$label <- label

# ----------------------------------------------------------------------------------------
# Domain
# ----------------------------------------------------------------------------------------
tomo$domain$alt_hi               <- 1300 * 10^3
tomo$domain$alt_lo               <- 0 * 10^3
tomo$domain$alt_res          	 <- 50 * 10^3
tomo$domain$alt_res_bndry_lo 	 <- 50 * 10^3
tomo$domain$alt_res_bndry_hi 	 <- 100 * 10^3
tomo$domain$alt_bndry_hi     	 <- 600 * 10^3
tomo$domain$long_west 			 <- -25
tomo$domain$long_east 			 <- 70
tomo$domain$lat_north 			 <- 85
tomo$domain$lat_south 			 <- 35
tomo$domain$long_bndry_west 	 <- tomo$domain$long_west + 2
tomo$domain$long_bndry_east 	 <- tomo$domain$long_east - 2
tomo$domain$long_res_bndry  	 <- 10
tomo$domain$long_res        	 <- 1
tomo$domain$lat_bndry_north 	 <- tomo$domain$lat_north - 2
tomo$domain$lat_bndry_south 	 <- tomo$domain$lat_south + 8
tomo$domain$lat_res_bndry   	 <- 2
tomo$domain$lat_res         	 <- 1

tomo$domain$limit_lat_north 	 <- tomo$domain$lat_north - 1
tomo$domain$limit_lat_south 	 <- tomo$domain$lat_south + 10
tomo$domain$limit_long_east 	 <- tomo$domain$long_east - 6
tomo$domain$limit_long_west 	 <- tomo$domain$long_west + 6
tomo$domain$lat_breaks           <- NULL
tomo$domain$long_breaks          <- NULL
tomo$domain$alt_breaks           <- NULL

# Below overrides resolutions set above
tomo$domain$lat_breaks   <- c(tomo$domain$lat_south, 50, 58, 70, 75, tomo$domain$lat_north)
tomo$domain$lat_res_vec  <- c(5, 1, 1, 1, 5)

tomo$domain$long_breaks  <- c(tomo$domain$long_west, 5, 10, 35, 40, tomo$domain$long_east)
tomo$domain$long_res_vec <- c(5, 1, 1, 1, 5)

tomo$domain$alt_breaks   <- c(tomo$domain$alt_lo, 50 * 10^3, 400 * 10^3, 600 * 10^3, tomo$domain$alt_hi)
tomo$domain$alt_res_vec   <- c(25, 25, 50, 100) * 10^3

# Distance (m) between points in raycasting approximation
tomo$domain$path_resolution <- 1000

tomo$domain$earth_rad <- 6371 * 10^3

# ----------------------------------------------------------------------------------------
# Beacon data
# ----------------------------------------------------------------------------------------
tomo$LEO$USE             <- FALSE
tomo$LEO$directory       <- LEO_directory
tomo$LEO$stations        <- c("SOD", "TAR", "TRO", "SVA", "KIR", "KUU", "IVA", "OUL", "KAU", "VES", "NUR", "JOM", "KEV", "LYC")
#tomo$LEO$satellites      <- c("CASSIOPE", "COSMOS_2463", "COSMOS_2407", "DMSP_5D-3_F15_(USA_147)")
tomo$LEO$satellites      <- c("COSMOS_2463", "COSMOS_2407")
tomo$LEO$averaging_t     <- 10 
tomo$LEO$limit_elevation <- tomo$domain$limit_elevation_LEO <- 10
tomo$LEO$modelling_error <- 2
tomo$LEO$level_to_zero   <- TRUE

# ----------------------------------------------------------------------------------------
# GNSS DATA
# ----------------------------------------------------------------------------------------
tomo$GNSS$directory       <- GNSS_directory
#tomo$GNSS$directory       <- paste(tomo$utils$data_directory, "gnss_data",  sep = "")
tomo$GNSS$averaging_t     <- tomo$minstep * 60
tomo$GNSS$average_std     <- FALSE
tomo$GNSS$rm_negat        <- FALSE # Only for precalibrated data
tomo$GNSS$limit_elevation <- 10
# Limit measurements used for calibration
tomo$GNSS$limit_elevation_calibration <- 10
tomo$GNSS$limit_pp_latitude_calibration <- NULL#c(40, 75)
tomo$GNSS$limit_pp_longitude_calibration <- NULL#c(-20, 60)
# Modelling error is (1 / sin(el / 180 * pi)) * tomo$GNSS$modelling_error
tomo$GNSS$modelling_error <- 0.1
tomo$GNSS$limit_tstart_calibration <- NULL#as.POSIXct("2018-11-09 06:00:00", tz = "UTC")
tomo$GNSS$limit_tstop_calibration <-  NULL#as.POSIXct("2018-11-09 14:00:00", tz = "UTC")

# ----------------------------------------------------------------------------------------
# RADIO OCCULTATION DATA
# ----------------------------------------------------------------------------------------
tomo$RO$USE             <- FALSE
tomo$RO$directory       <- RO_directory
tomo$RO$averaging_t     <- 30
tomo$RO$average_std     <- TRUE
tomo$RO$rm_negat        <- FALSE # Only for precalibrated data
tomo$RO$modelling_error <- tomo$GNSS$modelling_error
tomo$RO$limit_elevation <- -90
tomo$RO$LEO_remove      <- NULL
tomo$RO$GNSS_remove     <- NULL

# ----------------------------------------------------------------------------------------
# Ionosonde data
# ----------------------------------------------------------------------------------------
#tomo$IONOSONDE <-  list()
tomo$IONOSONDE$USE_AS_DIRECT_MEASUREMENTS <- TRUE
tomo$IONOSONDE$directory <- IONOSONDE_directory
# Stations read and used for validation  
tomo$IONOSONDE$station_all <- c("TR", "JR055")
# Stations in use if USE_AS_DIRECT_MEASUREMENTS = TRUE
tomo$IONOSONDE$station <- tomo$IONOSONDE$station_all#c("TR", "JR")
tomo$IONOSONDE$times  <- 15 #minutes backwards from t1
tomo$IONOSONDE$times_ahead <- 5
# tomo$IONOSONDE$paths <- IONOSONDE_paths
# tomo$IONOSONDE$station_lat  <- c(69.60, 54.6, 55.47, 51.7, 51.7, 50.0, 50.1, 40.8, 41.9, 40.6, 38.0)
# tomo$IONOSONDE$station_long <- c(19.2, 13.4, 37.3, -1.5, 14.6, 16.72, 4.6, 12.5, 17.8, 23.5)
# tomo$IONOSONDE$station_alt  <- c(86, 50, 50, 50, 50, 50, 50, 50, 50, 50)
tomo$IONOSONDE$sd_coeff <- 0.2
tomo$IONOSONDE$outlier_ne_limit_abs <- 20 * 10^11

tomo$IONOSONDE$topside <- TRUE
tomo$IONOSONDE$alt_limit <- 1000 * 10^3


# ----------------------------------------------------------------------------------------
# EISCAT data
# ----------------------------------------------------------------------------------------
tomo$ISR$USE_AS_DIRECT_MEASUREMENTS <- FALSE
# Stations in use if USE_AS_DIRECT_MEASUREMENTS = TRUE
tomo$ISR$station <- NA 
# Stations read and used for validation
tomo$ISR$station_all <- c("uhf", "vhf", "32m", "42m")
tomo$ISR$directory <- ISR_directory
tomo$ISR$paths <- ISR_paths
tomo$ISR$sd_coeff <- 0.1
tomo$ISR$alt_limit <- 350 * 10^3
tomo$ISR$range_limit <- 500 * 10^3

# ----------------------------------------------------------------------------------------
# NeQuick model parameters
# ----------------------------------------------------------------------------------------
tomo$prior$NeQuick$USE <- FALSE
tomo$prior$NeQuick$path <- NeQuick_directory
tomo$prior$NeQuick$flux <- 140 

# ----------------------------------------------------------------------------------------
# Background inosonde parameters
# ----------------------------------------------------------------------------------------
tomo$prior$ionosonde_avg_bg$USE <- FALSE
tomo$prior$ionosonde_avg_bg$profile <- tomo$IONOSONDE$station_all

# ----------------------------------------------------------------------------------------
# Precalculated tomoscand background
# ----------------------------------------------------------------------------------------
tomo$prior$tomoscand_bg$USE <- FALSE
tomo$prior$tomoscand_bg$path <- tomoscand_bg_path

# ----------------------------------------------------------------------------------------
# Prior parameters
# ----------------------------------------------------------------------------------------
tomo$prior$EF_boundary_alt <- 160 * 10^3

# Update prior variance according to available ionosonde measurements
# CALIBRATE = TRUE will override DYNAMIC_PRIOR_VAR
tomo$prior$DYNAMIC_PRIOR_VAR <- TRUE

tomo$prior$dynamic_profiles <- tomo$IONOSONDE$station_all

# Use previous posterior mean as prior
# Prior mean is multiplied with "mean_scale"
tomo$prior$RECURSIVE_PRIOR <- TRUE
tomo$prior$mean_parameters$mean_scale <- 0.9
# If previous data is older than t_max_mean
# The mean is then automatically recalibrated
tomo$prior$mean_parameters$t_max_diff <- 60


# When no ionosonnde measurements is available
# The dynamic prior will shift towards
# ionospheric_F_region_altitude_day   <- 220 * 10^3
# ionospheric_F_region_altitude_night <- 300 * 10^3
# Day and night is defined as
# day_starts <- 6 #UT
# day_ends   <- 15 #UT

# Prior profile
# Zero or ionosonde with exponential topside
# RECURSIVE_PRIOR = TRUE overrides and uses previous posterior mean
# tomo$prior$mean_parameters$type <- "ZERO"
tomo$prior$mean_parameters$type <- "ZERO"
tomo$prior$mean_parameters$peak_ne_offset_coeff <- 0
tomo$prior$mean_parameters$lat_trend_south  <- 1
tomo$prior$mean_parameters$lat_trend_north  <- 1
tomo$prior$mean_parameters$long_trend_west  <- 1
tomo$prior$mean_parameters$long_trend_east  <- 1

# Initial values if not taken from ionosonde
tomo$prior$peak_ne_F_default    <- 1 * 10^(-5)
tomo$prior$peak_alt_F_default   <- 300 * 10^3
tomo$prior$bottom_alt_F_default <- 180 * 10^3

tomo$prior$TEMPORAL_PRIOR_MODEL <- FALSE
tomo$prior$peak_temporal_profile <- NULL

# Uniform plasma parameters for GNSS measurements
#tomo$prior$plasma$mean <- 1 / (21 * 10^6)
#tomo$prior$plasma$var  <- 2#(1^2 / (21 * 10^6))
tomo$prior$plasma$mean <- 0
tomo$prior$plasma$var  <- 10^(-100)

tomo$prior$mean_parameters$scale_hi <- 100 * 10^3
tomo$prior$mean_parameters$scale_lo <- 60 * 10^3
tomo$prior$mean_parameters$hilevel  <- 0


# Prior covariance (F layer)

# Prior sd profile type if 
# TRUE: peak is extended tomo$prior$peak_alt_default to tomo$prior$peak_alt altitude
# FALSE: (default) peak is moved to tomo$prior$peak_alt altitude
tomo$prior$covariance_parameters$dynamic_width <- FALSE

tomo$prior$covariance_parameters$power_perc      <- 0.1
tomo$prior$covariance_parameters$sd_peak_width   <- 0 * 10^3
tomo$prior$covariance_parameters$peak_alt_offset <- 0 * 10^3
tomo$prior$covariance_parameters$scaleH          <- 120 * 10^3
tomo$prior$covariance_parameters$scaleL          <- 100 * 10^3
tomo$prior$covariance_parameters$sd_hilevel      <- 0.05
tomo$prior$covariance_parameters$lat_trend_south <- 1
tomo$prior$covariance_parameters$lat_trend_north <- 1
tomo$prior$covariance_parameters$long_trend_west <- 1
tomo$prior$covariance_parameters$long_trend_east <- 1

tomo$prior$covariance_parameters$l_vert 	     <- 50 * 10^3
tomo$prior$covariance_parameters$l_hor  	     <- 2.5 
tomo$prior$covariance_parameters$l_long 	     <- 2.5 
tomo$prior$covariance_parameters$a_scale         <- 0.2014112

tomo$prior$covariance_parameters$l_vert 	     <- 50 * 10^3
tomo$prior$covariance_parameters$l_hor  	     <- 5
tomo$prior$covariance_parameters$l_long 	     <- 5
tomo$prior$covariance_parameters$a_scale         <- 0.2085711

# tomo$prior$covariance_parameters$l_vert          <- 100 * 10^3
# tomo$prior$covariance_parameters$l_hor           <- 10
# tomo$prior$covariance_parameters$l_long          <- 10
# tomo$prior$covariance_parameters$a_scale         <- 0.1712219


# tomo$prior$covariance_parameters$l_vert 		 <- 500 * 10^3
# tomo$prior$covariance_parameters$l_hor     		 <- 40
# tomo$prior$covariance_parameters$l_long 		 <- 40
# tomo$prior$covariance_parameters$a_scale   		 <- 0.05523687


# tomo$prior$covariance_parameters$l_vert 	     <- 400 * 10^3
# tomo$prior$covariance_parameters$l_hor  	     <- 50 
# tomo$prior$covariance_parameters$l_long 	     <- 50 
# tomo$prior$covariance_parameters$a_scale         <- 0.2014112


# tomo$prior$covariance_parameters$l_vert 	       <- 150 * 10^3
# tomo$prior$covariance_parameters$l_hor  	       <- 5
# tomo$prior$covariance_parameters$l_long 	       <- 5
# tomo$prior$covariance_parameters$a_scale         <- 0.1407634

# tomo$prior$covariance_parameters$l_vert 	       <- 100 * 10^3
# tomo$prior$covariance_parameters$l_hor  	       <- 5
# tomo$prior$covariance_parameters$l_long 	       <- 5
# tomo$prior$covariance_parameters$a_scale         <- 0.1977257



# Change prior parameters for forenoon
tomo$prior$solar$SOLAR_VARIATION <- FALSE
tomo$prior$solar$lat <- 60
tomo$prior$solar$long <- 25
tomo$prior$solar$dawn_transition <- 30
tomo$prior$solar$dusk_transition <- 120  
tomo$prior$solar$mean_scale <- 0.9

# tomo$prior$covariance_parameters$l_vert_sunrise      <- 250 * 10^3
# tomo$prior$covariance_parameters$l_hor_sunrise       <- 20 #800 / 111
# tomo$prior$covariance_parameters$l_long_sunrise      <- 20 #1000 / 40 
# tomo$prior$covariance_parameters$a_scale_sunrise     <- 0.1995142
# tomo$prior$covariance_parameters$E_layer$use_sunrise <- FALSE

# tomo$prior$covariance_parameters$l_vert_sunrise      <- 250 * 10^3
# tomo$prior$covariance_parameters$l_hor_sunrise       <- 10 #800 / 111
# tomo$prior$covariance_parameters$l_long_sunrise      <- 10 #1000 / 40 
# tomo$prior$covariance_parameters$a_scale_sunrise     <- 0.09531429
# tomo$prior$covariance_parameters$E_layer$use_sunrise <- FALSE

# tomo$prior$covariance_parameters$l_vert_sunrise      <- 150 * 10^3
# tomo$prior$covariance_parameters$l_hor_sunrise       <- 10 #800 / 111
# tomo$prior$covariance_parameters$l_long_sunrise      <- 10 #1000 / 40 
# tomo$prior$covariance_parameters$a_scale_sunrise     <- 0.1391242
# tomo$prior$covariance_parameters$E_layer$use_sunrise <- FALSE

# tomo$prior$covariance_parameters$l_vert_sunrise      <- 100 * 10^3
# tomo$prior$covariance_parameters$l_hor_sunrise       <- 10 #800 / 111
# tomo$prior$covariance_parameters$l_long_sunrise      <- 10 #1000 / 40 
# tomo$prior$covariance_parameters$a_scale_sunrise     <- 0.1712219
# tomo$prior$covariance_parameters$E_layer$use_sunrise <- FALSE

# tomo$prior$covariance_parameters$l_vert_sunrise 	   <- 150 * 10^3
# tomo$prior$covariance_parameters$l_hor_sunrise  	   <- 5
# tomo$prior$covariance_parameters$l_long_sunrise  	   <- 5 #800 / 111
# tomo$prior$covariance_parameters$a_scale_sunrise     <- 0.1407634
# tomo$prior$covariance_parameters$E_layer$use_sunrise <- FALSE

# tomo$prior$covariance_parameters$l_vert_sunrise      <- 50 * 10^3
# tomo$prior$covariance_parameters$l_hor_sunrise       <- 2.5 #800 / 111
# tomo$prior$covariance_parameters$l_long_sunrise      <- 2.5 #1000 / 40 
# tomo$prior$covariance_parameters$a_scale_sunrise     <- 0.2014112
# tomo$prior$covariance_parameters$E_layer$use_sunrise <- FALSE

tomo$prior$covariance_parameters$l_hor_sunrise  	 <- tomo$prior$covariance_parameters$l_hor  
tomo$prior$covariance_parameters$l_long_sunrise  	 <- tomo$prior$covariance_parameters$l_long
tomo$prior$covariance_parameters$l_vert_sunrise 	 <- tomo$prior$covariance_parameters$l_vert 
tomo$prior$covariance_parameters$a_scale_sunrise     <- tomo$prior$covariance_parameters$a_scale
tomo$prior$covariance_parameters$E_layer$use_sunrise <- FALSE



# Prior covariance (E layer)
# Initial values if not taken from ionosonde
tomo$prior$peak_ne_E_default  <- 0.1 * 10^(-5)
tomo$prior$peak_alt_E_default <- 100 * 10^3
tomo$prior$covariance_parameters$E_layer$power_perc <- 0.1

tomo$prior$covariance_parameters$E_layer$a_scale <- 0.1265157

tomo$prior$covariance_parameters$E_layer$use           <- TRUE
tomo$prior$covariance_parameters$E_layer$peak_ne       <- tomo$prior$peak_ne_E_default
tomo$prior$covariance_parameters$E_layer$peak_alt      <- tomo$prior$peak_alt_E_default
tomo$prior$covariance_parameters$E_layer$scaleH        <- 60 * 10^3
tomo$prior$covariance_parameters$E_layer$scaleL        <- 40 * 10^3
tomo$prior$covariance_parameters$E_layer$sd_hilevel    <- 0
tomo$prior$covariance_parameters$E_layer$sd_peak_width <- 0
tomo$prior$covariance_parameters$E_layer$lat_mid       <- 69.600 
tomo$prior$covariance_parameters$E_layer$lat_spread    <- 8
tomo$prior$covariance_parameters$E_layer$long_mid      <- 19.200
tomo$prior$covariance_parameters$E_layer$long_spread   <- 100

# ----------------------------------------------------------------------------------------
# Bias and calibration
# ----------------------------------------------------------------------------------------

# When tomo$prior$bias$UPDATE = TRUE, all the station biases will be updated
# When tomo$prior$bias$UPDATE = FALSE, only new bias estimates are added
tomo$prior$bias$UPDATE <- FALSE
tomo$prior$bias$SHIFT_POSITIVE <- TRUE
# CALIBRATE = TRUE will override USE_SAVED_BIAS and UPDATE_BIAS

# If latest data is older than t_max_mean
# The mean is then automatically recalibrated
tomo$prior$bias$t_max_diff <- 60 * 24

# Re-estimate station biases
#CALIBRATION                             <- TRUE
tomo$prior$calibration$profile                     <- tomo$IONOSONDE$station_all
tomo$prior$calibration$E_layer_use                 <- FALSE
tomo$prior$calibration$sd_hilevel                  <- 0
tomo$prior$calibration$lat_trend_south             <- 1.5
tomo$prior$calibration$lat_trend_north             <- 0.5
tomo$prior$calibration$lat_trend_east              <- 1
tomo$prior$calibration$lat_trend_west              <- 1
tomo$prior$calibration$hilevel                     <- 0
tomo$prior$calibration$power_perc                  <- 1
tomo$prior$calibration$mean_type                   <- "IS_EXP"

tomo$prior$calibration$l_vert 	                   <- 800 * 10^3
tomo$prior$calibration$l_hor               	       <- 20 #800 / 111
tomo$prior$calibration$l_long 	                   <- 20 #1000 / 40 
tomo$prior$calibration$a_scale                     <- 0.1641364

# tomo$prior$calibration$l_vert 	                   <- 500 * 10^3
# tomo$prior$calibration$l_hor               	       <- 40
# tomo$prior$calibration$l_long 	                   <- 40
# tomo$prior$calibration$a_scale                     <- 0.05523687


# tomo$prior$calibration$l_vert 	                   <- tomo$prior$covariance_parameters$l_vert 
# tomo$prior$calibration$l_hor               	       <- tomo$prior$covariance_parameters$l_hor  
# tomo$prior$calibration$l_long 	                   <- tomo$prior$covariance_parameters$l_long 
# tomo$prior$calibration$a_scale                     <- tomo$prior$covariance_parameters$a_scale 


# tomo$prior$calibration$l_vert                      <- 800 * 10^3
# tomo$prior$calibration$l_hor                       <- 10
# tomo$prior$calibration$l_long                      <- 10 
# tomo$prior$calibration$a_scale                     <- 0.1075899

# tomo$prior$calibration$l_vert 	                   <- 250 * 10^3
# tomo$prior$calibration$l_hor               	       <- 20 #800 / 111
# tomo$prior$calibration$l_long 	                   <- 20 #1000 / 40 
# tomo$prior$calibration$a_scale                     <- 0.1995142

# tomo$prior$calibration$l_vert                      <- 100 * 10^3
# tomo$prior$calibration$l_hor                       <- 10
# tomo$prior$calibration$l_long                      <- 10
# tomo$prior$calibration$a_scale                     <- 0.1712219
 

# Prior parameters for station bias variances
# "station_bias_var_unestimated" is used when a new bias_subject is introduced
# #station_bias_var_estimated" is used for biases once estimated
station_bias_mean_unestimated <- 0
station_bias_var_estimated    <- 2^2
station_bias_var_unestimated  <- 100^2 #300^2
tomo$prior$bias$station_default    <- data.frame(
 "station_bias_mean_unestimated" = c(10, rep(station_bias_mean_unestimated, 5)),
 "station_bias_var_estimated"    = c(10, station_bias_var_estimated, 
                                         station_bias_var_estimated, 
                                         station_bias_var_estimated, 
                                         station_bias_var_estimated, 
                                         station_bias_var_estimated), 
 "station_bias_var_unestimated"  = c(10, rep(station_bias_var_unestimated, 5)),                            
 "satellite_type" = c("LEO", "GPS", "GALILEO", "GLONASS", "BEIDOU", "SWARM")   
)

#tomo$prior$bias$station_default[tomo$prior$bias$station_default$satellite_type == "GLONASS", "station_bias_var_estimated"] <- 4

satellite_bias_mean_unestimated <- 0
satellite_bias_var_estimated    <- 0.000001^2
satellite_bias_var_unestimated  <- 0.000001^2
tomo$prior$bias$satellite_default   <- data.frame(
   "satellite_bias_mean_unestimated" = rep(satellite_bias_mean_unestimated, 6),
   "satellite_bias_var_estimated"    = rep(satellite_bias_var_estimated, 6), 
   "satellite_bias_var_unestimated"  = rep(satellite_bias_var_unestimated, 6),                            
   "satellite_type" = c("LEO", "GPS", "GALILEO", "GLONASS", "BEIDOU", "SWARM")
)

# Empty dataframes that will be updated with bias estimates
tomo$prior$bias$station <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tomo$prior$bias$station) <- c( "station_bias_mean", "station_bias_var", "satellite_type", "bias_subject")

tomo$prior$bias$satellite <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(tomo$prior$bias$satellite) <- c( "satellite_bias_mean", "satellite_bias_var", "satellite_type", "satellite")

# ----------------------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------------------

tomo$style$SAVE_PLOTS <- TRUE

# Plot stuff during run
# Doesn't affect "SAVE_PLOTS"
tomo$style$PLOT_STUFF <- FALSE

# image/plot parameters
tomo$style$nelim   <- c(0, 20)
tomo$style$teclim  <- c(0, 50)
tomo$style$altlim  <- c(0, 800)
# tomo$style$longlim <- c(18, 32)
# tomo$style$latlim  <- c(58, 72)
tomo$style$longlim <- c(5, 40)
tomo$style$latlim  <- c(55, 80)
# tomo$style$longlim <- c(125, 150)
# tomo$style$latlim  <- c(25, 55)

tomo$style$lat  <- 65
tomo$style$lat2 <- 69
tomo$style$long <- 25
tomo$style$long2 <- 19

tomo$style$lat_3d  <- NULL
tomo$style$long_3d <- 25

tomo$style$col_palette     <- viridis::viridis_pal()(50)
tomo$style$col_palette_tec <- viridis::viridis_pal()(50)

# figure sizes

tomo$style$png_height    <- 800
tomo$style$png_width     <- 1400
tomo$style$png_pointsize <- 25

tomo$style$png_one_col <- 800
tomo$style$png_one_row <- 800

tomo$style$png_col <- 600
tomo$style$png_row <- 500

tomo$style$png_height_3d <- 1000
tomo$style$png_width_3d  <- 1000

tomo$style$png_height_tec <- 1000
tomo$style$png_width_tec  <- 1000

tomo$style$png_height_pp <- tomo$style$png_height_tec
tomo$style$png_width_pp  <- tomo$style$png_width_tec 

tomo$style$interpolation_resolution <- 350

tomo$style$interpolation_resolution_TEC <- 1000

tomo$style$interp <- FALSE


tomo$style$TEC_lab <- parse(text=paste("TEC~~units~~(10^16/m^2)", sep=" "))
tomo$style$long_lab <- "Longitude (°)"
tomo$style$lat_lab <- "Latitude (°)"
tomo$style$el_lab <- "Elevation (°)"
tomo$style$ne_lab <- parse(text=paste("Ne~~(10^11/m^3)", sep=" "))
tomo$style$alt_lab <- "Alt (km)"

tomo$style$legend.args <-list(text=parse(text=paste("Electron~~density~~(10^11/m^3)", sep=" ")), col=1, cex = 1, side = 4, line = 2)
tomo$style$legend_tec.args <-list(text=tomo$style$TEC_lab, col=1, cex = 1, side = 4, line = 2)
#tomo$style$main <- paste("Reconstruction: ", format(domain$t1, format = "%Y-%m-%d %H:%M"), " UTC, ", mainlab, " ", coord, " ", u, sep = "")
