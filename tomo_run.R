#!/usr/bin/env Rscript
#library(ionotomor)
#export OMP_NUM_THREADS=1

#devtools::load_all("/Users/norberg/Dropbox/R_packages/ionotomor")
#devtools::load_all("/home/ubuntu/ionotomor")
#setwd("/Users/norberg/Dropbox/Projects/2011_03_tohoku")

# ------------------------------------------------------------------------------
# Read configurations and GNSS data
# ------------------------------------------------------------------------------

library(ionotomor)
source("tomo_paths.R")
source("tomo_parameters.R")

cat("\nReading GNSS data\n")

GNSS_data_time0 <- as.POSIXct("2011-03-09 00:00:00", tz = "UTC")
GNSS_data_time1 <- as.POSIXct("2011-03-10 00:00:00", tz = "UTC")

GNSS_DATA_ALL <- ionotomor::read_madrigal_gnss_los(h5_path = tomo$GNSS$directory,
  sat_modeling_alt = 1200 * 10^3,
  t_start = as.numeric(GNSS_data_time0),
  t_stop  = as.numeric(GNSS_data_time1),
  lat_filt = c(tomo$domain$lat_south, tomo$domain$lat_north), 
  long_filt = c(tomo$domain$long_west, tomo$domain$long_east),
  return_as_datatable = TRUE,
  verbose = TRUE,
  rm_n = 50)

GNSS_DATA_ALL <- GNSS_DATA_ALL[stec > 0, ]
# GNSS_DATA_ALL <- GNSS_DATA_ALL[stec > 0 & stec < 150 & stec_std < 10, ]

verbose <- FALSE

shiny_path <- "/srv/shiny-server/www"

get_latest <- function(path) {
  paths <- dir(path, full.names=TRUE)
  return(paths[which.max(file.info(paths)$ctime)])
}
online <- TRUE


# tmp <- GNSS_DATA_ALL[!duplicated(GNSS_DATA_ALL$station),]
# plot(y = tmp$station_lat, x = tmp$station_long, xlim = c(tomo$domain$long_west, tomo$domain$long_east))
# maps::map("world", add = TRUE)
# 
# tmp <- GNSS_DATA
# points(y = tmp$lat, x = tmp$long, xlim = c(tomo$domain$long_west, tomo$domain$long_east))
# maps::map("world", add = TRUE)


# ------------------------------------------------------------------------------
# Bias  calibration run
# ------------------------------------------------------------------------------
cat("\nStarting bias calibration\n")

# setwd("/Users/norberg/Dropbox/Projects/2011_03_tohoku")
source("tomo_paths.R")
source("tomo_parameters.R")
# Read parameters from the simulation case

t0 <- as.POSIXct("2011-03-09 00:00:00", tz = "UTC")
t1 <- as.POSIXct("2011-03-09 00:30:00", tz = "UTC")
t0u <- as.numeric(t0)
t1u <- as.numeric(t1)
tomo$domain$t0 <- t0
tomo$domain$t1 <- t1
GNSS_DATA <- as.data.frame(GNSS_DATA_ALL[(t >= t0u) & (t <= t1u),])
tomo$GNSS$averaging_t <- as.numeric(t1 - t0) * 60

CALIBRATE_MEAN <- TRUE
ESTIMATE_DCB <- TRUE

tomo$IONOSONDE$topside <- TRUE
tomo$IONOSONDE$alt_limit <- 800 * 10^3

tomo$utils$label <- "_1"

ionotomor::tomo_main(tomo, GNSS_DATA, CALIBRATE_MEAN, ESTIMATE_DCB, verbose)

if (online) {
  orig_imgs <- vector()
	orig_imgs[1] <- get_latest(paste0(tomo$utils$results_directory, "/plots/piercepoints"))
	orig_imgs[2] <- get_latest(paste0(tomo$utils$results_directory, "/plots/tec_map"))
	orig_imgs[3] <- get_latest(paste0(tomo$utils$results_directory, "/plots/alt_lat"))
	orig_imgs[4] <- get_latest(paste0(tomo$utils$results_directory, "/plots/alt_long"))
 	orig_imgs[5] <- get_latest(paste0(tomo$utils$results_directory, "/plots/RAD_profiles"))

  shiny_imgs <- vector()
	shiny_imgs[1] <- paste0(shiny_path, "/pp.png")
	shiny_imgs[2] <- paste0(shiny_path, "/tec.png")
	shiny_imgs[3] <- paste0(shiny_path, "/lat.png")
	shiny_imgs[4] <- paste0(shiny_path, "/long.png")
 	shiny_imgs[5] <- paste0(shiny_path, "/RAD_prof.png")

  file.copy(orig_imgs, shiny_imgs, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
cat("\nStarting actual tomo run\n")
source("tomo_paths.R")
source("tomo_parameters.R")
# t1 <- as.POSIXct("2011-03-09 23:55:00", tz = "UTC")

for( tomo_i in 1 : 1000) {
  print(tomo_i)
  t0 <- t1
  t1 <- t1 + tomo$minstep * 60  
  t0u <- as.numeric(t0)
  t1u <- as.numeric(t1)
  tomo$domain$t0 <- t0
  tomo$domain$t1 <- t1


  if (t1 >= GNSS_data_time1) {
	GNSS_data_time0 <- GNSS_data_time0 + 24 * 60 * 60
	GNSS_data_time1 <- GNSS_data_time1 + 24 * 60 * 60

    tomo$GNSS$directory <- paste0(tomo$utils$data_directory, "/GNSS_data/los_", format(GNSS_data_time0, format = "%Y%m%d"), ".001.h5")

    GNSS_DATA0 <- as.data.frame(GNSS_DATA_ALL[(GNSS_DATA_ALL$t >= t0u) & (GNSS_DATA_ALL$t <= t1u),])
	GNSS_DATA_ALL <- ionotomor::read_madrigal_gnss_los(h5_path = tomo$GNSS$directory,
	  sat_modeling_alt = 1200 * 10^3,
	  t_start = as.numeric(GNSS_data_time0),
	  t_stop  = as.numeric(GNSS_data_time1),
	  lat_filt = c(tomo$domain$lat_south, tomo$domain$lat_north), 
	  long_filt = c(tomo$domain$long_west, tomo$domain$long_east),
	  return_as_datatable = TRUE,
	  verbose = TRUE,
	  rm_n = 50)

    GNSS_DATA_ALL <- rbind(GNSS_DATA0, GNSS_DATA_ALL)
    GNSS_DATA_ALL <- GNSS_DATA_ALL[GNSS_DATA_ALL$stec > 0, ]

  }

  GNSS_DATA <- as.data.frame(GNSS_DATA_ALL[(GNSS_DATA_ALL$t >= t0u) & (GNSS_DATA_ALL$t <= t1u),])

  CALIBRATE_MEAN <- FALSE
  ESTIMATE_DCB <- FALSE
  verbose <- FALSE
  

# # For some reason this helps and MUMPS doesn't get killed
# perc_rmv <- 0.001
# n <- dim(GNSS_DATA)[1]
# n_rmv <- max(10, round(perc_rmv * n))
# rnd_vec <- sample(n, size = n_rmv, replace = FALSE)
# GNSS_DATA_rmv <- GNSS_DATA[-rnd_vec,]
# 
# ionotomor::tomo_main(tomo, GNSS_DATA_rmv, CALIBRATE_MEAN, ESTIMATE_DCB, verbose)  

  ionotomor::tomo_main(tomo, GNSS_DATA, CALIBRATE_MEAN, ESTIMATE_DCB, verbose)

  check_n <- 1
  cat("\nChecking quality #", check_n, "\n")
  paths <- dir(paste0(tomo$utils$results_directory, "/tomo"), full.names = TRUE)
  tomo_check_path <- paths[which.max(file.info(paths)$ctime)]
  tomo_check <- readRDS(file = tomo_check_path)
  # Try again without some random measurements
  perc_rmv <- 0.001
  n <- dim(GNSS_DATA)[1]
  while ( max(abs(tomo_check$posterior$mean))  > 5 * 10^(-3) & perc_rmv < 1 ) { 
    cat("\nBAD SOLUTION. STARTING RE-ANALYSING.\n")
    Sys.sleep(5)
    file.remove(tomo_check_path)
    
    n_rmv <- max(10, round(perc_rmv * n))
    rnd_vec <- sample(n, size = n_rmv, replace = FALSE)
    GNSS_DATA_rmv <- GNSS_DATA[-rnd_vec,]
    
    ionotomor::tomo_main(tomo, GNSS_DATA_rmv, CALIBRATE_MEAN, ESTIMATE_DCB, verbose)  
    tomo_check <- readRDS(file = paths[which.max(file.info(paths)$ctime)])
    check_n <- 1 + check_n
    # If solution is still not found, increase the number of removed measurements
    perc_rmv <- perc_rmv * 5
  }
  if (max(abs(tomo_check$posterior$mean))  > 5 * 10^(-3) & perc_rmv < 1 ){
    cat("\nCOULD NOT GET A PROPER SOLUTION. CHECK PREVIOUS SOLUTIONS\n")
    break
  } else {
    cat("\nStarting new time timestep\n\n")
    Sys.sleep(2)
  }
  
  if (online) {
	# Copy images to Shiny server
	orig_imgs <- vector()
	  orig_imgs[1] <- get_latest(paste0(tomo$utils$results_directory, "/plots/piercepoints"))
	  orig_imgs[2] <- get_latest(paste0(tomo$utils$results_directory, "/plots/tec_map"))
	  orig_imgs[3] <- get_latest(paste0(tomo$utils$results_directory, "/plots/alt_lat"))
	  orig_imgs[4] <- get_latest(paste0(tomo$utils$results_directory, "/plots/alt_long"))
 	  orig_imgs[5] <- get_latest(paste0(tomo$utils$results_directory, "/plots/RAD_profiles"))

	shiny_imgs <- vector()
	  shiny_imgs[1] <- paste0(shiny_path, "/pp.png")
	  shiny_imgs[2] <- paste0(shiny_path, "/tec.png")
	  shiny_imgs[3] <- paste0(shiny_path, "/lat.png")
	  shiny_imgs[4] <- paste0(shiny_path, "/long.png")
 	  shiny_imgs[5] <- paste0(shiny_path, "/RAD_prof.png")

	file.copy(orig_imgs, shiny_imgs, overwrite = TRUE)
  }

  cat("\nStarting new time timestep\n")
  Sys.sleep(2)
}


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
