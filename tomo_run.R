#!/usr/bin/env Rscript
#library(ionotomor)
#export OMP_NUM_THREADS=1

#devtools::load_all("/Users/norberg/Dropbox/R_packages/ionotomor")
#devtools::load_all("/home/ubuntu/ionotomor")
#setwd("/Users/norberg/Dropbox/Projects/2023_01_21-26")
#setwd("/home/ubuntu/2023_01_21-26")

# ------------------------------------------------------------------------------
# Read configurations and GNSS data
# ------------------------------------------------------------------------------

library(ionotomor)
source("tomo_paths.R")
source("tomo_parameters.R")
source("tomo_parameters_remote.R")

cat("\nReading GNSS data\n")

files <- list.files(tomo$GNSS$directory)
paths <- list.files(tomo$GNSS$directory, full.names = TRUE)

file_time_char <- format(substr(files, start = 8, stop = 16), format = "%Y")
file_time <- strptime(file_time_char, "%y%j%H%M", tz = "UTC")

verbose <- FALSE

shiny_path <- "/srv/shiny-server/www"

get_latest <- function(path) {
  paths <- dir(path, full.names=TRUE)
  return(paths[which.max(file.info(paths)$ctime)])
}
online <- TRUE


# ------------------------------------------------------------------------------
# Bias  calibration run
# ------------------------------------------------------------------------------
cat("\nStarting bias calibration\n")

source("tomo_paths.R")
source("tomo_parameters.R")
source("tomo_parameters_remote.R")

t0 <- as.POSIXct("2023-01-20 10:00:00", tz = "UTC")
t1 <- as.POSIXct("2023-01-20 10:20:00", tz = "UTC")
t0u <- as.numeric(t0)
t1u <- as.numeric(t1)
tomo$domain$t0 <- t0
tomo$domain$t1 <- t1

file_is <- c(which((file_time > t0) & (file_time <= t1)), head(which(file_time > t1), 2))

GNSS_DATA_ALL <- NULL
for ( file_i in file_is) {
  h5_path <- paths[file_i]

  GNSS_DATA_ALL0 <- ionotomor::read_omotomo(h5_path,
	earth_rad = 6371 * 10^3,
	el_filt = c(0, 90),
	return_as_datatable = TRUE,
	verbose = FALSE,
	rm_n = 5)
  
  GNSS_DATA_ALL <- rbind(GNSS_DATA_ALL, GNSS_DATA_ALL0)
}

GNSS_DATA <- as.data.frame(GNSS_DATA_ALL[(t >= t0u) & (t <= t1u),])

tomo$GNSS$averaging_t <- (as.numeric(t1) - as.numeric(t0))

CALIBRATE_MEAN <- TRUE
ESTIMATE_DCB <- TRUE
tomo$GNSS$modelling_error <- 5

tomo$IONOSONDE$USE_AS_DIRECT_MEASUREMENTS <- TRUE
tomo$IONOSONDE$topside <- TRUE
tomo$IONOSONDE$alt_limit <- 800 * 10^3
tomo$IONOSONDE$times  <- 20 #minutes backwards from t1
tomo$IONOSONDE$times_ahead <- 0

tomo$utils$label <- "_1"
# undebug(solve_posterior)
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
source("tomo_parameters_remote.R")
# t1 <- as.POSIXct("2023-01-21 11:30:00", tz = "UTC")

for( tomo_i in 1 : 1000) {
  print(tomo_i)
  t0 <- t1
  t1 <- t1 + tomo$minstep * 60  
  t0u <- as.numeric(t0)
  t1u <- as.numeric(t1)
  tomo$domain$t0 <- t0
  tomo$domain$t1 <- t1

  if (t1u > max(GNSS_DATA_ALL$t)) {
	GNSS_DATA_ALL <- as.data.frame(GNSS_DATA_ALL[(t >= t0u) & (t <= t1u),])
   
	while (t1u > max(GNSS_DATA_ALL$t)) {
	  #file_i  <-  max(file_is)
	  file_i <- file_i + 1
	  h5_path <- paths[file_i]

	  GNSS_DATA_ALL0 <- ionotomor::read_omotomo(h5_path,
		earth_rad = 6371 * 10^3,
		el_filt = c(0, 90),
		return_as_datatable = TRUE,
		verbose = FALSE,
		rm_n = 5)
  
	  GNSS_DATA_ALL <- data.table::rbindlist(list(GNSS_DATA_ALL, GNSS_DATA_ALL0))
	}
  }

  GNSS_DATA <- as.data.frame(GNSS_DATA_ALL[(t >= t0u) & (t <= t1u),])

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
  res_paths <- dir(paste0(tomo$utils$results_directory, "/tomo"), full.names = TRUE)
  tomo_check_path <- res_paths[which.max(file.info(res_paths)$ctime)]
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
    tomo_check <- readRDS(file = res_paths[which.max(file.info(res_paths)$ctime)])
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
