# Read biases from latest tomo file from the run 
library("ggplot2")
library("ggmap")
library("viridis")  
library("data.table")
library("dplyr")

setwd("/Users/norberg/Dropbox/Projects/2023_01_21-26/")
file_paths <- list.files("results_10/tomo", full.names = TRUE)
i <- 248
tomo <- readRDS(file_paths[i])
 
 
files <- list.files(tomo$GNSS$directory)
paths <- list.files(tomo$GNSS$directory, full.names = TRUE)

file_time_char <- format(substr(files, start = 8, stop = 16), format = "%Y")
file_time <- strptime(file_time_char, "%y%j%H%M", tz = "UTC")

verbose <- FALSE

# ------------------------------------------------------------------------------
#
# ------------------------------------------------------------------------------

t0 <- as.POSIXct("2023-01-20 10:00:00", tz = "UTC")
t1 <- as.POSIXct("2023-01-20 10:20:00", tz = "UTC")
t0u <- as.numeric(t0)
t1u <- as.numeric(t1)

earth_rad <- 6371 * 10^3

file_is <- c(which((file_time > t0) & (file_time <= t1)), head(which(file_time > t1), 2))

GNSS_DATA_ALL <- NULL
for ( file_i in file_is) {
  h5_path <- paths[file_i]

  GNSS_DATA_ALL0 <- ionotomor::read_omotomo(h5_path,
	earth_rad = earth_rad,
	el_filt = c(0, 90),
	return_as_datatable = TRUE,
	verbose = FALSE,
	rm_n = 5)
  
  GNSS_DATA_ALL <- rbind(GNSS_DATA_ALL, GNSS_DATA_ALL0[!(satellite_type == "GLONASS")])
}

GNSS_DATA <- GNSS_DATA_ALL %>% 
  dplyr::left_join(x = tomo$posterior$bias$station[,c( "station_bias_mean", "station_bias_var", "bias_subject")], 
                   by = c("bias_subject")) %>%
  dplyr::mutate(stec_bias_rm = stec - station_bias_mean)

h_pp <- 350 * 10^3

khi <- asin(earth_rad  / (earth_rad  + 225 * 10^3) * cos(GNSS_DATA$el / 180 * pi))
GNSS_DATA$vtec <- cos(khi) * GNSS_DATA$stec_bias_rm

E <- GNSS_DATA$el / 180 * pi
A <- GNSS_DATA$az / 180 * pi
lat_r <- GNSS_DATA$station_lat / 180 * pi
long_r <- GNSS_DATA$station_long / 180 * pi

PSI_pp <- pi / 2 - E - asin((earth_rad / (earth_rad + h_pp)) * cos(E))

GNSS_DATA$pp_lat <- asin(sin(lat_r) * cos(PSI_pp) + cos(lat_r) * sin(PSI_pp) * cos(A)) / pi * 180
GNSS_DATA$pp_alt <- h_pp
GNSS_DATA$pp_long <- (long_r + asin((sin(PSI_pp) * sin(A)) / cos(GNSS_DATA$pp_lat / 180 * pi))) / pi * 180


longlim <- c(0, 40)
latlim <- c(55, 75)

maks <- 50
li <- c(0, maks)
la <- c(seq(0, maks, 10), paste0(">", maks))
br <- c(seq(0,  maks, 10), 100)
minstep <- 4
verbose <- FALSE
world <- map_data("world")

  if( format(t0, "%d") == format(t1, "%d")) {
    tlab <- paste(format(t0, "%Y-%m-%d %H:%M"), "--", format(t1, "%H:%M"))
  } else {
    tlab <- paste(format(t0, "%Y-%m-%d %H:%M"), "--", format(t1, "%m-%d %H:%M"))
  }

    
  label <- paste0("GNSS VTEC ", h_pp / 1000, " (km) | ", tlab)


  print(
	ggplot() + 
	  geom_point(
		data = GNSS_DATA,
		aes(x = pp_long, y = pp_lat, color = vtec),
		alpha = 0.7
	  ) +
# 	  geom_point(aes(x = epiclong, y = epiclat), shape = 21, size = 5, stroke = 2) +	  
	  geom_map(
		data = world, map = world,
		aes(map_id = region),
		color = "black", fill = "lightgray", linewidth = 0.5,
		alpha = 0.0
	  ) +  
# 	  scale_x_continuous(limits = c(-145, -55)) +
# 	  scale_y_continuous(limits = c(20, 55)) +
 	  scale_x_continuous(limits = c(0, 40)) +
 	  scale_y_continuous(limits = c(55, 75)) +
	  theme_minimal() +
	  scale_color_viridis(name = "VTEC", limits = li, breaks = br, labels = la, na.value = "red") +
	  ggtitle(label) +
	  xlab("Longitude (°)") +
	  ylab("Latitude (°)")
  )







t1 <- as.POSIXct("2023-01-20 10:20:00", tz = "UTC")
minstep <- 5
tomo_i <- 1
for( tomo_i in 1 : 5000) {
  print(tomo_i)
  t0 <- t1
  t1 <- t1 + minstep * 60  
  t0u <- as.numeric(t0)
  t1u <- as.numeric(t1)

  if( format(t0, "%d") == format(t1, "%d")) {
    tlab <- paste(format(t0, "%Y-%m-%d %H:%M"), "--", format(t1, "%H:%M"))
  } else {
    tlab <- paste(format(t0, "%Y-%m-%d %H:%M"), "--", format(t1, "%m-%d %H:%M"))
  }
    
  label <- paste0("GNSS VTEC ", h_pp / 1000, " (km) | ", tlab)

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


  GNSS_DATA <- GNSS_DATA_ALL %>% 
				 filter(satellite_type != "GLONASS" & (t >= t0u) & (t <= t1u)) %>%
				   dplyr::left_join(x = tomo$posterior$bias$station[,c( "station_bias_mean", "station_bias_var", "bias_subject")], 
					 by = c("bias_subject")) %>%
				  dplyr::mutate(stec_bias_rm = stec - station_bias_mean)

  h_pp <- 350 * 10^3

  khi <- asin(earth_rad  / (earth_rad  + 225 * 10^3) * cos(GNSS_DATA$el / 180 * pi))
  GNSS_DATA$vtec <- cos(khi) * GNSS_DATA$stec_bias_rm

  E <- GNSS_DATA$el / 180 * pi
  A <- GNSS_DATA$az / 180 * pi
  lat_r <- GNSS_DATA$station_lat / 180 * pi
  long_r <- GNSS_DATA$station_long / 180 * pi

  PSI_pp <- pi / 2 - E - asin((earth_rad / (earth_rad + h_pp)) * cos(E))

  GNSS_DATA$pp_lat <- asin(sin(lat_r) * cos(PSI_pp) + cos(lat_r) * sin(PSI_pp) * cos(A)) / pi * 180
  GNSS_DATA$pp_alt <- h_pp
  GNSS_DATA$pp_long <- (long_r + asin((sin(PSI_pp) * sin(A)) / cos(GNSS_DATA$pp_lat / 180 * pi))) / pi * 180

  png(paste0("vtec_plots/vtec_", format(t1, "%Y%m%d%H%M"), ".png"), width = 600, height = 600, pointsize = 25)
  print(
	ggplot() + 
	  geom_point(
		data = GNSS_DATA,
		aes(x = pp_long, y = pp_lat, color = vtec),
		alpha = 0.7
	  ) +
	  geom_map(
		data = world, map = world,
		aes(map_id = region),
		color = "black", fill = "lightgray", linewidth = 0.5,
		alpha = 0.0
	  ) +  
 	  scale_x_continuous(limits = longlim) +
 	  scale_y_continuous(limits = latlim) +
	  theme_minimal() +
	  scale_color_viridis(name = "VTEC", limits = li, breaks = br, labels = la, na.value = "red") +
	  ggtitle(label) +
	  xlab("Longitude (°)") +
	  ylab("Latitude (°)")
  )
  dev.off()
}


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


