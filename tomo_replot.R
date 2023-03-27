#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# ------------------------------------------------------------------------------
# devtools::load_all("/Users/norberg/Dropbox/R_packages/ionotomor")
library(ionotomor)
library(dplyr)
parameters_from_here <- TRUE

# setwd("/Users/norberg/Dropbox/Projects/2023_01_21-26")
# Real data tomo reconstruction files
results_main_dir <- "."
results_dir   <- args[1] #"results_remote_9"
# results_dir   <- "results_remote_9"

#par_def <- par()

par_ne <- function(cex = png_cex, mar = c(5, 5, 3, 2)) {
  par(
    mar      = mar,
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = cex,
    cex.lab  = cex,
    cex.main = cex
  )
}

par_tec <- function(cex = png_cex, mar = c(5, 5, 4, 3)) {
  par(
    mar      = mar,
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = cex,
    cex.lab  = cex,
    cex.main = cex
  )
}

par_rad <- function(cex = png_cex, mar = c(5, 5, 4, 2)) {
  par(
    mar      = mar,
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = cex,
    cex.lab  = cex,
    cex.main = 1
  )
}
par_prof <- function(cex = png_cex, mar = c(5, 5, 3, 2), mfrow = c(1, 1)) {
  par(
	mfrow    = mfrow,
	mar      = mar,
	xaxs     = "i",
	yaxs     = "i",
	cex.axis = cex,
	cex.lab  = cex,
	cex.main = cex
  )
}


nelim  <- c(0, 20)
teclim  <- c(0, 40)
altlim <- c(0, 800)
latlim <- c(54, 80)
longlim <- c(4, 35)
longitude <- 23
latitude <- 65
col_palette <- viridis::viridis_pal()(nelim[2] * 4 + 1)
col_palette_tec <- viridis::viridis_pal()(teclim[2] * 4 + 1)
#col_palette <- viridis::viridis_pal()(100)
nelim_diff <- c(-max(nelim), max(nelim))
ne_step_diff <- 1
col_palette_diff <- colorRampPalette(c("#FFE3E2", "#E8A8B4", "#DE93B3", "#D07FB5", "#BC6EB9", "#A261BB",  "#440154FF", "#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF", "#8FD744FF", "#FDE725FF"))(abs(diff(nelim_diff)) * (1 / ne_step_diff) + 1)

teclim_diff <- c(-max(teclim), max(teclim))
tec_step_diff <- 1
col_palette_diff_tec <- colorRampPalette(c("#FFE3E2", "#E8A8B4", "#DE93B3", "#D07FB5", "#BC6EB9", "#A261BB",  "#440154FF", "#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF", "#8FD744FF", "#FDE725FF"))(abs(diff(teclim_diff)) * (1 / tec_step_diff) + 1)

#col_palette_diff <- c(head(colorspace::heat_hcl(abs(nelim_diff[1]) + 1, h = c(0, -100), l = c(75, 40), c = c(40, 80), power = 1), -1), viridis::viridis_pal()(nelim_diff[2] + 1 ))
#colorspace::diverge_hcl(diff(nelim_diff) + 1) 
#RColorBrewer::brewer.pal(13, "Spectral")

png_width <- 800
png_height <- 400
png_units <- "px"
png_res <- 200
png_pointsize <- 6

#png_width_tec <- 2850
png_width_tec <- 500
png_height_tec <-  800

png_width_panel <-  600
png_height_panel <- 300

png_width_rad <- 360
png_height_rad <- 400

png_width_pp <- 800

lwd_rad <- 1.5
lwd_rad2 <- 2

png_cex <- 1

legend_cex <- 1

interp <- FALSE
direct_plot <- c("TR", "JR055")
altlab <- "Altitude (km)"
latlab_at_long <- paste("Latitude (°) | Longitude ", longitude, "°", sep = "")
latlab <- paste("Latitude (°)")
longlab <- "Longitude (°)"
longlab_at_lat <- paste("Longitude (°) | Latitude ", latitude, "°", sep = "")
deg <- "°"
#nelab <- parse(text=paste("Electron~~density~~(10^{11} / m^{3})", sep = " ")) #parse(text=paste("Ne~~(10^11/m^3)"))
nelab <- expression(paste("Ne (", 10^{11}, "/", m^{3}, ")"))
legendargs <- list(text = nelab, col = 1, cex = legend_cex, side = 4, line = 2)
#teclab <- parse(text=paste("TEC~~units~~(10^16/m^2)", sep=" "))
teclab <- expression(paste("TEC units (", 10^16, "/", m^2, ")"))
legendargs_tec <- list(text = teclab, col=1, cex = legend_cex, side = 4, line = 2)


nelab_panel <- expression(paste("Ne (", 10^{11}, "/", m^{3}, ")"))
legendargs_panel <- list(text = nelab_panel, col = 1, cex = legend_cex, side = 4, line = 2)



# Save destination
save_main_dir <- results_main_dir
#save_dir <- "plots_new"

save_dir <- "plots/panels"

#additional_label <- paste("_", strsplit(results_dir,  split = "_")[[1]][3], sep = "")
additional_label <- ""


# ------------------------------------------------------------------------------
getwd()
results_path <- file.path(results_main_dir, results_dir, "tomo")
# Save destination
save_path <- file.path(results_main_dir, results_dir, save_dir)
dir.create(save_path, showWarnings = TRUE, recursive = TRUE)
# setwd(save_path) # for system commands at the end



# ------------------------------------------------------------------------------
# Radar beams for TEC plots
beam_col <- 2
beam_col2 <- "orange"
beam_lwd <- 2

# UHF
earth_rad <- 6371 * 10^3
h        <- seq( 100000, 800000, 100000 )
E        <- 35 / 180 * pi
A        <- 145 / 180 * pi

uhf_lat <- 69+58/60
uhf_long <- 19+23/60
uhf_lat_r    <- uhf_lat / 180 * pi 
uhf_long_r   <- uhf_long / 180 * pi

pp <- function(E, A, lat_r, long_r, earth_rad = 6371 * 10^3, h = seq( 100000, 800000, 80000)) {
  PSI <- pi/2 - E - asin( ( earth_rad / ( earth_rad + h)) * cos(E) )
  pplat_rad <- asin( sin( lat_r ) * cos( PSI ) + cos( lat_r ) * sin( PSI ) * cos( A )) 
  pplong_rad <- (long_r + asin( ( sin( PSI ) * sin( A ) ) / cos( pplat_rad ) ) )
  pplat  <- pplat_rad / pi * 180
  pplong <- pplong_rad / pi * 180
  return(list("pplat" = pplat, "pplong" = pplong))
}




# ------------------------------------------------------------------------------
# Checking radar scan angles and selecting suitable ones

# E_vec_uhf  <- vector()
# A_vec_uhf  <- vector()
# E_vec_esr32  <- vector()
# A_vec_esr32  <- vector()
# 
# 
# for (i in 1 : length(times)) {
# 
#   tt <- as.POSIXct(times[i], tz = "UTC", format =  "%Y-%m-%d")
# 
#   print(tt)
#   t_stamp <- format(tt, "%j_%H%M")
# 
# 
#   paths <- dir(results_path, full.names = TRUE)
#   files <- dir(results_path, full.names = FALSE)
# 
#   tomo_times <- as.POSIXct(strptime(substr(files, start = 6, stop = 17),"%Y%m%d%H%M", tz = "UTC"))
# 
#   tomo <- readRDS(file = paths[which.min(abs(tomo_times - tt))])
# 
#   E_vec_uhf <- c(E_vec_uhf, (tomo$profiles_original$el[tomo$profiles_original$station == "uhf"]))
#   A_vec_uhf <- c(A_vec_uhf, (tomo$profiles_original$az[tomo$profiles_original$station == "uhf"]))
# 
#   E_vec_esr32 <- c(E_vec_esr32, (tomo$profiles_original$el[tomo$profiles_original$station == "32ma"]))
#   A_vec_esr32 <- c(A_vec_esr32, (tomo$profiles_original$az[tomo$profiles_original$station == "32ma"]))
# 
# }
# 
# E_uhf <- round(E_vec_uhf / 10) * 10
# A_uhf <- (round(A_vec_uhf / 10) * 10 ) %% 360
# A_uhf[E_uhf == 90] <- 180
# mlh_ind <- !(E_uhf < 90 & (A_uhf < 90 | A_uhf > 270))
# E_uhf <- E_uhf[mlh_ind]
# A_uhf <- A_uhf[mlh_ind]
# 
# mlh_pairs <- tibble(E_uhf, A_uhf) %>% group_by(E_uhf, A_uhf) %>% mutate(n = n()) %>% filter(n > 1000) %>% unique() %>%  mutate(lab  =  paste(E_uhf, A_uhf))

# unique(cbind(profiles_mlh$el,profiles_mlh$az))


# range(E_esr32)
# 
# E_esr32 <- round(E_vec_esr32 / 10) * 10
# A_esr32[E_esr32 > 90] <- (A_esr32[E_esr32 > 90] + 180) %% 360
# E_esr32[E_esr32 > 90] <- E_esr32[E_esr32 > 90] - 90
# A_esr32 <- (round(A_vec_esr32 / 10) * 10 ) %% 360
# A_esr32[E_esr32 == 90] <- 0
# esr32_ind <- !(E_esr32 < 90 & (A_esr32 < 110 | A_esr32 > 250))
# E_esr32 <- E_esr32[esr32_ind]
# A_esr32 <- A_esr32[esr32_ind]
# 
# unique(paste(E_esr32, A_esr32))

# ------------------------------------------------------------------------------
# Times for individual plots
# times <- as.POSIXct(c("2018-11-09 16:16:00", "2018-11-09 17:48:00"), tz = "UTC")

files <- list.files(file.path(results_path), ".RDS")
paths <- dir(file.path(results_path), ".RDS", full.names = TRUE)

times <- as.POSIXct(strptime(substr(files, start = 6, stop = 17),"%Y%m%d%H%M", tz = "UTC"))
# i_start <-  length(list.files("/Users/norberg/Dropbox/Projects/2021_ISSI_HSS//panels"))
# tt <- tail(list.files("/Users/norberg/Library/Mobile Documents/com~apple~CloudDocs/Projects/2021_ISSI_HSS/results_MIT_RO_15/panels"), 1)
# which.min(abs(strptime(tt,"%j_%H%M", tz = "UTC" ) - times))


# ------------------------------------------------------------------------------
# iiii <- 1
for (iiii in 1 : length(times)) {

  tt <- as.POSIXct(times[iiii], tz = "UTC", format =  "%Y-%m-%d")

  print(paste0(iiii, ": ", tt))
  t_stamp <- format(tt, "%j_%H%M")

  tt_lab <- format(tt, format =  "%Y-%m-%d %H:%M")

  # ------------------------------------------------------------------------------
  # Plot RECONSTRUCTED ionosphere

  paths <- dir(results_path, full.names = TRUE)
  files <- dir(results_path, full.names = FALSE)

  tomo_times <- as.POSIXct(strptime(substr(files, start = 6, stop = 17),"%Y%m%d%H%M", tz = "UTC"))

  tomo <- readRDS(file = paths[which.min(abs(tomo_times - tt))])

  profiles_mlh <- tomo$profiles_original[tomo$profiles_original$station %in% "mlh",]
  profiles_mlh$el <- round(profiles_mlh$el / 10) * 10
  profiles_mlh$az <-  round(profiles_mlh$az / 10) * 10
  profiles_mlh$az <- profiles_mlh$az %% 360
  profiles_mlh$az[profiles_mlh$el > 90] <- (profiles_mlh$az[profiles_mlh$el > 90] + 180) %% 360
  profiles_mlh$el[profiles_mlh$el > 90] <- 180 - profiles_mlh$el[profiles_mlh$el > 90]
  profiles_mlh$az[profiles_mlh$el == 90] <- 0
  profiles_mlh$az[profiles_mlh$az == 360] <- 0
  profiles_mlh <- profiles_mlh %>% mutate(lab = paste(round(el / 10) * 10, round(az / 10) * 10)) %>% filter(lab %in% paste(E_mlh, A_mlh))
  mlh_available <- unique(data.frame(el = round(profiles_mlh$el / 10) * 10, az = round(profiles_mlh$az / 10 ) * 10))


  latitudes <- tomo$domain$lat_ax
  altitudes  <- tomo$domain$alt_ax
  longitudes <- tomo$domain$long_ax

  ne_field <- ionotomor::get_ne_tomo(tomo, alt_m = NULL, lat_deg = NULL, long_deg = longitude)
  ne_field[ne_field < 0] <- 0

  png(paste(save_path, "/", t_stamp, "_reconst_lat.png", sep = ""), width = png_width, height = png_height, pointsize = png_pointsize, res = png_res, units = png_units)
  par_ne()
  fields::image.plot(z = t(ne_field)[, nrow(ne_field):1] * 10^5, y = tomo$domain$alt_ax / 1000, x = tomo$domain$lat_ax, 
                       zlim = nelim, xlim = latlim, ylim = altlim,  
                       ylab = altlab, xlab = latlab_at_long, 
                       col = col_palette, 
                       main = paste("TomoScand electron density |", tt_lab, "UTC"), 
                       axes = FALSE,
                       legend.args = legendargs)
  axis(1)
  axis(2)
#   points(x = 67.9, y = 200, pch = 4, col = 2, cex = 2)  
#   points(x = 74.89, y = 310, pch = 3, col = 2, cex = 2)    
#   legend("topright", pch = c(4, 3), legend = c("UHF beam location", "ESR32 beam location"), col = c(2, 2), bg = "white", pt.cex  = 2)

#   points(x  = 54.5, y = 750, pch = 22, col = 1, bg = "white", cex = 4)
#   text(labels = "A", x  = 54.5, y = 750, pch = 22, col = 1)  
  
  dev.off()

  ne_field <- ionotomor::get_ne_tomo(tomo, alt_m = NULL, lat_deg = latitude, long_deg = NULL)
  ne_field[ne_field < 0] <- 0


  png(paste(save_path, "/", t_stamp, "_reconst_long.png", sep = ""), width = png_width, height = png_height, pointsize = png_pointsize, res = png_res, units = png_units)
  par_ne()
  fields::image.plot(z = t(ne_field)[, nrow(ne_field):1] * 10^5, y = tomo$domain$alt_ax / 1000, x = tomo$domain$long_ax, 
                       zlim = nelim, xlim = longlim, ylim = altlim,  
                       ylab = altlab, xlab = longlab_at_lat, 
                       col = col_palette, 
                       main = paste("TomoScand electron density |", tt_lab, "UTC"), 
                       axes = FALSE,
                       legend.args = legendargs)
  axis(1)
  axis(2)
#   points(x = 67.9, y = 200, pch = 4, col = 2, cex = 2)  
#   points(x = 74.89, y = 310, pch = 3, col = 2, cex = 2)    
#   legend("topright", pch = c(4, 3), legend = c("UHF beam location", "ESR32 beam location"), col = c(2, 2), bg = "white", pt.cex  = 2)

#   points(x  = 54.5, y = 750, pch = 22, col = 1, bg = "white", cex = 4)
#   text(labels = "A", x  = 54.5, y = 750, pch = 22, col = 1)  
  
  dev.off()


  png(paste(save_path, "/", t_stamp, "_pp.png", sep = ""), width = png_width_pp, height = png_height_tec, pointsize = png_pointsize, res = png_res, units = png_units)
  par_tec(cex = png_cex)
  ionotomor::plot_tec(tomo, direct_measurements = tomo$profiles_original[tomo$profiles_original$station %in% direct_plot,],
           tec = FALSE, interp = interp, pp = TRUE,  pp_legend_style = 4, 
           teclim = teclim, rm.negat = TRUE,
           grid = FALSE, main_label_below = FALSE, 
           main = paste("Piercepoints 350 km \n", tt_lab, "UTC"),             
           map = TRUE, rectangle = FALSE,
           longlim_rect = longlim, latlim_rect = latlim,
           ylab = latlab, xlab = longlab, col = col_palette_tec)


  if (any(tomo$profiles_original$station == "mlh")) {
	text(labels = "mlh", x = mlh_long + 10, y = mlh_lat, col = 2, cex = png_cex)
  }

  for (iono_stat in unique(tomo$profiles_original$station[tomo$profiles_original$instrument == "ionosonde"])) {
	iono_i <- which(ionosondes$station == iono_stat)
	points(x = as.numeric(ionosondes$station_long[iono_i]), y = as.numeric(ionosondes$station_lat[iono_i]), col = 2, cex = png_cex)
	text(labels = ionosondes$station[iono_i], x = as.numeric(ionosondes$station_long[iono_i]) + 8, y = as.numeric(ionosondes$station_lat[iono_i]), col = 2, cex = png_cex)
  }

  RO <- tomo$satellite_measurements[!is.na(tomo$satellite_measurements$instrument) & tomo$satellite_measurements$instrument == "COSMIC",]

  points(x = RO$station_long, y = RO$station_lat, bg = "orange", pch = 21, cex = 2)
  points(x = RO$long, y = RO$lat, pch = "+", col = "orange", cex = 2)

  h <- seq( 100000, 550000, 100000 )
  E <- mlh_available$el / 180 * pi
  A <- mlh_available$az / 180 * pi

  for (ii in 1 : length(E)) {
    pp_coord <- pp(E[ii], A[ii], lat_r = mlh_lat_r, long_r = mlh_long_r, h = h)
    lines( x = c( mlh_long, pp_coord$pplong ), y = c( mlh_lat, pp_coord$pplat ), col = beam_col, lw = beam_lwd, cex = 1 )
    points( x = c(mlh_long, pp_coord$pplong ), y = c( mlh_lat, pp_coord$pplat ), col = beam_col, pch = 16, cex = png_cex )
  }
  
  graphics::rect(xleft = longlim[1], ybottom = latlim[1], xright = longlim[2], ytop = latlim[2], lwd = 2, lty = 2)
  graphics::rect(xleft = longlim[1] + 0.2 , ybottom = latlim[1] + 0.2, xright = longlim[2] - 0.2, ytop = latlim[2] - 0.2, lwd = 1, lty = 2, border = "white")

  lines(x = c(longitude, longitude), y = latlim, lty = 2, xpd = TRUE) 
  lines(y = c(latitude, latitude), x = longlim, lty = 2, xpd = TRUE) 

  dev.off()


  tec_rec <- ionotomor::get_tec_tomo(tomo)
  tec_rec[tec_rec < 0] <- 0
  
    png(paste(save_path, "/", t_stamp, "_tec.png", sep = ""), width = png_width_tec, height = png_height_tec, pointsize = png_pointsize, res = png_res, units = png_units)
  par_tec(cex = png_cex)
  fields::image.plot(
      x = longitudes,
      y = latitudes, 
      z = t(tec_rec)[, nrow(tec_rec) : 1],
      main = paste("TomoScand TEC \n", tt_lab, "UTC"),                             
      xlab = longlab,
      ylab = latlab,
      ylim = latlim,
      xlim = longlim,
      zlim = teclim,
      col = col_palette_tec,
      legend.args = legendargs_tec    
    )
  maps::map("world", xlim = longlim, ylim = latlim, add = TRUE)

  lines(x = c(longitude, longitude), y = latlim, lty = 2, xpd = TRUE) 
  lines(y = c(latitude, latitude), x = longlim, lty = 2, xpd = TRUE) 

  dev.off()



  # ------------------------------------------------------------------------------
  # Radar profiles
  # UHF

#   include_prof <- c("uhf")
#   profiles_mlh <- tomo$profiles_original[tomo$profiles_original$station %in% include_prof,]
#   profiles_mlh <- profiles_mlh %>% mutate(lab = paste(round(el / 10) * 10, round(az / 10) * 10))


# Plotting function with some global parameters
plot_isr <- function(RAD_name, E, A, lab = NULL, profiles, t_stamp = t_stamp, lat, long) {

  png(paste(save_path, "/", t_stamp, "_prof_", RAD_name, "_", E, "_", A, ".png", sep = ""), 
			width = png_width_rad, height = png_height_rad, 
			pointsize = png_pointsize, res = png_res, units = png_units)

  
  par_rad()
  plot(NULL, main = paste(lab, RAD_name, "|el:", E, "°, az:", A, "°", sep = ""), xlim = nelim, ylim = altlim, xlab = nelab, ylab = altlab, cex.main = png_cex, axes = FALSE)
  
  profs <- profiles %>% filter(lab %in% paste(E, A))
  profs <- split(profs, profs$lab)
  if (length(profs) > 0) {
    E <- profs[[1]]$el[1]
    A <- profs[[1]]$az[1]
    for (prof in profs) {
      lines(y = prof$alt / 10^3, x = prof$ne / 10^11, col = 2, lwd = lwd_rad2)
    }
  }
  axis(1)
  axis(2, labels = TRUE)
  box()

  location_at_alt1 <- pp_lat_long(lat_r = lat/ 180 * pi, 
									   long_r = long/ 180 * pi, 
									   E / 180 * pi, 
									   A / 180 * pi, 
									   h = alt1)

  raytrace_indices <- getTheoryMatrixCpp(
	start_lat = lat,
	start_long = long,
	start_alt = 100,
	stop_lat = location_at_alt1$lat,
	stop_long = location_at_alt1$long,
	stop_alt = alt1,
	lat_grid = tomo$domain$lat_grid,
	long_grid = tomo$domain$long_grid,
	alt_grid = tomo$domain$alt_grid,
	path_resolution = 100)
  
  # Plot corresponding profiles from prior and posterior mean
  alt_ax_vec <- rep(rep(rev(tomo$domain$alt_ax), 
							each = tomo$domain$lat_N), tomo$domain$long_N)
  lines(y = alt_ax_vec[raytrace_indices$j] / 1000, 
		x = tomo$prior$mean[raytrace_indices$j] * 10^5, 
		col = 5, lwd = lwd_rad)
  lines(y = alt_ax_vec[raytrace_indices$j] / 1000, 
		x = (tomo$prior$mean[raytrace_indices$j] - 1.96 * 
			tomo$prio$covariance_parameters$sd[raytrace_indices$j]) * 10^5, 
		col = 5, lwd = lwd_rad / 2, lty = 2)
  lines(y = alt_ax_vec[raytrace_indices$j] / 1000, 
		x = (tomo$prior$mean[raytrace_indices$j] + 1.96 * 
			tomo$prio$covariance_parameters$sd[raytrace_indices$j]) * 10^5, 
		col = 5, lwd = lwd_rad / 2, lty = 2)
  lines(y = alt_ax_vec[raytrace_indices$j] / 1000, 
		x = tomo$posterior$mean[raytrace_indices$j] * 10^5, 
		col = 1, lwd = lwd_rad)

  legend("topright", 
   legend = c("TomoScand", "Prior mean", "+/- 2*Prior SD", RAD_name), 
  lty = c(1, 1, 2, 1), lwd = c(lwd_rad, lwd_rad, lwd_rad / 2, lwd_rad), col = c(1, 5, 5, 2), cex = 1)

  dev.off()
}



plot_ionosonde <- function(ionosonde_name) {

  iono_i <- which(ionosondes$station == ionosonde_name)

  profiles_iono <- tomo$profiles_original[tomo$profiles_original$station %in% ionosonde_name,]
  iono_el <- 90
  iono_az <- 0
  iono_lat <- as.numeric(ionosondes$station_lat[ionosondes$station == ionosonde_name])
  iono_long <- as.numeric(ionosondes$station_long[ionosondes$station == ionosonde_name])
  iono_alt <- 100
  alt1 <- 1000 * 10^3

  png(paste(save_path, "/", t_stamp, "_prof_iono_", ionosonde_name, ".png", sep = ""), 
			width = png_width_rad, height = png_height_rad, 
			pointsize = png_pointsize, res = png_res, units = png_units)

  par_rad()
  plot(NULL, main = paste(ionosonde_name, "|el:", iono_el, "°, az:", iono_az, "°",  sep = ""), xlim = nelim, ylim = altlim, xlab = nelab, ylab = altlab, cex.main = png_cex, axes = FALSE)
  profs_iono <- split(profiles_iono, profiles_iono$t)
  for (prof in profs_iono) {
	lines(y = prof$alt / 10^3, x = prof$ne / 10^11, col = 2, lwd = lwd_rad2)
  }
  axis(1)
  axis(2, labels = TRUE)
  box()

  location_at_alt1 <- pp_lat_long(lat_r = iono_lat / 180 * pi, 
									   long_r = iono_long / 180 * pi, 
									   iono_el / 180 * pi, 
									   iono_az / 180 * pi, 
									   h = alt1)

  raytrace_indices_iono <- getTheoryMatrixCpp(
	start_lat = iono_lat,
	start_long = iono_long,
	start_alt = iono_alt,
	stop_lat = location_at_alt1$lat,
	stop_long = location_at_alt1$long,
	stop_alt = alt1,
	lat_grid = tomo$domain$lat_grid,
	long_grid = tomo$domain$long_grid,
	alt_grid = tomo$domain$alt_grid,
	path_resolution = 100)
  
  # Plot corresponding profiles from prior and posterior mean
  alt_ax_vec <- rep(rep(rev(tomo$domain$alt_ax), 
							each = tomo$domain$lat_N), tomo$domain$long_N)
  lines(y = alt_ax_vec[raytrace_indices_iono$j] / 1000, 
		x = tomo$prior$mean[raytrace_indices_iono$j] * 10^5, 
		col = 5, lwd = lwd_rad)
  lines(y = alt_ax_vec[raytrace_indices_iono$j] / 1000, 
		x = (tomo$prior$mean[raytrace_indices_iono$j] - 1.96 * 
			tomo$prio$covariance_parameters$sd[raytrace_indices_iono$j]) * 10^5, 
		col = 5, lwd = lwd_rad / 2, lty = 2)
  lines(y = alt_ax_vec[raytrace_indices_iono$j] / 1000, 
		x = (tomo$prior$mean[raytrace_indices_iono$j] + 1.96 * 
			tomo$prio$covariance_parameters$sd[raytrace_indices_iono$j]) * 10^5, 
		col = 5, lwd = lwd_rad / 2, lty = 2)
  lines(y = alt_ax_vec[raytrace_indices_iono$j] / 1000, 
		x = tomo$posterior$mean[raytrace_indices_iono$j] * 10^5, 
		col = 1, lwd = lwd_rad)


  legend("topright", 
   legend = c("TomoScand", "Prior mean", "+/- 2*Prior SD", ionosonde_name), 
  lty = c(1, 1, 2, 1), lwd = c(lwd_rad, lwd_rad, lwd_rad / 2, lwd_rad), col = c(1, 5, 5, 2), cex = 1)

  dev.off()
}


# for ( iii in 1 : length(E_mlh) ) {
#   plot_isr(RAD_name = "mlh", E = E_mlh[iii], A = A_mlh[iii], lab = NULL, profiles = profiles_mlh, t_stamp = t_stamp, lat = mlh_lat, long = mlh_long)
# }

for (ionosonde_name in direct_plot) {
 plot_ionosonde(ionosonde_name)
}



rad_ionos1 <- paste(save_path, "/", t_stamp , "_prof_iono_TR.png", sep = "")  
rad_ionos2 <- paste(save_path, "/", t_stamp , "_prof_iono_JR055.png", sep = "")  


piercepoints <- paste(save_path, "/", t_stamp , "_pp.png", sep = "")

tec_rec <- paste(save_path, "/", t_stamp, "_tec.png", sep = "")
ne_rec_lat <- paste(save_path, "/", t_stamp, "_reconst_lat.png", sep = "")
ne_rec_long <- paste(save_path, "/", t_stamp, "_reconst_long.png", sep = "")

append1 <- paste0(save_path, "/rad_append1.png")
ne_rec <- paste0(save_path, "/ne_rec.png")
final <- paste0(save_path, "/", t_stamp, ".png")
tmp <- paste0(save_path, "/tmp.png")

  system(paste("convert ", rad_ionos1," -crop 345x353+0+0  ",  rad_ionos1,   sep = ""))
  system(paste("convert ", rad_ionos2,  " -crop 345x390+0+0 ",  rad_ionos2,   sep = ""))

  
  system(paste("convert -append", rad_ionos1, rad_ionos2, "-append", append1, sep = " "))  


  system(paste("convert -append", ne_rec_lat, ne_rec_long, "-append", ne_rec, sep = " "))  
  system(paste("convert +append", ne_rec, tec_rec, "+append", ne_rec, sep = " "))  

  system(paste("convert +append", append1, piercepoints, "+append", tmp, sep = " "))  

  system(paste("convert -append", ne_rec, tmp, "-append ", final, sep = " "))  

  system(paste("rm",rad_mlh1, rad_mlh2, rad_mlh3, rad_mlh4, rad_ionos1, 
         rad_ionos2, rad_ionos3, rad_ionos4, rad_ionos5, rad_ionos6, rad_ionos7, rad_ionos8, rad_ionos9, tec_rec, piercepoints,
         ne_rec_lat, ne_rec_long, append1, ne_rec, tmp, sep = " "))

}