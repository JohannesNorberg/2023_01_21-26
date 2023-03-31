#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# devtools::load_all("/Users/norberg/Dropbox/R_packages/ionotomor")
library(ionotomor)
library(dplyr)
parameters_from_here <- TRUE

setwd("/Users/norberg/Dropbox/Projects/2023_01_21-26")
# Real data tomo reconstruction files
results_main_dir <- "."
# results_dir   <- args[1] #"results_remote_9"
results_dir   <- "results_remote_9"


col_palette <- viridis::viridis_pal()(10 * 4 + 1)
col_palette_diff <- viridis::viridis_pal()(10 * 2 + 1)

legend_cex <- 1 

interp <- FALSE
direct_plot <- c("TR", "JR055")
altlab <- "Altitude (km)"
latlab <- paste("Latitude (°)")
longlab <- "Longitude (°)"
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

# ------------------------------------------------------------------------------

files <- list.files(file.path(results_path), ".RDS")
paths <- dir(file.path(results_path), ".RDS", full.names = TRUE)

times <- as.POSIXct(strptime(substr(files, start = 6, stop = 17),"%Y%m%d%H%M", tz = "UTC"))


# ------------------------------------------------------------------------------

tomo <- readRDS(file = paths[1])

refs <- matrix(ncol = 2, byrow = TRUE,
               c(69 + 34/60, 26 + 42/60,
                 67 + 22/60, 26 + 38/60,
                 62.24,      25.75,
                 60 + 12/60, 24 + 57/60))

refs_index <- matrix(ncol = 2, nrow = nrow(refs))

for (loc_i in 1 : nrow(refs)) {
  refs_index[loc_i, 1] <- which.min(abs(rev(tomo$domain$lat_ax)  - refs[loc_i, 1]))
  refs_index[loc_i, 2] <- which.min(abs(tomo$domain$long_ax - refs[loc_i, 2]))
}

TEC <- matrix( ncol = nrow(refs), nrow = 0)

# iiii <- 13
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

  tec_rec <- ionotomor::get_tec_tomo(tomo)
  tec_rec[tec_rec < 0] <- 0
  
  TEC_ROW <- vector()
  for (refs_i_row in 1 : nrow(refs_index)) {
    TEC_ROW[refs_i_row] <- tec_rec[refs_index[refs_i_row, 1], refs_index[refs_i_row, 2]]  
  }
  TEC <- rbind(TEC, TEC_ROW)
  
}

png(paste0(results_dir, "/VTEC_ts.png"), height = 500, width = 1500)
plot(y = TEC[, 1], x = times, type = 'l', axes = FALSE, 
     ylim = c(0, 35), xlim = c(min(times), max(times) + 60 * 60 * 8),
     main = "VTEC", ylab = "TECU (10^16/m^2)", xlab = "")
lines(y = TEC[, 2], x = times, col = 2)
lines(y = TEC[, 3], x = times, col = 3)
lines(y = TEC[, 4], x = times, col = 4)
box()
axis(side = 2, at = seq(0, 35, 10))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "day"), round(tail(times, 1), "day") , by = "day"), format = "%b %d")
legend("topright", legend = rev(c("Utsjoki", "Sodankylä", "Jyväskylä", "Helsinki")), lty = 1, col = 4:1, cex = 1.5)
dev.off()





ref_long <- 25
ref_long_ind <- which.min(abs(tomo$domain$long_ax - ref_long))

TEC_lat <- matrix( ncol = length(times), nrow = tomo$domain$lat_N)


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

  tec_rec <- ionotomor::get_tec_tomo(tomo)
  tec_rec[tec_rec < 0] <- 0
  
  TEC_lat[, iiii] <- tec_rec[,ref_long_ind]
  
}

#write.table(TEC_lat, paste0("/Users/norberg/Dropbox/Projects/2023_01_21-26/data/TEC_data/vtec_lat_time_20-26.txt"))

png(paste0(results_dir, "/VTEC_lat_time.png"), height = 500, width = 1500)
fields::image.plot(main = "VTEC",
                   xlab = "Time", ylab = "Latitude (°)",z = t(TEC_lat)[,nrow(TEC_lat):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(0,40), col = col_palette)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "day"), round(tail(times, 1), "day") , by = "day"), format = "%b %d")
abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()

TEC_lat_mean <- matrix(ncol = length(unique(times_ind)), nrow = nrow(TEC_lat))
times_ind <- dense_rank(format(times, format = "%H%M"))

for (iiii in 1 : length(unique(times_ind))) {
  TEC_lat_mean[,iiii] <- rowMeans(TEC_lat[, which(times_ind == iiii)])
}

TEC_lat_scaled <-  TEC_lat - TEC_lat_mean[, times_ind]

png(paste0(results_dir, "/VTEC_lat_time-avg.png"), height = 500, width = 1500)
fields::image.plot(main = "VTEC - average day",
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_scaled)[,nrow(TEC_lat_scaled):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(-10, 10), col = col_palette_diff)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "day"), round(tail(times, 1), "day") , by = "day"), format = "%b %d")
#abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()



png(paste0(results_dir, "/VTEC_lat_time-avg_21.png"), height = 500, width = 1500)
fields::image.plot(main = "Jan 21 | VTEC - average day",
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_scaled)[,nrow(TEC_lat_scaled):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(-10, 10), 
                   xlim = as.POSIXct(c("2023-01-21 00:00:00 UTC", "2023-01-22 00:00:00 UTC"), tz = "UTC"),
                   col = col_palette_diff)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "hour"), round(tail(times, 1), "hour") , by = "hour"), format = "%H:%M")
abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()

png(paste0(results_dir, "/VTEC_lat_time-avg_22.png"), height = 500, width = 1500)
fields::image.plot(main = "Jan 22 | VTEC - average day",
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_scaled)[,nrow(TEC_lat_scaled):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(-10, 10), 
                   xlim = as.POSIXct(c("2023-01-22 00:00:00 UTC", "2023-01-23 00:00:00 UTC"), tz = "UTC"),
                   col = col_palette_diff)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "hour"), round(tail(times, 1), "hour") , by = "hour"), format = "%H:%M")
abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()

png(paste0(results_dir, "/VTEC_lat_time-avg_23.png"), height = 500, width = 1500)
fields::image.plot(main = "Jan 23 | VTEC - average day",
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_scaled)[,nrow(TEC_lat_scaled):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(-10, 10), 
                   xlim = as.POSIXct(c("2023-01-23 00:00:00 UTC", "2023-01-23 23:59:59 UTC"), tz = "UTC"),
                   col = col_palette_diff)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "hour"), round(tail(times, 1), "hour") , by = "hour"), format = "%H:%M")
abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()

png(paste0(results_dir, "/VTEC_lat_time-avg_24.png"), height = 500, width = 1500)
fields::image.plot(main = "Jan 24 | VTEC - average day",
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_scaled)[,nrow(TEC_lat_scaled):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(-10, 10), 
                   xlim = as.POSIXct(c("2023-01-24 00:00:00 UTC", "2023-01-24 23:59:59 UTC"), tz = "UTC"),
                   col = col_palette_diff)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "hour"), round(tail(times, 1), "hour") , by = "hour"), format = "%H:%M")
abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()

png(paste0(results_dir, "/VTEC_lat_time-avg_25.png"), height = 500, width = 1500)
fields::image.plot(main = "Jan 25 | VTEC - average day",
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_scaled)[,nrow(TEC_lat_scaled):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(-10, 10), 
                   xlim = as.POSIXct(c("2023-01-25 00:00:00 UTC", "2023-01-25 23:59:59 UTC"), tz = "UTC"),
                   col = col_palette_diff)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "hour"), round(tail(times, 1), "hour") , by = "hour"), format = "%H:%M")
abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()

png(paste0(results_dir, "/VTEC_lat_time-avg_26.png"), height = 500, width = 1500)
fields::image.plot(main = "Jan 26 | VTEC - average day",
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_scaled)[,nrow(TEC_lat_scaled):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(-10, 10), 
                   xlim = as.POSIXct(c("2023-01-26 00:00:00 UTC", "2023-01-26 23:59:59 UTC"), tz = "UTC"),
                   col = col_palette_diff)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "hour"), round(tail(times, 1), "hour") , by = "hour"), format = "%H:%M")
abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()








library(fields)

filter_width <- 180 #min
# filtfilt <- rep(1, filter_width /  5)/(filter_width / 5)
filtfilt <- dnorm(1 : (filter_width / 5) - (filter_width / 5)/2, sd = 5)

TEC_smooth <- t(apply(TEC_lat, FUN = function(x){stats::filter(x, filter = filtfilt)}, MARGIN = 1))
TEC_lat_smoothed <-  TEC_lat - TEC_smooth

png(paste0(results_dir, "/VTEC_lat_time-smooth_norm.png"), height = 500, width = 1500)
#par(mfrow = c(2, 1))
fields::image.plot(main = "VTEC - smoothed VTEC",
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_smoothed)[,nrow(TEC_lat_smoothed):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), 
                   zlim = c(-3, 3), 
                   col = col_palette_diff)
box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "day"), round(tail(times, 1), "day") , by = "day"), format = "%b %d")
#abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)
dev.off()

png(paste0(results_dir, "/VTEC_lat_time-smooth_65ts.png"), height = 500, width = 1500)
lat_ind <- which.min(abs(rev(tomo$domain$lat_ax) - 65))
plot(y = TEC_lat[lat_ind,], x = times, type = 'l', main = "VTEC and smoothed VTEC along lat 65°", ylab = "VTEC", xlab = "Time")
lines(y = TEC_smooth[lat_ind,], x= times, col = 2)
dev.off()

DAY <- 20
DAY <- DAY + 1

png(paste0(results_dir, "/VTEC_lat_time-smooth", "day_panels1",".png"), height = 800, width = 1500)
par(mfrow = c(3, 2))
for (DAY in 21:23) {

fields::image.plot(main = paste0("Jan ", DAY," | VTEC - smoothed VTEC"),
                   xlab = "Time", ylab = "Latitude (°)",
                   z = t(TEC_lat_smoothed)[,nrow(TEC_lat_smoothed):1], 
                   x = times, y = tomo$domain$lat_ax, 
                   axes = FALSE, ylim = c(55, 75), zlim = c(-3, 3), 
                   xlim = as.POSIXct(c(paste0("2023-01-", DAY, " 00:00:00 UTC"), paste0("2023-01-", DAY, " 23:59:50 UTC")), tz = "UTC"),
                   col = col_palette_diff)



box()
axis(side = 2, at = seq(55, 75, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "hour"), round(tail(times, 1), "hour") , by = "hour"), format = "%H:%M")
abline(v  = times[format(times, format = "%M") %in% c("00")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c("0000")], lwd = 0.5)
# abline(v  = times[format(times, format = "%H%M") %in% c( "1200")], lty = 2, lwd = 0.5)

lat_ind <- which.min(abs(rev(tomo$domain$lat_ax) - 65))
plot(y = TEC_lat[lat_ind,], x = times, xlim = as.POSIXct(c(paste0("2023-01-", DAY, " 00:00:00 UTC"), paste0("2023-01-", DAY,  " 23:59:59 UTC")), tz = "UTC"), type = 'l', main = "VTEC and smoothed VTEC along lat 65°", ylab = "VTEC", xlab = "Time", axes = FALSE)
lines(y = TEC_smooth[lat_ind,], x = times, col = 2)
axis.POSIXct(side = 1, at = seq(trunc(times[1], "hour"), round(tail(times, 1), "hour") , by = "hour"), format = "%H:%M")
axis(side = 2, at = seq(0, 30, 5))

}
dev.off()
