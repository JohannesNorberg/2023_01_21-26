file_paths <- list.files("/Users/norberg/Dropbox/Projects/2023_01_21-26/data/ionosonde_data/TR", full.names = TRUE)
file_names <- list.files("/Users/norberg/Dropbox/Projects/2023_01_21-26/data/ionosonde_data/TR")
times <- strptime(as.vector(unlist(sapply(file_names, strsplit, ".txt"))), "%Y%m%d%H%M%S", tz = "UTC")

# times <- as.POSIXct(TEC[, 1], origin="1970-01-01", tz = "UTC")


iono_mat <- matrix(ncol = 3, nrow = length(times))

# iono_i <- 1
for (iono_i in 1 : length(times)) {
  print(iono_i)
  iono_file <- read.table(file_paths[iono_i])
  ne <- 1.24 * 10^10 * iono_file[, 2]^2
  alt <- iono_file[, 1]
  peak_ind <- which.max(ne)

  iono_mat[iono_i, 1] <- as.numeric(times[iono_i])
  iono_mat[iono_i, 2] <- alt[peak_ind]
  iono_mat[iono_i, 3] <- ne[peak_ind]
}



peak_ne <- iono_mat[, 3] / 10^11
times <- as.POSIXct(iono_mat[, 1], origin="1970-01-01", tz = "UTC")

times <- times[1 : 4601]
peak_ne <- peak_ne[1 : 4601]

png(paste0("ionosonde_ne_ts.png"), height = 500, width = 1500)
plot(y = peak_ne, x = times, type = 'l', axes = FALSE, 
     ylim = c(0, 20), xlim = c(min(times), max(times)),
     main = "Tromsø ionosonde peak Ne", ylab = "Ne (10^11/m^3)", xlab = "", lwd = 0.5)
box()
axis(side = 2, at = seq(0, 20, 5))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "day"), round(tail(times, 1), "day") , by = "day"), format = "%b %d")

filt_w <- 15
filt <- rep(1/filt_w, filt_w, sdes = 2)
peak_ne_filt_15 <- stats::filter(peak_ne, filter = filt)
lines(y = peak_ne_filt_15, x = times, col = 2, lwd = 2)

filt_w <- 30
filt <- rep(1/filt_w, filt_w)
peak_ne_filt_30 <- stats::filter(peak_ne, filter = filt)
lines(y = peak_ne_filt_30, x = times, col = 3, lwd = 2)

filt_w <- 60
filt <- rep(1/filt_w, filt_w)
peak_ne_filt_60 <- stats::filter(peak_ne, filter = filt)
lines(y = peak_ne_filt_60, x = times, col = 4, lwd = 2)
legend("topleft", legend = c("2 min data", "30 min moving average", "60 min moving average", "120 min moving average"), col = 1:4, lty = 1)
dev.off()


ne_diff <- peak_ne - peak_ne_filt_30
# ne_diff <- peak_ne_filt_15 - peak_ne_filt_30
png(paste0("ionosonde_ne_ts_diff.png"), height = 500, width = 1500)
plot(y = ne_diff, x = times, type = 'l', axes = FALSE, 
     ylim = c(-5, 5), xlim = c(min(times), max(times)),
     main = "Tromsø ionosonde peak Ne", ylab = "Ne (10^11/m^3)", xlab = "")
box()
axis(side = 2, at = seq(-5, 5, 1))
axis.POSIXct(side = 1, at = seq(trunc(times[1], "day"), round(tail(times, 1), "day") , by = "day"), format = "%b %d")

sd_vec <- vector()
for (i_sd in 0 : (length(times) - filt_w)) {
  sd_vec <- c(sd_vec, sd((ne_diff)[i_sd : (i_sd + filt_w)]))
}
lines(y = sd_vec, x = times[1:length(sd_vec)], col = 2, lwd = 2)
legend("topleft", legend = c("2 min data - 60 min moving average", "Standard deviation"), col = 1:2, lty = 1)

dev.off()



