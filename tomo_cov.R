library(Matrix)
library(rmumps)
library(ionotomor)

tomo <- readRDS(file = "/Users/norberg/Dropbox/Projects/2023_01_21-26/results_4/tomo/tomo_202301201025.RDS")
# tomo <- solve_posterior_sd(tomo)

longlim <- c(5, 40)
latlim  <- c(55, 80)

lat  <- 65
lat2 <- 69
long <- 25
long2 <- 19

sdlim <- 0.5# tomo$prior$calibration$power_perc * 10


col_palette     <- viridis::viridis_pal()(200)


names(tomo$prior$covariance_parameters)

L <- tomo$prior$L_mat
prec_mat <- Matrix::t(L) %*% L

# Sample
set.seed(1)
xi <- rnorm(dim(L)[1])

Lxi <- Matrix::t(L)%*%xi

system.time(smple <- rmumps::mumps.solve(mat = prec_mat, rhs = Lxi, np = 1, sym = 0, host.involved = 1, pivot.order = 7, out.of.core = 1))

smple_array_3D <- aperm(array(smple, dim = c(tomo$domain$lat_N, tomo$domain$alt_N, tomo$domain$long_N)), perm = c(2, 1, 3))

mean_array_3D <- aperm(array(tomo$prior$mean, dim = c(tomo$domain$lat_N, tomo$domain$alt_N, tomo$domain$long_N)), perm = c(2, 1, 3)) 

par(mfrow  = c(2, 2), mar = c(5, 4, 4, 4))
#lat-alt
plot_ne(tomo, array_3D = smple_array_3D, interp = FALSE, alt_m = NULL, lat_deg = NULL, long_deg = long, xlab = NULL, ylab = NULL, legend.args = NULL, main = NULL, add = FALSE, zlim = c(-sdlim, sdlim))
plot_ne(tomo, array_3D = smple_array_3D + mean_array_3D, interp = FALSE, alt_m = NULL, lat_deg = NULL, long_deg = long, xlab = NULL, ylab = NULL, legend.args = NULL, main = NULL, add = FALSE, zlim = c(0, max(smple_array_3D + mean_array_3D) * 10^5))
range(smple_array_3D + mean_array_3D)
#long-alt
plot_ne(tomo, array_3D = smple_array_3D, interp = FALSE, alt_m = NULL, lat_deg = lat, long_deg = NULL, xlab = NULL, ylab = NULL, legend.args = NULL, main = NULL, add = FALSE, zlim = c(-sdlim, sdlim))
plot_ne(tomo, array_3D = smple_array_3D + mean_array_3D, interp = FALSE, alt_m = NULL, lat_deg = lat, long_deg = NULL, xlab = NULL, ylab = NULL, legend.args = NULL, main = NULL, add = FALSE, zlim = c(0, max(smple_array_3D + mean_array_3D) * 10^5))

range(mean_array_3D)
range(smple_array_3D)
# plot_ne(tomo, array_3D = mean_array_3D, interp = FALSE, alt_m = NULL, lat_deg = lat, long_deg = NULL, xlab = NULL, ylab = NULL, legend.args = NULL, main = NULL, add = FALSE, zlim = c(0, max(smple_array_3D + mean_array_3D) * 10^5))

# Scale prior covariance
system.time(prior_var <- rmumps::mumps.calc.inverse.diagonal(mat = prec_mat, np = 1, sym = 2))

prior_sd <- sqrt(prior_var[ 1 : tomo$domain$N ])
prior_sd_array_3D <- aperm(array(prior_sd, dim = c(tomo$domain$lat_N, tomo$domain$alt_N, tomo$domain$long_N)), perm = c(2, 1, 3))
# Only prior$covariance_parameters$alpha in get_prior_cov.R is scaled with "a_scale"
# Hence none of the below profiles is affected by "a_scale". 
# tomo$prior$covariance_parameters$sd and thus initial_sd_array_3D are scaled correctly
initial_sd_array_3D <- aperm(array(tomo$prior$covariance_parameters$sd, dim = c(tomo$domain$lat_N, tomo$domain$alt_N, tomo$domain$long_N)), perm = c(2, 1, 3)) 
a_scale <- (max(prior_sd) / 
            (tomo$prior$peak_ne * tomo$prior$covariance_parameters$power_perc) * 
            sqrt(tomo$prior$covariance_parameters$a_scale))^2


E_F <- c(tomo$prior$covariance_parameters$E_layer$peak_ne, tomo$prior$peak_ne)
E_F_power <- c(tomo$prior$covariance_parameters$E_layer$power_perc, 
               tomo$prior$covariance_parameters$power_perc)
E_F_max <- which.max(E_F)

a_scale <- (max(prior_sd) / 
            (E_F[E_F_max] * E_F_power[E_F_max]) * 
            sqrt(tomo$prior$covariance_parameters$a_scale))^2


tomo$prior$covariance_parameters$a_scale


# lat2 <- 69
# long <- 19#25
par(mfrow  = c(1, 1), mar = c(5, 4, 4, 4))
# Initial prior profile
plot(y = rev(tomo$domain$alt_ax) / 10^3 , x = initial_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))] * 10^5, col = 5, xlim = c(0, max(initial_sd_array_3D) * 10^5 + 0.2))
lines(y = rev(tomo$domain$alt_ax) / 10^3 , x = initial_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))] * 10^5, type = 'l', col = 5)

# Solved prior profile
#plot(y = rev(tomo$domain$alt_ax) / 10^3 , x = prior_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))]*10^5, type = 'l', col = 2)
lines(y = rev(tomo$domain$alt_ax) / 10^3 , x = prior_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))] * 10^5, type = 'l', col = 3)
points(y = rev(tomo$domain$alt_ax) / 10^3 , x = prior_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))] * 10^5, col = 3)
abline(h = 100)
tomo$prior$covariance_parameters$E_layer$peak_alt / 10 ^3

initial_F_sd_array_3D <- aperm(array(tomo$prior$covariance_parameters$sd_F, dim = c(tomo$domain$lat_N, tomo$domain$alt_N, tomo$domain$long_N)), perm = c(2, 1, 3)) 
initial_E_sd_array_3D <- aperm(array(tomo$prior$covariance_parameters$E_layer$sd, dim = c(tomo$domain$lat_N, tomo$domain$alt_N, tomo$domain$long_N)), perm = c(2, 1, 3)) 
lines(y = rev(tomo$domain$alt_ax) / 10^3 , x = initial_F_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))] * 10^5 * tomo$prior$covariance_parameters$power_perc, type = 'l', col = 3)
lines(y = rev(tomo$domain$alt_ax) / 10^3 , x = initial_E_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))] * 10^5 * tomo$prior$covariance_parameters$E_layer$power_perc , type = 'l', col = 2)

lines(y = rev(tomo$domain$alt_ax) / 10^3 , 
      x = (initial_F_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), 
                                   which.min(abs(tomo$domain$long_ax - long))] * 
          tomo$prior$covariance_parameters$power_perc + 
          initial_E_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), 
                                  which.min(abs(tomo$domain$long_ax - long))] * 
          tomo$prior$covariance_parameters$E_layer$power_perc) * 10^5, 
      type = 'l', col = 4)




points(y = rev(tomo$domain$alt_ax) / 10^3, x = tomo$prior$covariance_parameters$sd_prof * 10^5 * tomo$prior$covariance_parameters$a_scale * , col =  2)
lines(y = rev(tomo$domain$alt_ax) / 10^3, x = tomo$prior$covariance_parameters$sd_prof * 10^5 * tomo$prior$covariance_parameters$a_scale, col =  2)



# Prior covariance profile

alt_ind <- which.min(abs(rev(tomo$domain$alt_ax) - tomo$prior$peak_alt)) - 1
#floor(tomo$domain$alt_N / 2)
middlepixel <- floor(tomo$domain$long_N / 2) * tomo$domain$lat_N * tomo$domain$alt_N + tomo$domain$lat_N * alt_ind + ceiling(tomo$domain$lat_N / 2)

lat_inds <- middlepixel + (1 : tomo$domain$lat_N - ceiling(tomo$domain$lat_N / 2))
lat_mask <- spMatrix(nrow = tomo$domain$N, ncol = tomo$domain$N, i = rep(middlepixel, length(lat_inds)), j = lat_inds, x = rep(1, length(lat_inds)))

alt_inds <- middlepixel + (1 : tomo$domain$alt_N - alt_ind) * tomo$domain$lat_N
#c(middlepixel - ((floor(tomo$domain$alt_N / 2) : 1) * tomo$domain$lat_N), middlepixel + ((0 : (ceiling(tomo$domain$alt_N / 2) - 1)) * tomo$domain$lat_N))
alt_mask <- spMatrix(nrow = tomo$domain$N, ncol = tomo$domain$N, i = rep(middlepixel, length(alt_inds)), j = alt_inds, x = rep(1, length(alt_inds)))

# masks can be combined and solved at once
mask <- alt_mask + lat_mask
mask[middlepixel, middlepixel] <- 1

system.time(prior_cov_prof <- rmumps::mumps.calc.inverse.elements(mat = prec_mat, mask = mask, np = 4, sym = 2))

par(mfrow  = c(3, 2), mar = c(5, 4, 4, 4))
#correlation
plot(y = prior_cov_prof[middlepixel, lat_inds] / (prior_sd[lat_inds] * prior_sd[middlepixel]), x = tomo$domain$lat_ax - tomo$domain$lat_ax[ceiling(tomo$domain$lat_N / 2)])
abline(h = 0.1, v = c(-tomo$prior$covariance_parameters$l_hor , tomo$prior$covariance_parameters$l_hor ) )

plot(y = prior_cov_prof[middlepixel, alt_inds] / (prior_sd[alt_inds] * prior_sd[middlepixel]), x = rev(tomo$domain$alt_ax - rev(tomo$domain$alt_ax)[alt_ind])/10^3)
abline(h = 0.1, v = c(-tomo$prior$covariance_parameters$l_vert, tomo$prior$covariance_parameters$l_vert ) / 10^3 )

lines(y = prior_cov_prof[middlepixel, alt_inds] / (prior_sd[alt_inds] * prior_sd[middlepixel]), x = rev(tomo$domain$alt_ax - rev(tomo$domain$alt_ax)[alt_ind])/10^3, col = 1)
points(y = prior_cov_prof[middlepixel, alt_inds] / (prior_sd[alt_inds] * prior_sd[middlepixel]), x = rev(tomo$domain$alt_ax - rev(tomo$domain$alt_ax)[alt_ind])/10^3, col = 1)
abline(h = 0.1 )




# Posterior covariance
tomo$IONOSONDE$USE_AS_DIRECT_MEASUREMENTS <- TRUE

tomo <- solve_posterior_sd(tomo)

tomo <- tomo1
tomo <- tomo2

post_sd <- sqrt(tomo$posterior$sd)
post_sd_array_3D <- aperm(array(post_sd, dim = c(tomo$domain$lat_N, tomo$domain$alt_N, tomo$domain$long_N)), perm = c(2, 1, 3))

lines(y = rev(tomo$domain$alt_ax) / 10^3 , x = post_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))] * 10^5, type = 'l', col = 6)
points(y = rev(tomo$domain$alt_ax) / 10^3 , x = post_sd_array_3D[, which.min(abs(tomo$domain$lat_ax - lat2)), which.min(abs(tomo$domain$long_ax - long))] * 10^5, col = 6)
abline(h = 100)

 

long3 <- long2
par(mfrow  = c(2, 1), mar = c(5, 4, 4, 4))

plot_ne(tomo, array_3D = prior_sd_array_3D, long_deg = long3)

plot_ne(tomo, array_3D = post_sd_array_3D, long_deg = long3)


plot_ne(tomo, array_3D = (prior_sd_array_3D - post_sd_array_3D) / prior_sd_array_3D * 10^(-5), long_deg = long3)

abline(h = 495)