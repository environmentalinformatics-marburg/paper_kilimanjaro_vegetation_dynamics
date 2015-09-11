### environmental stuff

## working directory
library(Orcs)
setwdOS()

## packages
lib <- c("reshape2", "raster", "remote", "lubridate", "doParallel", 
         "RColorBrewer", "Rsenal", "ggplot2", "grid", "rgdal", "rasterVis")
jnk <- sapply(lib, function(i) library(i, character.only = TRUE))

## functions
ch_dir_repo <- "repositories/magic/dfg_for_kilimanjaro/ndvi_kilimanjaro/src/gimms/"
source("repositories/magic/dfg_for_kilimanjaro/fire_ndvi/src/panel.smoothconts.R")
source(paste0(ch_dir_repo, "importOni.R"))
source(paste0(ch_dir_repo, "importDmi.R"))
source(paste0(ch_dir_repo, "extentEnsoSeason.R"))

## parallelization
cl <- makeCluster(4)
registerDoParallel(cl)

## input filepath
ch_dir_fig <- "publications/paper/detsch_et_al__ndvi_dynamics/figures/"
ch_dir_data <- paste0(ch_dir_fig, "data/")


### data processing

## start and end
st <- "1982-01"
nd <- "2011-12"

## oni, dmi
oni_mlt <- importOni(file = "/media/fdetsch/XChange/kilimanjaro/gimms3g/gimms3g/data/oni/enso_and_iod.csv")
dat_dmi <- importDmi(file = "/media/fdetsch/XChange/kilimanjaro/gimms3g/gimms3g/data/dmi/dmi.dat")
dat_oni_dmi <- merge(oni_mlt, dat_dmi, all = TRUE, by = 1)

# dat_dmi <- read.table("data/dmi/dmi.dat")
# names(dat_dmi) <- c("Year", "Month", "DMI")
# 
# dat_dmi[, 2] <- formatC(dat_dmi[, 2], width = 2, flag = "0")
# dates <- paste(dat_dmi[, 1], dat_dmi[, 2], "01", sep = "-")
# dat_dmi$Date <- as.Date(dates)
# 
# id_st <- grep(st, dat_dmi$Date)
# dat_dmi <- dat_dmi[id_st:nrow(dat_dmi), c("Date", "DMI")]
# 
# # merge oni, dmi
# dat_oni_dmi <- merge(oni_mlt, dat_dmi, all = TRUE, by = 1)

## ndvi
fls_ndvi <- paste0(ch_dir_data, "gimms_ndvi3g_dwnscl_8211.tif")
rst_ndvi <- stack(fls_ndvi)
mat_ndvi <- as.matrix(rst_ndvi)

ndvi_date <- seq(as.Date("1982-01-01"), as.Date("2011-12-01"), "month")

# exclusion of areas with low ndvi, cf brown et al. (2006), based on upper quartile
fls_ndvi_quan75 <- paste0(ch_dir_data, "ndvi_quan75.tif")
rst_ndvi_quan75 <- raster(fls_ndvi_quan75)
val_ndvi_quan75 <- getValues(rst_ndvi_quan75)

id_lndvi <- which(val_ndvi_quan75 <= .3)
rst_ndvi[id_lndvi] <- NA

# areas with low variance, based on 10-90% iqr
fls_ndvi_quan10 <- paste0(ch_dir_data, "ndvi_quan10.tif")
rst_ndvi_quan10 <- raster(fls_ndvi_quan10)
val_ndvi_quan10 <- getValues(rst_ndvi_quan10)

fls_ndvi_quan90 <- paste0(ch_dir_data, "ndvi_quan90.tif")
rst_ndvi_quan90 <- raster(fls_ndvi_quan90)
val_ndvi_quan90 <- getValues(rst_ndvi_quan90)

id_lvar <- which((val_ndvi_quan90 - val_ndvi_quan10) <= .1)

# tmp <- rst_ndvi[[1]]
# tmp[id_lvar] <- NA
# plot(tmp)

# deseason
rst_ndvi_dsn <- deseason(rst_ndvi)
mat_ndvi_dsn <- as.matrix(rst_ndvi_dsn)

# long-term means, cf. anyamba et al. (2002)
fls_ndvi_ltm <- paste0(ch_dir_data, "gimms_ndvi3g_dwnscl_8211_ltm.tif")
rst_ndvi_ltm <- stack(fls_ndvi_ltm)
mat_ndvi_ltm <- as.matrix(rst_ndvi_ltm)
med_ndvi_ltm <- apply(mat_ndvi_ltm, 2, function(...) median(..., na.rm = TRUE))
df_med_ndvi_ltm <- data.frame(x = 1:12, y = med_ndvi_ltm[c(7:12, 1:6)])

rst_ndvi_anom <- (rst_ndvi/rst_ndvi_ltm - 1) * 100

# exclusion of areas with low ndvi, cf brown et al. (2006)
rst_ndvi[id_lndvi] <- NA
rst_ndvi_dsn[id_lndvi] <- NA
rst_ndvi_anom[id_lndvi] <- NA

mat_ndvi_dsn <- as.matrix(rst_ndvi_dsn)
mat_ndvi_anom <- as.matrix(rst_ndvi_anom)

# remove all np pixels
np <- readOGR(dsn = paste0(ch_dir_data, "shp/"), layer = "fdetsch-kilimanjaro-1420532792846", 
              p4s = "+init=epsg:4326")
np <- spTransform(np, CRS = "+init=epsg:21037")

rst_ndvi_anom_rmnp <- rst_ndvi_anom
rst_ndvi_anom_rmnp <- mask(rst_ndvi_anom_rmnp, np, inverse = TRUE)
mat_ndvi_anom_rmnp <- as.matrix(rst_ndvi_anom_rmnp)


# ### correlation
# 
# ## oni ~ ndvi
# 
# # el nino
# groups <- list("all El Ninos", "pure El Ninos", "pure m/s El Ninos", 
#                "El Ninos w IOD+", "m/s El Ninos w IOD+", "purest IOD+")
# 
# groups_nino <- list(c("WE", "ME", "SE"),
#                     c("WE", "ME", "SE"),
#                     c("ME", "SE"), 
#                     c("WE", "ME", "SE"), 
#                     c("ME", "SE"), 
#                     "WE")
# 
# groups_iod <- list(c("M", "", "P"), 
#                    "", 
#                    "", 
#                    "P", 
#                    "P", 
#                    "P")
# 
# ndvi_sp <- foreach(nino = groups_nino, iod = groups_iod, group = groups, 
#                    span = rep(8, length(groups)), .packages = lib, 
#                    .combine = "rbind") %dopar% {
#   oni_mlt_sub <- subset(oni_mlt, (Type %in% nino) & (IOD %in% iod))
# 
#   # remove incomplete data (e.g. late 2011)
#   seasons <- split(oni_mlt_sub, as.factor(oni_mlt_sub$Season))
#   seasons_len <- sapply(seasons, nrow)
#   seasons <- seasons[seasons_len == 12]
#   
#   oni_mlt_sub <- lapply(seasons, function(i) extentEnsoSeason(i, ndvi = rst_ndvi, span = span))
#   oni_mlt_sub <- do.call("rbind", oni_mlt_sub)
#   
#   sp <- smooth.spline(oni_mlt_sub$Month, oni_mlt_sub$NDVI, spar = .01)
#   sp_pred <- predict(sp, seq(1, nlevels(oni_mlt_sub$Month), .01))
#   
#   mat_sp_pred <- do.call("cbind", sp_pred)
#   df_sp_pred <- data.frame(group = group, mat_sp_pred)
# #   df_sp_pred$x <- unique(oni_mlt_sub$Month)
#   
#   return(df_sp_pred)
# }
# 
# red <- brewer.pal(4, "Reds")
# blue <- brewer.pal(4, "Blues")
# 
# plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen")
# names(plot.colors) <- groups
# 
# # extend ltm
# # df_med_ndvi_ltm$x <- factor(df_med_ndvi_ltm$x)
# ltm_ext <- merge(data.frame(x = 1:19), df_med_ndvi_ltm, by = "x", 
#                  all = TRUE)
# ltm_ext[is.na(ltm_ext$y), 2] <- ltm_ext[1:sum(is.na(ltm_ext$y)), 2]
# # ltm_ext$x <- as.numeric(as.character(ltm_ext[, 1]))
# 
# # x-axis labels
# lbl <- rep(c(7:12, 1:6, 7:12, 1:2))
# names(lbl) <- 1:(12+span-1)
# 
# # df_med_ndvi_ltm$x <- factor(df_med_ndvi_ltm$x)
# ltm_ext <- merge(data.frame(x = 1:20), df_med_ndvi_ltm, by = "x", 
#                  all = TRUE)
# ltm_ext[is.na(ltm_ext$y), 2] <- ltm_ext[1:sum(is.na(ltm_ext$y)), 2]
# # ltm_ext$x <- as.numeric(as.character(ltm_ext[, 1]))
# 
# plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen")
# names(plot.colors) <- groups
# 
# p_nino <- ggplot(aes(x, y, group = group, colour = group), data = ndvi_sp) + 
#   geom_line() +
#   geom_line(aes(x, y), data = ltm_ext, colour = "grey65", linetype = 2) + 
#   scale_colour_manual("", values = plot.colors) + 
# #   scale_x_discrete("\nMonth", labels = lbl) + 
#   scale_x_continuous("\nMonth", breaks = 1:(12+span-1), labels = lbl) + 
#   labs(x = "\nMonth", y = expression(atop(NDVI[median], "\n"))) + 
#   theme_bw()
# 
# png("vis/cor_ndvi_oni/ts_nino_ndvi.png", width = 30, height = 12, units = "cm", 
#     pointsize = 18, res = 300)
# print(p_nino)
# dev.off()
# 
# # occurrence of each group
# ssn_cnt_nino <- foreach(nino = groups_nino, iod = groups_iod, group = groups, 
#                         span = rep(8, length(groups)), .combine = "rbind") %do% {
#   oni_mlt_sub <- subset(oni_mlt, (Type %in% nino) & (IOD %in% iod))
#   ssn_cnt <- length(unique(oni_mlt_sub$Season))
#   df_ssn_cnt <- data.frame(Group = group, Count = ssn_cnt)
#   return(df_ssn_cnt)
# }
# 
# # la nina
# groups <- list("all La Ninas", "pure La Ninas", "pure m/s La Ninas", 
#                "La Ninas w IOD+", "m/s La Ninas w IOD+", "purest IOD+")
# 
# groups_nina <- list(c("WL", "ML", "SL"),
#                     c("WL", "ML", "SL"),
#                     c("ML", "SL"), 
#                     c("WL", "ML", "SL"), 
#                     c("ML", "SL"), 
#                     "WL")
# 
# groups_iod <- list(c("M", "", "P"), 
#                    "", 
#                    "", 
#                    "P", 
#                    "P", 
#                    "P")
# 
# nina <- groups_nina[[6]]
# iod <- groups_iod[[6]]
# group <- groups[[6]]
# span <- 8
# 
# ndvi_sp <- foreach(nina = groups_nina, iod = groups_iod, group = groups, 
#                    span = rep(8, length(groups)), .packages = lib, 
#                    .combine = "rbind") %dopar% {
#   oni_mlt_sub <- subset(oni_mlt, (Type %in% nina) & (IOD %in% iod))
#   
#   # remove incomplete data (e.g. late 2011)
#   seasons <- split(oni_mlt_sub, as.factor(oni_mlt_sub$Season))
#   seasons_len <- sapply(seasons, nrow)
#   seasons <- seasons[seasons_len == 12]
#   
# #   if (span > 0) {
# #     dates_wspan <- lapply(seasons, function(i, span = span) {
# #       dates <- i$Date
# #       
# #       dates_st <- dates[1]
# #       dates_nd <- dates[length(dates)]
# #       month(dates_nd) <- month(dates_nd) + span
# #       
# #       dates <- seq(dates_st, dates_nd, "month")
# #       
# #       return(dates)
# #     })
# #     dates_wspan <- do.call("c", dates_wspan)
# #     dates_wspan <- dates_wspan[-which(duplicated(dates_wspan))]
# #     
# #     oni_mlt_sub <- oni_mlt[oni_mlt$Date %in% dates_wspan, ]
# #   }
# #   
# #   dates_sub <- oni_mlt_sub$Date
# #   ndvi_id <- which(ndvi_date %in% dates_sub)
# #   rst_ndvi_sub <- rst_ndvi[[ndvi_id]]
# #   mat_ndvi_sub <- as.matrix(rst_ndvi_sub)
# #   
# #   num_ndvi_sub_median <- apply(mat_ndvi_sub, 2, 
# #                                function(...) median(..., na.rm = TRUE))
# #   
# #   oni_mlt_sub$NDVI <- num_ndvi_sub_median
# #   oni_mlt_sub$Month <- as.numeric(substr(oni_mlt_sub$Date, 6, 7))
# 
# oni_mlt_sub <- lapply(seasons, function(i) extentEnsoSeason(i, span = span))
# oni_mlt_sub <- do.call("rbind", oni_mlt_sub)
# 
#   sp <- smooth.spline(oni_mlt_sub$Month, oni_mlt_sub$NDVI, spar = .01)
#   sp_pred <- predict(sp, seq(1, nlevels(oni_mlt_sub$Month), .01))
#   
#   mat_sp_pred <- do.call("cbind", sp_pred)
#   df_sp_pred <- data.frame(group = group, mat_sp_pred)
# #   df_sp_pred$x <- unique(oni_mlt_sub$Month)
#   return(df_sp_pred)
# }
# 
# red <- brewer.pal(4, "Reds")
# blue <- brewer.pal(4, "Blues")
# 
# plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen")
# names(plot.colors) <- groups
# 
# # extend ltm
# # df_med_ndvi_ltm$x <- factor(df_med_ndvi_ltm$x)
# ltm_ext <- merge(data.frame(x = 1:20), df_med_ndvi_ltm, by = "x", 
#                  all = TRUE)
# ltm_ext[is.na(ltm_ext$y), 2] <- ltm_ext[1:sum(is.na(ltm_ext$y)), 2]
# # ltm_ext$x <- as.numeric(as.character(ltm_ext[, 1]))
# 
# # x-axis labels
# lbl <- rep(c(7:12, 1:6, 7:12, 1:2))
# names(lbl) <- 1:(12+span-1)
# 
# p_nina <- ggplot(aes(x, y, group = group, colour = group), data = ndvi_sp) + 
#   geom_line() +
#   geom_line(aes(x, y), data = ltm_ext, colour = "grey65", 
#             linetype = 2) + 
#   scale_colour_manual("", values = plot.colors) + 
# #   scale_x_discrete("\nMonth", labels = lbl) + 
#   scale_x_continuous("\nMonth", breaks = 1:(12+span-1), labels = lbl) + 
#   labs(x = "\nMonth", y = expression(atop(NDVI[median], "\n"))) + 
#   theme_bw()
# 
# png("vis/cor_ndvi_oni/ts_nina_ndvi.png", width = 30, height = 12, units = "cm", 
#     pointsize = 18, res = 300)
# print(p_nina)
# dev.off()
# 
# # occurrence of each group
# ssn_cnt_nina <- foreach(nina = groups_nina, iod = groups_iod, group = groups, 
#                         span = rep(8, length(groups)), .combine = "rbind") %do% {
#   oni_mlt_sub <- subset(oni_mlt, (Type %in% nina) & (IOD %in% iod))
#   ssn_cnt <- length(unique(oni_mlt_sub$Season))
#   df_ssn_cnt <- data.frame(Group = group, Count = ssn_cnt)
#   return(df_ssn_cnt)
# }
#                      

## ccf

# remove low-variace areas (evergreen forest)
rst_ndvi_anom_rmf <- rst_ndvi_anom
rst_ndvi_anom_rmf[id_lvar] <- NA
mat_ndvi_anom_rmf <- as.matrix(rst_ndvi_anom_rmf)

# lag with maximum correlation
Find_Max_CCF <- function(a, b, lag.max = NULL) {
  d <- ccf(a, b, plot = FALSE, lag.max = lag.max)
  cor = d$acf[,,1]
  lag = d$lag[,,1]
  res = data.frame(cor,lag)
  res_max = res[which.max(res$cor),]
  return(res_max)
} 

# ccf
ls_ccf <- foreach(id = c("ONI", "DMI")) %do% {
  dat_ccf <- foreach(i = 1:nrow(mat_ndvi_anom_rmnp), .combine = "rbind") %dopar% {
    
    # pixel ts
    val <- mat_ndvi_anom_rmnp[i, ]

    # skip current iteration if pixel was previously masked 
    if (all(is.na(val))) {
      return(NA)
    }
    
    dat <- cbind(dat_oni_dmi[, id], val)
    dat <- na.omit(dat)
    
    cor_lag <- Find_Max_CCF(dat[, 1], dat[, 2], lag.max = 8)
    
    if (cor_lag$lag >= 0) {
      mod <- lm(dat[, 1] ~ lag(dat[, 2], cor_lag$lag))
    } else {
      mod <- lm(lag(dat[, 1], abs(cor_lag$lag)) ~ dat[, 2])
    }
    
    p <- summary(mod)$coefficients[2, 4]
    cor_lag_p <- cbind(cor_lag, p)
    
    return(cor_lag_p)
  }
  
  rst_ccf <- rst_lag <- rst_ndvi[[1]]
  rst_ccf[] <- dat_ccf[, 1]
  rst_lag[] <- dat_ccf[, 2]
  rst_ccf_lag <- stack(rst_ccf, rst_lag)
  
  rst_ccf_lag[dat_ccf$p >= .01] <- NA
  
  return(rst_ccf_lag)
}

reds <- colorRampPalette(brewer.pal(7, "OrRd"))
blues <- brewer.pal(7, "YlGnBu")
cols_div <- colorRampPalette(rev(brewer.pal(5, "PuOr")))

# dem contours
dem <- raster("kilimanjaro/coordinates/coords/DEM_ARC1960_30m_Hemp.tif")
dem <- trim(projectRaster(dem, crs = "+init=epsg:4326"))
dem_flipped <- flip(dem, "y")
x <- coordinates(dem_flipped)[, 1]
y <- coordinates(dem_flipped)[, 2]
z <- dem_flipped[]

p_dem <- levelplot(z ~ x * y, colorkey = FALSE, at = seq(1000, 6000, 1000), 
                   panel = function(...) {
                     panel.smoothconts(zlevs.conts = seq(1000, 5500, 500), 
                                       labels = c(1000, "", 2000, "", 3000, "", 4000, "", 5000, ""),
                                       col = "grey50", ...)
                   })

# oni
p_ndvi_oni_lag <- spplot(trim(projectRaster(ls_ccf[[1]][[2]], crs = "+init=epsg:4326")), 
                         col.regions = cols_div(100), at = -6.5:6.5, 
                         scales = list(draw = TRUE), main = list("Lag (months)"), 
                         sp.layout = list("sp.text", loc = c(37.025, -3.35), 
                                          txt = "a)", font = 2, cex = 1.2))
p_ndvi_oni_lag_dem <- p_ndvi_oni_lag + as.layer(p_dem)
p_ndvi_oni_lag_dem_env <- envinmrRasterPlot(p_ndvi_oni_lag_dem, rot = 0)
p_ndvi_oni_r <- spplot(trim(projectRaster(ls_ccf[[1]][[1]], crs = "+init=epsg:4326")), 
                       col.regions = reds(100), at = seq(.025, .525, .05), 
                       scales = list(draw = TRUE), 
                       sp.layout = list("sp.text", loc = c(37.025, -3.35), 
                                        txt = "b)", font = 2, cex = 1.2))
p_ndvi_oni_r_dem <- p_ndvi_oni_r + as.layer(p_dem)
p_ndvi_oni_r_dem_env <- envinmrRasterPlot(p_ndvi_oni_r_dem, rot = 0)

# dmi
p_ndvi_dmi_lag <- spplot(trim(projectRaster(ls_ccf[[2]][[2]], crs = "+init=epsg:4326")), 
                         col.regions = cols_div(100), at = -6.5:6.5, 
                         scales = list(draw = TRUE), 
                         sp.layout = list("sp.text", loc = c(37.025, -3.35), 
                                          txt = "c)", font = 2, cex = 1.2))
p_ndvi_dmi_lag_dem <- p_ndvi_dmi_lag + as.layer(p_dem)
p_ndvi_dmi_lag_dem_env <- envinmrRasterPlot(p_ndvi_dmi_lag_dem, rot = 0)
p_ndvi_dmi_r <- spplot(trim(projectRaster(ls_ccf[[2]][[1]], crs = "+init=epsg:4326")), 
                       col.regions = reds(100), at = seq(.025, .525, .05), 
                       scales = list(draw = TRUE), 
                       sp.layout = list("sp.text", loc = c(37.025, -3.35), 
                                        txt = "d)", font = 2, cex = 1.2))
p_ndvi_dmi_r_dem <- p_ndvi_dmi_r + as.layer(p_dem)
p_ndvi_dmi_r_dem_env <- envinmrRasterPlot(p_ndvi_dmi_r_dem, rot = 0)

p_comb <- latticeCombineGrid(list(p_ndvi_oni_lag_dem_env, p_ndvi_dmi_lag_dem_env, 
                                  p_ndvi_oni_r_dem_env, p_ndvi_dmi_r_dem_env))

## manuscript version
png(paste0(ch_dir_data, "vis/fig08__ccf_oni_dmi.png"), width = 26 * .75, 
    height = 30 * .75, units = "cm", pointsize = 15, res = 300)
plot.new()

vp0 <- viewport(x = 0, y = .05, width = 1, height = .95, 
                just = c("left", "bottom"), name = "figure.vp")
pushViewport(vp0)
print(p_comb, newpage = FALSE)


vp1 <- viewport(x = .5, y = .065,
                height = 0.07, width = .8,
                just = c("center", "bottom"),
                name = "key.vp")
pushViewport(vp1)
draw.colorkey(key = list(col = reds(100), width = 1, height = .5,
                         #                          at = -1:1, labels = c("-", "0", "+"), 
                         at = seq(.025, .525, .05), 
                         space = "bottom"), draw = TRUE)
grid.text("r", x = 0.5, y = -.05, just = c("centre", "top"), 
          gp = gpar(font = 2, cex = 1))

dev.off()

## standalone version
tiff(paste0(ch_dir_data, "vis/figure_08.tiff"), width = 26 * .75, 
    height = 30 * .75, units = "cm", pointsize = 15, res = 500, compression = "lzw")
plot.new()

vp0 <- viewport(x = 0, y = .05, width = 1, height = .95, 
                just = c("left", "bottom"), name = "figure.vp")
pushViewport(vp0)
print(p_comb, newpage = FALSE)


vp1 <- viewport(x = .5, y = .065,
                height = 0.07, width = .8,
                just = c("center", "bottom"),
                name = "key.vp")
pushViewport(vp1)
draw.colorkey(key = list(col = reds(100), width = 1, height = .5,
                         #                          at = -1:1, labels = c("-", "0", "+"), 
                         at = seq(.025, .525, .05), 
                         space = "bottom"), draw = TRUE)
grid.text("r", x = 0.5, y = -.05, just = c("centre", "top"), 
          gp = gpar(font = 2, cex = 1))

dev.off()
