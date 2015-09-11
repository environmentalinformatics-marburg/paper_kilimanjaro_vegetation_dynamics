### environmental stuff

## working directory
library(Orcs)
setwdOS()

## packages
lib <- c("raster", "rgdal", "Kendall", "doParallel", "ggplot2", "rgeos", 
         "reshape2", "plyr", "gridExtra", "remote", "grid")
jnk <- sapply(lib, function(x) library(x, character.only = TRUE, quietly = TRUE))

## functions
source("repositories/magic/dfg_for_kilimanjaro/fire_ndvi/src/uniqueFires.R")

## parallelization
registerDoParallel(cl <- makeCluster(4))

## input filepath
ch_dir_fig <- "publications/paper/detsch_et_al__ndvi_dynamics/figures/"
ch_dir_data <- paste0(ch_dir_fig, "data/")
ch_dir_ext <- "/media/fdetsch/XChange/kilimanjaro/ndvi/"


### data processing

## modified gimms ndvi

# import data
fls_ndvi1km <- paste0(ch_dir_data, "gimms_ndvi3g_dwnscl_8211.tif")
rst_ndvi1km <- stack(fls_ndvi1km)
# (deseasoned)
fls_ndvi1km_dsn <- paste0(ch_dir_data, "gimms_ndvi3g_dwnscl_8211_dsn.tif")
rst_ndvi1km_dsn <- stack(fls_ndvi1km_dsn)

# corresponding months
ch_months <- formatC(1:12, width = 2, flag = "0")
ls_ndvi1km_months <- lapply(1982:2011, function(i) {
  paste0(i, ch_months)
})
ch_ndvi1km_months <- do.call("c", ls_ndvi1km_months)

# prepare 82-00 matrix for subsequent analysis
int_id_dec00 <- grep("200012", ch_ndvi1km_months)
rst_ndvi1km_8200 <- rst_ndvi1km[[1:int_id_dec00]]
mat_ndvi1km_8200 <- as.matrix(rst_ndvi1km_8200)
# (deseasoned)
rst_ndvi1km_dsn_8200 <- rst_ndvi1km_dsn[[1:int_id_dec00]] 
mat_ndvi1km_dsn_8200 <- as.matrix(rst_ndvi1km_dsn_8200)


## modis fire
# monthly fire data
fls_agg1m <- list.files(paste0(ch_dir_ext, "data/md14a1/low/aggregated"), 
                        pattern = "^aggsum_md14a1", full.names = TRUE)
rst_agg1m <- stack(fls_agg1m)

# dates
fire_months <- substr(basename(fls_agg1m), 15, 20)
fire_years <- substr(fire_months, 1, 4)

# unique fire cells
unique_fires <- uniqueFires(rst_agg1m, n_cores = 3)
unique_fires$month <- fire_months[unique_fires$id]
unique_fires$year <- fire_years[unique_fires$id]

# limit to 82-11
int_id_11 <- grep("2011", unique_fires$month)
int_id_11 <- int_id_11[length(int_id_11)]
unique_fires_8211 <- unique_fires[1:int_id_11, ]


## np

# initial np boundaries 
spy_np <- readOGR(dsn = paste0(ch_dir_data, "shp/"), 
                  layer = "fdetsch-kilimanjaro-national-park-1420535670531", 
                  p4s = "+init=epsg:4326")
spy_np <- spTransform(spy_np, CRS = "+init=epsg:21037")

# 1-km modis fire cells inside np
rst_tmp <- rst_agg1m[[1]]
df_np_cells <- extract(rst_tmp, spy_np, cellnumbers = TRUE, df = TRUE)
df_np_cells <- df_np_cells[, c(2, 3)]

# remove cells intersecting np boundary
# spy_np_lines <- as(spy_np, "SpatialLines")
# log_intersects <- foreach(i =1:ncell(rst_tmp), .packages = lib, 
#                           .combine = "c") %dopar% {
#   tmp_rst <- rst_tmp
#   tmp_rst[][-i] <- NA
#   tmp_shp <- rasterToPolygons(tmp_rst)
#   return(gIntersects(spy_np_lines, tmp_shp))
# }
# int_intersects <- which(log_intersects)
# save(int_intersects, file = "data/RData/gIntersects_md14a1_npold.RData")
load(paste0(ch_dir_ext, "data/RData/gIntersects_md14a1_npold.RData"))

log_np_cells_intersect <- df_np_cells[, 1] %in% int_intersects
df_np_cells <- df_np_cells[!log_np_cells_intersect, ]

# rst_tmp[df_np_cells[, 1]] <- 100
# plot(rst_tmp)
# plot(spy_np, col = "transparent", add = TRUE)


## fires inside np

log_np_cells_fire <- unique_fires_8211$cell %in% df_np_cells[, 1]
unique_fires_np <- unique_fires_8211[log_np_cells_fire, ]

# pixel-based mk before and after fire
# df_mk_stats <- foreach(i = 1:nrow(unique_fires_np), .packages = lib, 
#                        .combine = "rbind") %dopar% {
#   int_fire_cell <- unique_fires_np[i, "cell"]
#   ch_fire_month <- unique_fires_np[i, "month"]
#   
#   int_id_fire <- grep(ch_fire_month, ch_ndvi1km_months)
#   
#   # pre fire
#   int_id_prefire <- int_id_fire - 1
#   rst_ndvi1km_prefire <- rst_ndvi1km[[(int_id_dec00+1):int_id_prefire]]
#   mat_ndvi1km_prefire <- as.matrix(rst_ndvi1km_prefire)
#   mat_ndvi1km_prefire <- cbind(mat_ndvi1km_8200, mat_ndvi1km_prefire)
#   num_ndvi1km_prefire <- mat_ndvi1km_prefire[int_fire_cell, ]
#   num_ndvi1km_prefire_mu <- mean(num_ndvi1km_prefire, na.rm = TRUE)
#   num_ndvi1km_prefire_md <- median(num_ndvi1km_prefire, na.rm = TRUE)
#   num_ndvi1km_prefire_mx <- max(num_ndvi1km_prefire, na.rm = TRUE)
#   # (deseasoned)
#   rst_ndvi1km_dsn_prefire <- rst_ndvi1km_dsn[[(int_id_dec00+1):int_id_prefire]]
#   mat_ndvi1km_dsn_prefire <- as.matrix(rst_ndvi1km_dsn_prefire)
#   mat_ndvi1km_dsn_prefire <- cbind(mat_ndvi1km_dsn_8200, mat_ndvi1km_dsn_prefire)
#   num_ndvi1km_dsn_prefire <- mat_ndvi1km_dsn_prefire[int_fire_cell, ]
#   mk_ndvi1km_dsn_prefire <- MannKendall(num_ndvi1km_dsn_prefire)
# 
#   # post fire
#   int_id_postfire <- int_id_fire + 1
#   rst_ndvi1km_postfire <- rst_ndvi1km[[int_id_postfire:nlayers(rst_ndvi1km)]]
#   mat_ndvi1km_postfire <- as.matrix(rst_ndvi1km_postfire)
#   num_ndvi1km_postfire <- mat_ndvi1km_postfire[int_fire_cell, ]
#   num_ndvi1km_postfire_mu <- mean(num_ndvi1km_postfire, na.rm = TRUE)
#   num_ndvi1km_postfire_md <- median(num_ndvi1km_postfire, na.rm = TRUE)
#   num_ndvi1km_postfire_mx <- max(num_ndvi1km_postfire, na.rm = TRUE)
#   # (deseasoned)
#   rst_ndvi1km_dsn_postfire <- rst_ndvi1km_dsn[[(int_id_dec00+1):int_id_postfire]]
#   mat_ndvi1km_dsn_postfire <- as.matrix(rst_ndvi1km_dsn_postfire)
#   mat_ndvi1km_dsn_postfire <- cbind(mat_ndvi1km_dsn_8200, mat_ndvi1km_dsn_postfire)
#   num_ndvi1km_dsn_postfire <- mat_ndvi1km_dsn_postfire[int_fire_cell, ]
#   mk_ndvi1km_dsn_postfire <- MannKendall(num_ndvi1km_dsn_postfire)
#   
#   # merge
#   data.frame(month = ch_fire_month, cell = int_fire_cell, 
#              tau_pre = mk_ndvi1km_dsn_prefire$tau, p_pre = mk_ndvi1km_dsn_prefire$sl, 
#              mu_pre = num_ndvi1km_prefire_mu, md_pre = num_ndvi1km_prefire_md, 
#              mx_pre = num_ndvi1km_prefire_mx, 
#              tau_post = mk_ndvi1km_dsn_postfire$tau, p_post = mk_ndvi1km_dsn_postfire$sl, 
#              mu_post = num_ndvi1km_postfire_mu, md_post = num_ndvi1km_postfire_md, 
#              mx_post = num_ndvi1km_postfire_mx)
# }
# 
# save(df_mk_stats, file = "data/RData/Stats_per_fire_inside_np_dsn.RData")
load(paste0(ch_dir_ext, "data/RData/Stats_per_fire_inside_np_dsn.RData"))

# visualization
df_mk_stats_mlt <- melt(df_mk_stats, id.vars = c(1, 2))
df_mk_stats_mlt$param <- sapply(strsplit(as.character(df_mk_stats_mlt$variable), "_"), "[[", 1)
df_mk_stats_mlt$period <- sapply(strsplit(as.character(df_mk_stats_mlt$variable), "_"), "[[", 2)
df_mk_stats_mlt$period <- factor(df_mk_stats_mlt$period, levels = c("pre", "post"))
levels(df_mk_stats_mlt$period) <- c("pre-fire", "post-fire")

prm <- c("mx", "tau")
lbl <- list(expression("NDVI"[EOTmax]), expression("Kendall\'s " ~ tau))
ls_p <- list()

for (i in 1:2) {
  ls_p[[i]] <- 
    ggplot(aes(x = period, y = value), 
           data = subset(df_mk_stats_mlt, param == prm[i])) + 
      geom_boxplot(notch = TRUE, fill = "grey75", lwd = 1) + 
      labs(x = "", y = lbl[[i]]) + 
      theme_bw() + 
      theme(text = element_text(size = 18), panel.grid = element_blank())

  if (i == 2)
    ls_p[[i]] <- ls_p[[i]] + 
      geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "grey75")    
}

## manuscript version
png(paste0(ch_dir_data, "vis/fig07__vegetation_response.png"), width = 16, 
    height = 20, units = "cm", pointsize = 18, res = 300)
grid.newpage()
print(do.call(function(...) grid.arrange(..., ncol = 2), ls_p), newpage = FALSE)
grid.text("a)", x = .185, y = .9, gp = gpar(font = 2, cex = 1.2))
grid.text("b)", x = .71, y = .9, gp = gpar(font = 2, cex = 1.2))
dev.off()

## standalone version
setEPS()
postscript(paste0(ch_dir_data, "vis/figure_07.eps"), width = 16*.3937, 
           height = 20*.3937, pointsize = 18)
grid.newpage()
print(do.call(function(...) grid.arrange(..., ncol = 2), ls_p), newpage = FALSE)
grid.text("a)", x = .185, y = .9, gp = gpar(font = 2, cex = 1.2))
grid.text("b)", x = .71, y = .9, gp = gpar(font = 2, cex = 1.2))
dev.off()
