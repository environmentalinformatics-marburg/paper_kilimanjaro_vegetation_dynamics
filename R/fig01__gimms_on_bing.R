## working directory
setwd("/media/permanent/publications/paper/detsch_et_al__ndvi_dynamics/figures/data")

## packages
lib <- c("raster", "rgdal", "plyr", "Rsenal", "OpenStreetMap", "ggplot2")
sapply(lib, function(x) library(x, character.only = TRUE))

## gimms 8-km
fls_gimms <- "geo198107-VI3g_crp_utm_wht_aggmax.tif"
rst_gimms <- raster(fls_gimms)
template <- rasterToPolygons(rst_gimms)

template@data$id <- rownames(template@data)
template_points <- fortify(template, region = "id")

## modis 1-km
fls_prd <- list.files("MODIS_ARC/PROCESSED/", pattern = "^MOD14A1", 
                      full.names = TRUE, recursive = TRUE)
rst_prd <- raster(fls_prd)
rst_prd <- crop(rst_prd, kiliAerial(rasterize = TRUE))

rst_prd_crp <- crop(rst_prd, template[c(8, 9, 17, 18), ], snap = "out")
template_1km_crp <- rasterToPolygons(rst_prd_crp)
template_1km_crp <- crop(template_1km_crp, template[c(8, 9, 17, 18), ])
template_1km_crp@data$id <- rownames(template_1km_crp@data)
template_1km_crp_points <- fortify(template_1km_crp, region = "id")

## bing aerial image
ext_gimms <- projectExtent(template, "+init=epsg:4326")

kili.map <- openproj(openmap(upperLeft = c(ymax(ext_gimms), xmin(ext_gimms)), 
                             lowerRight = c(ymin(ext_gimms), xmax(ext_gimms)), 
                             type = "bing", minNumTiles = 12L), 
                             projection = "+init=epsg:21037")

# quadrant margins
ext_tmp <- extent(template)
num_cntr_x <- xmin(ext_tmp) + (xmax(ext_tmp) - xmin(ext_tmp)) / 2
num_cntr_y <- ymin(ext_tmp) + (ymax(ext_tmp) - ymin(ext_tmp)) / 2

## np borders
np_old <- readOGR(dsn = "shp/", 
                  layer = "fdetsch-kilimanjaro-national-park-1420535670531", 
                  p4s = "+init=epsg:4326")
np_old_utm <- spTransform(np_old, CRS("+init=epsg:21037"))

np_old_utm@data$id <- rownames(np_old_utm@data) #join id column to data slot on SpatialLinesDataFrame
np_old_utm_df <- fortify(np_old_utm,region="id") #create data frame from SpatialLinesDataFrame
np_old_utm_df <- join(np_old_utm_df, np_old_utm@data, by="id") #add Turbity information to the data frame object

np_new <- readOGR(dsn = "shp/", layer = "fdetsch-kilimanjaro-1420532792846", 
                  p4s = "+init=epsg:4326")
np_new_utm <- spTransform(np_new, CRS("+init=epsg:21037"))

np_new_utm@data$id <- rownames(np_new_utm@data) #join id column to data slot on SpatialLinesDataFrame
np_new_utm_df <- fortify(np_new_utm,region="id") #create data frame from SpatialLinesDataFrame
np_new_utm_df <- join(np_new_utm_df, np_new_utm@data, by="id") #add Turbity information to the data frame object

## visualization
png("vis/fig01__map_incl_borders.png", units = "cm", width = 20, 
    height = 16, res = 300, pointsize = 18)
autoplot(kili.map) + 
  geom_polygon(aes(long, lat), data = np_new_utm_df, colour = "grey75", 
               lwd = 1.1, fill = "transparent") + 
  geom_polygon(aes(long, lat), data = np_old_utm_df, colour = "grey75", 
               lwd = 1, linetype = "dotted", fill = "transparent") + 
  geom_polygon(aes(long, lat, group = group), template_points, lwd = .5, 
               fill = "transparent", colour = "black") + 
  geom_polygon(aes(long, lat, group = group), template_1km_crp_points, lwd = .3,
               fill = "transparent", colour = "black") + 
  geom_hline(aes(yintercept = num_cntr_y), colour = "black", lty = "dashed", 
             lwd = .7) +
  geom_vline(aes(xintercept = num_cntr_x), colour = "black", lty = "dashed", 
             lwd = .7) + 
  scale_x_continuous(breaks = seq(290000, 350000, 20000), expand = c(.001, 0)) + 
  theme(axis.title.x = element_text(size = rel(1.1)), 
        axis.text.x = element_text(size = rel(.9), colour = "black"), 
        axis.title.y = element_text(size = rel(1.1)), 
        axis.text.y = element_text(size = rel(.9), colour = "black"), 
        text = element_text(family = "Arial", colour = "black"))
dev.off()
