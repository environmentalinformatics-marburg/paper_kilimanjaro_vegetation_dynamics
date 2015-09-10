### environmental stuff

## working directory
library(Orcs)
setwdOS()

## packages
library(GSODTools)
library(Rsenal)
library(reshape2)
library(foreach)
library(ggplot2)

## input filepath
ch_dir_fig <- "publications/paper/detsch_et_al__ndvi_dynamics/figures/"
ch_dir_data <- paste0(ch_dir_fig, "data/")


### data processing

# start and end
st <- "1982-01"
nd <- "2011-12"

# metoffice data
prcp_kia_moshi <- read.table(paste0(ch_dir_data, "metoffice_1973-2013.csv"), 
                             sep = " ", header = TRUE)
prcp_kia_moshi[, 1] <- as.Date(prcp_kia_moshi[, 1])

id_st <- grep(st, prcp_kia_moshi[, 1])
id_nd <- grep(nd, prcp_kia_moshi[, 1])
prcp_kia_moshi <- prcp_kia_moshi[id_st:id_nd, ]

harm_prcp_sth <- lapply(c(6), function(i) {
  tmp_harm_prcp <- foreach(j = c(1982), k = c(2011), 
                           .combine = "rbind") %do% {
  
  tmp_id_st <- grep(paste0(j, "-01"), prcp_kia_moshi[, 1])
  tmp_id_nd <- grep(paste0(k, "-01"), prcp_kia_moshi[, 1])
                             
  x <- vectorHarmonics(prcp_kia_moshi[tmp_id_st:tmp_id_nd, i], st = c(j, 1), nd = c(k, 12))
  x[x < 0] <- 0
  df_x <- data.frame(station = i, 
                     period = paste(j, k, sep = "-"), 
                     month = month.abb, x)
  
  return(df_x)
}})
harm_prcp_sth <- do.call("rbind", harm_prcp_sth)
harm_prcp_sth$station <- "mean(KIA, Moshi)"

# gsod data
data(gsodstations)
dir_gsod <- paste0(ch_dir_data, "gsod")
harm_prcp_nth <- lapply(c("MAKINDU"), function(i) {
  tmp_harm_prcp <- foreach(j = c(1982), k = c(2011), 
                           .combine = "rbind") %do% {
    gsod_sub <- gsodstations[gsodstations$STATION.NAME == i, ]
    gsod_nai <- dlGsodStations(usaf = gsod_sub$USAF, start_year = 1982, end_year = 2011, 
                               dsn = dir_gsod, unzip = TRUE, rm_gz = TRUE)
    
    # inches to mm conversion
    gsod_nai$PRCP <- gsod_nai$PRCP * 25.4
    
    # monthly aggregation
    prcp_nai_mnth <- aggregate(gsod_nai$PRCP, FUN = function(...) sum(..., na.rm = TRUE),
                               by = list(substr(gsod_nai$YEARMODA, 1, 6)))
    names(prcp_nai_mnth)[1] <- "YEARMODA"
    prcp_nai_mnth[, 1] <- paste0(prcp_nai_mnth[, 1], 01)
    prcp_nai_mnth[, 1] <- as.Date(prcp_nai_mnth[, 1], format = "%Y%m%d")
   
    # output storage
    
    fls_out <- paste0("prcp_mnth_", paste(gsod_sub$USAF, j, k, sep = "_"), ".csv")
    write.csv(prcp_nai_mnth, paste0(dir_gsod, "/", fls_out), row.names = FALSE)
    
    # seasonality
    tmp_id_st <- grep(paste0(j, "-01-01"), prcp_nai_mnth[, 1])
    tmp_id_nd <- grep(paste0(k, "-12-01"), prcp_nai_mnth[, 1])
    
    x <- vectorHarmonics(prcp_nai_mnth[tmp_id_st:tmp_id_nd, 2], 
                         st = c(j, 1), nd = c(k, 12))
    x[x < 0] <- 0
    df_x <- data.frame(station = i, 
                       period = paste(j, k, sep = "-"), 
                       month = month.abb, x)
    return(df_x)
  }
})
harm_prcp_nth <- do.call("rbind", harm_prcp_nth)

# merge data
harm_prcp <- rbind(harm_prcp_nth, harm_prcp_sth)
harm_prcp$month <- factor(harm_prcp$month, levels = month.abb)

# visualize
levels(harm_prcp[, 1]) <- c("Makindu", 
                            "Southern Kilimanjaro region")

cols <- c("grey50", "black")
names(cols) <- c("Makindu", "Southern Kilimanjaro region")

p_ssn <- ggplot(aes(x = month, y = x, group = station, colour = station, 
                    fill = station), data = harm_prcp) + 
  # geom_line(lwd = 1.5) + 
  geom_histogram(stat = "identity", position = "dodge") + 
  scale_colour_manual("", values = cols) + 
  scale_fill_manual("", values = cols) + 
  labs(x = "\nMonth", y = "Precipitation (mm)\n") + 
  theme_bw() + 
  theme(legend.position = c(1, 1), legend.justification = c(1, 1))

## manuscript version
png(paste0(ch_dir_data, "vis/fig06__rainfall_seasonality.png"), width = 16.8, 
    height = 12, units = "cm", pointsize = 15, res = 300)
print(p_ssn)
dev.off()

## standalone version
setEPS()
postscript(paste0(ch_dir_data, "vis/figure_06.eps"), width = 16.8 * .3937, 
           height = 12 * .3937, pointsize = 15)
print(p_ssn)
dev.off()
