### environmental stuff

## working directory
Orcs::setwdOS()

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

gsod_sub <- gsodstations[gsodstations$STATION.NAME == "MAKINDU", ]
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

fls_out <- paste0("prcp_mnth_", paste(gsod_sub$USAF, "1982_2011", sep = "_"), ".csv")
write.csv(prcp_nai_mnth, paste0(dir_gsod, "/", fls_out), row.names = FALSE)

# seasonality
tmp_id_st <- grep("1982-01-01", prcp_nai_mnth[, 1])
tmp_id_nd <- grep("2011-12-01", prcp_nai_mnth[, 1])

x <- vectorHarmonics(prcp_nai_mnth[tmp_id_st:tmp_id_nd, 2], 
                     st = c(j, 1), nd = c(k, 12))
x[x < 0] <- 0
harm_prcp_nth <- data.frame(station = "MAKINDU", 
                            period = "1982-2011", 
                            month = month.abb, x)

# merge data
harm_prcp <- rbind(harm_prcp_nth, harm_prcp_sth)
harm_prcp$month <- factor(harm_prcp$month, levels = month.abb)

# visualize
levels(harm_prcp[, 1]) <- c("Makindu", 
                            "Southern Kilimanjaro region")

cols <- c("grey60", "black")
names(cols) <- c("Makindu", "Southern Kilimanjaro region")

p_ssn <- ggplot(aes(x = month, y = x, group = station, colour = station, 
                    fill = station), data = harm_prcp) + 
  # geom_line(lwd = 1.5) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_colour_manual("", values = cols) + 
  scale_fill_manual("", values = cols) + 
  labs(x = "Month", y = "Precipitation (mm)") + 
  theme_bw() + 
  theme(legend.position = c(1, 1), legend.justification = c(.95, .75), 
        legend.text = element_text(size = 8), legend.key.size = unit(4, "mm"),
        axis.title = element_text(size = 10), 
        axis.text = element_text(size = 8), 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        legend.background = element_rect(fill = "transparent"))

## manuscript version
png(paste0(ch_dir_data, "vis/fig06__rainfall_seasonality.png"), width = 10, 
    height = 9, units = "cm", res = 500)
print(p_ssn)
dev.off()

## standalone version
setEPS()
postscript(paste0(ch_dir_data, "vis/figure_06.eps"), width = 10 * .3937, 
           height = 9 * .3937)
print(p_ssn)
dev.off()
