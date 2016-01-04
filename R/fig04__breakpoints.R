### environmental stuff

## working directory
library(Orcs)
setwdOS()

## packages
library(zoo)
library(bfast)
library(latticeExtra)
library(reshape2)

## functions
source("repositories/magic/dfg_for_kilimanjaro/ndvi_kilimanjaro/src/gimms/mergeBfast.R")
source("repositories/magic/dfg_for_kilimanjaro/ndvi_kilimanjaro/src/gimms/importOni.R")

## parallelization
library(doParallel)
registerDoParallel(cl <- makeCluster(3))

## input filepath
ch_dir_fig <- "publications/paper/detsch_et_al__ndvi_dynamics/figures/"
ch_dir_data <- paste0(ch_dir_fig, "data/")


### data processing

st <- as.Date("1982-01-01")
nd <- as.Date("2011-12-01")

# upper-left quadrant
fls_prd_ul <- paste0(ch_dir_data, "gimms_ndvi3g_dwnscl_8211_0_0.tif")
rst_prd_ul <- stack(fls_prd_ul)

fls_prd_mk_ul <- paste0(ch_dir_data, "gimms_ndvi3g_dwnscl_8211_mk001_0_0.tif")
rst_prd_mk_ul <- raster(fls_prd_mk_ul)

id_rm <- which((is.na(rst_prd_mk_ul[])) | (rst_prd_mk_ul[] > -.25))
rst_prd_ul[id_rm] <- NA

df_prd_ul <- as.data.frame(rst_prd_ul)

md_prd_ul <- apply(df_prd_ul, 2, FUN = median, na.rm = TRUE)
ts_prd_ul <- ts(md_prd_ul, start = c(1982, 1), end = c(2011, 12), 
                frequency = 12)

bfast_prd_ul <- bfast(ts_prd_ul, max.iter = 10, season = "harmonic", 
                      hpc = "foreach")

bp_id <- bfast_prd_ul$output[[2]]$bp.Vt$breakpoints
bp_months <- seq(st, nd, "month")[bp_id]

df_bfast_prd_ul <- mergeBfast(bfast_prd_ul)
df_bfast_prd_ul$time <- seq(st, nd, "month")
mlt_bfast_prd_ul <- melt(df_bfast_prd_ul, id.vars = "time", 
                         variable.name = "component")
levels(mlt_bfast_prd_ul$component) <- c("Input time series", paste(c("Seasonal", 
                                        "Trend", "Remainder"), "component"))

# enso data
df_oni <- importOni(file = "/media/fdetsch/XChange/kilimanjaro/gimms3g/gimms3g/data/oni/enso_and_iod.csv")
df_oni_high <- subset(df_oni, ONI >= 1)
ls_oni_high <- split(df_oni_high, as.factor(df_oni_high$Season))

ls_oni_max <- lapply(ls_oni_high, function(i) {
  int_id_max <- which.max(i$ONI)
  return(i[int_id_max, ])
})
df_oni_max <- do.call("rbind", ls_oni_max)
df_oni_max <- df_oni_max[c(6, 7, 9), ]

# set back date to july
ls_oni_max_ssn_st <- strsplit(df_oni_max$Season, "-")
ch_oni_max_ssn_st <- sapply(ls_oni_max_ssn_st, "[[", 1)
ch_oni_max_ssn_st <- paste0(ch_oni_max_ssn_st, "-07-01")
dt_oni_max_ssn_st <- as.Date(ch_oni_max_ssn_st)

df_oni_max$Date <- dt_oni_max_ssn_st

# labels
ls_oni_max_lbl <- strsplit(df_oni_max$Season, "-")
mat_oni_max_lbl <- sapply(1:2, function(i) {
  ch_yr4 <- sapply(ls_oni_max_lbl, "[[", i)
  ch_yr2 <- substr(ch_yr4, 3, 4)
  return(ch_yr2)
})
ch_oni_max_lbl <- paste(mat_oni_max_lbl[, 1], mat_oni_max_lbl[, 2], sep = "/")

## manuscript version
old_theme <- lattice.options()
old_theme$layout.widths$left.padding
old_theme$layout.widths$left.padding$units <- "points"
old_theme$layout.widths$right.padding$x <- 0
old_theme$layout.widths$right.padding$units <- "points"
old_theme$layout.heights$top.padding$x <- 0
old_theme$layout.heights$top.padding$units <- "points"
old_theme$layout.heights$bottom.padding$x <- 0
old_theme$layout.heights$bottom.padding$units <- "points"

png(paste0(ch_dir_data, "vis/fig04__breakpoints.png"), width = 14.5, height = 12, 
    units = "cm", res = 500)
xyplot(value ~ time | component, data = mlt_bfast_prd_ul, layout = c(1, 4),
         xlab = "Time (months)", ylab = NULL, 
         as.table = TRUE, scales = list(y = list(relation = "free", rot = 0)), 
         panel = function(x, y) {
           panel.xyplot(x, y, col = "black", type = "l", lwd = 2)
         }, par.settings = list(strip.background = list(col = "lightgrey"), 
                                axis.text = list(cex = 0.65),
                                par.xlab.text = list(cex = 0.8),
                                par.ylab.text = list(cex = 0.8)), 
       lattice.options = old_theme, 
       strip = strip.custom(par.strip.text = list(cex = .7)))

trellis.focus(name = "panel", column = 1, row = 3)
panel.abline(v = bp_months, lty = 3, col = "red", lwd = 2.5)
panel.text(x = bp_months, y = c(.5, .5, .625), labels = paste0("\n", as.yearmon(bp_months)), 
           srt = 90, col = "red", cex = .5)
trellis.unfocus()

for (i in c(2, 4)) {
  trellis.focus(name = "panel", column = 1, row = i)
  panel.abline(h = 0, lty = 2, col = "grey50")
  trellis.unfocus()
}

jnk <- foreach (i = 1:4, label = c("a)", "b)", "c)", "d)")) %do%  {
  trellis.focus(name = "panel", column = 1, row = i)
  y_min <- current.panel.limits()$ylim[1]
  y_max <- current.panel.limits()$ylim[2]
  y_pos <- y_max - 0.15 * (y_max - y_min)
  panel.text(x = as.Date("1981-01-01"), y = y_pos, labels = label, cex = .75, 
             fontface = "bold")
  trellis.unfocus()
}

# major enso events
trellis.focus(name = "panel", column = 1, row = 1)
panel.abline(v = df_oni_max$Date, lty = 2, col = "darkgreen", lwd = 2.5)
panel.text(x = df_oni_max$Date, y = rep(.415, length(df_oni_max$Date)), cex = .5,
           labels = paste0("\n", ch_oni_max_lbl), srt = 90, col = "darkgreen")
trellis.unfocus()

dev.off()

## standalone version
setEPS()
postscript(paste0(ch_dir_data, "vis/figure_04.eps"), width = 14.5*.3937, 
           height = 12*.3937)

xyplot(value ~ time | component, data = mlt_bfast_prd_ul, layout = c(1, 4),
       xlab = "Time (months)", ylab = NULL, 
       as.table = TRUE, scales = list(y = list(relation = "free", rot = 0)), 
       panel = function(x, y) {
         panel.xyplot(x, y, col = "black", type = "l", lwd = 2)
       }, par.settings = list(strip.background = list(col = "lightgrey"), 
                              axis.text = list(cex = 0.65),
                              par.xlab.text = list(cex = 0.8),
                              par.ylab.text = list(cex = 0.8)), 
       lattice.options = old_theme, 
       strip = strip.custom(par.strip.text = list(cex = .7)))

trellis.focus(name = "panel", column = 1, row = 3)
panel.abline(v = bp_months, lty = 3, col = "red", lwd = 2.5)
panel.text(x = bp_months, y = c(.5, .5, .625), labels = paste0("\n", as.yearmon(bp_months)), 
           srt = 90, col = "red", cex = .5)
trellis.unfocus()

for (i in c(2, 4)) {
  trellis.focus(name = "panel", column = 1, row = i)
  panel.abline(h = 0, lty = 2, col = "grey50")
  trellis.unfocus()
}

jnk <- foreach (i = 1:4, label = c("a)", "b)", "c)", "d)")) %do%  {
  trellis.focus(name = "panel", column = 1, row = i)
  y_min <- current.panel.limits()$ylim[1]
  y_max <- current.panel.limits()$ylim[2]
  y_pos <- y_max - 0.15 * (y_max - y_min)
  panel.text(x = as.Date("1981-01-01"), y = y_pos, labels = label, cex = .75, 
             fontface = "bold")
  trellis.unfocus()
}

# major enso events
trellis.focus(name = "panel", column = 1, row = 1)
panel.abline(v = df_oni_max$Date, lty = 2, col = "darkgreen", lwd = 2.5)
panel.text(x = df_oni_max$Date, y = rep(.415, length(df_oni_max$Date)), cex = .5,
           labels = paste0("\n", ch_oni_max_lbl), srt = 90, col = "darkgreen")
trellis.unfocus()

dev.off()

## deregister parallel backend
stopCluster(cl)
