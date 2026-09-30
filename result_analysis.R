library(data.table)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(sf)

rm_fct <- function(object) {
  rm(list = setdiff(ls(envir = .GlobalEnv), c(object, "rm_fct")), envir = .GlobalEnv)
  invisible(gc())
}
load_sing_chem_fct <- function(cn) {
  file_path <- file.path("results/sing_sub", cn)
  folder_files <- list.files(file_path, pattern = "\\.rds", full.names = TRUE)
  
  info_list <- list()
  
  for (f in folder_files) {
    obj <- readRDS(f)                 # load object
    name <- tools::file_path_sans_ext(basename(f))  # create object name
    info_list[[name]] <- obj
  }
  return(info_list)
}

setwd("C:/Users/leonz/Documents/Ecotox_Studium/Master_Thesis/Data")
fra <- readRDS("Tables_prepared/fra.rds")

conc_col <- c("conc", "ld", "lq")
phychem_col <- c("sw", "kow", "pka", "p", "h", "dr_soil", "dr_wat_sed", "koc", "gus")
pred_names <- names(fra)[names(fra) != "site.id"]
scale_rm <- grepl("_100m|_1km", pred_names)
scale_rm <- pred_names[scale_rm]
rm_col <- c("site.id", "sample.d", "si_mean", scale_rm)

cols_cas <-  setdiff(names(fra), c(rm_col, conc_col, phychem_col))
cols_phychem <- setdiff(names(fra), c(rm_col, conc_col, "cas"))

omni <- readRDS("results/omnibus.rds")
ss <- readRDS("results/sing_sub.rds")
gc()
phychem <- readRDS("results/phychem.rds")
block <- readRDS("results/block_by_station/blocked/rfs.rds")
not_block <- readRDS("results/block_by_station/normal/rfs.rds")
gc()

# Compare Singular with Omnibus Model ---------------------------------
fra <- readRDS("Tables_prepared/fra.rds")
chem <- readRDS("Raw/chemical.processed.rds")
omni <- readRDS("results/omnibus.rds")
all_mod <- readRDS("results/sing_sub.rds")
substances <- unique(fra$cas)
chem_name <- chem[cas %in% substances]$chemical


all_comp <- setNames(vector("list", length(chem_name)), chem_name)

for (cas_nr in substances) {
  
  cn <- chem[cas == cas_nr]$chemical
  # Data Preperation
  data_cas <- fra[complete.cases(fra[, ..cols_cas]), ..cols_cas][cas == cas_nr]
  set.seed(25)
  partition_subset <- createDataPartition(data_cas$det,
                                          p = 0.7,
                                          list = FALSE,
                                          times = 1)
  sub <- data_cas[partition_subset]
  data_test <- data_cas[-partition_subset, ]
  
  set.seed(25)
  partition <- createDataPartition(sub$det, 
                                   p = 0.8,
                                   list = FALSE,
                                   times = 1)
  data_train <- sub[partition]
  data_val <- sub[-partition]
  
  opt_res_list <- IR_eval_fct(omni$models$opt$model, data_val, data_test, 10)
  cf_opt <- confusionMatrix(as.factor(opt_res_list$pred),
                            data_test$det, 
                            positive = "D")
  
  comparison <- list("cf_sing" = all_mod[[cn]]$models$opt$confu_mat,
                     "cf_omni" = cf_opt)
  
  all_comp[[cn]] <- comparison
}

saveRDS(all_comp, "results/omni_ss_comp.rds")


# Features Omni vs. SS ----------------------------------------------------
feat_grps <- as.data.table(read.xlsx("parameter_grouping.xlsx"))
omni <- readRDS("results/omnibus.rds")
ss <- readRDS("results/sing_sub.rds")
feat_dt <- data.table(feature = colnames(fra[, ..cols_cas]))

for (cn in names(ss)) {
  feat <- ss[[cn]]$models$opt$features
  feat_dt[, (cn) := as.integer(feature %in% feat)]
}

feat_dt[, total := rowSums(.SD), .SDcols = names(ss)]
feat_dt[, percent := round(total/30, 2)]
feat_dt[, omni := as.integer(feature %in% omni$rfs$remaining_features)]

total_row <- feat_dt[, lapply(.SD, sum), .SDcols = is.numeric]
total_row[, feature := "Total"]
feat_dt <- rbind(feat_dt, total_row, fill = TRUE)
feat_dt <- feat_grps[feat_dt, on = "feature"]

feat_long <- melt(feat_dt[nrow(feat_dt), 3:32], variable.name = "chem", value.name = "n_feat")


# Compare Relative Occurrence to Sensitivity
chem <- readRDS("Raw/chemical.processed.rds")
substances <- unique(fra$cas)
chem_name <- chem[cas %in% substances]$chemical
confma_all <- readRDS("results/omni_ss_comp.rds")

dt <- data.table(chem = names(confma_all),
                 total_n = sapply(chem_name, function(cn){sum(confma_all[[cn]]$cf_omni$table)}),
                 rel_ab = unname(sapply(chem_name, function(cn){confma_all[[cn]]$cf_omni$byClass[8]})),
                 se_omni = unname(sapply(chem_name, function(cn){confma_all[[cn]]$cf_omni$byClass[1]})),
                 se_sing = unname(sapply(chem_name, function(cn){confma_all[[cn]]$cf_sing$byClass[1]})),
                 sp_omni = unname(sapply(chem_name, function(cn){confma_all[[cn]]$cf_omni$byClass[2]})),
                 sp_sing = unname(sapply(chem_name, function(cn){confma_all[[cn]]$cf_sing$byClass[2]}))
)
num_cols <- names(dt)[sapply(dt, is.numeric)]
dt[, (num_cols) := lapply(.SD, round, digits = 3), .SDcols = num_cols]
dt[chem, p.type := i.p.type, on = c("chem" = "chemical")]
setorder(dt, rel_ab)
val_ord <- dt$chem

dt_long <- melt(dt, 
                id.vars = c("chem", "total_n", "rel_ab", "p.type"), 
                variable.name = "model",
                value.name = "sens")
# dt_long <- feat_long[dt_long, on = "chem"]

# Dumbell Plot
ggplot(dt_long[model %in% c("se_omni", "se_sing")], aes(sens, reorder(chem, rel_ab))) +
  geom_line() +
  geom_point(aes(color = model))

# Relative Frequency Plot
ggplot(dt_long, aes(rel_ab, sens, colour = model)) +
  geom_point() +
  geom_smooth(se = FALSE)

# Grouped Features Plot
feat_grps <- melt(feat_dt[, c(2:32, 35)], id.vars = "group", variable.name = "chem")
feat_grps <- feat_grps[complete.cases(feat_grps) & value == 1]

ggplot(feat_grps, aes(chem, fill = group)) +
  geom_bar() +
  coord_flip()





# Weekly Predictions ------------------------------------------------------
site_sub <- CJ(unique(fra$site.id), unique(fra$cas))
setnames(site_sub, c("V1", "V2"), c("site.id", "cas"))
site_sub <- unique(fra[, -c(2:4, 7:10, 55:76, 102)])[site_sub, on = c("site.id", "cas")]

week_vec <- seq(from = as.Date("2013-01-01"),
                to = as.Date("2023-12-31"),
                by = "week")
week_dt <- CJ(unique(fra$site.id), unique(fra$cas), week_vec)
setnames(week_dt, c("V1", "V2", "week_vec"), c("site.id", "cas", "sample.d"))

# Meteorology Data -------------------------------------------------
met_8x8 <- readRDS("C:/Users/leonz/Documents/Ecotox_Studium/4. Semester/RPC/1. RPC_data_and_calc/1. RPC_raw_data/new_raw_data/met_france8x8grid_0523.rds")
met_8x8[, ":="(LAMBX = LAMBX * 100, 
               LAMBY = LAMBY * 100)]

met_8x8_coord <- unique(met_8x8[, .(LAMBX, LAMBY)])
met_8x8_coord[, id := seq(1, nrow(met_8x8_coord))]
met_8x8_coord_sf <- st_as_sf(met_8x8_coord, coords = c("LAMBX", "LAMBY"), crs = 27572)

met_8x8 <- met_8x8_coord[met_8x8, on = .(LAMBX, LAMBY)]
setnames(met_8x8, 
         old = c("PRELIQ", "T", "FF", "HU", "EVAP", "ETP", "SWI", "RUNC"),
         new = c("tot_precip", "mean_temp", "windspeed", "mean_humid", "tot_evapot", "pot_evapot", "soil_moist", "runoff"))
setorder(met_8x8, id, DATE)

# Create Lagged Variables
vars <- c("tot_precip", "mean_temp", "windspeed")
lag_dt <- met_8x8[, as.data.table(unlist(
  c(lapply(vars, \(v) setNames(data.table::shift(get(v), n = 1:3, type = "lag"), paste0(v, "_lag", 1:3))),
    lapply(vars, \(v) setNames(data.table::shift(get(v), n = 1:3, type = "lead"), paste0(v, "_lead", 1:3)))),
  recursive = FALSE)),
  by = id]
met_8x8[, tot_precip_7d := frollmean(tot_precip, n = 7, na.rm = TRUE)]
met_8x8[, tot_precip_14d := frollmean(tot_precip, n = 14, na.rm = TRUE)]
met_8x8[, tot_precip_30d := frollmean(tot_precip, n = 30, na.rm = TRUE)]

lag_dt[, id := NULL]
met_8x8 <- cbind(met_8x8, lag_dt)

site_coord <- readRDS("Raw/stations.processed.rds")
site_coord <- st_as_sf(site_coord, coords = c("coord.x", "coord.y"), crs = 2154)
met_8x8_coord_sf <- st_transform(met_8x8_coord_sf, st_crs(site_coord))

nearest_id <- st_nearest_feature(site_coord, met_8x8_coord_sf)
site_coord <- st_drop_geometry(site_coord)
site_coord[, grid_id := nearest_id]

week_dt <- week_dt[site_coord, on = "site.id"]
week_dt <- merge.data.table(week_dt, met_8x8[, .(id, DATE, tot_precip, mean_temp, windspeed, mean_humid, tot_evapot,
                                               pot_evapot, soil_moist, runoff, tot_precip_7d, tot_precip_14d, tot_precip_30d,
                                               tot_precip_lag1, tot_precip_lag2, tot_precip_lag3, tot_precip_lead1, tot_precip_lead2, tot_precip_lead3)], 
                           by.x = c("sample.d", "grid_id"), by.y = c("DATE", "id"),
                           all.x = TRUE)
week_dt <- week_dt[, -c(5:10)]
rm_fct(c("week_dt", "fra"))
try <- site_sub[week_dt, on = c("site.id", "cas")]

# 1.4 Substance Sales Data --------------------------------------------
france_postal <- read_sf("Raw/France_PostalCode/codes_postaux_region.shp")
france_postal <- st_transform(france_postal, crs = 2154) |> 
  filter(DEP != "20")

obs_area <- st_read("Raw/GIS/Results/zonal_slope_layers.gpkg", layer = "10kmb400")
obs_area$area_obs <- st_area(obs_area)

intsec <- st_intersection(obs_area, france_postal)
intsec$area_insec <- st_area(intsec)
intsec$overlay <- as.numeric((intsec$area_insec/intsec$area_obs) * 100)
intsec <- intsec[c("site.id", "ID", "overlay")] |> 
  st_drop_geometry() |> 
  setDT()

sales_dt <- readRDS("Raw/sales_data.rds")
sales_dt <- sales_dt[, .(annee, code_postal_acheteur, cas, quantite_substance)]
setnames(sales_dt,
         c("annee", "code_postal_acheteur", "quantite_substance"),
         c("year", "ID", "amount"))
sales_dt[, ":="(ID = as.character(ID),
                amount = as.numeric(amount))]
sales_dt[, amount := sum(amount), by = .(year, ID, cas)]
sales_dt <- unique(sales_dt)

week_dt <- merge.data.table(week_dt, intsec,
                           by = "site.id",
                           all.x = TRUE,
                           allow.cartesian = TRUE)
week_dt <- merge.data.table(week_dt, sales_dt,
                           by = c("year", "ID", "cas"),
                           all.x = TRUE)
week_dt[, yr_sales := weighted.mean(amount, overlay, na.rm = TRUE), 
       by = .(year, cas, site.id)]
week_dt <- unique(week_dt, by = c("sample.d", "cas", "site.id"))
week_dt[is.na(yr_sales), yr_sales := 0]
week_dt[, c("overlay", "amount", "ID") := NULL]

rm_fct(c("france", "fra"))
