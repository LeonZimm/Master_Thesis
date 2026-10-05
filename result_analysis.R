setwd("C:/Users/leonz/Documents/Ecotox_Studium/Master_Thesis/Data")

library(data.table)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)
library(sf)

rm_fct <- function(keep, env = parent.frame()) {
  rm(list = setdiff(ls(envir = env), keep), envir = env)
  invisible(gc())
}
table_fct <- function(data, phychem_extra = FALSE) {
  # Substance Sales Data --------------------------------------------
  cat("SALES", "\n\n")
  data_postal <- read_sf("Raw/France_PostalCode/codes_postaux_region.shp")
  data_postal <- st_transform(data_postal, crs = 2154) |> 
    filter(DEP != "20")
  
  obs_area <- st_read("Raw/GIS/Results/zonal_slope_layers.gpkg", layer = "10kmb400")
  obs_area$area_obs <- st_area(obs_area)
  
  intsec <- st_intersection(obs_area, data_postal)
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
  
  data <- data[intsec, on = "site.id", allow.cartesian = TRUE]
  idx <- sales_dt[data, on = c("year", "ID", "cas"), which = TRUE]
  data[, amount := sales_dt[idx, .(amount)]]
  data[, yr_sales := weighted.mean(amount, overlay, na.rm = TRUE), 
       by = .(year, cas, site.id)]
  data <- unique(data, by = c("sample.d", "cas", "site.id"))
  data[is.na(yr_sales), yr_sales := 0]
  data[, c("overlay", "amount", "ID") := NULL]
  data[, month := month(as.Date(sample.d))]
  
  
  rm_fct("data")
  
  # Meteorology Data -------------------------------------------------
  cat("METEROLOGY", "\n\n")
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
  met_8x8[, tot_precip_30d := frollmean(tot_precip, n = 30, na.rm = TRUE), 
          by = id]
  
  lag_dt[, id := NULL]
  met_8x8 <- cbind(met_8x8, lag_dt)
  
  site_coord <- readRDS("Raw/stations.processed.rds")
  site_coord <- st_as_sf(site_coord, coords = c("coord.x", "coord.y"), crs = 2154)
  met_8x8_coord_sf <- st_transform(met_8x8_coord_sf, st_crs(site_coord))
  
  nearest_id <- st_nearest_feature(site_coord, met_8x8_coord_sf)
  site_coord <- st_drop_geometry(site_coord)
  site_coord[, grid_id := nearest_id]
  
  data[site_coord, grid_id := i.grid_id, on = "site.id"]
  rm_fct(c("data", "met_8x8"))
  
  met_cols <- c("id", "DATE", "tot_precip", "mean_temp", "windspeed", "mean_humid", "tot_evapot",
                "pot_evapot", "soil_moist", "runoff", "tot_precip_7d", "tot_precip_14d", "tot_precip_30d",
                "tot_precip_lag1", "tot_precip_lag2", "tot_precip_lag3", "tot_precip_lead1", "tot_precip_lead2", "tot_precip_lead3")
  
  idx <- met_8x8[data, on = c("DATE" = "sample.d", "id" = "grid_id"), which = TRUE]
  data[, (met_cols) := met_8x8[idx, ..met_cols]]
  data[, c("grid_id", "id", "DATE") := NULL]
  
  rm_fct("data")
  
  
  # WWTP --------------------------------------------------------------------
  cat("WWTP", "\n\n")
  wwtp <- st_read("Raw/GIS/Results/wwtp_dist.gpkg") |> 
    st_drop_geometry() |> 
    drop_na() |> 
    setDT()
  setnames(wwtp, old = c("origin_id", "destination_id", "network_cost"),
           new = c("wwtp_id", "site.id", "distance"))
  wwtp <- wwtp[, .(site.id, wwtp_id, distance)]
  wwtp_days <- unique(data[, .(site.id, sample.d)][wwtp, on = .NATURAL])
  
  uwwtp <- st_read("Raw/UWWTD_TreatmentPlants.gpkg") |> 
    st_drop_geometry() |> 
    setDT()
  uwwtp <- uwwtp[, .(OBJECTID, uwwDateClosing, uwwBeginLife, uwwLoadEnteringUWWTP, uwwWasteWaterTreated, uwwPrimaryTreatment, uwwSecondaryTreatment, uwwOtherTreatment)]
  uwwtp[, wwtp_level := ifelse(uwwOtherTreatment == 1, 3, ifelse(uwwSecondaryTreatment == 1, 2, 1))]
  
  uwwtp <- uwwtp[wwtp_days, on = c("OBJECTID" = "wwtp_id")]
  wwtp_dis <- uwwtp[, .SD[uwwBeginLife <= sample.d | is.na(uwwBeginLife)][which.min(distance)],
                    by = .(site.id, sample.d)]
  nr_wwtp <- uwwtp[, .(nr_wwtp = sum(uwwBeginLife <= sample.d | is.na(uwwBeginLife))),
                   by = .(site.id, sample.d)]
  wwtp_para <- wwtp_dis[nr_wwtp, on = c("site.id", "sample.d")]
  uwwtp <- uwwtp[wwtp_para, on = .NATURAL, nomatch = NULL]
  setnames(uwwtp, c("uwwLoadEnteringUWWTP", "uwwWasteWaterTreated", "distance"), c("ww_load", "ww_treated", "wwtp_dis"))
  
  uwwtp <- uwwtp[, .(site.id, sample.d, nr_wwtp, wwtp_level, wwtp_dis, ww_treated, ww_load)]
  data <- uwwtp[data, on = .(site.id, sample.d)]
  
  na_fill <- c(wwtp_dis = -1, wwtp_level = 0, nr_wwtp = 0, ww_load = 0, ww_treated = 0)
  for (col in names(na_fill)) {
    set(data, which(is.na(data[[col]])), col, na_fill[[col]])
  }
  
  rm_fct("data")
  
  
  # PhysicoChemcial Data ----------------------------------------------------
  cat("PHYSICOCHEMICAL", "\n")
  if (phyhem_extry == TRUE) {
    phychem_extra <- read.xlsx("Raw/PhysicoChemical.xlsx", sep = ";", sheet = "Tabelle3") |> 
      select(!c("dr_soil_typ", "name", "rel_ab", "ges", "p.type")) |> 
      setDT()
    data <- phychem_extra[data, on = "cas"]
  } else {
    phychem <- read.csv("Raw/PhysicoChemical.csv", sep = ";") |> 
      select(!c("name", "dr_soil_typ")) |> 
      setDT()
    data <- phychem[data, on = "cas"]
  }
  
  rm_fct("data")
  
  # Environmental Data --------------------------------------------------
  cat("ENVIRONMENT", "\n\n")
  env_sa10b400 <- readRDS("Raw/spatial.sa10.b400.processed.rds")
  names_env <- names(env_sa10b400)
  rm_names <- grepl("rz.", names_env)
  rm_names <- names_env[rm_names]
  env_sa10b400[, c(rm_names, "fid", "type", "start") := NULL]
  
  data[, (names(env_sa10b400)) := env_sa10b400[.SD, on = c("site_id" = "site.id"), .SD, .SDcols =  names(env_sa10b400)]]
  data[, site_id := NULL]
  
  rm_fct("data")
  
  # Line Slope --------------------------------------------------------------
  cat("LINE SLOPE", "\n\n")
  linslp_layer <- st_layers("Raw/GIS/Results/line_slope_layers.gpkg")$name
  linslp_lst <- lapply(linslp_layer, function(l) {
    
    dt <- st_read("Raw/GIS/Results/line_slope_layers.gpkg", layer = l) |>
      st_drop_geometry() |>
      setDT()
    
    dt[, slope_ratio := slope/length]  
    min_length <- 10
    floor_noise <- 0.7 * sqrt(2)                           
    dt <- dt[slope_ratio <= floor_noise/length & length >= min_length]          
    
    dt <- dt[, {
      lw <- sum(length * slope_ratio) / sum(length)               
      sd <- sqrt(sum(length * (slope_ratio - lw)^2) / sum(length))
      .(slp = lw, sd = sd)
    }, by = site.id]
    
    setnames(dt, old = c("slp", "sd"), new = c(paste0("linslp_", l),
                                               paste0("linslp_", l, "_sd")))
    
  })
  
  linslp <- mergelist(linslp_lst, on = "site.id", how = "full")
  data <- linslp[data, on = "site.id"]
  
  rm_fct("data")
  
  
  # Riparian Area Slope --------------------------------------------------------------
  cat("RIPARiAN AREA SLOPE", "\n\n")
  ripslp_layer <- st_layers("Raw/GIS/Results/rip_slope_layers.gpkg")$name
  ripslp <- lapply(ripslp_layer, function(l) {
    
    dt <- st_read("Raw/GIS/Results/rip_slope_layers.gpkg", layer = l) |>
      st_drop_geometry() |>
      setDT()
    
    dt[, slope_ratio := slope/100]  
    min_length <- 7.5
    dt <- dt[length >= min_length]
    
    dt <- dt[, {
      lw <- sum(length * slope_ratio) / sum(length)               
      sd <- sqrt(sum(length * (slope_ratio - lw)^2) / sum(length))
      .(slp = lw, sd = sd)
    }, by = site.id]
    
    setnames(dt, old = c("slp", "sd"), new = c(paste0("ripslp_", l),
                                               paste0("ripslp_", l, "_sd")))
    
  })
  
  ripslp <- mergelist(ripslp, on = "site.id", how = "full")
  data <- ripslp[data, on = "site.id"]
  
  rm_fct("data")
  
  
  # Zonal Slope --------------------------------------------------------------
  cat("ZONAL SLOPE", "\n\n")
  zonslp_layer <- st_layers("Raw/GIS/Results/zonal_slope_layers.gpkg")$name
  zonslp <- lapply(zonslp_layer, function(l) {
    
    dt <- st_read("Raw/GIS/Results/zonal_slope_layers.gpkg", layer = l) |>
      st_drop_geometry() |>
      setDT()
    dt <- dt[, .(site.id, X_mean, X_stdev)]
    
    setnames(dt, old = c("X_mean", "X_stdev"), new = c(paste0("zonslp_", l),
                                                       paste0("zonslp_", l, "_sd")))
  })
  
  zonslp <- mergelist(zonslp, on = "site.id", how = "full")
  data <- zonslp[data, on = "site.id"]
  
  rm_fct("data")
  
  
  
  # Riparian Land Use -------------------------------------------------------
  cat("RIPARIAN LAND USE", "\n\n")
  ripagri_layer <- st_layers("Raw/GIS/Results/rip_stats_layers.gpkg")$name
  ripagri <- lapply(ripagri_layer, function(l) {
    
    dt <- st_read("Raw/GIS/Results/rip_stats_layers.gpkg", layer = l) |>
      st_drop_geometry() |>
      setDT()
    
    dt[, perc := HISTO_1/(HISTO_1 + HISTO_0)]
    dt <- dt[, .(site.id, perc)]
    
    setnames(dt, old = "perc", new = paste0("rip_agri_", l))
  })
  
  ripagri <- mergelist(ripagri, on = "site.id", how = "full")
  data <- ripagri[data, on = "site.id"]
  
  rm_fct("data")
  
  
  # Width -------------------------------------------------------------------
  cat("WIDTH", "\n\n")
  width_layer <- st_layers("Raw/GIS/Results/line_layers_topage.gpkg")$name
  width_layer <- width_layer[grepl("split_", width_layer) & !grepl("_topagePoints", width_layer)]
  
  width_lst <- lapply(width_layer, function(l) {
    
    dt <- st_read("Raw/GIS/Results/line_layers_topage.gpkg", layer = l) |>
      st_drop_geometry() |>
      setDT()
    
    dt <- dt[, {
      wm <- weighted.mean(width, length)              
      sd <- sqrt(sum(length * (width - wm)^2) / sum(length))
      .(wid = wm, sd = sd)
    }, by = site.id]
    
    l_part <- sub(".*_", "", l)
    setnames(dt, old = c("wid", "sd"), new = c(paste0("wid_", l_part),
                                               paste0("wid_", l_part, "_sd")))
    
  })
  
  width <- mergelist(width_lst, on = "site.id", how = "full")
  
  site_wid <- st_read("Raw/GIS/Results/line_layers_topage.gpkg", layer = "width_at_site") |> 
    select(site.id, width) |> 
    st_drop_geometry() |>
    setDT()
  width <- site_wid[width, on = "site.id"]
  data <- width[data, on = "site.id"]
  
  rm_fct("data")
  
  
  # Length ------------------------------------------------------------------
  cat("LENGTH", "\n\n")
  site_len <- st_read("Raw/GIS/Results/line_layers_topage.gpkg", layer = "merged_10km") |> 
    select(site.id, length) |> 
    st_drop_geometry() |>
    setDT()
  data <- site_len[data, on = "site.id"]
  
  rm_fct("data")
  
  return(data)
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

# Substance Sales Data --------------------------------------------
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



# Prediction on Unseen Substances -----------------------------------------
phychem_extra <- read.xlsx("Raw/PhysicoChemical.xlsx", sep = ";", sheet = "Tabelle3") |> 
  select(!c("dr_soil_typ")) |> 
  setDT()
parameter <- names(phychem_extra[, 6:ncol(phychem_extra)])
p_type <- c("H", "F", "I")
res <- phychem_extra[, .(name, p.type)]

for (p in parameter) {
  out     <- rep("", nrow(phychem_extra))
  rng_gen <- range(fra[[p]], na.rm = TRUE)   # general range, all p.types
  
  for (t in unique(phychem_extra$p.type)) {
    idx <- which(phychem_extra$p.type == t)
    rng <- range(fra[p.type == t][[p]], na.rm = TRUE)   # type-specific range
    v   <- phychem_extra[[p]][idx]
    
    out[idx] <- fcase(
      is.na(v),                          "NA",
      v < rng_gen[1] | v > rng_gen[2],   "X",   # outside general range (and therefore also the type range)
      v < rng[1]     | v > rng[2],       ".",   # outside type range only
      default = ""
    )
  }
  
  res[, (p) := out]
}
print(res, nrows = Inf)



france <- readRDS("Raw/data.france.rds")
france[, c("sample.id", "chem.code", "chemical", "c.class") := NULL]
france <-  france[year %in% 2010:2023 & p.type %in% c("H", "I", "F"),]
france[ndq == "ND", det := "ND"]
france[ndq != "ND", det := "D"]
france[, det := as.factor(det)]
france <- france[, -c("ndq")]
france <- france[cas %in% phychem_extra$cas]
france <- table_fct(france, phychem)

phychem <- readRDS("results/phychem.rds")
phychem_model <- phychem$models$opt$model

