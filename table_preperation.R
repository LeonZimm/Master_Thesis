setwd("C:/Users/leonz/Documents/Ecotox_Studium/Master_Thesis/Data")

library(sf)
library(RSQLite)
library(dplyr)
library(tidyr)
library(terra)
library(tmap)
library(data.table)

rm_fct <- function(keep = character(), env = parent.frame()) {
  caller_args <- character()
  if (sys.parent() > 0) {
    caller <- sys.function(sys.parent())
    if (!is.primitive(caller)) caller_args <- names(formals(caller))
  }
  
  keep_all <- union(c(keep, caller_args), "rm_fct")
  rm(list = setdiff(ls(envir = env), keep_all), envir = env)
  invisible(gc())
}
table_fct <- function(data, 
                      phychem_extra = FALSE, 
                      res_anal = FALSE,
                      features = character()) {
  # Substance Sales Data --------------------------------------------
  cat("\n", "SALES", "\n\n")
  data_postal <- read_sf("table_fct/codes_postaux_region.shp")
  data_postal <- st_transform(data_postal, crs = 2154) |> 
    filter(DEP != "20")
  
  obs_area <- st_read("table_fct/zonal_slope_layers.gpkg", layer = "10kmb400")
  obs_area$area_obs <- st_area(obs_area)
  
  intsec <- st_intersection(obs_area, data_postal)
  intsec$area_insec <- st_area(intsec)
  intsec$overlay <- as.numeric((intsec$area_insec/intsec$area_obs) * 100)
  intsec <- intsec[c("site.id", "ID", "overlay")] |> 
    st_drop_geometry() |> 
    setDT()
  
  sales_dt <- readRDS("table_fct/sales_data.rds")
  sales_dt <- sales_dt[cas %in% unique(data$cas)]
  sales_dt <- sales_dt[, .(annee, code_postal_acheteur, cas, quantite_substance)]
  setnames(sales_dt,
           c("annee", "code_postal_acheteur", "quantite_substance"),
           c("year", "ID", "amount"))
  sales_dt[, ":="(ID = as.character(ID),
                  amount = as.numeric(amount))]
  sales_dt[, amount := sum(amount), by = .(year, ID, cas)]
  sales_dt <- unique(sales_dt)
  
  red_sales <- unique(data[, .(site.id, year, cas)]) 
  red_sales <- red_sales[intsec, on = "site.id", allow.cartesian = TRUE, nomatch = NULL]
  red_sales[sales_dt, amount := i.amount, on = .(ID, year, cas)]
  red_sales[, yr_sales := weighted.mean(amount, overlay, na.rm = TRUE), 
       by = .(year, cas, site.id)]
  red_sales <- unique(red_sales, by = c("year", "cas", "site.id"))
  data[, yr_sales := red_sales[.SD, on = .(site.id, year, cas), x.yr_sales]]
  data[is.na(yr_sales), yr_sales := 0]
  data[, month := month(as.Date(sample.d))]
  
  rm_fct()
  
  # Meteorology Data -------------------------------------------------
  cat("\n", "METEROLOGY", "\n\n")
  met_8x8 <- readRDS("table_fct/met_france8x8grid_0523.rds")
  rng <- range(data$sample.d, na.rm = TRUE) + c(-30, 3)
  met_8x8 <- met_8x8[DATE >= rng[1] & DATE <= rng[2]]
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
  vars <- c("tot_precip")
  lag_dt <- met_8x8[, as.data.table(unlist(
    c(lapply(vars, \(v) setNames(data.table::shift(get(v), n = 1:3, type = "lag"), paste0(v, "_lag", 1:3))),
      lapply(vars, \(v) setNames(data.table::shift(get(v), n = 1:3, type = "lead"), paste0(v, "_lead", 1:3)))),
    recursive = FALSE)),
    by = id]
  met_8x8[, tot_precip_7d := frollmean(tot_precip, n = 7, na.rm = TRUE), by = id]
  met_8x8[, tot_precip_14d := frollmean(tot_precip, n = 14, na.rm = TRUE), by = id]
  met_8x8[, tot_precip_30d := frollmean(tot_precip, n = 30, na.rm = TRUE), by = id]
  
  lag_dt[, id := NULL]
  met_8x8 <- cbind(met_8x8, lag_dt)
  
  site_coord <- readRDS("table_fct/stations.processed.rds")
  site_coord <- st_as_sf(site_coord, coords = c("coord.x", "coord.y"), crs = 2154)
  met_8x8_coord_sf <- st_transform(met_8x8_coord_sf, st_crs(site_coord))
  
  nearest_id <- st_nearest_feature(site_coord, met_8x8_coord_sf)
  site_coord <- st_drop_geometry(site_coord)
  site_coord[, grid_id := nearest_id]
  
  data[site_coord, grid_id := i.grid_id, on = "site.id"]
  rm_fct("met_8x8")
  
  met_cols <- c("id", "DATE", "tot_precip", "mean_temp", "windspeed", "mean_humid", "tot_evapot",
                "pot_evapot", "soil_moist", "runoff", "tot_precip_7d", "tot_precip_14d", "tot_precip_30d",
                "tot_precip_lag1", "tot_precip_lag2", "tot_precip_lag3", "tot_precip_lead1", "tot_precip_lead2", "tot_precip_lead3")
  
  idx <- met_8x8[data, on = c("DATE" = "sample.d", "id" = "grid_id"), which = TRUE]
  data[, (met_cols) := met_8x8[idx, ..met_cols]]
  data[, c("grid_id", "id", "DATE") := NULL]
  
  rm_fct()
  
  
  # WWTP --------------------------------------------------------------------
  cat("\n", "WWTP", "\n\n")
  wwtp <- st_read("table_fct/wwtp_dist.gpkg") |> 
    st_drop_geometry() |> 
    drop_na() |> 
    setDT()
  setnames(wwtp, old = c("origin_id", "destination_id", "network_cost"),
           new = c("wwtp_id", "site.id", "distance"))
  wwtp <- wwtp[, .(site.id, wwtp_id, distance)]
  wwtp_days <- unique(data[, .(site.id, sample.d)][wwtp, on = .NATURAL])
  
  uwwtp <- st_read("table_fct/UWWTD_TreatmentPlants.gpkg") |> 
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
  
  rm_fct()
  
  
  # PhysicoChemcial Data ----------------------------------------------------
  cat("\n", "PHYSICOCHEMICAL", "\n")
  if (phychem_extra) {
    phychem <- read.xlsx("table_fct/PhysicoChemical.xlsx", sheet = "Tabelle3") |> 
      select(!c("dr_soil_typ", "name", "rel_ab", "ges", "p.type")) |> 
      setDT()
    data <- phychem[data, on = "cas"]
  } else {
    phychem <- read.xlsx("table_fct/PhysicoChemical.xlsx", sheet = "Tabelle1") |> 
      select(!c("name", "dr_soil_typ")) |> 
      setDT()
    data <- phychem[data, on = "cas"]
  }
  
  rm_fct()
  
  
  
  # Site Joins --------------------------------------------------------------
  if (file.exists("results/site_features.rds")) {
    cat("\n", "SITE TABLE EXISTS", "\n\n")
    site_dt <- readRDS("results/site_features.rds")
  } else{
    # Environmental --------------------------------------------------
    site_dt <- unique(data[, .(site.id)])
    
    cat("\n", "ENVIRONMENT", "\n\n")
    env_sa10b400 <- readRDS("table_fct/spatial.sa10.b400.processed.rds")
    names_env <- names(env_sa10b400)
    rm_names <- grepl("rz.", names_env)
    rm_names <- names_env[rm_names]
    env_sa10b400[, c(rm_names, "fid", "type", "start") := NULL]
    
    cols <- setdiff(names(env_sa10b400), "site_id")
    site_dt[, (cols) := env_sa10b400[.SD, on = c("site_id" = "site.id"), .SD, .SDcols = cols]]
    
    rm_fct("site_dt")
    
    # Line Slope --------------------------------------------------------------
    cat("\n", "LINE SLOPE", "\n\n")
    linslp_layer <- st_layers("table_fct/line_slope_layers.gpkg")$name
    linslp_lst <- lapply(linslp_layer, function(l) {
      
      dt <- st_read("table_fct/line_slope_layers.gpkg", layer = l) |>
        st_drop_geometry() |>
        setDT()
      
      dt[, slope_ratio := slope/length]  
      min_length <- 10
      floor_noise <- 0.7 * sqrt(2)                           
      dt <- dt[slope_ratio <= floor_noise/length & length >= min_length]          
      
      dt <- dt[, {
        lw <- sum(length * slope_ratio, na.rm = TRUE) / sum(length, na.rm = TRUE)               
        sd <- sqrt(sum(length * (slope_ratio - lw)^2, na.rm = TRUE) / sum(length, na.rm = TRUE))
        .(slp = lw, sd = sd)
      }, by = site.id]
      
      setnames(dt, old = c("slp", "sd"), new = c(paste0("linslp_", l),
                                                 paste0("linslp_", l, "_sd")))
      
    })
    
    linslp <- mergelist(linslp_lst, on = "site.id", how = "full")
    cols <- setdiff(names(linslp), "site.id")
    site_dt[, (cols) := linslp[.SD, on = "site.id", .SD, .SDcols = cols]]
    
    rm_fct("site_dt")
    
    
    # Riparian Area Slope --------------------------------------------------------------
    cat("\n", "RIPARiAN AREA SLOPE", "\n\n")
    ripslp_layer <- st_layers("table_fct/rip_slope_layers.gpkg")$name
    ripslp <- lapply(ripslp_layer, function(l) {
      
      dt <- st_read("table_fct/rip_slope_layers.gpkg", layer = l) |>
        st_drop_geometry() |>
        setDT()
      
      dt[, slope_ratio := slope/100]  
      min_length <- 7.5
      dt <- dt[length >= min_length]
      
      dt <- dt[, {
        lw <- sum(length * slope_ratio, na.rm = TRUE) / sum(length, na.rm = TRUE)               
        sd <- sqrt(sum(length * (slope_ratio - lw)^2, na.rm = TRUE) / sum(length, na.rm = TRUE))
        .(slp = lw, sd = sd)
      }, by = site.id]
      
      setnames(dt, old = c("slp", "sd"), new = c(paste0("ripslp_", l),
                                                 paste0("ripslp_", l, "_sd")))
      
    })
    
    ripslp <- mergelist(ripslp, on = "site.id", how = "full")
    cols <- setdiff(names(ripslp), "site.id")
    site_dt[, (cols) := ripslp[.SD, on = "site.id", .SD, .SDcols = cols]]
    
    rm_fct("site_dt")
    
    
    # Zonal Slope --------------------------------------------------------------
    cat("\n", "ZONAL SLOPE", "\n\n")
    zonslp_layer <- st_layers("table_fct/zonal_slope_layers.gpkg")$name
    zonslp <- lapply(zonslp_layer, function(l) {
      
      dt <- st_read("table_fct/zonal_slope_layers.gpkg", layer = l) |>
        st_drop_geometry() |>
        setDT()
      dt <- dt[, .(site.id, X_mean, X_stdev)]
      
      setnames(dt, old = c("X_mean", "X_stdev"), new = c(paste0("zonslp_", l),
                                                         paste0("zonslp_", l, "_sd")))
    })
    
    zonslp <- mergelist(zonslp, on = "site.id", how = "full")
    cols <- setdiff(names(zonslp), "site.id")
    site_dt[, (cols) := zonslp[.SD, on = "site.id", .SD, .SDcols = cols]]
    
    rm_fct("site_dt")
    
    
    # Riparian Land Use -------------------------------------------------------
    cat("\n", "RIPARIAN LAND USE", "\n\n")
    ripagri_layer <- st_layers("table_fct/rip_stats_layers.gpkg")$name
    ripagri <- lapply(ripagri_layer, function(l) {
      
      dt <- st_read("table_fct/rip_stats_layers.gpkg", layer = l) |>
        st_drop_geometry() |>
        setDT()
      
      dt[, perc := HISTO_1/(HISTO_1 + HISTO_0)]
      dt <- dt[, .(site.id, perc)]
      
      setnames(dt, old = "perc", new = paste0("rip_agri_", l))
    })
    
    ripagri <- mergelist(ripagri, on = "site.id", how = "full")
    cols <- setdiff(names(ripagri), "site.id")
    site_dt[, (cols) := ripagri[.SD, on = "site.id", .SD, .SDcols = cols]]
    
    rm_fct("site_dt")
    
    
    # Width -------------------------------------------------------------------
    cat("\n", "WIDTH", "\n\n")
    
    width_layer <- st_layers("table_fct/line_layers_topage.gpkg")$name
    width_layer <- width_layer[grepl("split_", width_layer) & !grepl("_topagePoints", width_layer)]
    
    width_lst <- lapply(width_layer, function(l) {
      
      dt <- st_read("table_fct/line_layers_topage.gpkg", layer = l) |>
        st_drop_geometry() |>
        setDT()
      
      dt <- dt[, {
        wm <- weighted.mean(width, length, na.rm = TRUE)              
        sd <- sqrt(sum(length * (width - wm)^2, na.rm = TRUE) / sum(length, na.rm = TRUE))
        .(wid = wm, sd = sd)
      }, by = site.id]
      
      l_part <- sub(".*_", "", l)
      setnames(dt, old = c("wid", "sd"), new = c(paste0("wid_", l_part),
                                                 paste0("wid_", l_part, "_sd")))
      
    })
    
    width <- mergelist(width_lst, on = "site.id", how = "full")
    
    site_wid <- st_read("table_fct/line_layers_topage.gpkg", layer = "width_at_site") |> 
      select(site.id, width) |> 
      st_drop_geometry() |>
      setDT()
    width <- site_wid[width, on = "site.id"]
    cols <- setdiff(names(width), "site.id")
    site_dt[, (cols) := width[.SD, on = "site.id", mult = "first", .SD, .SDcols = cols]]
    
    rm_fct("site_dt")
    
    
    # Length ------------------------------------------------------------------
    cat("\n", "LENGTH", "\n\n")
    site_len <- st_read("table_fct/line_layers_topage.gpkg", layer = "merged_10km") |> 
      select(site.id, length) |> 
      st_drop_geometry() |>
      setDT()
    cols <- setdiff(names(site_len), "site.id")
    site_dt[, (cols) := site_len[.SD, on = "site.id", .SD, .SDcols = cols]]
    
    rm_fct("site_dt")
    
    saveRDS(site_dt, "results/site_features.rds")
  }
  
  # Final Join --------------------------------------------------------------
  if (res_anal) {
    cols <- setdiff(names(site_dt[, -c("site.id")]), features)
    site_dt[, (cols) := NULL]
    
    cols <- setdiff(names(data[, -c("site.id")]), features)
    data[, (cols) := NULL]
  } 
  
  cols <- setdiff(names(site_dt), "site.id")
  data[, (cols) := site_dt[.SD, on = "site.id", .SD, .SDcols = cols]]
  
  return(data)
}



# Data Preparation ---------------------------------------------------
# France Data
france <- readRDS("table_fct/data.france.rds")
france <-  france[year %in% 2010:2023 & p.type %in% c("H", "I", "F"),]
france[ndq == "ND", det := "ND"]
france[ndq != "ND", det := "D"]
france[, det := as.factor(det)]
france <- france[, -c("ndq")]

# Most abundant pesticides (10 x H, I, F)
sales_dt <- readRDS("table_fct/sales_data.rds")
pest_dt <- france[, .N, by = .(cas, det, p.type)][, rel_ab := N/sum(N), by = cas]
pest_inters <- pest_dt[cas %in% unique(sales_dt$cas)]
pest_inters <- pest_inters[, ges := sum(N), by = cas][det == "D"]
most_pes <- pest_inters[ges >= 100000][order(-rel_ab), .SD[1:10], by = p.type]
france <- france[cas %in% most_pes$cas]
france[, c("sample.id", "chem.code", "chemical", "c.class") := NULL]

sites <- read_sf("table_fct/line_layers_topage.gpkg", 
                 layer = "merged_10km")
france <- france[site.id %in% sites$site.id]
france <- readRDS("Tables_prepared/france_reduced.rds")
france <- table_fct(france) 

# Compare Measured to Estimated Width -------------------------------------
library(epiR)
library(hydroGOF)

site_wid <- st_read("Raw/GIS/Results/line_layers_topage.gpkg", layer = "width_at_site") |> 
  select(site.id, width, RANG) |> 
  st_drop_geometry() |>
  setDT()
nai_wid <- read.csv2("Raw/operation.csv") |> 
  select(CdStationMesureEauxSurface, LargeurPleinBord) |> 
  rename("measured" = "LargeurPleinBord") |> 
  setDT()
wid_comp <- site_wid[nai_wid, on = c("site.id" = "CdStationMesureEauxSurface")]
wid_comp <- wid_comp[, measured := as.numeric(measured)][complete.cases(wid_comp) & measured != 999 & RANG > 0]
setnames(wid_comp, old = c("width", "RANG"), new = c("estimated", "strahler"))

# Remove probable Outliers
whisker_mult <- 3

wid_comp <- wid_comp |> 
  group_by(strahler) |>
  mutate(q1 = quantile(measured, 0.25, na.rm = TRUE),
         q3 = quantile(measured, 0.75, na.rm = TRUE),
         iqr = q3 - q1,
         upper = q3 + whisker_mult * iqr,
         lower = q1 - whisker_mult * iqr,
         is_outlier = measured > upper | measured < lower) |>
  ungroup()

# Inspect Flagged Observations
wid_comp |>
  filter(is_outlier) |>
  arrange(strahler, desc(measured)) |>
  select(site.id, strahler, measured, estimated) |>
  print(n = 50)

wid_comp_clean <- wid_comp |> 
  filter(!is_outlier)

# Validation Metrics
recompute_metrics <- function(data) {
  ccc_res   <- epi.ccc(data$estimated, data$measured, ci = "z-transform", conf.level = 0.95)
  ccc_log   <- epi.ccc(log(data$estimated), log(data$measured), ci = "z-transform", conf.level = 0.95)
  spearman  <- cor(data$estimated, data$measured, method = "spearman")
  rmse_val  <- rmse(data$estimated, data$measured)
  mae_val   <- mae(data$estimated, data$measured)
  pbias_val <- pbias(data$estimated, data$measured)
  nse_val   <- NSE(data$estimated, data$measured)
  
  list(n = nrow(data),
       CCC = ccc_res$rho.c$est,
       CCC_log = ccc_log$rho.c$est,
       Spearman = spearman,
       RMSE = rmse_val,
       MAE = mae_val,
       PBIAS = pbias_val,
       NSE = nse_val)
}
recompute_metrics(wid_comp_clean)


wid_comp_clean |>
  group_by(strahler) |>
  summarise(n = n(),
            RMSE = rmse(estimated, measured),
            MAE = mae(estimated, measured),
            Bias = mean(estimated - measured),
            PBIAS = pbias(estimated, measured),
            .groups = "drop")


# Compare Per-point/Mean Relative Rrror against PBIAS - Test for Bias Homogeneity
wid_comp_clean$rel_error <- (wid_comp_clean$estimated - wid_comp_clean$measured) / wid_comp_clean$measured * 100

median(wid_comp_clean$rel_error, na.rm = TRUE)
pbias(wid_comp_clean$estimated, wid_comp_clean$measured)

wid_comp_clean %>%
  group_by(strahler) %>%
  summarise(
    n = n(),
    median_rel_error = median(rel_error, na.rm = TRUE),
    mean_rel_error   = mean(rel_error, na.rm = TRUE),
    PBIAS = pbias(estimated, measured),
    .groups = "drop"
  )

lims <- range(c(wid_comp_clean$measured, wid_comp_clean$estimated))

ggplot(wid_comp_clean, aes(x = measured, y = estimated)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  coord_equal(xlim = lims, ylim = lims) +
  theme_minimal()

ggplot(wid_comp_clean, aes(x = as.factor(strahler), y = measured)) +
  geom_boxplot() +
  theme_minimal()

cor_table_fct <- function(dataset) {
  library(foreign)
  library(dplyr)
  
  dataset <- dataset |> 
    select(where(is.numeric))
  cor <- cor(dataset) 
  cor[lower.tri(cor, diag = TRUE)] <- NA 
  cor[abs(cor) < 0.70] <- NA 
  
  relevant <- which(abs(cor) >= .70, arr.ind = T)
  vals <- abs(as.matrix(cor)[relevant])
  
  final <- tibble(var1 = rownames(cor)[relevant[, "row"]],
                  var2 = colnames(cor)[relevant[, "col"]],
                  correlation = vals)
  final <- arrange(final, -correlation)
  return(final)
}
print(cor_table_fct(fra), n = 100)
