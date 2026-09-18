library(data.table)
library(dplyr)
library(ggplot2)
library(openxlsx)
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

load_fct <- function(cn) {
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
dt_long <- feat_long[dt_long, on = "chem"]

# Dumbell Plot
ggplot(dt_long, aes(sens, reorder(chem, rel_ab))) +
  geom_line() +
  geom_point(aes(color = model))

# Relative Frequency Plot
ggplot(dt_long[model %in% c("se_omni", "se_sing")], aes(p.type, sens, colour = model)) +
  geom_point() +
  geom_smooth(se = FALSE)

# Grouped Features Plot
feat_grps <- melt(feat_dt[, c(2:32, 35)], id.vars = "group", variable.name = "chem")
feat_grps <- feat_grps[complete.cases(feat_grps) & value == 1]

ggplot(feat_grps, aes(chem, fill = group)) +
  geom_bar() +
  coord_flip()



