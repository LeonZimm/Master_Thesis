# RPC - France River Monitoring
setwd("C:/Users/leonz/Documents/Ecotox_Studium/Master_Thesis/Data")

# required packages
library(tidyverse)
library(caret)
library(data.table)
library(sf)
library(ranger)
library(foreach)
library(doParallel)
library(pROC)
library(MLmetrics)
library(mltools)
library(precrec)
library(tmap)
library(cowplot)
library(GGally)
library(maptiles)
library(terra)
library(grid)


scale_stability_fct <- function(data, 
                                n_reps = 20,
                                cores = 7) {
  results <- list()
  sites <- unique(fra$site.id) 
  
  for (i in 1:n_reps) {
    cat("Iteration", i, "/", n_reps,
        "| Time:", format(Sys.time(), "%H:%M:%S"),
        "\n")
    set.seed(i)
    boot_sites <- sample(sites, 
                         length(sites),
                         replace = TRUE)
    boot_data <- data[site.id %in% boot_sites]
    boot_data <- boot_data[, setdiff(names(boot_data), c(conc_col, phychem_col, "site.id", "sample.d", "si_mean")), with = FALSE]
    
    fit <- ranger(det ~ ., 
                  data = boot_data, 
                  num.trees = 120,
                  importance = "permutation",
                  probability = TRUE,
                  num.threads = cores,
                  seed = i)
    
    imp <- as.data.table(fit$variable.importance,
                         keep.rownames = TRUE)
    setorder(imp, -V2)
    imp[, ":="(rep = i,
               order = seq(nrow(imp)))]
    results[[i]] <- imp
  }
  rbindlist(results)
}
IR_eval_fct <- function(model, validation, test, n_cores) {
  
  # Prediction
  prob_val <- predict(model, validation)$predictions[, "D"]
  
  thresholds <- seq(0.01, 0.99, 0.02)
  
  # Set up parallelization
  n_cores <- n_cores
  clus <- makeCluster(n_cores)
  on.exit(stopCluster(clus), add = TRUE)
  registerDoParallel(clus)
  
  thre_res <- foreach(t = thresholds,
                      .combine = rbind,
                      .packages = c("MLmetrics", "mltools", "data.table")) %dopar% {
                        
                        pred <- ifelse(prob_val >= t, "D", "ND")
                        
                        # MCC
                        mcc <- mltools::mcc(as.factor(pred),
                                            validation$det)
                        
                        data.table(threshold = t, mcc = mcc)
                      }
  
  best_t <- thre_res[which.max(mcc), threshold]
  
  prob_test <- predict(model, test)$predictions[, "D"]
  pred <- ifelse(prob_test >= best_t, "D", "ND")
  pred_num <- ifelse(pred == "D", 1, 0)
  test_num <- ifelse(test$det == "D", 1, 0)
  # AUC
  mm <- mmdata(prob_test,
               test_num,
               posclass = 1)
  evalmod <- evalmod(mm)
  
  auc <- auc(evalmod)[1, 4]
  pr_auc <- auc(evalmod)[2, 4]
  # F1
  f1 <- F1_Score(y_pred = pred,
                 y_true = test$det,
                 positive = "D")
  # MCC
  mcc <- mltools::mcc(as.factor(pred),
                      test$det)
  
  list("auc" = auc,
       "pr_auc" = pr_auc, 
       "f1" = f1,
       "mcc" = mcc,
       "threshold" = best_t,
       "pred" = pred,
       "evalmod" = as.data.table(evalmod))
}
evaluation_fct <- function(model_type, model, test) {
  # Prediction
  predi <- predict(model, test)$predictions
  
  if (model_type == "CONCENTRATION") {
    # MSE
    mse <- mean((test$conc - predi)^2)
    # RMSE
    rmse <- sqrt(mse)
    # MAE
    mae <- mean(abs(test$conc - predi))
    # R squared
    r2 <- 1 - sum((test$conc - predi)^2) / sum((test$conc - mean(test$conc))^2)
    
    list("mse" = mse,
         "rmse" = rmse,
         "mae" = mae, 
         "r2" = r2,
         "pred" = predi)
  } else {
  # AUC
  mm <- mmdata(predi[, 1],
               test$det,
               posclass = "D")
  evalmod <- evalmod(mm)
  auc <- auc(evalmod)[1, 4]
  pr_auc <- auc(evalmod)[2, 4]
  # F1
  pred <- as.factor(ifelse(predi[, 1] >= 0.5, "D", "ND"))
  f1 <- F1_Score(y_pred = pred,
                 y_true = test$det,
                 positive = "D")
  # MCC
  mcc <- mltools::mcc(as.factor(pred),
                      test$det)
  
  list("auc" = auc,
       "pr_auc" = pr_auc,
       "f1" = f1, 
       "mcc" = mcc,
       "pred" = pred,
       "evalmod" = as.data.table(evalmod))
  }
}
hp_tun_fct <- function(data,
                       param_name,
                       param_values,
                       model_type,
                       fixed_params = list(),
                       kfolds,
                       cores,
                       seed = 25) {
  
  # Table for results
  res_log <- data.table()
  # Create folds
  if (model_type == "CONCENTRATION") {
    set.seed(seed)
    folds <- createFolds(data$conc, 
                         k = kfolds,
                         list = TRUE,
                         returnTrain = FALSE)
    
  } else {
    set.seed(seed)
    folds <- createFolds(data$det, 
                         k = kfolds,
                         list = TRUE,
                         returnTrain = FALSE)
  }
  for (val in param_values) {
    cat(model_type, "\n",
        "Tuning", param_name, "=", val,
        "| Time:", format(Sys.time(), "%H:%M:%S"),
        "\n")
    
    
    for (k in 1:kfolds) {
      start <- proc.time()[[3]]
      cat(" Fold", k,
          " | Time:", format(Sys.time(), "%H:%M:%S"),
          "\n")
      
      # Prepare data
      idx <- setdiff(seq_len(nrow(data)), folds[[k]])
      train <- data[idx]
      test <- data[folds[[k]]]
      
      # Build dynamic parameter list
      if (model_type == "CONCENTRATION") {
        param_list <- list(formula = conc ~ .,
                           data = train,
                           respect.unordered.factors = "order",
                           num.threads = cores,
                           seed = seed)
      } else {
      param_list <- list(formula = det ~ .,
                         data = train,
                         respect.unordered.factors = "order",
                         probability = TRUE,
                         num.threads = cores,
                         seed = seed)
      }
      # Add tuning parameter
      param_list[[param_name]] <- val
      
      # Add fixed parameters
      param_list <- c(param_list, fixed_params)
      
      # Train model
      model <- do.call(ranger, param_list)
      
      # Predict, calibrate probabilities (isotonic regression), compute AUC/PR AUC and threshold 
      res_list <- evaluation_fct(model_type = model_type, 
                                 model, 
                                 test)
      
      end <- proc.time()[[3]]
      # Error log
      if (model_type == "CONCENTRATION") {
        res_log <- rbind(res_log, data.table(param_value = val,
                                             mse = res_list$mse,
                                             rmse = res_list$rmse,
                                             mae = res_list$mae,
                                             r2 =  res_list$r2,
                                             time = (end - start)/60))
      } else {
        res_log <- rbind(res_log, data.table(param_value = val,
                                             auc = res_list$auc,
                                             pr_auc = res_list$pr_auc,
                                             f1 = res_list$f1,
                                             mcc =  res_list$mcc,
                                             time = (end - start)/60))
      }
    }
    
    
  }
  
  setnames(res_log, "param_value", param_name)
  if (model_type == "CONCENTRATION") {
    res <- res_log[, .(mse = mean(mse),
                       rmse = mean(rmse),
                       mae = mean(mae),
                       r2 = mean(r2),
                       time = mean(time)), 
                   by = param_name]
  } else {
    res <- res_log[, .(auc = mean(auc),
                       pr_auc = mean(pr_auc),
                       f1 = mean(f1),
                       mcc = mean(mcc),
                       time = mean(time)), 
                   by = param_name]
    }  
    
  
  if (class(param_values) == "numeric") {
    if (param_name == "min.node.size") {
      tunval_plot <- ndsz_num_hypa_fct(model_type = model_type,
                                       res,
                                       param_name)
    } else{
      tunval_plot <- num_hypa_fct(model_type = model_type,
                                  res,
                                  param_name)  
    }
  } else {
    tunval_plot <- char_hypa_fct(model_type = model_type,
                                 res, 
                                 param_name)
  }
  
  output <- list("res_log" = res_log[],
                 "res" = res, 
                 "tunval" = tunval_plot[[1]],
                 "plot" = tunval_plot[[2]])
  
  return(output)
}
num_hypa_fct <- function(model_type, res, para){
  if (model_type == "CONCENTRATION") {
    best <- min(res$mse)
    best_zone <- best + 0.01 * best
    
    if (!any((res$mse < (best_zone)) & (res[[para]] < res[[para]][which.min(res$mse)]))) {
      para_tun = res[[para]][which.min(res$mse)]
    } else {
      para_tun = res[[para]][res[[para]] == min(res[[para]][(res$mse < (best_zone)) & (res[[para]] < res[[para]][which.max(res$mse)])])]
    }
    
    plot <- ggplot(res, aes(x = !!sym(para), y = mse, group = 1)) +
      geom_point() +
      geom_line() +
      geom_hline(yintercept = best_zone, color = "#cc0000")
    
    output <- list("para_tun" = para_tun, "plot" = plot)
    return(output)
    
  } else {
    best <- max(res$mcc)
    best_zone <- best - 0.01 * best
    
    if (!any((res$mcc > (best_zone)) & (res[[para]] < res[[para]][which.max(res$mcc)]))) {
      para_tun = res[[para]][which.max(res$mcc)]
    } else {
      para_tun = res[[para]][res[[para]] == min(res[[para]][(res$mcc > (best_zone)) & (res[[para]] < res[[para]][which.max(res$mcc)])])]
    }
    
    plot <- ggplot(res, aes(x = !!sym(para), y = mcc, group = 1)) +
      geom_point() +
      geom_line() +
      geom_hline(yintercept = best_zone, color = "#cc0000")
    
    output <- list("para_tun" = para_tun, "plot" = plot)
    return(output)  
  }
  
}
ndsz_num_hypa_fct <- function(model_type, res, para){
  if (model_type == "CONCENTRATION") {
    best <- min(res$mse)
    best_zone <- best + 0.01 * best
    
    if (!any((res$mse < (best_zone)) & (res[[para]] < res[[para]][which.min(res$mse)]))) {
      para_tun = res[[para]][which.min(res$mse)]
    } else {
      para_tun = res[[para]][res[[para]] == max(res[[para]][(res$mse < (best_zone)) & (res[[para]] > res[[para]][which.min(res$mse)])])]
    }
    
    plot <- ggplot(res, aes(x = !!sym(para), y = mse, group = 1)) +
      geom_point() +
      geom_line() +
      geom_hline(yintercept = best_zone, color = "#cc0000")
    
    output <- list("para_tun" = para_tun, "plot" = plot)
    return(output)
  } else {
    best <- max(res$mcc)
    best_zone <- best - 0.01 * best
    
    if (!any((res$mcc > (best_zone)) & (res[[para]] > res[[para]][which.max(res$mcc)]))) {
      para_tun = res[[para]][which.max(res$mcc)]
    } else {
      para_tun = res[[para]][res[[para]] == max(res[[para]][(res$mcc > (best_zone)) & (res[[para]] > res[[para]][which.max(res$mcc)])])]
    }
    
    plot <- ggplot(res, aes(x = !!sym(para), y = mcc, group = 1)) +
      geom_point() +
      geom_line() +
      geom_hline(yintercept = best_zone, color = "#cc0000")
    
    output <- list("para_tun" = para_tun, "plot" = plot)
    return(output)
  }
}
char_hypa_fct <- function(model_type, res, para) {
  if (model_type == "CONCENTRATION") {
    best <- min(res$mse)
    best_zone <- best + 0.01 * best
    
    if (!any((res$mse < (best_zone)) & (res$time < res$time[which.min(res$mse)]))) {
      para_tun = res[[para]][which.min(res$mse)]
    } else {
      para_tun = res[[para]][res[[para]] == min(res[[para]][(res$mse < (best_zone)) & (res$time < res$time[which.min(res$mse)])])]
    }
    
    plot <- ggplot(res, aes(x = !!sym(para), y = mcc, group = 1)) +
      geom_point() +
      geom_line() +
      geom_hline(yintercept = best_zone, color = "#cc0000")
    
    output <- list("para_tun" = para_tun, "plot" = plot)
    return(output)
  } else {
  best <- max(res$mcc)
  best_zone <- best - 0.01 * best
  
  if (!any((res$mcc > (best_zone)) & (res$time < res$time[which.max(res$mcc)]))) {
    para_tun = res[[para]][which.max(res$mcc)]
  } else {
    para_tun = res[[para]][res[[para]] == max(res[[para]][(res$mcc > (best_zone)) & (res$time < res$time[which.max(res$mcc)])])]
  }
  
  plot <- ggplot(res, aes(x = !!sym(para), y = mcc, group = 1)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = best_zone, color = "#cc0000")
  
  output <- list("para_tun" = para_tun, "plot" = plot)
  return(output)
  }
}
feat_rm_steps_fct <- function(n_features) {
  if (n_features > 40) {
    ceiling(0.10 * n_features)   
  } else if (n_features > 20) {
    ceiling(0.05 * n_features)  
  } else {
    1                             
  }
}
rf_fct <- function(data, 
                   cas_nr = NULL, 
                   model_type,
                   n_trees = 120,
                   kfolds = 5,
                   cores = 7,
                   ntree_val = c(50, 100, 150, 200)) {
  
  if (model_type == "SINGLE SUBSTANCE") {
    chem_name <- chem[cas == cas_nr]$chemical
    file_path <- file.path("results/sing_sub", chem_name)
    folder <- dir.create(file_path)
    data <- data[cas == cas_nr]
  } else if (model_type == "OMNIBUS") {
    file_path <- file.path("results/omni")
  } else if (model_type == "PHYSICO CHEMICAL") {
    file_path <- file.path("results/phychem")
  } else {
    file_path <- file.path("results/conc")
  }
  
  if (model_type == "CONCENTRATION") {
    set.seed(25)
    partition <- createDataPartition(data$conc,
                                     p = 0.7,
                                     list = FALSE,
                                     times = 1)
    data_train <- data[partition]
    data_test <- data[-partition]
  } else {
    set.seed(25)
    partition_subset <- createDataPartition(data$det,
                                            p = 0.7,
                                            list = FALSE,
                                            times = 1)
    sub <- data[partition_subset]
    data_test <- data[-partition_subset]
    
    set.seed(25)
    partition <- createDataPartition(sub$det, 
                                     p = 0.8,
                                     list = FALSE,
                                     times = 1)
    data_train <- sub[partition]
    data_val <- sub[-partition]
  }
  
  # 1. Feature Selection -------------------------------------------------------
  
  t_max <- ncol(data_train) - 1
  eval_log_path <- paste0(file_path, "/eval_log.rds")
  feat_log_path <- paste0(file_path, "/feat_log.rds")
  
  # Store results for analysis
  if (file.exists(eval_log_path) && file.exists(feat_log_path)) {
    
    rfs_eval_log <- readRDS(eval_log_path)
    rfs_feat_log <- readRDS(feat_log_path)
    beg <- max(rfs_feat_log$iter) + 1
    rfs_data <- data_train[, setdiff(names(data_train), rfs_feat_log$feature_rm), with = FALSE]
    
  } else {
    
    rfs_eval_log <- data.table()
    rfs_feat_log <- data.table()
    beg <- 1
    rfs_data <- copy(data_train)
  }
  
  
  
  if (model_type == "CONCENTRATION") {
    for (t in beg:t_max) {
      cat(model_type, "\n",
          "Variables left", ncol(rfs_data), 
          "| Iteration", t,
          "| Time:", format(Sys.time(), "%H:%M:%S"),
          "| Nrow", nrow(rfs_data),
          "\n")
      
      # Store k-fold variable importance measurements
      impu_list <- list()
      perm_list <- list()
      
      # Create k-folds 
      set.seed(25)
      folds_rfs <- createFolds(rfs_data$conc, 
                               k = kfolds, 
                               list = TRUE, 
                               returnTrain = FALSE)
      
      start <- proc.time()[[3]]
      for (k in 1:kfolds) {
        cat(" Fold", k,
            " | Time:", format(Sys.time(), "%H:%M:%S"),
            "\n")
        
        # Prepare data
        rfs_index <- folds_rfs[[k]]
        rfs_train <- rfs_data[-rfs_index]
        rfs_test <- rfs_data[rfs_index]
        
        # RF for impurity importance
        rfs_impu <- ranger(conc ~ .,
                           data = rfs_train,
                           num.trees = n_trees,
                           importance = "impurity",
                           respect.unordered.factors = "order",
                           num.threads = cores,
                           seed = 25)
        
        # RF for permutation importance
        rfs_perm <- ranger(conc ~ .,
                           data = rfs_train,
                           num.trees = n_trees,
                           importance = "permutation",
                           respect.unordered.factors = "order",
                           num.threads = cores,
                           seed = 25)
        
        # Predict and Calculate Performance Metrics 
        res_list <- evaluation_fct(model_type = model_type,
                                   rfs_perm,
                                   rfs_test)
        
        # Variable importance table + error log
        impu_list[[k]] <- as.data.table(rfs_impu$variable.importance, keep.rownames = TRUE)
        perm_list[[k]] <- as.data.table(rfs_perm$variable.importance, keep.rownames = TRUE)
        
        rfs_eval_log <- rbind(rfs_eval_log, data.table(iter = t,
                                                       mse = res_list$mse,
                                                       rmse = res_list$rmse,
                                                       mae = res_list$mae,
                                                       r2 = res_list$r2))
      }
      end <- proc.time()[[3]]
      time_rfs <- end - start                           
      
      saveRDS(rfs_eval_log, paste0(file_path, "/eval_log.rds"))
      
      # Calculate mean of cross validated importance measurements
      impu_dt <- rbindlist(impu_list)
      perm_dt <- rbindlist(perm_list)
      
      impu_mean <- impu_dt[, .(impu = mean(V2)), by = V1]
      perm_mean <- perm_dt[, .(perm = mean(V2)), by = V1]
      
      va_impo_dt <- merge(impu_mean, perm_mean, by = "V1") 
      
      # Normalizing importance measurements and scoring features
      z_score <- function(x) (x - mean(x))/sd(x)
      
      va_impo_dt[, impu_z := z_score(impu)]
      va_impo_dt[, perm_z := z_score(perm)]
      va_impo_dt[, score := 0.5*impu_z + 0.5*perm_z] 
      va_impo_dt <- va_impo_dt[order(score)]
      
      # Removing 5% of the worst features
      feat_rem <- va_impo_dt$V1[1:feat_rm_steps_fct(nrow(va_impo_dt))]
      rfs_data <- rfs_data[, (feat_rem) := NULL]
      
      # Storing results
      rfs_feat_log <- rbind(rfs_feat_log, data.table(iter = t,
                                                     n_features = nrow(va_impo_dt),
                                                     feature_rm = feat_rem,
                                                     time = time_rfs))
      saveRDS(rfs_feat_log, paste0(file_path, "/feat_log.rds"))
      
      if (ncol(rfs_data) <= 1) {
        break
      }
      
      print(t)
    }
    
    rfs_res <- rfs_eval_log[, .(mse = mean(mse),
                                rmse = mean(rmse),
                                mae = mean(mae),
                                r2 = mean(r2)),
                            by = iter]
    
    rfs_res <- merge(rfs_res, rfs_feat_log, 
                     by = "iter", 
                     all.x = TRUE)
    
    rfs_plot <- ggplot(rfs_res, aes(as.factor(n_features), mse)) +
      geom_point() +
      theme_bw() +
      geom_hline(yintercept = min(rfs_res$mse) + 0.01 * min(rfs_res$mse), color = "#cc0000")
    
    # Identify the columns to remove
    cols_to_remove <- rfs_res$feature_rm[1:(which.min(rfs_res$mse) - 1)]
    # Get remaining column names
    rfs_features <- setdiff(colnames(data), cols_to_remove)
    
    rfs <- list("eval_log" = rfs_eval_log,
                "feat_log" = rfs_feat_log,
                "res" = rfs_res,
                "plot" = rfs_plot,
                "remaining_features" = rfs_features)
    saveRDS(rfs, file.path(file_path, "rfs.rds"))
    
  } else {
    for (t in beg:t_max) {
      cat(model_type, "\n",
          "Variables left", ncol(rfs_data), 
          "| Iteration", t,
          "| Time:", format(Sys.time(), "%H:%M:%S"),
          "| Nrow", nrow(rfs_data),
          "\n")
      
      # Store k-fold variable importance measurements
      impu_list <- list()
      perm_list <- list()
      
      # Create k-folds 
      set.seed(25)
      folds_rfs <- createFolds(rfs_data$det, 
                               k = kfolds, 
                               list = TRUE, 
                               returnTrain = FALSE)
      
      start <- proc.time()[[3]]
      for (k in 1:kfolds) {
        cat(" Fold", k,
            " | Time:", format(Sys.time(), "%H:%M:%S"),
            "\n")
        
        # Prepare data
        rfs_index <- folds_rfs[[k]]
        rfs_train <- rfs_data[-rfs_index]
        rfs_test <- rfs_data[rfs_index]
        
        # RF for impurity importance
        rfs_impu <- ranger(det ~ .,
                           data = rfs_train,
                           num.trees = n_trees,
                           importance = "impurity",
                           respect.unordered.factors = "order",
                           num.threads = cores,
                           seed = 25)
        
        # RF for permutation importance
        rfs_perm <- ranger(det ~ .,
                           data = rfs_train,
                           num.trees = n_trees,
                           importance = "permutation",
                           respect.unordered.factors = "order",
                           num.threads = cores,
                           probability = TRUE,
                           seed = 25)
        
        # Predict and Calculate Performance Metrics 
        res_list <- evaluation_fct(model_type = model_type,
                                   rfs_perm, 
                                   rfs_test)
        
        # Variable importance table + error log
        impu_list[[k]] <- as.data.table(rfs_impu$variable.importance, keep.rownames = TRUE)
        perm_list[[k]] <- as.data.table(rfs_perm$variable.importance, keep.rownames = TRUE)
        
        rfs_eval_log <- rbind(rfs_eval_log, data.table(iter = t,
                                                       auc = res_list$auc,
                                                       pr_auc = res_list$pr_auc,
                                                       f1 = res_list$f1,
                                                       mcc = res_list$mcc))
      }
      end <- proc.time()[[3]]
      time_rfs <- end - start                           
      
      saveRDS(rfs_eval_log, paste0(file_path, "/eval_log.rds"))
      
      # Calculate mean of cross validated importance measurements
      impu_dt <- rbindlist(impu_list)
      perm_dt <- rbindlist(perm_list)
      
      impu_mean <- impu_dt[, .(impu = mean(V2)), by = V1]
      perm_mean <- perm_dt[, .(perm = mean(V2)), by = V1]
      
      va_impo_dt <- merge(impu_mean, perm_mean, by = "V1") 
      
      # Normalizing importance measurements and scoring features
      z_score <- function(x) (x - mean(x))/sd(x)
      
      va_impo_dt[, impu_z := z_score(impu)]
      va_impo_dt[, perm_z := z_score(perm)]
      va_impo_dt[, score := 0.5*impu_z + 0.5*perm_z] 
      va_impo_dt <- va_impo_dt[order(score)]
      
      # Removing 5% of the worst features
      feat_rem <- va_impo_dt$V1[1:feat_rm_steps_fct(nrow(va_impo_dt))]
      rfs_data <- rfs_data[, (feat_rem) := NULL]
      
      # Storing results
      rfs_feat_log <- rbind(rfs_feat_log, data.table(iter = t,
                                                     n_features = nrow(va_impo_dt),
                                                     feature_rm = feat_rem,
                                                     time = time_rfs))
      saveRDS(rfs_feat_log, paste0(file_path, "/feat_log.rds"))
      
      if (ncol(rfs_data) <= 1) {
        break
      }
      
      print(t)
    }
    
    rfs_res <- rfs_eval_log[, .(auc = mean(auc),
                                pr_auc = mean(pr_auc),
                                f1 = mean(f1),
                                mcc = mean(mcc)),
                            by = iter]
    
    rfs_res <- merge(rfs_res, rfs_feat_log, by = "iter", all.x = TRUE)
    
    rfs_plot <- ggplot(rfs_res, aes(as.factor(n_features), mcc)) +
      geom_point() +
      theme_bw() +
      geom_hline(yintercept = max(rfs_res$mcc) - 0.01 * max(rfs_res$mcc), color = "#cc0000")
    
    # Identify the columns to remove
    cols_to_remove <- rfs_res$feature_rm[1:(which.max(rfs_res$mcc) - 1)]
    # Get remaining column names
    rfs_features <- setdiff(colnames(data), cols_to_remove)
    
    rfs <- list("eval_log" = rfs_eval_log,
                "feat_log" = rfs_feat_log,
                "res" = rfs_res,
                "plot" = rfs_plot,
                "remaining_features" = rfs_features)
    saveRDS(rfs, file.path(file_path, "rfs.rds"))
  }
  
  # 2. Hyperparameter Optimization using Sequential Method ------------------------------------------------------
  hp_tun_data <- copy(data_train[, ..rfs_features])
  
  
  # 2.1 Number of trees -----------------------------------------------------
  # ntree_val <- c(50, 100, 120, 150, 200, 500)
  
  ntree <- hp_tun_fct(hp_tun_data,
                      "num.trees",
                      ntree_val,
                      fixed_params = list(),
                      kfolds = kfolds,
                      cores = cores,
                      model_type = model_type)
  
  
  # 2.2 Mtry ----------------------------------------------------------------
  mtry_val <- seq(2, ncol(hp_tun_data) - 1, by = 2)
  
  mtry <- hp_tun_fct(hp_tun_data,
                     "mtry",
                     mtry_val,
                     fixed_params = list(num.trees = ntree$tunval),
                     kfolds = kfolds,
                     cores = cores,
                     model_type = model_type)
  
  
  # 2.3 Sample fraction ---------------------------------------------------------
  safr_val <- seq(1, 0.5, by = -0.1)
  
  safr <- hp_tun_fct(hp_tun_data,
                     "sample.fraction",
                     safr_val,
                     fixed_params = list(num.trees = ntree$tunval,
                                         mtry = mtry$tunval),
                     kfolds = kfolds,
                     cores = cores,
                     model_type = model_type)
  
  
  # 2.4 Minimal node size ---------------------------------------------------
  ndsz_val <- seq(1, 6, 1)
  
  ndsz <- hp_tun_fct(hp_tun_data,
                     "min.node.size",
                     ndsz_val,
                     fixed_params = list(num.trees = ntree$tunval,
                                         mtry = mtry$tunval,
                                         sample.fraction = safr$tunval),
                     kfolds = kfolds,
                     cores = cores,
                     model_type = model_type)
  
  
  # 2.5 Splitting rule ------------------------------------------------------
  if (model_type == "CONCENTRATION") {
    sprl_val <- c("variance", "extratrees", "poisson")
  } else {
    sprl_val <- c("gini", "extratrees", "hellinger")
  }
  
  sprl <- hp_tun_fct(hp_tun_data,
                     "splitrule",
                     sprl_val,
                     fixed_params = list(num.trees = ntree$tunval,
                                         mtry = mtry$tunval,
                                         sample.fraction = safr$tunval,
                                         min.node.size = ndsz$tunval),
                     kfolds = kfolds,
                     cores = cores,
                     model_type = model_type)
  
  hypa <- list("ntree" = ntree,
               "mtry" = mtry,
               "safr" = safr,
               "ndsz" = ndsz,
               "sprl" = sprl)
  saveRDS(hypa, file.path(file_path, "hypa.rds"))
  
  
  # 3. RF-Models --------------------------------------------------
  if (model_type == "CONCENTRATION") {
    basli_model_ranger <- ranger(conc ~ .,
                                 data = data_train,
                                 num.trees = 120,
                                 respect.unordered.factors = "order",
                                 importance = "permutation",
                                 num.threads = cores,
                                 seed = 25)
    
    basli_res_list <- evaluation_fct(model_type = model_type, 
                                     basli_model_ranger, 
                                     data_test)
    
    basli_model <- list("res" = basli_res_list,
                        "model" = basli_model_ranger)
    
    opt_model_ranger <- ranger(conc ~ .,
                               data = hp_tun_data,
                               num.trees = ntree$tunval,
                               mtry = mtry$tunval,
                               sample.fraction = safr$tunval,
                               min.node.size = ndsz$tunval,
                               splitrule = sprl$tunval,
                               importance = "permutation",
                               respect.unordered.factors = "order",
                               num.threads = cores,
                               seed = 25)
    
    opt_res_list <- evaluation_fct(model_type = model_type,
                                   opt_model_ranger,
                                   data_test)
    
    opt_model <- list("res" = opt_res_list,
                      "model" = opt_model_ranger,
                      "features" = rfs_features)
    
    model_list <- list("basli" = basli_model,
                       "opt" = opt_model)
    saveRDS(model_list, file.path(file_path, "models.rds"))
  } else {
    basli_model_ranger <- ranger(det ~ .,
                                 data = data_train,
                                 num.trees = 120,
                                 respect.unordered.factors = "order",
                                 importance = "permutation",
                                 num.threads = cores,
                                 probability = TRUE,
                                 seed = 25)
    
    basli_res_list <- evaluation_fct(model_type = model_type, 
                                     basli_model_ranger,
                                     data_test)
    cf_basli <- confusionMatrix(as.factor(basli_res_list$pred),
                                data_test$det, 
                                positive = "D")
    
    basli_model <- list("res" = basli_res_list,
                        "confu_mat" = cf_basli,
                        "model" = basli_model_ranger)
    
    opt_model_ranger <- ranger(det ~ .,
                               data = hp_tun_data,
                               num.trees = ntree$tunval,
                               mtry = mtry$tunval,
                               sample.fraction = safr$tunval,
                               min.node.size = ndsz$tunval,
                               splitrule = sprl$tunval,
                               importance = "permutation",
                               probability = TRUE,
                               respect.unordered.factors = "order",
                               num.threads = cores,
                               seed = 25)
    
    opt_res_list <- IR_eval_fct(opt_model_ranger,
                                data_val, 
                                data_test, 
                                10)
    cf_opt <- confusionMatrix(as.factor(opt_res_list$pred),
                              data_test$det, 
                              positive = "D")
    
    opt_model <- list("res" = opt_res_list,
                      "confu_mat" = cf_opt,
                      "model" = opt_model_ranger,
                      "features" = rfs_features)
    
    model_list <- list("basli" = basli_model,
                       "opt" = opt_model)
    saveRDS(model_list, file.path(file_path, "models.rds"))
  }
}
model_results_fct <- function(model) {
  file_path <- list.files(path = paste0("results/", model),
                          pattern = "\\.rds$",
                          full.names = TRUE)
  file_names <- gsub(pattern = "\\.rds$", 
                     replacement = "", 
                     x = basename(file_path))
  results <- lapply(file_path, 
                    readRDS)
  names(results) <- file_names
  return(results)
}


# Data Preparation ----------------------------------------------
fra <- readRDS("Tables_prepared/fra.rds")

# Columns to be Removed
# Concentration
conc_col <- c("conc", "ld", "lq")
# Physicochemcial
phychem_col <- c("sw", "kow", "pka", "p", "h", "dr_soil", "dr_wat_sed", "koc", "gus")

# Compare Catchment Scale Predictors 
# scale_res <- scale_stability_fct(fra)
# scale_res <- scale_res[, ":="(avg_order = mean(order), sd_order = sd(order)), by = V1]
# scale_res <- scale_res[, avg_order := mean(order), by = V1]
# saveRDS(scale_res, "results/catch_scale.rds")
# 10km catchment continuously outperforms 1km/100m predictors -> remove 1km/100m predictors
pred_names <- names(fra)[names(fra) != "site.id"]
scale_rm <- grepl("_100m|_1km", pred_names)
scale_rm <- pred_names[scale_rm]
rm_col <- c("site.id", "sample.d", "si_mean", scale_rm) # si_mean has cor>0.8 with sa_mean and k_mean


# Omnibus Model -----------------------------------------------------------
try_cas <- fra_cas[, .N, by = .(cas, det, p.type)][, rel_ab := N/sum(N), by = cas]
try_cas <- try_cas[, .SD[which.max(rel_ab)], by = p.type]

cols_cas <-  setdiff(names(fra), c(rm_col, conc_col, phychem_col))

rf_fct(fra[complete.cases(fra[, ..cols_cas]), ..cols_cas],
       model_type = "OMNIBUS",
       cores = 5)
omnibus <- model_results_fct("omni")
saveRDS(omnibus, "results/omnibus.rds")

# Single-Substance Model --------------------------------------------------
chem <- readRDS("Raw/chemical.processed.rds")
substances <- unique(fra$cas)

for (cas_nr in substances) {
  chem_name <- chem[cas == cas_nr]$chemical
  cat(chem_name, "\n")
  
  start <- proc.time()[[3]]
  
  rf_fct(fra[complete.cases(fra[, ..cols_cas]), ..cols_cas], 
         cas_nr = cas_nr,
         model_type = "SINGLE SUBSTANCE",
         cores = 5)
  
  end <- proc.time()[[3]]
  time <- (end - start)/60/60
  
  cat(chem_name, 
      "\nTime taken:", round(time, 2), "hours",
      "\n\n")
}

chem_name <- chem[cas %in% substances]$chemical
all_mod <- setNames(vector("list", length(chem_name)), chem_name)
for (cn in chem_name) {
  results <- model_results_fct(paste0("sing_sub/", cn))
  all_mod[[cn]] <- results
}
saveRDS(all_mod, "results/sing_sub.rds")


# Physicochemical Model -----------------------------------------------------------
cols_phychem <- setdiff(names(fra), c(rm_col, conc_col, "cas"))
rf_fct(fra[complete.cases(fra[, ..cols_phychem]), ..cols_phychem],
       model_type = "PHYSICO CHEMICAL",
       cores = 5)

phychem <- model_results_fct("phychem")
saveRDS(phychem, "results/phychem.rds")










# Concentration Model -----------------------------------------------------
cols_conc <-  setdiff(names(fra), c(rm_col, phychem_col, "det"))
conc_dt <- fra[det == "D"]
conc_dt <- conc_dt[complete.cases(conc_dt[, ..cols_conc]), ..cols_conc]

rf_fct_conc(conc_dt,
            model_type = "CONCENTRATION",
            cores = 7)


