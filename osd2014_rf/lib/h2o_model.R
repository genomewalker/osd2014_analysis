library(tidyverse)

ohi_median2class <- function(X){

  case_when(X < 1.4 ~ "Very low impact", #64DD17
            X >= 1.4 & X < 4.95 ~ "Low impact", #AEEA00
            X >= 4.95 & X < 8.47 ~ "Medium impact", #FFD600
            X >= 8.47 & X < 12 ~ "Medium high impact", #FFAB00
            X >= 12 & X < 15.52 ~ "High impact", #FF6D00
            X >= 15.52 ~ "Very high impact") #DD2C00

}

ohi_median2class_short <- function(X){

  case_when(X < 4.95 ~ "Low impact", #AEEA00
            X >= 4.95 ~ "Impacted") #DD2C00

}


get_error_ensemble <- function(X){
  pred <- predict(X$ensemble, newdata = as.h2o(X$test))
  X$test %>%
    as_tibble() %>% select(median) %>% #mutate(median = median^2) %>%
    mutate(pred = pred[,1] %>% as.vector())  %>%
    dplyr::rename(actual = median) %>%
    mutate(
      error     = actual - pred,
      error_pct = error / actual,
      class_actual = ohi_median2class_short(actual),
      class_pred = ohi_median2class_short(pred),
      classif = ifelse(class_actual == class_pred, "correct", "wrong"),
      class_actual = as.factor(class_actual),
      class_pred = as.factor(class_pred),
      class_actual = fct_expand(class_actual, levels(class_pred)),
      class_pred = fct_expand(class_pred, levels(class_actual)),
      class_actual = fct_relevel(class_actual, ohi_thrs),
      class_pred = fct_relevel(class_pred, ohi_thrs)
    )
}

get_error <- function(X, data = data){
  rep <- X
  data <- data[[X]]
  pred <- predict(data$aml@leader, newdata = as.h2o(data$test))
  data$test %>%
    as_tibble() %>% dplyr::select(median) %>% #mutate(median = median^2) %>%
    mutate(pred = pred[,1] %>% as.vector())  %>%
    dplyr::rename(actual = median) %>%
    mutate(
      error     = actual - pred,
      error_pct = error / actual,
      class_actual = ohi_median2class_short(actual),
      class_pred = ohi_median2class_short(pred),
      classif = ifelse(class_actual == class_pred, "correct", "wrong"),
      class_actual = as.factor(class_actual),
      class_pred = as.factor(class_pred),
      class_actual = fct_expand(class_actual, levels(class_pred)),
      class_pred = fct_expand(class_pred, levels(class_actual)),
      class_actual = fct_relevel(class_actual, ohi_thrs),
      class_pred = fct_relevel(class_pred, ohi_thrs),
      model_id = data$aml@leader@model_id,
      algo = data$aml@leader@algorithm,
      rep = rep
    )
}

summarise_errors <- function(X){
  X %>%
    summarise(
      me   = mean(error),
      rmse = mean(error^2)^0.5,
      mae  = mean(abs(error)),
      MAE = median(abs(error)),
      mape = mean(abs(error_pct)),
      smape = round(mean(2 * abs(actual - pred)/(abs(actual) + abs(pred)), na.rm = TRUE), digits = 3),
      mpe  = mean(error_pct)
    )
}

plot_error_point <- function(X){
  X %>%
    ggplot(aes(pred, actual)) +
    geom_abline(linetype = 2, alpha = 0.6) +
    geom_smooth(method=glm) +
    geom_point(color = "darkred", alpha = 0.5) +
    expand_limits(x = 0, y = 0) +
    ggtitle(paste0(unique(X$rep),": ", unique(X$algo)))
}


get_ensemble <- function(models){
  get_leaders <- function(X) {
    X@leader
  }

  base_models <- unique(lapply(models$aml, get_leaders))

  ensemble <- h2o.stackedEnsemble(x = setdiff(names(as.h2o(models$train)), "median"),
                                  y = "median",
                                  training_frame = as.h2o(models$train),
                                  base_models = base_models)


  var_importance <- map_df(base_models, h2o.varimp, .id = "base_model")
  var_importance_aggregated <- var_importance %>%
    group_by(variable) %>%
    summarise(relative_importance = mean(relative_importance)) %>%
    mutate(scaled_importance = relative_importance/max(relative_importance)) %>%
    arrange(desc(relative_importance)) %>%
    filter(relative_importance > 0)

  list(ensemble = ensemble, train = models$train, test = models$test, base_models = base_models, var_importance = var_importance, var_importance_aggregated = var_importance_aggregated)

}


save_model_data <- function(X, type = type, local_path = local_path, remote_path = remote_path, n_fold = n_fold){
  folder <- as.Date(as.POSIXct(Sys.time()))
  path <- file.path(remote_path, type, folder)
  path_local <- file.path(local_path, type, folder)
  if (dir.exists(path_local)){
    unlink(path_local)
    cat("Removing existing folder:", path_local)
  }else{
    dir.create(path_local, recursive = TRUE)
  }
  for (i in paste0("rep_",1:n_fold)){
    map(X[i], function(Y){
      h2o.saveModel(h2o.getModel(Y$aml), path = file.path(path, i), force = TRUE)
      h2o.saveModelDetails(h2o.getModel(Y@model_id), path = file.path(path, i), force = TRUE)
      h2o.download_mojo(h2o.getModel(Y@model_id), path = path_local, get_genmodel_jar = TRUE)
    })
  }
}

save_ensemble_model_data <- function(X, type = type, local_path = local_path, remote_path = remote_path){
  folder <- as.Date(as.POSIXct(Sys.time()))
  path <- file.path(remote_path, type, folder)
  path_local <- file.path(local_path, type, folder)
  if (dir.exists(path_local)){
    unlink(path_local)
    cat("Removing existing folder:", path_local)
  }else{
    dir.create(path_local, recursive = TRUE)
  }
  for (i in paste0("rep_",1:5)){
    map(ensembles[[i]]$base_models, function(X){
      h2o.saveModel(h2o.getModel(X@model_id), path = file.path(path, i), force = TRUE)
      h2o.saveModelDetails(h2o.getModel(X@model_id), path = file.path(path, i), force = TRUE)
      h2o.download_mojo(h2o.getModel(X@model_id), path = path_local, get_genmodel_jar = TRUE)
    })
    h2o.saveModel(h2o.getModel(ensembles[[i]]$ensemble@model_id), path = file.path(path, i), force = TRUE)
    h2o.saveModel(h2o.getModel(ensembles[[i]]$ensemble@model_id), path = file.path(path, i), force = TRUE)
    h2o.download_mojo(h2o.getModel(ensembles[[i]]$ensemble@model_id), path = path_local, get_genmodel_jar = TRUE)
  }
}



load_model_data <- function(X, type = type, remote_path = remote_path){
  folder <- as.Date(as.POSIXct(Sys.time()))
  path <- file.path(remote_path, type, folder)
  path_local <- file.path(local_path, type, folder)
  if (dir.exists(path_local)){
    cat("Found folder:", path_local)
  }else{
    stop("Folder not found")
  }

  for (i in paste0("rep_",1:5)){
    map(ensembles[[i]]$base_models, function(X){
      h2o.saveModel(h2o.getModel(X@model_id), path = file.path(path, i), force = TRUE)
      h2o.saveModelDetails(h2o.getModel(X@model_id), path = file.path(path, i), force = TRUE)
      h2o.download_mojo(h2o.getModel(X@model_id), path = path_local, get_genmodel_jar = TRUE)
    })
    h2o.saveModel(h2o.getModel(ensembles[[i]]$ensemble@model_id), path = file.path(path, i), force = TRUE)
    h2o.saveModel(h2o.getModel(ensembles[[i]]$ensemble@model_id), path = file.path(path, i), force = TRUE)
    h2o.download_mojo(h2o.getModel(ensembles[[i]]$ensemble@model_id), path = path_local, get_genmodel_jar = TRUE)
  }
}



plot_residuals_point <- function(X){
  X %>%
    ggplot(aes(actual, error)) +
    geom_abline(linetype = 2, alpha = 0.6, slope = 0, intercept = 0) +
    geom_point(color = "darkred", alpha = 0.5) +
    ggtitle(paste0(unique(X$rep),": ", unique(X$algo)))
}



