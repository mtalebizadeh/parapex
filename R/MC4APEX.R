#' MC4APEX
#' @description Performs Monte Carlo simulations using parallel computation.
#' @param Input A list containing analysis setting.
#'
#' @return A GSA_Object that is used for sensitivity analysis.
#' @export
#'
#' @examples GSA <- MC4APEX(Input)
#'
MC4APEX <- function(Input) {
  sample_size <- Input$sample_size
  coreNumber <- Input$coreNumber
  GSA_Type = Input$GSA_Type
  caption_var_sim <- Input$caption_var_sim
  caption_var_obs <- Input$caption_var_obs
  start_date <- Input$start_date
  end_date <- Input$end_date
  label_APEX_exe <- Input$label_APEX_exe
  label_watershed_param <- Input$label_watershed_param
  label_control_param <- Input$label_control_param
  label_output_variable_AWP <- Input$label_output_variable_AWP
  label_output_variable_ACY <- Input$label_output_variable_ACY
  label_output_variable_DWS <- Input$label_output_variable_DWS
  label_prior_dist <- Input$label_prior_dist
  label_observed_var.txt <- Input$label_observed_var.txt
  folder_path_project <- Input$folder_path_project
  folder_path_R_codes <- Input$folder_path_R_codes
  folder_path_observed <- Input$folder_path_observed
  folder_path_GSA_Outputs <- Input$folder_path_GSA_Outputs
  setwd(folder_path_R_codes)
  Back_Up_PARM0806.dat <- Input$Back_Up_PARM0806.dat
  Back_Up_APEXCONT.dat <- Input$Back_Up_APEXCONT.dat
  Back_Up_SUB0027.SUB <- Input$Back_Up_SUB0027.SUB
  store_folder_path_watershed <- Input$store_folder_path_watershed
  store_folder_path_control <- Input$store_folder_path_control
  Calculated_output_folder_AWP <- Input$Calculated_output_folder_AWP
  Calculated_output_folder_ACY <- Input$Calculated_output_folder_ACY
  Calculated_output_folder_DWS <- Input$Calculated_output_folder_DWS
  all_basin_parameters <- Input$APEX_PARM
  all_control_parameters <- Input$APEX_CONT
  uncertain_basin_parameters <- param_extract(all_basin_parameters)
  uncertain_control_parameters <- param_extract(all_control_parameters)

  if (ncol(uncertain_basin_parameters) != 0
      && ncol(uncertain_control_parameters) != 0) {
    all_uncertain_model_parameters <- cbind(uncertain_basin_parameters,
                                            uncertain_control_parameters)

  } else if (ncol(uncertain_basin_parameters) !=0
             && ncol(uncertain_control_parameters) == 0) {
    all_uncertain_model_parameters <- uncertain_basin_parameters
  } else if (ncol(uncertain_basin_parameters) == 0
             && ncol(uncertain_control_parameters) != 0) {
    all_uncertain_model_parameters <- uncertain_control_parameters
  } else {
    stop("No uncertain parameter was found!")
  }
  uncertain_basin_parameters <- param_extract(all_basin_parameters)
  uncertain_control_parameters <- param_extract(all_control_parameters)
  all_uncertain_model_parameters <- cbind(uncertain_basin_parameters,
                                          uncertain_control_parameters)
  lower_limits <-as.vector(all_uncertain_model_parameters[1,],
                           mode = "numeric")
  upper_limits <-as.vector(all_uncertain_model_parameters[2,],
                           mode = "numeric")
  if (Input$GSA_Type == "MORRIS") {
    GSA_Object <- sensitivity::morris(
      model = NULL,
      factor = names(all_uncertain_model_parameters),
      r = as.numeric(Input$SA_Parms$morris_r),
      design = list(type = "oat",
                    levels = as.numeric(Input$SA_Parms$morris_levels),
                    grid.jump = 1),
     binf = lower_limits,
     bsup = upper_limits, scale = TRUE)
    GSA_Object$X <- as.data.frame(GSA_Object$X);
  } else {
    GSA_Object <- gsa_object(sample_size,
                             lower_limits,
                             upper_limits,
                             all_uncertain_model_parameters,
                             GSA_Type)
  }
  save.image(paste(folder_path_GSA_Outputs,
                   "/", GSA_Type, "_", "Xmat",
                   "_", format(Sys.time(), "%Y_%B_%d_%H_%M"),
                   ".RData", sep = ""))
  write.table(GSA_Object$X,
              paste(folder_path_GSA_Outputs,
                    "/", label_prior_dist,
                    sep = ""),
              sep = "\t")
  sample_size <- nrow(GSA_Object$X)
  X_list <- list(i = 1)
  X_list$sample_size <- sample_size
  X_list$folder_path_project <- folder_path_project
  X_list$label_watershed_param <- label_watershed_param
  X_list$GSA_Object <- GSA_Object
  X_list$uncertain_basin_parameters <- uncertain_basin_parameters
  X_list$Back_Up_PARM0806.dat <- Back_Up_PARM0806.dat
  X_list$label_control_param <- label_control_param
  X_list$uncertain_control_parameters <- uncertain_control_parameters
  X_list$Back_Up_APEXCONT.dat <- Back_Up_APEXCONT.dat
  X_list$label_APEX_exe <- label_APEX_exe
  X_list$folder_path_R_codes <- folder_path_R_codes
  X_list$store_folder_path_watershed <- store_folder_path_watershed
  X_list$store_folder_path_control <- store_folder_path_control
  X_list$label_output_variable_AWP <- label_output_variable_AWP
  X_list$Calculated_output_folder_AWP <- Calculated_output_folder_AWP
  X_list$label_output_variable_ACY <- label_output_variable_ACY
  X_list$Calculated_output_folder_ACY <- Calculated_output_folder_ACY
  X_list$label_output_variable_DWS <- label_output_variable_DWS
  X_list$Calculated_output_folder_DWS <- Calculated_output_folder_DWS
  parallel_runer(X_list, coreNumber)
  GSA_Object <- list(GSA_Object = GSA_Object, Input=Input)
}
#
# Helper functions
#####
parallel_runer <- function(X_list, Number_cores = parallel::detectCores() - 1) {
  # Setups inputs for parallel computation.

  XX_list <- seq_2_parallel(X_list, Number_cores)
  cl <- parallel::makeCluster(Number_cores, type = "PSOCK")
  parallel::parLapply(cl,XX_list, output_model_run_Parallel)
  parallel::stopCluster(cl)
  #Removing clone folders...
  i = 1
  folder_name_colones <- 1
  while (i <= Number_cores) {
    folder_name_colones <- paste(X_list$folder_path_project,
                                 "_colon", toString(i), sep = "")
    unlink(folder_name_colones, recursive = TRUE)
    i = i + 1
  }
}
#
#####
seq_2_parallel <- function(X_list, Number_cores) {
  # Prepares inputs for parallel computation.
  batch_size_average <- floor(X_list$sample_size/Number_cores)
  solid_length <- batch_size_average*Number_cores
  batch_start_idx <- seq(1, solid_length, by = batch_size_average)
  batch_end_idx <- batch_start_idx - 1
  batch_end_idx <- batch_end_idx[-1]
  batch_end_idx[length(batch_end_idx) + 1] <- X_list$sample_size
  XX_list <- list(X_list)

  j = 1
  while (j <= length(batch_start_idx)) {
    XX_list[[j]] <- X_list
    XX_list[[j]]["i"] <- batch_start_idx[j]
    XX_list[[j]]["sample_size"] <- batch_end_idx[j]
    j = j + 1
  }
  # Creating multiple project folders for parallel computation
  folder_name_colones <- X_list$folder_path_project
  file_path_from <- paste(X_list$folder_path_project,
                          "/", list.files(X_list$folder_path_project),
                          sep = "")
  i = 1
  while (i <= Number_cores) {
    folder_name_colones[i] <- paste(X_list$folder_path_project,
                                    "_colon",
                                    toString(i),
                                    sep="")
    file_path_to <- paste(X_list$folder_path_project,
                          "_", "colon", toString(i),
                          "/", list.files(X_list$folder_path_project),
                          sep = "")
    dir.create(folder_name_colones[i])
    file.copy(file_path_from,
              file_path_to,
              overwrite = TRUE)
    XX_list[[i]]["folder_path_project"] <- folder_name_colones[i]
    i = i + 1
  }
  return(XX_list)
}
#
#####
output_model_run_Parallel <- function(X_list) {
  # Used by 'parallel_runer' function.

  i <- X_list$i
  sample_size <- X_list$sample_size
  folder_path_project <- X_list$folder_path_project
  label_watershed_param <- X_list$label_watershed_param
  GSA_Object <- X_list$GSA_Object
  uncertain_basin_parameters <- X_list$uncertain_basin_parameters
  Back_Up_PARM0806.dat <- X_list$Back_Up_PARM0806.dat
  label_control_param <- X_list$label_control_param
  uncertain_control_parameters <- X_list$uncertain_control_parameters
  Back_Up_APEXCONT.dat <- X_list$Back_Up_APEXCONT.dat
  label_APEX_exe <- X_list$label_APEX_exe
  folder_path_R_codes <- X_list$folder_path_R_codes
  store_folder_path_watershed <- X_list$store_folder_path_watershed
  store_folder_path_control <- X_list$store_folder_path_control
  label_output_variable_AWP <- X_list$label_output_variable_AWP
  Calculated_output_folder_AWP <- X_list$Calculated_output_folder_AWP
  label_output_variable_ACY <- X_list$label_output_variable_ACY
  Calculated_output_folder_ACY <- X_list$Calculated_output_folder_ACY
  label_output_variable_DWS <- X_list$label_output_variable_DWS
  Calculated_output_folder_DWS <- X_list$Calculated_output_folder_DWS

  while (i <= sample_size[1]) {
    #Wrting basin parameters....
    temp_watershed_file_path <- paste(folder_path_project,
                                      "/", label_watershed_param,
                                      ".dat", sep = "")
    write_param_basin(GSA_Object$X[i,1:length(uncertain_basin_parameters)],
                      temp_watershed_file_path,
                      Back_Up_PARM0806.dat)
    #Writing control parameters...
    temp_control_file_path <- paste(folder_path_project,
                                    "/", label_control_param,
                                    ".dat", sep = "")
    if (length(uncertain_control_parameters) != 0) {
      write_param_cont(GSA_Object$X[i,(length(uncertain_basin_parameters) + 1)
                                    :(length(uncertain_control_parameters) +
                                        length(uncertain_basin_parameters))]
                       ,temp_control_file_path, Back_Up_APEXCONT.dat)
    }
    setwd(folder_path_project)
    system(paste(label_APEX_exe, ".exe", sep = ""))
    setwd(folder_path_R_codes)
    # Parm files...
    file.copy(temp_watershed_file_path,
              paste(store_folder_path_watershed,
                    "/", label_watershed_param,
                    toString(i), ".dat", sep = ""),
              recursive = TRUE, overwrite = TRUE)
    # APEXCONT files...
    file.copy(temp_control_file_path,
              paste(store_folder_path_control,
                    "/", label_control_param,
                    toString(i), ".dat", sep = ""),
              recursive = TRUE, overwrite = TRUE)
    #Output files...
    #.AWP file...
    temp_output_file_path_AWP <- paste(folder_path_project,
                                       "/", label_output_variable_AWP,
                                       ".AWP", sep = "")
    file.copy(temp_output_file_path_AWP,
              paste(Calculated_output_folder_AWP,
                    "/", label_output_variable_AWP,
                    toString(i), ".AWP", sep = ""),
              recursive = TRUE, overwrite = TRUE)
    #.ACY file...
    temp_output_file_path_ACY <- paste(folder_path_project,
                                       "/",label_output_variable_ACY,
                                       ".ACY", sep = "")
    file.copy(temp_output_file_path_ACY,
              paste(Calculated_output_folder_ACY,
                    "/", label_output_variable_ACY,
                    toString(i), ".ACY", sep = "")
              ,recursive = TRUE, overwrite = TRUE)
    #.DWS files...
    temp_output_file_path_DWS <- paste(folder_path_project,
                                       "/", label_output_variable_DWS,
                                       ".DWS", sep = "")
    file.copy(temp_output_file_path_DWS,
              paste(Calculated_output_folder_DWS,
                    "/", label_output_variable_DWS,
                    toString(i),".DWS", sep = ""),
              recursive = TRUE, overwrite = TRUE)
    file.copy(temp_output_file_path_DWS,
              paste(Calculated_output_folder_DWS, "/",
                    "TimeStamp/", as.numeric(Sys.time()),
                    ".DWS", as.character(runif(1)) ,sep = ""),
              recursive = TRUE, overwrite = TRUE)
    i = i + 1
  }
}
#
#####
param_extract <- function(param_range) {
  # Extracts uncertain parameters from a data frame cantaining lower and upper limits.
  a<-dim(param_range)
  numb_param <-a[2]
  i <- 1
  while (i <= numb_param) {
    if (param_range[[i]][1]==param_range[[i]][2]) {
      param_range[i]=NULL
      i=max(1, i - 1)
      numb_param = numb_param - 1
    } else
      i = i + 1
    }
  param_range
}
#
#####
gsa_object <- function(
  sample_size,
  lower_limits,
  upper_limits,
  all_uncertain_model_parameters,
  type) {
  # Creates a GSA_Object containing sensitivity analysis settings.
  if (type == "MORRIS") {
    GSA_Object <- sensitivity::morris(
      model = NULL,
      factor = names(all_uncertain_model_parameters),
      r = as.numeric(SA_Parms$morris_r),
      design = list(
        type = "oat",
        levels = as.numeric(SA_Parms$morris_levels),
        grid.jump = 1),
      binf = lower_limits,
      bsup = upper_limits,
      scale = TRUE)
    GSA_Object$X <- as.data.frame(GSA_Object$X);
    GSA_Object
  } else if (type == "SRC") {
    X<-all_uncertain_model_parameters[-1,]
    X[1,] = array(NaN,dim=c(1,length(X)))
    X <- data_frame_aug(X,sample_size)
    x_size <- dim(X)
    i = 1
    while(i <= x_size[2]) {
      X[i] <- runif(
        sample_size,
        min=lower_limits[i],
        max=upper_limits[i])
      i = i + 1
      }
    GSA_Object = list(X=X)
  } else if (type == "SRRC") {
    X <- all_uncertain_model_parameters[-1,]
    X[1,] = array(NaN,dim=c(1,length(X)))
    X <- data_frame_aug(X,sample_size)
    x_size <- dim(X)
    i = 1
    while(i<=x_size[2]) {
      X[i] <- runif(
        sample_size,
        min=lower_limits[i],
        max=upper_limits[i])
      i = i + 1
    }
    GSA_Object = list(X=X)
  } else if (type=="FAST99") {
    lower_limits <- as.vector(all_uncertain_model_parameters[1,],
                              mode = "numeric")
    upper_limits <- as.vector(all_uncertain_model_parameters[2,],
                              mode = "numeric")
    limit_list <- list()
    i = 1
    while (i <= ncol(all_uncertain_model_parameters)) {
      limit_list[[i]] <- list(
        min = lower_limits[i],
        max = upper_limits[i])
      i = i + 1
    }
    GSA_Object <- sensitivity::fast99(
      model = NULL,
      factor = names(all_uncertain_model_parameters),
      n = sample_size,
      M = 4,
      q = rep("qunif",
              ncol(all_uncertain_model_parameters)),
      q.arg = limit_list)

  } else {
    stop("Please enter a valid method")
  }
  return(GSA_Object)
}
#
#####
data_frame_aug <- function(data_frame, rep_numb) {
  dim_x <- dim(data_frame)
  data_frame_augment <- data_frame
  i = 1
  while (i < rep_numb) {
    data_frame_augment <-rbind(
      data_frame_augment,
      data_frame,
      make.row.names = FALSE)
    i = i + 1
  }
  data_frame_augment
}




