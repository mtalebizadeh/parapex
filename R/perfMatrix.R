#' perfMatrix
#' @description Calculates model performance or simulated mean.
#' @param Input
#'
#' @return A data.frame object containing model performance.
#' @export
#'
#' @examples perfMat <- perfMatrix(Input)
#'
#'
perfMatrix <- function(Input) {
  folder_path_observed <- Input$folder_path_observed
  Calculated_output_folder_DWS <- Input$Calculated_output_folder_DWS
  folder_path_GSA_Outputs  <- Input$folder_path_GSA_Outputs
  folder_path_GSA_Outputs <- Input$folder_path_GSA_Outputs
  i_start <- 1
  i_last <- length(list.files(path = Input$Calculated_output_folder_DWS,
                              pattern = "*.DWS"))
  start_date <- gsub(pattern = " ",
                     replacement = "-",
                     x = Input$start_date)
  end_date <- gsub(pattern = " ",
                   replacement = "-",
                   x = Input$end_date)
  window_length <- Input$window_length
  caption_var_obs <- Input$caption_var_obs
  caption_var_sim <- Input$caption_var_sim
  label_observed_var.txt <- Input$label_observed_var.txt
  label_output_variable_DWS <- Input$label_output_variable_DWS
  func_type <- Input$func_type
  TW <- Input$TW
  number_cores <- Input$coreNumber
  # Calculating model performance...
  Perf_Matrix_Final <- perf_matrix_dync_proTW(
    i_start,
    i_last,
    start_date,
    end_date,
    window_length,
    caption_var_obs,
    caption_var_sim,
    folder_path_observed,
    Calculated_output_folder_DWS,
    label_observed_var.txt,
    label_output_variable_DWS,
    func_type = func_type,
    TW = TW,
    number_cores = number_cores)
  # Writing performance matrix
  write.table(
    x = Perf_Matrix_Final,
    file = paste(
      folder_path_GSA_Outputs,
      "/", "Perf_mat", "_",
      window_length,"_","TW",
      "_",TW,".txt",sep = ""),
    sep = "\t",row.names = FALSE)

  Perf_Matrix_Final
}
#
#### Helper functions
perf_matrix_dync_proTW <- function(
  i_start,
  i_last,
  start_date,
  end_date,
  window_length,
  caption_var_obs,
  caption_var_sim,
  folder_path_observed,
  Calculated_output_folder_DWS,
  label_observed_var.txt,
  label_output_variable_DWS,
  func_type,
  TW,
  number_cores) {
  # Calculates dynamic performance matrix fo the analysis period.
  #
  # Args:
  #   i_start: An integer indicating the first file index to be processed.
  #   i_last:  An integer indicating the last file index to be processed.
  #   start_date: A charachter indicating the start of analysis period.
  #       The accepted format is: "YYYY MM DD".
  #   end_date: A charachter indicating the end of analysis period.
  #       The accepted format is: "YYYY MM DD".
  #   window_length: A charachter indicating the analysis sub-periods. Can take
  #       "day", "week", "month", "year" or any multiples of them, e.g. "3 day".
  #   caption_var_obs: Observed variable name as it appears inside the observed file.
  #   caption_var_sim: Simulated variable name as it appears inside the .DWS file.
  #   folder_path_observed: Folder path to where the observed file resides.
  #   Calculated_output_folder_DWS: Folder path to where the simulated .DWS files
  #       are stored.
  #   label_observed_var.txt: File name containing observed data.
  #   label_output_variable_DWS: .DWS file name.
  #   func_type: Function type used for calculating model performance. Can take:
  #       "MEAN", "RMSE", "NASH", "PBIAS".
  #   TW: lenght of time aggregation period. Can take
  #       "day", "week", "month", "year" or any multiples of them, e.g. "3 day".
  #   number_cores: An integer representing number of computation nodes.
  #
  # Returns:
  #     A data.frame containing performance values during the analysis period.
  x_lister <- function(
    i_start,
    i_last,
    start_date,
    end_date,
    window_length,
    caption_var_obs,
    caption_var_sim,
    folder_path_observed,
    Calculated_output_folder_DWS,
    label_observed_var.txt,
    label_output_variable_DWS,
    func_type,
    TW) {
    # Prepares inputs for 'dynamic_performance_matrix' function.
    X_list=list(i_start=i_start)
    X_list$i_last <- i_last
    X_list$start_date <- start_date
    X_list$end_date <- end_date
    X_list$window_length <- window_length
    X_list$caption_var_obs <- caption_var_obs
    X_list$caption_var_sim <- caption_var_sim
    X_list$folder_path_observed <- folder_path_observed
    X_list$Calculated_output_folder_DWS <- Calculated_output_folder_DWS
    X_list$label_observed_var.txt <- label_observed_var.txt
    X_list$label_output_variable_DWS <- label_output_variable_DWS
    X_list$func_type <- func_type
    X_list$TW <-TW
    return(X_list)
  }
  dynamic_performance_matrix <- function(X_list) {
    # Calculates performance matrix using input list.
    #
    # Args:
    #   X_list: A list containing analysis inputs.
    #
    # Returns:
    #   A data.frame containing model performance values.
    date_seq <- as.character(
      seq.Date(from = as.Date(X_list$start_date),
               to = as.Date(X_list$end_date),
               by = X_list$window_length))

    date_seq_start <- date_seq[1:(length(date_seq)-1)]
    date_seq_daily <- as.character(
      seq.Date(from = as.Date(X_list$start_date),
               to=as.Date(end_date),
               by="day"))
    date_idx <- 1:length(date_seq_daily)
    date_seq_end <- date_seq[-1]

    i = 1
    while (i <= length(date_seq_end)) {
      end_date_idx_temp <- date_idx[date_seq_daily == date_seq_end[i]]
      end_date_idx_temp <- end_date_idx_temp-1
      date_seq_end[i] <- date_seq_daily[end_date_idx_temp]
      i = i + 1
    }
    # Generating an empty performance matrix
    performance_dynamic <- as.data.frame(
      array(NaN,
            dim = c(X_list$i_last,
                    length(date_seq_start))))
    names(performance_dynamic) <- date_seq_start
    # Filling the matrix
    i <- X_list$i_start
    while (i <= X_list$i_last) {
      j = 1
      while(j <= ncol(performance_dynamic)) {
        temp_simulated_output <- dws2timeseries(
          paste(X_list$Calculated_output_folder_DWS,
                "/", X_list$label_output_variable_DWS,
                toString(i),".DWS",
                sep = ""),
          X_list$caption_var_sim,
          date_seq_start[j], date_seq_end[j]
        )
        if (func_type=="MEAN") {
          temp_observed_output=-999
        } else {
          temp_observed_output <- obs2timeseries(
            paste(X_list$folder_path_observed,
                  "/", X_list$label_observed_var.txt,
                  sep = ""),
            X_list$caption_var_obs,
            date_seq_start[j],
            date_seq_end[j])
        }
      performance_dynamic[i,j] <- performance_func_proTW(
        temp_simulated_output,
        temp_observed_output,
        date_seq_start[j],
        date_seq_end[j],
        TW = X_list$TW,
        func_type = X_list$func_type)
      j = j + 1
      }
    i  = i + 1
    }
  return(performance_dynamic)
  }
#
#####
  seq_2_par_dyn_performance  <-  function(X_list,number_cores) {
    # Converts inputs for sequential simulation to parallel.
    #
    # Args:
    #   X_list: List of simulation inputs.
    #   number_cores: Number of cores to be used for parallel simulation.
    #
    # Returns: List of slimulation inputs to be used for computation nodes.
    #
    batch_size_average <- floor(X_list$i_last/number_cores)
    solid_length <- batch_size_average * number_cores
    batch_start_idx <- seq(1,
                           solid_length,
                           by = batch_size_average)
    batch_end_idx <- batch_start_idx - 1
    batch_end_idx <- batch_end_idx[-1]
    batch_end_idx[length(batch_end_idx) + 1] <- X_list$i_last
    XX_list <- list(X_list)
    j = 1
    while (j <= length(batch_start_idx)) {
      XX_list[[j]] <- X_list
      XX_list[[j]]["i_start"] <- batch_start_idx[j]
      XX_list[[j]]["i_last"] <- batch_end_idx[j]
      j = j + 1
      }

    return(XX_list)
  }

  X_list <- x_lister(
    i_start,
    i_last,
    start_date,
    end_date,
    window_length,
    caption_var_obs,
    caption_var_sim,
    folder_path_observed,
    Calculated_output_folder_DWS,
    label_observed_var.txt,
    label_output_variable_DWS,
    func_type,
    TW)

  XX_list <- seq_2_par_dyn_performance(
    X_list,
    number_cores)

  cl <- parallel::makeCluster(number_cores,
                              type = "PSOCK")
  Mat_list= parallel::parLapply(cl,
                                XX_list,
                                dynamic_performance_matrix)
  parallel::stopCluster(cl)
  i = 1
  performance_matrix_final <- na.omit(Mat_list[[i]])
  while(i < number_cores) {
    performance_matrix_final <- rbind(
      performance_matrix_final,
      na.omit(Mat_list[[i + 1]]))
    i = i + 1
  }
  return(performance_matrix_final)
}
#
#####
performance_func_proTW <- function(
  sim_series,
  obs_series,
  start_date,
  end_date,
  TW = "day",
  func_type) {
  # Claculates model performance using simulated and observed data.
  #
  # Args:
  #   sim_series: A dataframe containing simulated variable.
  #   obs_series: A dataframe containing observed variable.
  #   start_date: A character representing start date of analysis.
  #   end_date: A character representing end date of analysis.
  #   TW: lenght of time aggregation period. Can take
  #       "day", "week", "month", "year" or any multiples of them, e.g. "3 day".
  #   func_type: Function type used for calculating model performance. Can take:
  #       "MEAN", "RMSE", "NASH", "PBIAS".
  switch (func_type,
          "RMSE"=return(output_RMSE_proTW(sim_series,
                                          obs_series,
                                          start_date,
                                          end_date,
                                          TW)),
          "NASH"=return(output_NASH_proTW(sim_series,
                                          obs_series,
                                          start_date,
                                          end_date,
                                          TW)),
          "PBIAS"=return(output_PBIAS_proTW(sim_series,
                                            obs_series,
                                            start_date,
                                            end_date,
                                            TW)),
          "MEAN"= return(output_MEAN_proTW(sim_series,
                                           start_date,
                                           end_date,
                                           TW)))
}
#
##### Helper functions
output_RMSE_proTW <- function(
  sim_series,
  obs_series,
  start_date,
  end_date, TW="day") {
  # Calculates RMSE
  sim_series_agg <- agg4timeseries(
    time_series = sim_series,
    start_date = start_date,
    end_date = end_date,
    TW = TW)
  obs_series_agg <- agg4timeseries(
    time_series = obs_series,
    start_date = start_date,
    end_date = end_date,
    TW = TW)
  RMSE_series <- sqrt(mean((sim_series_agg[[2]] -
                              obs_series_agg[[2]]) ^ 2,
                           na.rm = TRUE))
  return(RMSE_series)
}
#
#####
output_NASH_proTW <- function(
  sim_series,
  obs_series,
  start_date,
  end_date,
  TW) {
  # Calculates Nash-Sutcliffe performance measure.
  sim_series_agg <- agg4timeseries(
    time_series = sim_series,
    start_date = start_date,
    end_date = end_date,
    TW = TW)
  obs_series_agg <- agg4timeseries(
    time_series = obs_series,
    start_date = start_date,
    end_date = end_date, TW = TW)
  NS_Numerator <- sum((sim_series_agg[[2]] -
                         obs_series_agg[[2]]) ^ 2,
                      na.rm = TRUE)
  NS_Denominator <- sum((obs_series_agg[[2]] -
                           mean(obs_series_agg[[2]],
                                na.rm = TRUE))^2,
                        na.rm = TRUE)
  NASH <- 1 - (NS_Numerator/NS_Denominator)
  return(NASH)
}
#
#####
output_PBIAS_proTW <- function(
  sim_series,
  obs_series,
  start_date,
  end_date,
  TW) {
  # Calculates PBIAS performance measure
  sim_series_agg <- agg4timeseries(
    time_series = sim_series,
    start_date = start_date,
    end_date = end_date,
    TW = TW)
  obs_series_agg <- agg4timeseries(
    time_series = obs_series,
    start_date = start_date,
    end_date = end_date,
    TW = TW)

  PB_Numerator <- sum(obs_series_agg[[2]] -
                        sim_series_agg[[2]],
                      na.rm = TRUE) * 100
  PB_Denominator <- sum(obs_series_agg[[2]],
                        na.rm = TRUE)
  PBIAS <- PB_Numerator/PB_Denominator
  return(PBIAS)
}
output_MEAN_proTW <- function(
  sim_series,
  start_date,
  end_date,
  TW="day") {
  # Calculates simulated mean.
  sim_series_agg <- agg4timeseries(
    time_series = sim_series,
    start_date = start_date,
    end_date = end_date,
    TW = TW)
  MEAN_series <- mean(sim_series_agg[[2]] ,
                      na.rm = TRUE)
  return(MEAN_series)
}
#
#####
agg4timeseries <- function(
  time_series,
  start_date,
  end_date,
  TW ="day") {
  # Aggregates simulated or observed time series.
  #
  # Args:
  #   time_series: A dataframe containing simulated/observed data.
  #   start_date: Start date of analysis.
  #   end_date: End date of analysis.
  #   TW: Length of time aggregation.
  #
  # Returns:
  #   A dataframe containing aggregated time series.
  #
  if (TW == "day") {
    return(time_series)
  } else {
    agg_date_seq <- as.character(
      seq.Date(from = as.Date(start_date),
               to = as.Date(end_date), by = TW))
    idx_seq <- 1:nrow(time_series)
    idx_seq_start <- idx_seq[time_series$Date %in%
                               agg_date_seq[1:(length(agg_date_seq) - 1)]]
    idx_seq_end <- (idx_seq[time_series$Date %in%
                              agg_date_seq][-1]) - 1
    idx_seq_end[length(idx_seq_end)] <- idx_seq_end[length(idx_seq_end)] + 1
    agg_timeseries <- as.data.frame(array(
      data = NaN,
      dim = c(length(idx_seq_start),2)))
    names(agg_timeseries) <- names(time_series)
    agg_timeseries$Date <- agg_date_seq[1:length(agg_date_seq)-1]
    i = 1
    while (i <= nrow(agg_timeseries)) {
      agg_timeseries[[2]][i] <- mean(c(time_series[[2]][idx_seq_start[i]:
                                                          idx_seq_end[i]]))
      i = i + 1
    }
    return(agg_timeseries)
  }
}
#
#####
obs2timeseries <- function(
  file_name,
  caption_obs_var,
  start_date,
  end_date) {
  # Creates time series of observed data.
  #
  # Args:
  #   File name containing observed data.
  #   caption_obs_var: Column name of observed variable.
  #   start_date: Start date of analysis.
  #   end_date: End date of analysis.
  #
  # Returns:
  #   A dataframe containing observed data.
  #
  A_temp <- read.table(file_name, header = TRUE)
  first_dws_date <- paste(
    substr(x = A_temp$Date[1], start = 1,stop = 4), "-",
    substr(x = A_temp$Date[1], start = 5,stop = 6), "-",
    substr(x = A_temp$Date[1],start = 7,stop = 8),
    sep = "")
  last_dws_date <- paste(
    substr(x = A_temp$Date[nrow(A_temp)], start = 1,stop = 4), "-",
    substr(x= A_temp$Date[nrow(A_temp)], start = 5, stop = 6), "-",
    substr(x= A_temp$Date[nrow(A_temp)], start = 7, stop = 8),
    sep = "")
  A_temp$Date <- as.character(seq.Date(
    as.Date(first_dws_date),
    as.Date(last_dws_date),
    by = "day"))
  date_seq <- as.character(seq.Date(
    as.Date(start_date),
    as.Date(end_date),
    by = "day"))
  idx_all_dates <- 1:nrow(A_temp)
  idx_start_date <- idx_all_dates[A_temp$Date == start_date]
  idx_end_date <- idx_all_dates[A_temp$Date == end_date]
  time_series_obs <- array(NaN,
                           dim = c(length(date_seq),
                                   length(caption_obs_var) + 1))
  time_series_obs[, 1] <- date_seq
  col_names <- c("Date", caption_obs_var)
  colnames(time_series_obs) <- col_names
  time_series_obs <- as.data.frame(time_series_obs)
  i = 1
  while (i <= length(col_names) - 1) {
    time_series_obs[col_names[i + 1]] <-
      A_temp[caption_obs_var[i]][idx_start_date:idx_end_date,]
    i = i + 1
  }
  return(time_series_obs)
}
#
#####
dws2timeseries <- function(
  file_name,
  caption_dws_var,
  start_date,
  end_date) {
  # Creates simulation time series from .DWS output files.
  #
  # Args:
  #   file_name: DWS utput file name containing simulated variables.
  #   caption_dws_var: Variable name as it appears inside the DWS files.
  #   start_date: Start date of analysis.
  #   end_date: End date of analysis.
  #
  # Returns:
  #   A dataframe containing simulated time series.
  A_temp = read.table(file_name,
                      header = TRUE,
                      skip = 8)
  first_dws_date <- paste(toString(A_temp$Y[1]), "-",
                          toString(A_temp$M[1]), "-",
                          toString(A_temp$D[1]), sep = "")
  last_dws_date <- paste(toString(A_temp$Y[nrow(A_temp)]), "-",
                         toString(A_temp$M[nrow(A_temp)]), "-",
                         toString(A_temp$D[nrow(A_temp)]),
                         sep = "")
  A_temp$Date <- as.character(seq.Date(
    as.Date(first_dws_date),
    as.Date(last_dws_date),
    by = "day"))
  if (caption_dws_var == "TN") {
    A_temp$TN <- A_temp$QN + A_temp$YN +
      A_temp$QDRN + A_temp$RSFN +
      A_temp$QRFN + A_temp$SSFN
  }
  if (caption_dws_var == "TP") {
    A_temp$TP <- A_temp$QP + A_temp$YP +
      A_temp$QDRP + A_temp$QRFP
  }

  date_seq <- as.character(seq.Date(
    as.Date(start_date),
    as.Date(end_date),
    by = "day"))
  idx_all_dates <- 1:nrow(A_temp)
  idx_start_date <- idx_all_dates[A_temp$Date == start_date]
  idx_end_date <- idx_all_dates[A_temp$Date == end_date]
  time_series_dws <- array(NaN,
                           dim = c(length(date_seq),
                                   length(caption_dws_var) + 1))
  time_series_dws[, 1] <- date_seq
  col_names <- c("Date", caption_dws_var)
  colnames(time_series_dws) <- col_names
  time_series_dws <- as.data.frame(time_series_dws)
  i = 1
  while (i <= length(col_names) - 1) {
    time_series_dws[col_names[i + 1]] <-
      A_temp[caption_dws_var[i]][idx_start_date:idx_end_date,]
    i = i + 1
  }
  return(time_series_dws)
}
