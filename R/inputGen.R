
#' inputGen
#' @description Creates a list object for setting inputs.
#' @return A list object containing analysis inputs.
#' @export
#'
#' @examples Input <- inputGen()
inputGen <- function() {
  Simulation_Inputs <- list(sample_size = 100)
  Simulation_Inputs$caption_var_sim = "ET"
  Simulation_Inputs$caption_var_obs <- "ET"
  Simulation_Inputs$start_date <- "2002 01 01"
  Simulation_Inputs$end_date <- "2003 01 01"
  Simulation_Inputs$label_APEX_exe <- "APEX0806"
  Simulation_Inputs$label_watershed_param <- "PARM"
  Simulation_Inputs$label_control_param <- "Apexcont"
  Simulation_Inputs$label_output_variable_AWP <- "APEX001"
  Simulation_Inputs$label_output_variable_ACY <- "APEX001"
  Simulation_Inputs$label_output_variable_DWS <- "APEX001"
  Simulation_Inputs$label_prior_dist <- "prior_dist.txt"
  Simulation_Inputs$label_observed_var.txt <- "Observed.txt"
  Simulation_Inputs$folder_path_project <- "Example/APEX"
  Simulation_Inputs$Back_Up_PARM0806.dat <- "Example/Back_Up/PARM.dat"
  Simulation_Inputs$Back_Up_APEXCONT.dat <- "Example/Back_Up/Apexcont.dat"
  Simulation_Inputs$folder_path_R_codes <- getwd()
  Simulation_Inputs$folder_path_observed <- "Example/Observed_Files"
  Simulation_Inputs$folder_path_GSA_Outputs <- "Example/GSA_Outputs"
  Simulation_Inputs$store_folder_path_watershed <- "Example/Generated_Inputs/PARM"
  Simulation_Inputs$store_folder_path_control <- "Example/Generated_Inputs/CONT"
  Simulation_Inputs$Calculated_output_folder_AWP <- "Example/Calculated_Outputs/AWP"
  Simulation_Inputs$Calculated_output_folder_ACY <- "Example/Calculated_Outputs/ACY"
  Simulation_Inputs$Calculated_output_folder_DWS <- "Example/Calculated_Outputs/DWS"
  Simulation_Inputs$GSA_Type <- "FAST99"
  Simulation_Inputs$coreNumber <- parallel::detectCores()
  Simulation_Inputs$window_length <- "year"
  Simulation_Inputs$TW = "3 day"
  Simulation_Inputs$func_type <- "MEAN"
  #
  ##SA method-specific parameters
  #
  #Morris
  Simulation_Inputs$SA_Parms$morris_r <- 15
  Simulation_Inputs$SA_Parms$morris_levels <- 8
  Simulation_Inputs$APEX_PARM <- empty_parm_range_file()
  Simulation_Inputs$APEX_CONT <- empty_cont_range_file()
  Simulation_Inputs
}

# Helper functions
# Creates a dataframe containing 'APEXCONT' parameters.
empty_cont_range_file <- function () {

  watershed_cont_dataframe <- data.frame(RFN = array(-1, dim = c(2,1)))
  watershed_cont_dataframe$CO2 <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$CQN <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$PSTX <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$YWI <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$BTA <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$EXPK <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$QG <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$QCF <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$CHSO <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$BWD <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$FCW <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$FPSC <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$GWSO <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$RFTO <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$RFPO <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$SATO <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$FL <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$FW <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$ANG <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$UXP <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$DIAM <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$ACW <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$GZL0 <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$RTN0 <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$BXCT <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$BYCT <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$DTHY <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$QTH <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$STND <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$DRV <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$PCO0 <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$RCC0 <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$CSLT <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$BUS1 <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$BUS2 <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe$BUS3 <- array(-1, dim = c(2, 1))
  watershed_cont_dataframe
}

# Creates a dataframe containing 'PARM' parameters.
#
empty_parm_range_file <- function(){
  watershed_param <- data.frame(Crop_canopy_PET = array(-1, dim = c(2, 1)))
  watershed_param$Root_growth_soil = array(-1, dim = c(2, 1))
  watershed_param$Water_stress_harvest = array(-1, dim = c(2,1))
  watershed_param$Water_storage_N = array(-1, dim = c(2, 1))
  watershed_param$Soil_water_limit = array(-1, dim = c(2, 1))
  watershed_param$Winter_dormancy = array(-1, dim = c(2, 1))
  watershed_param$N_fixation = array(-1, dim = c(2, 1))
  watershed_param$Soluble_P_runoff = array(-1, dim = c(2, 1))
  watershed_param$Pest_damage_moisture = array(-1, dim = c(2,1))
  watershed_param$Pest_damage_cover = array(-1, dim = c(2,1))
  watershed_param$Moisture_req_seed_germ = array(-1, dim = c(2,1))
  watershed_param$Soil_evap_coeff = array(-1, dim = c(2, 1))
  watershed_param$Wind_erod_coeff = array(-1, dim = c(2, 1))
  watershed_param$Nitrate_leac_ratio = array(-1, dim = c(2,1))
  watershed_param$Runoff_CN_Adj_parm = array(-1, dim = c(2,1))
  watershed_param$Expand_CN_ret_parm = array(-1, dim = c(2,1))
  watershed_param$Soil_evap_plant_cover = array(-1, dim = c(2,1))
  watershed_param$Sedim_rout_exponent = array(-1, dim = c(2,1))
  watershed_param$Sedim_rout_coeff = array(-1, dim = c(2, 1))
  watershed_param$Runoff_CN_int_abs = array(-1, dim = c(2,1))
  watershed_param$Soluble_C_adsorp_coeff = array(-1, dim = c(2,1))
  watershed_param$CN_retention_frozen_soil = array(-1, dim = c(2,1))
  watershed_param$Harg_equation_parm = array(-1, dim = c(2,1))
  watershed_param$Pest_leach_ratio = array(-1, dim = c(2, 1))
  watershed_param$Expo_coeff_rainfall = array(-1, dim = c(2,1))
  watershed_param$Matur_frac_spring = array(-1, dim = c(2,1))
  watershed_param$CEC_effect_nitrification = array(-1, dim = c(2,1))
  watershed_param$N_fixation_limit = array(-1, dim = c(2, 1))
  watershed_param$Biological_mix_efficiency = array(-1, dim = c(2,1))
  watershed_param$Soluble_P_exponent = array(-1, dim = c(2,1))
  watershed_param$Max_depth_bio_mixing = array(-1, dim = c(2,1))
  watershed_param$OrgP_loss_exponent = array(-1, dim = c(2,1))
  watershed_param$MUST_coeff = array(-1, dim = c(2, 1))
  watershed_param$Harg_PET_exponent = array(-1, dim = c(2,1))
  watershed_param$Denit_soil_threshold = array(-1, dim = c(2,1))
  watershed_param$Daily_denit_limit = array(-1, dim = c(2,1))
  watershed_param$SWAT_delivery_ratio_exponent = array(-1, dim = c(2,1))
  watershed_param$Water_stress_coeff = array(-1, dim = c(2, 1))
  watershed_param$Puddling_sat_conduct = array(-1, dim = c(2,1))
  watershed_param$Groundwater_stor_threshold = array(-1, dim = c(2,1))
  watershed_param$Root_temp_stress_exponent = array(-1, dim = c(2,1))
  watershed_param$SCS_index_coeff = array(-1, dim = c(2,1))
  watershed_param$Plow_depth = array(-1, dim = c(2,1))
  watershed_param$CN_retention_param = array(-1, dim = c(2,1))
  watershed_param$sediment_rout_travel_coeff = array(-1, dim = c(2,1))
  watershed_param$RUSLE_c_factor_res = array(-1, dim = c(2,1))
  watershed_param$RUSLE_c_factor_height = array(-1, dim = c(2,1))
  watershed_param$Climate_stress_factor = array(-1, dim = c(2,1))
  watershed_param$Max_rain_intercept = array(-1, dim = c(2,1))
  watershed_param$Rain_intercept_coeff = array(-1, dim = c(2,1))
  watershed_param$Water_stor_residue_coeff = array(-1, dim = c(2,1))
  watershed_param$Tillage_residue_decay_rate_coeff = array(-1,dim = c(2, 1))
  watershed_param$Microbial_soil_depth_coeff = array(-1, dim = c(2,1))
  watershed_param$N_enrich_coeff = array(-1, dim = c(2, 1))
  watershed_param$N_enrich_rout_exponent = array(-1, dim = c(2,1))
  watershed_param$Fraction_destroyed_burn = array(-1, dim = c(2,1))
  watershed_param$P_enrich_rout_coeff = array(-1, dim = c(2,1))
  watershed_param$P_enrich_rout_exponent = array(-1, dim = c(2,1))
  watershed_param$P_move_evap_coeff = array(-1, dim = c(2,1))
  watershed_param$Max_days_grazed_rotation = array(-1, dim = c(2,1))
  watershed_param$Soil_water_up_flow_limit = array(-1, dim = c(2,1))
  watershed_param$Manure_erosion_equation_coeff = array(-1,dim = c(2, 1))
  watershed_param$N_enrich_ratio_delivery = array(-1, dim = c(2, 1))
  watershed_param$Dust_distribution_coeff = array(-1, dim = c(2,1))
  watershed_param$RUSLE2_trans_capacity = array(-1, dim = c(2,1))
  watershed_param$RUSLE2_trans_capacity_threshold = array(-1,dim = c(2, 1))
  watershed_param$Dust_distribution_exponent = array(-1, dim = c(2,1))
  watershed_param$Manure_erosion_exponent = array(-1, dim = c(2,1))
  watershed_param$Microbial_top_soil_coeff = array(-1, dim = c(2,1))
  watershed_param$Microbial_decay_coeff = array(-1, dim = c(2,1))
  watershed_param$Manure_erosion_coeff = array(-1, dim = c(2,1))
  watershed_param$Volt_nitrification_partition_coeff = array(-1,dim = c(2, 1))
  watershed_param$Hydrograph_dev_param = array(-1, dim = c(2,1))
  watershed_param$Partition_N_flow_groundwater = array(-1,dim = c(2,1))
  watershed_param$P_enrich_ratio_deliver_SWAT = array(-1, dim = c(2,1))
  watershed_param$Stand_dead_fall_rate_coeff = array(-1, dim = c(2,1))
  watershed_param$Runoff_2_delay_pest = array(-1, dim = c(2,1))
  watershed_param$Soil_water_2_delay_tillage = array(-1, dim = c(2,1))
  watershed_param$Auto_mov_lower_limit = array(-1, dim = c(2,1))
  watershed_param$Nitrification_vol_upper_limit = array(-1,dim = c(2, 1))
  watershed_param$Tech_coeff = array(-1, dim = c(2, 1))
  watershed_param$Drainage_lateral_conduct = array(-1, dim = c(2,1))
  watershed_param$P_flux_labile_active_coeff = array(-1, dim = c(2,1))
  watershed_param$P_flux_active_stable_coeff = array(-1, dim = c(2,1))
  watershed_param$N_salt_evap_coeff = array(-1, dim = c(2,1))
  watershed_param$Water_table_recession_coeff = array(-1, dim = c(2,1))
  watershed_param$Water_table_move_limit = array(-1, dim = c(2,1))
  watershed_param$Water_table_recession_exponent = array(-1,dim = c(2, 1))
  watershed_param$Subsurface_flow_factor = array(-1, dim = c(2,1))
  watershed_param$Flood_evap_limit = array(-1, dim = c(2, 1))
  watershed_param$Runoff_adj_link = array(-1, dim = c(2, 1))
  watershed_param$Water_erosion_threshold = array(-1, dim = c(2,1))
  watershed_param$Wind_erosion_threshold = array(-1, dim = c(2,1))
  watershed_param$Crop_stress_temp_exponent = array(-1, dim = c(2,1))
  watershed_param$Soluble_P_leach_KD = array(-1, dim = c(2,1))
  watershed_param$Unknown1 = array(-1, dim = c(2, 1))
  watershed_param$Unknown2 = array(-1, dim = c(2, 1))
  watershed_param$Unknown3 = array(-1, dim = c(2, 1))
  watershed_param$Irrigation_cost = array(-1, dim = c(2, 1))
  watershed_param$Lime_cost = array(-1, dim = c(2, 1))
  watershed_param$Fuel_cost = array(-1, dim = c(2, 1))
  watershed_param$Labor_cost = array(-1, dim = c(2, 1))
  watershed_param$Unknown4 = array(-1, dim = c(2, 1))
  watershed_param
}






