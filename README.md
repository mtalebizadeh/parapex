parapex quick start manual
================


Introduction
------------

Parapex is a package in R for performing dynamic sensitivity and model performance analysis for the APEX model.The package was primarily developed for USDA-ARS Grazinglands Research Laboratory to assist in parametrization of the APEX model. For more information on the structure of parapex and citation, please refer to: [Talebizadeh et al. (2018)](https://onlinelibrary.wiley.com/doi/full/10.1111/1752-1688.12686)

Example folder
--------------

An example folder containing an APEX project and other inputs is available for users to test the package. The rest of this manual provides details of implementing an SA project using the accompanying example folder which can be created through a call to:

``` r
# Creating a copy of tutorial folder inside the working directory
  parapex::getExampleFolder()
```

Steps for performing SA
-----------------------

After loading parapex and generating a copy of the example folder, the following four steps, described in the next sections should be followed for performing SA.

``` r
# 1) Generating a list object with a predefined structure compatible to APEXSENSUN:
     Input <- parapex::inputGen()

# 2) Setting the required inputs (e.g. uncertainty boubds, SA method, sample size, ...)
  #
  # Setting uncertainty bounds:
    Input$APEX_PARM$Root_growth_soil[1]=1.1
    Input$APEX_PARM$Root_growth_soil[2]=1.20

    Input$APEX_PARM$Soil_water_limit[1]=0
    Input$APEX_PARM$Soil_water_limit[2]=1

    Input$APEX_PARM$Soil_evap_coeff[1]=1.5
    Input$APEX_PARM$Soil_evap_coeff[2]=2.5

    Input$APEX_PARM$Soil_evap_plant_cover[1]=0
    Input$APEX_PARM$Soil_evap_plant_cover[2]=0.5

    Input$APEX_PARM$Runoff_CN_int_abs[1]=0.05
    Input$APEX_PARM$Runoff_CN_int_abs[2]=0.4

    Input$APEX_PARM$Max_rain_intercept[1]=0
    Input$APEX_PARM$Max_rain_intercept[2]=15

    Input$APEX_PARM$Rain_intercept_coeff[1]=0.05
    Input$APEX_PARM$Rain_intercept_coeff[2]=0.3

    Input$APEX_PARM$Microbial_top_soil_coeff[1]=0.1
    Input$APEX_PARM$Microbial_top_soil_coeff[2]=1

    Input$APEX_PARM$Microbial_decay_coeff[1]=0.5
    Input$APEX_PARM$Microbial_decay_coeff[2]=1.5
  
  # SA method and sample size:
    Input$GSA_Type <- "FAST99"
    Input$sample_size = 1000

# 3) Performing Monte Carlo simulation using the setting Input:
     GSA <- parapex::MC4APEX(Input)
    
# 4) Calculation of performance matrix:
     Input$func_type <- "MEAN"
     perfMat <- parapex::perfMatrix(Input)
     
# 5) Calculation of sensitivity indices:
     SA <- parapex::SA4APEX(GSA, perfMat)
```

APENDIX
-------

This section provides 4 main tables containing the name of different parameters and their description.

| Parameters                      | Description                                                          |
|:--------------------------------|:---------------------------------------------------------------------|
| sample\_size                    | Sample size or length of the discretization of parameter space.      |
| caption\_var\_sim               | Simulated variable name as it appears inside .DWS APEX file.         |
| caption\_var\_obs               | Observed variable names as appears inside observed file.             |
| start\_date                     | Start date for analysis with format: YYYY MM DD e.g., "2002 01 01"   |
| end\_date                       | End date for analysis with format: YYYY MM DD e.g., "2003 01 01"     |
| label\_APEX\_exe                | APEX executable file name excluding file's extension.                |
| label\_watershed\_param         | APEX PARM file name excluding file's extension.                      |
| label\_control\_param           | APEXCONT file name excluding file's extension.                       |
| label\_output\_variable\_AWP    | APEX .AWP file name excluding file's extension.                      |
| label\_output\_variable\_ACY    | APEX .ACY file name excluding file's extension.                      |
| label\_output\_variable\_DWS    | APEX .DWS file name excluding file's extension.                      |
| label\_prior\_dist              | File name containing sampled parameters for MC simulation.           |
| label\_observed\_var.txt        | Observed file name containing observed time series.                  |
| folder\_path\_project           | Path to original APEX file containing PARM parameters.               |
| Back\_Up\_PARM0806.dat          | Path to original APEX file containing APEXCONT parameters.           |
| Back\_Up\_APEXCONT.dat          | Path to folder containing APEX model.                                |
| folder\_path\_R\_codes          | No need to set!                                                      |
| folder\_path\_observed          | Path to folder containing observed data file.                        |
| folder\_path\_GSA\_Outputs      | Path to folder storing SA results.                                   |
| store\_folder\_path\_watershed  | Path to folder storing generated PARM files for Monte Carlo runs.    |
| store\_folder\_path\_control    | Path to folder storing generated APEXCONT file for Monte Carlo runs. |
| Calculated\_output\_folder\_AWP | Path to folder storing calculated .AWP files for Monte Carlo runs.   |
| Calculated\_output\_folder\_ACY | Path to folder storing calculated .ACY files for Monte Carlo runs.   |
| Calculated\_output\_folder\_DWS | Path to folder storing calculated .DWS files for Monte Carlo runs.   |
| GSA\_Type                       | Type of SA method: MORRIS, SRC, SRRC, FAST99.                        |
| coreNumber                      | Number of computation cores.                                         |
| window\_length                  | lenght of analysis sub-periods                                       |
| TW                              | length of time aggregation window                                    |
| func\_type                      | Type of model performance measure: "MEAN", "PBIAS", "NASH", "RMSE".  |

</br>

| Parameters     | Description                                                                                                  |
|:---------------|:-------------------------------------------------------------------------------------------------------------|
| morris\_r      | an integer representing design repetition number (i.e. the number of elementary effect computed per factor). |
| morris\_levels | an integer specifying the number of levels of the design in OAT (Once At a Time) design.                     |

</br>

| Parameters                            | Description                                                                         |
|:--------------------------------------|:------------------------------------------------------------------------------------|
| Crop\_canopy\_PET                     | Crop canopy-PET                                                                     |
| Root\_growth\_soil                    | Root growth-soil strength                                                           |
| Water\_stress\_harvest                | Water stress-harvest index                                                          |
| Water\_storage\_N                     | Water storage N leaching                                                            |
| Soil\_water\_limit                    | Soil water lower limit                                                              |
| Winter\_dormancy                      | Winter dormancy                                                                     |
| N\_fixation                           | N fixation                                                                          |
| Soluble\_P\_runoff                    | Soluble phosphorus runoff coefficient                                               |
| Pest\_damage\_moisture                | Pest damage moisture threshold                                                      |
| Pest\_damage\_cover                   | Pest damage cover threshold                                                         |
| Moisture\_req\_seed\_germ             | Moisture required for seed germination                                              |
| Soil\_evap\_coeff                     | Soil evaporation coefficient                                                        |
| Wind\_erod\_coeff                     | Wind erodibility coefficient                                                        |
| Nitrate\_leac\_ratio                  | Nitrate leaching ratio                                                              |
| Runoff\_CN\_Adj\_parm                 | Runoff CN Residue Adjustment Parameter                                              |
| Expand\_CN\_ret\_parm                 | Expands CN retention parameter                                                      |
| Soil\_evap\_plant\_cover              | Soil evaporation â<U+0080><U+0093> plant cover factor                               |
| Sedim\_rout\_exponent                 | Sediment routing exponent                                                           |
| Sedim\_rout\_coeff                    | Sediment routing coefficient                                                        |
| Runoff\_CN\_int\_abs                  | Runoff curve number initial abstraction                                             |
| Soluble\_C\_adsorp\_coeff             | Soluble Carbon adsorption Coefficient                                               |
| CN\_retention\_frozen\_soil           | Reduces NRCS Runoff CN Retention Parameter for Frozen Soil                          |
| Harg\_equation\_parm                  | Hargreaves PET equation coefficient                                                 |
| Pest\_leach\_ratio                    | Pesticide leaching ratio                                                            |
| Expo\_coeff\_rainfall                 | Exponential coefficient used to account for rainfall intensity on curve number      |
| Matur\_frac\_spring                   | Fraction of maturity at spring growth initiation                                    |
| CEC\_effect\_nitrification            | CEC effect on nitrification & volatilization                                        |
| N\_fixation\_limit                    | Upper Nitrogen Fixation limit                                                       |
| Biological\_mix\_efficiency           | Biological mixing efficiency                                                        |
| Soluble\_P\_exponent                  | Soluble phosphorus runoff exponent                                                  |
| Max\_depth\_bio\_mixing               | Maximum depth for biological mixing                                                 |
| OrgP\_loss\_exponent                  | Organic P loss exponent                                                             |
| MUST\_coeff                           | Coefficient in MUST EQ                                                              |
| Harg\_PET\_exponent                   | Hargreaves PET equation exponent                                                    |
| Denit\_soil\_threshold                | Denitrification soil-water threshold                                                |
| Daily\_denit\_limit                   | Upper Limit of Daily Denitrification rate                                           |
| SWAT\_delivery\_ratio\_exponent       | Exponent in Delivery Ratio for SWAT Output                                          |
| Water\_stress\_coeff                  | Water stress weighting coefficient                                                  |
| Puddling\_sat\_conduct                | Puddling Saturated conductivity                                                     |
| Groundwater\_stor\_threshold          | Groundwater storage threshold                                                       |
| Root\_temp\_stress\_exponent          | Plant root temperature stress exponent                                              |
| SCS\_index\_coeff                     | SCS curve number index coefficient                                                  |
| Plow\_depth                           | Plow layer depth                                                                    |
| CN\_retention\_param                  | Upper Limit of Curve Number Retention Parameter                                     |
| sediment\_rout\_travel\_coeff         | Sediment routing travel time coefficient                                            |
| RUSLE\_c\_factor\_res                 | RUSLE C-factor coefficient                                                          |
| RUSLE\_c\_factor\_height              | RUSLE C-factor coefficient                                                          |
| Climate\_stress\_factor               | Adjusts climatic stress factor                                                      |
| Max\_rain\_intercept                  | Maximum rainfall interception by plant canopy                                       |
| Rain\_intercept\_coeff                | Rainfall interception coefficient                                                   |
| Water\_stor\_residue\_coeff           | Water stored in litter (residue) coefficient                                        |
| Tillage\_residue\_decay\_rate\_coeff  | Exponential coefficient in EQUATION expressing tillage effect on residue decay rate |
| Microbial\_soil\_depth\_coeff         | Coefficient in oxygen EQUATION used in modifying microbial activity with soil depth |
| N\_enrich\_coeff                      | N enrichment ratio coefficient for routing                                          |
| N\_enrich\_rout\_exponent             | N enrichment ratio exponent for routing                                             |
| Fraction\_destroyed\_burn             | Fraction destroyed by burn operation                                                |
| P\_enrich\_rout\_coeff                | P enrichment ratio coefficient for routing                                          |
| P\_enrich\_rout\_exponent             | P enrichment ratio exponent for routing                                             |
| P\_move\_evap\_coeff                  | P upward movement by evaporation coefficient                                        |
| Max\_days\_grazed\_rotation           | Maximum number of days a pasture is grazed before rotation                          |
| Soil\_water\_up\_flow\_limit          | Soil water Upward Flow Limit                                                        |
| Manure\_erosion\_equation\_coeff      | Manure erosion equation coefficient                                                 |
| N\_enrich\_ratio\_delivery            | N Enrichment Ratio for Delivery to SWAT                                             |
| Dust\_distribution\_coeff             | Dust distribution coefficient                                                       |
| RUSLE2\_trans\_capacity               | RUSLE2 transport capacity parameter                                                 |
| RUSLE2\_trans\_capacity\_threshold    | RUSLE2 threshold transport capacity coefficient                                     |
| Dust\_distribution\_exponent          | Dust distribution dispersion exponent                                               |
| Manure\_erosion\_exponent             | Manure erosion exponent                                                             |
| Microbial\_top\_soil\_coeff           | Coefficient adjusts microbial activity function in the top soil layer               |
| Microbial\_decay\_coeff               | Microbial decay rate coefficient                                                    |
| Manure\_erosion\_coeff                | Manure erosion coefficient                                                          |
| Volt\_nitrification\_partition\_coeff | Volatilization/nitrification partitioning coefficient                               |
| Hydrograph\_dev\_param                | Hydrograph development parameter                                                    |
| Partition\_N\_flow\_groundwater       | Partitions Nitrogen flow from groundwater                                           |
| P\_enrich\_ratio\_deliver\_SWAT       | P Enrichment Ratio for Delivery to SWAT                                             |
| Stand\_dead\_fall\_rate\_coeff        | Standing Dead fall rate coefficient                                                 |
| Runoff\_2\_delay\_pest                | Runoff amount to delay pest application                                             |
| Soil\_water\_2\_delay\_tillage        | Soil water value to delay tillage                                                   |
| Auto\_mov\_lower\_limit               | Auto mow lower limit                                                                |
| Nitrification\_vol\_upper\_limit      | Upper Limit of Nitrification-Volatilization                                         |
| Tech\_coeff                           | Technology Coefficient                                                              |
| Drainage\_lateral\_conduct            | Estimates drainage system lateral hydraulic conductivity                            |
| P\_flux\_labile\_active\_coeff        | Coefficient regulating P flux between labile and active pool                        |
| P\_flux\_active\_stable\_coeff        | Coefficient regulating P flux between active and stable pool                        |
| N\_salt\_evap\_coeff                  | Nitrogen and Salt Upward movement by evaporation coefficient                        |
| Water\_table\_recession\_coeff        | Water table recession coefficient                                                   |
| Water\_table\_move\_limit             | Limits daily water table movement                                                   |
| Water\_table\_recession\_exponent     | Water table recession                                                               |
| Subsurface\_flow\_factor              | Subsurface flow factor                                                              |
| Flood\_evap\_limit                    | Flood Evaporation Limit                                                             |
| Runoff\_adj\_link                     | Runoff Volume Adjustment for Direct Link                                            |
| Water\_erosion\_threshold             | Water Erosion Threshold                                                             |
| Wind\_erosion\_threshold              | Wind Erosion Threshold                                                              |
| Crop\_stress\_temp\_exponent          | Exponent of Crop Stress Temperature function                                        |
| Soluble\_P\_leach\_KD                 | Soluble Phosphorus Leaching KD value                                                |
| Unknown1                              | ---                                                                                 |
| Unknown2                              | ---                                                                                 |
| Unknown3                              | ---                                                                                 |
| Irrigation\_cost                      | Cost of Irrigation Water                                                            |
| Lime\_cost                            | Cost of Lime                                                                        |
| Fuel\_cost                            | Cost of Fuel                                                                        |
| Labor\_cost                           | Cost of Labor                                                                       |
| Unknown4                              | ---                                                                                 |

</br>

| Parameters | Description                                                                                              |
|:-----------|:---------------------------------------------------------------------------------------------------------|
| RFN        | Average concentration of nitrogen in rainfall                                                            |
| CO2        | Carbon dioxide concentration in atmosphere                                                               |
| CQN        | Concentration of NO3-N in irrigation water in ppm                                                        |
| PSTX       | Pest damage scaling factor                                                                               |
| YWI        | Number years of maximum monthly 0.5 hour rainfall available                                              |
| BTA        | COEF (0-1) governing wet-dry probabilities given days of rain                                            |
| EXPK       | Parameter used to modify exponential rainfall amount distribution                                        |
| QG         | Channel Capacity Flow Rate                                                                               |
| QCF        | Exponent in watershed area flow rate equation                                                            |
| CHSO       | Average upland slope (m/m) in watershed                                                                  |
| BWD        | Channel bottom width/depth in m/m; Channel flow rate (QG) &gt; 0                                         |
| FCW        | Floodplain width/channel width in m/m                                                                    |
| FPSC       | Floodplain saturated hydraulic conductivity in mm/h                                                      |
| GWSO       | Maximum ground water storage in mm                                                                       |
| RFTO       | Ground water residence time in days                                                                      |
| RFPO       | Return Flow / (Return Flow + Deep Percolation)                                                           |
| SATO       | Saturated Conductivity adjustment factor                                                                 |
| FL         | Field length (if wind erosion is to be considered) in kilometers                                         |
| FW         | Field width (if wind erosion is to be considered) in kilometers                                          |
| ANG        | Clockwise angle of field length from north (if wind erosion is to be considered)                         |
| UXP        | Power Parameter of Modified Exponential Distribution of Wind Speed (if wind erosion is to be considered) |
| DIAM       | Soil Particle Diameter(if wind erosion is to be considered)                                              |
| ACW        | Wind Erosion Adjustment Factor                                                                           |
| GZL0       | Grazing limit                                                                                            |
| RTN0       | Number of years of cultivation at start of simulation                                                    |
| BXCT       | Linear coefficient of change in rainfall from east to west (PI/PO/KM)                                    |
| BYCT       | Linear coefficient of change in rainfall from south to north (PI/PO/KM)                                  |
| DTHY       | Time interval for flood routing (hours)                                                                  |
| QTH        | Routing Threshold (mm) â<U+0080><U+0093> VSC routing used when QVOL&gt;QTH                               |
| STND       | VSC Routing used when reach storage &gt; STND                                                            |
| DRV        | Equation for Water Erosion                                                                               |
| PCO0       | Fraction of subareas controlled by ponds                                                                 |
| RCC0       | USLE Crop Management Channel Factor                                                                      |
| CSLT       | Salt Concentration in Irrigation Water                                                                   |
| BUS1       | Exponents of MUSI equation                                                                               |
| BUS2       | Exponents of MUSI equation                                                                               |
| BUS3       | Exponents of MUSI equation                                                                               |
