write_param_cont <-function(uncertain_cont_param,
                            APEXCONT.dat,
                            Back_Up_APEXCONT.dat) {
  #Writes parameter values inside APEXCONT file.
  #
  # Args:
  #   uncertain_watershed_param: dataframe containing APEX basin parameters.
  #   param_file.dat: APEXCONT file name.
  #   Back_Up_param_file.dat: APEXCONT file name containing default APEX parameters.
  #
  # Returns:
  #   void (writes sampled parameter values inside the APEXCONT file).
  #
  conn_source <- file(description = Back_Up_APEXCONT.dat,open = "r+",
                      blocking = TRUE)
  conn_target <- file(description = APEXCONT.dat,open = "r+",
                      blocking = TRUE )
  #Reading from back up file into a readLines object...
  readlines_obj <- readLines(con = conn_source,n = -1)
  i = 1
  col_numbers <-ncol(uncertain_cont_param)
  while(i <= col_numbers) {
    readlines_obj <- cont2replace(readlines_obj, uncertain_cont_param[i])
    i = i + 1
  }
  ###Writing the final readLine object...
  writeLines(text = readlines_obj,con = conn_target)
  close(conn_source)
  close(conn_target)
}
#
###Sub-functionn
  cont2loc <- function(APEXCONT) {
    switch(names(APEXCONT),

      #Line 3...
           "RFN"= return(c(3,1,8,"%8.2f")),
           "CO2"= return(c(3,9,16,"%8.2f")),
           "CQN"= return(c(3,17,24,"%8.2f")),
           "PSTX"= return(c(3,25,32,"%8.2f")),
           "YWI"= return(c(3,33,40,"%8.2f")),
           "BTA"= return(c(3,41,48,"%8.2f")),
           "EXPK"= return(c(3,49,56,"%8.2f")),
           "QG"= return(c(3,57,64,"%8.2f")),
           "QCF"= return(c(3,65,72,"%8.2f")),
           "CHSO"= return(c(3,73,80,"%8.2f")),

       #Line 4...
            "BWD"=return(c(4,1,8,"%8.2f")),
            "FCW"=return(c(4,9,16,"%8.2f")),
            "FPSC"=return(c(4,17,24,"%8.2f")),
            "GWSO"=return(c(4,25,32,"%8.2f")),
            "RFTO"=return(c(4,33,40,"%8.2f")),
            "RFPO"=return(c(4,41,48,"%8.2f")),
            "SATO"=return(c(4,49,56,"%8.2f")),
            "FL"=return(c(4,57,64,"%8.2f")),
            "FW"=return(c(4,65,72,"%8.2f")),
            "ANG"=return(c(4,73,80,"%8.2f")),

      #Line 5...
          "UXP"= return(c(5,1,8,"%8.2f")),
          "DIAM"= return(c(5,9,16,"%8.2f")),
          "ACW"= return(c(5,17,24,"%8.2f")),
          "GZL0"= return(c(5,25,32,"%8.2f")),
          "RTN0"=return(c(5,33,40,"%8.2f")),
          "BXCT"= return(c(5,41,48,"%8.2f")),
          "BYCT"=return(c(5,49,56,"%8.2f")),
          "DTHY"=return(c(5,57,64,"%8.2f")),
          "QTH"=return(c(5,65,72,"%8.2f")),
          "STND"=return(c(5,73,80,"%8.2f")),

      #Line 6...
          "DRV"=return(c(6,1,8,"%8.2f")),
          "PCO0"=return(c(6,9,16,"%8.2f")),
          "RCC0"=return(c(6,17,24,"%8.2f")),
          "CSLT"=return(c(6,25,32,"%8.2f")),
          "BUS1"=return(c(6,33,40,"%8.2f")),
          "BUS2"=return(c(6,41,48,"%8.2f")),
          "BUS3"=return(c(6,49,56,"%8.2f"))

    )
  }

####Sub-Function....
    cont2replace <- function(readlines_obj,APEXCONT) {
      loc_vec <- cont2loc(APEXCONT)
      substr(x = readlines_obj[strtoi(loc_vec[1])], start = strtoi(loc_vec[2]),
             stop = strtoi(loc_vec[3])) <- sprintf(fmt = loc_vec[4],APEXCONT[[1]])
      return(readlines_obj)

    }





