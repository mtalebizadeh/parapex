#' SA4APEX
#' @description Calculates sensitivty indices.
#' @param GSA_Object An object containing analysis settings (i.e., the object
#'  returning from a call to MC4APEX)
#' @param performanceMatrix An object containing performance matrix (i.e. the
#'  object returning from a call to perfMatrix)
#'
#' @return A list object contanining sensitivity indices
#' @export
#'
#' @examples SA <- SA4APEX(GSA_Object = GSA,performanceMatrix = perfMat)

SA4APEX <- function(GSA_Object,performanceMatrix) {

  SA_Coeff_Matrix_main <- array(NaN, c(ncol(GSA_Object$GSA_Object$X),
                                      ncol(performanceMatrix)))
  SA_Coeff_Matrix_total <- array(NaN, c(ncol(GSA_Object$GSA_Object$X),
                                       ncol(performanceMatrix)))
  row.names(SA_Coeff_Matrix_main) <- names(GSA_Object$GSA_Object$X)
  colnames(SA_Coeff_Matrix_main) <- names(performanceMatrix)
  row.names(SA_Coeff_Matrix_total) <- names(GSA_Object$GSA_Object$X)
  colnames(SA_Coeff_Matrix_total) <- names(performanceMatrix)
  if (GSA_Object$Input$GSA_Type %in% c("SRC","SRRC")) {
    j = 1
    while(j<=ncol(SA_Coeff_Matrix_main)) {
      vector_SRC <- gsa_analysis_dync(GSA_Object = GSA_Object$GSA_Object,
                                      Y = performanceMatrix[,j],
                                      GSA_Type = GSA_Object$Input$GSA_Type)
      #First column for main effect
      #Second column for total effect
      SA_Coeff_Matrix_main[,j] <- vector_SRC[[1]]
      j = j + 1
    }

  } else if (GSA_Object$Input$GSA_Type %in% c("MORRIS")) {
    j = 1
    while(j <= ncol(SA_Coeff_Matrix_main)) {
      vector_SRC <- gsa_analysis_dync(GSA_Object = GSA_Object$GSA_Object,
                                      Y = performanceMatrix[,j],
                                      GSA_Type = GSA_Object$Input$GSA_Type)
      #First column for main effect
      #Second column for total effect
      SA_Coeff_Matrix_main[,j] <- vector_SRC[[2]]
      SA_Coeff_Matrix_total[,j] <- vector_SRC[[3]]
      j = j + 1
    }

  } else {
    j = 1
    while (j <= ncol(SA_Coeff_Matrix_main)) {
    vector_SRC <- gsa_analysis_dync(GSA_Object = GSA_Object$GSA_Object,
                                    Y = performanceMatrix[,j],
                                    GSA_Type = GSA_Object$Input$GSA_Type)
    #First column for main effect
    #Second column for total effect
    SA_Coeff_Matrix_main[,j] <- vector_SRC[[1]]
    SA_Coeff_Matrix_total[,j] <- vector_SRC[[2]]
    j = j + 1
    }
  }
  SA_Coeff_Matrix <- list(mainEffect = SA_Coeff_Matrix_main,
                          TotalEffect = SA_Coeff_Matrix_total)
   if (GSA_Object$Input$GSA_Type %in% c("MORRIS")) {
     names(SA_Coeff_Matrix) <- c("MuStar","Sig")
   } else if (GSA_Object$Input$GSA_Type %in% c("SRC")) {
     SA_Coeff_Matrix <- list(SRC=SA_Coeff_Matrix_main)
   }
   else if (GSA_Object$Input$GSA_Type %in% c("SRRC")) {
     SA_Coeff_Matrix <- list(SRRC=SA_Coeff_Matrix_main)
   }
  SA_Coeff_Matrix
}
#
#####
gsa_analysis_dync <- function(GSA_Object,Y,GSA_Type) {
  # Calculates sensitivity indices.

  # Args:
  #  GSA_object: An object containing sensitivity analysis settings.
  #  Y: A vector containing model outputs.
  #  GSA_Type: A character indicating sensitivity analysis method,
  #            possible values are: "MORRIS", "SRC", "SRRC", "FAST99".
  #
  # Retunrs:
  #  A data.frame containing sensitivity indices.

  if (GSA_Type == "MORRIS" ) {
    sensitivity::tell(GSA_Object ,y = Y)
    Mu_Sigma <- print(GSA_Object)
    SA_Coeff <- Mu_Sigma
  } else if (GSA_Type == "SRC") {
    GSA_Object <- src(GSA_Object$X,y=Y, rank = F,
                      nboot = 0, conf = 0.95)
    SRC_Coeff <- print(GSA_Object)
    SA_Coeff <- SRC_Coeff
  } else if (GSA_Type == "SRRC") {
    GSA_Object <- src(GSA_Object$X,y=Y, rank = TRUE,
                      nboot = 0, conf = 0.95)
    SRC_Coeff <- print(GSA_Object)
    SA_Coeff <- SRC_Coeff

  } else if (GSA_Type=="FAST99") {
    GSA_Object$X <- scale(x = GSA_Object$X, center = TRUE,
                          scale = TRUE)
    Y = scale(x = Y,center = TRUE,scale = TRUE)
    GSA_Object <- sensitivity::tell(GSA_Object, y = Y)
    FAST99_Idx <- as.data.frame(sensitivity:::print.fast99(GSA_Object))
    colnames(FAST99_Idx) <- c("Main Efffect","Total Effect")
    SA_Coeff <- FAST99_Idx
  } else {
    stop("Please Enter a Valid GSA Method")
  }
  return(SA_Coeff)
}
