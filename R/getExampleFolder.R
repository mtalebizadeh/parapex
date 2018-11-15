
#' getExampleFolder
#' @description Creates a copy of the Example forlder.
#' @param folderName Folder name containing example project.
#'
#' @return void
#' @export
#'
#' @examples getExampleFolder()
getExampleFolder <- function(folderName="Example") {

  if (file.exists(folderName)) {
    stop(paste(folderName),"The folder already exists! Please select a different name for example folder.")}

  else {
    folder_location <- system.file("extdata","Exampleee",package = "APEXSENSUN")
    if (file.exists("Exampleee")) {stop("Can not copy the example folder!")}
    file.copy(from = folder_location,to = getwd(),overwrite = T,recursive = T)
    file.rename(from = "Exampleee",folderName)
  }

}
