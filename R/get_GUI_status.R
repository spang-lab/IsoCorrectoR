#' Function to get activation status of IsoCorrectoR GUI
#' 
#' The function \code{getGUIstatus} can be used in scripts to get the activation status of the IsoCorrectoR GUI (logical, FALSE for closed, TRUE for active). 
#' It is required if the GUI is to be started e.g. via a bash script without manually starting an R session.
#' 
#' @export
#' 
#' @examples 
#' 
#' # this function is never called directly but BiocCheck() requires a runnable example.
#' 
#' status<-getGUIstatus()
#' 
#' @return Returns the activation status of the IsoCorrectoR GUI. FALSE represents closed and TRUE represents active state. Returns NULL if GUI has not been started.
getGUIstatus <- function() {
    baseEnvRef <- baseenv()
    return(baseEnvRef$continueIsoCorrection)
}
