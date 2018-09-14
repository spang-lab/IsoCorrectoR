# Extraction of molecule parameters

#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stringr str_extract
#' @importFrom stringr str_extract_all
#' @importFrom magrittr '%>%'
#'  
MoleculeInfoExtraction <- function(MoleculeData, ElementInfo, UltraHighRes, CorrectTracerImpurity, MoleculesFound, logEnvironment, verbose) {
  
  if(verbose){message(date(), " :: processing molecule file ...")}
  
  #Find molecules from measurement file in molecule data
  
  rownames(MoleculeData) <- as.character(MoleculeData[, 1])
  
  MoleculeLocation.vec <- vector()
  MoleculeLocationLabel.vec <- vector()
  
  for (MoleculeFound in MoleculesFound) {
    
    loc <- grep(paste0("^", MoleculeFound, "$"), rownames(MoleculeData))
    
    if (length(loc) == 1) {
      MoleculeLocation.vec <- c(MoleculeLocation.vec, loc)  # simply add current MoleculeName to vector of all valid MoleculeNames
      MoleculeLocationLabel.vec <- c(MoleculeLocationLabel.vec, MoleculeFound)
    } else {
      notification <- stringr::str_c("Molecule '", MoleculeFound, "' from the measurement data was not found in the molecule file.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
    }
  }  #MoleculeFound
  
  names(MoleculeLocation.vec) <- MoleculeLocationLabel.vec
  
  MoleculeList <- list()  # big list containing all necessary information on a molecule
  
  # MoleculeLocation.vec contains location indices of found molecules !
  
  #Analyze imported molecule formulas with regex
  
  for (MoleculeNo in seq_len(length(MoleculeLocation.vec))) {
    
    MoleculeLocation <- MoleculeLocation.vec[MoleculeNo]
    
    MoleculeInformation <- MoleculeData[MoleculeLocation,]
    
    if(verbose){message(date(), " :: :: found molecule: #", MoleculeNo, " [", names(MoleculeLocation.vec)[MoleculeNo], "]")}
    # Check if MoleculeFile contains entries in both the product ion/MS ion 
    # and the neutral loss column or only in the product ion/MS ion column.  
    # This yields NumberFragments for a given molecule and determines whether 
    # MS or MS/MS correction will be applied in the following.
    
    FragmentList <- ParseMoleculeInformation(MoleculeInformation=MoleculeInformation,
                                             ElementInfo=ElementInfo,UltraHighRes=UltraHighRes,
                                             verbose=verbose)
    
    MoleculeList[[MoleculeNo]] <- FragmentList
    
  }  # MoleculeNo
  
  names(MoleculeList) <- names(MoleculeLocation.vec)
  
  #Check logic of extracted molecule information
  
  checkMoleculeDataLogic(MoleculeList, UltraHighRes, CorrectTracerImpurity, ElementInfo, logEnvironment, verbose=verbose)
  
  #Provide log information on molecules and tracers
  
  logEnvironment$param$molecules <- names(MoleculeList)
  
  tracerLogInfo <- list()
  tracerNames <- character()
  
  #Get tracerNames, a vector of all tracer elements used
  
  for (molecule in names(MoleculeList)) {
    
    for (fragmentNo in seq_len(sum(stringr::str_detect(names(MoleculeList[[molecule]]), "Fragment")))) {
      
      fragment <- paste0("Fragment_", fragmentNo)
      
      tracers <- MoleculeList[[molecule]][[fragment]]$Tracer
      
      if (length(tracers) > 0) {
        
        tracerNames <- c(tracerNames, names(tracers))
        
      }
      
    }
    
  }
  
  tracerNames <- unique(tracerNames)
  
  for (tracer in tracerNames) {
    
    if (CorrectTracerImpurity == TRUE) {
      
      tracerLogInfo[[tracer]] <- ElementInfo[[tracer]]$Tracer.purity
      
    } else {
      
      tracerLogInfo[[tracer]] <- NA
      
    }
    
  }
  
  logEnvironment$param$tracers <- tracerLogInfo
  
  if(verbose){message(date(), " :: processing molecule file [OK]\n")}
  
  return(MoleculeInfo = MoleculeList)
}  # function()