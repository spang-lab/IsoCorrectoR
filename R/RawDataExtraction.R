
#' @importFrom magrittr '%>%'

RawDataExtraction <- function(data, MoleculeArray, logEnvironment, verbose) {
    
    if(verbose){message(date(), " :: processing raw data ...")}
    
    for (i in 2:ncol(data)) {
        data[, i] <- as.numeric(data[, i])
    }
    
    # remove rows that are ALL NA (=> don't remove individual NAs !)
    data <- data[!apply(apply(data, 1, is.na), 2, sum) == ncol(data), ]
    # remove columns that are all NA
    data <- data[, !apply(apply(data, 2, is.na), 2, sum) == nrow(data)]
    
    rownames(data) <- as.character(data[, "Measurements/Samples"])
    data <- as.matrix(data[, -which(colnames(data) %in% "Measurements/Samples"), drop = FALSE])
    txt_extracted <- rownames(data)
    
    # Extraction of measurement locations from the imported data to check for presence.
    
    MoleculesTotal <- length(MoleculeArray)
    MoleculeName <- names(MoleculeArray)
    
    MoleculeLocationList <- list()
    NumberTransitionsList <- list()
    TransitionLocationExpectedList <- list()
    TransitionLocationList <- list()
    MissingMoleculesList <- list()
    MissingTransitionsList <- list()
    
    for (MoleculeNo in seq_len(MoleculesTotal)) {
        
        MoleculeData <- MoleculeArray[[MoleculeNo]]  # entire 'MoleculeInfo' data for current molecule
        
        TransitionsExpected <- MoleculeData[["TransitionsExpected"]]  # from CalculateTransitions()
        NumberTransitions <- nrow(TransitionsExpected)
        
        # Check: Do the measurement labels in the first column of the data match the theoretically expected measurement labels?
        
        MissingTransitions.vec <- vector()
        TransitionLocationExpected.vec <- vector()
        
        for (TransitionNo in seq_len(NumberTransitions)) {
            
            tmpIdx <- stringr::str_detect(txt_extracted, paste0("^", rownames(TransitionsExpected)[TransitionNo], "$")) %>% which
            
            MissingTransitions.vec[TransitionNo] <- 0
            
            if (length(tmpIdx) == 0) {
                TransitionLocationExpected.vec[TransitionNo] <- 0
                MissingTransitions.vec[TransitionNo] <- 1
                
                notification <- stringr::str_c("In measurement data file: The expected measurement ID '", rownames(TransitionsExpected)[TransitionNo], "' could not be found in the measurement data file.\nCorrection will be performed, however the results may be less accurate. Please check your input file for typos.", 
                  "\nBe especially careful when considering fraction and mean enrichment values from molecules with missing measurements.")
                errorHandler(notification, logEnvironment, "warning", verbose=verbose)
                
            } else {
                
                TransitionLocationExpected.vec[TransitionNo] <- tmpIdx
                
            }
            
        }  #TransitionNo
        
        if (sum(TransitionLocationExpected.vec) == 0) {
            
            notification <- stringr::str_c("In measurement data file: None of the expected measurement IDs found for molecule '", names(MoleculeArray[MoleculeNo]), 
                "'.")
            errorHandler(notification, logEnvironment, "error", verbose=verbose)
            
        }
        
        names(TransitionLocationExpected.vec) <- rownames(TransitionsExpected)
        names(MissingTransitions.vec) <- rownames(TransitionsExpected)
        
        # If measurement IDs are missing, a new list of Measurement IDs is produced containing only those that were actually found in the input data.  This list
        # is then used for further computations.
        
        Transitions <- TransitionsExpected[MissingTransitions.vec == 0, , drop = FALSE]
        TransitionLocation <- TransitionLocationExpected.vec[MissingTransitions.vec == 0]
        
        TransitionLocationExpectedList[[MoleculeNo]] <- TransitionLocationExpected.vec
        TransitionLocationList[[MoleculeNo]] <- TransitionLocation
        
        MoleculeArray[[MoleculeNo]][["Transitions"]] <- Transitions
        MoleculeArray[[MoleculeNo]][["MissingTransitions"]] <- MissingTransitions.vec
        MoleculeArray[[MoleculeNo]][["TransitionLocationExpected"]] <- TransitionLocationExpected.vec
        MoleculeArray[[MoleculeNo]][["TransitionLocation"]] <- TransitionLocation
        
    }  #MoleculeNo
    
    if(verbose){message(date(), " :: processing raw data [OK]\n")}
    return(list(MoleculeInfo = MoleculeArray, dataRaw = data))
    
}  #RawDataExtraction
