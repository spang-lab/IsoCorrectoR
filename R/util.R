# Supporting functions

#' @importFrom methods is
#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom readxl read_excel
#' @importFrom readr read_csv
#' @importFrom tools file_ext

#' 
# Function to check parameters of IsoCorrection() function call

checkIsoCorrectionParameters <- function(MeasurementFile, ElementFile, MoleculeFile, CorrectTracerImpurity, CorrectTracerElementCore, CalculateMeanEnrichment, UltraHighRes, FileOut, FileOutFormat, ReturnResultsObject, CorrectAlsoMonoisotopic, CalculationThreshold, CalculationThreshold_UHR, logEnvironment, Testmode, verbose) {
  
  # Check filenames
  
  filenames <- list()
  
  filenames[["MoleculeFile"]] <- MoleculeFile
  filenames[["ElementFile"]] <- ElementFile
  filenames[["MeasurementFile"]] <- MeasurementFile
  filenames[["FileOut"]] <- FileOut
  
  for (name in names(filenames)) {
    
    if (is.na(filenames[[name]]) || filenames[[name]] == "") {
      
      notification <- paste0("No or unsuitable value provided for parameter ", name, ".")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  }
  
  # Check boolean values
  
  booleans <- list()
  
  booleans[["CorrectTracerImpurity"]] <- CorrectTracerImpurity
  booleans[["CorrectTracerElementCore"]] <- CorrectTracerElementCore
  booleans[["CorrectAlsoMonoisotopic"]] <- CorrectAlsoMonoisotopic
  booleans[["UltraHighRes"]] <- UltraHighRes
  booleans[["CalculateMeanEnrichment"]] <- CalculateMeanEnrichment
  booleans[["ReturnResultsObject"]] <- ReturnResultsObject
  booleans[["Testmode"]] <- Testmode
  
  for (name in names(booleans)) {
    
    if (is.na(booleans[[name]]) || (!(booleans[[name]] == TRUE) && !(booleans[[name]] == FALSE) && !(booleans[[name]] == 0) && !(booleans[[name]] == 1))) {
      
      notification <- paste0("No or unsuitable value provided for parameter ", name, ".")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  }
  
  # Check format of output file
  isValid<-FileOutFormat%in%c("xls","csv")
  if(!isValid){
    notification <- paste0("The output format for the result file should be either \"xls\" or \"csv\". Yours was \"", FileOutFormat, "\".")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
  }
  
  # Check CalculationThreshold values
  
  if (!UltraHighRes) {
    
    if (CalculationThreshold > (10^-2) || CalculationThreshold < 0) {
      
      notification <- paste0("Input parameter 'CalculationThreshold' must be between 0 and 0.01 in normal resolution correction. If unsure, use the default value of 10^-8.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  } else if (UltraHighRes) {
    
    if (CalculationThreshold_UHR < 0 || (CalculationThreshold_UHR != round(CalculationThreshold_UHR))) {
      
      notification <- paste0("Input parameter 'CalculationThreshold_UHR' must be integer and non-negative in ultra high resolution correction. If unsure, use the default value of 8.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  }
  
}

# Function to write paths/filenames, correction parameters and errors/warnings to a log file

writeLog <- function(logEnvironment) {
  
  logFile <- file.path(logEnvironment$param$OutputDirectory, logEnvironment$param$LogfileName)
  
  disclaimer <- "IsoCorrectoR is free software and comes without any warranty."
  
  write(paste0("ISOCORRECTOR LOG-FILE - ", logEnvironment$param$Timestamp, "\nIsoCorrectoR-version: ", logEnvironment$param$version, "\n\n", disclaimer, 
               "\n"), file = logFile)
  
  # ERRORS AND WARNINGS
  
  if (length(logEnvironment$error) > 0) {
    
    write(paste0("AN ERROR HAS OCCURED, THE CORRECTION PROCESS WAS ABORTED. ERROR:\n\n", paste(logEnvironment$error, collapse = "\nError in function: "), 
                 "\n"), file = logFile, append = TRUE)
    
  }
  
  if (length(logEnvironment$warning) > 0) {
    
    write(paste0("WARNINGS HAVE OCCURRED, THE CORRECTION MAY BE FAULTY. WARNING(s):\n\n", paste(logEnvironment$warning, collapse = "\n\n"), "\n"), file = logFile, 
          append = TRUE)
    
  }
  
  # FILES AND DIRECTORIES
  
  write(paste0("FILES AND DIRECTORIES\n\nMeasurement File: ", logEnvironment$param$MeasurementFile, "\nMolecule File: ", logEnvironment$param$MoleculeFile, 
               "\nElement File: ", logEnvironment$param$ElementFile, "\nWorking Directory: ", logEnvironment$param$WorkingDirectory, "\nOutput Directory: ", logEnvironment$param$OutputDirectory, 
               "\nResult File(s): ", logEnvironment$param$FileOut, "\n"), file = logFile, append = TRUE)
  
  # PARAMETERS
  
  write(paste0("PARAMETERS\n\nCorrection for tracer impurity: ", logEnvironment$param$CorrectTracerImpurity, "\nCorrection of tracer element in core molecule: ", 
               logEnvironment$param$CorrectTracerElementCore, "\nUltra-High-Resolution correction mode: ", logEnvironment$param$UltraHighRes, "\nCalculation of mean enrichment: ", 
               logEnvironment$param$CalculateMeanEnrichment, "\nMolecules in the correction process: ", paste(logEnvironment$param$molecules, collapse = ", "), "\nTracer(s): ", 
               paste(names(logEnvironment$param$tracers), collapse = ", "), "\nTracer purity: ", paste(logEnvironment$param$tracers, collapse = ", "), "\nLimit value for calculation: ", 
               logEnvironment$param$threshold, "\n"), file = logFile, append = TRUE)
  
}

# Function that is called in the case of anticipated exceptions or warnings for a clean exit and logging of warnings

errorHandler <- function(notification, logEnvironment, type, verbose) {
  
  if (type == "warning" || type == "general") {
    
    logEnvironment[[type]] <- c(logEnvironment[[type]], notification)
    
    if (type == "warning") {
      
      warning(notification)
      
    } else if (type == "general") {
      
      if(verbose){message(date(), " ", notification)}
      
    }
  } else if (type == "error") {
    
    stop(notification)
    
  } else {
    
    notification <- paste0("Invalid error-type ", type, " supplied to errorHandler function")
    
    stop(notification)
    
  }
  
}

# Function to check the structure of measurement data

checkRawData <- function(inputFile, logEnvironment, verbose) {
  
  if(verbose){message(date(), " :: checking raw data ...")}
  
  location <- "In measurement data file: "
  
  data <- as.data.frame(readFileByExt(inputFile, verbose=verbose))
  
  if (ncol(data) < 2) {
    
    notification <- paste0(location, ncol(data), " columns detected but at least 2 are required.")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  if (nrow(data) < 2) {
    
    notification <- paste0(location, nrow(data), " rows of measurements detected but at least 2 are required.")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  if (names(data[1]) != "Measurements/Samples") {
    
    notification <- paste0(location, "Entry in column 1, row 1 must be 'Measurements/Samples'.")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  duplicateRow <- anyDuplicated(data[[1]])
  
  if (duplicateRow) {
    
    notification <- paste0(location, "Duplicate measurement identifier in row ", duplicateRow, " of column 1.")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  duplicateCol <- anyDuplicated(colnames(data))
  
  if (duplicateCol) {
    
    notification <- paste0(location, "Duplicate sample name in column ", duplicateCol, ".")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  rawDataMoleculeNames <- as.character(data[, 1])
  
  checkRes <- checkRawDataMoleculeNames(rawDataMoleculeNames, splitString = "_")
  
  MoleculesFound <- checkRes$moleculesFound
  NoMoleculeSyntax <- checkRes$invalidMolecules
  
  # there are some invalid molecules
  if (length(NoMoleculeSyntax) != 0) {
    
    notification <- stringr::str_c(location, paste0(length(NoMoleculeSyntax), " measurement identifiers in column 1 of the measurement data did not show proper syntax (name of molecule followed by _): ", 
                                                    paste0(NoMoleculeSyntax, collapse = ", ")))
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
  }
  
  for (col in 2:ncol(data)) {
    
    #if (is.element(FALSE, sapply(data[[col]], function(x) (is.numeric(x) || is.na(x))))) {
    if (is.element(FALSE, vapply(data[[col]], function(x) (is.numeric(x) || is.na(x)),TRUE))) {
      notification <- stringr::str_c(location, "Invalid value (not numeric and not NA) in column ", col, " of measurement file.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  }
  
  MoleculesFound <- unique(MoleculesFound)
  
  if(verbose){message(date(), " :: checking raw data [OK]\n")}
  
  return(list(RawData = data, MoleculesFound = MoleculesFound))
  
}

checkRawDataMoleculeNames <- function(moleculeNames, splitString) {
  tmp.split <- strsplit(moleculeNames, split = splitString)
  wrongSyntax.index <- which(lapply(tmp.split, function(x) length(x) == 1) == TRUE)  #no split
  if (length(wrongSyntax.index) != 0) {
    tmp.split <- tmp.split[-wrongSyntax.index]
  }
  
  return(list(moleculesFound = unlist(lapply(lapply(tmp.split, function(x) x[seq_len(length(x) - 1)]), function(x) if (length(x) == 1) {
    x
  } else {
    paste0(x, collapse = splitString)
  })), invalidMolecules = moleculeNames[wrongSyntax.index]))
}

# Function to check the structure of the supplied molecule file

checkMoleculeDataStructure <- function(data, logEnvironment, verbose) {
  if(verbose){message(date(), " :: checking molecule data structure ...")}
  
  # Check number of columns
  
  if (ncol(data) != 3) {
    
    notification <- paste0(ncol(data), " columns detected in Molecule File but 3 are required.")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  # Check column names
  
  colnamesMoleculeFileExp <- c("molecule", "ms ion or ms/ms product ion", "ms/ms neutral loss")
  colnamesMoleculeFile <- tolower(colnames(data))
  
  for (col in seq_len(ncol(data))) {
    
    if (colnamesMoleculeFileExp[col] != colnamesMoleculeFile[col]) {
      
      notification <- paste0("Column names '", paste0(colnamesMoleculeFile, collapse = ", "), "' found in the molecule file do not match those expected for the file: '", 
                             paste0(colnamesMoleculeFileExp, collapse = ", "), "'.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
  }
  
  # Check if there are rows
  
  if (nrow(data) == 0) {
    
    notification <- "No rows found in molecule file."
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  # Use regular expressions to check column data
  
  moleculeRegex <- "^((|Lab)[A-Z][a-z]?[1-9][0-9]?)+$"
  
  invalidSyntax <- " doesn't match the required syntax."
  
  for (row in seq_len(nrow(data))) {
    
    moleculeName <- data[[1]][[row]]
    
    # Check if Molecule names are NA
    
    if (is.na(moleculeName)) {
      
      notification <- paste0("In molecule file row ", row, ": No molecule name specified.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
    productIon <- data[[2]][[row]]
    neutralLoss <- data[[3]][[row]]
    
    substr1 <- "In molecule file column "
    substr2 <- paste0(", row ", row, ", molecule ", moleculeName, ": Molecule formula ")
    
    locationSpecifierProduct <- paste0(substr1, "2", substr2)
    locationSpecifierNeutral <- paste0(substr1, "3", substr2)
    
    if (is.na(productIon) == FALSE && is.na(neutralLoss) == FALSE) {
      
      if (stringr::str_detect(productIon, moleculeRegex) == FALSE) {
        
        notification <- paste0(locationSpecifierProduct, productIon, invalidSyntax)
        errorHandler(notification, logEnvironment, "error", verbose=verbose)
        
      } else if (stringr::str_detect(neutralLoss, moleculeRegex) == FALSE) {
        
        notification <- paste0(locationSpecifierNeutral, neutralLoss, invalidSyntax)
        errorHandler(notification, logEnvironment, "error", verbose=verbose)
        
      }
      
    } else if (is.na(productIon) == FALSE) {
      
      if (stringr::str_detect(productIon, moleculeRegex) == FALSE) {
        
        notification <- paste0(locationSpecifierProduct, productIon, invalidSyntax)
        errorHandler(notification, logEnvironment, "error", verbose=verbose)
        
      }
    } else {
      
      notification <- paste0(locationSpecifierProduct, "for ", colnamesMoleculeFileExp[2], " must be provided.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
    duplicate <- anyDuplicated(data[[1]])
    
    if (duplicate != 0) {
      
      notification <- paste0("In molecule file: Molecule ", data[[1]][duplicate], " present in duplicate.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
  }
  if(verbose){message(date(), " :: checking molecule data structure [OK]\n")}
  
  return(data.frame(data))
}

# Check if molecule data makes sense

checkMoleculeDataLogic <- function(extractedData, UltraHighRes, CorrectTracerImpurity, ElementList, logEnvironment, verbose) {
  if(verbose){message(date(), " :: :: checking molecule data logic ...")}
  molecules <- names(extractedData)
  
  numberFragments <- list()
  
  tracers <- list()
  
  for (molecule in molecules) {
    
    numberFragments[[molecule]] <- sum(stringr::str_detect(names(extractedData[[molecule]]), "Fragment"))
    
    tracers[[molecule]] <- c()
    
    tracers[[molecule]][1] <- length(extractedData[[molecule]]$Fragment_1$Tracer)
    
    if (numberFragments[[molecule]] > 1) {
      
      tracers[[molecule]][2] <- length(extractedData[[molecule]]$Fragment_2$Tracer)
      
    } else {
      
      tracers[[molecule]][2] <- 0
      
    }
    
    if (sum(tracers[[molecule]]) == 0) {
      
      notification <- paste0("In molecule file: No tracer element provided for molecule ", molecule, ".")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  }
  
  if (UltraHighRes == TRUE) {
    
    for (molecule in molecules) {
      
      if (numberFragments[[molecule]] != 1) {
        
        notification <- paste0("In molecule file: More than one fragment provided for molecule ", molecule, ", but ultra high resolution correction is not compatible with MS/MS data.")
        errorHandler(notification, logEnvironment, "error", verbose=verbose)
        
      }
      
    }
    
  } else {
    
    for (molecule in molecules) {
      
      if ((tracers[[molecule]][1] > 1) || (tracers[[molecule]][2] > 1)) {
        
        notification <- paste0("In molecule file: More than one tracer element provided for ms ion/product Ion and/or neutral loss of molecule ", molecule, ", but normal resolution correction is only compatible with a single tracer element.")
        errorHandler(notification, logEnvironment, "error", verbose=verbose)
        
      }
      
      if ((tracers[[molecule]][1] > 0) && (tracers[[molecule]][2] > 0)) {
        
        if (extractedData[[molecule]]$Fragment_1$IDTracer != extractedData[[molecule]]$Fragment_2$IDTracer) {
          
          notification <- paste0("In molecule file: Tracer elements in product ion and neutral loss of molecule ", molecule, " do not match.")
          errorHandler(notification, logEnvironment, "error", verbose=verbose)
          
        }
      }
    }
  }
  
  # Check for duplicate elements in the molecule formula
  
  for (molecule in molecules) {
    
    duplicatesPresent <- 0
    
    for (fragmentNo in seq_len(numberFragments[[molecule]])) {
      
      fragment <- paste0("Fragment_", fragmentNo)
      
      duplicatesPresent <- duplicatesPresent + anyDuplicated(names(extractedData[[molecule]][[fragment]]$Element)) + anyDuplicated(names(extractedData[[molecule]][[fragment]]$Tracer))
      
    }
    
    if (duplicatesPresent) {
      
      notification <- paste0("In molecule file: Molecule ", molecule, " contains duplicate elements.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  }
  
  elementsProvided <- names(ElementList)
  
  for (molecule in molecules) {
    
    for (fragmentNo in seq_len(numberFragments[[molecule]])) {
      
      fragment <- paste0("Fragment_", fragmentNo)
      
      elementsFragment <- extractedData[[molecule]][[fragment]]$Element
      
      if (fragmentNo == 1) {
        fragmentString = "ms ion/product ion"
      } else {
        fragmentString = "neutral loss"
      }
      
      commonString <- paste0("In molecule file: In ", fragmentString, " of molecule ", molecule, ": ")
      
      # Check if elements in the molecule formula are present in the element file
      
      if (anyNA(match(names(elementsFragment), elementsProvided))) {
          
        notification <- paste0(commonString, "One or more elements in the provided formula were not found in the element file.")
        errorHandler(notification, logEnvironment, "error", verbose=verbose)
          
      }
  
      # Do this only if the fragment contains a tracer
      
      if (tracers[[molecule]][[fragmentNo]] > 0) {
        
        tracersFragment <- extractedData[[molecule]][[fragment]]$Tracer
        
        # Check if tracer elements in the molecule formula are also present as non-tracer elements (which gives the total number of atoms, not just the tracer
        # atoms)
        
        if (anyNA(match(names(tracersFragment), names(elementsFragment)))) {
            
          notification <- paste0(commonString, "Information on total atom count (labelled + unlabelled) missing for tracer element(s).")
          errorHandler(notification, logEnvironment, "error", verbose=verbose)
            
        }
          
        # For each tracer element found in the molecule formula, check the number of tracer atoms exceeds the total number of atoms of that element
        
        for (tracer in names(tracersFragment)) {
          
          if (tracersFragment[[tracer]] > elementsFragment[[tracer]]) {
            
            notification <- paste0(commonString, "Number of tracer isotopes of element ", tracer, " exceeds number of atoms of that element present in total.")
            errorHandler(notification, logEnvironment, "error")
            
          }
          
          # For each tracer element, check if a mass shift is provided in the element file
          
          TracerMassShift <- ElementList[[tracer]]$Tracer.isotope.mass.shift
          
          if (TracerMassShift == 0) {
            
            notification <- paste0(commonString, "No mass shift information for tracer element ", tracer, " provided in element file.")
            errorHandler(notification, logEnvironment, "error", verbose=verbose)
            
          }
          
        }
          
        # Given correction for tracer impurity is active, for each tracer element, check if purity is provided in the element file
          
        if (CorrectTracerImpurity == TRUE) {
          
          for (tracer in names(tracersFragment)) {
            
            if (ElementList[[tracer]]$Tracer.purity == 0) {
              
              notification <- paste0(commonString, "No tracer purity information for tracer element ", tracer, " provided in element file, but correction for tracer impurity is active.")
              errorHandler(notification, logEnvironment, "error", verbose=verbose)
              
            }
            
          }
          
        }
          
      } 
      
    }
    
  }
  if(verbose){message(date(), " :: :: checking molecule data logic [OK]")}
}

checkElementDataStructure <- function(data, logEnvironment, verbose) {
  if(verbose){message(date(), " :: checking element data structure ...")}
  # Check number of columns
  
  if (ncol(data) != 4) {
    
    notification <- paste0(ncol(data), " columns detected in Element File but 4 are required.")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  # Check column names
  
  colnamesElementFileExp <- c("element", "isotope abundance_mass shift", "tracer isotope mass shift", "tracer purity")
  colnamesElementFile <- tolower(colnames(data))
  
  for (col in seq_len(ncol(data))) {
    
    if (colnamesElementFileExp[col] != colnamesElementFile[col]) {
      
      notification <- paste0("Column names '", paste0(colnamesElementFile, collapse = ", "), "' found in the element file do not match those expected for the file: '", 
                             paste0(colnamesElementFileExp, collapse = ", "), "'.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  }
  
  # Check if there are rows
  
  if (nrow(data) == 0) {
    
    notification <- "No rows found in element file."
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  # Use regular expressions to check column data
  
  regexVector <- c()
  
  regexVector[1] <- "^[A-Z][a-z]?$"
  
  regexVector[2] <- "(((^1|^1\\.[0]+|^[0]\\.[0-9]+)_(|-)[0-9])(/((1|1\\.[0]+|[0]\\.[0-9]+)_(|-)[0-9]))*)$"
  
  regexVector[3] <- "^(|-)[1-9]$"
  regexVector[4] <- "^1$|^1\\.[0]+$|^[0]\\.[0-9]+$"
  
  notificationVector <- c()
  
  notificationVector[1] <- " is no valid element symbol."
  notificationVector[2] <- " doesn't show appropriate isotope abundance/mass shift definition syntax."
  notificationVector[3] <- " is no valid tracer isotope mass shift definition."
  notificationVector[4] <- " is no valid tracer purity definition."
  
  for (col in seq_len(ncol(data))) {
    
    currentCol <- data[[col]]
    
    # This conversion of NA to empty string for columns 1 and 2 is done because NA values in columns 3 and 4 must not throw an error
    
    for (row in seq_len(nrow(data))) {
      
      if (is.na(currentCol[row]) && col == (1 || 2)) {
        
        currentCol[row] <- ""
        
      }
      
      if (is.na(currentCol[row]) == FALSE) {
        
        if (stringr::str_detect(currentCol[row], regexVector[col]) == FALSE) {
          
          locationSpecifier <- paste0("In element file column ", col, ", row ", row, ": ")
          
          notification <- paste0(locationSpecifier, currentCol[row], notificationVector[col])
          
          errorHandler(notification, logEnvironment, "error", verbose=verbose)
          
        }
      }
    }
  }
  
  # Check for duplicate elements
  
  duplicate <- anyDuplicated(data[[1]])
  
  if (duplicate != 0) {
    
    notification <- paste0("In element file: Element ", data[[1]][duplicate], " present in duplicate.")
    errorHandler(notification, logEnvironment, "error", verbose=verbose)
    
  }
  
  if(verbose){message(date(), " :: checking element data structure [OK]\n")}
  
  return(data.frame(data))
}

# Check if element data makes sense

checkElementDataLogic <- function(extractedData, UltraHighRes, logEnvironment, verbose) {
  #if(verbose){message(date(), " [checkElementDataLogic]")}
  if(verbose){message(date(), " :: :: checking element data logic ...")}
  for (element in names(extractedData)) {
    
    if (any(duplicated(as.data.frame(extractedData[[element]]$Isotopes)$MassShift))) {
      notification <- paste0("In element file: Duplicated MassShift for element ", element, ".")
      errorHandler(notification, logEnvironment, "error")
    }
    
    
    # Without rounding, Si doesn't work as the sum is 1.000001. Howver, should we really round here?
    if (round(sum(data.frame(extractedData[[element]][["Isotopes"]])[, 1]), 5) > 1) {
      
      notification <- paste0("In element file: Sum of isotope abundances of element ", element, " is > 1.")
      errorHandler(notification, logEnvironment, "error", verbose=verbose)
      
    }
    
  }
  
  # In case of UHR mode, check if the mass shift of the tracers is also found in the tracer elements natural isotopes. If not, probabilities cannot be
  # calculated
  
  if (UltraHighRes) {
    
    for (element in names(extractedData)) {
      
      TracerMassShift <- extractedData[[element]]$Tracer.isotope.mass.shift
      
      if (TracerMassShift > 0 && is.na(match(TracerMassShift, as.data.frame(extractedData[[element]]$Isotopes)$MassShift))) {
        
        notification <- paste0("In element file: Mass shift '", TracerMassShift, "' of tracer element ", element, " is not found in the natural isotope information of that element.")
        errorHandler(notification, logEnvironment, "error", verbose=verbose)
        
      }
      
    }
    
  }
  
  if(verbose){message(date(), " :: :: checking element data logic [OK]")}
  
}

readFileByExt <- function(file, verbose) {

  fileExt <- tools::file_ext(file)
  if (fileExt == "") {
    stop(date(), " :: no file extension for file: ", file, "\n")
  } else {
    if(verbose){message(date(), " :: :: detected file format: ", fileExt)}
    if (fileExt %in% c("xls", "xlsx")) {
      fileData <- readxl::read_excel(file)
    } else if (fileExt == "csv") {
      fileData <- readr::read_csv(file)
    } else {
      stop(date(), " :: unsupported file format: ", fileExt, "\n")
    }
  }
  return(fileData)
}

convertEnrichmentList2matrix <- function(input, verbose) {
  if(verbose){message(date(), " :: :: converting enrichment results into a matrix ...")}
  # matrix with sampleIds as columns and mean enrichment for molecule on row
  m <- matrix(ncol = length(input), nrow = length(input[[1]]))
  colnames(m) <- names(input)
  rownames(m) <- names(input[[1]])
  
  # TODO: quick and dirty loop. find more elegant solution !!
  for (i in seq_len(length(input))) {
    for (i2 in seq_len(length(input[[i]]))) {
      m[i2, i] <- as.numeric(input[[i]][[i2]])
    }  #i2
  }  #i
  if(verbose){message(date(), " :: :: converting enrichment results into a matrix [OK]")}
  return(m)
  
}


convert2matrix <- function(input) {
  vec.transitions <- vector()
  # we need to loop over all samples to really get all transitions
  for (sample in names(input)) {
    for (molecule in names(input[[sample]])) {
      vec.transitions <- c(vec.transitions, names(input[[sample]][[molecule]][[1]]))
    }  #molecule
  }  #sample
  transitions.unique <- unique(vec.transitions)
  
  templateMat <- matrix(NA, ncol = length(names(input)), nrow = length(transitions.unique))
  colnames(templateMat) <- names(input)
  rownames(templateMat) <- transitions.unique
  
  # now populate the four matrices
  list.corrMat <- list()
  for (corrMat in names(input[[1]][[1]])) {
    
    templateMat2 <- templateMat
    for (sample in names(input)) {
      for (molecule in names(input[[sample]])) {
        for (transition in names(input[[sample]][[molecule]][[corrMat]])) {
          templateMat2[transition, sample] <- input[[sample]][[molecule]][[corrMat]][[transition]]
        }  #transition
      }  #molecule
    }  #sample
    list.corrMat[[corrMat]] <- templateMat2
  }  #corrMat
  
  return(list.corrMat)
}  #convert2matrix

defineOutputFilenames <- function(logEnvironment,FileOut,FileOutFormat,
                                  CalculateMeanEnrichment,CorrectAlsoMonoisotopic) {
  
  #DEFINE NAMES OF OUTPUT FILES
  
  logEnvironment$param$LogfileName <-
    paste0(logEnvironment$param$NamePrefix, "_", FileOut, ".log")
  
  CorrectedDataNamesBase <-
    c("Corrected",
      "CorrectedFractions",
      "Residuals",
      "RelativeResiduals")
  
  if (CorrectAlsoMonoisotopic) {
    CorrectedDataNamesBase <-
      c(
        CorrectedDataNamesBase,
        "CorrectedMonoisotopic",
        "CorrectedMonoisotopicFractions"
      )
    
  }
  
  if (CalculateMeanEnrichment) {
    CorrectedDataNames <-
      c(CorrectedDataNamesBase, "MeanEnrichment", "RawData")
    
  }
  else {
    CorrectedDataNames <- c(CorrectedDataNamesBase, "RawData")
    
  }
  
  if (FileOutFormat == "csv") {
    CorrectedDataFileNames <-
      vector(mode = "character",
             length = length(CorrectedDataNames))
    names(CorrectedDataFileNames) <- CorrectedDataNames
    
    for (name in CorrectedDataNames) {
      CorrectedDataFileNames[name] <-
        paste0(logEnvironment$param$NamePrefix,
               "_",
               FileOut,
               "_",
               name,
               ".",
               FileOutFormat)
      
    }
    
    logEnvironment$param$FileOut <-
      paste0(CorrectedDataFileNames, collapse = ", ")
    
  }
  else if (FileOutFormat == "xls") {
    logEnvironment$param$FileOut <-
      paste0(logEnvironment$param$NamePrefix,
             "_",
             FileOut,
             ".",
             FileOutFormat)
    
    CorrectedDataFileNames <- NA
    
  }
  
  return(list("CorrectedDataNames" = CorrectedDataNames,
              "CorrectedDataNamesBase" = CorrectedDataNamesBase,
              "CorrectedDataFileNames" = CorrectedDataFileNames))
  
}

#WRITE CORRECTED DATA TO FILE

writeCorrectionResultsToFile <- function(CorrectedDataOutput,NamesCorrectedOutput,logEnvironment, verbose) {

  if(verbose){message(date()," :: writing results to ",logEnvironment$param$FileOut," ...")}
  if (logEnvironment$param$FileOutFormat == "xls") {
    WriteXLS::WriteXLS(
      CorrectedDataOutput,
      file.path(logEnvironment$param$OutputDirectory, logEnvironment$param$FileOut),
      perl = "perl",
      row.names = TRUE
    )
    
  } else if (logEnvironment$param$FileOutFormat == "csv") {
    for (name in NamesCorrectedOutput$CorrectedDataNames) {
      utils::write.csv(
        CorrectedDataOutput[[name]],
        file = file.path(
          logEnvironment$param$OutputDirectory,
          NamesCorrectedOutput$CorrectedDataFileNames[name]
        )
      )
      
    }
  }
  if(verbose){message(date()," :: writing results to ",logEnvironment$param$FileOut," [OK]\n")}
}
