# Extraction of molecule parameters

#' @importFrom stringr str_c
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stringr str_extract
#' @importFrom stringr str_extract_all
#' @importFrom magrittr '%>%'
#'  
MoleculeInfoExtraction <- function(MoleculeData, ElementArray, UltraHighRes, CorrectTracerImpurity, MoleculesFound, logEnvironment) {
  
  message(date(), " [MoleculeInfoExtraction()]")
  message(date(), " :: reading MoleculeFile")
  
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
      errorHandler(notification, logEnvironment, "error")
    }
  }  #MoleculeFound
  
  names(MoleculeLocation.vec) <- MoleculeLocationLabel.vec
  
  # The strings imported from MoleculeFile containing information on the number 
  # of atoms per element are analyzed using regular expressions.  The regular
  # expression Characters extracts all characters, the regular expression 
  # Numbers extracts all numbers.  TracerElement looks for the tracer element
  # identifier.
  
  Characters <- "([A-Za-z]*)"
  Numbers <- "([0-9]+[0-9]|[0-9])"
  TracerElement <- "(Lab)"
  
  MoleculeList <- list()  # big list containing all necessary information on a molecule
  
  # MoleculeLocation.vec contains location indices of found molecules !
  
  for (MoleculeNo in seq_len(length(MoleculeLocation.vec))) {
    
    MoleculeLocation <- MoleculeLocation.vec[MoleculeNo]
    
    message(date(), " :: molecule: ", MoleculeNo, " [", names(MoleculeLocation.vec)[MoleculeNo], "]")
    # Check if MoleculeFile contains entries in both the product ion/MS ion 
    # and the neutral loss column or only in the product ion/MS ion column.  
    # This yields NumberFragments for a given molecule and determines whether 
    # MS or MS/MS correction will be applied in the following.
    
    if (is.na(MoleculeData[MoleculeLocation, 3])) {
      NumberFragments <- 1
    } else {
      NumberFragments <- 2
    }
    
    # collect features for each individual fragment
    tmpFragmentList <- list()
    
    # FragmentList contains one or two tmpFragmentList()s
    FragmentList <- list()
    
    # The following is done for each Molecule and each fragment 
    # of a given molecule
    for (Fragment in seq_len(NumberFragments)) {
      
      # For each molecule(-fragment), the string containinig the element count 
      # information is split into the individual elements using a regex
      
      MoleculeFormulaRegex <- "(|Lab)[A-Z][a-z]?[1-9][0-9]?"
      
      Elements <- stringr::str_extract_all(MoleculeData[MoleculeLocation, Fragment + 1], MoleculeFormulaRegex)[[1]]
      
      # Now each element in Elements, consisting of the element ID and the 
      # element count, is again separated into ID and element count using 
      # the regular expressions constructed above. 
      # Additionally, a check is made to identify the 
      # tracer element information ('Lab').
      
      ElementCount.vec <- vector()
      Element.vec <- vector()
      TracerCount.vec <- vector()
      Tracer.vec <- vector()
      
      NumberElementsandTracer <- length(Elements)
      
      for (Element in seq_len(NumberElementsandTracer)) {
        
        # TRUE => current Element is Tracer // FALSE => current Element is no Tracer
        if (Elements[Element] %>% stringr::str_detect(TracerElement)) {
          
          TracerCount.vec <- c(TracerCount.vec, Elements[Element] %>% stringr::str_replace(TracerElement, "") %>% stringr::str_extract(Numbers) %>% 
                                 as.numeric)
          Tracer.vec <- c(Tracer.vec, Elements[Element] %>% stringr::str_replace(TracerElement, "") %>% stringr::str_extract(Characters))
        } else {
          
          ElementCount.vec <- c(ElementCount.vec, Elements[Element] %>% stringr::str_extract(Numbers) %>% as.numeric)
          Element.vec <- c(Element.vec, Elements[Element] %>% stringr::str_extract(Characters))
        }  #else
      }  #Element    
      names(ElementCount.vec) <- Element.vec
      names(TracerCount.vec) <- Tracer.vec
      
      tmpFragmentList <- list(Element = ElementCount.vec, Tracer = TracerCount.vec)
      
      NumberElements <- length(ElementCount.vec)
      NumberTracers <- length(TracerCount.vec)
      
      # 'NumberElements' - 'NumberTracers' gives the number of non-tracer elements NumberElementsNonTracer.
      NumberElementsNonTracer <- NumberElements - NumberTracers
      
      if (!UltraHighRes) {
        
        # If a tracer is present in the molecule(-fragment) considered, 
        # the tracer parameters are extracted.  'MaxLabel' is the maximum amount
        # of tracer isotope that is expected to be found in the 
        # molecule(-fragment) due to metabolism while 'nTracerMax' is 
        # the maximum amount of tracer element (labelled or
        # unlabelled) in that same species. 'IDTracer' is the tracer elements ID.
        
        if (NumberTracers > 0) {
          
          MaxLabel <- max(TracerCount.vec)
          IDTracer <- TracerCount.vec %>% which.max %>% names
          
          # This part deals with the possibility of having a tracer element 
          # that shows isotopes with a negative mass shift.  
          # In this case IsoCombinationsMaster()
          # has to calculate `PlacesToAssign` differently for the tracer.
          
          if (sum(ElementArray[[IDTracer]][["Isotopes"]][[1]][["MassShift"]] < 0) > 0) {
            NegIsoTracer <- 1  # => TRUE
          } else {
            NegIsoTracer <- 2  # => FALSE
          }
          
          names(NegIsoTracer) <- IDTracer
          
          tmp.which <- which(names(ElementCount.vec) %in% names(TracerCount.vec))
          nTracerMax <- as.numeric(ElementCount.vec[tmp.which])
          
          if (length(nTracerMax) == NumberTracers) {
            
            names(nTracerMax) <- names(TracerCount.vec)
            
          }
          
          names(MaxLabel) <- IDTracer
          
          tmpFragmentList[["MaxLabel"]] <- MaxLabel
          tmpFragmentList[["IDTracer"]] <- IDTracer
          tmpFragmentList[["nTracerMax"]] <- nTracerMax
          tmpFragmentList[["NegIsoTracer"]] <- NegIsoTracer
        } else {
          tmpFragmentList[["MaxLabel"]] <- NA
          tmpFragmentList[["IDTracer"]] <- NA
          tmpFragmentList[["nTracerMax"]] <- NA
          tmpFragmentList[["NegIsoTracer"]] <- NA
        }  #if(NumberTracers>0)
        
        # In this section, the non-tracer element parameters are extracted in 
        # the same way as for the tracer element.
        
        ElementsNonTracerList <- list()
        
        if (NumberElementsNonTracer > 0) {
          
          NonTracer <- 1
          
          tmpElementNonTracer.vec <- vector()
          tmpElementNonTracerCount.vec <- vector()
          tmpElementZeroTracer.vec <- vector()
          tmpElementZeroTracerCount.vec <- vector()
          
          for (Element in seq_len(NumberElements)) {
            # tracers exist
            if (NumberTracers > 0) {
              # current element is no tracer
              if (names(ElementCount.vec)[Element] != IDTracer) {
                tmpElementNonTracer.vec <- c(tmpElementNonTracer.vec, names(ElementCount.vec)[Element])
                tmpElementNonTracerCount.vec <- c(tmpElementNonTracerCount.vec, as.numeric(ElementCount.vec[Element]))
              }
            } else {
              tmpElementZeroTracer.vec <- c(tmpElementZeroTracer.vec, names(ElementCount.vec)[Element])
              tmpElementZeroTracerCount.vec <- c(tmpElementZeroTracerCount.vec, as.numeric(ElementCount.vec[[Element]]))
              
            }  #NumberTracers
          }  #Element
          
          names(tmpElementNonTracerCount.vec) <- tmpElementNonTracer.vec
          names(tmpElementZeroTracerCount.vec) <- tmpElementZeroTracer.vec
          
          tmpFragmentList[["NonTracer"]] <- tmpElementNonTracerCount.vec
          tmpFragmentList[["ZeroTracer"]] <- tmpElementZeroTracerCount.vec
          
        } else {
          # NumberElementsNonTracer>0
          
          tmpFragmentList[["NonTracer"]] <- c()
          tmpFragmentList[["ZeroTracer"]] <- c()
          
        }
        
      } else if (UltraHighRes) 
      {
        
        if (NumberTracers > 0) {
          
          # For multiple tracer correction, the gathering of tracer parameters 
          # has to loop through the number of tracers present in a 
          # given molecule (-fragment).
          
          MaxLabel.vec <- vector()
          IDTracer.vec <- vector()
          nTracerMax.vec <- vector()
          
          for (TracerNo in seq_len(NumberTracers)) {
            MaxLabel <- as.numeric(TracerCount.vec[TracerNo])
            IDTracer <- names(TracerCount.vec)[TracerNo]
            
            MaxLabel.vec <- c(MaxLabel.vec, MaxLabel)
            IDTracer.vec <- c(IDTracer.vec, IDTracer)
            
            for (Element in seq_len(NumberElements)) {
              # if current element is a tracer
              if (names(tmpFragmentList[["Element"]][Element]) == IDTracer) 
              {
                nTracerMax <- tmpFragmentList[["Element"]][Element] %>% as.numeric
                nTracerMax.vec <- c(nTracerMax.vec, nTracerMax)
              }  #if
            }  #Element
          }  #TracerNo
          names(MaxLabel.vec) <- IDTracer.vec
          
          if (length(nTracerMax.vec) == NumberTracers) {
            
            names(nTracerMax.vec) <- IDTracer.vec
            
          }
          
          tmpFragmentList[["MaxLabel"]] <- MaxLabel.vec
          tmpFragmentList[["IDTracer"]] <- IDTracer.vec
          tmpFragmentList[["nTracerMax"]] <- nTracerMax.vec
          
          # Gaining natural abundance information associated with the 
          # tracer isotopes for probability calculations
          
          # Number of isotopes per Element
          NumberIso <- unlist(lapply(ElementArray, function(x) nrow(data.frame(x[["Isotopes"]]))))
          
          ### store information for each individual Tracer Element
          NatAbuTracerList <- list()
          NatAbuBaseList <- list()
          MassShiftTracerList <- list()
          
          MassShiftTracer.vec <- vector()
          
          for (TracerNo in names(TracerCount.vec)) {
            MassShiftTracer <- ElementArray[[TracerNo]][[2]]
            
            # for each Tracer isotope
            NatAbuTracer.vec <- vector()
            NatAbuBase.vec <- vector()
            
            for (IsotopeNo in seq_len(NumberIso[TracerNo])) {
              
              if (data.frame(ElementArray[[TracerNo]][[1]])[IsotopeNo, 2] == MassShiftTracer) {
                
                NatAbuTracer.vec <- c(NatAbuTracer.vec, data.frame(ElementArray[[TracerNo]][[1]])[IsotopeNo, 1])  # 1=>IsotopeAbundance
                
              } else if (data.frame(ElementArray[[TracerNo]][[1]])[IsotopeNo, 2] == 0) {
                
                NatAbuBase.vec <- c(NatAbuBase.vec, data.frame(ElementArray[[TracerNo]][[1]])[IsotopeNo, 1])
                
              } else {
                
                message(date(), " :: skipping IsotopeNo ", IsotopeNo, " for Tracer ", TracerNo)
                
              }
            }  #IsotopeNo
            
            NatAbuTracerList[[TracerNo]] <- NatAbuTracer.vec
            NatAbuBaseList[[TracerNo]] <- NatAbuBase.vec
            MassShiftTracerList[[TracerNo]] <- MassShiftTracer
            
          }  #TracerNo
          
          tmpFragmentList[["NatAbuTracer"]] <- unlist(NatAbuTracerList)
          tmpFragmentList[["NatAbuBase"]] <- unlist(NatAbuBaseList)
          tmpFragmentList[["MassShiftTracer"]] <- unlist(MassShiftTracerList)
        } else {
          message(date(), " :: NO TRACER FOR MOLECULE #", MoleculeNo, " AND FRAGMENT #", Fragment)
        }  #NumberTracers>0
      }  #if(UltraHighRes==1)
      
      FragmentList[[Fragment]] <- tmpFragmentList
    }  # Fragment
    
    names(FragmentList) <- stringr::str_c("Fragment_", rep(seq_len(NumberFragments)))
    
    MoleculeList[[MoleculeNo]] <- FragmentList
    
  }  # MoleculeNo
  
  names(MoleculeList) <- names(MoleculeLocation.vec)
  
  checkMoleculeDataLogic(MoleculeList, UltraHighRes, CorrectTracerImpurity, ElementArray, logEnvironment)
  
  logEnvironment$param$molecules <- names(MoleculeList)
  
  tracerLogInfo <- list()
  tracerNames <- character()
  
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
      
      tracerLogInfo[[tracer]] <- ElementArray[[tracer]]$Tracer.purity
      
    } else {
      
      tracerLogInfo[[tracer]] <- NA
      
    }
    
  }
  
  logEnvironment$param$tracers <- tracerLogInfo
  
  message(date(), " [MoleculeInfoExtraction()] [OK]")
  
  return(MoleculeInfo = MoleculeList)
}  # function()
