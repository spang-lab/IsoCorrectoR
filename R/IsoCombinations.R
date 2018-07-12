#' @importFrom dplyr last
#' @importFrom dplyr arrange_
#' @importFrom stringr str_c
#' @importFrom tibble as_tibble
#' @importFrom magrittr '%>%'

IsoCombinations <- function(MoleculeNo, Fragment, x, Element, ElementArray, AvailablePlacesTotal, IsoCluster, Prob_threshold_isotopes) {
  
  # NumberIso is the number of isotopes per element. 
  # The number of nested loops used depends on the number of 
  # isotopes of the current element.
  
  NumberIso <- nrow(ElementArray[[Element]][["Isotopes"]][[1]])
  
  Combinations <- vector()
  Combinations[1] <- 1
  
  # Calculation of Combinations[IsotopeNo+1], the number of isotope combinations
  # when considering IsotopeNo different isotopes (also including combinations
  # where the sum of isotopes exceeds the number of available places IsoCluster).
  
  for (IsotopeNo in seq_len(NumberIso)) {
    Combinations[IsotopeNo + 1] <- Combinations[IsotopeNo] * ((IsoCluster) + 1)
  }  #IsotopeNo
  
  TotalCombinations <- dplyr::last(Combinations)  # get last element
  
  # Computation of the isotope combinations list `CombinationsArrayPreList`.  
  # At this point, combinations where the sum of the different isotopes exceeds
  # IsoCluster are still considered, leading to invalid combinations that are 
  # later removed. In its columns, CombinationsArrayPreList contains the number
  # of occurences of a given isotope in a given combination. 
  # Each row corresponds to a combination.
  
  CombinationsArrayPreList <- list()
  
  for (IsotopeNo in seq_len(NumberIso)) {
    
    p <- 1
    tmpCombinationsArrayPre <- vector()
    
    # first, for each column (isotope) the repetitive column element is computed
    
    for (IsoCount in 0:IsoCluster) {
      tmpCombinationsArrayPre[p:(p + Combinations[IsotopeNo] - 1)] <- IsoCount
      p <- p + Combinations[IsotopeNo]
    }  #IsoCount
    
    # Then, the repetitive column element is repeatedly inserted until the 
    # number of rows matches `TotalCombinations`
    
    for (Repeat in seq_len((TotalCombinations/Combinations[IsotopeNo + 1])) - 1) {
      
      if (Repeat > 0) {
        tmpCombinationsArrayPre[p:(p + Combinations[IsotopeNo + 1] - 1)] <- tmpCombinationsArrayPre[seq_len(Combinations[IsotopeNo + 1])]
        p <- p + Combinations[IsotopeNo + 1]
      }
      
    }  #Repeat
    
    CombinationsArrayPreList[[IsotopeNo]] <- tmpCombinationsArrayPre
  }  #IsotopeNo
  names(CombinationsArrayPreList) <- stringr::str_c("Isotope_", seq_len(NumberIso))
  
  #CombinationsArrayPreList is sorted to get the most probable combinations on top
  
  CombinationsArrayPreList2 <- CombinationsArrayPreList %>% as_tibble %>% dplyr::arrange_(.dots = names(CombinationsArrayPreList)) %>% as.list
  
  CombinationsArrayPre.df <- as.data.frame(lapply(CombinationsArrayPreList2, function(x) rev(x)))
  
  RowSum <- apply(CombinationsArrayPre.df, 1, sum)
  j <- 1
  
  # CombinationsArray is derived from CombinationsArrayPreList by removing 
  # all combinations that do not have a row sum of IsoCluster
  
  CombinationsArray <- NULL
  
  for (i in seq_len(TotalCombinations)) {
    
    if (RowSum[i] == IsoCluster) 
    {
      
      tmpA <- vector()
      
      tmpA <- as.numeric(CombinationsArrayPre.df[i, seq_len(NumberIso)])
      tmpA[1] <- tmpA[1] + (AvailablePlacesTotal - IsoCluster)
      CombinationsArray <- rbind(CombinationsArray, tmpA)
      
      j <- j + 1
      
    }  #RowSum[i]==IsoCluster
  }  #i
  
  # The number of combinations in CombinationsArray is derived, 
  # CombinationsArray is stored in Isotopes
  
  NumberIsoCombMax <- nrow(CombinationsArray)
  
  Isotopes <- CombinationsArray
  
  # The probabilities and mass shifts of the isotope combinations in 
  # CombinationsArray are computed using IsoCombinationsProbMassShift()
  
  IsoComb <- 1
  ProbArray.vec <- vector()
  MassShiftArray.vec <- vector()
  NumberIsoComb.vec <- vector()
  
  for (IsoCombPre in seq_len(NumberIsoCombMax)) {
    
    
    icpms <- IsoCombinationsProbMassShift(Element = Element, NumberIso = NumberIso, ElementArray = ElementArray, CombinationsArray = CombinationsArray, 
                                          IsoCombPre = IsoCombPre, AvailablePlacesTotal = AvailablePlacesTotal)
    
    # Due to Prob_threshold_isotopes which limits calculations to probability 
    # values large enough to be relevant, probability and mass shift are 
    # usually not computed for all isotope combinations.
    
    if (icpms[["ProbIsoComb"]] > Prob_threshold_isotopes) {
      
      ProbArray.vec[IsoComb] <- icpms[["ProbIsoComb"]]
      MassShiftArray.vec[IsoComb] <- icpms[["MassShiftIsoComb"]]
      NumberIsoComb <- IsoComb
      IsoComb <- IsoComb + 1
    } else {
      if (IsoComb == 1) {
        NumberIsoComb <- 0
      }
      
      IsoCombPre <- NumberIsoCombMax + 1
    }
    
  }  #IsoCombPre
  
  return(list(Isotopes = Isotopes, ProbArray = ProbArray.vec, MassShiftArray = MassShiftArray.vec, NumberIsoComb = NumberIsoComb))
  
  
}  #IsoCombinations
