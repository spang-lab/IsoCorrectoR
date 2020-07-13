# Calculation of probability and mass shift associated with an isotope combination

IsoCombinationsProbMassShift <- function(Element, NumberIso, ElementInfo, CombinationsArray, IsoCombPre, AvailablePlacesTotal) {
  
  ProbFactorialTermLog <- lfactorial(AvailablePlacesTotal)
  ProbExpTerm <- 1
  MassShiftIsoComb <- 0
  ElementTable <- data.frame(ElementInfo[[Element]][[1]])
  
  for (IsotopeNo in seq_len(NumberIso)) {
    
    IsotopeCount <- as.numeric(CombinationsArray[IsoCombPre, IsotopeNo])
    IsotopeAbundance <- ElementTable[IsotopeNo, 1]
    IsotopeMassShift <- ElementTable[IsotopeNo, 2]
    
    # Probability factorial term (log factorial as factorial can lead to numbers that are too big)
    ProbFactorialTermLog <- ProbFactorialTermLog - lfactorial(IsotopeCount)
    
    # Probability exponential term
    ProbExpTerm <- ProbExpTerm * (IsotopeAbundance^IsotopeCount)
    
    # mass shift
    MassShiftIsoComb <- MassShiftIsoComb + IsotopeMassShift * IsotopeCount
  }  #IsotopeNo
  
  ProbIsoComb <- exp(ProbFactorialTermLog) * ProbExpTerm
  
  return(list(ProbIsoComb = ProbIsoComb, MassShiftIsoComb = MassShiftIsoComb))
  
}
