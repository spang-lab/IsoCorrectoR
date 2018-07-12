# Calculation of probability and mass shift associated with an isotope combination

IsoCombinationsProbMassShift <- function(Element, NumberIso, ElementArray, CombinationsArray, IsoCombPre, AvailablePlacesTotal) {
  
  VariableTerm <- 1
  MassShiftIsoComb <- 0
  
  for (IsotopeNo in seq_len(NumberIso)) {
    
    # abundance
    VariableTerm <- VariableTerm * (data.frame(ElementArray[[Element]][[1]])[IsotopeNo, 1])^as.numeric(CombinationsArray[IsoCombPre, IsotopeNo])/(factorial(as.numeric(CombinationsArray[IsoCombPre, IsotopeNo])))
    
    # mass shift
    MassShiftIsoComb <- MassShiftIsoComb + data.frame(ElementArray[[Element]][[1]])[IsotopeNo, 2] * as.numeric(CombinationsArray[IsoCombPre, IsotopeNo])
  }  #IsotopeNo
  
  ProbIsoComb <- VariableTerm * (factorial(AvailablePlacesTotal))
  
  return(list(ProbIsoComb = ProbIsoComb, MassShiftIsoComb = MassShiftIsoComb))
  
}
