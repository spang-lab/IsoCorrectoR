# Computation of isotope combinations, associated probabilities and mass shifts 
# for an element

#' @importFrom stringr str_c

IsoCombinationsMaster <- function(MoleculeData, MoleculeName, MoleculeNo, Fragment, ElementInfo, CalculationThreshold, CorrectTracerImpurity, CorrectTracerElementCore) {
  
  MoleculeFragmentData <- MoleculeData[[stringr::str_c("Fragment_", Fragment)]]  # get individual fragment data
  
  Transitions <- MoleculeData[["Transitions"]]
  
  NumberTransitions <- nrow(Transitions)
  
  NumberTracers <- length(MoleculeFragmentData[["Tracer"]])
  
  # Isotope combinations for tracer element
  
  IsoCombinationsTracerResult <- IsoCombinationsTracer(MoleculeFragmentData=MoleculeFragmentData,
                                                       Fragment=Fragment,
                                                       ElementInfo = ElementInfo,
                                                       Transitions=Transitions,
                                                       NumberTracers=NumberTracers,
                                                       CorrectTracerImpurity=CorrectTracerImpurity, 
                                                       CorrectTracerElementCore=CorrectTracerElementCore,
                                                       CalculationThreshold=CalculationThreshold)
  
  ProbTracerList <- IsoCombinationsTracerResult$ProbTracerList
  MassShiftTracerList <- IsoCombinationsTracerResult$MassShiftTracerList
  NumberIsoCombTracerEff <- IsoCombinationsTracerResult$NumberIsoCombTracerEff
  IsotopesTracer <- IsoCombinationsTracerResult$IsotopesTracer
  LabelPresent <- IsoCombinationsTracerResult$LabelPresent
  
  # Isotope combinations for non-tracer elements
  
  if (NumberTracers > 0) {
    
    IDTracer <- MoleculeFragmentData[["IDTracer"]] 
    
  } else {
    
    IDTracer <- NA
    
  }
  
  IsoCombinationsNonTracerResult <- IsoCombinationsNonTracer(MoleculeFragmentData=MoleculeFragmentData, 
                                                            Fragment=Fragment,
                                                            ElementInfo=ElementInfo,
                                                            NumberTracers=NumberTracers,
                                                            IDTracer=IDTracer,
                                                            CalculationThreshold=CalculationThreshold)
  
  ProbElemList <- IsoCombinationsNonTracerResult$ProbElemList
  MassShiftElemList <- IsoCombinationsNonTracerResult$MassShiftElemList
  NumberIsoCombEff <- IsoCombinationsNonTracerResult$NumberIsoCombEff
  IsotopesElem <- IsoCombinationsNonTracerResult$IsotopesElem

  # Compute tracer impurity combinations
  
  ProbTracerImpurityList <- list()
  MassShiftTracerImpurityList <- list()
  NumberImpurityCombEff <- list()
  
  if (CorrectTracerImpurity) {
    
    if (NumberTracers > 0) {
      
      IsoCombinationsImpurityResult <- IsoCombinationsTracerImpurity(MoleculeFragmentData=MoleculeFragmentData,
                                                                     ElementInfo=ElementInfo,
                                                                     Transitions=Transitions,
                                                                     NumberTransitions=NumberTransitions,
                                                                     IDTracer=IDTracer,
                                                                     LabelPresent=LabelPresent,
                                                                     CalculationThreshold=CalculationThreshold)
      
      ProbTracerImpurityList <- IsoCombinationsImpurityResult$ProbTracerImpurityList
      MassShiftTracerImpurityList <- IsoCombinationsImpurityResult$MassShiftTracerImpurityList
      NumberImpurityCombEff <- IsoCombinationsImpurityResult$NumberImpurityCombEff
      
    } 
  } else {
    ProbTracerImpurityList <- "No Correction for Tracer Impurity performed"
    MassShiftTracerImpurityList <- "No Correction for Tracer Impurity performed"
    NumberImpurityCombEff <- "No Correction for Tracer Impurity performed"
  } #CorrectTracerImpurity=TRUE
  
  IsoCombinationsResult <- list(ProbTracerList = ProbTracerList, ProbElemList = ProbElemList, ProbTracerImpurityList = ProbTracerImpurityList, MassShiftTracerList = MassShiftTracerList, 
                     MassShiftElemList = MassShiftElemList, MassShiftTracerImpurityList = MassShiftTracerImpurityList, NumberIsoCombEff = unlist(NumberIsoCombEff), 
                     NumberImpurityCombEff = NumberImpurityCombEff, NumberIsoCombTracerEff = unlist(NumberIsoCombTracerEff), IsotopesTracer = IsotopesTracer, IsotopesElem = IsotopesElem, 
                     LabelPresent = LabelPresent)
  
  return(IsoCombinationsResult)
  
}  #IsoCombinationsMaster()
