# Calculation of the probability matrix for Ultra High Resolution (UHR) 
# data correction

ProbUltraHighRes <- function(MoleculeInfo,MoleculesTotal,ElementInfo,CorrectTracerElementCore,CorrectTracerImpurity,
                             CalculationThreshold_UHR, verbose) {
  
  if(verbose){message(date(), " :: ProbUltraHighRes")
  message(date(), " :: CorrectTracerElementCore: ", CorrectTracerElementCore)
  message(date(), " :: CorrectTracerImpurity: ", CorrectTracerImpurity)
  message(date(), " :: CalculationThreshold_UHR: ", CalculationThreshold_UHR)}
  
  CombinedProbList <- list()
  
  for (MoleculeNo in seq_len(MoleculesTotal)) {
    
    MoleculeData <- MoleculeInfo[[MoleculeNo]]
    
    NatAbuImpurityList <- IsoCombinationsUHR(MoleculeData=MoleculeData,
                                                 ElementInfo=ElementInfo,
                                                 CorrectTracerElementCore=CorrectTracerElementCore,
                                                 CorrectTracerImpurity=CorrectTracerImpurity,
                                                 verbose=verbose)
    
    CumProbList <- TracerCountProbabilitiesUHR(MoleculeData=MoleculeData,
                                               NatAbuImpurityList=NatAbuImpurityList,
                                               verbose=verbose)
    
    CombinedProb <- TracersCombinedProbabilityUHR(MoleculeData=MoleculeData,
                                                  CumProbList=CumProbList,
                                                  CalculationThreshold_UHR=CalculationThreshold_UHR)
    
    CombinedProbList[[MoleculeNo]] <- CombinedProb
    
  }

  names(CombinedProbList) <- names(MoleculeInfo)
  
  if(verbose){message(date(), " :: leaving ProbUltraHighRes()")}
  
  return(CombinedProbList)
  
}
