# CALCULATION OF PROBABILITY MATRIX
#
# In this section, the probability matrix used for data correction is produced.
# To this end, for each molecule to be corrected, molecule/fragment
# information from the molecule information input file is used
# to calculate the probability that a certain transition or full MS ion x
# is derived from another transition/full MS ion y due to
# natural stable isotope abundance/tracer purity.

ProbNormalRes <- function(MoleculeInfo,MoleculesTotal,ElementInfo,CorrectTracerElementCore,
                          CorrectTracerImpurity,CalculationThreshold, verbose) {

  if(verbose){message(date(), " :: calculating probability matrix with parameters: ")
  message(date(), " :: :: CorrectTracerElementCore: ", CorrectTracerElementCore)
  message(date(), " :: :: CorrectTracerImpurity: ", CorrectTracerImpurity)
  message(date(), " :: :: CalculationThreshold: ", CalculationThreshold)}
  
  MaxMassShiftList <- list()
  icmResultList <- list()
  ecResultList <- list()
  msprobResultList <- list()
  mstResultList <- list()
  
  CombinedProbList <- list()

  for (MoleculeNo in seq_len(MoleculesTotal)) {
    
    MoleculeData <- MoleculeInfo[[MoleculeNo]]
    MoleculeName <- names(MoleculeInfo[MoleculeNo])
    
    NumberFragments <-
      sum(stringr::str_detect(names(MoleculeData), "Fragment"))
    Transitions <- MoleculeData[["Transitions"]]
    
    tmpIcmResultList <- list()
    tmpEcResultList <- list()
    tmpMsprobResultList <- list()
    tmpMstResultList <- list()
    tmpMaxMassShiftList <- list()
    
    # The following functions are run for both fragments of a molecule
    # in the case of MS/MS measurements or just for one molecule(-fragment) for
    # full-MS measurements.
    
    for (Fragment in seq_len(NumberFragments)) {
      # For each molecule, the expected transitions are checked and MaxMassShift
      # gives the mass shifts of the labelled molecule(-fragment(s)) with the highest mass.
      
      MaxMassShift <- apply(Transitions, 2, max)[Fragment]
      tmpMaxMassShiftList[[Fragment]] <- MaxMassShift
      
      # CALCULATE ISOTOPE COMBINATION PROBABILITIES AND MASS SHIFTS
      
      IsoCombinationsResult <-
        IsoCombinationsMaster(
          MoleculeData = MoleculeData,
          MoleculeName = MoleculeName,
          MoleculeNo = MoleculeNo,
          Fragment = Fragment,
          ElementInfo = ElementInfo,
          CorrectTracerElementCore = CorrectTracerElementCore,
          CorrectTracerImpurity = CorrectTracerImpurity,
          CalculationThreshold = CalculationThreshold
        )
      
      tmpIcmResultList[[Fragment]] <- IsoCombinationsResult
      
      # CALCULATE ELEMENT COMBINATION PROBABILITIES AND MASS SHIFTS
      
      ElementCombinationsResult <-
        ElementCombinations(
          MoleculeInfo = MoleculeInfo,
          MoleculeNo = MoleculeNo,
          Fragment = Fragment,
          IsoCombinationsResult = IsoCombinationsResult,
          CalculationThreshold = CalculationThreshold,
          CorrectTracerImpurity = CorrectTracerImpurity
        )
      
      tmpEcResultList[[Fragment]] <- ElementCombinationsResult
      
      #CALCULATE CUMULATIVE PROBABILITIES OF MASS SHIFTS
      
      CumProbList <-
        MassShiftProbabilities(
          MoleculeInfo = MoleculeInfo,
          MoleculeNo = MoleculeNo,
          Fragment = Fragment,
          ElementCombinations = ElementCombinationsResult,
          MaxMassShift = MaxMassShift,
          CorrectTracerImpurity = CorrectTracerImpurity
        )
      
      tmpMsprobResultList[[Fragment]] <- CumProbList
      
      #GET MASS SHIFTS OF MEASUREMENTS/TRANSITIONS
      
      mstResult <-
        MassShiftTransitions(
          MoleculeInfo = MoleculeInfo,
          MoleculeNo = MoleculeNo,
          Fragment = Fragment
        )
      
      tmpMstResultList[[Fragment]] <- mstResult
      
    }#Fragment
    names(tmpIcmResultList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    names(tmpEcResultList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    names(tmpMsprobResultList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    names(tmpMstResultList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    names(tmpMaxMassShiftList) <-
      paste0("Fragment_", seq(seq_len(NumberFragments)))
    
    # CALCULATE FINAL PROBABILTIY MATRIX CombFragmentProb
    
    CombFragmentProb <-
      FragmentsCombinedProbability(
        MoleculeInfo = MoleculeInfo,
        MoleculeNo = MoleculeNo,
        CumProbList = tmpMsprobResultList,
        FragmentMassShift =
          tmpMstResultList
      )
    
    CombinedProbList[[MoleculeNo]] <-
      CombFragmentProb
    
    MaxMassShiftList[[MoleculeNo]] <- tmpMaxMassShiftList
    icmResultList[[MoleculeNo]] <- tmpIcmResultList
    ecResultList[[MoleculeNo]] <- tmpEcResultList
    msprobResultList[[MoleculeNo]] <- tmpMsprobResultList
    mstResultList[[MoleculeNo]] <- tmpMstResultList
    
  }#MoleculeNo  
  
  names(MaxMassShiftList) <- names(MoleculeInfo)
  names(icmResultList) <- names(MoleculeInfo)
  names(ecResultList) <- names(MoleculeInfo)
  names(msprobResultList) <- names(MoleculeInfo) # CumProbList
  names(mstResultList) <- names(MoleculeInfo) # FragmentMassShift
  names(CombinedProbList) <-
    names(MoleculeInfo) # CombFragmentProb
  
  if(verbose){message(date(), " :: calculating probability matrix [OK]")}
  
  return(CombinedProbList)
  
}