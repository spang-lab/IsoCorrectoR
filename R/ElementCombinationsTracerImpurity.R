# COMBINATION WITH TRACER IMPURITY STATES

# If correction for tracer impurity is switched on, the arrays TracerElemProbList and TracerElemMassShiftList are additonally combined with the
# probabilites and mass shifts associated with different numbers of 'impure' tracer atoms in the molecule(-fragment). This yields the arrays
# TracerImpurityCombProbList and TracerImpurityCombMassShiftList. They give the probability and mass shift associated with defined elemental isotope
# combinations and a defined number of 'impure' tracer atoms.

ElementCombinationsTracerImpurity <- function(ProbTracerImpurityList,MassShiftTracerImpurityList,NumberImpurityCombEff,
                                              TracerElemProbList,TracerElemMassShiftList,NumberTracerElemComb,
                                              Transitions,NumberTransitions,CalculationThreshold) {
  
  TracerImpurityCombProbList <- list()
  TracerImpurityCombMassShiftList <- list()
  NumberTracerImpurityComb <- vector()
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    
    tmpTracerImpurityCombProb_Array <- vector()
    tmpTracerImpurityCombMassShift_Array <- vector()
    ImpurityComb <- 1
    for (ImpureTracer in 0:NumberImpurityCombEff[[TransitionNo]]) {
      
      for (EleComb in seq_len(NumberTracerElemComb[TransitionNo])) {
        
        TracerImpurityCombProb <- ProbTracerImpurityList[[TransitionNo]][ImpureTracer + 1] * TracerElemProbList[[TransitionNo]][EleComb]
        
        if (TracerImpurityCombProb > CalculationThreshold) {
          
          tmpTracerImpurityCombProb_Array[ImpurityComb] <- TracerImpurityCombProb
          tmpTracerImpurityCombMassShift_Array[ImpurityComb] <- MassShiftTracerImpurityList[[TransitionNo]][ImpureTracer + 1] + TracerElemMassShiftList[[TransitionNo]][EleComb]
          
          NumberTracerImpurityComb[TransitionNo] <- ImpurityComb
          ImpurityComb <- ImpurityComb + 1
        } else {
          EleComb <- NumberTracerElemComb[TransitionNo] + 1
        }  #TracerImpurityCombProb>CalculationThreshold
      }  #EleComb
    }  #ImpureTracer
    TracerImpurityCombProbList[[TransitionNo]] <- tmpTracerImpurityCombProb_Array
    TracerImpurityCombMassShiftList[[TransitionNo]] <- tmpTracerImpurityCombMassShift_Array
  }  #TransitionNo
  
  names(TracerImpurityCombProbList) <- rownames(Transitions)
  names(TracerImpurityCombMassShiftList) <- rownames(Transitions)
  names(NumberTracerImpurityComb) <- rownames(Transitions)
  
  return(list("TracerImpurityCombProbList"=TracerImpurityCombProbList,
              "TracerImpurityCombMassShiftList"=TracerImpurityCombMassShiftList,
              "NumberTracerImpurityComb"=NumberTracerImpurityComb))
  
}