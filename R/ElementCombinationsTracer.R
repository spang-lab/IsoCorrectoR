# CALCULATION OF THE TRACER ELEMENT - NON-TRACER ELEMENT COMBINATIONS

# Given that there are more elements than just the tracer element, the element combinations from CombElemProbList are now multiplied with the isotope
# combinations of the tracer element, for all labelling states. This is analogous to the previous algorithms from this function.  However, in the
# computation of the mass shift array, the intrinsic mass shift (due to tracer incorporation) associated with the respective labelling state is added in
# addtion.  In the end, the probability array TracerElemProbList and the mass shift array TracerElemMassShiftList are derived. They contain the
# probabilities and mass shifts of all relevant elemental combinations of isotope combinations for all labelling states of a given molecule(-fragment).

ElementCombinationsTracer <- function(ProbTracerList, MassShiftTracerList, NumberIsoCombTracerEff, 
                                         CombElemProbList, CombElemMassShiftList, NumberEleComb,
                                         Fragment, NumberElementsNonTracer, NumberTracers, Transitions,
                                         NumberTransitions, CalculationThreshold) {

  TracerElemProbList <- list()
  TracerElemMassShiftList <- list()
  NumberTracerElemComb <- vector()
  
  if (NumberTracers > 0) {
    
    if (NumberElementsNonTracer > 0) {
      
      for (TransitionNo in seq_len(NumberTransitions)) {
        
        tmpTracerElemProb_Array <- vector()
        tmpTracerElemMassShift_Array <- vector()
        EleComb <- 1
        
        for (ElementCombination in seq_len(NumberEleComb[[NumberElementsNonTracer]])) {
          
          for (IsoCombination in seq_len(NumberIsoCombTracerEff[[TransitionNo]])) {
            
            TracerElemProb <- CombElemProbList[[NumberElementsNonTracer]][[ElementCombination]] * ProbTracerList[[TransitionNo]][[IsoCombination]]
            if (TracerElemProb > CalculationThreshold) {
              tmpTracerElemProb_Array[EleComb] <- TracerElemProb
              
              tmpTracerElemMassShift_Array[EleComb] <- CombElemMassShiftList[[NumberElementsNonTracer]][ElementCombination] + MassShiftTracerList[[TransitionNo]][IsoCombination] + 
                Transitions[TransitionNo, Fragment]
              
              NumberTracerElemComb[TransitionNo] <- EleComb
              
              EleComb <- EleComb + 1
              
            } else {
              IsoCombination <- NumberIsoCombTracerEff[[TransitionNo]] + 1
            }  #TracerElemProb > CalculationThreshold
          }  #IsoCombination
        }  #ElementCombination
        
        TracerElemProbList[[TransitionNo]] <- tmpTracerElemProb_Array
        TracerElemMassShiftList[[TransitionNo]] <- tmpTracerElemMassShift_Array
        
      }  #TransitionNo
    } else {
      
      # If the only element considered is the tracer element, the probability- and mass shift arrays from 'IsoCombinationsResult' are directly fed into the
      # probability array TracerElemProbList and the mass shift array TracerElemMassShiftList.  However, for each transition its intrinsic mass shift is
      # added to the mass shift array in addition.
      
      for (i in seq_len(NumberTransitions)) {
        NumberTracerElemComb[i] <- NumberIsoCombTracerEff[[i]]
      }  #i
      
      for (TransitionNo in seq_len(NumberTransitions)) {
        tmpTracerElemProb_Array <- vector()
        tmpTracerElemMassShift_Array <- vector()
        
        for (IsoCombination in seq_len(NumberIsoCombTracerEff[[TransitionNo]])) {
          
          tmpTracerElemProb_Array[IsoCombination] <- ProbTracerList[[TransitionNo]][IsoCombination]
          
          tmpTracerElemMassShift_Array[IsoCombination] <- MassShiftTracerList[[TransitionNo]][IsoCombination] + Transitions[TransitionNo, Fragment]
          
        }  #IsoCombination
        TracerElemProbList[[TransitionNo]] <- tmpTracerElemProb_Array
        TracerElemMassShiftList[[TransitionNo]] <- tmpTracerElemMassShift_Array
      }  #TransitionNo
    }  #NumberElementsNonTracer>0
  } else {
    
    # If the molecule(-fragment) in question contains no tracer element, the probability array TracerElemProbList and the mass shift array
    # TracerElemMassShiftList are equal to CombElemProbList and CombElemMassShiftList for all transitions of this molecule(-fragment).
    
    for (TransitionNo in seq_len(NumberTransitions)) {
      
      tmpTracerElemProb_Array <- vector()
      tmpTracerElemMassShift_Array <- vector()
      
      for (i in seq_len(NumberEleComb[[NumberElementsNonTracer]])) {
        
        tmpTracerElemProb_Array[i] <- CombElemProbList[[NumberElementsNonTracer]][i]
        tmpTracerElemMassShift_Array[i] <- CombElemMassShiftList[[NumberElementsNonTracer]][i]
      }  #i
      
      TracerElemProbList[[TransitionNo]] <- tmpTracerElemProb_Array
      TracerElemMassShiftList[[TransitionNo]] <- tmpTracerElemMassShift_Array
      NumberTracerElemComb[TransitionNo] <- NumberEleComb[[NumberElementsNonTracer]]
    }  #TransitionNo
    
  }  #NumberTracers>0
  
  names(TracerElemProbList) <- rownames(Transitions)
  names(TracerElemMassShiftList) <- rownames(Transitions)
  names(NumberTracerElemComb) <- rownames(Transitions)
  
  return(list("TracerElemProbList"=TracerElemProbList,
              "TracerElemMassShiftList"=TracerElemMassShiftList,
              "NumberTracerElemComb"=NumberTracerElemComb))
  
}