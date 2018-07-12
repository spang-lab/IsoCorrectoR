
ElementCombinations <- function(MoleculeArray, MoleculeNo, Fragment, IsoCombinationsMaster, Prob_threshold_elements, Prob_threshold_tracer, Prob_threshold_impurity, 
    CorrectTracerImpurity) {
    
    MoleculeData <- MoleculeArray[[MoleculeNo]][[str_c("Fragment_", Fragment)]]
    ProbTracer_Array <- IsoCombinationsMaster[["ProbTracer_Array"]]
    MassShiftTracer_Array <- IsoCombinationsMaster[["MassShiftTracer_Array"]]
    NumberIsoCombTracerEff <- IsoCombinationsMaster[["NumberIsoCombTracerEff"]]
    ProbElem_Array <- IsoCombinationsMaster[["ProbElem_Array"]]
    MassShiftElem_Array <- IsoCombinationsMaster[["MassShiftElem_Array"]]
    NumberIsoCombEff <- IsoCombinationsMaster[["NumberIsoCombEff"]]
    ProbTracerImpurity_Array <- IsoCombinationsMaster[["ProbTracerImpurity_Array"]]
    MassShiftTracerImpurity_Array <- IsoCombinationsMaster[["MassShiftTracerImpurity_Array"]]
    NumberImpurityCombEff <- IsoCombinationsMaster[["NumberImpurityCombEff"]]
    ElementsNonTracer <- MoleculeData[["NonTracer"]]
    
    if (length(ElementsNonTracer) == 0) {
        ElementsNonTracer <- MoleculeData[["ZeroTracer"]]
    }
    
    NumberElementsNonTracer <- length(ElementsNonTracer)
    NumberTracers <- length(MoleculeData[["Tracer"]])
    
    Transitions <- MoleculeArray[[MoleculeNo]][["Transitions"]]
    NumberTransitions <- nrow(Transitions)
    
    # a)CALCULATION OF THE NON-TRACER ELEMENT COMBINATIONS If the number of non-tracer elements is > 1, all isotope combination probabilities of non-tracer
    # elements 1 and 2 of a molecule(-fragment) are multiplied with each other and the corresponding mass shifts are added.  Like in 'IsoCombinations' there
    # is a probability threshold which resets to the outer loop if the probability calculated is below.
    
    CombElemProb_Array <- list()
    tmpCombElemProb_Array <- vector()
    
    CombElemMassShift_Array <- list()
    tmpCombElemMassShift_Array <- vector()
    
    NumberEleComb <- list()
    tmpNumberEleComb <- vector()
    
    if (NumberElementsNonTracer > 1) {
        
        EleComb <- 1
        
        for (IsoCombinationA in seq_len(as.numeric(NumberIsoCombEff[1]))) {
            
            for (IsoCombinationB in seq_len(as.numeric(NumberIsoCombEff[[2]]))) {
                
                tmpCombElemProb <- ProbElem_Array[[1]][IsoCombinationA] * ProbElem_Array[[2]][IsoCombinationB]
                
                if (tmpCombElemProb > Prob_threshold_elements) {
                  
                  tmpCombElemProb_Array[EleComb] <- tmpCombElemProb
                  tmpCombElemMassShift_Array[EleComb] <- MassShiftElem_Array[[1]][IsoCombinationA] + MassShiftElem_Array[[2]][IsoCombinationB]
                  tmpNumberEleComb <- EleComb
                  
                  EleComb <- EleComb + 1
                } else {
                  IsoCombinationB <- NumberIsoCombEff[[2]] + 1
                }  #tmpCombElemProb > Prob_threshold_elements
            }  #IsoCombinationB
        }  #IsoCombinationA
        
        CombElemProb_Array[[2]] <- tmpCombElemProb_Array
        CombElemMassShift_Array[[2]] <- tmpCombElemMassShift_Array
        NumberEleComb[[2]] <- tmpNumberEleComb
        
        rm(tmpCombElemProb_Array)
        rm(tmpCombElemMassShift_Array)
        rm(tmpNumberEleComb)
        
    } else if (NumberElementsNonTracer > 0) 
        {
            
            # If the number of non-tracer elements is 1, the probability and mass shift arrays are equal to those generated in 'IsoCombinationsMaster'.
            
            tmpCombElemProb_Array[seq_len(NumberIsoCombEff[[1]])] <- ProbElem_Array[[1]][seq_len(NumberIsoCombEff[[1]])]
            tmpCombElemMassShift_Array[seq_len(NumberIsoCombEff[[1]])] <- MassShiftElem_Array[[1]][seq_len(NumberIsoCombEff[[1]])]
            tmpNumberEleComb <- NumberIsoCombEff[[1]]
            
            CombElemProb_Array[[1]] <- tmpCombElemProb_Array
            CombElemMassShift_Array[[1]] <- tmpCombElemMassShift_Array
            NumberEleComb[[1]] <- tmpNumberEleComb
            
            rm(tmpCombElemProb_Array)
            rm(tmpCombElemMassShift_Array)
            rm(tmpNumberEleComb)
            
        }  #NumberElementsNonTracer>1
    
    # If there are more than 2 non-tracer elements, the element combinations are calculated in iterations of the following outer loop until all elements have
    # been covered. The algorithm uses the CombElemProb_Array covering the previously combined element probabilities and multiplies these with the isotope
    # combination probabilities of the next element (ProbElem_Array).  This is fed back into CombElemProb_Array at the current element index of the loop). The
    # same is done additively for the mass shift array CombElemMassShift_Array.
    
    if (NumberElementsNonTracer > 2) 
        {
            
            tmpCombElemProb_Array <- vector()
            tmpCombElemMassShift_Array <- vector()
            tmpNumberEleComb <- vector()
            
            for (NonTracer in 3:NumberElementsNonTracer) {
                EleComb <- 1
                for (ElementCombination in seq_len(NumberEleComb[[NonTracer - 1]])) {
                  for (IsoCombination in seq_len(NumberIsoCombEff[[NonTracer]])) {
                    
                    tmpCombElemProb <- CombElemProb_Array[[NonTracer - 1]][[ElementCombination]] * ProbElem_Array[[NonTracer]][[IsoCombination]]
                    if (tmpCombElemProb > Prob_threshold_elements) {
                      
                      tmpCombElemProb_Array[EleComb] <- tmpCombElemProb
                      tmpCombElemMassShift_Array[EleComb] <- CombElemMassShift_Array[[NonTracer - 1]][[ElementCombination]] + MassShiftElem_Array[[NonTracer]][[IsoCombination]]
                      tmpNumberEleComb <- EleComb
                      EleComb <- EleComb + 1
                      
                    } else {
                      
                      IsoCombination <- NumberIsoCombEff[[NonTracer]] + 1
                      
                    }  #tmpCombElemProb > Prob_threshold_elements
                  }  #IsoCombination
                  
                }  #ElementCombination
                CombElemProb_Array[[NonTracer]] <- tmpCombElemProb_Array
                CombElemMassShift_Array[[NonTracer]] <- tmpCombElemMassShift_Array
                NumberEleComb[[NonTracer]] <- tmpNumberEleComb
                
            }  #NonTracer
        }  #NumberElementsNonTracer>2
    
    # b) CALCULATION OF THE TRACER ELEMENT - NON-TRACER ELEMENT COMBINATIONS
    
    # Given that there are more elements than just the tracer element, the element combinations from CombElemProb_Array are now multiplied with the isotope
    # combinations of the tracer element, for all labelling states. This is analogous to the previous algorithms from this function.  However, in the
    # computation of the mass shift array, the intrinsic mass shift (due to tracer incorporation) associated with the respective labelling state is added in
    # addtion.  In the end, the probability array TracerElemProb_Array and the mass shift array TracerElemMassShift_Array are derived. They contain the
    # probabilities and mass shifts of all relevant elemental combinations of isotope combinations for all labelling states of a given molecule(-fragment).
    
    TracerElemProb_Array <- list()
    TracerElemMassShift_Array <- list()
    
    NumberTracerElemComb <- vector()
    if (NumberTracers > 0) {
        
        if (NumberElementsNonTracer > 0) {
            
            for (TransitionNo in seq_len(NumberTransitions)) {
                
                tmpTracerElemProb_Array <- vector()
                tmpTracerElemMassShift_Array <- vector()
                EleComb <- 1
                
                for (ElementCombination in seq_len(NumberEleComb[[NumberElementsNonTracer]])) {
                  
                  for (IsoCombination in seq_len(NumberIsoCombTracerEff[[TransitionNo]])) {
                    
                    TracerElemProb <- CombElemProb_Array[[NumberElementsNonTracer]][[ElementCombination]] * ProbTracer_Array[[TransitionNo]][[IsoCombination]]
                    if (TracerElemProb > Prob_threshold_tracer) {
                      tmpTracerElemProb_Array[EleComb] <- TracerElemProb
                      
                      tmpTracerElemMassShift_Array[EleComb] <- CombElemMassShift_Array[[NumberElementsNonTracer]][ElementCombination] + MassShiftTracer_Array[[TransitionNo]][IsoCombination] + 
                        Transitions[TransitionNo, Fragment]
                      
                      NumberTracerElemComb[TransitionNo] <- EleComb
                      
                      EleComb <- EleComb + 1
                      
                    } else {
                      IsoCombination <- NumberIsoCombTracerEff[[TransitionNo]] + 1
                    }  #TracerElemProb > Prob_threshold_tracer
                  }  #IsoCombination
                }  #ElementCombination
                
                TracerElemProb_Array[[TransitionNo]] <- tmpTracerElemProb_Array
                TracerElemMassShift_Array[[TransitionNo]] <- tmpTracerElemMassShift_Array
                
            }  #TransitionNo
        } else {
            
            # If the only element considered is the tracer element, the probability- and mass shift arrays from 'IsoCombinationsMaster' are directly fed into the
            # probability array TracerElemProb_Array and the mass shift array TracerElemMassShift_Array.  However, for each transition its intrinsic mass shift is
            # added to the mass shift array in addition.
            
            for (i in seq_len(NumberTransitions)) {
                NumberTracerElemComb[i] <- NumberIsoCombTracerEff[[i]]
            }  #i
            
            for (TransitionNo in seq_len(NumberTransitions)) {
                tmpTracerElemProb_Array <- vector()
                tmpTracerElemMassShift_Array <- vector()
                
                for (IsoCombination in seq_len(NumberIsoCombTracerEff[[TransitionNo]])) {
                  
                  tmpTracerElemProb_Array[IsoCombination] <- ProbTracer_Array[[TransitionNo]][IsoCombination]
                  
                  tmpTracerElemMassShift_Array[IsoCombination] <- MassShiftTracer_Array[[TransitionNo]][IsoCombination] + Transitions[TransitionNo, Fragment]
                  
                }  #IsoCombination
                TracerElemProb_Array[[TransitionNo]] <- tmpTracerElemProb_Array
                TracerElemMassShift_Array[[TransitionNo]] <- tmpTracerElemMassShift_Array
            }  #TransitionNo
        }  #NumberElementsNonTracer>0
    } else {
        
        # If the molecule(-fragment) in question contains no tracer element, the probability array TracerElemProb_Array and the mass shift array
        # TracerElemMassShift_Array are equal to CombElemProb_Array and CombElemMassShift_Array for all transitions of this molecule(-fragment).
        
        for (TransitionNo in seq_len(NumberTransitions)) {
            
            tmpTracerElemProb_Array <- vector()
            tmpTracerElemMassShift_Array <- vector()
            
            for (i in seq_len(NumberEleComb[[NumberElementsNonTracer]])) {
                
                tmpTracerElemProb_Array[i] <- CombElemProb_Array[[NumberElementsNonTracer]][i]
                tmpTracerElemMassShift_Array[i] <- CombElemMassShift_Array[[NumberElementsNonTracer]][i]
            }  #i
            
            TracerElemProb_Array[[TransitionNo]] <- tmpTracerElemProb_Array
            TracerElemMassShift_Array[[TransitionNo]] <- tmpTracerElemMassShift_Array
            NumberTracerElemComb[TransitionNo] <- NumberEleComb[[NumberElementsNonTracer]]
        }  #TransitionNo
        
    }  #NumberTracers>0
    
    # c)COMBINATION WITH TRACER IMPURITY STATES
    
    # If correction for tracer impurity is switched on, the arrays TracerElemProb_Array and TracerElemMassShift_Array are additonally combined with the
    # probabilites and mass shifts associated with different numbers of 'impure' tracer atoms in the molecule(-fragment). This yields the arrays
    # TracerImpurityCombProb_Array and TracerImpurityCombMassShift_Array. They give the probability and mass shift associated with defined elemental isotope
    # combinations and a defined number of 'impure' tracer atoms.
    
    TracerImpurityCombProb_Array <- list()
    TracerImpurityCombMassShift_Array <- list()
    
    NumberTracerImpurityComb <- vector()
    
    if (CorrectTracerImpurity) 
        {
            
            if (NumberTracers > 0) 
                {
                  
                  
                  for (TransitionNo in seq_len(NumberTransitions)) {
                    
                    tmpTracerImpurityCombProb_Array <- vector()
                    tmpTracerImpurityCombMassShift_Array <- vector()
                    ImpurityComb <- 1
                    for (ImpureTracer in 0:NumberImpurityCombEff[[TransitionNo]]) {
                      
                      for (EleComb in seq_len(NumberTracerElemComb[TransitionNo])) {
                        
                        TracerImpurityCombProb <- ProbTracerImpurity_Array[[TransitionNo]][ImpureTracer + 1] * TracerElemProb_Array[[TransitionNo]][EleComb]
                        
                        if (TracerImpurityCombProb > Prob_threshold_impurity) {
                          
                          tmpTracerImpurityCombProb_Array[ImpurityComb] <- TracerImpurityCombProb
                          tmpTracerImpurityCombMassShift_Array[ImpurityComb] <- MassShiftTracerImpurity_Array[[TransitionNo]][ImpureTracer + 1] + TracerElemMassShift_Array[[TransitionNo]][EleComb]
                          
                          NumberTracerImpurityComb[TransitionNo] <- ImpurityComb
                          ImpurityComb <- ImpurityComb + 1
                        } else {
                          EleComb <- NumberTracerElemComb[TransitionNo] + 1
                        }  #TracerImpurityCombProb>Prob_threshold_impurity
                      }  #EleComb
                    }  #ImpureTracer
                    TracerImpurityCombProb_Array[[TransitionNo]] <- tmpTracerImpurityCombProb_Array
                    TracerImpurityCombMassShift_Array[[TransitionNo]] <- tmpTracerImpurityCombMassShift_Array
                  }  #TransitionNo
                  names(TracerImpurityCombProb_Array) <- rownames(Transitions)
                  names(TracerImpurityCombMassShift_Array) <- rownames(Transitions)
                  names(NumberTracerImpurityComb) <- rownames(Transitions)
                  
                }  #NumberTracers>0
        }  #CorrectTracerImpurity==TRUE
    
    names(TracerElemProb_Array) <- rownames(Transitions)
    names(TracerElemMassShift_Array) <- rownames(Transitions)
    names(NumberTracerElemComb) <- rownames(Transitions)
    
    returnList <- list(TracerElemProb_Array = TracerElemProb_Array, TracerElemMassShift_Array = TracerElemMassShift_Array, TracerImpurityCombProb_Array = TracerImpurityCombProb_Array, 
        TracerImpurityCombMassShift_Array = TracerImpurityCombMassShift_Array, NumberTracerImpurityComb = NumberTracerImpurityComb, NumberTracerElemComb = NumberTracerElemComb)
    
    
    return(returnList)
    
}  #ElementCombinations()
