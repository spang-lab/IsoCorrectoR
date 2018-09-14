# COMPUTE TRACER IMPURITY COMBINATIONS

# Tracer purity is the probability that a tracer element atom in the 
# tracer substrate that should be labelled actually is labelled. 
# E.g.  a 99.0% isotopic purity 1,2-13C-Glucose has a 99.0% chance for 
# each of its carbon atom positions 1 and 2 of containing a 13C. 
# Thus, there is a 1.0% chance for each of those carbons that they are not 13C.
# Consequently, molecules that contain the tracer
# isotope due to metabolic transfer from the tracer substrate and not due to 
# natural abundance have a certain chance of contributing to labelling states
# with less tracer incorporated. This is due to the decrease in mass shift 
# associated with tracer impurity (e.g. 12C instead of 13C at a carbon position).
# In the following, for each labelling state the probability of having a 
# certain number of 'impure' tracer positions as well as the associated 
# (negative) mass shift are is calculated. 
# The maximum amount of impure tracer positions possible depends on the number
# of label present in a given transition.

IsoCombinationsTracerImpurity <- function(MoleculeFragmentData, ElementInfo, Transitions, NumberTransitions,
                                          IDTracer, LabelPresent, CalculationThreshold) { 

  ProbTracerImpurityList <- list()
  MassShiftTracerImpurityList <- list()
  NumberImpurityCombEff <- list()
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    
    tmpProbTracerImpurity_Array <- vector()
    tmpMassShiftTracerImpurity_Array <- vector()
    tmpNumberImpurityCombEff <- vector()
    
    for (ImpureTracer in 0:LabelPresent[TransitionNo]) {
      
      ProbImpurity <- choose(LabelPresent[TransitionNo], ImpureTracer) * (ElementInfo[[IDTracer]][["Tracer.purity"]]^(LabelPresent[TransitionNo] - 
                                                                                                                        ImpureTracer)) * ((1 - ElementInfo[[IDTracer]][["Tracer.purity"]])^ImpureTracer)
      
      if (ProbImpurity > CalculationThreshold) {
        tmpProbTracerImpurity_Array[ImpureTracer + 1] <- ProbImpurity
        tmpMassShiftTracerImpurity_Array[ImpureTracer + 1] <- (-1) * ImpureTracer * ElementInfo[[IDTracer]][["Tracer.isotope.mass.shift"]]
        NumberImpurityCombEff[TransitionNo] <- ImpureTracer
      } else {
        ImpureTracer <- LabelPresent + 1
      }  #ProbImpurity > CalculationThreshold
    }  #ImpureTracer
    
    ProbTracerImpurityList[[TransitionNo]] <- tmpProbTracerImpurity_Array
    MassShiftTracerImpurityList[[TransitionNo]] <- tmpMassShiftTracerImpurity_Array
  }  #TransitionNo
  
  names(ProbTracerImpurityList) <- rownames(Transitions)
  names(MassShiftTracerImpurityList) <- rownames(Transitions)
  names(NumberImpurityCombEff) <- rownames(Transitions)
  
  return(list("ProbTracerImpurityList"=ProbTracerImpurityList,
              "MassShiftTracerImpurityList"=MassShiftTracerImpurityList,
              "NumberImpurityCombEff"=NumberImpurityCombEff))
 
}