# Calculate probabilites and mass shifts of element combinations

# With the probabilities and mass shifts of the different
# isotope combinations calculated, the next step is to combine
# the isotope combination probabilities and mass shifts of
# the different elements (and the tracer impurity, if checked) using the function
# 'ElementCombinations'. For a given molecule(-fragment) the function combinatorically
# multiplies all isotope combination probabilities of all non-tracer elements,
# the tracer element and the tracer impurity with each other
# In the same way, the isotope combination mass shifts are added to each
# other. In the end, the probability that a molecule(-fragment) contains a combination
# of specific element isotope combinations for its various elements
# (and a certain number of 'impure' tracer atoms) is yielded,
# together with the mass shift associated with
# such a combination. Because combinatorically multiplying all
# isotope combinations is resource intensive, the
# calculation stops at certain probability thresholds defined by the
# parameter CalculationThreshold. In those cases the probability of a combination is so low
# that it does not affect the correction.

ElementCombinations <- function(MoleculeInfo, MoleculeNo, Fragment, IsoCombinationsResult, CalculationThreshold,
    CorrectTracerImpurity) {
    
    MoleculeFragmentData <- MoleculeInfo[[MoleculeNo]][[str_c("Fragment_", Fragment)]]
    ProbTracerList <- IsoCombinationsResult[["ProbTracerList"]]
    MassShiftTracerList <- IsoCombinationsResult[["MassShiftTracerList"]]
    NumberIsoCombTracerEff <- IsoCombinationsResult[["NumberIsoCombTracerEff"]]
    ProbElemList <- IsoCombinationsResult[["ProbElemList"]]
    MassShiftElemList <- IsoCombinationsResult[["MassShiftElemList"]]
    NumberIsoCombEff <- IsoCombinationsResult[["NumberIsoCombEff"]]
    ProbTracerImpurityList <- IsoCombinationsResult[["ProbTracerImpurityList"]]
    MassShiftTracerImpurityList <- IsoCombinationsResult[["MassShiftTracerImpurityList"]]
    NumberImpurityCombEff <- IsoCombinationsResult[["NumberImpurityCombEff"]]
    ElementsNonTracer <- MoleculeFragmentData[["NonTracer"]]
    
    if (length(ElementsNonTracer) == 0) {
        ElementsNonTracer <- MoleculeFragmentData[["ZeroTracer"]]
    }
    
    NumberElementsNonTracer <- length(ElementsNonTracer)
    
    # CALCULATION OF THE NON-TRACER ELEMENT COMBINATIONS
    
    ElementCombinationsNonTracerResult <- ElementCombinationsNonTracer(ProbElemList=ProbElemList,
                                 MassShiftElemList=MassShiftElemList,
                                 NumberIsoCombEff=NumberIsoCombEff,
                                 ElementsNonTracer=ElementsNonTracer,
                                 NumberElementsNonTracer = NumberElementsNonTracer,
                                 CalculationThreshold=CalculationThreshold)
    
    CombElemProbList <- ElementCombinationsNonTracerResult$CombElemProbList
    CombElemMassShiftList <- ElementCombinationsNonTracerResult$CombElemMassShiftList
    NumberEleComb <- ElementCombinationsNonTracerResult$NumberEleComb
    
    # CALCULATION OF THE TRACER ELEMENT - NON-TRACER ELEMENT COMBINATIONS
    
    NumberTracers <- length(MoleculeFragmentData[["Tracer"]])
    Transitions <- MoleculeInfo[[MoleculeNo]][["Transitions"]]
    NumberTransitions <- nrow(Transitions)
    
    ElementCombinationsTracerResult <- ElementCombinationsTracer(ProbTracerList=ProbTracerList,
                                                                       MassShiftTracerList=MassShiftTracerList, 
                                                                       NumberIsoCombTracerEff=NumberIsoCombTracerEff, 
                                                                       CombElemProbList=CombElemProbList, 
                                                                       CombElemMassShiftList=CombElemMassShiftList,
                                                                       NumberEleComb=NumberEleComb,
                                                                       Fragment=Fragment,
                                                                       NumberElementsNonTracer=NumberElementsNonTracer,
                                                                       NumberTracers=NumberTracers, 
                                                                       Transitions=Transitions,
                                                                       NumberTransitions=NumberTransitions, 
                                                                       CalculationThreshold=CalculationThreshold)
    
    TracerElemProbList <- ElementCombinationsTracerResult$TracerElemProbList
    TracerElemMassShiftList <- ElementCombinationsTracerResult$TracerElemMassShiftList
    NumberTracerElemComb <- ElementCombinationsTracerResult$NumberTracerElemComb
    
    # COMBINATION WITH TRACER IMPURITY STATES
    
    if (CorrectTracerImpurity && (NumberTracers > 0)) {
      
      ElementCombinationsTracerImpurityResult <- ElementCombinationsTracerImpurity(ProbTracerImpurityList=ProbTracerImpurityList,
                                                                                   MassShiftTracerImpurityList=MassShiftTracerImpurityList,
                                                                                   NumberImpurityCombEff=NumberImpurityCombEff,
                                                                                   TracerElemProbList=TracerElemProbList,
                                                                                   TracerElemMassShiftList=TracerElemMassShiftList,
                                                                                   NumberTracerElemComb=NumberTracerElemComb,
                                                                                   Transitions=Transitions,
                                                                                   NumberTransitions=NumberTransitions,
                                                                                   CalculationThreshold=CalculationThreshold) 
      
      TracerImpurityCombProbList <- ElementCombinationsTracerImpurityResult$TracerImpurityCombProbList
      TracerImpurityCombMassShiftList <- ElementCombinationsTracerImpurityResult$TracerImpurityCombMassShiftList
      NumberTracerImpurityComb <- ElementCombinationsTracerImpurityResult$NumberTracerImpurityComb

    } else {
      
      TracerImpurityCombProbList <- list()
      TracerImpurityCombMassShiftList <- list()
      NumberTracerImpurityComb <- vector()
      
    } 

    return(list(TracerElemProbList = TracerElemProbList, TracerElemMassShiftList = TracerElemMassShiftList, TracerImpurityCombProbList = TracerImpurityCombProbList, 
        TracerImpurityCombMassShiftList = TracerImpurityCombMassShiftList, NumberTracerImpurityComb = NumberTracerImpurityComb, NumberTracerElemComb = NumberTracerElemComb))
    
}  #ElementCombinations()
