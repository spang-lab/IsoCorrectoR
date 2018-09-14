# Calculation of probability that a given labelling state produces the mass shift of another labelling state

# In the function 'FragmentsCombinedProbability' the final
# probability matrix CombFragmentProb is produced. In this matrix, for each
# transition y of a molecule the probability that it produces a transition x
# of higher mass through natural isotope abundance is given.
# This is done by matching the mass shifts from
# FragmentMassShift and CumProbList, as described before for the
# function 'MassShiftTransitions'. For MS/MS data, the
# probability that the product ion y produces the product ion
# mass shift of transition x and and the probability that neutral
# loss y produces the neutral loss mass shift of transition x
# are multiplied. This way the overall probability that
# transition y produces transition x through natural isotope abundance
# is yielded.

#' @importFrom stringr str_detect
#' @importFrom magrittr '%>%'
#' 
FragmentsCombinedProbability <- function(MoleculeInfo, MoleculeNo, CumProbList, FragmentMassShift) {
  
  NumberFragments <- stringr::str_detect(names(MoleculeInfo[[MoleculeNo]]), "Fragment") %>% sum
  Transitions <- MoleculeInfo[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  
  CombinedProb <- matrix(NA, nrow = NumberTransitions, ncol = NumberTransitions)
  colnames(CombinedProb) <- rownames(Transitions)
  rownames(CombinedProb) <- rownames(Transitions)
  
  if (NumberFragments == 2) {
    for (x in seq_len(NumberTransitions)) {
      for (y in seq_len(NumberTransitions)) {
        CombinedProb[x, y] <- as.data.frame(CumProbList[[1]])[y, FragmentMassShift[[1]][x] + 1] * as.data.frame(CumProbList[[2]])[y, FragmentMassShift[[2]][x] + 1]
      }  #y
    }  #x
  } else if (NumberFragments == 1) {
    for (x in seq_len(NumberTransitions)) {
      for (y in seq_len(NumberTransitions)) {
        CombinedProb[x, y] <- as.data.frame(CumProbList[[1]])[y, FragmentMassShift[[1]][x] + 1]
      }  #y
    }  #x
  }
  
  return(CombinedProb)
  
}
