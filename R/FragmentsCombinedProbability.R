# Calculation of probability that a given labelling state produces the mass shift of another labelling state

#' @importFrom stringr str_detect
#' @importFrom magrittr '%>%'
#' 
FragmentsCombinedProbability <- function(MoleculeArray, MoleculeNo, CumProbList, FragmentMassShift) {
  
  NumberFragments <- stringr::str_detect(names(MoleculeArray[[MoleculeNo]]), "Fragment") %>% sum
  Transitions <- MoleculeArray[[MoleculeNo]][["Transitions"]]
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
