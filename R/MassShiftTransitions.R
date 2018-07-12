# Generation of a mass shift array of labelling states

MassShiftTransitions <- function(MoleculeArray, MoleculeNo, Fragment) {
  
  Transitions <- MoleculeArray[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  
  FragmentMassShift <- vector()
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    FragmentMassShift[TransitionNo] <- Transitions[TransitionNo, Fragment]
  }  #TransitionNo
  
  return(FragmentMassShift)
  
}  #MassShiftTransitions
