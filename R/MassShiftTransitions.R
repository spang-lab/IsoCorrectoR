# Generation of a mass shift vector of labelling states

# The purpose of the function 'MassShiftTransitions' is
# to generate a vector (FragmentMassShift) which contains
# the mass shift of a molecule(-fragment) of a given
# transition in relation to the same molecule(-fragment)
# containing no isotopes of higher mass. This is required to
# be able to yield the probability that a certain transition
# or full MS ion x is derived from another transition/full MS ion
# y due to natural stable isotope abundance. This probability is yielded
# by matching the mass shift k of transition x in FragmentMassShift
# with the mass shift k of transition y in CumProbList
# (this is the mass shift 'produced' by y through natural abundance).
# The probability is found associated with mass shift k in CumProbList.

MassShiftTransitions <- function(MoleculeInfo, MoleculeNo, Fragment) {
  
  Transitions <- MoleculeInfo[[MoleculeNo]][["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  
  FragmentMassShift <- vector()
  
  for (TransitionNo in seq_len(NumberTransitions)) {
    FragmentMassShift[TransitionNo] <- Transitions[TransitionNo, Fragment]
  }  #TransitionNo
  
  return(FragmentMassShift)
  
}  #MassShiftTransitions
