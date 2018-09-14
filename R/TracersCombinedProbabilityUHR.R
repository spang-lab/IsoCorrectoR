# Calculate combined tracer probabilites for UHR

# In a way similar to the function FragmentsCombinedProbability in the normal 
# resolution mode, this section gives the probability that a given labeling state
# y produces the same isotope pattern as a labeling state x by natural abundance 
# (and tracer impurity). To this end, the probability that the number of tracer isotopes in 
# transition x is produced by transition y for each tracer element is 
# multiplied with the other tracer element probabilities. 
# The variable CalculationThreshold_UHR constraints
# the transitions y for which the contribution to x is calculated to 
# CalculationThreshold_UHR which is +/- the total amount of 
# additional/substracted label
# (all tracers summed up) compared to x.

TracersCombinedProbabilityUHR <- function(MoleculeData,CumProbList,CalculationThreshold_UHR) {

  NumberTracers <- length(MoleculeData[[1]][["Tracer"]])
  
  Transitions <- MoleculeData[["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  
  CombinedProb <- matrix(ncol = NumberTransitions, nrow = NumberTransitions)
  
  for (x in seq_len(NumberTransitions)) {
    if (CalculationThreshold_UHR > 0) {
      MultipleTracerProbVector <- vector()
      for (y in seq_len(NumberTransitions)) {
        if (((Transitions[x, NumberTracers + 1] - Transitions[y, NumberTracers + 1])^2)^0.5 <= CalculationThreshold_UHR) {
          Prob <- 1
          for (TracerNo in seq_len(NumberTracers)) {
            Prob <- Prob * CumProbList[[TracerNo]][[Transitions[y, TracerNo] + 1]][Transitions[x, TracerNo] + 1]
          }  #TracerNo
        } else {
          Prob <- 0
        }
        MultipleTracerProbVector[y] <- Prob
      }  #y
      CombinedProb[x, ] <- MultipleTracerProbVector
      rm(MultipleTracerProbVector)
    } else {
      for (y in seq_len(NumberTransitions)) {
        Prob <- 1
        for (TracerNo in seq_len(NumberTracers)) {
          Prob <- Prob * CumProbList[[TracerNo]][[Transitions[y, TracerNo] + 1]][Transitions[x, TracerNo] + 1]
        }  #TracerNo
        MultipleTracerProbVector[y] <- Prob
      }  #y
      CombinedProb[x, ] <- MultipleTracerProbVector
      rm(MultipleTracerProbVector)
    }  #if(CalculationThreshold_UHR>0)
  }  #x
  colnames(CombinedProb) <- rownames(Transitions)
  rownames(CombinedProb) <- rownames(Transitions)
 
  return(CombinedProb)
   
}