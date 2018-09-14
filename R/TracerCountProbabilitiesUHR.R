# Calculate probabilites for net tracer counts for UHR

# Similarly to the function MassShiftProbabilities in the normal resolution 
# mode, this part of the code finds all probabilites of a certain tracer in 
# the probability list NatAbuImpurityProbList that are associated with a total 
# number of tracer isotope TotalLabel.
# It then sums them up to yield CumProbList.

TracerCountProbabilitiesUHR <- function(MoleculeData,NatAbuImpurityList, verbose) {

  NumberTracers <- length(MoleculeData[[1]][["Tracer"]])
  MaxLabel <- MoleculeData[[1]][["MaxLabel"]]
  
  NatAbuImpurityShiftList <- NatAbuImpurityList$NatAbuImpurityShiftList
  NatAbuImpurityProbList <- NatAbuImpurityList$NatAbuImpurityProbList
  
  CumProbList <- list()
  for (TracerNo in seq_len(NumberTracers)) {
    if(verbose){message(date(), " :: [TracerNo] ", TracerNo)}
    
    tmpCumProbList <- list()
    for (IntrinsicLabel in 0:MaxLabel[TracerNo]) {
      if(verbose){message(date(), " :: [IntrinsicLabel] ", IntrinsicLabel)}
      tmpCumProb.vec <- vector()
      for (TotalLabel in 0:MaxLabel[TracerNo]) {
        if(verbose){message(date(), " :: [TotalLabel] ", TotalLabel)}
        
        tmp.idx <- lapply(NatAbuImpurityShiftList[[TracerNo]][[IntrinsicLabel + 1]], function(x) x == TotalLabel)
        
        tmpProb <- unlist(NatAbuImpurityProbList[[TracerNo]][[IntrinsicLabel + 1]])[unlist(tmp.idx)]
        
        tmpCumProb.vec[TotalLabel + 1] <- sum(unlist(lapply(tmpProb, sum)))
        
      }  #TotalLabel
      
      tmpCumProbList[[IntrinsicLabel + 1]] <- tmpCumProb.vec
    }  #IntrinsicLabel
    
    CumProbList[[TracerNo]] <- tmpCumProbList
    
  }  #TracerNo
  
  return(CumProbList)
  
}