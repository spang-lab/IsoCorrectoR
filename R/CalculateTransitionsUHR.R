# Automated Computation of Transitions

#' @importFrom stringr str_detect
#' @importFrom stringr str_c
#' @importFrom magrittr '%>%'
#' 
  CalculateTransitionsUHR <- function(MoleculeInfo, ElementInfo, verbose) {
  
  if(verbose){message(date(), " :: calculating transitions ...")}
  MoleculesTotal <- length(MoleculeInfo)
  MoleculesName <- names(MoleculeInfo)
    
  # CalculateTransitions() in the case of UltraHighRes option active for n different tracers.
  
  TotalCombinationsList <- list()
  TransitionsExpectedList <- list()
  for (MoleculeNo in seq_len(MoleculesTotal)) {
    
    MoleculeData <- MoleculeInfo[[MoleculeNo]]  # entire 'MoleculeInfo' data for current molecule
    NumberTracers <- unlist(lapply(MoleculeData, function(x) length(x[["Tracer"]])))[[1]]  # number of tracers of First Fragment only !
    
    MaxLabel <- lapply(MoleculeData, function(x) x[["MaxLabel"]])[[1]]
    IDTracer <- lapply(MoleculeData, function(x) x[["IDTracer"]])[[1]]
    
    tmpTotalCombinationsList <- list()
    
    # SubCombinations.vec[TracerNo] contains the number of possible combined labelling states that can arise when combining the labelling of the
    # [1...TracerNo-1] tracer elements that occur in the list of tracer elements prior to the tracer element with number TracerNo.
    
    SubCombinations.vec <- 1
    
    if (NumberTracers > 0) 
    {
      tmpTransitionsExpectedList <- list()
      # number of tracers
      
      for (TracerNo in seq_len(NumberTracers)) {
        SubCombinations.vec[TracerNo + 1] <- c(SubCombinations.vec[TracerNo] * (MaxLabel[TracerNo] + 1))
      }  #TracerNo 
      
      # TotalCombinations is the total number of labelling states that is possible with the set of tracer elements
      
      TotalCombinations <- SubCombinations.vec[NumberTracers[[1]] + 1]
      
      # Now for each expected labelling state the number of each tracer element is calculated
      
      TransitionsExpected.vec <- vector()
      for (TracerNo in seq_len(NumberTracers)) {
        
        LabelSpacer <- 1
        
        # In this loop, each number of label of a tracer is repeated 
        # SubCombinations.vec[TracerNo] times in succession in the 
        # TransitionsExpected.vec. LabelSpacer
        # then provides the first index of the next number of label in the vector
        
        for (Label in 0:MaxLabel[TracerNo]) {
          
          TransitionsExpected.vec[LabelSpacer:(LabelSpacer + SubCombinations.vec[TracerNo] - 1)] <- Label
          LabelSpacer <- LabelSpacer + SubCombinations.vec[TracerNo]
          
        }  #Label
        
        tmpTransitionsExpectedList[[TracerNo]] <- TransitionsExpected.vec
        
        RepeatSpacer <- LabelSpacer
        
        # Now, the block generated in TransitionsExpected.vec in the previous
        # loop is copied until the length of the vector matches 
        # TotalCombinations. Then, the
        # labelling with the given tracer is present in 
        # TransitionsExpected.vec for all labelling states of the molecule.
        
        for (Repeat in seq_len((TotalCombinations/SubCombinations.vec[TracerNo + 1])) - 1) {
          
          if (Repeat > 0) {
            TransitionsExpected.vec[RepeatSpacer:(RepeatSpacer + SubCombinations.vec[TracerNo + 1] - 1)] <- tmpTransitionsExpectedList[[TracerNo]][seq_len(SubCombinations.vec[TracerNo + 1])]
            
            RepeatSpacer <- RepeatSpacer + SubCombinations.vec[TracerNo + 1]
          }
        }  #Repeat
        
        # This procedure is performed for all tracers to get a list of 
        # expected labelling states
        
        tmpTransitionsExpectedList[[TracerNo]] <- TransitionsExpected.vec
      }  #TracerNo
      
      names(tmpTransitionsExpectedList) <- MoleculeData[[1]][["IDTracer"]]
      
      NumberTransitions <- length(tmpTransitionsExpectedList[[1]])
      
      tmpTransitionsExpectedList[["Sum"]] <- apply(as.data.frame(tmpTransitionsExpectedList), 1, sum)
      
      TransitionsExpected.df <- as.data.frame(tmpTransitionsExpectedList)
      
      # Compute the expected names of the expected labelling states
      
      Prefix <- stringr::str_c(MoleculesName[MoleculeNo], "_")
      
      for (TransitionNo in seq_len(NumberTransitions)) {
        PreName <- Prefix
        
        for (TracerNo in seq_len(NumberTracers)) {
          
          PreName <- stringr::str_c(PreName, IDTracer[TracerNo], as.numeric(TransitionsExpected.df[TransitionNo, TracerNo]))
          if (TracerNo < NumberTracers) 
          {
            PreName <- stringr::str_c(PreName, ".")
          }  #if
          if (TracerNo == NumberTracers) 
          {
            rownames(TransitionsExpected.df)[TransitionNo] <- PreName
          }  #if
          
        }  #TracerNo
        
      }  #TransitionNo
      
      TransitionsExpected.df <- TransitionsExpected.df[order(rownames(TransitionsExpected.df)), ]
      
      MoleculeInfo[[MoleculeNo]][["TransitionsExpected"]] <- TransitionsExpected.df
      MoleculeInfo[[MoleculeNo]][["TotalCombinations"]] <- TotalCombinations
      
    }  #NumberTracers[[1]]>0
    
  }  #MoleculeNo
  
  if(verbose){message(date(), " :: calculating transitions [OK]\n")}
  return(MoleculeInfo)

}
