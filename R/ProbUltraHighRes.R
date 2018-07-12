# Calculation of the probability matrix for Ultra High Resolution (UHR) 
# data correction

ProbUltraHighRes <- function(MoleculeData, ElementArray, LabelDiffConstraint, CorrectTracerElementCore, CorrectTracerImpurity) {
  
  message(date(), " :: ProbUltraHighRes")
  message(date(), " :: CorrectTracerElementCore: ", CorrectTracerElementCore)
  message(date(), " :: CorrectTracerImpurity: ", CorrectTracerImpurity)
  message(date(), " :: LabelDiffConstraint: ", LabelDiffConstraint)
  
  NumberTracers <- length(MoleculeData[[1]][["Tracer"]])
  nTracerMax <- MoleculeData[[1]][["nTracerMax"]]
  IDTracer <- MoleculeData[[1]][["IDTracer"]]
  MaxLabel <- MoleculeData[[1]][["MaxLabel"]]
  NatAbuTracer <- MoleculeData[[1]][["NatAbuTracer"]]
  NatAbuBase <- MoleculeData[[1]][["NatAbuBase"]]
  AvailableTracerPlaces <- vector()
  MaxTracerPlaces <- vector()
  variableTracerMaximum <- list()
  
  Transitions <- MoleculeData[["Transitions"]]
  NumberTransitions <- nrow(Transitions)
  
  NatAbuProbList <- list()
  NatAbuShiftList <- list()
  for (TracerNo in seq_len(NumberTracers)) {
    
    if (CorrectTracerElementCore) {
      AvailableTracerPlaces[TracerNo] <- nTracerMax[TracerNo]
      MaxTracerPlaces[TracerNo] <- MaxLabel[TracerNo]
    } else if (!CorrectTracerElementCore) {
      AvailableTracerPlaces[TracerNo] <- nTracerMax[TracerNo] - MaxLabel[TracerNo]
      if (MaxLabel[TracerNo] <= AvailableTracerPlaces[TracerNo]) {
        MaxTracerPlaces[TracerNo] <- MaxLabel[TracerNo]
      } else {
        MaxTracerPlaces[TracerNo] <- AvailableTracerPlaces[TracerNo]
      }
    } else {
      stop(date(), " :: [STOP] invalid `CorrectTracerElementCore` value: ", CorrectTracerElementCore)
    }
    tmpNatAbuProbList <- list()
    tmpNatAbuShiftList <- list()
    
    tmp_variableTracerMaximum.vec <- vector(mode = "numeric", length = MaxLabel[TracerNo] + 1)
    
    for (IntrinsicLabel in 0:MaxLabel[TracerNo]) {
      
      if (CorrectTracerElementCore) {
        x <- IntrinsicLabel
        
      } else if (!CorrectTracerElementCore) {
        x <- 0
      }
      NatAbuProb <- vector()
      NatAbuShift <- vector()
      
      if (!CorrectTracerImpurity) {
        tmp_variableTracerMaximum <- MaxTracerPlaces[TracerNo] - x
      } else if (CorrectTracerImpurity && (MaxTracerPlaces[TracerNo] <= (AvailableTracerPlaces[TracerNo] - x))) {
        tmp_variableTracerMaximum = MaxTracerPlaces[TracerNo]
      } else {
        tmp_variableTracerMaximum = AvailableTracerPlaces[TracerNo] - x
      }
      
      for (NatAbuLabel in 0:tmp_variableTracerMaximum) {
        
        NatAbuProb[NatAbuLabel + 1] <- choose(AvailableTracerPlaces[TracerNo] - x, NatAbuLabel) * NatAbuTracer[TracerNo]^NatAbuLabel * NatAbuBase[TracerNo]^(AvailableTracerPlaces[TracerNo] - 
                                                                                                                                                               NatAbuLabel - x)
        
        NatAbuShift[NatAbuLabel + 1] <- IntrinsicLabel + NatAbuLabel
      }  #NatAbuLabel
      
      tmpNatAbuProbList[[IntrinsicLabel + 1]] <- NatAbuProb
      tmpNatAbuShiftList[[IntrinsicLabel + 1]] <- NatAbuShift
      
      tmp_variableTracerMaximum.vec[x + 1] <- tmp_variableTracerMaximum
      
    }  #IntrinsicLabel
    
    variableTracerMaximum[[TracerNo]] <- tmp_variableTracerMaximum.vec
    
    NatAbuProbList[[TracerNo]] <- tmpNatAbuProbList
    NatAbuShiftList[[TracerNo]] <- tmpNatAbuShiftList
    
  }  #TracerNo
  
  # As in the case of the normal resolution approach, a correction for tracer 
  # impurity can be performed, in this case for all tracers present in the
  # molecule. The amount of impure tracer that may be found depends on the 
  # amount of intrinsic label of this tracer in the molecule.
  
  if (CorrectTracerImpurity) {
    ImpureTracerProbList <- list()
    for (TracerNo in seq_len(NumberTracers)) {
      tmpImpureTracerProbList <- list()
      for (IntrinsicLabel in 0:MaxLabel[TracerNo]) {
        ImpureTracerProb <- vector()
        for (ImpureTracer in 0:IntrinsicLabel) {
          ImpureTracerProb[ImpureTracer + 1] <- choose(IntrinsicLabel, ImpureTracer) * (1 - ElementArray[[IDTracer[TracerNo]]][[3]])^ImpureTracer * 
            ElementArray[[IDTracer[TracerNo]]][[3]]^(IntrinsicLabel - ImpureTracer)
        }  #ImpureTracer
        tmpImpureTracerProbList[[IntrinsicLabel + 1]] <- ImpureTracerProb
      }  #IntrinsicLabel
      ImpureTracerProbList[[TracerNo]] <- tmpImpureTracerProbList
    }  #TracerNo
    
    # In the following section, natural abundance and tracer impurity 
    # probabilities are combined multiplicatively. 
    # The total shift in tracer number associated
    # with natural abundance and tracer impurity (negative shift) is derived by 
    # adding/substracting all individual shifts.
    
    NatAbuImpurityProbList <- list()
    NatAbuImpurityShiftList <- list()
    for (TracerNo in seq_len(NumberTracers)) {
      
      tmpIntrinsicLabelList <- list()
      tmpIntrinsicLabelShiftList <- list()
      tmp_variableTracerMaximum.vec <- variableTracerMaximum[[TracerNo]]
      
      for (IntrinsicLabel in 0:MaxLabel[TracerNo]) {
        message(date(), " :: IntrinsicLabel:", IntrinsicLabel)
        if (CorrectTracerElementCore) {
          x <- IntrinsicLabel
        } else if (!CorrectTracerElementCore) {
          x <- 0
        }
        tmpNatAbuLabelList <- list()
        tmpNatAbuLabelShiftList <- list()
        
        tmp_variableTracerMaximum <- tmp_variableTracerMaximum.vec[x + 1]
        
        for (NatAbuLabel in 0:tmp_variableTracerMaximum) {
          message(date(), " :: NatAbuLabel:", NatAbuLabel)
          NatAbuImpurityProb <- vector()
          NatAbuImpurityShift <- vector()
          for (ImpureTracer in 0:IntrinsicLabel) {
            message(date(), " :: ImpureTracer:", ImpureTracer)
            NatAbuImpurityProb[ImpureTracer + 1] <- NatAbuProbList[[TracerNo]][[IntrinsicLabel + 1]][NatAbuLabel + 1] * ImpureTracerProbList[[TracerNo]][[IntrinsicLabel + 
                                                                                                                                                            1]][ImpureTracer + 1]
            NatAbuImpurityShift[ImpureTracer + 1] <- IntrinsicLabel + NatAbuLabel - ImpureTracer
          }
          tmpNatAbuLabelList[[NatAbuLabel + 1]] <- NatAbuImpurityProb
          tmpNatAbuLabelShiftList[[NatAbuLabel + 1]] <- NatAbuImpurityShift
        }  #NatAbuLabel
        tmpIntrinsicLabelList[[IntrinsicLabel + 1]] <- tmpNatAbuLabelList
        tmpIntrinsicLabelShiftList[[IntrinsicLabel + 1]] <- tmpNatAbuLabelShiftList
      }  #IntrinsicLabel
      
      NatAbuImpurityProbList[[TracerNo]] <- tmpIntrinsicLabelList
      NatAbuImpurityShiftList[[TracerNo]] <- tmpIntrinsicLabelShiftList
      
    }  #TracerNo
  } else {
    NatAbuImpurityProbList <- NatAbuProbList
    NatAbuImpurityShiftList <- NatAbuShiftList
  }  #CorrectTracerImpurity==TRUE
  
  # Similarly to the function MassShiftProbabilities in the normal resolution 
  # mode, this part of the code finds all probabilites of a certain tracer in 
  # the probability list NatAbuImpurityProbList that are associated with a total 
  # number of tracer isotope TotalLabel.
  # It then sums them up to yield CumProbList.
  
  CumProbList <- list()
  for (TracerNo in seq_len(NumberTracers)) {
    message(date(), " :: [TracerNo] ", TracerNo)
    
    tmpCumProbList <- list()
    for (IntrinsicLabel in 0:MaxLabel[TracerNo]) {
      message(date(), " :: [IntrinsicLabel] ", IntrinsicLabel)
      tmpCumProb.vec <- vector()
      for (TotalLabel in 0:MaxLabel[TracerNo]) {
        message(date(), " :: [TotalLabel] ", TotalLabel)
        
        tmp.idx <- lapply(NatAbuImpurityShiftList[[TracerNo]][[IntrinsicLabel + 1]], function(x) x == TotalLabel)
        
        tmpProb <- unlist(NatAbuImpurityProbList[[TracerNo]][[IntrinsicLabel + 1]])[unlist(tmp.idx)]
        
        tmpCumProb.vec[TotalLabel + 1] <- sum(unlist(lapply(tmpProb, sum)))
        
      }  #TotalLabel
      
      tmpCumProbList[[IntrinsicLabel + 1]] <- tmpCumProb.vec
    }  #IntrinsicLabel
    
    CumProbList[[TracerNo]] <- tmpCumProbList
    
  }  #TracerNo
  
  # In a way similar to the function FragmentsCombinedProbability in the normal 
  # resolution mode, this section gives the probability that a given transition
  # y produces a transition x by natural abundance (and tracer impurity). 
  # To this end, the probability that the number of tracer isotopes in 
  # transition x is produced by transition y for each tracer element is 
  # multiplied with the other tracer element probabilities. 
  # The variable LabelDiffConstraint constraints
  # the transitions y for which the contribution to x is calculated to 
  # LabelDiffConstraint which is +/- the total amount of 
  # additional/substracted label
  # (all tracers summed up) compared to x.
  
  CombinedProb <- matrix(ncol = NumberTransitions, nrow = NumberTransitions)
  for (x in seq_len(NumberTransitions)) {
    if (LabelDiffConstraint > 0) {
      MultipleTracerProbVector <- vector()
      for (y in seq_len(NumberTransitions)) {
        if (((Transitions[x, NumberTracers + 1] - Transitions[y, NumberTracers + 1])^2)^0.5 <= LabelDiffConstraint) {
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
    }  #if(LabelDiffConstraint>0)
  }  #x
  colnames(CombinedProb) <- rownames(Transitions)
  rownames(CombinedProb) <- rownames(Transitions)
  
  message(date(), " :: leaving ProbUltraHighRes()")
  return(CombinedProb)
}  # function
