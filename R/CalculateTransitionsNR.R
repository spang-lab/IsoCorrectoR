CalculateTransitionsNR <- function(MoleculeInfo, ElementInfo, verbose) {
  
  if(verbose){message(date(), " :: calculating transitions ...")}
  MoleculesTotal <- length(MoleculeInfo)
  MoleculesName <- names(MoleculeInfo)
  
  # For all molecules, determine the maximum amount of tracer isotopes expected
  # in product ion or neutral loss.
  for (MoleculeNo in seq_len(MoleculesTotal)) {
    
    ### assemble required data
    MoleculeData <- MoleculeInfo[[MoleculeNo]]  # entire 'MoleculeInfo' data for current molecule
    
    NumberFragments <- sum(stringr::str_detect(names(MoleculeData), "Fragment"))  # %>%sum
    NumberTracers <- unlist(lapply(MoleculeData, function(x) length(x[["Tracer"]])))
    MaxLabel <- unlist(lapply(MoleculeData, function(x) x[["MaxLabel"]]))
    IDTracer <- unlist(lapply(MoleculeData, function(x) x[["IDTracer"]]))
    
    if (NumberFragments == 2) {
      
      #if(verbose){message(date(), " :: [", MoleculesName[MoleculeNo], "] // NumberFragments==", NumberFragments)}
      if(verbose){message(date(), " :: :: molecule [", MoleculesName[MoleculeNo], "] with ",NumberFragments," fragments")}
      if (NumberTracers[[1]] > 0) {
        MaxLabelProduct <- MaxLabel[[1]]
      } else {
        MaxLabelProduct <- 0
      }
      
      if (NumberTracers[[2]] > 0) {
        MaxLabelNeutralLoss <- MaxLabel[[2]]
      } else {
        MaxLabelNeutralLoss <- 0
      }
      if (NumberTracers[[1]] > 0) {
        Tracer <- IDTracer[[1]]
      } else if (NumberTracers[[2]] > 0) {
        Tracer <- IDTracer[[2]]
      } else {
        stop(date(), " :: [CalculateTransitions()] No tracer element specified for molecule ", MoleculesName[MoleculeNo], "!")
      }
      
      MaxLabelPrecursor <- MaxLabelProduct + MaxLabelNeutralLoss
      
      # CALCULATION OF MS/MS TRANSITIONS
      
      TransitionNo <- 0
      LabelNeutralLoss <- 0
      ProductIon_SiteSaturation <- 0
      SiteSaturationCorrection <- 0
      
      # This part of the function computes how label in the precursor ion can be distributed among product ion and neutral loss, within the constraints given by
      # MaxLabelProduct and MaxLabelNeutralLoss. This yields the expected MS/MS transitions.
      
      tmpTransitionsExpectedList <- list()
      tmpTransitionsExpectedListNames.vec <- vector()
      
      for (LabelPrecursor in 0:MaxLabelPrecursor) {
        # Until the label in the precursor has exceeded the amount of MaxLabelNeutralLoss, all label is considered to be in the neutral loss at this stage.
        
        if (LabelNeutralLoss < MaxLabelNeutralLoss) {
          LabelNeutralLoss <- LabelPrecursor
        }
        
        # ProductIon_SiteSaturation is a means to determine by how much the label in the precursor exceeds the maximum capacity of the product ion given in
        # MaxLabelProduct.
        
        ProductIon_SiteSaturation <- (MaxLabelProduct - LabelPrecursor) * (-1)
        if (ProductIon_SiteSaturation > 0) {
          SiteSaturationCorrection <- ProductIon_SiteSaturation
        }
        
        # At this stage, the MS/MS transitions are calculated and written into the molecule information list.  The vector named Precursor contains the mass shift
        # (in relation to M+0) associated with the precursor (Number of label in the precursor multiplied with the tracer isotope mass shift given in
        # ElementList).  The vector named NeutralLoss contains the mass shift of the neutral loss. This is determined in a loop which starts at the maximum
        # possible amount of current LabelPrecursor from the parent loop in the neutral loss. It then decreases the amount of LabelNeutralLoss by 'shifting' label
        # from the neutral loss to the product ion (NLLabel_in_ProductIon). Consequently, the product ion mass shift in the vector named ProductIon is given as
        # the difference between precursor mass shift and neutral loss mass shift of the current loop iteration.  This way all combinations of product ion and
        # neutral loss labelling for a given precursor labelling state are derived. The shifting just described has limits given by MaxLabelProduct.  This is
        # where SiteSaturationCorrection becomes relevant.  If the amount of label in the precursor exceeds the maximum amount of labelling possible in the
        # product ion, it will reduce the number of iterations of the shifting loop by just this difference. In addition, the loop generates specific name tags
        # for each transition of a given molecule by combining the molecule name from the molecules information file with '_T' for transition, followed by the
        # precuror label x and product ion label y as x.y (Name_x.y).
        
        # TransitionNo<-0
        for (NLLabel_in_ProductIon in 0:(LabelNeutralLoss - SiteSaturationCorrection)) {
          
          TransitionNo <- TransitionNo + 1
          
          tmpTrans1 <- LabelPrecursor * ElementInfo[[Tracer]][[2]]
          tmpTrans2 <- (LabelNeutralLoss - NLLabel_in_ProductIon) * ElementInfo[[Tracer]][[2]]
          tmpTrans3 <- tmpTrans1 - tmpTrans2
          
          tmpTransitionsExpectedList[[TransitionNo]] <- list(ProductIon = tmpTrans3, NeutralLoss = tmpTrans2, Precursor = tmpTrans1)
          
          tmpTransitionsExpectedListNames.vec <- c(tmpTransitionsExpectedListNames.vec, str_c(MoleculesName[MoleculeNo], "_", tmpTrans1, ".", tmpTrans3))
          
        }  # NLLabel
        
      }  # LabelPrecursor
      
      names(tmpTransitionsExpectedList) <- tmpTransitionsExpectedListNames.vec
      
      TransitionsExpected.df <- data.frame(ProductIon = as.numeric(unlist(lapply(tmpTransitionsExpectedList, function(x) x["ProductIon"]))), NeutralLoss = as.numeric(unlist(lapply(tmpTransitionsExpectedList, 
                                                                                                                                                                                    function(x) x["NeutralLoss"]))), Precursor = as.numeric(unlist(lapply(tmpTransitionsExpectedList, function(x) x["Precursor"]))))
      
      rownames(TransitionsExpected.df) <- tmpTransitionsExpectedListNames.vec
      
      TransitionsExpected <- TransitionsExpected.df
      
      # In this section the expected MS (not MS/MS) measurements and their name tags are generated for each molecule An MS name tag has the structure 'Name_x'
      
    } else if (NumberFragments == 1) 
    {
      
      if (NumberTracers[[1]] > 0) {
        MaxLabelProduct <- MaxLabel[[1]]
        Tracer <- IDTracer[[1]]
      } else {
        MaxLabelProduct <- 0
        warning(date(), " [WARNING] No tracer element specified for Molecule #", MoleculeNo, ".")
      }
      
      TransitionNo <- 0
      tmpTransitionsExpectedList <- list()
      tmpTransitionsExpectedListNames.vec <- vector()
      for (LabelIon in 0:MaxLabelProduct) {
        TransitionNo <- TransitionNo + 1
        tmpTrans1 <- LabelIon * ElementInfo[[Tracer]][[2]]
        
        tmpTransitionsExpectedList[[TransitionNo]] <- list(ProductIon = tmpTrans1, Precursor = tmpTrans1)
        tmpTransitionsExpectedListNames.vec <- c(tmpTransitionsExpectedListNames.vec, str_c(MoleculesName[MoleculeNo], "_", tmpTrans1))
        
      }  #LabelIon
      names(tmpTransitionsExpectedList) <- tmpTransitionsExpectedListNames.vec
      
      TransitionsExpected.df <- data.frame(ProductIon = as.numeric(unlist(lapply(tmpTransitionsExpectedList, function(x) x["ProductIon"]))), Precursor = as.numeric(unlist(lapply(tmpTransitionsExpectedList, 
                                                                                                                                                                                  function(x) x["Precursor"]))))
      
      
      rownames(TransitionsExpected.df) <- tmpTransitionsExpectedListNames.vec
      
      TransitionsExpected <- TransitionsExpected.df
      
    }  # NumberFragments==1
    MoleculeInfo[[MoleculeNo]][["TransitionsExpected"]] <- TransitionsExpected
  }  #  MoleculeNo
  
  if(verbose){message(date(), " :: calculating transitions [OK]\n")}
  return(MoleculeInfo)
  
}