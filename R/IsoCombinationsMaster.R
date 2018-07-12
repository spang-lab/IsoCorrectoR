# Computation of isotope combinations, associated probabilities and mass shifts 
# for an element

#' @importFrom stringr str_c

IsoCombinationsMaster <- function(MoleculeArray, MoleculeNo, Fragment, ElementArray, Prob_threshold_isotopes, CorrectTracerImpurity, CorrectTracerElementCore) {
  
  MoleculeData <- MoleculeArray[[MoleculeNo]][[stringr::str_c("Fragment_", Fragment)]]  # get individual fragment data
  
  ElementNonTracer_label <- "NonTracer"
  
  ElementsNonTracer <- MoleculeData[[ElementNonTracer_label]]
  
  if (length(ElementsNonTracer) == 0) {
    ElementsNonTracer <- MoleculeData[["ZeroTracer"]]
    ElementNonTracer_label <- "ZeroTracer"
  }
  
  NumberElementsNonTracer <- length(ElementsNonTracer)
  NumberTracers <- length(MoleculeData[["Tracer"]])
  
  Transitions <- MoleculeArray[[MoleculeNo]][["Transitions"]]
  
  NumberTransitions <- nrow(Transitions)
  
  NegIsoTracer <- MoleculeData[["NegIsoTracer"]]
  nTracerMax <- MoleculeData[["nTracerMax"]]
  MaxLabel <- MoleculeData[["MaxLabel"]]
  
  ProbTracer_Array <- list()
  MassShiftTracer_Array <- list()
  NumberIsoCombTracerEff <- list()
  IsotopesTracer <- list()
  ProbElem_Array <- list()
  MassShiftElem_Array <- list()
  NumberIsoCombEff <- list()
  IsotopesElem <- list()
  ProbTracerImpurity_Array <- list()
  MassShiftTracerImpurity_Array <- list()
  NumberImpurityCombEff <- list()
  
  # a) TRACER ISOTOPE COMBINATIONS
  
  LabelPresent <- rep(0, nrow(Transitions))
  
  if (NumberTracers > 0) {
    
    IDTracer <- MoleculeData[["IDTracer"]]
    Element <- IDTracer
    LabelPresent <- vector()
    nTracer <- vector()
    
    for (TransitionNo in seq_len(NumberTransitions)) {
      
      # x is passed to the daughter function 'IsoCombinations'
      
      x <- TransitionNo
      
      # ElementArray[[Element]][[2]] is the tracer isotope mass shift.  
      # Dividing by this transforms the mass shift of a transition into 
      # the number of label present.
      
      LabelPresent[x] <- Transitions[x, Fragment]/ElementArray[[Element]][[2]]
      
      nTracer[x] <- as.numeric(nTracerMax - LabelPresent[x])
      
      # Do not constrain PlacesToAssign if tracer purity correction is active, 
      # as tracer purity produces negative mass shifts which leads to similar 
      # effects as NegIsoTracer.
      
      if (CorrectTracerImpurity || NegIsoTracer == 1) {
        PlacesToAssign <- as.numeric(ElementArray[[Element]][[2]]) * MaxLabel
      } else {
        PlacesToAssign <- as.numeric(ElementArray[[Element]][[2]] * (MaxLabel - LabelPresent[x]))
      }
      
      # Correct molecule core or not
      
      if (CorrectTracerElementCore) {
        AvailablePlacesTotal <- nTracer[x]
      } else if (!CorrectTracerElementCore) {
        AvailablePlacesTotal <- as.numeric(nTracerMax) - as.numeric(MaxLabel)
      } else {
        stop(date(), " :: Invalid value for `CorrectTracerElementCore`: ", CorrectTracerElementCore)
      }
      
      # Definition of the IsoCluster variable that is transferred to the 
      # daughter function IsoCombinations() that computes the actual 
      # isotope combinations
      
      if (PlacesToAssign <= AvailablePlacesTotal) {
        IsoCluster <- PlacesToAssign
      } else {
        IsoCluster <- AvailablePlacesTotal
      }
      
      IsoCombinationsResultList <- IsoCombinations(MoleculeNo = MoleculeNo, Fragment = Fragment, x = x, Element = Element, ElementArray = ElementArray, 
                                                   AvailablePlacesTotal = AvailablePlacesTotal, IsoCluster = IsoCluster, Prob_threshold_isotopes = Prob_threshold_isotopes)
      NumberIsoComb <- IsoCombinationsResultList[["NumberIsoComb"]]
      tmpProbTracer_Array <- vector()
      tmpMassShiftTracer_Array <- vector()
      tmpIsotopesTracer <- list()
      for (i in seq_len(NumberIsoComb)) {
        
        tmpProbTracer_Array[i] <- IsoCombinationsResultList[["ProbArray"]][[i]]
        tmpMassShiftTracer_Array[i] <- IsoCombinationsResultList[["MassShiftArray"]][[i]]
        
      }  #i
      
      ProbTracer_Array[[x]] <- tmpProbTracer_Array
      MassShiftTracer_Array[[x]] <- tmpMassShiftTracer_Array
      NumberIsoCombTracerEff[[x]] <- NumberIsoComb
      IsotopesTracer[[x]] <- IsoCombinationsResultList[["Isotopes"]]
      
    }  #TransitionsNo (represented by x)
    
    names(ProbTracer_Array) <- rownames(Transitions)
    names(MassShiftTracer_Array) <- rownames(Transitions)
    names(NumberIsoCombTracerEff) <- rownames(Transitions)
    names(IsotopesTracer) <- rownames(Transitions)
  } else {
    
    LabelPresent <- rep(0, NumberTransitions)
    
  }
  
  # b)NON-TRACER ELEMENTS Analogous to the tracer element, the isotope 
  # combinations, probabilites and mass shifts of the other elements are 
  # computed, given that not only the tracer element is provided in 
  # the molcule file. As with the tracer element, 
  # IsoCombinationsMaster() feeds into IsoCombinations() also for the other elements.
  
  if (NumberElementsNonTracer > 0) 
  {
    
    for (NonTracer in seq_len(NumberElementsNonTracer)) {
      
      # Definition of parameters for the non-tracer elements
      x <- NonTracer
      
      Element <- names(ElementsNonTracer)[x]
      
      if (NumberTracers > 0) {
        
        IDTracer <- MoleculeData[["IDTracer"]]
        
        PlacesToAssign <- as.numeric(ElementArray[[IDTracer]][[2]] * MoleculeData[["MaxLabel"]])
        
      } else {
        PlacesToAssign <- 0
      }  #NumberTracers > 0
      
      AvailablePlacesTotal <- as.numeric(MoleculeData[[ElementNonTracer_label]][Element])
      
      if (PlacesToAssign <= AvailablePlacesTotal) {
        IsoCluster <- PlacesToAssign
      } else {
        IsoCluster <- AvailablePlacesTotal
      }
      
      # Calling of IsoCombinations()
      IsoCombinationsResultList <- IsoCombinations(MoleculeNo = MoleculeNo, Fragment = Fragment, x = x, Element = Element, ElementArray = ElementArray, 
                                                   AvailablePlacesTotal = AvailablePlacesTotal, IsoCluster = IsoCluster, Prob_threshold_isotopes = Prob_threshold_isotopes)
      
      NumberIsoComb <- IsoCombinationsResultList[["NumberIsoComb"]]
      tmpProbElem_Array <- vector()
      tmpMassShiftElem_Array <- vector()
      
      tmpIsotopesElem <- list()
      
      for (i in seq_len(NumberIsoComb)) {
        tmpProbElem_Array[i] <- IsoCombinationsResultList[["ProbArray"]][i]
        tmpMassShiftElem_Array[i] <- IsoCombinationsResultList[["MassShiftArray"]][i]
      }  #i
      
      ProbElem_Array[[x]] <- tmpProbElem_Array
      MassShiftElem_Array[[x]] <- tmpMassShiftElem_Array
      NumberIsoCombEff[[x]] <- NumberIsoComb
      IsotopesElem[[x]] <- IsoCombinationsResultList[["Isotopes"]]
      
    }  #NonTracer
    
    names(ProbElem_Array) <- names(ElementsNonTracer)
    names(MassShiftElem_Array) <- names(ElementsNonTracer)
    names(NumberIsoCombEff) <- names(ElementsNonTracer)
    names(IsotopesElem) <- names(ElementsNonTracer)
    
  }  # NumberElementsNonTracer > 0
  
  # c) PROBABILITIES AND MASS SHIFTS ASSOCIATED WITH TRACER IMPURITY 
  # Tracer purity is the probability that a tracer element atom in the 
  # tracer substrate that should be labelled actually is labelled. 
  # E.g.  a 99.0% isotopic purity 1,2-13C-Glucose has a 99.0% chance for 
  # each of its carbon atom positions 1 and 2 of containing a 13C. 
  # Thus, there is a 1.0% chance for each of those carbons that they are not 13C.
  # Consequently, molecules that contain the tracer
  # isotope due to metabolic transfer from the tracer substrate and not due to 
  # natural abundance have a certain chance of contributing to labelling states
  # with less tracer incorporated. This is due to the decrease in mass shift 
  # associated with tracer impurity (e.g. 12C instead of 13C at a carbon position).
  # In the following, for each labelling state the probability of having a 
  # certain number of 'impure' tracer positions as well as the associated 
  # (negative) mass shift are is calculated. 
  # The maximum amount of impure tracer positions possible depends on the number
  # of label present in a given transition.
  
  if (CorrectTracerImpurity) {
    if (NumberTracers > 0) 
    {
      for (TransitionNo in seq_len(NumberTransitions)) {
        
        tmpProbTracerImpurity_Array <- vector()
        tmpMassShiftTracerImpurity_Array <- vector()
        tmpNumberImpurityCombEff <- vector()
        
        for (ImpureTracer in 0:LabelPresent[TransitionNo]) {
          
          ProbImpurity <- choose(LabelPresent[TransitionNo], ImpureTracer) * (ElementArray[[IDTracer]][["Tracer.purity"]]^(LabelPresent[TransitionNo] - 
                                                                                                                             ImpureTracer)) * ((1 - ElementArray[[IDTracer]][["Tracer.purity"]])^ImpureTracer)
          
          if (ProbImpurity > Prob_threshold_isotopes) {
            tmpProbTracerImpurity_Array[ImpureTracer + 1] <- ProbImpurity
            tmpMassShiftTracerImpurity_Array[ImpureTracer + 1] <- (-1) * ImpureTracer * ElementArray[[IDTracer]][["Tracer.isotope.mass.shift"]]
            NumberImpurityCombEff[TransitionNo] <- ImpureTracer
          } else {
            ImpureTracer <- LabelPresent + 1
          }  #ProbImpurity > Prob_threshold_isotopes
        }  #ImpureTracer
        
        ProbTracerImpurity_Array[[TransitionNo]] <- tmpProbTracerImpurity_Array
        MassShiftTracerImpurity_Array[[TransitionNo]] <- tmpMassShiftTracerImpurity_Array
      }  #TransitionNo
      
      names(ProbTracerImpurity_Array) <- rownames(Transitions)
      names(MassShiftTracerImpurity_Array) <- rownames(Transitions)
      names(NumberImpurityCombEff) <- rownames(Transitions)
    }  #NumberTracers > 0
  } else {
    # message(date(),' :: no CorrectTracerImpurity performed ...')
    ProbTracerImpurity_Array <- "No Correction for Tracer Impurity performed"
    MassShiftTracerImpurity_Array <- "No Correction for Tracer Impurity performed"
    NumberImpurityCombEff <- "No Correction for Tracer Impurity performed"
    
  }  #CorrectTracerImpurity=TRUE
  
  returnList <- list(ProbTracer_Array = ProbTracer_Array, ProbElem_Array = ProbElem_Array, ProbTracerImpurity_Array = ProbTracerImpurity_Array, MassShiftTracer_Array = MassShiftTracer_Array, 
                     MassShiftElem_Array = MassShiftElem_Array, MassShiftTracerImpurity_Array = MassShiftTracerImpurity_Array, NumberIsoCombEff = unlist(NumberIsoCombEff), 
                     NumberImpurityCombEff = NumberImpurityCombEff, NumberIsoCombTracerEff = unlist(NumberIsoCombTracerEff), IsotopesTracer = IsotopesTracer, IsotopesElem = IsotopesElem, 
                     LabelPresent = LabelPresent)
  
  return(returnList)
  
}  #IsoCombinationsMaster()
