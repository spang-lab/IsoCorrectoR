# COMPUTE ISOTOPE COMBINATIONS OF TRACER ELEMENT

IsoCombinationsTracer <- function(MoleculeFragmentData, Fragment, ElementInfo, Transitions, NumberTracers, CorrectTracerImpurity, CorrectTracerElementCore, CalculationThreshold) {
  ProbTracerList <- list()
  MassShiftTracerList <- list()
  NumberIsoCombTracerEff <- list()
  IsotopesTracer <- list()

  NumberTransitions <- nrow(Transitions)
  NegIsoTracer <- MoleculeFragmentData[["NegIsoTracer"]]
  nTracerMax <- MoleculeFragmentData[["nTracerMax"]]
  MaxLabel <- MoleculeFragmentData[["MaxLabel"]]

  LabelPresent <- rep(0, nrow(Transitions))

  if (NumberTracers > 0) {
    IDTracer <- MoleculeFragmentData[["IDTracer"]]
    Element <- IDTracer
    LabelPresent <- vector()
    nTracer <- vector()

    for (TransitionNo in seq_len(NumberTransitions)) {

      # x is passed to the daughter function 'IsoCombinations'

      x <- TransitionNo

      # ElementInfo[[Element]][[2]] is the tracer isotope mass shift.
      # Dividing by this transforms the mass shift of a transition into
      # the number of label present.

      LabelPresent[x] <- Transitions[x, Fragment] / ElementInfo[[Element]][[2]]

      nTracer[x] <- as.numeric(nTracerMax - LabelPresent[x])

      # Do not constrain PlacesToAssign if tracer purity correction is active,
      # as tracer purity produces negative mass shifts which leads to similar
      # effects as NegIsoTracer.

      if (CorrectTracerImpurity || NegIsoTracer == 1) {
        PlacesToAssign <- as.numeric(ElementInfo[[Element]][[2]]) * MaxLabel
      } else {
        PlacesToAssign <- as.numeric(ElementInfo[[Element]][[2]] * (MaxLabel - LabelPresent[x]))
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

      IsoCombinationsResultList <- IsoCombinations(Fragment = Fragment, x = x, Element = Element, ElementInfo = ElementInfo, AvailablePlacesTotal = AvailablePlacesTotal, IsoCluster = IsoCluster, CalculationThreshold = CalculationThreshold)
      NumberIsoComb <- IsoCombinationsResultList[["NumberIsoComb"]]
      tmpProbTracer_Array <- vector()
      tmpMassShiftTracer_Array <- vector()
      tmpIsotopesTracer <- list()
      for (i in seq_len(NumberIsoComb)) {
        tmpProbTracer_Array[i] <- IsoCombinationsResultList[["ProbArray"]][[i]]
        tmpMassShiftTracer_Array[i] <- IsoCombinationsResultList[["MassShiftArray"]][[i]]
      } # i

      ProbTracerList[[x]] <- tmpProbTracer_Array
      MassShiftTracerList[[x]] <- tmpMassShiftTracer_Array
      NumberIsoCombTracerEff[[x]] <- NumberIsoComb
      IsotopesTracer[[x]] <- IsoCombinationsResultList[["Isotopes"]]
    } # TransitionsNo (represented by x)

    names(ProbTracerList) <- rownames(Transitions)
    names(MassShiftTracerList) <- rownames(Transitions)
    names(NumberIsoCombTracerEff) <- rownames(Transitions)
    names(IsotopesTracer) <- rownames(Transitions)
  } else {
    LabelPresent <- rep(0, NumberTransitions)
  }

  return(list(
    "ProbTracerList" = ProbTracerList,
    "MassShiftTracerList" = MassShiftTracerList,
    "NumberIsoCombTracerEff" = NumberIsoCombTracerEff,
    "IsotopesTracer" = IsotopesTracer,
    "LabelPresent" = LabelPresent
  ))
}
