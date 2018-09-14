# COMPUTE ISOTOPE COMBINATIONS OF NON-TRACER ELEMENT

# Combinations, probabilites and mass shifts of the non-tracer elements are
# computed, given that not only the tracer element is provided in
# the molcule file. As with the tracer element,
# IsoCombinationsMaster() feeds into IsoCombinations() also for the other elements.

IsoCombinationsNonTracer <- function(MoleculeFragmentData, Fragment, ElementInfo,
                                     NumberTracers, IDTracer, CalculationThreshold) {
  ProbElemList <- list()
  MassShiftElemList <- list()
  NumberIsoCombEff <- list()
  IsotopesElem <- list()

  ElementNonTracer_label <- "NonTracer"
  ElementsNonTracer <- MoleculeFragmentData[[ElementNonTracer_label]]

  if (length(ElementsNonTracer) == 0) {
    ElementsNonTracer <- MoleculeFragmentData[["ZeroTracer"]]
    ElementNonTracer_label <- "ZeroTracer"
  }

  NumberElementsNonTracer <- length(ElementsNonTracer)

  if (NumberElementsNonTracer > 0) {
    for (NonTracer in seq_len(NumberElementsNonTracer)) {

      # Definition of parameters for the non-tracer elements
      x <- NonTracer

      Element <- names(ElementsNonTracer)[x]

      if (NumberTracers > 0) {
        PlacesToAssign <- as.numeric(ElementInfo[[IDTracer]][[2]] * MoleculeFragmentData[["MaxLabel"]])
      } else {
        PlacesToAssign <- 0
      } # NumberTracers > 0

      AvailablePlacesTotal <- as.numeric(MoleculeFragmentData[[ElementNonTracer_label]][Element])

      if (PlacesToAssign <= AvailablePlacesTotal) {
        IsoCluster <- PlacesToAssign
      } else {
        IsoCluster <- AvailablePlacesTotal
      }

      # Calling of IsoCombinations()
      IsoCombinationsResultList <- IsoCombinations(
        Fragment = Fragment, x = x, Element = Element, ElementInfo = ElementInfo,
        AvailablePlacesTotal = AvailablePlacesTotal, IsoCluster = IsoCluster, CalculationThreshold = CalculationThreshold
      )
      NumberIsoComb <- IsoCombinationsResultList[["NumberIsoComb"]]
      tmpProbElem_Array <- vector()
      tmpMassShiftElem_Array <- vector()

      # tmpIsotopesElem <- list()

      for (i in seq_len(NumberIsoComb)) {
        tmpProbElem_Array[i] <- IsoCombinationsResultList[["ProbArray"]][i]
        tmpMassShiftElem_Array[i] <- IsoCombinationsResultList[["MassShiftArray"]][i]
      } # i

      ProbElemList[[x]] <- tmpProbElem_Array
      MassShiftElemList[[x]] <- tmpMassShiftElem_Array
      NumberIsoCombEff[[x]] <- NumberIsoComb
      IsotopesElem[[x]] <- IsoCombinationsResultList[["Isotopes"]]
    } # NonTracer

    names(ProbElemList) <- names(ElementsNonTracer)
    names(MassShiftElemList) <- names(ElementsNonTracer)
    names(NumberIsoCombEff) <- names(ElementsNonTracer)
    names(IsotopesElem) <- names(ElementsNonTracer)
  } # NumberElementsNonTracer > 0

  return(list(
    "ProbElemList" = ProbElemList,
    "MassShiftElemList" = MassShiftElemList,
    "NumberIsoCombEff" = NumberIsoCombEff,
    "IsotopesElem" = IsotopesElem
  ))
}
