# Data correction using the probability matrix

# The function 'RawDataCorrection' performs the actual
# data correction. It uses the values extracted from the measurement file
# and corrects molecule by molecule for each sample. This is done by
# numerically solving a linear equation system. Here,
# the uncorrected value of a given transition is a linear combination
# of the corrected transition values with their
# probability of contributing to the uncorrected value
# as their coefficients (derived from the probability matrix,
# ProbMatrixComplete). The algorithm works with the
# constraint that a solution of the linear equation
# system cannot be < 0.

#' @importFrom pracma lsqlincon
RawDataCorrection <- function(UncorrectedData, MoleculeData, MoleculeName,
                              ProbMatrix, MoleculeNo, SampleNo, SampleName,
                              roundDigits, logEnvironment, verbose) {
  NumberTransitions <- nrow(MoleculeData[["Transitions"]])

  ProbMatrixComplete <- ProbMatrix

  # In the following, the linear equation system is solved with linear
  # inequality constraints: ConstraintMatrix*SolutionVector <= ConstraintVector.
  # The ConstraintVector is 0, the ConstraintMatrix is -1 at each
  # ConstraintMatrix(TransitionNo, TransitionNo) position.
  # It is thereby assured that the solutions cannot be < 0.
  # If a solution is < 0, ConstraintMatrix*SolutionVector becomes positive and
  # the constraint ConstraintMatrix*SolutionVector <=
  # ConstraintVector (which is 0) is not fulfilled.

  ConstraintVector <- rep(0, NumberTransitions)
  ConstraintMatrix <- matrix(0, nrow = NumberTransitions, ncol = NumberTransitions)

  for (TransitionNo in seq_len(NumberTransitions)) {
    ConstraintMatrix[TransitionNo, TransitionNo] <- -1
  }

  # Check the vector of uncorrected values for values that are not a number
  # (e.g. NA) and thus missing.
  CheckForNaN <- is.na(UncorrectedData)

  # If there are no missing values, perform the correction for the
  # current vector of uncorrected values UncorrectedData.
  # This is done by numerically solving a linear equation system where each
  # uncorrected value is a linear combination of unknown corrected values with
  # their probabilites of contributing to the uncorrected value as the constants
  # (from ProbMatrix).  ||ProbMatrix*CorrTransitionsPlus - UncorrectedData||^2 is minimized in a linear least squares
  # approach with the constraint that the corrected values (CorrTransitionsPlus)
  # must not be < 0. In addition to CorrTransitionsPlus, also the residuals of
  # the solving process are given (CorrResiduals).
  
  if (sum(CheckForNaN) == 0) {
    CorrTransitionsPlus <- pracma::lsqlincon(C = ProbMatrix, d = as.numeric(UncorrectedData), A = ConstraintMatrix, b = ConstraintVector)
    ManResiduals <- ProbMatrix %*% CorrTransitionsPlus - UncorrectedData
  } else {
    notification <- paste0(
      "In measurement data file: Measurement data of molecule ", MoleculeName, " in sample ", SampleName, " contains NA values. The correction performed may be less acurate.",
      "\nBe especially careful when considering fraction and mean enrichment values from samples with missing values."
    )
    errorHandler(notification, logEnvironment, "warning", verbose = verbose)

    CorrTransitionsPlus <- vector()
    ManResiduals <- vector()

    # Get indices of missing values. The rows with these indices removed from
    # UncorrectedData and the ConstraintVector. Row and column with these
    # indices are removed from ProbMatrix and ConstraintMatrix.
    # This way, the missing values are completely removed from the
    # linear equation system.

    # get index of NaN values
    nan.index <- as.numeric(which(CheckForNaN))

    # remove respective rows/columns
    UncorrectedData.nan <- UncorrectedData[-nan.index]
    ProbMatrix.nan <- ProbMatrix[-nan.index, -nan.index]
    ConstraintMatrix.nan <- ConstraintMatrix[-nan.index, -nan.index]
    ConstraintVector.nan <- ConstraintVector[-nan.index]

    # Here, the linear equation system is solved with a set of equations that is reduced because of the missing values. If the missing values are expected to
    # be high, the effect on the accuracy of the correction is more pronounced than if they are expected to be close to 0.

    if (length(nan.index) < NumberTransitions) {
      CorrTransitionsPlusNaN <- pracma::lsqlincon(C = ProbMatrix.nan, d = as.numeric(UncorrectedData.nan), A = ConstraintMatrix.nan, b = ConstraintVector.nan)
      ManResidualsNaN <- ProbMatrix.nan %*% CorrTransitionsPlusNaN - UncorrectedData.nan
    } else {
      notification <- paste0(
        "In measurement data file: All measurements are NA for molecule ", MoleculeName, " in sample ", SampleName ,
        "."
      )
      errorHandler(notification, logEnvironment, "warning", verbose)
    }

    # To be able to correctly assign corrected values to their measurement tag when missing values are present, a full size corrected value vector
    # CorrTransitionsPlus is made from the reduced size corrected values vector CorrTransitionsPlusNaN.  This is done by adding the missing values as NA at
    # the associated vector index.  In the same way a full size CorrResiduals vector is produced. A warning is given in the WarningsLog for each missing
    # value.

    i <- 1

    # NumberTransitions<-ncol(ProbMatrix)

    for (TransitionNo in seq_len(NumberTransitions)) {
      if (!CheckForNaN[TransitionNo]) {
        # 2017-08-28

        CorrTransitionsPlus[TransitionNo] <- CorrTransitionsPlusNaN[i]
        ManResiduals[TransitionNo] <- ManResidualsNaN[i]

        i <- i + 1
      } else {
        CorrTransitionsPlus[TransitionNo] <- NA
        ManResiduals[TransitionNo] <- NA
      } # if
    } # TransitionNo
  } # sum(CheckForNaN)==0

  CorrTransitionsPlus <- round(CorrTransitionsPlus, digits = roundDigits)
  ManResiduals <- round(ManResiduals, digits = roundDigits)

  # The corrected values in CorrTransitionsPlus are values that correspond to
  # the full integral of the isotopologue abundance distribution of a given
  # labeled species. To get the value (CorrTransitions) that corresponds to the
  # species in the distribution that contains no natural isotopes of higher
  # mass, the corresponding CorrTransitionsPlus value has to be multiplied
  # with the probability that the species contains no isotopes of higher mass
  # due to natural abundance

  CorrTransitions <- vector()
  CorrTransitionsFractions <- vector()
  CorrTransitionsPlusFractions <- vector()
  RelativeResiduals <- vector()

  # Calculation of residuals relative to corrected data

  RelativeResiduals <- round(ManResiduals / CorrTransitionsPlus, digits = roundDigits)

  for (TransitionNo in seq_len(NumberTransitions)) {
    if (is.na(CorrTransitionsPlus[TransitionNo]) == FALSE) {
      CorrTransitions[TransitionNo] <- CorrTransitionsPlus[TransitionNo] * ProbMatrixComplete[TransitionNo, TransitionNo]
    } else {
      CorrTransitions[TransitionNo] <- NA
    }
  } # TransitionNo

  round(CorrTransitions, digits = roundDigits)

  # To be able to calculate fractions, the vectors of the corrected values are summed.
  CorrTransitionsPlusSum <- sum(CorrTransitionsPlus[which(CheckForNaN == FALSE)])
  CorrTransitionsSum <- sum(CorrTransitions[which(CheckForNaN == FALSE)])

  # Calculation of fractions

  CorrTransitionsFractions <- vector()
  CorrTransitionsPlusFractions <- vector()

  for (TransitionNo in seq_len(NumberTransitions)) {
    CorrTransitionsFractions[TransitionNo] <- CorrTransitions[TransitionNo] / CorrTransitionsSum
    CorrTransitionsPlusFractions[TransitionNo] <- CorrTransitionsPlus[TransitionNo] / CorrTransitionsPlusSum
  } # TransitionNo

  round(CorrTransitionsPlusFractions, digits = roundDigits)
  round(CorrTransitionsFractions, digits = roundDigits)

  names(CorrTransitions) <- colnames(ProbMatrix)
  names(CorrTransitionsPlus) <- colnames(ProbMatrix)
  names(ManResiduals) <- colnames(ProbMatrix)
  names(RelativeResiduals) <- colnames(ProbMatrix)
  names(CorrTransitionsFractions) <- colnames(ProbMatrix)
  names(CorrTransitionsPlusFractions) <- colnames(ProbMatrix)

  returnList <- list(
    CorrectedMonoisotopic = CorrTransitions, Corrected = CorrTransitionsPlus, CorrectedMonoisotopicFractions = CorrTransitionsFractions,
    CorrectedFractions = CorrTransitionsPlusFractions, Residuals = ManResiduals, RelativeResiduals = RelativeResiduals
  )

  return(returnList)
}
