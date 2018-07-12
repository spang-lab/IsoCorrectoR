## ----setup, include=FALSE--------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readxl)

## ---- echo=FALSE, eval=TRUE, message=FALSE---------------------------------
# not nice, but we need the IsoCorrectoR data to generate to tool_features table !

# load package
library(IsoCorrectoR)

# load IsoCorrectoR example data
data(IsoCorrectoR)

## ----toolFeaturesTable, echo=FALSE, eval=TRUE------------------------------
knitr::kable(
  IsoCorrectoR$tool_features
)

## ---- echo=TRUE, eval=FALSE------------------------------------------------
#  IsoCorrection(MeasurementFile=NA, ElementFile=NA, MoleculeFile=NA,
#                CorrectTracerImpurity=FALSE, CorrectTracerElementCore=TRUE,
#                CalculateMeanEnrichment=TRUE, UltraHighRes=FALSE,
#                DirOut='.', FileOut='result', FileOutFormat='csv',
#                ReturnResultsObject=FALSE, CorrectAlsoMonoisotopic=FALSE,
#                CalculationThreshold=10^-8, CalculationThreshold_UHR=8,
#                Testmode=FALSE)

## ----moleculeFileExampleNormalRes, echo=FALSE, eval=TRUE-------------------

fileExample <- IsoCorrectoR[["normal_resolution"]][["molecule_file"]]

fileExample[4:7, 3] <- ""

knitr::kable(
  fileExample, align = "l", caption="Molecule information for normal resolution data"
)

## ----moleculeFileExampleHighRes, echo=FALSE, eval=TRUE---------------------

fileExample <- IsoCorrectoR[["high_resolution"]][["molecule_file"]]

fileExample[is.na(fileExample)] <- ""

knitr::kable(
  fileExample, align = "l", caption="Molecule information for high resolution data"
)

## ----measurementFileExampleNormRes, echo=FALSE, eval=TRUE------------------
fileExample <- IsoCorrectoR[["normal_resolution"]][["measurement_file"]]

#Subset and adjust example to illustrate explanations

fileExample <- fileExample[c(1:6, 40:43),1:6]

fileExample[7:10,1] <- gsub('.{2}$', '', fileExample[7:10,1])

fileExample[c(4,6), "Sample3"] <- ""
fileExample[4, "Sample5"] <- ""
fileExample[9, "Sample1"] <- ""

knitr::kable(
  fileExample, align = "l", row.names = FALSE, caption="Measurement information for normal resolution data"
)

## ----measurementFileExampleHighRes, echo=FALSE, eval=TRUE------------------
fileExample <- IsoCorrectoR[["high_resolution"]][["measurement_file"]]

#Subset and adjust example to illustrate explanations

fileExample <- fileExample[1:21,1:6]

fileExample[20,"Sample3"] <- ""
fileExample[21,"Sample3"] <- "0"
fileExample[16,"Sample3"] <- ""

knitr::kable(
  fileExample, align = "l", row.names = FALSE, caption="Measurement information for high resolution data"
)

## ----elementFileExample, echo=FALSE, eval=TRUE-----------------------------
fileExample <- IsoCorrectoR$element_file

fileExample[is.na(fileExample)] <- ""

knitr::kable(
  fileExample, align = "l", caption="Element information (resolution independent)"
)

## ----sessionInfo, echo=FALSE, eval=TRUE, messages=FALSE--------------------
sessionInfo()

