# Calculates the mean isotopic enrichment for the corrected values.


CalcMeanEnrichment <- function(MoleculeInfo, MoleculesTotal, SamplesTotal, roundDigits, CorrectionResultList, UltraHighRes, verbose) {
    
    if(verbose){message("\n", date(), " :: calculating mean enrichment ...")}
    
    MeanEnrichmentResultList <- list()
    
    for (SampleNo in seq_len(SamplesTotal)) {
        
        MeanEnrichmentSample <- list()
        
        for (MoleculeNo in seq_len(MoleculesTotal)) {
            
            MoleculeData <- MoleculeInfo[[MoleculeNo]]
            MoleculeName <- names(MoleculeInfo[MoleculeNo])
            CorrectionData <- CorrectionResultList[[SampleNo]][[MoleculeNo]]
            Transitions <- MoleculeData[["Transitions"]]
            NumberTransitions <- nrow(Transitions)
            
            if (UltraHighRes == FALSE) {
                
                MeanEnrichmentTmp <- 0
                
                for (TransitionNo in seq_len(NumberTransitions)) {
                  
                  if (is.na(CorrectionData[["CorrectedFractions"]][TransitionNo]) == FALSE) {
                    
                    # It is ok to use the mass shift values and not neccessarily the label values here, because dividing by TotalMassShiftMax will in the end lead to correct
                    # values
                    
                    TotalMassShift <- MoleculeData$Transitions$Precursor[[TransitionNo]]
                    
                    MeanEnrichmentTmp <- MeanEnrichmentTmp + as.numeric(CorrectionData[["CorrectedFractions"]][TransitionNo] * TotalMassShift)
                    
                  }
                  
                }  #TransitionNo
                
                TotalMassShiftMax <- max(MoleculeData$TransitionsExpected$Precursor)
                
                MeanEnrichmentSample[[MoleculeName]] <- round(MeanEnrichmentTmp/TotalMassShiftMax, digits = roundDigits)
                
            } else if (UltraHighRes) {
                
                Tracers <- MoleculeData$Fragment_1$Tracer
                
                for (Tracer in names(Tracers)) {
                  
                  MoleculeNamePlusTracer <- paste0(MoleculeName, "_", Tracer)
                  Label <- MoleculeData$Transitions[[Tracer]]
                  MeanEnrichmentTmp <- 0
                  
                  for (TransitionNo in seq_len(NumberTransitions)) {
                    
                    if (is.na(CorrectionData[["CorrectedFractions"]][TransitionNo]) == FALSE) {
                      
                      MeanEnrichmentTmp <- MeanEnrichmentTmp + CorrectionData[["CorrectedFractions"]][TransitionNo] * Label[[TransitionNo]]
                      
                    }
                    
                  }
                  
                  MeanEnrichmentSample[[MoleculeNamePlusTracer]] <- round(MeanEnrichmentTmp/Tracers[[Tracer]], digits = roundDigits)
                  
                }
                
            }
            
        }
        
        MeanEnrichmentResultList[[SampleNo]] <- MeanEnrichmentSample
        
    }
    
    names(MeanEnrichmentResultList) <- names(CorrectionResultList)
    
    MeanEnrichmentDataOutput <-
      as.data.frame(convertEnrichmentList2matrix(input = MeanEnrichmentResultList, verbose=verbose))
    
    if(verbose){message(date(), " :: calculating mean enrichment [OK]\n")}
    return(MeanEnrichmentDataOutput)
    
}  #MeanEnrichment
