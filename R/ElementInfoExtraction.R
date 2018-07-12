# Element Parameter Extraction

#' @importFrom stringr str_split
#' @importFrom tibble as_tibble
#' @importFrom magrittr '%>%'
#' 
ElementInfoExtraction <- function(ElementData, UltraHighRes, logEnvironment) {
  message(date(), " :: reading ElementFile")
  
  # Using two nested loops, the isotope information of each element in column 2 is extracted by splitting it using the delimiters '/'
  # (isotopes1/isotope2...) and '_' (abundance-mass shift of given isotope).
  
  ElementList <- list()
  
  for (Element in seq_len(nrow(ElementData))) {
    
    Isotopes <- stringr::str_split(ElementData[Element, 2], pattern = "/", simplify = TRUE)
    
    ElementArray <- matrix(NA, nrow = length(Isotopes), ncol = 2)
    colnames(ElementArray) <- c("IsotopeAbundance", "MassShift")
    for (Isotope in seq_len(length(Isotopes))) {
      
      IsoAbundanceMassShift <- stringr::str_split(Isotopes[Isotope], pattern = "_", simplify = TRUE)
      
      for (i in seq_len(ncol(ElementArray))) {
        ElementArray[Isotope, i] <- as.numeric(IsoAbundanceMassShift[i])
      }
      
    }  #Isotope
    
    # Sort Isotopes of each element in ElementArray ascending by probability. This is required to be able to later calculate probabilities in descending
    # order, making it possible to stop the calculation at a defined threshold without losing higher probability values.
    
    ElementArray_tbl <- tibble::as_tibble(ElementArray)
    
    # make R CMD check happy
    ElementArray_tbl <- tibble::as_tibble(ElementArray[order(ElementArray[, "IsotopeAbundance"], decreasing = TRUE), ])
    
    # makes R CMD check unhappy ElementArray_tbl<-ElementArray_tbl %>% arrange(desc(IsotopeAbundance))
    
    
    tmp.list <- list()
    tmp2.list <- list()
    tmp.list[[1]] <- ElementArray_tbl
    tmp2.list[[1]] <- tmp.list
    
    # now add rest of ElementData
    for (i in 3:ncol(ElementData)) {
      tmp <- ElementData[Element, i]
      
      # Making NAs of ElementArray 0
      
      if (is.na(tmp)) {
        tmp <- 0
      }
      tmp2.list[[i - 1]] <- tmp
    }  #i
    
    names(tmp2.list) <- c("Isotopes", colnames(ElementData)[3:ncol(ElementData)])
    ElementList[[Element]] <- tmp2.list
    rm(tmp.list)
    rm(tmp2.list)
    
  }  #Element
  
  names(ElementList) <- ElementData[[1]]
  
  checkElementDataLogic(ElementList, UltraHighRes, logEnvironment)
  
  message(date(), " :: Leaving ElementInfoExtraction()")
  return(ElementList)
}  # end of function
