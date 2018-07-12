#' Graphical User Interface for IsoCorrectoR
#'
#' @import tcltk
#' @importFrom tcltk2 tk2tip
#' @importFrom tcltk2 tk2frame
#' @importFrom utils browseURL
#' 
#' @export
#' 
#' @examples
#' 
#'  IsoCorrectionGUI()
#' 
#' @return Calls internal function to display IsoCorrectoR's Graphical User Interface
IsoCorrectionGUI <- function() {
    
    initGUI()
    
}

initGUI <- function(advancedOptions = FALSE, filenameStartvalue = "", elementfileStartvalue = "", moleculefileStartvalue = "", DirOutStartvalue = "", FileOutStartvalue = "result", 
    UHRStartvalue = FALSE, CorrectTracerImpurityStartvalue = FALSE, CorrectTracerElementCoreStartvalue = TRUE, CalculateMeanEnrichmentStartvalue = TRUE, Calculation_thresholdStartvalue = 1e-08, 
    CalcTreshUHRStartvalue = 8, CalculateMonoisotopicProbStartvalue = FALSE, FileFormatStartvalue = "csv") {
    
    # Set this base environment variable as a flag to determine if execution should continue or be aborted in the case of directly starting from a batch file
    
    baseEnvRef <- baseenv()
    
    baseEnvRef$continueIsoCorrection <- TRUE
    
    # Specify styles and fonts
    .Tcl("ttk::style theme use default")
    fontTextLabel <- tcltk::tkfont.create(family = "font", size = 12)
    fontHelptext <- tcltk::tkfont.create(family = "font", size = 13)
    fontErrortext <- tcltk::tkfont.create(family = "font", size = 12)
    fontSmall <- tcltk::tkfont.create(family = "font", size = 4)
    tcltk::tcl("ttk::style", "configure", "But.TButton", "-background", "cornflower blue", "-borderwidth", "5", font = "helvetica 12", padding = 0, focusthickness = 5)
    tcltk::tcl("ttk::style", "configure", "Select.TButton", "-background", "cornflower blue", "-borderwidth", "2", font = "helvetica 10", padding = 0, focusthickness = 5)
    tcltk::tcl("ttk::style", "configure", "Check.TCheckbutton", relief = "raised", shiftrelief = 0, padding = 0, background = "cornflower blue", indicatorcolor = "white", 
        indicatorrelief = "sunken", indicatordiameter = 20, indicatormargin = 0, borderwidth = 4, focusthickness = 0, text = "test", justify = "left", foreground = "red")
    tcltk::tcl("ttk::style", "configure", "Test.TFrame", "-background", "white", "-borderwidth", "10", relief = "raised")
    tcltk::tcl("ttk::style", "configure", "Gray.TButton", "-background", "lightgray", "-borderwidth", "5", font = "helvetica 12", padding = 0, focusthickness = 5)
    tcltk::tcl("ttk::style", "configure", "Rad.TRadiobutton", "-background", "SystemMenu", fg = "cornflower blue", shiftrelief = 0, padding = 0, bordercolor = "cornflower blue", 
        indicatorcolor = "white", indicatorrelief = "sunken", indicatordiameter = 20, indicatormargin = 0, borderwidth = 3, focusthickness = 5, font = "helvetica 12")
    
    # Function of the 'Advanced Options' button
    ShowAdvancedOptions <- function() {
        tcltk::tkdestroy(win)
        initGUI(advancedOptions = TRUE, filenameStartvalue = tclvalue(filename.global), elementfileStartvalue = tclvalue(elementfile.global), moleculefileStartvalue = tclvalue(moleculefile.global), 
            DirOutStartvalue = tclvalue(DirOut.global), FileOutStartvalue = tclvalue(FileOut.global), UHRStartvalue = as.logical(as.integer(tclvalue(UltraHighResvalue))), 
            CorrectTracerImpurityStartvalue = as.logical(as.integer(tclvalue(CorrectTracerImpurityvalue))), CorrectTracerElementCoreStartvalue = as.logical(as.integer(tclvalue(CorrectTracerElementCorevalue))), 
            CalculateMeanEnrichmentStartvalue = as.logical(as.integer(tclvalue(CalculateMeanEnrichmentvalue))), Calculation_thresholdStartvalue = as.double(tclvalue(Calculation_thresholdvalue)), 
            CalcTreshUHRStartvalue = as.double(tclvalue(CalcTreshUHRvalue)), CalculateMonoisotopicProbStartvalue = as.logical(as.integer(tclvalue(CalcMonoisotopicProbabilitiesvalue))), 
            FileFormatStartvalue = tclvalue(FileFormatValue))
    }
    
    # Function of the 'Hide Advanced Options' button
    HideAdvancedOptions <- function() {
        tcltk::tkdestroy(win)
        initGUI(advancedOptions = FALSE, filenameStartvalue = tclvalue(filename.global), elementfileStartvalue = tclvalue(elementfile.global), moleculefileStartvalue = tclvalue(moleculefile.global), 
            DirOutStartvalue = tclvalue(DirOut.global), FileOutStartvalue = tclvalue(FileOut.global), UHRStartvalue = as.logical(as.integer(tclvalue(UltraHighResvalue))), 
            CorrectTracerImpurityStartvalue = as.logical(as.integer(tclvalue(CorrectTracerImpurityvalue))), CorrectTracerElementCoreStartvalue = as.logical(as.integer(tclvalue(CorrectTracerElementCorevalue))), 
            CalculateMeanEnrichmentStartvalue = as.logical(as.integer(tclvalue(CalculateMeanEnrichmentvalue))), Calculation_thresholdStartvalue = as.double(tclvalue(Calculation_thresholdvalue)), 
            CalcTreshUHRStartvalue = as.double(tclvalue(CalcTreshUHRvalue)), CalculateMonoisotopicProbStartvalue = as.logical(as.integer(tclvalue(CalcMonoisotopicProbabilitiesvalue))), 
            FileFormatStartvalue = tclvalue(FileFormatValue))
    }
    
    # Function for creating the error-pop-up-windows
    errorwindow <- function(text) {
        
        # Function of the 'OK' button in error-pop-up-windows
        destroy <- function() {
            tcltk::tkdestroy(errwin)
        }
        
        errwin <- tcltk::tktoplevel()
        tcltk::tkwm.title(errwin, "Error")
        errtext <- tcltk::tklabel(errwin, text = text, font = fontErrortext)
        tcltk::tkgrid(errtext)
        OK.button <- tcltk::ttkbutton(errwin, text = "OK", command = destroy)
        tcltk::tkgrid(OK.button)
    }
    
    # Function creating the window after hitting the Correct-Button
    finish <- function(results) {
        continuewindow <- tcltk::tktoplevel()
        tcltk::tkwm.title(continuewindow, "Continue")
        tcltk::tkfocus(continuewindow)
        
        contdestroy <- function() {
            tcltk::tkdestroy(continuewindow)
            baseEnvRef$continueIsoCorrection <- FALSE
        }
        
        # Rewrite window manager function to call contdestroy()
        
        tcltk::tcl("wm", "protocol", continuewindow, "WM_DELETE_WINDOW", contdestroy)
        
        # Function to start a completely new correction
        newcorrection <- function() {
            tcltk::tkdestroy(continuewindow)
            initGUI()
        }
        
        if (results$success == TRUE) {
            Headtext <- tcltk::tklabel(continuewindow, text = "Correction successful!", font = fontErrortext)
            tcltk::tkgrid(Headtext)
            newcorrection.button <- tcltk::ttkbutton(continuewindow, text = "Start new correction", command = newcorrection)
            tcltk::tkgrid(newcorrection.button)
        } else {
            Headtext <- tcltk::tklabel(continuewindow, text = paste0("Correction was aborted because an error has occured. Error:\n\n", results$error, "\n"), 
                font = fontErrortext)
            tcltk::tkgrid(Headtext)
            newcorrection.button <- tcltk::ttkbutton(continuewindow, text = "Start new correction", command = newcorrection)
            tcltk::tkgrid(newcorrection.button)
        }
        blankspace <- tcltk::tklabel(continuewindow, text = " ", font = fontSmall)
        tcltk::tkgrid(blankspace)
        OK.button <- tcltk::ttkbutton(continuewindow, text = "Close", command = contdestroy)
        tcltk::tkgrid(OK.button)
        
    }
    
    # Function of the 'Correct' button
    correction <- function() {
        
        # Changing from tclvalues to R classes
        filenameVal <- tcltk::tclvalue(filename.global)
        elementfileVal <- tcltk::tclvalue(elementfile.global)
        moleculefileVal <- tcltk::tclvalue(moleculefile.global)
        FileOutVal <- tcltk::tclvalue(FileOut.global)
        UltraHighResVal <- as.logical(as.integer(tcltk::tclvalue(UltraHighResvalue)))
        CorrectTracerImpurityVal <- as.logical(as.integer(tcltk::tclvalue(CorrectTracerImpurityvalue)))
        CorrectTracerElementCoreVal <- as.logical(as.integer(tcltk::tclvalue(CorrectTracerElementCorevalue)))
        CalculateMeanEnrichmentVal <- as.logical(as.integer(tcltk::tclvalue(CalculateMeanEnrichmentvalue)))
        Calculation_thresholdVal <- as.double(tcltk::tclvalue(Calculation_thresholdvalue))
        CalcTreshUHRVal <- as.integer(tcltk::tclvalue(CalcTreshUHRvalue))
        CalcMonoisotopicProbabilitiesVal <- as.logical(as.integer(tcltk::tclvalue(CalcMonoisotopicProbabilitiesvalue)))
        DirOutVal <- tcltk::tclvalue(DirOut.global)
        FileFormatVal <- tcltk::tclvalue(FileFormatValue)
        
        # Intercept input errors
        if (filenameVal != "") {
            if (elementfileVal != "") {
                if (moleculefileVal != "") {
                  if (!is.na(Calculation_thresholdVal)) {
                    if (DirOutVal != "") {
                      if (file.exists(filenameVal)) {
                        if (file.exists(elementfileVal)) {
                          if (file.exists(moleculefileVal)) {
                            if (dir.exists(DirOutVal)) {
                              if (FileOutVal != "") {
                                if (!is.na(CalcTreshUHRVal)) {
                                  print("Correction")
                                  tcltk::tkdestroy(win)
                                  results <- IsoCorrection(MeasurementFile = filenameVal, ElementFile = elementfileVal, FileOutFormat = FileFormatVal, MoleculeFile = moleculefileVal, 
                                    UltraHighRes = UltraHighResVal, CorrectTracerImpurity = CorrectTracerImpurityVal, CorrectTracerElementCore = CorrectTracerElementCoreVal, 
                                    CalculateMeanEnrichment = CalculateMeanEnrichmentVal, CorrectAlsoMonoisotopic = CalcMonoisotopicProbabilitiesVal, DirOut = DirOutVal, 
                                    CalculationThreshold = Calculation_thresholdVal, CalculationThreshold_UHR = CalcTreshUHRVal, FileOut = FileOutVal)
                                  finish(results)
                                } else {
                                  errorwindow("Please enter numeric value for Calculation threshold UHR")
                                }
                              } else {
                                errorwindow("Please enter name of the output file")
                              }
                            } else {
                              errorwindow("Output directory does not exist")
                            }
                          } else {
                            errorwindow("Molecule File does not exist")
                          }
                        } else {
                          errorwindow("Element File does not exist")
                        }
                      } else {
                        errorwindow("Measurement File does not exist")
                      }
                    } else {
                      errorwindow("Please choose output directory")
                    }
                  } else {
                    errorwindow("Please enter numeric value for Calculation threshold")
                  }
                } else {
                  errorwindow("Please choose molecule file")
                }
            } else {
                errorwindow("Please choose element file")
            }
        } else {
            errorwindow("Please choose measurement file")
        }
    }
    # Function of the different categories in the help menu
    help <- function(text) {
        # Function of the 'OK' button in helpmenu window
        destroy <- function() {
            tcltk::tkdestroy(helpwin)
        }
        helpwin <- tcltk::tktoplevel(background = "white")
        tcltk::tkwm.title(helpwin, "Help")
        helptext <- tcltk::tklabel(helpwin, text = text, justify = "left", font = fontHelptext, background = "white")
        tcltk::tkgrid(helptext)
        OK.button <- tcltk::ttkbutton(helpwin, text = "OK", command = destroy)
        tcltk::tkgrid(OK.button)
    }
    
    # Function to select file via filechooser/browser
    filebrowser.filename <- function() {
        filename.global <<- tcltk::tclVar(tk_choose.files())
        tcltk::tkconfigure(entry.filename, textvariable = filename.global)
    }
    
    # Function to select elementfile via filechooser/browser
    filebrowser.elementfile <- function() {
        elementfile.global <<- tcltk::tclVar(tk_choose.files())
        tcltk::tkconfigure(entry.elementfile, textvariable = elementfile.global)
    }
    
    # Function to select moleculefile via filechooser/browser
    filebrowser.moleculefile <- function() {
        moleculefile.global <<- tcltk::tclVar(tk_choose.files())
        tcltk::tkconfigure(entry.moleculefile, textvariable = moleculefile.global)
    }
    
    # Function to select Output Directory via directorychooser/browser
    dirbrowser.DirOut <- function() {
        DirOut.global <<- tcltk::tclVar(tk_choose.dir())
        tcltk::tkconfigure(entry.DirOut, textvariable = DirOut.global)
    }
    
    # Function to open general help html-file
    openhelp <- function() {
        utils::browseURL(system.file("doc/IsoCorrectoR.html", package = "IsoCorrectoR"))
    }
    
    # Function to show information on IsoCorrectoR
    
    showAbout <- function() {
        
        # Function of the 'OK' button in error-pop-up-windows
        destroy <- function() {
            tcltk::tkdestroy(aboutwin)
        }
        
        aboutwin <- tcltk::tktoplevel()
        tcltk::tkwm.title(aboutwin, "About")
        abouttext <- tcltk::tklabel(aboutwin, text = paste0("IsoCorrectoR version ", packageVersion("IsoCorrectoR"), "\n", "Institute of Functional Genomics , University of Regensburg\n", 
            "IsoCorrectoR is licensed under ", packageDescription("IsoCorrectoR", fields = "License"), ". It is free software and comes without any warranty.\n Please run 'citation(\"IsoCorrectoR\")' for correct citation when using IsoCorrectoR for your research."))
        tcltk::tkgrid(abouttext)
        OK.button <- tcltk::ttkbutton(aboutwin, text = "OK", command = destroy)
        tcltk::tkgrid(OK.button)
        
    }
    
    # Create window
    win <- tcltk::tktoplevel()
    
    windestroy <- function() {
        tcltk::tkdestroy(win)
        baseEnvRef$continueIsoCorrection <- FALSE
    }
    
    # Rewrite window manager function to call windestroy()
    
    tcltk::tcl("wm", "protocol", win, "WM_DELETE_WINDOW", windestroy)
    
    tcltk::tkwm.title(win, "IsoCorrectoR")
    tcltk::tkfocus(win)
    
    # Add help menu
    topMenu <- tcltk::tkmenu(win)
    tcltk::tkconfigure(win, menu = topMenu, bd = 0)
    helpMenu <- tcltk::tkmenu(topMenu, tearoff = FALSE, activebackground = "cornflower blue", bg = "white", bd = 0)
    tcltk::tkadd(helpMenu, "command", label = "How to use IsoCorrectoR", command = function() openhelp())
    tcltk::tkadd(helpMenu, "command", label = "About", command = function() showAbout())
    tcltk::tkadd(topMenu, "cascade", label = "Help", menu = helpMenu)
    
    # Enter measurement file
    filename.global <- tcltk::tclVar(filenameStartvalue)
    filename.tip <- tcltk::tklabel(win, text = "Measurement File:", font = fontTextLabel, justify = "left")
    entry.filename <- tcltk::tkentry(win, width = "60", textvariable = filename.global, bg = "white")
    browse.button.filename <- tcltk::ttkbutton(win, text = "Select", command = filebrowser.filename)
    tcltk::tkgrid(filename.tip, row = 0, column = 0, sticky = "w")
    tcltk::tkgrid(entry.filename, row = 0, column = 2)
    tcltk::tkgrid(browse.button.filename, row = 0, column = 3)
    tcltk2::tk2tip(filename.tip, "Please enter or choose file with measured data")
    
    # Enter Molecule File
    moleculefile.global <- tcltk::tclVar(moleculefileStartvalue)
    moleculefile.tip <- tcltk::tklabel(win, text = "Molecule File:", font = fontTextLabel, justify = "left")
    entry.moleculefile <- tcltk::tkentry(win, width = "60", textvariable = moleculefile.global, bg = "white")
    browse.button.moleculefile <- tcltk::ttkbutton(win, text = "Select", command = filebrowser.moleculefile)
    tcltk::tkgrid(moleculefile.tip, row = 1, column = 0, sticky = "w")
    tcltk::tkgrid(entry.moleculefile, row = 1, column = 2)
    tcltk::tkgrid(browse.button.moleculefile, row = 1, column = 3)
    tcltk2::tk2tip(moleculefile.tip, "Please enter or choose file with molecule data")
    
    # Enter Element File
    elementfile.global <- tcltk::tclVar(elementfileStartvalue)
    elementfile.tip <- tcltk::tklabel(win, text = "Element File:", font = fontTextLabel, justify = "left")
    entry.elementfile <- tcltk::tkentry(win, width = "60", textvariable = elementfile.global, bg = "white")
    browse.button.elementfile <- tcltk::ttkbutton(win, text = "Select", command = filebrowser.elementfile)
    tcltk::tkgrid(elementfile.tip, row = 2, column = 0, sticky = "w")
    tcltk::tkgrid(entry.elementfile, row = 2, column = 2)
    tcltk::tkgrid(browse.button.elementfile, row = 2, column = 3)
    tcltk2::tk2tip(elementfile.tip, "Please enter or choose file with element data")
    
    # Enter DirOut
    
    DirOut.global <- tcltk::tclVar(DirOutStartvalue)
    DirOut.tip <- tcltk::tklabel(win, text = "Output Directory:", font = fontTextLabel, justify = "left")
    entry.DirOut <- tcltk::tkentry(win, width = "60", textvariable = DirOut.global, bg = "white")
    browse.button.DirOut <- tcltk::ttkbutton(win, text = "Select", command = dirbrowser.DirOut)
    tcltk::tkgrid(DirOut.tip, row = 3, column = 0, sticky = "w")
    tcltk::tkgrid(entry.DirOut, row = 3, column = 2)
    tcltk::tkgrid(browse.button.DirOut, row = 3, column = 3)
    tcltk2::tk2tip(DirOut.tip, "Please enter or choose directory path to write results to")
    
    # Name Output file
    if (FileOutStartvalue == "") {
        FileOutStartvalue <- "result"
    }
    FileOut.global <- tcltk::tclVar(FileOutStartvalue)
    FileOut.tip <- tcltk::tklabel(win, text = "Name Outputfile:", font = fontTextLabel, justify = "left")
    entry.FileOut <- tcltk::tkentry(win, width = "60", textvariable = FileOut.global, bg = "white")
    tcltk::tkgrid(FileOut.tip, row = 4, column = 0, sticky = "w")
    tcltk::tkgrid(entry.FileOut, row = 4, column = 2)
    tcltk2::tk2tip(FileOut.tip, "Please enter name for your output file")
    
    # Choose output file format via radiobuttons
    if (FileFormatStartvalue == "") {
        FileFormatStartvalue <- "csv"
    }
    FileFormatValue <- tcltk::tclVar(FileFormatStartvalue)
    FileFormat.tip <- tcltk::tklabel(win, text = "Format Outputfile:", font = fontTextLabel, justify = "left")
    frame <- tcltk2::tk2frame(win, borderwidth = 0)
    tcltk::tkgrid(frame, row = 5, column = 2)
    win$env$rb1 <- tcltk::ttkradiobutton(frame, text = ".csv  ")
    win$env$rb2 <- tcltk::ttkradiobutton(frame, text = ".xls")
    tcltk::tkconfigure(win$env$rb1, variable = FileFormatValue, value = "csv")
    tcltk::tkconfigure(win$env$rb2, variable = FileFormatValue, value = "xls")
    tcltk::tkgrid(FileFormat.tip, row = 5, column = 0, sticky = "w")
    tcltk::tkgrid(win$env$rb1, win$env$rb2)
    
    UltraHighResvalue <- tcltk::tclVar(as.character(as.integer(UHRStartvalue)))
    CorrectTracerImpurityvalue <- tcltk::tclVar(as.character(as.integer(CorrectTracerImpurityStartvalue)))
    CorrectTracerElementCorevalue <- tcltk::tclVar(as.character(as.integer(CorrectTracerElementCoreStartvalue)))
    CalculateMeanEnrichmentvalue <- tcltk::tclVar(as.character(as.integer(CalculateMeanEnrichmentStartvalue)))
    Calculation_thresholdvalue <- tcltk::tclVar(as.character(Calculation_thresholdStartvalue))
    CalcTreshUHRvalue <- tcltk::tclVar(as.character(CalcTreshUHRStartvalue))
    CalcMonoisotopicProbabilitiesvalue <- tcltk::tclVar(as.character(as.integer(CalculateMonoisotopicProbStartvalue)))
    
    # Check Button CorrectTracerImpurity
    check.CorrectTracerImpurity <- tcltk::ttkcheckbutton(win)
    CorrectTracerImpurity.tip <- tcltk::tklabel(win, text = "Correct Tracer Impurity:", font = fontTextLabel, justify = "left")
    tcltk::tkgrid(CorrectTracerImpurity.tip, sticky = "w")
    tcltk::tkgrid(check.CorrectTracerImpurity, row = 6, column = 1)
    tcltk::tkconfigure(check.CorrectTracerImpurity, variable = CorrectTracerImpurityvalue)
    tcltk2::tk2tip(CorrectTracerImpurity.tip, "Correct for isotopic impurity of the tracer substrate?")
    
    # Check Button CorrectTracerElementCore
    check.CorrectTracerElementCore <- tcltk::ttkcheckbutton(win)
    CorrectTracerElementCore.tip <- tcltk::tklabel(win, text = "Corr. Tracer Element Core:", font = fontTextLabel, justify = "left")
    tcltk::tkgrid(CorrectTracerElementCore.tip, sticky = "w")
    tcltk::tkgrid(check.CorrectTracerElementCore, row = 7, column = 1)
    tcltk::tkconfigure(check.CorrectTracerElementCore, variable = CorrectTracerElementCorevalue)
    tcltk2::tk2tip(CorrectTracerElementCore.tip, "Take into account the natural isotope abundance of the tracer element atoms in the core molecule when correcting?")
    
    # Check Button CalculateMeanEnrichment
    check.CalculateMeanEnrichment <- tcltk::ttkcheckbutton(win)
    CalculateMeanEnrichment.tip <- tcltk::tklabel(win, text = "Calculate Mean Enrichment:   ", font = fontTextLabel, justify = "left")
    tcltk::tkgrid(CalculateMeanEnrichment.tip, sticky = "w")
    tcltk::tkgrid(check.CalculateMeanEnrichment, row = 8, column = 1)
    tcltk::tkconfigure(check.CalculateMeanEnrichment, variable = CalculateMeanEnrichmentvalue)
    tcltk2::tk2tip(CalculateMeanEnrichment.tip, "Calculate mean enrichment?")
    
    # Check Button UltraHighRes
    check.UltraHighRes <- tcltk::ttkcheckbutton(win)
    UltraHighRes.tip <- tcltk::tklabel(win, text = "High Resolution Mode:", font = fontTextLabel, justify = "left")
    tcltk::tkgrid(UltraHighRes.tip, sticky = "w")
    tcltk::tkgrid(check.UltraHighRes, row = 9, column = 1)
    tcltk::tkconfigure(check.UltraHighRes, variable = UltraHighResvalue)
    tcltk2::tk2tip(UltraHighRes.tip, "Perform correction of high resolution MS data (high resolution data allows handling of multiple different tracer elements in the same measurement)?")
    
    # Advanced options
    
    if (advancedOptions == FALSE) {
        AO.button <- tcltk::ttkbutton(win, text = "Advanced Options", command = ShowAdvancedOptions)
        tcltk::tkgrid(AO.button, row = 9, column = 2)
    } else {
        # Hide Advanced Options Button
        HAO.button <- tcltk::ttkbutton(win, text = "Hide Advanced Options", command = HideAdvancedOptions)
        tcltk::tkgrid(HAO.button, row = 9, column = 2)
    }
    
    if (advancedOptions == TRUE) {
        # Check Monoisotopic Probability
        check.CalcMonoisotopicProbabilities <- tcltk::ttkcheckbutton(win)
        CalcMonoisotopicProbabilities.tip <- tcltk::tklabel(win, text = "Monoisotopic Results:", font = fontTextLabel, justify = "left")
        tcltk::tkgrid(CalcMonoisotopicProbabilities.tip, sticky = "w")
        tcltk::tkgrid(check.CalcMonoisotopicProbabilities, row = 10, column = 1)
        tcltk::tkconfigure(check.CalcMonoisotopicProbabilities, variable = CalcMonoisotopicProbabilitiesvalue)
        tcltk2::tk2tip(CalcMonoisotopicProbabilities.tip, "Get monoisotopic correction results in addition to normal results, see Help for more information.")
        
        # Enter Calculation threshold
        Calculation_threshold.tip <- tcltk::tklabel(win, text = "Calculation Threshold:", font = fontTextLabel, justify = "left")
        entry.Calculation_threshold <- tcltk::tkentry(win, width = "20", textvariable = Calculation_thresholdvalue, bg = "white")
        tcltk::tkgrid(Calculation_threshold.tip, sticky = "w")
        tcltk::tkgrid(entry.Calculation_threshold, row = 11, column = 2)
        tcltk2::tk2tip(Calculation_threshold.tip, "Usually, this parameter should not be changed. See Help for further information")
        
        # Enter Calculation threshold regarding UHR
        CalcTreshUHR.tip <- tcltk::tklabel(win, text = "Calc. Thresh. UHR:", font = fontTextLabel, justify = "left")
        entry.CalcTreshUHR <- tcltk::tkentry(win, width = "20", textvariable = CalcTreshUHRvalue, bg = "white")
        tcltk::tkgrid(CalcTreshUHR.tip, sticky = "w")
        tcltk::tkgrid(entry.CalcTreshUHR, row = 12, column = 2)
        tcltk2::tk2tip(CalcTreshUHR.tip, "Usually, this parameter should not be changed. See Help for further information")
        
    }
    
    blankspace <- tcltk::tklabel(win, text = " ", font = fontSmall)
    tcltk::tkgrid(blankspace)
    Correct.button <- tcltk::ttkbutton(win, text = "Start Correction", command = correction)
    tcltk::tkgrid(Correct.button, column = 2)
    
    tcltk::tkfocus(win)
    
}
