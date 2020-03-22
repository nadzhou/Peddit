#
# $Id: fftk_GeomOpt.tcl,v 1.10 2019/08/27 22:31:22 johns Exp $
#

#======================================================
namespace eval ::ForceFieldToolKit::GeomOpt:: {

    variable pdb
    variable com
    
    variable qmProc
    variable qmMem
    variable qmCharge
    variable qmMult
    variable qmRoute

    variable logFile
    variable optPdb

    variable qmSoft $::ForceFieldToolKit::qmSoft
}
#======================================================
proc ::ForceFieldToolKit::GeomOpt::init {} {
    
    # localize variables
    variable pdb
    variable com
    
    variable qmProc
    variable qmMem
    variable qmCharge
    variable qmMult
    variable qmRoute

    variable logFile
    variable optPdb

    variable qmSoft $::ForceFieldToolKit::qmSoft
    
    # Set Variables to Initial value
    set pdb {}
    set com {}
    
    ::ForceFieldToolKit::${qmSoft}::resetDefaultsGeomOpt

    set logFile {}
    set optPdb {}
    
}
#======================================================
proc ::ForceFieldToolKit::GeomOpt::sanityCheck {} {
    # checks to see that appropriate information is set prior to running
    
    # returns 1 if all input is sane
    # returns 0 if there is a problem
    
    # localize relevant GeomOpt variables
    
    variable pdb
    variable com
    
    variable qmProc
    variable qmMem
    variable qmCharge
    variable qmMult
    variable qmRoute
    
    # local variables
    set errorList {}
    set errorText ""
    
    # checks
    # make sure that pdb is entered and exists
    if { $pdb eq "" } {
        lappend errorList "No PDB file was specified."
    } else {
        if { ![file exists $pdb] } { lappend errorList "Cannot find PDB file." }
    }
    
    # make sure that com is enetered and exists
    if { $com eq "" } {
        lappend errorList "No output path was specified."
    } else {
        if { ![file writable [file dirname $com]] } { lappend errorList "Cannot write to output path." }
    }
     
    # validate gaussian settings (not particularly vigorous validation)
    # qmProc (processors)
    if { $qmProc eq "" } { lappend errorList "No processors were specified." }
    if { $qmProc <= 0 || $qmProc != [expr int($qmProc)] } { lappend errorList "Number of processors must be a positive integer." }
    # qmMem (memory)
    if { $qmMem eq "" } { lappend errorList "No memory was specified." }
    if { $qmMem <= 0 || $qmMem != [expr int($qmMem)]} { lappend errorList "Memory must be a postive integer." }
    # qmCharge (charge)
    if { $qmCharge eq "" } { lappend errorList "No charge was specified." }
    if { $qmCharge != [expr int($qmCharge)] } { lappend errorList "Charge must be an integer." }
    # qmMult (multiplicity)
    if { $qmMult eq "" } { lappend errorList "No multiplicity was specified." }
    if { $qmMult < 0 || $qmMult != [expr int($qmMult)] } { lappend errorList "Multiplicity must be a positive integer." }
    # qmRoute (route card for gaussian; just make sure it isn't empty)
    if { $qmRoute eq "" } { lappend errorList "Route card is empty." }
    
    # if there is an error, tell the user about it
    # return -1 to tell the calling proc that there is a problem
    if { [llength $errorList] > 0 } {
        foreach ele $errorList {
            set errorText [concat $errorText\n$ele]
        }
        tk_messageBox \
            -type ok \
            -icon warning \
            -message "Application halting due to the following errors:" \
            -detail $errorText
        
        # there are errors, return the error response
        return 0
    }

    # if you've made it this far, there are no errors
    return 1    

}
#======================================================
proc ::ForceFieldToolKit::GeomOpt::writeComFile {} {
    # writes the QM input file for the geometry optimization
    
    # localize relevant variables
    variable pdb
    variable com
    variable qmProc
    variable qmMem
    variable qmCharge
    variable qmMult
    variable qmRoute
    variable qmSoft $::ForceFieldToolKit::qmSoft

    # sanity check
    if { ![::ForceFieldToolKit::GeomOpt::sanityCheck] } { return }

    # call procedure to write input file for QM geometry optimization
    ::ForceFieldToolKit::${qmSoft}::writeComGeomOpt $pdb $com $qmProc $qmMem $qmCharge $qmMult $qmRoute 

}
#======================================================
proc ::ForceFieldToolKit::GeomOpt::loadLogFile {} {
    # loads the output file from the geometry optimization
    
    # localize relevant variables
    variable pdb
    variable logFile
    variable qmSoft $::ForceFieldToolKit::qmSoft

    # check to makes sure that pdb is set
    if { $pdb eq "" || ![file exists $pdb] } {
        tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "PDB file was not specified or could not be found."
        return
    }

    # make sure that logFile is set
    if { $logFile eq "" || ![file exists $logFile] } {
        tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "Cannot find optimization output file."
        return
    }

    # make sure qmSoft variable is set to the right value for the selected output file
    if {[::ForceFieldToolKit::SharedFcns::checkWhichQM $logFile]} {return}

    # call procedure to load output file from QM geometry optimization
    ::ForceFieldToolKit::${qmSoft}::readOutGeomOpt $pdb $logFile

    # message the console
    ::ForceFieldToolKit::gui::consoleMessage "Geometry optimization output file loaded"
}
#======================================================
proc ::ForceFieldToolKit::GeomOpt::writeOptPDB {} {
    # writes a new pdb file with coordinates for optimized geometry
    
    # localize relevant variables
    variable pdb
    variable logFile
    variable optPdb
    variable qmSoft $::ForceFieldToolKit::qmSoft

    # check to makes sure that pdb is set
    if { $pdb eq "" || ![file exists $pdb] } {
        tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "PDB file was not specified or could not be found."
        return
    }

    # make sure that logFile is set
    if { $logFile eq "" || ![file exists $logFile] } {
        tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "Cannot find optimization output file."
        return
    }

    # make sure qmSoft variable is set to the right value for the selected output file
    if {[::ForceFieldToolKit::SharedFcns::checkWhichQM $logFile]} {return}

    # make sure that optPdb is set
    if { ![file writable [file dirname $optPdb]] } {
        tk_messageBox -type ok -icon warning -message "Action halded on error!" -detail "Cannot write to output directory."
        return
    }
    
    # call procedure to write the optimized file as a new PDB file
    ::ForceFieldToolKit::${qmSoft}::writePDBGeomOpt $pdb $logFile $optPdb
    
    # message the console
    ::ForceFieldToolKit::gui::consoleMessage "Optimized geometry written to PDB file"
    
}
#======================================================

