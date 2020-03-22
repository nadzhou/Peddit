#
# $Id: fftk_GenZMatrix.tcl,v 1.16 2019/08/27 22:31:22 johns Exp $
#
#======================================================
namespace eval ::ForceFieldToolKit::GenZMatrix:: {

    variable psfPath
    #variable pdbPath
    variable outFolderPath
    variable basename
    
    variable donList
    variable accList
    variable atomLabels
    variable vizSpheresDon
    variable vizSpheresAcc
    
    variable qmProc
    variable qmMem
    variable qmRoute
    variable qmCharge
    variable qmMult

    variable qmSoft $::ForceFieldToolKit::qmSoft
}
#======================================================
proc ::ForceFieldToolKit::GenZMatrix::init {} {

    # IO variables
    variable psfPath
    #variable pdbPath
    variable outFolderPath
    variable basename

    # Hydrogen Bonding Variables    
    variable donList
    variable accList

    # QM Input File Variables   
    variable qmProc
    variable qmMem
    variable qmRoute
    variable qmCharge
    variable qmMult
    
    # Initialize
    set psfPath {}
    set pdbPath {}
    set outFolderPath {}
    set basename {}
    set donList {}
    set accList {}
    variable qmSoft $::ForceFieldToolKit::qmSoft

    # initialize the gaussian defaults
    ::ForceFieldToolKit::${qmSoft}::resetDefaultsGenZMatrix
}
#======================================================
proc ::ForceFieldToolKit::GenZMatrix::sanityCheck {} {
    # checks to see that appropriate information is set prior to running
    
    # returns 1 if all input is sane
    # returns 0 if there is a problem
    
    # localize relevant GenZMatrix variables
    variable psfPath
    #variable pdbPath
    set pdbPath $::ForceFieldToolKit::Configuration::geomOptPDB
    variable outFolderPath
    variable basename
    
    variable donList
    variable accList
    
    variable qmProc
    variable qmMem
    variable qmRoute
    variable qmCharge
    variable qmMult
    
    # local variables
    set errorList {}
    set errorText ""
    
    # checks
    # make sure that psfPath is entered and exists
    if { $psfPath eq "" } {
        lappend errorList "No PSF file was specified."
    } else {
        if { ![file exists $psfPath] } { lappend errorList "Cannot find PSF file." }
    }
    
    # make sure that pdbPath is entered and exists
    if { $pdbPath eq "" } {
        lappend errorList "No PDB file was specified."
    } else {
        if { ![file exists $pdbPath] } { lappend errorList "Cannot find PDB file." }
    }
    
    # make sure that outFolderPath is specified and writable
    if { $outFolderPath eq "" } {
        lappend errorList "No output path was specified."
    } else {
        if { ![file writable $outFolderPath] } { lappend errorList "Cannot write to output path." }
    }
    
    # make sure that basename is not empty
    if { $basename eq "" } { lappend errorList "No basename was specified." }
    
    # it's OK if donor and/or acceptor lists are emtpy, nothing will be written
    
    # validate gaussian settings (not particularly vigorous validation)
    # qmProc (processors)
    if { $qmProc eq "" } { lappend errorList "No processors were specified." }
    if { $qmProc <= 0 || $qmProc != [expr int($qmProc)] } { lappend errorList "Number of processors must be a positive integer." }
    # qmMem (memory)
    if { $qmMem eq "" } { lappend errorList "No memory was specified." }
    if { $qmMem <= 0 || $qmMem != [expr int($qmMem)]} { lappend errorList "Memory must be a positive integer." }
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
proc ::ForceFieldToolKit::GenZMatrix::genZmatrix {} {
    # computes the geometry-dependent position of water and
    # writes QM input files for optimizing water interaction

    # initialize some variables
    variable outFolderPath
    variable basename
    variable donList
    variable accList
    variable qmProc
    variable qmMem
    variable qmRoute
    variable qmCharge
    variable qmMult
    variable qmSoft $::ForceFieldToolKit::qmSoft

    # run sanity check
    if { ![::ForceFieldToolKit::GenZMatrix::sanityCheck] } { return }

    ::ForceFieldToolKit::${qmSoft}::genZmatrix $outFolderPath $basename $donList $accList $qmProc $qmMem $qmRoute $qmCharge $qmMult 
    
}
#==========================================================================
# proc ::ForceFieldToolKit::GenZMatrix::writeZmat moved to QMGaussian file
#=========================================================================
# proc ::ForceFieldToolKit::GenZMatrix::placeprobe moved to QMGaussian file
#=========================================================================
proc ::ForceFieldToolKit::GenZMatrix::writeExceptionZMats { aName aInd gnames atom_info } {
	# checks if C=O, P=O, or S=O which require two additional interaction
	# files at alternate water positions (120 degrees off axis)

	# localize some required variables
	variable outFolderPath
	variable basename
	variable qmProc
	variable qmMem
	variable qmRoute
	variable qmCharge
	variable qmMult
        variable qmSoft $::ForceFieldToolKit::qmSoft

	# lookup info from atoms A and B
	set aSel [atomselect top "index $aInd"]
	set bondlistA [lindex [$aSel getbonds] 0]
	set aGname [lindex $gnames $aInd]
	set aElem [string index $aGname 0]

	set bInd [lindex $bondlistA 0]
	set bSel [atomselect top "index $bInd"]
	set bondlistB [lindex [$bSel getbonds] 0]
	set bGname [lindex $gnames $bInd]
	set bElem [string index $bGname 0]

	# find a valid C atom (and associated information)
	set cInd {}
	foreach ele $bondlistB {
		# make sure that C != A, and is non-linear
		if { $ele == $aInd } {
			continue
		} elseif { [expr {abs([measure angle [list $aInd $bInd $ele]])}] >= 178.0 } {
			continue
		} else {
			set cInd $ele
		}
	}
	if { ![llength $cInd] } { return }
	set cSel [atomselect top "index $cInd"]
	set cGname [lindex $gnames $cInd]

	# check if exception case of X=O
	if { $aElem eq "O" && ($bElem eq "C" || $bElem eq "P" || $bElem eq "S" || $bElem eq "N") && [llength $bondlistA] == 1 } {

                # call procedure to write single point input files for water interaction
                ::ForceFieldToolKit::${qmSoft}::write120filesWI $outFolderPath $basename $aName $aInd $atom_info $aGname $bGname $cGname $qmProc $qmMem $qmRoute $qmCharge $qmMult 

		# clean up
		$cSel delete
	}

	# clean up and return
	$aSel delete; $bSel delete
	return
}
#======================================================
proc ::ForceFieldToolKit::GenZMatrix::writeSPfiles {} {
    # writes single point energy files required for charge optimization
    # hard coded for HF/6-31G* and MP2/6-31G*

    # localize some variables
    variable psfPath
    #variable pdbPath
    set pdbPath $::ForceFieldToolKit::Configuration::geomOptPDB
    variable outFolderPath
    variable basename
    variable qmProc
    variable qmMem
    variable qmCharge
    variable qmMult
    variable qmSoft $::ForceFieldToolKit::qmSoft

    # write compound sp GAU file
    
    # assign atom names and gather x,y,z for output com file
    mol new $psfPath
    mol addfile $pdbPath
    set Gnames {}
    set atom_info {}

    # call procedure to write single point input files for water interaction
    ::ForceFieldToolKit::${qmSoft}::writeSPfilesWI $outFolderPath $basename $qmProc $qmMem $qmCharge $qmMult 

}
#======================================================
proc ::ForceFieldToolKit::GenZMatrix::loadCOMFile { comfile } {
    # New QM Input loader
    variable qmSoft $::ForceFieldToolKit::qmSoft
    ::ForceFieldToolKit::${qmSoft}::loadCOMFile $comfile
}
#======================================================
proc ::ForceFieldToolKit::GenZMatrix::loadLOGFile { logfile } {
    # New QM Output loader
    variable qmSoft $::ForceFieldToolKit::qmSoft
    ::ForceFieldToolKit::${qmSoft}::loadLOGFile $logfile
}
#======================================================
