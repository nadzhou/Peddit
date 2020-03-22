#===========================================================================================================
#===========================================================================================================
# ffTK QM Module Support for ORCA
#===========================================================================================================
#===========================================================================================================
namespace eval ::ForceFieldToolKit::ORCA {

}
#===========================================================================================================
# GEOMETRY OPTIMIZATION
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::writeComGeomOpt { pdb com qmProc qmMem qmCharge qmMult qmRoute } {

    mol new $pdb

    # NEW 02/01/2019: Checking "element" is defined
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    set Gnames {}
    set atom_info {}

    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element]
        $temp delete
    }

    # open the output com file
    set outfile [open $com w]
    # write the header ORCA
    puts $outfile "${qmRoute}"
    ::ForceFieldToolKit::ORCA::setMemProc $outfile $qmProc $qmMem
    puts $outfile "\%output"
    puts $outfile "  PrintLevel Mini"
    puts $outfile "  Print \[ P_Mulliken \] 1"
    puts $outfile "  Print \[ P_AtCharges_M \] 1"
    puts $outfile "end"
    puts $outfile ""
    puts $outfile "\%coords"
    puts $outfile "  CTyp xyz"
    puts $outfile "  Charge $qmCharge"
    puts $outfile "  Mult $qmMult"
    puts $outfile "  Units Angs"
    puts $outfile "  coords"
    puts $outfile ""
    # write the coordinates
    foreach atom_entry $atom_info {
       puts $outfile "[lindex $atom_entry 0] [format %16.8f [lindex $atom_entry 1]] [format %16.8f [lindex $atom_entry 2]] [format %16.8f [lindex $atom_entry 3]]"
    }
    # end lines to terminate
    puts $outfile "  end"
    puts $outfile "end"

    # clean up
    close $outfile
    mol delete top

}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::readOutGeomOpt { pdb logFile {final "false"}} {
    # load output file from geometry optimization

    # load the pdb file, followed by the coordinates from gaussian log
    set molId [mol new $pdb]
    set inFile [open $logFile r]

    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "CARTESIAN COORDINATES (ANGSTROEM)" } {
            # jump to coordinates
            gets $inFile
            # read coordinates
            set coords {}
            while { ![string match "" [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 1 3]
            }
            # add a new frame, set the coords 
            mol addfile $pdb
            for {set i 0} {$i < [llength $coords]} {incr i} {
                set temp [atomselect $molId "index $i"]
                $temp set x [lindex $coords $i 0]
                $temp set y [lindex $coords $i 1]
                $temp set z [lindex $coords $i 2]
                $temp delete
            }
            unset coords
        } else {
            continue
        }
    }

    # clean up
    close $inFile

    # Additional functionality, returns molId with final frame
    if { $final } {
        mol top $molId
        animate delete beg 0 end [expr {[molinfo $molId get numframes] - 2}]
        return $molId
    }

}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::writePDBGeomOpt { pdb logFile optPdb } {
    # write the optimized structure to as new PDB file

    # load the pdb and load the coords from the log file
    set molId [mol new $pdb]

    # NEW 02/01/2019: Checking "element" is defined
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    # parse the geometry optimization output file
    set inFile [open $logFile r]

    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "CARTESIAN COORDINATES (ANGSTROEM)" } {
            # jump to coordinates
            gets $inFile
            # read coordinates
            set coords {}
            while { ![string match "" [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 1 3]
            }
            # add a new frame, set the coords 
            mol addfile $pdb
            for {set i 0} {$i < [llength $coords]} {incr i} {
                set temp [atomselect $molId "index $i"]
                $temp set x [lindex $coords $i 0]
                $temp set y [lindex $coords $i 1]
                $temp set z [lindex $coords $i 2]
                $temp delete
            }
            unset coords
        } else {
            continue
        }
    }

    # write the new coords to file
    [atomselect $molId all] writepdb $optPdb

    # clean up 
    close $inFile
    mol delete $molId

}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::resetDefaultsGeomOpt {} {
    # resets the gaussian settings to default

    # reset to default
    set ::ForceFieldToolKit::GeomOpt::qmProc 1
    set ::ForceFieldToolKit::GeomOpt::qmMem 1
    set ::ForceFieldToolKit::GeomOpt::qmCharge 0
    set ::ForceFieldToolKit::GeomOpt::qmMult 1
    set ::ForceFieldToolKit::GeomOpt::qmRoute "\! MP2 6-31G* TightSCF opt"; # for ORCA

}

#===========================================================================================================
# CHARGE PARAMETRIZATION: WATER INTERACTION
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::genZmatrix {outFolderPath basename donList accList qmProc qmMem qmRoute qmCharge qmMult } {
 
    # First call the gaussian proc to create the files with the aux flag  
    set ORCAFLAG "true"
    set gauRoute "# HF/6-31G"
    set gauName "$basename-aux"
    ::ForceFieldToolKit::Gaussian::genZmatrix $outFolderPath $gauName $donList $accList $qmProc $qmMem $gauRoute $qmCharge $qmMult  

    set numAtoms [molinfo top get numatoms]
    #========#
    # DONORS #
    #========#
    foreach donorAtom $donList {
       
        set atomName [[atomselect top "index $donorAtom"] get name]     
        set outname [file join $outFolderPath ${basename}-DON-${atomName}.inp]
        set outfile [open $outname w]  

        puts $outfile "$qmRoute"
        ::ForceFieldToolKit::ORCA::setMemProc $outfile $qmProc $qmMem
        puts $outfile ""
        puts $outfile "\%geom"
        puts $outfile " ConnectFragments" ;# testing new contrains
        puts $outfile "   \{ 1 2 O $donorAtom $numAtoms \}"  ;# index of Ow is equal to numatoms
        puts $outfile " end"
        puts $outfile " Constraints"
        puts $outfile "   \{ B $donorAtom $numAtoms C \}"  ;# index of Ow is equal to numatoms
        if { [[atomselect top "index $donorAtom"] get element] eq "H" } {
             puts $outfile "   \{ D * $numAtoms [[atomselect top "index $donorAtom"] getbonds] * C \}"       ;# all the dihedrals *-Ow-
        } elseif { [[atomselect top "index $donorAtom"] get element] eq "C" } {
             puts $outfile "   \{ D * $numAtoms $donorAtom * C \}"       ;# all the dihedrals *-Ow-
        }
        puts $outfile " end"
        puts $outfile " invertConstraints true"
        puts $outfile "end"
        puts $outfile ""
        puts $outfile "\%coords"
        puts $outfile "  CTyp xyz"
        puts $outfile "  Charge $qmCharge"
        puts $outfile "  Mult $qmMult"
        puts $outfile "  Units Angs"
        puts $outfile "  coords"
        puts $outfile ""
    	close $outfile

        ::ForceFieldToolKit::ORCA::auxXYZwriter [file join $outFolderPath ${gauName}-DON-${atomName}.gau] $outname    

        file delete [file join $outFolderPath ${gauName}-DON-${atomName}.gau] ;# Delete the Gaussian files created
    }

    #===========#
    # ACCEPTORS #
    #===========#
    foreach accAtom $accList {

        set atomName [[atomselect top "index $accAtom"] get name]
        set outname [file join $outFolderPath ${basename}-ACC-${atomName}.inp]
        set outfile [open $outname w]

        puts $outfile "$qmRoute"
        ::ForceFieldToolKit::ORCA::setMemProc $outfile $qmProc $qmMem
        puts $outfile ""
        puts $outfile "\%geom"
        puts $outfile " ConnectFragments" ;# testing new contrains
        puts $outfile "   \{ 1 2 O $accAtom $numAtoms \}"  ;# index of Hw is equal to numatoms
        puts $outfile " end"
        puts $outfile " Constraints"
        puts $outfile "   \{ B $accAtom $numAtoms C \}"  ;# index of Hw is equal to numatoms
        puts $outfile "   \{ D * [expr {$numAtoms + 1}] $accAtom * C \}"       ;# all the dihedrals *-Ow-AccAtom-*
        puts $outfile "    end"
        puts $outfile "  invertConstraints true"
        puts $outfile "end"
        puts $outfile ""
        puts $outfile "\%coords"
        puts $outfile "  CTyp xyz"
        puts $outfile "  Charge $qmCharge"
        puts $outfile "  Mult $qmMult"
        puts $outfile "  Units Angs"
        puts $outfile "  coords"
        puts $outfile ""
        close $outfile

        ::ForceFieldToolKit::ORCA::auxXYZwriter [file join $outFolderPath ${gauName}-ACC-${atomName}.gau] $outname    

        file delete [file join $outFolderPath ${gauName}-ACC-${atomName}.gau] ;# Delete the Gaussian files created
    }

}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::auxXYZwriter { gauFile outname } {
        
    set outfile [open $outname a]  
    # Read gaussian file with QMTools
    ::QMtool::load_file $gauFile com

    # File loaded to VMD
    set xyzList [[atomselect top "not atomicnumber 0"] get index]
    set xyzSize [llength $xyzList] 
    foreach ind [lrange $xyzList 0 [expr {$xyzSize - 4} ]] { ;# Loop only the molecule index. Fragment (1).
        set auxAtom [atomselect top "index $ind"]
        puts $outfile "[$auxAtom get element](1) [format "%7s  %7s  %7s" [$auxAtom get x] [$auxAtom get y] [$auxAtom get z]]"
    }
    foreach ind [lrange $xyzList [expr {$xyzSize - 3} ] end] { ;# Loop for the H2O index. Fragment (2).
        set auxAtom [atomselect top "index $ind"]
        puts $outfile "[$auxAtom get element](2) [format "%7s  %7s  %7s" [$auxAtom get x] [$auxAtom get y] [$auxAtom get z]]"
    }
    puts $outfile "  end"
    puts $outfile "end"
    mol delete top
    close $outfile
}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::write120filesWI { outFolderPath basename aName aInd atom_info aGname bGname cGname qmProc qmMem qmRoute qmCharge qmMult } {
    # writes single point energy files required for charge optimization
    # first step is to call the Gaussian counterpart with ORCAFLAG true
    set ORCAFLAG "true"
    set gauRoute "# HF/6-31G"
    set gauName "$basename-aux"
    ::ForceFieldToolKit::Gaussian::write120filesWI $outFolderPath $gauName $aName $aInd $atom_info $aGname $bGname $cGname $qmProc $qmMem $gauRoute $qmCharge $qmMult
    
    set numAtoms [molinfo top get numatoms]
    
    # write two slightly different files
    foreach altPos {"a" "b"} dihed {180 0} {

        set outname [file join $outFolderPath ${basename}-ACC-${aName}-120${altPos}.inp]
        set outfile [open $outname w]
        
        puts $outfile "$qmRoute"
        ::ForceFieldToolKit::ORCA::setMemProc $outfile $qmProc $qmMem
        puts $outfile ""
        puts $outfile "\%geom"
        puts $outfile " ConnectFragments" ;# testing new contrains
        puts $outfile "   \{ 1 2 O $aInd $numAtoms \}"  ;# index of Hw is equal to numatoms
        puts $outfile " end"
        puts $outfile " Constraints"
        puts $outfile "   \{ B $aInd $numAtoms C \}"  ;# index of Hw is equal to numatoms
        puts $outfile "   \{ D * [expr {$numAtoms + 1}] $aInd * C \}"       ;# all the dihedrals *-Ow-AccAtom-*
        puts $outfile " end"
        puts $outfile " invertConstraints true"
        puts $outfile "end"
        puts $outfile ""
        puts $outfile "\%coords"
        puts $outfile "  CTyp xyz"
        puts $outfile "  Charge $qmCharge"
        puts $outfile "  Mult $qmMult"
        puts $outfile "  Units Angs"
        puts $outfile "  coords"
        puts $outfile ""
    	close $outfile
       
        ::ForceFieldToolKit::ORCA::auxXYZwriter [file join $outFolderPath ${gauName}-ACC-${aName}-120${altPos}.gau] $outname    

        file delete [file join $outFolderPath ${gauName}-ACC-${aName}-120${altPos}.gau] ;# Delete the Gaussian files created
    }
}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::writeSPfilesWI { outFolderPath basename qmProc qmMem qmCharge qmMult } {
    # writes single point energy files required for charge optimization
    # hard coded for HF/6-31G* and MP2/6-31G*
   
    # NEW 02/01/2019: Checking "element" is defined
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    } 
    mol delete top

    # Write the HF Single Point File
    # open output file
    set outfile [open [file join $outFolderPath ${basename}-sp-HF.inp] w]

    # write the header
    puts $outfile "\! HF 6-31G* TightSCF"
    ::ForceFieldToolKit::ORCA::setMemProc $outfile $qmProc $qmMem
    puts $outfile ""
    puts $outfile "\%coords"
    puts $outfile "  CTyp xyz"
    puts $outfile "  Charge $qmCharge"
    puts $outfile "  Mult $qmMult"
    puts $outfile "  Units Angs"
    puts $outfile "  coords"
    puts $outfile ""
    # write the cartesian coords
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
    }
    puts $outfile "  end"
    puts $outfile "end"

    close $outfile

    # Write the MP2 Single Point File
    # open output file
    set outfile [open [file join $outFolderPath ${basename}-sp-MP2.inp] w]

    # write the header
    puts $outfile "\! MP2 6-31G* TightSCF"
    ::ForceFieldToolKit::ORCA::setMemProc $outfile $qmProc $qmMem
    puts $outfile ""
    puts $outfile "\%coords"
    puts $outfile "  CTyp xyz"
    puts $outfile "  Charge $qmCharge"
    puts $outfile "  Mult $qmMult"
    puts $outfile "  Units Angs"
    puts $outfile "  coords"
    puts $outfile ""
    # write the cartesian coords
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
    }
    puts $outfile "  end"
    puts $outfile "end"

    close $outfile

    # Write TIP3P Single Point File
    set outfile [open [file join $outFolderPath wat-sp.inp] w]
    puts $outfile "\! HF 6-31G* TightSCF"
    puts $outfile ""
    puts $outfile "* gzmt $qmCharge $qmMult"
    puts $outfile "O"
    puts $outfile "H 1 0.9572"
    puts $outfile "H 1 0.9572 2 104.52"
    puts $outfile "*"

    close $outfile

}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::loadCOMFile { comfile } {
    # New QM Input loader

    set inFile [open $comfile r]
    set tempFile [open temp.xyz w]
    set route [string trim [gets $inFile]]

    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "%coords" } {
            # read CTyp
            set coordType [lindex [string trim [gets $inFile]] 1]
            # read coordinates
            set coords {}
            gets $inFile  ;# this would be the charge
            gets $inFile  ;# this would be the multiplicity
            gets $inFile  ;# this would be the units
            gets $inFile  ;# extra coords menu
            gets $inFile  ;# extra empty line
            if { $coordType eq "xyz" } {
                while { [regexp {[A-Z]} [set inLine [string trim [gets $inFile]]]] } {
                    lappend coords [lrange $inLine 0 3]
                }
            } elseif { $coordType eq "inp" } {  ;# to do
            } elseif { $coordType eq "gzmt" } { ;# to do
            }
        } elseif { [lrange $inLine 0 1] eq "* gzmt" } { ;# only for the wat-sp.inp
             set coords {}
             # first atom at 0,0,0
             lappend coords "O 0.0000 0.0000 0.0000" 
             gets $inFile ;# jump first atom
             set firstbond [lindex [string trim [gets $inFile]] 2]
             lappend coords "H $firstbond 0.0000 0.0000"
             set nextline [string trim [gets $inFile]]
             set yval [ expr { [lindex $nextline 2] * sin([lindex $nextline 4]*3.14159265359/180) } ]
             set xval [ expr { [lindex $nextline 2] * cos([lindex $nextline 4]*3.14159265359/180) } ]
             lappend coords "H $xval $yval 0.0000" 
        }
    }
    set numAtoms [llength $coords]
    puts $tempFile $numAtoms
    puts $tempFile ""

    for {set i 0} {$i < $numAtoms} {incr i} {
         puts $tempFile [lindex $coords $i]
    }

    close $tempFile 
    set molId [mol new temp.xyz] ;# load the temporary xyz file in VMD
    file delete temp.xyz ;# remove the temporary xyz file 
    close $inFile
    return $molId

}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::loadLOGFile { logfile } {
    # New QM Output loader
    set inFile [open $logfile r]

    ## Create a temporary xyz file from the ORCA input
    set firstXYZ 1
    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "CARTESIAN COORDINATES (ANGSTROEM)" } {
            # jump to the coordinates
            gets $inFile
            # read coordinates
            set coords {}
            while { [regexp {[A-Z]} [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 0 3]
            }
            if { $firstXYZ eq 1 } { ;# for the first geom opt iteration, create xyz, fill and load it in VMD
                set tempFile [open temp.xyz w] ;# create the temp.xyz
                set numAtoms [llength $coords]
                puts $tempFile $numAtoms
                puts $tempFile ""
                for {set i 0} {$i < $numAtoms} {incr i} {
                     puts $tempFile [lindex $coords $i]
                }
                close $tempFile 
                set molId [mol new temp.xyz] ;# load the temporary xyz file in VMD
                set firstXYZ 0
            } elseif { $firstXYZ eq 0 } { ;# for the following geom opt iterations addfiles
                # add a new frame, set the coords 
                mol addfile temp.xyz
                for {set i 0} {$i < [llength $coords]} {incr i} {
                    set temp [atomselect $molId "index $i"]
                    $temp set x [lindex $coords $i 1]
                    $temp set y [lindex $coords $i 2]
                    $temp set z [lindex $coords $i 3]
                    $temp delete
                }
            }
            unset coords 
        }
    }
    file delete temp.xyz ;# remove the temporary xyz file 
###########################
# Testing refiting if ORCA
###########################
    set sel2 [atomselect $molId "index < [expr { $numAtoms - 3 }]"]
    set all2 [atomselect $molId all]
    set newid [mol new $::ForceFieldToolKit::Configuration::geomOptPDB]
    set sel1 [atomselect $newid "index < [expr { $numAtoms - 3 }]"]
    $all2 move [measure fit $sel2 $sel1]
    mol delete $newid
###########################
    close $inFile
    return $molId
}

#===========================================================================================================
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::getscf_ChargeOpt { file simtype } {
    # Get SCF energies from ORCA Output 
   set scfenergies {}

   set fid [open $file r]

   set hart_kcal 1.041308e-21; # hartree in kcal
   set mol 6.02214e23;

   set num 0
   set ori 0
   set tmpscf {}
   set optstep 0
   set scanpoint 0

   while {![eof $fid]} {
      set line [string trim [gets $fid]]

      # Stop reading on errors
      # if {[string match "Error termination*" $line]} { puts $line; return $scfenergies }

      # We only read Link0
      if {[string match "****ORCA TERMINATED NORMALLY****" $line]} { variable normalterm 1; break }
            
      if {[string match "Total Energy*" $line]} {
         set scf [lindex $line 3]
         set scfkcal [expr {$scf*$hart_kcal*$mol}]
         set tmpscf [list $num $scfkcal]
         if {[llength $tmpscf]} { lappend scfenergies $tmpscf; set tmpscf {} }
         incr num
      }
   }
   close $fid
   if {[llength $tmpscf]} { lappend scfenergies $tmpscf }
   
   return $scfenergies

}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::getMolCoords_ChargeOpt { file numMolAtoms } {

   set fid [open $file r]
   set coordlist {}
   while { ![eof $fid] } {
      set inLine [string trim [gets $fid]]
      if { $inLine eq "*** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***" } {
         # jump to coordinates
         gets $fid
         gets $fid
         gets $fid
         gets $fid
         gets $fid
         # read coordinates
         while { [regexp {[A-Z]} [set inLine [string trim [gets $fid]]]] } {
            lappend coordlist [lrange $inLine 1 3]
         }
      } else {
          continue
      } 
  }

  close $fid

  return [lrange $coordlist 0 [expr $numMolAtoms - 1]]

}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::getWatCoords_ChargeOpt { file } {

   set fid [open $file r]
   set coordlist {}
   set atomnames {}
   while { ![eof $fid] } {
      set inLine [string trim [gets $fid]]
      if { $inLine eq "*** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***" } {
         # jump to coordinates
         gets $fid
         gets $fid
         gets $fid
         gets $fid
         gets $fid
         # read coordinates
         while { [regexp {[A-Z]} [set inLine [string trim [gets $fid]]]] } {
            lappend atomnames [lindex $inLine 0]
            lappend coordlist [lrange $inLine 1 3]
         }
      } else {
          continue
      } 
  }
  set numAtoms [llength $atomnames]
  set Hcount 0
  for {set i [expr $numAtoms - 3]} {$i < $numAtoms} {incr i} {
      set name [lindex $atomnames $i]
      if { [string match "O*" $name] } {
         set Ocoord [lindex $coordlist $i]
      } elseif { [string match "H*" $name] && $Hcount == 1} {
         set H2coord [lindex $coordlist $i]
         set Hcount 2
      } elseif { [string match "H*" $name] && $Hcount == 0} {
         set H1coord [lindex $coordlist $i]
         set Hcount 1
      }
  }

  close $fid

  set coords [list $Ocoord $H1coord $H2coord]
  return $coords
}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::getDipoleData_ChargeOpt { file } {
    # initialize some variables
    set coords {}
    set qmVec {}
    set qmMag {}
    # open the file for reading
    set inFile [open $file r]

    # parse the output file
    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "CARTESIAN COORDINATES (ANGSTROEM)" } {
            gets $inFile
            # read coordinates
            while { [regexp {[A-Z]} [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 1 3]
            }
        } elseif { $inLine eq "DIPOLE MOMENT" } {
            for {set i 0} {$i < 5} {incr i} { gets $inFile }
            set inLine [string trim [gets $inFile]]
            set qmVec [list [lindex $inLine 4] [lindex $inLine 5] [lindex $inLine 6]]
            gets $inFile
            gets $inFile
            set inLine [string trim [gets $inFile]]
            set qmMag [lindex $inLine 3]
        } else {
            continue
        }
    }
    close $inFile
    return [list $coords $qmVec $qmMag]
}

#======================================================
proc ::ForceFieldToolKit::ORCA::resetDefaultsGenZMatrix {} {
    # resets ORCA settings to default for water interaction tab

    # reset
    set ::ForceFieldToolKit::GenZMatrix::qmProc 1
    set ::ForceFieldToolKit::GenZMatrix::qmMem 1
    set ::ForceFieldToolKit::GenZMatrix::qmCharge 0
    set ::ForceFieldToolKit::GenZMatrix::qmMult 1
    set ::ForceFieldToolKit::GenZMatrix::qmRoute "\! HF 6-31G* TightSCF opt"

}
#======================================================
proc ::ForceFieldToolKit::ORCA::auxFit_ChargeOpt { refmolid resName } {
    set sel2 [atomselect $refmolid "resname $resName"]
    set all2 [atomselect $refmolid all]
    set newid [mol new $::ForceFieldToolKit::Configuration::geomOptPDB]
    set sel1 [atomselect $newid "resname $resName"]
    $all2 move [measure fit $sel2 $sel1]
    mol delete $newid
}

#===========================================================================================================
# CHARGE PARAMETRIZATION: ELECTROSTATIC POTENTIAL
#===========================================================================================================

proc ::ForceFieldToolKit::ORCA::writeGauFile_ESP { chk qmProc qmMem qmCharge qmMult qmRoute gau } {

    variable atom_info
    variable atom_mas
    variable atom_rad

    ### make sure that the chk file is from ORCA or that a pdb file was provided ###
    set noORCA 1
    set inFile [open $chk r]
    # read and match chk file
    while { ![eof $inFile] } {
        set line [string trim [string map { \" {} } [gets $inFile]]]
        # Check if ORCA
        if { [string match -nocase "*O   R   C   A*" $line] } {
               set noORCA 0
               break
        }
    }
    close $inFile
    # Give error message if no suitable file was found
    if { $noORCA && [file extension $chk] ne ".pdb" } {
	tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "The structure file was not recognized.\
		 Please provide an ORCA output file or a PDB file (.pdb)."
	return
    }

    ### Get atomic coordinates ###
    if { [file extension $chk] eq ".pdb" } {
    # if a pdb file was provided extract the atomic information from there.
        mol new $chk
        # NEW 02/01/2019: Checking "element" is defined
        ::ForceFieldToolKit::SharedFcns::checkElementPDB
        # Get coordinates from VMD 
        set atom_info {}
        for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
            set temp [atomselect top "index $i"]
            lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
            $temp delete
        }
        # remove loaded molecule
        mol delete top
    } else {
    # otherwise extract the information from ORCA output file
        set inFile [open $chk r]
        while { ![eof $inFile] } {
            set inLine [string trim [gets $inFile]]
            # jump to the coordinates of only the final optimized geometry
            if {[string first "FINAL ENERGY EVALUATION AT THE STATIONARY POINT" $inLine] ne -1 } {
                # Skip 5 lines
                for {set i 1} {$i <= 5} {incr i} {gets $inFile}
                # read atom and coordinates
                set atom_info {}
                while { [regexp {[A-Z]} [set inLine [string trim [gets $inFile]]]] } {
                    lappend atom_info [lrange $inLine 0 3]
                }
            } else {
                continue
            }
        }
        close $inFile
    }

    ### Files preparation ###
    # Generate 4 files. First file is the SP .inp ORCA file. Second file is a temporary one.
    # Other two are scripts to run the SP and ESP calculations.
    set SPinp [join [list [file dirname $gau] "/1-SP."  [file tail $gau]] ""] 
    set tmp   [join [list [file dirname $gau] "/TMP."   [file tail [file rootname $gau]] ".xyz"] ""]
    set SPsh  [join [list [file dirname $gau] "/1-SP."  [file tail [file rootname $gau]] ".sh"] ""]
    set ESPsh [join [list [file dirname $gau] "/2-ESP." [file tail [file rootname $gau]] ".sh"] ""]

    set gauFile [open $SPinp w]
    set tmpFile [open $tmp w]
    puts $tmpFile [llength $atom_info]
    puts $tmpFile "Temporary xyz file"
    puts "Preparing file for ESP calculation in ORCA:"
    puts ""
    # write input for single point calculation in ORCA
    puts $gauFile "${qmRoute}"
    ::ForceFieldToolKit::ORCA::setMemProc $gauFile $qmProc $qmMem
    puts $gauFile ""
    puts $gauFile "\%coords"
    puts $gauFile "  CTyp xyz"
    puts $gauFile "  Charge $qmCharge"
    puts $gauFile "  Mult $qmMult"
    puts $gauFile "  Units Angs"
    puts $gauFile "  coords"
    puts $gauFile ""
    # write the coordinates for ORCA and for a temporary .xyz file needed to get the atomic radii.
    foreach atom_entry $atom_info {
       puts $gauFile "  [lindex $atom_entry 0]\t[format %10.6f [lindex $atom_entry 1]]\t[format %10.6f [lindex $atom_entry 2]]\t[format  %10.6f [lindex $atom_entry 3]]"
       puts $tmpFile "  [lindex $atom_entry 0]\t[format %10.6f [lindex $atom_entry 1]]\t[format %10.6f [lindex $atom_entry 2]]\t[format  %10.6f [lindex $atom_entry 3]]"
    }
    # end lines to terminate
    puts $gauFile "  end"
    puts $gauFile "end"
    close $gauFile
    close $tmpFile
    puts " Input for single point energy calculation written in: $SPinp"

    # Get the atomic radii needed for the generation of the ESP grid
    mol new $tmp
    set atom_rad {}
    set atom_mas {}
    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_rad [$temp get radius]
        lappend atom_mas [$temp get  mass ]
        $temp delete
    }
    mol delete top

    # Set names for .sh files
    set job [join [list "1-SP."  [file tail [file rootname $gau]]] ""]
    set jobesp [join [list "2-ESP."  [file tail [file rootname $gau]]] ""]

    # write input to calculate the SP energy in ORCA: call "orca" script
    set ext [file extension $gau];# use same extension as the one provided by the user for the iput file
    set command "\${DIR_ORCA}/orca \${DIR}/\${JOB}${ext} > \${DIR}/\${JOB}.out &"
    set comment "Single point job name"
    ::ForceFieldToolKit::ORCA::makesh_ESP $SPsh $command $job $comment
    puts "                                                       $SPsh"

    # write input to calculate the ESP in ORCA: call "orca_vpot" script
    set command "\${DIR_ORCA}/orca_vpot \${DIR}/\${JOB}.gbw \${DIR}/\${JOB}.scfp \${DIR}/\${JOB2}_grid.xyz \${DIR}/\${JOB2}_ESP.out"
    set comment "Previous single point job name"
    set commentesp "Electrostatic potential calculation job name"
    ::ForceFieldToolKit::ORCA::makesh_ESP $ESPsh $command $job $comment $jobesp $commentesp
    puts " Input for electrostatic potential calculation written in: $ESPsh"
    # Generate the grid at which to calculate the potential. Needed by ORCA.
    #::ForceFieldToolKit::ORCA::genGridFile_ESP_msms $gau $atom_info $atom_rad $atom_mas   ;# this procedure uses MSMS, which may not be installed in all workstations
    ::ForceFieldToolKit::ORCA::genGridFile_ESP $gau $atom_info $atom_rad $atom_mas  ;# this procedure uses "measure sasa", which is internal in VMD but it is less accurate that MSMS

    file delete $tmp
    puts "Done."

}
#===========================================================================================================

proc ::ForceFieldToolKit::ORCA::writeDatFile_ESP { inputName gauLog } {

        set count 1
        set auScale 0.52917720


        # open output file and read the number of fit centers
        set logFile [open $gauLog r]
        set line [gets $logFile]
	# check if a proper file was provided (not so robust check).
	# ORCA ESP output contains the number of fitting centers as first line. Check if first line is an integer.
	if { ![string is integer -strict $line] } {
            tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "The following QM output file was not recognized: ${gauLog}"
	    close $logFile
	    return
	}
	# Get number of atoms
	mol new $::ForceFieldToolKit::Configuration::geomOptPDB	
	set nAtoms [molinfo top get numatoms]
	mol delete top
	# Calculate number of fitting centers
        set nFitCenters [expr $line - $nAtoms]


	# Check if second line in ORCA ESP output is exactly a list of 4 items (to avoid that users load the grid file instead of the real output)
        set line [gets $logFile]
	if {[llength $line] != 4} {
            tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "The following QM output file was not recognized: ${gauLog}"
            close $logFile
            return
        }

        # name the file based on the inputName. Open file and start writing
        set datName ${inputName}.dat
        set datFile [open $datName w]
        puts $datFile "  $nAtoms  $nFitCenters"

	# Write coordinates of atoms first (no fit value)
        set formatStr "                   %s   %s   %s"
	for {set i 0} {$i<$nAtoms} {incr i} {
        	set line2 [lrange $line 0 2] 
		# The coordinates in the ORCA output are alreadz in AU units, no need to convert
        	puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $line2 0] ] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $line2 1] ] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $line2 2]]]
        set line [gets $logFile]
	}

        set formatStr "   %s   %s   %s   %s"
        for {set i 0} { $i < $nFitCenters } {incr i} {
                set line2 [lrange [string map { \" {} } $line] 0 3]         
                # The coordinates in the ORCA output are alreadz in AU units, no need to convert
                puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $line2 3] ] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $line2 0] ] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $line2 1] ] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $line2 2]]]
        set line [gets $logFile]
        }

	close $logFile
	close $datFile
}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::genGridFile_ESP_msms { gau atom_info atom_rad atom_mas } {
	# Generate a grid of points that will be used by ORCA to generate an electrostatic potential at each point.
	# This procedure uses MSMS software, which must be installed in the workstation.

	# Set specific variables
	set density 25
	set layers  10
	set increm  0.09 ;# in Ang
	set PI 3.141592653589
	set auScale 0.52917720
	set ScaleProbeRad 1.0
	set ScaleAtomRad 1.0
	set atom_rad_scaled {}
	set totalpoints 0
        set incrval 0.7 ;# do not start from VDW surface (incrval=0), otherwise values will greatly differ from the ones obtained with Gaussian

	# Set name of output grid file and open it for writing
	set outName [join [list [file dirname $gau] "/2-ESP." [file tail [file rootname $gau]] "_grid.xyz"] ""]
	set tmpName [join [list [file dirname $gau] "/2-ESP." [file tail [file rootname $gau]] "_grid.tmp"] ""]
	set outGrid [open $outName w]
	set tmpGrid [open $tmpName w]
	
	# Start writing temporary file with xyz coordinates
	puts -nonewline " Generating grid for ESP calculation...WAIT"
	# First add atomic coordinates
        foreach atom_entry $atom_info {
        	puts $tmpGrid "[format %10.6f [expr [lindex $atom_entry 1]/$auScale]]\t[format %10.6f [expr [lindex $atom_entry 2]/$auScale]]\t[format  %10.6f [expr [lindex $atom_entry 3]/$auScale]]"
		set totalpoints [expr $totalpoints + 1]
        }

        # Calculate weighted average radius of molecule
        set sum_mass 0
	set av_rad 0
        for {set n 0} {$n < [llength $atom_rad]} {incr n} {
        	set sum_mass [expr $sum_mass + [lindex $atom_mas $n]]
        	set av_rad   [expr $av_rad + [lindex $atom_rad $n]*[lindex $atom_mas $n]]
        }
	set av_rad [expr $av_rad/$sum_mass]
	# Scale the atomic radii
	for {set n 0} {$n < [llength $atom_rad]} {incr n} {
        	set tmprad [expr [lindex $atom_rad $n]*$ScaleAtomRad]
        	lappend atom_rad_scaled $tmprad
        }

	# Calculate the grid points using MSMS
	set tmpInpFile  [join [list [file dirname $gau] "/tmp_msms.xyzr" ] ""]
	set tmpOutFile  [join [list [file dirname $gau] "/tmp_msms_OUT" ] ""]
	# Loop over all layer requested
	for {set i 0} {$i < $layers} {incr i} {
	puts -nonewline "."
	set inFile [open $tmpInpFile w]
	foreach xyz $atom_info r $atom_rad_scaled {
		puts $inFile "[format %10.6f [lindex $xyz 1]]\t[format %10.6f [lindex $xyz 2]]\t[format %10.6f [lindex $xyz 3]]\t[format %8.4f [expr $r + $incrval]]"
	}
	close $inFile
	exec msms -if $tmpInpFile -of $tmpOutFile -no_area -probe_radius [expr $av_rad*$ScaleProbeRad] -density $density -one_cavity 1 1

	# Write coordinates in temporary file and convert to AU.
	set rMSMS [open $tmpOutFile.vert r]
	for {set j 0} {$j<3} {incr j} {set inLine [string trim [gets $rMSMS]]}
 	set Npoints [lindex $inLine 0]	
	for {set count 0} {$count < $Npoints} {incr count} {	
		set inLine [string trim [gets $rMSMS]]
        	puts $tmpGrid "[format %10.6f [expr [lindex $inLine 0]/$auScale]]\t[format %10.6f [expr [lindex $inLine 1]/$auScale]]\t[format  %10.6f [expr [lindex $inLine 2]/$auScale]]"
		set totalpoints [expr $totalpoints + 1]
	}
	set incrval [expr $incrval + $increm]

	}
	close $tmpGrid

	# Start writing in file. First write total number of points.
	puts $outGrid $totalpoints
	# Then add atomic coordinates
	set xyzFile [open $tmpName r]
        for {set i 0} {$i<$totalpoints} {incr i} {
		set inLine [gets $xyzFile]
        	puts $outGrid $inLine
	}
	close $xyzFile
	close $outGrid
	puts "Done."
	puts ""
	puts "Clean up..."
	file delete $tmpName
	file delete $tmpInpFile
	file delete ${tmpOutFile}.face
	file delete ${tmpOutFile}.vert

}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::genGridFile_ESP { gau atom_info atom_rad atom_mas } {
        # Generate a grid of points that will be used by ORCA to generate an electrostatic potential at each point.
        # This procedure uses "measure sasa", which is in VMD.

        # Set specific variables
        set density 1200
        set layers  10
        set increm  0.09 ;# in Ang
        set PI 3.141592653589
        set auScale 0.52917720
        set ScaleProbeRad 1.0
        set ScaleAtomRad 1.0
        set atom_rad_scaled {}
        set totalpoints 0
        set incrval 0.6 ;# do not start from VDW surface (incrval=0), otherwise values will greatly differ from the ones obtained with Gaussian

        # Set name of output grid file and open it for writing
        set outName [join [list [file dirname $gau] "/2-ESP." [file tail [file rootname $gau]] "_grid.xyz"] ""]
        set tmpName [join [list [file dirname $gau] "/2-ESP." [file tail [file rootname $gau]] "_grid.tmp"] ""]
        set outGrid [open $outName w]
        set tmpGrid [open $tmpName w]

	# Write a temporary xyz file to load the molecule. Necessary to use 'measure sasa'.
        set tmpXYZ [join [list [file dirname $gau] "/2-ESP." [file tail [file rootname $gau]] "_coor.xyz"] ""]
        set xyzFile [open $tmpXYZ w]
	puts $xyzFile [llength $atom_info]
	puts $xyzFile "temporary file"
        foreach xyz $atom_info {
                puts $xyzFile "[format %5s [lindex $xyz 0]]\t[format %10.6f [lindex $xyz 1]]\t[format %10.6f [lindex $xyz 2]]\t[format %10.6f [lindex $xyz 3]]"
        }
        close $xyzFile
	mol new $tmpXYZ

        # Start writing temporary file with xyz coordinates
        puts -nonewline " Generating grid for ESP calculation...WAIT"
        # First add atomic coordinates
        foreach atom_entry $atom_info {
                puts $tmpGrid "[format %10.6f [expr [lindex $atom_entry 1]/$auScale]]\t[format %10.6f [expr [lindex $atom_entry 2]/$auScale]]\t[format  %10.6f [expr [lindex $atom_entry 3]/$auScale]]"
                set totalpoints [expr $totalpoints + 1]
        }


        # Calculate the grid points using SASA
        # Make selection of all atoms in molecule.
        set all [atomselect top "all"]
        set TotPoints 0
        for {set i 0} {$i<$layers} {incr i} {
            puts -nonewline "."
            set sasaval [measure sasa $incrval $all -points sasapoints -samples $density]
            display update ui
            set incrval [expr $incrval + $increm]
            set density [expr int( $density +  $density*($incrval/0.2)*(0.01)) ]
            set Npoints [llength $sasapoints]
            # Write coordinates in temporary file and convert to AU.
            for {set count 0} {$count < $Npoints} {incr count} {
                    puts $tmpGrid "[format %10.6f [expr [lindex $sasapoints $count 0]/$auScale]]\t[format %10.6f [expr [lindex $sasapoints $count 1]/$auScale]]\t[format  %10.6f [expr [lindex $sasapoints $count 2]/$auScale]]"
                    incr totalpoints 
            }
            set incrval [expr $incrval + $increm]

            set TotPoints [expr $TotPoints + [llength $sasapoints] ]
        }
        close $tmpGrid


        # Start writing in file. First write total number of points.
        puts $outGrid $totalpoints
        # Then add atomic coordinates
        set xyzFile [open $tmpName r]
        for {set i 0} {$i<$totalpoints} {incr i} {
                set inLine [gets $xyzFile]
                puts $outGrid $inLine
        }
	close $xyzFile
        close $outGrid
        puts "Done."
        puts ""
        puts "Clean up..."
	mol delete top
        file delete $tmpName
        file delete $tmpXYZ

}
#===========================================================================================================

proc ::ForceFieldToolKit::ORCA::makesh_ESP { wFile cmd job1 comment1 args } {
# Write a .sh file to facilitate running ORCA jobs.

set shFile [open $wFile w]

# Modify directory for ORCA
puts $shFile "#!/bin/bash
# Please add the full path to your ORCA (DIR_ORCA) installation if it is not already in your \$PATH.
DIR_ORCA=/Scr/mariano/QM_Software/orca_4_1_1_linux_x86-64_openmpi303/
"
puts $shFile "# ${comment1}
JOB=${job1}"

if {[string trimleft $args] ne ""} {
puts $shFile "# [lindex ${args} 1]
JOB2=[lindex ${args} 0]"
}

puts $shFile "# Current directory
DIR=\$PWD
# Actual command to run the calculation
${cmd}"

close $shFile
exec chmod +x $wFile
}
#======================================================
proc ::ForceFieldToolKit::ORCA::resetDefaultsESP {} {

    # resets the ORCA Settings to default values
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmProc 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmMem 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmCharge 0
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmMult 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmRoute "\! HF 6-31G* TightSCF Grid4 keepdens"

}

#===========================================================================================================
# BOND/ANGLE PARAMETRIZATION
#===========================================================================================================

proc ::ForceFieldToolKit::ORCA::WriteComFile_GenBonded { geomCHK com qmProc qmMem qmRoute qmCharge qmMult psf pdb lbThresh } {
    # make sure qmSoft variable is set to the right value for the selected output file
    if {[::ForceFieldToolKit::SharedFcns::checkWhichQM $geomCHK]} {return}

    # NEW 02/26/2019: Use readOutGeomOpt to get optimized geometry
    set newid [::ForceFieldToolKit::ORCA::readOutGeomOpt $pdb $geomCHK "true"]

    #Load optimized pdb
    # set newid [mol new $pdb]

    # NEW 02/01/2019:
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    #Set auxiliar list for xyz coordinates
    set xyzList [[atomselect $newid all] get index] 

    #Open the output inp file
    # set com $::ForceFieldToolKit::GenBonded::com
    set inFile [open $com w]
   
    puts $inFile "$qmRoute"
    ::ForceFieldToolKit::ORCA::setMemProc $inFile $qmProc $qmMem
    puts $inFile "\%base\"initial\""
    puts $inFile ""
    puts $inFile "\%coords"
    puts $inFile "  CTyp xyz"
    puts $inFile "  Charge $qmCharge"
    puts $inFile "  Mult $qmMult"
    puts $inFile "  Units Angs"
    puts $inFile "  coords"
    puts $inFile ""
    foreach ind $xyzList {
        set auxAtom [atomselect $newid "index $ind"]
        puts $inFile [format "%3s  %7s  %7s  %7s" [$auxAtom get element] [$auxAtom get x] [$auxAtom get y] [$auxAtom get z]]
    }
    puts $inFile "  end"
    puts $inFile "end"

    # Second part of the input converts the hessian to internal coordinates
    
    puts $inFile ""
    puts $inFile "\$new_job" ;# This line adds another job to the input file, it runs after the first one.
    puts $inFile ""
    puts $inFile "$qmRoute"
    ::ForceFieldToolKit::ORCA::setMemProc $inFile $qmProc $qmMem
    #puts $inFile "\%base\"final\""
    puts $inFile ""
    puts $inFile "\%geom inhess read"
    puts $inFile "      inhessname \"initial.hess\" "
    puts $inFile "      PRINTINTERNALHESS true"
    puts $inFile "end"
    puts $inFile ""
    puts $inFile "\%coords"
    puts $inFile "  CTyp xyz"
    puts $inFile "  Charge $qmCharge"
    puts $inFile "  Mult $qmMult"
    puts $inFile "  Units Angs"
    puts $inFile "  coords"
    puts $inFile ""
    foreach ind $xyzList {
        set auxAtom [atomselect $newid "index $ind"]
        puts $inFile [format "%3s  %7s  %7s  %7s" [$auxAtom get element] [$auxAtom get x] [$auxAtom get y] [$auxAtom get z]]
    }
    puts $inFile "  end"
    puts $inFile "end"
    
    close $inFile

    # Write a .sh file to facilitate running ORCA jobs.

    set wFile [join [list [file rootname $com] ".sh"] "" ]
    set shFile [open $wFile w]
    set nodirfile [lindex [file split [file rootname $com]] end]
    
    puts $shFile "#!/bin/bash"
    puts $shFile ""
    puts $shFile "# Please add the path to your ORCA (DIR_ORCA) installation if it is not already in your \$PATH."
    puts $shFile "DIR_ORCA=\"/Scr/mariano/QM_Software/orca_4_1_1_linux_x86-64_openmpi313/\""
    puts $shFile ""
    puts $shFile "# Name of the files" 
    puts $shFile "Project=$nodirfile"
    puts $shFile "# Current directory"
    puts $shFile "DIR=\$PWD"
    puts $shFile "# Actual command to run the calculation"
    puts $shFile "\$DIR_ORCA/orca \$DIR/\$Project.inp > \$DIR/\$Project.out &"
    
    close $shFile
    file attributes $wFile -permissions ugo+x
    return
}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::zmatqm_BondAngleOpt { debug debugLog hessLogID hessLog {asktypeList 0} } {

    # localize forcefieldtoolkit debugging variables
#    variable debug
#    variable debugLog

    # set hessLogID [mol new $psf]
    set natoms [molinfo $hessLogID get numatoms]
    set nbonds 0
    set nangles 0
    set ndiheds 0
    set nimprops 0
    set zmatqm {0}
    set flag "Q"
    set scan {}
    set readOptParam 0
    set typeList {} ;# only for the get_inthessian_kcal proc

    set inFile [open "[file rootname $hessLog].out"]
    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "***        THE OPTIMIZATION HAS CONVERGED     ***" } {
            # jump to the coordinates
            for { set i 0 } { $i < 11 } { incr i } { gets $inFile }
            while { [regexp {[0-9]} [set inLine [string trim [gets $inFile]]]] } {
                 set typePar [string index [lindex $inLine 1] 0]
                 lappend $typeList $typePar
                 switch $typePar {
                   B { set typePar bond; incr nbonds }
                   A { set typePar angle; incr nangles }
                   D { set typePar dihed; incr ndiheds }
                   T { set typePar dihed; incr ndiheds }
                   O { set typePar imprp; incr nimprops }
                   X { continue }
                   Y { continue }
                   Z { continue }
                 }
                 if { [string match "*bond" $typePar] } {
                     set indexlist [list [string index [lindex $inLine 3] 0] [string index [lindex $inLine 2] 0]]
                     set val [lindex $inLine 7]
                     lappend zmatqm [list "R$nbonds" $typePar $indexlist $val {{} {}} $flag $scan]

                 } elseif { [string equal "angle" $typePar] } {
                     set indexlist [list [string index [lindex $inLine 2] 0] [string index [lindex $inLine 3] 0] [string index [lindex $inLine 4] 0]]
                     set val [lindex $inLine 8]
                     lappend zmatqm [list "A$nangles" $typePar $indexlist $val {{} {} {} {}} $flag $scan]

                 } elseif {[string equal "dihed" $typePar]} {
                     set indexlist [list [string index [lindex $inLine 2] 0] [string index [lindex $inLine 3] 0] [string index [lindex $inLine 4] 0] [string index [lindex $inLine 5] 0]]
                     set val [lindex $inLine 9]
                     lappend zmatqm [list "D$ndiheds" $typePar $indexlist $val {{} {} {}} $flag $scan]
                 }
            }
            set ncoords [expr { $nbonds + $nangles + $ndiheds + $nimprops }]
            set havepar 1
            set havefc 1 ;# not sure where this comes from
            lset zmatqm 0 [list $natoms $ncoords $nbonds $nangles $ndiheds $nimprops $havepar $havefc]
            foreach ele $zmatqm { puts $ele }
        }
    }

    if { $asktypeList } { 
        return $typeList ;# info required for get_inthessian_kcal proc 
    } else {
        return $zmatqm ;# normal return for the proc 
        #return [list $zmatqm $hessLogID]
    }

}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::get_inthessian_kcal_BondAngleOpt { hessLogID hessLog } { 

    set fid [open "[file rootname $hessLog]_internal.hess" r]
    set dimHess 0
    set hess {}
    set inthessian_kcal {}

    while {![eof $fid]} {
        set line [gets $fid]
        # jump to $hessian line
        if {[string first "\$hessian" $line]>=0} {
            set dimHess [lindex [set line [gets $fid]] 0]
            for { set j 0 } { $j <= [expr { $dimHess / 5 }] } { incr j } {
                gets $fid ;# skip first line of each matrix block
                for { set i 0 } { $i < $dimHess } { incr i } {
                    set line [gets $fid]
                    set row [lindex $line 0]
                    set rowdata [lrange $line 1 end]
                    set currow [lindex $hess $row]
                    set currow [concat $currow $rowdata]
                    if { $j == 0 } {
                        lappend hess $currow
                    } else {
                        lset hess $row $currow
                    }
                }
            }
        }

    }

    # calling the zmatqm proc to obtain the typeList for the appropiate dim values 
    set debug 0
    set debugLog ""
    set asktypeList 1
    set typeList [ zmatqm_BondAngleOpt $debug $debugLog $hessLogID $hessLog $asktypeList ]

    # convert hess from hartree*bohr2 to kcal*A2
    for { set i 0 } { $i < $dimHess } { incr i } {
        set dim1 1
        if { [lindex $typeList $i] eq "B" } { set dim1 0.529177249 }
        set rowdata [lrange [lindex $hess $i] 0 $i]
        set rowdata_kcal {}
        for { set j 0 } { $j < [llength $rowdata] } { incr j } {
            set dim2 1
            if { [lindex $typeList $j] eq "B" } { set dim2 0.529177249 }
            set kcal_data [expr { [lindex $rowdata $j]*0.5*1.041308e-21*6.02214e23*0.943*0.943/($dim1*$dim2) }]
            lappend rowdata_kcal $kcal_data
        }
        lappend inthessian_kcal $rowdata_kcal
    }

    return $inthessian_kcal
}
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::resetDefaultsGenBonded {} {
    # resets the QM settings to the default values
    set ::ForceFieldToolKit::GenBonded::qmProc 1
    set ::ForceFieldToolKit::GenBonded::qmMem 1
    set ::ForceFieldToolKit::GenBonded::qmRoute "\! MP2 6-31G* TightSCF opt NumFreq"

    # Reset name of output QM file.
    set ::ForceFieldToolKit::GenBonded::com "hess.inp"

}
#===========================================================================================================


#===========================================================================================================
# DIHEDRALS PARAMETRIZATION
#===========================================================================================================

proc ::ForceFieldToolKit::ORCA::buildFiles_GenDihScan { dihData outPath basename qmProc qmCharge qmMem qmMult qmRoute psf pdb } {

    mol new $psf; mol addfile $pdb
    # NEW 02/01/2019:
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    set Gnames {}
    set atom_info {}
    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element]
        $temp delete
    }

    # cycle through each dihedral to scan
    set scanCount 1
    foreach dih $dihData {
        # In ORCA the index are 0-based, differently than in Gaussian. Therefor here use the 0-based
        set zeroInds [lindex $dih 0]

        # Measure current value of dihedral
        set dihVal [format "%.6f" [measure dihed $zeroInds] ]

        # negative scan
        # open the output file
        set outfile [open ${outPath}/${basename}.scan${scanCount}.neg.inp w]
        
        # Write route section
        puts $outfile "$qmRoute"
        ::ForceFieldToolKit::ORCA::setMemProc $outfile $qmProc $qmMem
        puts $outfile ""
        # Write dihedral variable to scan
        puts $outfile "\%geom Scan"
        # The number of steps is increasing by 1 to match with the Gaussian output, where the initial energy is considered as well as a step.
        puts $outfile "         D $zeroInds = [format "%10.6f" $dihVal], [format "%10.6f" [expr $dihVal - [lindex $dih 1]]], [expr int([lindex $dih 1]/[lindex $dih 2] + 1)]" 
        puts $outfile "         end"
        puts $outfile "      end"
        puts $outfile ""
        # Write other settings
        puts $outfile "\%coords"
        puts $outfile "  CTyp xyz"
        puts $outfile "  Charge $qmCharge"
        puts $outfile "  Mult $qmMult"
        puts $outfile "  Units Angs"
        puts $outfile "  coords"
        puts $outfile ""
        # write coords
        foreach atom_entry $atom_info {
            puts $outfile "[lindex $atom_entry 0] [format "%10.6f" [lindex $atom_entry 1]] [format "%10.6f" [lindex $atom_entry 2]] [format "%10.6f" [lindex $atom_entry 3]]"
        }
        puts $outfile "  end"
        puts $outfile "end"

        close $outfile

        # positive scan
        # open the output file
        set outfile [open ${outPath}/${basename}.scan${scanCount}.pos.inp w]

        # Write route section
        puts $outfile "$qmRoute"
        ::ForceFieldToolKit::ORCA::setMemProc $outfile $qmProc $qmMem
        puts $outfile ""
        # Write dihedral variable to scan
        puts $outfile "\%geom Scan"
        # The number of steps is increasing by 1 to match with the Gaussian output, where the initial energy is considered as well as a step.
        puts $outfile "         D $zeroInds = [format "%10.6f" $dihVal], [format "%10.6f" [expr $dihVal + [lindex $dih 1]]], [expr int([lindex $dih 1]/[lindex $dih 2] + 1)]"
        puts $outfile "         end"
        puts $outfile "      end"
        puts $outfile ""
        # Write other settings
        puts $outfile "\%coords"
        puts $outfile "  CTyp xyz"
        puts $outfile "  Charge $qmCharge"
        puts $outfile "  Mult $qmMult"
        puts $outfile "  Units Angs"
        puts $outfile "  coords"
        puts $outfile ""
        # write coords
        foreach atom_entry $atom_info {
            puts $outfile "[lindex $atom_entry 0] [format "%10.6f" [lindex $atom_entry 1]] [format "%10.6f" [lindex $atom_entry 2]] [format "%10.6f" [lindex $atom_entry 3]]"
        }
        puts $outfile "  end"
        puts $outfile "end"

        close $outfile

        incr scanCount
    }

    # clean up
    mol delete top
}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::readLog_TorExplor { log } {

    # initialize some variables
    set indDef {}; set stepEn {}; set stepCoords {}

    # open the log file for reading
    set infile [open $log r]

    # Read the ORCA .out file until finding the Relaxed Surface Scan section. Always check if EOF.
    while { ! [ regexp {(Relaxed Surface Scan)} [set inline [gets $infile]] ] && ! [eof $infile] } {
           continue
    }
    ::ForceFieldToolKit::ORCA::checkEOF $log $infile "relaxed surface scan"
    # keep reading until finding the dihedral being scanned
    while { [lindex [set inline [gets $infile]] 0] ne "Dihedral" && ! [eof $infile] } {
           continue
    }
    ::ForceFieldToolKit::ORCA::checkEOF $log $infile "dihedral indices"
    # Get number of steps for a better control.
    set scanDihSteps [lindex $inline end]
     # Get dihedral definition (already 0-based in ORCA)
    regexp {\(.*\)} $inline match1             ;# string between parentheses
    regsub -all {[(),]} $match1 " " match2   ;# substitute parenteses and comma with a space
    # Final list with only one space between numbers. Reverse list to match same output as Gaussian.
    # For some reason Gaussian give the inverted order in the Log with respect to the input file.
    set indDef [lreverse [regsub -all {\s+} [string trim $match2] " "] ] ;# final list with only one space between numbers

    # Continue reading through ORCA output file but in a while loop, to get all the optimized steps in the scan.
    for { set scanstep 1 } {$scanstep<=$scanDihSteps} {incr scanstep} {
        # Go to the beginning of one scan step
        while { ! [ regexp {(RELAXED SURFACE SCAN STEP)} [set inline [gets $infile]] ] && ! [eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $log $infile "a relaxed surface scan step"
        # Go to number of atoms (there are multiple matches as the number of optimization steps)
        while { ! [ regexp {(Number of atoms)} [set inline [gets $infile]] ] && ! [eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $log $infile "number of atoms"
        # Get the number of atoms to have a better control later when getting the coordinates
        set numAtoms [lindex $inline end]
        # Go to end of ione scan step optimization
        while { ! [ regexp {(THE OPTIMIZATION HAS CONVERGED)} [set inline [gets $infile]] ] && ! [eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $log $infile "optimization convergence"
        # keep reading until finding the optimized geometry in Cartesian coordinates
        while { [set inline [gets $infile]] ne "CARTESIAN COORDINATES (ANGSTROEM)" && ! [ eof $infile ] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $log $infile "coordinates"
        # clear any existing coordinates
        set currCoords {}
        # burn  1 lines 
        gets $infile
        # parse coordinates
        for {set i 0} {$i<$numAtoms} {incr i} {
             set inline [gets $infile]
             lappend currCoords [lrange $inline 1 3]
        }

        while { ! [ regexp {(FINAL SINGLE POINT ENERGY)} [set inline [gets $infile]] ] && ! [eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $log $infile "final energy"
        # Parse E(RHF) energy; convert hartrees to kcal/mol
        set currEnergy [expr {[lindex $inline end] * 627.5095}]
        # NOTE: this value corresponds always to the correct energy (HF, MP2 or DFT)

        # We've reached the end of one scan cycle.
        # Add the collected information to the master list
        lappend stepEn $currEnergy
        lappend stepCoords $currCoords 


    }; # end of cycling through Glog lines (FOR loop)

    close $infile

    # reverse data if it's a negative scan
    if { [regexp {\.neg\.} $log] } {
        set stepEn [lreverse $stepEn]
        set stepCoords [lreverse $stepCoords]
    }

    # return the parsed data
    return [list $indDef $stepEn $stepCoords]

}
#===========================================================================================================

proc ::ForceFieldToolKit::ORCA::parseGlog_DihOpt { debug debugLog GlogFile } {

    # localize forcefieldtoolkit debugging variables
    # variable debug
    # variable debugLog

    # initialize log-wide variables
    set currDihDef {}; set currDihVal {}; set currCoords {}; set currEnergy {}
    set infile [open $GlogFile r]
    set GlogData {}

    # Read the ORCA .out file until finding the Relaxed Surface Scan section. Always check if EOF.
    while { ! [ regexp {(Relaxed Surface Scan)} [set inline [gets $infile]] ] && ! [eof $infile] } {
           continue
    }
    ::ForceFieldToolKit::ORCA::checkEOF $GlogFile $infile "relaxed surface scan"
    # keep reading until finding the dihedral being scanned
    while { [lindex [set inline [gets $infile]] 0] ne "Dihedral" && ! [eof $infile] } {
           continue
    }
    ::ForceFieldToolKit::ORCA::checkEOF $GlogFile $infile "dihedral indices"
    # Get number of steps for a better control.
    set scanDihSteps [lindex $inline end]
     # Get dihedral definition (already 0-based in ORCA)
    regexp {\(.*\)} $inline match1             ;# string between parentheses
    regsub -all {[(),]} $match1 " " match2   ;# substitute parenteses and comma with a space
    # Final list with only one space between numbers. Reverse list to match same output as Gaussian.
    # For some reason Gaussian give the inverted order in the Log with respect to the input file.
    set scanDihInds [lreverse [regsub -all {\s+} [string trim $match2] " "] ] ;# final list with only one space between numbers

    if { $debug } {
        puts $debugLog "Scan dihedral FOUND:"
        puts $debugLog "\t0-based Indicies (VMD): $scanDihInds"
        puts $debugLog "--------------------------------------"
        puts $debugLog "Ind (VMD)\tCurrDih\tEnergy (kcal/mol)"
        puts $debugLog "--------------------------------------"
        flush $debugLog
    }
            
    # Continue reading through ORCA output file but in a while loop, to get all the optimized steps in the scan.
    for { set scanstep 1 } {$scanstep<=$scanDihSteps} {incr scanstep} {

        while { ! [ regexp {(RELAXED SURFACE SCAN STEP)} [set inline [gets $infile]] ] && ! [eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $GlogFile $infile "a relaxed surface scan step"
        # keep reading until finding the scanned dihedral
        while { [lindex [set inline [gets $infile]] 1] ne "Dihedral" && ! [ eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $GlogFile $infile "dihedral value"
        # parse out the current dihedral value; round to integer
        set currDihVal [expr { round([lindex [ regsub -all {\*} $inline " "] end]) }]
        while { ! [ regexp {(Number of atoms)} [set inline [gets $infile]] ] && ! [eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $GlogFile $infile "number of atoms"
        # Get the number of atoms to have a better control later when getting the coordinates
        set numAtoms [lindex $inline end]

        while { ! [ regexp {(THE OPTIMIZATION HAS CONVERGED)} [set inline [gets $infile]] ] && ! [eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $GlogFile $infile "optimization convergence"
        # keep reading until finding the optimized geometry in Cartesian coordinates
        while { [set inline [gets $infile]] ne "CARTESIAN COORDINATES (ANGSTROEM)" && ! [ eof $infile ] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $GlogFile $infile "coordinates"
        # clear any existing coordinates
        set currCoords {}
        # burn  1 lines 
        gets $infile
        # parse coordinates
        for {set i 0} {$i<$numAtoms} {incr i} {
     	 set inline [gets $infile]
             lappend currCoords [lrange $inline 1 3]
        }

        while { ! [ regexp {(FINAL SINGLE POINT ENERGY)} [set inline [gets $infile]] ] && ! [eof $infile] } {
               continue
        }
        ::ForceFieldToolKit::ORCA::checkEOF $GlogFile $infile "final energy"
        # Parse E(RHF) energy; convert hartrees to kcal/mol
        set currEnergy [expr {[lindex $inline end] * 627.5095}]
        # NOTE: this value corresponds always to the correct energy (HF, MP2 or DFT)

        # We've reached the end of one scan cycle.
        # Add the collected information to the master list
        lappend tempGlogData [list $scanDihInds $currDihVal $currEnergy $currCoords]

        if { $debug } {
            puts $debugLog "$scanDihInds\t$currDihVal\t$currEnergy"
            flush $debugLog
        }
    }; # end of cycling through Glog lines (FOR loop)

    # if the QM log file runs the scan in negative direction, reverse the order of entries
    # if not explicitely in the negative direction, preserve the order of entries
    if { [lsearch -exact [split $GlogFile \.] "neg"] != -1 } {
        foreach ele [lreverse $tempGlogData] { lappend GlogData $ele }
    } else {
        foreach ele $tempGlogData { lappend GlogData $ele }
    }

    # clean up
    close $infile

    return $GlogData
}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::checkEOF { inputname inputfile tag } {

   if { [eof ${inputfile}] } {
      tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "End of file. No match was found for ${tag} in: ${inputname}. No frames were loaded."
      return -level 2 -code 2
   }

}

#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::resetDefaultsGenDihScan {} {
    # reset QM settings for generation of dihedral scan to the default values
    #
    # set variables
    set ::ForceFieldToolKit::GenDihScan::qmProc 1
    set ::ForceFieldToolKit::GenDihScan::qmCharge 0
    set ::ForceFieldToolKit::GenDihScan::qmMem 1
    set ::ForceFieldToolKit::GenDihScan::qmMult 1
    set ::ForceFieldToolKit::GenDihScan::qmRoute "\! MP2 6-31G* TightSCF opt"

}
#===========================================================================================================

#===========================================================================================================
# SHARED PROCEDURES FOR ORCA
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::setMemProc { outfile qmProc qmMem } {
    set settings {}
    # In ORCA user specify the memory per core in MB
    set mem [expr ($qmMem*1000)/$qmProc]
    # Write memory information to file
    puts $outfile "%maxcore $mem"
    
    # Write information about number of processors if more then one core is requested
    if {$qmProc>1} {
        puts $outfile "%pal"
        puts $outfile "   nprocs $qmProc"
        puts $outfile "end"
    }

}

#===========================================================================================================
# TEMPORARY  PROCEDURES
#===========================================================================================================
proc ::ForceFieldToolKit::ORCA::tempORCAmessage {} {
     # Print message about the current implementation of ORCA
     tk_messageBox -type ok -icon warning -message "User info" -detail "ORCA is under implementation. Use only Gaussian software for the moment."
}
#===========================================================================================================
