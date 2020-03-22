#===========================================================================================================
#===========================================================================================================
# ffTK QM Module Support for Gaussian
#===========================================================================================================
#===========================================================================================================
namespace eval ::ForceFieldToolKit::Gaussian {

}
#===========================================================================================================
# GEOMETRY OPTIMIZATION
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::resetDefaultsGeomOpt {} {
    # resets the gaussian settings to default

    # reset to default
    set ::ForceFieldToolKit::GeomOpt::qmProc 1
    set ::ForceFieldToolKit::GeomOpt::qmMem 1
    set ::ForceFieldToolKit::GeomOpt::qmCharge 0
    set ::ForceFieldToolKit::GeomOpt::qmMult 1
    set ::ForceFieldToolKit::GeomOpt::qmRoute "\# MP2/6-31G* Opt=(Redundant) SCF=Tight Geom=PrintInputOrient"

}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::writeComGeomOpt { pdb com qmProc qmMem qmCharge qmMult qmRoute } {
    # write input file for geometry optimization

    # load structure
    mol new $pdb

    # NEW 02/01/2019: Checking "element" is defined
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    # assign atom names and gather x,y,z for output com file
    set Gnames {}
    set atom_info {}

    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element][expr $i+1] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    }
    # open the output com file
    set outfile [open $com w]

    # write the header
    puts $outfile "%chk=[file tail [file rootname $com]].chk"
    puts $outfile "%nproc=$qmProc"
    puts $outfile "%mem=${qmMem}GB"
    puts $outfile "$qmRoute"
    puts $outfile ""
    puts $outfile "<qmtool> simtype=\"Geometry optimization\" </qmtool>"
    puts $outfile ""
    puts $outfile "$qmCharge $qmMult"

    # write the coordinates
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [format %16.8f [lindex $atom_entry 1]] [format %16.8f [lindex $atom_entry 2]] [format %16.8f [lindex $atom_entry 3]]"
    }

    # empty line to terminate
    puts $outfile ""

    # clean up
    close $outfile
    mol delete top

}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::readOutGeomOpt { pdb logFile } {
    # load output file from geometry optimization

    # load the pdb file, followed by the coordinates from gaussian log
    set molId [mol new $pdb]
    set inFile [open $logFile r]

    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "Input orientation:" } {
            # burn the coord header
            for {set i 0} {$i < 4} {incr i} { gets $inFile }
            # read coordinates
            set coords {}
            while { ![regexp {^-*$} [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 3 5]
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

}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::writePDBGeomOpt { pdb logFile optPdb } {
    # write the optimized structure to as new PDB file
   
    # load the pdb and load the coords from the log file
    set molId [mol new $pdb]

    # NEW 02/01/2019: Checking "element" is defined
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    # parse the geometry optimization output file
    set inFile [open $logFile r]

    while { ![eof $inFile] } {
        set inLine [string trim [gets $inFile]]
        if { $inLine eq "Input orientation:" } {
            # burn the coord header
            for {set i 0} {$i < 4} {incr i} { gets $inFile }
            # read coordinates
            set coords {}
            while { ![regexp {^-*$} [set inLine [string trim [gets $inFile]]]] } {
                lappend coords [lrange $inLine 3 5]
            }
            # (re)set the coords 
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
proc ::ForceFieldToolKit::Gaussian::auxFit_ChargeOpt { refmolid resName } {
     # Empty proc, used for other QM softwares (i.e. ORCA)
}

#===========================================================================================================
# CHARGE PARAMETRIZATION: WATER INTERACTION
#===========================================================================================================

proc ::ForceFieldToolKit::Gaussian::genZmatrix {outFolderPath basename donList accList qmProc qmMem qmRoute qmCharge qmMult } {
    # NEW 02/01/2019: Checking "element" is defined
    ::ForceFieldToolKit::SharedFcns::checkElementPDB

    # computes the geometry-dependent position of water and
    # writes Gaussian input files for optimizing water interaction
    # assign Gaussian atom names and gather x,y,z for output gau file
    set Gnames {}
    set atom_info {}
    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element][expr $i+1] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    }

    #========#
    # DONORS #
    #========#
    foreach donorAtom $donList {

		set donorName [[atomselect top "index $donorAtom"] get name]

		# open output file
    	        set outname [file join $outFolderPath ${basename}-DON-${donorName}.gau]
    	        set outfile [open $outname w]

    	        # write the header
		puts $outfile "%chk=${basename}-DON-${donorName}.chk"
		puts $outfile "%nproc=$qmProc"
		puts $outfile "%mem=${qmMem}GB"
                puts $outfile "$qmRoute"
		puts $outfile ""
		puts $outfile "<qmtool> simtype=\"Geometry optimization\" </qmtool>"
		puts $outfile "${basename}-DON-${donorName}"
		puts $outfile ""
		puts $outfile "$qmCharge $qmMult"

		# write the cartesian coords for the molecule
		foreach atom_entry $atom_info {
		    puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
		}

		# build / write the zmatrix
		::ForceFieldToolKit::Gaussian::writeZmat $donorAtom donor $Gnames $outfile
		close $outfile
    }; # end DONORS

    #===========#
    # ACCEPTORS #
    #===========#
    foreach acceptorAtom $accList {
		set acceptorName [[atomselect top "index $acceptorAtom"] get name]

		# open output file
		set outname [file join $outFolderPath ${basename}-ACC-${acceptorName}.gau]
		set outfile [open $outname w]

		# write the header
		puts $outfile "%chk=${basename}-ACC-${acceptorName}.chk"
		puts $outfile "%nproc=$qmProc"
		puts $outfile "%mem=${qmMem}GB"
                puts $outfile "$qmRoute"
		puts $outfile ""
		puts $outfile "<qmtool> simtype=\"Geometry optimization\" </qmtool>"
		puts $outfile "${basename}-ACC-${acceptorName}"
		puts $outfile ""
		puts $outfile "$qmCharge $qmMult"

		# write the cartesian coords for the molecule
		foreach atom_entry $atom_info {
		    puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
		}

		# build / write the zmatrix
		::ForceFieldToolKit::Gaussian::writeZmat $acceptorAtom acceptor $Gnames $outfile
		close $outfile

		# CARBONYL EXCEPTIONS
		# there is a special exception for X=O cases (e.g. C=O, P=O, S=O)
		# where we need to write two additional files
		::ForceFieldToolKit::GenZMatrix::writeExceptionZMats $acceptorName $acceptorAtom $Gnames $atom_info

    }
}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::writeZmat { aInd class gnames outfile } {
	# builds and writes the z-matrix to file

	# passed:
	#		aInd = atom index for interacting atom
	#		class = acceptor / donor
	#		gnames = list of Gaussian-style names
	#		outfile = file handle where output is written
	
	# returns: nothing

	# make a selection for atom A (the interaction site), and find all attached atoms
	set aSel [atomselect top "index $aInd"]
	set bondlistA [lindex [$aSel getbonds] 0]

	# z-matrix will be different for n=1 and n>1 (each will have their own special cases as well)
	if { [llength $bondlistA] == 1} {
		#=======#
		# N = 1 #
		#=======#
		# when n=1 there are two cases in which only 1D scanning is possible:
		# a diatomic molecule, or when C--B--A angle is 180 degrees (i.e. linear)

		set diatomic 0; # defaults to false
		set linear 1; # defaults to true

		# check for a diatomic molecule
		if { [molinfo top get numatoms] == 2} { set diatomic 1 }
		# check if C--B--A is linear
		if { !$diatomic } {
			set bSel [atomselect top "index $bondlistA"]
			set bondlistB [lindex [$bSel getbonds] 0]
			foreach ele $bondlistB {
				# if C = A, skip
				if { $ele == $aInd } { continue }
				# check for a non-linear angle (+/- 2 degrees)
				if { [expr {abs([measure angle [list $aInd $bondlistA $ele]])}] <= 178.0 } {
					# found a non-linear C atom; unset the flag and stop looking
					set linear 0; break
				} else {
					# keep looking
					continue
				}
			}
			# clean up atom selections
			$bSel delete
		}

		if { $diatomic || $linear } {
			# if either diatomic or linear, build a zmatrix for a 1D scan
			
			# get the relevant gnames
			set aGname [lindex $gnames $aInd]
			set bGname [lindex $gnames $bondlistA]

			# for non-diatomic linear cases, we will define x in cartesian coords
			if { $linear } {
				set bSel [atomselect top "index $bondlistA"]
				set v1 [vecnorm [vecsub [measure center $aSel] [measure center $bSel]]]
				if { [lindex $v1 0] == 0 } {
					set v2 "1 0 0"
				} elseif { [lindex $v1 1] == 0 } {
					set v2 "0 1 0"
				} elseif { [lindex $v1 2] == 0 } {
					set v2 "0 0 1"
				} else {
					set v2 [list 1 [expr -1.0*[lindex $v1 0]/[lindex $v1 1]] 0]
				}
				set v2 [vecnorm $v2]
				set xPos [vecadd $v2 [measure center $aSel]]
			}

			# positioning of x is the same for donors and acceptors
			if { $diatomic } {
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s" x $aGname 1.0 $bGname 90.00]	
			} else {
				puts $outfile [format "%3s  %.4f  %.4f  %.4f" x [lindex $xPos 0] [lindex $xPos 1] [lindex $xPos 2]]
			}

			# write the rest of the z-matrix
			if { $class eq "donor" } {
				# donor
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s" Ow $aGname rAH x 90.00 $bGname 180.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s" H1w Ow 0.9572 $aGname 127.74 x 0.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s\n" H2w Ow 0.9572 $aGname 127.74 x 180.00]
				puts $outfile "rAH 2.0"
			} else {
				# acceptor
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s" Ow $aGname rAH x 90.00 $bGname 180.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s" H2w Ow 0.9572 $aGname 104.52 x   0.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s %7s  %6s\n" H1w Ow 0.9572 H2w 104.52 x 0.0]
				puts $outfile "rAH 2.0"
			}

			# clean up and return
			$aSel delete
			if { $linear } { $bSel delete }
			return
			# done with special n=1 cases (diatomic or linear)

		} else {
			# handle the n=1 case 'normally'

			# find some information about B atom
			set bInd $bondlistA
			set bSel [atomselect top "index $bInd"]
			set bondlistB [lindex [$bSel getbonds] 0]

			# find a valid C atom
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

			# make an atom selection of C atom
			set cSel [atomselect top "index $cInd"]

			# get gnames
			set aGname [lindex $gnames $aInd]
			set bGname [lindex $gnames $bInd]
			set cGname [lindex $gnames $cInd]

			if { $class eq "donor" } {
				# donor
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x $aGname 1.0 $bGname 90.00 $cGname dih]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow $aGname rAH x 90.00 $bGname 180.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w Ow 0.9572 $aGname 127.74 x 0.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9572 $aGname 127.74 x 180.00]
				puts $outfile "rAH  2.0"
				puts $outfile "dih  0.0"
			} else {
				# acceptor
				# call helper function to find probe position
				set probePos [::ForceFieldToolKit::Gaussian::placeProbe $aSel]
				# note that Gaussian doesn't like 180 degree angles, so we have to be a little clever
				# make some measurements of probe position
				set mAng [::QMtool::bond_angle $probePos [measure center $aSel] [measure center $cSel]]
				set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $cSel] [measure center $bSel]]

				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $cGname [format %3.2f $mAng] $bGname [format %3.2f $mDih]]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x H1w 1.0 $aGname 90.00 $cGname 0.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow H1w 0.9527 x 90.00 $aGname 180.00]
				puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9527 H1w 104.52 x dih]
				puts $outfile "rAH  2.0"
				puts $outfile "dih  0.0"
			}

			# clean up and return
			$aSel delete; $bSel delete; $cSel delete
			return			
		}; # end of 'normal' N=1

	} else {
		#=======#
		# N > 1 #
		#=======#

		# find some information about atom B
		set bInd [lindex $bondlistA 0]
		set bSel [atomselect top "index $bInd"]
		set cInd [lindex $bondlistA 1]
		set cSel [atomselect top "index $cInd"]

		# find gnames
		set aGname [lindex $gnames $aInd]
		set bGname [lindex $gnames $bInd]
		set cGname [lindex $gnames $cInd]

		# test if C is valid choice
		set abcAng [expr {abs([measure angle [list $aInd $bInd $cInd]])}]
		if { $abcAng < 2 || $abcAng > 178 } { set validC 0 } else {	set validC 1 }
		unset abcAng

		# find probe coords
		set probePos [::ForceFieldToolKit::Gaussian::placeProbe $aSel]
		set mAng [::QMtool::bond_angle $probePos [measure center $aSel] [measure center $bSel]]

		if { !$validC } {
			# if C is invalid, ABC are linear and we need a second dummy atom in lieu of the original C atom
			set cGname "x2"
			set x2Pos [coordtrans [trans center [measure center $aSel] axis [vecsub [measure center $cSel] [measure center $bSel]] 180.0 deg] $probePos]
			set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $bSel] $x2Pos]
			puts $outfile [format "%3s  %.4f  %.4f  %.4f" x2 [lindex $x2Pos 0] [lindex $x2Pos 1] [lindex $x2Pos 2]]
		} else {
			# C is valid, we can use it to define the dihedral
			set mDih [::QMtool::dihed_angle $probePos [measure center $aSel] [measure center $bSel] [measure center $cSel]]	
		}
		
		if { $class eq "donor" } {
			# donor
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow $aGname rAH $bGname [format %3.2f $mAng] $cGname [format %3.2f $mDih]]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x Ow 1.0 $aGname 90.00 $bGname dih]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w Ow 0.9572 $aGname 127.74 x 0.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9572 $aGname 127.74 x 180.00]
            puts $outfile "rAH  2.0"
            puts $outfile "dih  0.0"

		} else {
			# acceptor
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $bGname [format %3.2f $mAng] $cGname [format %3.2f $mDih]]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x H1w 1.0 $aGname 90.00 $bGname 0.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow H1w 0.9572 x 90.00 $aGname 180.00]
			puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9572 H1w 104.52 x dih]
			puts $outfile "rAH  2.0"
			puts $outfile "dih  0.0"
		}

		# clean up atomselections and return
		$aSel delete; $bSel delete; $cSel delete
		return
	}
}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::placeProbe { aSel } {
	# helper function for genZmatrix
	# computes the position of the first probe atom
	# (i.e. Oxygen for HB donors or Hydrogen for HB acceptors)

	# takes in atomselection for the interacting atom
	# returns x,y,z of the first probe atom (i.e. O or H)

	# setup some variables
	set mol [$aSel molid]
	set aCoord [measure center $aSel]

	set bonded [lindex [$aSel getbonds] 0]
	set numBonded [llength $bonded]

	set all [atomselect $mol "all"]
	set allCenter [measure center $all]


	set count 0
	set normavg {0 0 0}


	# compute the probe position based on geometry
	# n = 1, i.e. b---a
	if { $numBonded == 1 } {
		set temp1 [atomselect $mol "index [lindex $bonded 0]"]
		set normavg [vecsub $aCoord [measure center $temp1]]
		incr count
		$temp1 delete
	}


	# n = 2, i.e. b1---a---b2
	if { $numBonded == 2 } {
		set temp1 [atomselect $mol "index [lindex $bonded 0]"]
		set bondvec1 [vecnorm [vecsub $aCoord [measure center $temp1]]]
		set temp2 [atomselect $mol "index [lindex $bonded 1]"]
		set bondvec2 [vecnorm [vecsub $aCoord [measure center $temp2]]]

		# check for linearity, project point away from rest of molecule
        if { abs([vecdot $bondvec1 $bondvec2]) > 0.95 } {
            # check that center of molecule doesn't already lie on our line, or worse, our atom
            # if it does, we don't care about where we place the water relative to the rest of
            # molecule
            set flag 0
            if { [veclength [vecsub $allCenter $aCoord]] < 0.01 } {
                set flag 1
            } elseif { abs([vecdot $bondvec1 [vecnorm [vecsub $allCenter $aCoord]]]) > 0.95 } {
                set flag 1
            }
            if { $flag } {
                if { [lindex $bondvec1 0] == 0 } {
                   set probeCoord] "1 0 0"
                } elseif { [lindex $bondvec1 1] == 0 } {
                   set probeCoord "0 1 0"
                } elseif { [lindex $bondvec1 2] == 0 } {
                   set probeCoord "0 0 1"
                } else {
                   set probeCoord [list 1 [expr -1.0*[lindex $bondvec1 0]/[lindex $bondvec1 1]] 0]
                }
                set probeCoord [vecadd $aCoord [vecscale 2.0 [vecnorm $probeCoord]]]
                return $probeCoord
            }
            set alpha [vecdot $bondvec1 [vecsub $allCenter $aCoord]]
            # same side as mol center
            set probeCoord [vecsub $allCenter [vecscale $alpha $bondvec1]]
            # reflect, extend
            set probeCoord [vecadd $aCoord [vecscale 2.0 [vecnorm [vecsub $aCoord $probeCoord]]]]
            return $probeCoord
        } else {
            set normavg [vecadd $bondvec1 $bondvec2]
        }

		incr count
		$temp1 delete; $temp2 delete; $all delete
	}

	# n > 2, there are many cases; majority are n=3 and n=4
	if { $numBonded > 2 } {
		# cycle through to find all normals
	    for {set i 0} {$i <= [expr $numBonded-3]} {incr i} {
	    
	        set temp1 [atomselect $mol "index [lindex $bonded $i]"]
	        # normalize bond vectors first
	        set normPos1 [vecadd $aCoord [vecnorm [vecsub [measure center $temp1] $aCoord]]]
	    
	        for {set j [expr $i+1]} {$j <= [expr $numBonded-2]} {incr j} {
				set temp2 [atomselect $mol "index [lindex $bonded $j]"]
				set normPos2 [vecadd $aCoord [vecnorm [vecsub [measure center $temp2] $aCoord]]]

				for {set k [expr $j+1]} {$k <= [expr $numBonded-1]} {incr k} {
					set temp3 [atomselect $mol "index [lindex $bonded $k]"]
					set normPos3 [vecadd $aCoord [vecnorm [vecsub [measure center $temp3] $aCoord]]]

					# get the normal vector to the plane formed by the three atoms
					set vec1 [vecnorm [vecsub $normPos1 $normPos2]]
					set vec2 [vecnorm [vecsub $normPos2 $normPos3]]
					set norm [veccross $vec1 $vec2]

					# check that the normal vector and atom of interest are on the same side of the plane
					set d [expr -1.0*[vecdot $norm $normPos1]]
					if { [expr $d + [vecdot $norm $aCoord]] < 0 } {
						set norm [veccross $vec2 $vec1]
					}

					# will average normal vectors at end
					set normavg [vecadd $normavg $norm]
					incr count
					$temp3 delete
				}
				$temp2 delete
	        }
	        $temp1 delete
	    }
	}

	# finish up and return
	set normavg [vecscale [expr 1.0/$count] $normavg]
	set probeCoord [vecadd $aCoord [vecscale 2.0 [vecnorm $normavg]]]
	return $probeCoord
}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::write120filesWI { outFolderPath basename aName aInd atom_info aGname bGname cGname qmProc qmMem qmRoute qmCharge qmMult } {
    # writes single point energy files required for charge optimization
    # hard coded for HF/6-31G* and MP2/6-31G*
    
    # write two slightly different files
    foreach altPos {"a" "b"} dihed {180 0} {

    	# open output file
    	set outname [file join $outFolderPath ${basename}-ACC-${aName}-120${altPos}.gau]
    	set outfile [open $outname w]

    	# write the header
    	puts $outfile "%chk=${basename}-ACC-${aName}-120${altPos}.chk"
    	puts $outfile "%nproc=$qmProc"
    	puts $outfile "%mem=${qmMem}GB"
      puts $outfile "$qmRoute"
    	puts $outfile ""
    	puts $outfile "<qmtool> simtype=\"Geometry optimization\" </qmtool>"
    	puts $outfile "${basename}-ACC-${aName}-120${altPos}"
    	puts $outfile ""
    	puts $outfile "$qmCharge $qmMult"

    	# write the cartesian coords for the molecule
    	foreach atom_entry $atom_info {
    	    puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
    	}

    	# write custom zmatrix
    	set ang 120

    	puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" H1w $aGname rAH $bGname [format %3.2f $ang] $cGname [format %3.2f $dihed]]
    	puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" x H1w 1.0 $aGname 90.00 $bGname 0.00]
    	puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s" Ow H1w 0.9527 x 90.00 $aGname 180.00]
    	puts $outfile [format "%3s  %7s  %6s  %7s  %6s  %7s  %6s\n" H2w Ow 0.9527 H1w 104.52 x dih]
    	puts $outfile "rAH  2.0"
    	puts $outfile "dih  0.0"

    	# close up
    	close $outfile
    }
}
#===========================================================================================================

proc ::ForceFieldToolKit::Gaussian::writeSPfilesWI { outFolderPath basename qmProc qmMem qmCharge qmMult } {
    # writes single point energy files required for charge optimization
    # hard coded for HF/6-31G* and MP2/6-31G*

    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element][expr $i+1] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    } 
    mol delete top

    # Write the HF Single Point File
    # open output file
    set outfile [open [file join $outFolderPath ${basename}-sp-HF.gau] w]

    # write the header
    puts $outfile "%chk=${basename}-sp-HF.chk"
    puts $outfile "%nproc=$qmProc"
    puts $outfile "%mem=${qmMem}GB"
    puts $outfile "# HF/6-31G* SCF=Tight"
    puts $outfile ""
    puts $outfile "<qmtool> simtype=\"Single point calculation\" </qmtool>"
    puts $outfile "${basename}-sp-HF"
    puts $outfile ""
    puts $outfile "$qmCharge $qmMult"

    # write the cartesian coords
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
    }

    puts $outfile ""
    close $outfile

    # Write the MP2 Single Point File
    # open output file
    set outfile [open [file join $outFolderPath ${basename}-sp-MP2.gau] w]

    # write the header
    puts $outfile "%chk=${basename}-sp-MP2.chk"
    puts $outfile "%nproc=$qmProc"
    puts $outfile "%mem=${qmMem}GB"
    puts $outfile "# MP2/6-31G* SCF=Tight Density=Current"
    puts $outfile ""
    puts $outfile "<qmtool> simtype=\"Single point calculation\" </qmtool>"
    puts $outfile "${basename}-sp-MP2"
    puts $outfile ""
    puts $outfile "$qmCharge $qmMult"

    # write the cartesian coords
    foreach atom_entry $atom_info {
        puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
    }

    puts $outfile ""
    close $outfile

    # Write TIP3P Single Point File
    set outfile [open [file join $outFolderPath wat-sp.gau] w]
    puts $outfile "%chk=wat-sp.chk"
    puts $outfile "%nproc=1"
    puts $outfile "%mem=1GB"
    puts $outfile "# RHF/6-31G* SCF=Tight"
    puts $outfile ""
    puts $outfile "wat-sp"
    puts $outfile ""
    puts $outfile "0 1"
    puts $outfile "O1"
    puts $outfile "H2 1 0.9572"
    puts $outfile "H3 1 0.9572 2 104.52"
    puts $outfile ""
    close $outfile

}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::loadCOMFile { comfile } {
    # New QM Input loader

    set molId [mol new]
    ::QMtool::use_vmd_molecule $molId
    ::QMtool::read_gaussian_input $comfile $molId
    return $molId
}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::loadLOGFile { logfile } {
    # New QM Output loader
    
    set molId [mol new]
    ::QMtool::use_vmd_molecule $molId
    ::QMtool::read_gaussian_log $logfile $molId
    return $molId
}

#===========================================================================================================
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::getscf_ChargeOpt { file simtype } {
    # Get SCF energies from Gaussian Output o
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
      if {[string match "Error termination*" $line]} { puts $line; return $scfenergies }

      # We only read Link0
      if {[string match "Normal termination of Gaussian*" $line]} { variable normalterm 1; break }
            
      if {$simtype=="Relaxed potential scan"} {
         if {[string match "Step number * out of a maximum of * on scan point * out of *" $line]} {
            set optstep   [lindex $line 2]
            set scanpoint [lindex $line 12]
            set scansteps [lindex $line 15]
#            puts "SCAN: optstep $optstep on scan point $scanpoint out of $scansteps"
         }
      }
            
     if {[string match "SCF Done:*" $line] || [string match "Energy=* NIter=*" $line]} {
         if {[string match "SCF Done:*" $line]} {
            set scf [lindex $line 4]
         } else {
            set scf [lindex $line 1]
         }
         set scfkcal [expr {$scf*$hart_kcal*$mol}]
         if {$num==0} { set ori $scf }
         set scfkcalori [expr {($scf-$ori)*$hart_kcal*$mol}]
         # In case of a relaxed potential scan we replace the previous energy of the same scanstep,
         # otherwise we just append all new scf energies
         if {$optstep==1 || !($simtype=="Relaxed potential scan")} {
            if {[llength $tmpscf]} { lappend scfenergies $tmpscf; set tmpscf {} }
#            puts [format "%i: SCF = %f hart = %f kcal/mol; rel = %10.4f kcal/mol" $num $scf $scfkcal $scfkcalori]
         }
         set tmpscf [list $num $scfkcal]

         incr num
      }

   }
   close $fid
   if {[llength $tmpscf]} { lappend scfenergies $tmpscf }

   return $scfenergies
}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::getMolCoords_ChargeOpt { file numMolAtoms } {

   set fid [open $file r]
   ::QMtool::init_variables ::QMtool

   ::QMtool::read_gaussian_cartesians $fid qmtooltemppdb.pdb last
   file delete qmtooltemppdb.pdb
   set coordlist [lindex [::QMtool::get_cartesian_coordinates] 0]
    
   close $fid

   return [lrange $coordlist 0 [expr $numMolAtoms - 1]]
}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::getWatCoords_ChargeOpt { file } {

   set fid [open $file r]
   ::QMtool::init_variables ::QMtool

   ::QMtool::read_gaussian_cartesians $fid qmtooltemppdb.pdb last
   file delete qmtooltemppdb.pdb
   set coordlist [lindex [::QMtool::get_cartesian_coordinates] 0]
   set atomlist [::QMtool::get_atomproplist]
   set numAtoms [llength $atomlist]

   set Hcount 0
   for {set i [expr $numAtoms - 4]} {$i < $numAtoms} {incr i} {
      set name [lindex [lindex $atomlist $i] 1]
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
proc ::ForceFieldToolKit::Gaussian::getDipoleData_ChargeOpt { filename } {

    # initialize some variables
    set coords {}
    set qmVec {}
    set qmMag {}

    # open the file for reading
    set inFile [open $filename r]

    # parse the output file
    # burn lines until we find the std orientation
    while { [set inLine [string trim [gets $inFile]]] ne "Standard orientation:" } { continue }

    # once std orientation is found, burn header (4 lines)
    for {set i 0} {$i < 4} {incr i} { gets $inFile }

    # read in the coords
    while { ![regexp {^-*$} [set inLine [string trim [gets $inFile]]]] } {
        lappend coords [lrange $inLine 3 5]
    }

    # burn until find the dipole moment
    while { [set inLine [string trim [gets $inFile]]] ne "Dipole moment (field-independent basis, Debye):"} { continue }

    # parse the dipole moment
    set inLine [string trim [gets $inFile]]
    set qmVec [list [lindex $inLine 1] [lindex $inLine 3] [lindex $inLine 5]]
    set qmMag [lindex $inLine 7]

    # we're done with the output file
    unset inLine
    close $inFile
    
    return [list $coords $qmVec $qmMag]
}

#======================================================
proc ::ForceFieldToolKit::Gaussian::resetDefaultsGenZMatrix {} {
    # resets gaussian settings to default for water interaction tab

    # reset
    set ::ForceFieldToolKit::GenZMatrix::qmProc 1
    set ::ForceFieldToolKit::GenZMatrix::qmMem 1
    set ::ForceFieldToolKit::GenZMatrix::qmCharge 0
    set ::ForceFieldToolKit::GenZMatrix::qmMult 1
    set ::ForceFieldToolKit::GenZMatrix::qmRoute "# HF/6-31G* Opt=(Z-matrix,MaxCycles=100) Geom=PrintInputOrient"

}
#======================================================


#===========================================================================================================
# CHARGE PARAMETRIZATION: ELECTROSTATIC POTENTIAL
#===========================================================================================================

proc ::ForceFieldToolKit::Gaussian::writeGauFile_ESP { chk qmProc qmMem qmCharge qmMult qmRoute gau } {

    ### make sure that the chk file is from Gaussian or that a pdb file was provided ###
    # Give error message if no suitable file was found
    if { [file extension $chk] ne ".chk" && [file extension $chk] ne ".pdb" } {
        tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "The structure file was not recognized.\
                 Please provide a Gaussian Checkpoint file (.chk) or a PDB file (.pdb)."
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
        # write the .gau file
	set newCHKname "[file rootname $gau].chk"
        set gauFile [open $gau w]
        # write header
        puts $gauFile "%chk=[file tail $newCHKname]"
        puts $gauFile "%nproc=$qmProc"
        puts $gauFile "%mem=${qmMem}GB"
	# Modify route section to remove Geom=Checkpoint, which is not needed if pdb is used. Give message in tk.
	regsub -nocase {Geom=Checkpoint} $qmRoute {} qmRoute_mod
   	puts ""
        puts $gauFile "$qmRoute_mod"
        puts $gauFile ""
        puts $gauFile "<qmtool> simtype=\"ESP Calculation\" </qmtool>"
        puts $gauFile ""
        puts $gauFile "$qmCharge $qmMult"
        # write the coordinates
        foreach atom_entry $atom_info {
            puts $gauFile "[lindex $atom_entry 0] [format %16.8f [lindex $atom_entry 1]] [format %16.8f [lindex $atom_entry 2]] [format %16.8f [lindex $atom_entry 3]]"
        }
        # empty line to terminate
        puts $gauFile ""
        close $gauFile

    } else {
    # otherwise copy structural information from Gaussian checkpoint file
        # make a copy of the CHK file to prevent from overwriting the original
        set newCHKname "[file rootname $gau].chk"
        file copy -force $chk $newCHKname

         set gauFile [open $gau w]

         # write the .gau file
         puts $gauFile "%chk=[file tail $newCHKname]"
         puts $gauFile "%nproc=$qmProc"
         puts $gauFile "%mem=${qmMem}GB"
         puts $gauFile "$qmRoute"
         puts $gauFile ""
         puts $gauFile "<qmtool> simtype=\"ESP Calculation\" </qmtool>"
         puts $gauFile ""
         puts $gauFile "$qmCharge $qmMult"

         close $gauFile

    }

}
#===========================================================================================================

proc ::ForceFieldToolKit::Gaussian::writeDatFile_ESP { inputName gauLog } {
        set count 1
        set auScale 0.52917720

        # make sure qmSoft variable is set to the right value for the selected output file
        if {[::ForceFieldToolKit::SharedFcns::checkWhichQM $gauLog]} {return}

        # name the file based on the inputName
        set datName ${inputName}.dat
        set datFile [open $datName w]

        set logFile [open $gauLog r]

	# read the number of atoms
	while { [lindex [set line [string map { \" {} } [gets $logFile]] ] 0] ne "NAtoms=" } { 
		continue
        }

	set nAtoms [lindex $line 1]

	# read the number of fit centers
	while { [lrange [set line  [gets $logFile]] 1 8] ne "points will be used for fitting atomic charges" } {
		continue
        }

	set nFitCenters [lindex $line 0]
	
	puts $datFile "  $nAtoms  $nFitCenters"

	# go back to the beginning of the .log file
	seek $logFile 0
	
	# read the Atom Fit Centers
	while { [lrange [set line [string map { \" {} } [gets $logFile]]] 0 1] ne "Atomic Center" } { 
		continue
        }
   set formatStr "                   %s   %s   %s"
#	set temp [format $formatStr [lindex $line 5] [lindex $line 6] [lindex $line 7]]
   puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 5] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 6] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 7] / $auScale]]]
   while { [lrange [set line [gets $logFile]] 0 1] eq "Atomic Center" } {
       puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 5] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 6] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 7] / $auScale]]]
		continue
    }
	
	# go back to the beginning of the .log file
	seek $logFile 0
	
	# read the fit values
	while { [lrange [set line [string map { \" {} } [gets $logFile]]] 0 3] ne "Electrostatic Properties (Atomic Units)" } {
		continue
    }
	# G09 Revision B does not give the fit values, so there must be a way to check this at this point of the code and give an error.
	# Otherwise the code will loop forever.
	while { [lindex [set line [gets $logFile]] 1] ne "Fit" && ![eof $logFile] } {
		continue
    }
        if { [eof $logFile] } {
            tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "Cannot find fitting values in log file. This may happen because you are using g09 revision B. Try to use a later version."
            return
        }
	lappend fitValue [lindex $line 2]
	while { [lindex [set line [gets $logFile]] 1] eq "Fit" } {
		lappend fitValue [lindex $line 2]
		continue
    }
	
	# go back to the beginning of the .log file
	seek $logFile 0
	
	# write formatted fit values to the data file
	while { [lrange [set line [string map { \" {} } [gets $logFile]]] 0 2] ne "ESP Fit Center" } { 
		continue
   }
	set formatStr "   %s   %s   %s"
	puts -nonewline $datFile [format "   %s" [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $fitValue 0]]]
	puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 6] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 7] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 8] / $auScale]]]
	while { [lrange [set line [gets $logFile]] 0 2] eq "ESP Fit Center" } {
		puts -nonewline $datFile [format "   %s" [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [lindex $fitValue $count]]]
		puts $datFile [format $formatStr [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 6] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 7] / $auScale]] [::ForceFieldToolKit::ChargeOpt::ESP::formatRESP [expr [lindex $line 8] / $auScale]]]
		incr count
		continue
	}	
	close $logFile
	close $datFile
}
#======================================================
proc ::ForceFieldToolKit::Gaussian::resetDefaultsESP {} {

    # resets the Gaussian Settings to default values
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmProc 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmMem 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmCharge 0
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmMult 1
    set ::ForceFieldToolKit::ChargeOpt::ESP::qmRoute "#P HF/6-31G* SCF=Tight Geom=Checkpoint Pop=MK IOp(6/33=2,6/41=10,6/42=17)"

}

#===========================================================================================================
# BOND/ANGLE PARAMETRIZATION
#===========================================================================================================

proc ::ForceFieldToolKit::Gaussian::WriteComFile_GenBonded { geomCHK com qmProc qmMem qmRoute qmCharge qmMult psf pdb lbThresh } {
    # Check if a file with .chk extension was given 
    set fext [file extension $geomCHK]
    if { $fext ne ".chk" } {
        tk_messageBox -type ok -icon warning -message "Action halted on error!" -detail "File has not a \".chk\" extension. It might not be a Gaussian Checkpoint File."
        return
    }

    # load the molecule psf/pdb to get the internal coordinates
    set logID [mol new $psf]
    mol addfile $pdb $logID
    ::QMtool::use_vmd_molecule $logID
    set zmat [::QMtool::modredundant_zmat]
    
    # make a copy of the CHK file to prevent Gaussian from overwriting the original
    set newCHKname "[file rootname $com].chk"
    file copy $geomCHK $newCHKname
    
    # write the com file
    set outfile [open $com w]
    puts $outfile "%chk=[file tail $newCHKname]"
    puts $outfile "%nproc=$qmProc"
    puts $outfile "%mem=${qmMem}GB"
    puts $outfile "$qmRoute"
    puts $outfile ""

    # shamelessly stolen from qmtool
    # First delete all existing internal coordinates
    puts $outfile "B * * K"
    puts $outfile "A * * * K"
    puts $outfile "L * * * K"
    puts $outfile "D * * * * K"
    #puts $outfile "O * * * * R"
    
    set num 0
    set lbList {}
    foreach entry $zmat {
        # skip the qmtool zmat header (first line)
        if {$num==0} { incr num; continue }

        set indexes {}
        foreach ind [lindex $entry 2] {
            lappend indexes [expr {$ind+1}]
        }
        set type [string toupper [string index [lindex $entry 1] 0]]

        # check for linear angle
        if { $type eq "A" && [lindex $entry 3] > $lbThresh } {
            # angle qualifies as a "linear bend"
            set type "L"
            lappend lbList $indexes
        }

        # check if standard dihedrals are part of linear bend (undefined)
        set skipflag 0
        if { $type eq "D" && [llength $lbList] > 0 } {
            # test each linear bend in lbList against current dih indices definition
            foreach ang $lbList {
                # test forward and rev angle definitions
                if { [string match "*$ang*" $indexes] || [string match "*[lreverse $ang]*" $indexes] } {
                    # positive test -> leave this dihedral out
                    set skipflag 1
                    incr num
                    break
                }
            }
        }
        if { $skipflag } { continue }
      
        # impropers modeled as dihedrals because Gaussian ignores out-of-plane bends
        if {$type=="I"} { set type "D" }
        if {$type=="O"} { set type "D" }

        # write the entry to the input file
        puts $outfile "$type $indexes A [regsub {[QCRM]} [lindex $entry 5] {}]"
        #puts $outfile "$type $indexes $val [regsub {[QCRM]} [lindex $entry 5] {}]"
        incr num
    }

    puts $outfile ""
    close $outfile
    
}

#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::zmatqm_BondAngleOpt { debug debugLog hessLogID hessLog {asktypeList 0} } {

    # localize forcefieldtoolkit debugging variables
#    variable debug
#    variable debugLog
    
    # load the Gaussian Log files from the hessian calculation
    if { $debug } { puts -nonewline $debugLog "loading hessian log file..."; flush $debugLog }
    # set hessLogID [mol new $psf]

    ::QMtool::use_vmd_molecule $hessLogID
    ::QMtool::load_gaussian_log $hessLog $hessLogID

    if { $debug } { puts $debugLog "DONE"; flush $debugLog }

    # store internal coordinates from the hessian calculation
    set zmatqm [::QMtool::get_internal_coordinates]
    #return [list $zmatqm $hessLogID]
    return $zmatqm 

}
   
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::get_inthessian_kcal_BondAngleOpt { {hessLogID ""} {hessLog ""} } { 
   set inthessian_kcal [::QMtool::get_internal_hessian_kcal]
   return $inthessian_kcal
}
#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::resetDefaultsGenBonded {} {
    # resets the QM settings to the default values
    set ::ForceFieldToolKit::GenBonded::qmProc 1
    set ::ForceFieldToolKit::GenBonded::qmMem 1
    set ::ForceFieldToolKit::GenBonded::qmRoute "\# MP2/6-31G* Geom=(AllCheck,ModRedundant) Freq NoSymm IOp(7/33=1) SCF=Tight Guess=Read"

    # Reset name of output QM file.
    set ::ForceFieldToolKit::GenBonded::com "hess.gau"

}
#===========================================================================================================


#===========================================================================================================
# DIHEDRALS PARAMETRIZATION
#===========================================================================================================

proc ::ForceFieldToolKit::Gaussian::buildFiles_GenDihScan { dihData outPath basename qmProc qmCharge qmMem qmMult qmRoute psf pdb } {

    # assign Gaussian atom names and gather x,y,z for output com file
    mol new $psf; mol addfile $pdb

    # NEW 02/01/2019:
    ::ForceFieldToolKit::SharedFcns::checkElementPDB
    
    set Gnames {}
    set atom_info {}
    for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
        set temp [atomselect top "index $i"]
        lappend atom_info [list [$temp get element][expr $i+1] [$temp get x] [$temp get y] [$temp get z]]
        lappend Gnames [$temp get element][expr $i+1]
        $temp delete
    }

    # cycle through each dihedral to scan
    set scanCount 1
    foreach dih $dihData {
        # change 0-based indices to 1-based
        set zeroInds [lindex $dih 0]
        set oneInds {}
        foreach ind $zeroInds {
            lappend oneInds [expr {$ind + 1}]
        }

        # negative scan
        # open the output file
        set outfile [open ${outPath}/${basename}.scan${scanCount}.neg.gau w]

        # write the header
        puts $outfile "%chk=${basename}.scan${scanCount}.neg.chk"
        puts $outfile "%nproc=$qmProc"
        puts $outfile "%mem=${qmMem}GB"
        puts $outfile "$qmRoute"
        puts $outfile ""
        puts $outfile "$basename Dihedral Scan"
        puts $outfile ""
        puts $outfile "$qmCharge $qmMult"
        # write coords
       foreach atom_entry $atom_info {
           puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
       }
       # write scan
       puts $outfile ""
       puts $outfile "D $oneInds S [expr int([expr [lindex $dih 1]/[lindex $dih 2]])] [format "%.6f" [expr {-1*[lindex $dih 2]}]]"

       close $outfile

       # positive scan
       # open the output file
       set outfile [open ${outPath}/${basename}.scan${scanCount}.pos.gau w]

       # write the header
       puts $outfile "%chk=${basename}.scan${scanCount}.pos.chk"
       puts $outfile "%nproc=$qmProc"
       puts $outfile "%mem=${qmMem}GB"
       puts $outfile "$qmRoute"
       puts $outfile ""
       puts $outfile "$basename Dihedral Scan at MP2/6-31G*"
       puts $outfile ""
       puts $outfile "$qmCharge $qmMult"
       # write coords
      foreach atom_entry $atom_info {
          puts $outfile "[lindex $atom_entry 0] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]"
      }
      # write scan
      puts $outfile ""
      puts $outfile "D $oneInds S [expr int([expr [lindex $dih 1]/[lindex $dih 2]])] [format "%.6f" [lindex $dih 2]]"

      close $outfile

      incr scanCount

    }

}

#===========================================================================================================
proc ::ForceFieldToolKit::Gaussian::readLog_TorExplor { log } {

    # initialize some variables
    set indDef {}; set stepEn {}; set stepCoords {}

    # open the log file for reading
    set infile [open $log r]
    while { ![eof $infile] } {
        # read a line at a time
        set inline [string trim [gets $infile]]

        switch -regexp $inline {
            {Initial Parameters} {
                # keep reading until finding the dihedral being scanned
                # and parse out 1-based indices, convert to 0-based
                while { ![regexp {^\!.*D\(([0-9]+),([0-9]+),([0-9]+),([0-9]+)\).*Scan[ \t]+\!$} [string trim [gets $infile]] full_match ind1 ind2 ind3 ind4] && ![eof $infile] } { continue }
                if { [eof $infile] } { return }
                foreach ele [list $ind1 $ind2 $ind3 $ind4] { lappend indDef [expr {$ele - 1}] }
            }

            {Input orientation:} {
                # clear any existing coordinates
                set currCoords {}
                # burn the header
                for {set i 0} {$i<=3} {incr i} { gets $infile }
                # parse coordinates
                while { [string range [string trimleft [set line [gets $infile]] ] 0 0] ne "-" } { lappend currCoords [lrange $line 3 5] }
            }

            {SCF[ \t]*Done:} {
                # parse E(RHF) energy; convert hartrees to kcal/mol
                set currEnergy [expr {[lindex $inline 4] * 627.5095}]
                # NOTE: this value will be overridden if E(MP2) is also found
            }

            {E2.*EUMP2} {
                # convert from Gaussian notation in hartrees to scientific notation
                set currEnergy [expr {[join [split [lindex [string trim $inline] end] D] E] * 627.5095}]
                # NOTE: this overrides the E(RHF) parse from above
            }

            {Optimization completed} {
                # we've reached the optimized conformation
                lappend stepEn $currEnergy
                lappend stepCoords $currCoords
            }

            default {continue}
        }
    }

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
proc ::ForceFieldToolKit::Gaussian::parseGlog_DihOpt { debug debugLog GlogFile } {

    # localize forcefieldtoolkit debugging variables
    # variable debug
    # variable debugLog

    # initialize log-wide variables
    set currDihDef {}; set currDihVal {}; set currCoords {}; set currEnergy {}
    set infile [open $GlogFile r]
    set GlogData {}

    # read through Gaussian Log File (Glog)
    while {[eof $infile] != 1} {
        # read a line in
        set inline [gets $infile]
        # parse line
        switch -regexp $inline {
            {Initial Parameters} {
                # keep reading until finding the dihedral being scanned
                while { [lindex [set inline [gets $infile]] 4] ne "Scan" } {
                    continue
                }
                # parse out the dihedral definition (1-based indices)
                set scanDihDef [lindex $inline 2]
                # strip out four indices from D(#1,#2,#3,#4)
                set scanDihInds {}
                foreach ind [split [string range $scanDihDef 2 [expr [string length $scanDihDef] - 2]] ","] {
                   lappend scanDihInds [expr $ind - 1]
                }

                if { $debug } {
                    puts $debugLog "Scan dihedral FOUND:"
                    puts $debugLog "\tGaussian Indicies: $scanDihDef"
                    puts $debugLog "\t0-based Indicies (VMD): $scanDihInds"
                    puts $debugLog "--------------------------------------"
                    puts $debugLog "Ind (VMD)\tCurrDih\tEnergy (kcal/mol)"
                    puts $debugLog "--------------------------------------"
                    flush $debugLog
                }
            }

            {Input orientation:} {
                # clear any existing coordinates
                set currCoords {}
                # burn the header
                for {set i 0} {$i<=3} {incr i} {
                    gets $infile
                }
                # parse coordinates
                while { [string range [string trimleft [set line [gets $infile]] ] 0 0] ne "-" } {
                    lappend currCoords [lrange $line 3 5]
                }
            }

            {SCF[ \t]*Done:} {
                # parse E(RHF) energy; convert hartrees to kcal/mol
                set currEnergy [expr {[lindex $inline 4] * 627.5095}]
                # NOTE: this value will be overridden if E(MP2) is also found
            }

            {E2.*EUMP2} {
                # convert from Gaussian notation in hartrees to scientific notation
                set currEnergy [expr {[join [split [lindex [string trim $inline] end] D] E] * 627.5095}]
                # NOTE: this overrides the E(RHF) parse from above
            }

            {Optimization completed} {
                # we've reached an optimized conformation
                # keep reading until finding the scanned dihedral
                while { [lindex [set inline [gets $infile]] 2] ne $scanDihDef } {
                    continue
                }
                # parse out the current dihedral value; round to integer
                set currDihVal [expr { round([lindex $inline 3]) }]
                # add the collected information to the master list
                # lappend GlogData [list $scanDihInds $currDihVal $currEnergy $currCoords]
                lappend tempGlogData [list $scanDihInds $currDihVal $currEnergy $currCoords]

                if { $debug } {
                    puts $debugLog "$scanDihInds\t$currDihVal\t$currEnergy"
                    flush $debugLog
                }
            }
        }; # end of line parse (switch)
    }; # end of cycling through Glog lines (while)

    # if the Gaussian log file runs the scan in negative direction, reverse the order of entries
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
proc ::ForceFieldToolKit::Gaussian::resetDefaultsGenDihScan {} {
    # reset QM settings for generation of dihedral scan to the default values
    #
    # set variables
    set ::ForceFieldToolKit::GenDihScan::qmProc 1
    set ::ForceFieldToolKit::GenDihScan::qmCharge 0
    set ::ForceFieldToolKit::GenDihScan::qmMem 1
    set ::ForceFieldToolKit::GenDihScan::qmMult 1
    set ::ForceFieldToolKit::GenDihScan::qmRoute "\# opt=modredundant MP2/6-31g(d) Geom=PrintInputOrient"

}
#===========================================================================================================
#===========================================================================================================
