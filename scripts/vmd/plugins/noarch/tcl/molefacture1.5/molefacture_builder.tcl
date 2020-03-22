##
## $Id: molefacture_builder.tcl,v 1.41 2019/06/28 20:20:38 jribeiro Exp $
##

# This file contains all procs relating to molecule building, fragments, and the like

proc ::Molefacture::add_all_hydrogen { {runAutoPsf 1} } {
   variable tmpmolid
   variable protlasthyd
   variable picklist
   variable unknownRes
   variable selTop
   global env


   Molefacture::saveState
   mol reanalyze $tmpmolid
   update_openvalence


   set molList [molinfo list]

   set seltextorig "($Molefacture::notDumSel) and occupancy >= 0.8 and (not ions and not water)"

   if {[llength $picklist] > 0} {
    set seltextorig "($Molefacture::notDumSel) and index [join $picklist] and occupancy >= 0.8 and (not ions and not water)"
   }

   set seltext $seltextorig

   if {$runAutoPsf == 1 && [llength $picklist] == 0} {
    append seltext " and not (protein and not nucleic and not glycan)"
   }
   if {[llength $unknownRes] > 0} {
    append seltext " and resname [join $unknownRes]"
   }
   set selUnkAdd [atomselect $tmpmolid $seltext]
   set mindexlist [$selUnkAdd get index]
   set addHresname [lsort -unique -dictionary [$selUnkAdd get resname]]

   $selUnkAdd set $Molefacture::modFlagH 1

  $selUnkAdd delete
  set selUnkAdd ""
   set newatoms {}
   ### try again but this time with runAutoPsf 0
   ### in case a non-protein residue has a protein resname
   ### The former fails in the "and not protein" selection


   foreach mindex $mindexlist {
      set picklist $mindex
      set numhyd $Molefacture::openvalence($mindex)
      set j 0

      while {$j < $numhyd} {
        set index [add_atom_noupdate H]
        if {$index != -1} {
          lappend newatoms $index
        }

        incr j
      }

    }

    update_openvalence

    if {[llength $newatoms] > 1} {
      set newatoms [join $newatoms]
    }

    if {[llength $newatoms] > 0} {
      topo -molid "$tmpmolid" -sel "index $newatoms" guessatom radius element
      topo -molid "$tmpmolid" -sel "index $newatoms" guessatom mass element
    }
    # mol reanalyze $tmpmolid
    if {$runAutoPsf == 0} {
      ::Molefacture::updateAtomTable
    } elseif {$runAutoPsf == 1 && $addHresname != ""} {
      set selUnkAdd [atomselect $tmpmolid "resname [join $addHresname]"]
    }

    ### Is there anything left to do?
    if {[llength $addHresname] > 0 } {
      append seltext " and not resname [join $addHresname]"
    }
    if {[llength $unknownRes] > 0} {
      append seltext " and not resname [join $unknownRes]"
    }

    set seltext "($seltextorig) and not ($Molefacture::modFlagPSF == 1) and not ($Molefacture::modFlagH == 1) "

    set selKnAdd [atomselect $tmpmolid $seltext]

    if {[$selKnAdd num] == 0} {
      set runAutoPsf 0
      $selKnAdd delete
    }

    ### Use the topology information and autopsf to add the hydrogens
    if {$runAutoPsf == 1} {

      set pwd [pwd]

      cd $env(TMPDIR)
      if {[file exist molefacture_tmp] != 1} {
        file mkdir molefacture_tmp
      }
      cd molefacture_tmp
      ### use autopsf (psfgen to add the hydrogen to known residues)

      set resnames [lsort -unique -dictionary [$selKnAdd get resname]]

      if {[info exist Molefacture::deftopinfo] != 1 || $Molefacture::deftopinfo == ""} {
        Molefacture::loadDefaulTopo $resnames
      }

      if {[llength $selTop] > 0} {
        set molList [molinfo list]

        set autopsfLog ""
        $selKnAdd writepdb knownres.pdb
        $selKnAdd delete
        set mol [mol new knownres.pdb]

        set atpsfOk [catch {autopsf -mol ${mol} -inlcude "($seltext) or water" -prefix knownres -top $selTop -regen -noterm} autopsfLog ]

        if {[file exists knownres_formatted_autopsf.psf] != 1} {
          puts "Autpsf Error $autopsfLog"
          foreach mol [molinfo list] {
            if {[lsearch $molList $mol] == -1} { mol delete $mol }
          }
          cd $pwd
          ::Molefacture::add_all_hydrogen 0

          return
        } else {
          set autopsfDone 1
        }

        set newmol [mol new knownres_formatted_autopsf.psf]
        mol addfile knownres_formatted_autopsf.pdb waitfor all $newmol

        set knownsel [atomselect $newmol "all"]

        set dumSel [atomselect $tmpmolid "resname DUM or segname CAP"]
        set selist [list $knownsel  $dumSel]

        if {$selUnkAdd != ""} {
           lappend selist $selUnkAdd
        }

        set k 0
        foreach sel $selist {
          incr k
          $sel global
        }

        Molefacture::mergeSels $molList $selist

        foreach sel $selist {
          $sel delete
        }

        ::Molefacture::updateAtomTable
        Molefacture::updateprotlasth
        update_openvalence
        return
      } else {
        ::Molefacture::updateAtomTable
        $selKnAdd delete
        cd $pwd
      }
    }

    if {$runAutoPsf == 0} {
     ::Molefacture::updateAtomTable
    }


}


proc ::Molefacture::add_atom_noupdate { element } {
  #Same as add_hydrogen, but doesn't call update_openvalence, which
  #is a major timesuck for big structures
  variable tmpmolid
  
  variable atomaddlist
  variable picklist
  variable topGui
  
  set newatom ""
  ### do == 1; add isolated atom
   ### do == 2; add atom to the selected atom(s)
  set do 0
  if {[$Molefacture::atomTable size] > 0} {
    set do 2
  } else {
    set do 1
  }

   #Check to make sure we don't need more dummy hydrogens
   set dummies [atomselect $tmpmolid "$Molefacture::dumSel"]
   if {[$dummies num] == 0} {
     topo -molid "$tmpmolid" -sel "$Molefacture::notDumSel" guessatom radius element
     topo -molid "$tmpmolid" -sel "$Molefacture::notDumSel" guessatom mass element
     hrealloc "$Molefacture::notDumSel"
     $dummies delete
     set dummies [atomselect $tmpmolid "$Molefacture::dumSel"]
   }


   if {$do == 1} {


    set aindex [lindex [$dummies list] 0]

    set newatom $aindex
    set newatomsel  [atomselect $tmpmolid "index $aindex"]
     #Set the occupancy to indicate its new status
    if {$element == "H"} {
      $newatomsel set occupancy 0.5
    } else {
      $newatomsel set occupancy 0.8
    }

    $newatomsel set segid A
    $newatomsel set resid 1
    $newatomsel set resname TMP
    $newatomsel set chain A
    $newatomsel set atomicnumber [lsearch $Molefacture::periodic [string toupper $element]]
    $newatomsel set name ${element}

    $newatomsel set type $element
    $newatomsel set element $element

    $newatomsel delete
    
  } else {
    set incrdum 0
    foreach mindex $picklist {

      set aindex [lindex [$dummies list] $incrdum]
      lappend newatom $aindex
      set newatomsel  [atomselect $tmpmolid "index $aindex"]
      set mother [atomselect $tmpmolid "index $mindex"]
      set mcoor [join [$mother get {x y z}]]
      

      set numbonds [Molefacture::newbond $mindex]
      if {$numbonds <= -1} {
        return -1
      }
      ### Add bond an delete possible hydrogen atoms.
      ### In this case, this test is done before addDelBond beacuse the number of hydrogens
      ### will influence the geometry of the new added atom
      ### Also, only call update_openvalence if really needed
      if {$element != "H" && [info exist Molefacture::openvalence($mindex)] == 1} {
         ### Check if the openvalence is negative (more atoms than supposed)
         ## and delete hydrogens if possible to make it 0
         if {[expr $Molefacture::openvalence($mindex) -1] < 0} {
          #  set sel [atomselect $Molefacture::tmpmolid "index $atom and occupancy > 0.8"]
           set bonds [join [$mother getbonds]]
           if {[$mother num] > 0 && [llength $bonds] > 0} {
              set selbonds [atomselect $Molefacture::tmpmolid "index $bonds and occupancy 0.5"]
              if {[$selbonds num] > 0} {
                set picklistaux $picklist 
                set pick [list]
                foreach index [$selbonds get index] {
                  lappend pick $index
                  if {[expr $Molefacture::openvalence($mindex) + [llength $pick]] == 1} {
                    break
                  }
                }
                set picklist $pick
                ::Molefacture::del_atom_gui
                set picklist $picklistaux
                update_openvalence
                set numbonds [expr $numbonds - [llength $pick]]
              }
              $selbonds delete
           }
          
         } 
      } 


      $newatomsel set type [string toupper $element]
      $newatomsel set element $element
      $newatomsel set segid   [$mother get segid]
      $newatomsel set resid   [$mother get resid]
      $newatomsel set resname [$mother get resname]
      $newatomsel set chain [$mother get chain]

      $newatomsel set atomicnumber [lsearch $Molefacture::periodic [string toupper $element]]
      $newatomsel set name $element

      set name [::Molefacture::getAtomNewName $element $element [$mother get resid]]

      $newatomsel set name $name

      set newcoor [calc_geo $mindex $numbonds]

      set coorind  [::Molefacture::positionAvailable $newcoor $mindex [$newatomsel get index]]
      set movecoor [lindex $coorind 0]
      set avalableindex [lindex $coorind 1]

      $newatomsel moveto [::Molefacture::placeAtom $mindex $aindex $movecoor]

      vmd_addbond $mindex $aindex

      set motherbonds [atomselect $tmpmolid "index [join [$mother getbonds]] and occupancy 0.5"]

      if {$avalableindex < [llength $newcoor]} {
         set newcoor [lreplace $newcoor $avalableindex $avalableindex]
      }


      for {set i 0} {$i < [$motherbonds num]} {incr i} {
        if {[llength $newcoor] == 0} {break}
        set moveatom [atomselect $tmpmolid "index [lindex [$motherbonds get index] $i]"]
        set movecoor [lindex $newcoor 0]
        $moveatom moveto [::Molefacture::placeAtom $mindex [lindex [$motherbonds get index] $i] $movecoor]
        set newcoor [lreplace $newcoor 0 0]
      }
      $motherbonds delete
      #### reinforce the tetrahedral conformation
      set nbonds [llength [lindex [$mother getbonds] 0] ]
      if {$element == "H"} {
        set occupancy 0.5
      } else {
        set occupancy 0.8
      }
      $newatomsel set occupancy $occupancy
      if  {$nbonds > 0} {
        set bonds [lindex [$mother getbonds] 0]
        set selheavy [atomselect $tmpmolid "index [join $bonds] and occupancy >= 0.8"]
        ### Assign a 120 degree angle to the bonded heavy atoms (occupancy >= 0.8)
        if {[$selheavy num] == 2 && $occupancy >= 0.8} {
          set bondatms [$selheavy get index]
          set atm1 [lindex $bondatms 0]
          set atm2 $mindex
          set atm3 [lindex $bondatms 1]
          while {1} {
            set coords [$selheavy get {x y z}]
            # puts "DEBUG coordiantes order [$selheavy get index]"
            set curangle [measure angle [list $atm1 $atm2 $atm3] molid $tmpmolid]
            if {[expr abs($curangle)] >= 119.5 && [expr abs($curangle)] <= 120.5 } {
              break
            }
            set delta [expr (120.0 - $curangle)]
            # puts "DEBUG curangle $curangle|| delta $delta || angles [list [lindex $bonds 0] $mindex [lindex $bonds 1]]"
            #puts "DEBUG coord [lindex $coords 0] [lindex [$mother get {x y z}] 0] [lindex $coords 1]"
            set mat1 [trans angle [lindex $coords 0] [lindex [$mother get {x y z}] 0] [lindex $coords 1]  [expr (120.0 - $curangle)] deg]
            $newatomsel move $mat1
            # puts "DEBUG New Angle [measure angle [list $atm1 $atm2 $atm3] molid $tmpmolid] || atomindex [lindex $bonds 0] $mindex [lindex $bonds 1]"
          }
        }
        $selheavy delete
      }
      $newatomsel delete
      incr incrdum
   }
  }

  $dummies delete
  return $newatom

   # Restore the selection

}

######################################################################
### Place the atoms based on the vectors returned by the calc_geo  ###
### proc                                                           ###
######################################################################
proc ::Molefacture::placeAtom {mindex newatom bondvector} {
  variable tmpmolid

  set mother [atomselect $tmpmolid "index $mindex"]
  set natom [atomselect $tmpmolid "index $newatom"]
  set dist [::Molefacture::getBondDist [$mother get element] [$natom get element]]

  set newposition [vecadd [lindex [$mother get {x y z}] 0] [vecscale $bondvector $dist]]
  $mother delete
  $natom delete
  return $newposition
}
### Return a position available to move an atom to
proc ::Molefacture::positionAvailable {newcoor { mindex -1} {moveindex -1} } {
  variable tmpmolid
  set newcoorind 0
  set movecoor [Molefacture::placeAtom $mindex $moveindex [lindex $newcoor $newcoorind] ]

  set mother [atomselect $tmpmolid "index $mindex"]
  set natom [atomselect $tmpmolid "index $moveindex"]
  set mindist [expr [::Molefacture::getBondDist [$mother get element] [$natom get element]] * 0.5]
  $mother delete
  $natom delete
  ### Check if the new coordinates are already taken
  set selaux [atomselect $tmpmolid "sqrt( sqr(x - [lindex $movecoor 0]) + sqr( y - [lindex $movecoor 1]) + sqr( z - [lindex $movecoor 2]) ) < $mindist and ($Molefacture::notDumSel) and not index $mindex $moveindex"]

  set newlist $newcoor
  foreach coor $newcoor {
    lappend newlist [vecscale -1.0 $coor]
  }

  while {[$selaux num] != 0 && $newcoorind < [expr [llength $newlist] -1] } {
    $selaux delete
    incr newcoorind

    set movecoor [Molefacture::placeAtom $mindex $moveindex [lindex $newlist $newcoorind]]
    set selaux [atomselect $tmpmolid "(sqrt( sqr(x - [lindex $movecoor 0]) + sqr( y - [lindex $movecoor 1]) + sqr( z - [lindex $movecoor 2]) ) < $mindist and ($Molefacture::notDumSel)) and not index $mindex $moveindex"]

  }

  set coor [lindex $newlist $newcoorind]

  $selaux delete


  return [list $coor $newcoorind]
}

proc ::Molefacture::add_atom { element} {
   set newatoms [add_atom_noupdate $element]

   update_openvalence
  
   return $newatoms
}
#####################################################
# Add atoms from fragment file                      #
#####################################################
proc ::Molefacture::add_atom_from_frag {name type resname resid chain segname element x y z charge} {

  variable tmpmolid

   #Check to make sure we don't need more dummy hydrogens
   set dummies [atomselect $tmpmolid "$Molefacture::dumSel"]
   if {[expr [$dummies num] -[llength $name]] <= 0} {
     hrealloc "$Molefacture::notDumSel"
     $dummies delete
     set dummies [atomselect $tmpmolid "$Molefacture::dumSel"]
   }


  set hindex [lindex [$dummies list] 0]
  set newatom [atomselect $tmpmolid "index $hindex"]

  # Assign the proper characteristics
  set name [Molefacture::getAtomNewName $name $element $resid]

  $newatom set name $name
  $newatom set type $type
  $newatom set resname $resname
  $newatom set resid $resid
  $newatom set chain $chain
  $newatom set segname $segname
  $newatom set element $element
  $newatom set x $x
  $newatom set y $y
  $newatom set z $z
  $newatom set charge $charge

  #Set the occupancy to indicate its new status
  if {$element == "H"} {
    $newatom set occupancy 0.5
  } else {
    $newatom set occupancy 0.8
  }

  return $hindex
}

###########################################################
## Check if the atom's name already exist in the residue ##
## being edited                                          ##
###########################################################
proc ::Molefacture::getAtomNewName {name element resid } {
  variable tmpmolid
  set hselatms [atomselect $tmpmolid "($Molefacture::notDumSel) and element $element and resid $resid"]

  set atmnames [$hselatms get name]
  set ind 0
  set atmindex [lsearch -all $atmnames $name]
  # && [lsearch [lindex $atmresids $atmindex] $resid] > -1
  if {[llength $atmindex] > 1} {
    set name "${element}$ind"
    set atmindex [lsearch $atmnames $name]

    #&& [lsearch [lindex $atmresids $atmindex] $resid] > -1
    while {$atmindex != -1 } {
      set name "${element}$ind"
      if {[string length $name] > 4} {
        incr ind -1
        set name "${element}$ind"
        break
      }
      incr ind
      set atmindex [lsearch $atmnames $name]


    }
  }
  $hselatms delete
  return $name
}
proc ::Molefacture::read_fragmentlist {file} {
  #This proc will read a set of fragments in the molefacture fragment DB format
  #from a file, and then add them to the fragmentlist array

  variable fragmentlist
  variable periodic

  set infile [open $file r]

  while {[gets $infile line]>=0} {
#          puts $line
    if {[regexp {(^FRAGMENT\s*)(\w*)} $line match head fragname] > 0} {
      set myfrag [list]
      while {[gets $infile line]>=0 && [regexp {^END} $line]==0} {
        set i [regexp {\w+\s+\w+\s+(\w+)\s+\w+\s+\w+\s+\w+\s+(\-?\w+\.\w+)\s+(\-?\w+\.\w+)\s+(\-?\w+\.\w+)\s+\w+\.\w+\s+\w+\.\w+\s+\w+\s+(\w+)} $line match name x y z element]
        #set elementid [lsearch $periodic $element]
        set myatom [list $name $element $x $y $z]
        lappend myfrag $myatom
      }
      set fragmentlist($fragname) $myfrag
    }
  }
  close $infile
}

proc ::Molefacture::replace_atom_with_fragment {fragpath hydindex} {
  #This procedure will replace an atom
  #with a fragment

  variable tmpmolid
  # variable availabledummy

  set basesel [atomselect $tmpmolid "index $hydindex"]
  if {[llength [lindex [$basesel getbonds] 0] ]>1} {
      puts "You can only replace singly bonded atom!"
      return -1
  }
  Molefacture::saveState
  # Open the fragment file and read all the atoms to be added
  set fragfile [open $fragpath "r"]
  gets $fragfile line
  set names [list]
  set addx [list]
  set addy [list]
  set addz [list]
  set addelem [list]
  #First get the h-replacing atom
  gets $fragfile line
  set linearray [split $line]
  set linearray [noblanks $linearray]
  set repname [lindex $linearray 0]
  set repx [lindex $linearray 1]
  set repy [lindex $linearray 2]
  set repz [lindex $linearray 3]
  set repelem [lindex $linearray 4]

  #Now all of the others
  gets $fragfile line
  while {[string range $line 0 5] != "BONDS"} {
    set linearray [split $line]
    set linearray [noblanks $linearray]
    lappend names [lindex $linearray 0]
    lappend addx [lindex $linearray 1]
    lappend addy [lindex $linearray 2]
    lappend addz [lindex $linearray 3]
    lappend addelem [lindex $linearray 4]
    gets $fragfile line
  }

  #Finally, read the specified bonds
  set bondpart1 [list]
  set bondpart2 [list]
  set bondorders [list]
  while {[gets $fragfile line] != -1} {
    set linearray [split $line]
    set linearray [noblanks $linearray]
    lappend bondpart1 [lindex $linearray 0]
    lappend bondpart2 [lindex $linearray 1]
    lappend bondorders [lindex $linearray 2]
  }
  close $fragfile

   #Check to make sure we don't need more dummy hydrogens
   set dummies [atomselect $tmpmolid "$Molefacture::dumSel"]
   if {[expr [$dummies num] -[llength $names]] <= 0} {
     hrealloc "$Molefacture::notDumSel"
     $dummies delete
     set dummies [atomselect $tmpmolid "$Molefacture::dumSel"]
   }
  set basesel [atomselect $tmpmolid "index $hydindex"]
  $basesel set occupancy 0.8

  # Get residue/segment/etc information
  set myres [$basesel get resname]
  set myresid [$basesel get resid]
  set mchain [$basesel get chain]
  set mysegname [$basesel get segname]
  set mysegid [$basesel get segid]

  #Orient our coordinate system
  set mothers [$basesel getbonds]
  set mindex [lindex [lindex $mothers 0] 0]
  set msel [atomselect $tmpmolid "index $mindex"]
  set basecoor [lindex [$basesel get {x y z}] 0]
  set mcoor [lindex [$msel get {x y z}] 0]
  set myvec [vecsub $basecoor $mcoor]

  set tmpmolidaux $tmpmolid
  #Move the hydrogen to the proper place and replace it with a new atom
  $basesel moveto [list $repx 0 0]
  
  set name [::Molefacture::getAtomNewName $repname $repelem $myresid]
  $basesel set name $name

  $basesel set element $repelem
  $basesel set type $repname

  #List for processing bonds later
  set addedatoms [list]
  lappend addedatoms $hydindex


  #Now, the new piece is built near the origin, and then translated and rotated
  #to the proper position
  foreach name $names x $addx y $addy z $addz elem $addelem {
    set newid [add_atom_from_frag $name $name $myres $myresid $mchain $mysegname $elem $x $y $z 0.0]
    lappend addedatoms $newid
  }

  #### In case the addition of the new atoms implied the creation of a new molecule.
  if {$tmpmolidaux != $tmpmolid} {
    $msel delete
    set msel [atomselect $tmpmolid "index $mindex"]
  }
  # Adjust the order of the first bond, if desired
  if {[lindex $bondpart2 0] == -1} {
    set order [lindex $bondorders 0]
    set bondpart1 [lreplace $bondpart1 0 0]
    set bondpart2 [lreplace $bondpart2 0 0]
    set bondorders [lreplace $bondorders 0 0]
    set i 1
    while {$i < $order} {
      raise_bondorder [list [lindex [$msel get index] 0] [lindex [$basesel get index] 0]]
      incr i
    }
  }

  #Add the specified bonds
  foreach atom1 $bondpart1 atom2 $bondpart2 order $bondorders {
    vmd_addbond [lindex $addedatoms $atom1] [lindex $addedatoms $atom2] $order
  }

  #Properly orient and translate the new fragment
  set newsel [atomselect $tmpmolid "index $addedatoms"]
  $newsel move [transvec $myvec]
  $newsel moveby [list [$msel get x] [$msel get y] [$msel get z]]

  # Update molecule characteristics
  update_openvalence
  mol reanalyze $tmpmolid
  set newindices [$newsel get index]

  $newsel delete
  $msel delete
  $basesel delete
  set Molefacture::topocgenff 0
  
  return $newindices
  
}

proc ::Molefacture::new_independent_fragment {fragpath {dupselec 0}} {
  #This procedure creates a new molecule from a given base fragment
  #changes to previously edited molecules will be lost
  # The dupselec option indicates if the selected atoms are being duplicated as
  # a new fragment


  variable tmpmolid
  variable picklist

  if {$dupselec == 1 && [llength $picklist] == 0} {
    return
  }

  if {$tmpmolid == -1} {
   Molefacture::molefacture_gui_aux
  }

  set names [list]
  set addx [list]
  set addy [list]
  set addz [list]
  set addelem [list]
  set bondpart1 [list]
  set bondpart2 [list]
  set bondorders [list]
  # In the case of adding a fragment from file
  if {$dupselec == 0} {
    # Open the fragment file and read all the atoms to be added
    set fragfile [open $fragpath "r"]
    gets $fragfile line
    set linearray [split $line]
    set linearray [noblanks $linearray]
    set myresname [string range [lindex $linearray 1] 0 2]
    gets $fragfile line
    while {[string range $line 0 5] != "BONDS"} {
      set linearray [split $line]
      set linearray [noblanks $linearray]
      lappend names [lindex $linearray 0]
      lappend addx [lindex $linearray 1]
      lappend addy [lindex $linearray 2]
      lappend addz [lindex $linearray 3]
      lappend addelem [lindex $linearray 4]

      gets $fragfile line
    }

    #Finally, read the specified bonds

    while {[gets $fragfile line] != -1} {
      set linearray [split $line]
      set linearray [noblanks $linearray]
      lappend bondpart1 [lindex $linearray 0]
      lappend bondpart2 [lindex $linearray 1]
      lappend bondorders [lindex $linearray 2]
    }
    close $fragfile
  } else {
    set atomsel [atomselect $tmpmolid "index [join $picklist]"]
    set myresname [lindex [$atomsel get resname] 0]
    foreach name [$atomsel get name] coor [$atomsel get {x y z}] element [$atomsel get element] resid [$atomsel get resid] {
      set auxname "D$name"
      if {[string length $auxname] > 4} {
        set auxname [string range $auxname 0 3]
      }
      lappend names $auxname
      lappend addx [lindex $coor 0]
      lappend addy [lindex $coor 1]
      lappend addz [lindex $coor 2]
      lappend addelem $element
    }
    set atmindex [$atomsel get index]
    set atmbonds [$atomsel getbonds]
    set atmbondord [$atomsel getbondorders]
    set i 0

    foreach index $atmindex {
      set j 0
      foreach bond [lindex $atmbonds $i] {
        if {$index < $bond && [lsearch $atmindex $bond] > -1} {
          lappend bondpart1 [lsearch $atmindex $index]
          lappend bondpart2 [lsearch $atmindex $bond]
          lappend bondorders [lindex $atmbondord $i $j]
          incr newindex
        }
        incr j
      }
      incr i
    }
  }

  set dummies [atomselect $tmpmolid "$Molefacture::dumSel"]
  if {[expr [$dummies num] -[llength $names]] <= 0} {
     hrealloc "$Molefacture::notDumSel"
     $dummies delete
     set dummies [atomselect $tmpmolid "$Molefacture::dumSel"]
  }
  Molefacture::saveState
  ### get xmax (simple approach) to avoid overlapping the new frag with existing molecule
  set nondum [atomselect $tmpmolid "$Molefacture::notDumSel"]
  set xmax 0.0
  if {[$nondum num] > 0} {
    set xmax [expr [lindex [measure minmax $nondum] 1 0] + 1]
  }
  $nondum delete


  set addedatoms [list]


  set cursel [atomselect $tmpmolid "($Molefacture::notDumSel)"]
  set resid [lindex [lsort -unique -integer [$cursel get resid]] 0]

  if {$dupselec == 0} {
    set resid [expr $resid +1]
  }

  $cursel delete

  #Build the fragment itself
  foreach name $names x $addx y $addy z $addz elem $addelem {
    # puts "DEBUG: $name | $name | $myresname | 1 | A | 1 | $elem | $x | $y | $z"
    set newid [add_atom_from_frag $name $elem $myresname $resid "A" "A" $elem [expr $x + $xmax] $y $z 0.0]
    lappend addedatoms $newid
  }

  #Add the specified bonds
  foreach atom1 $bondpart1 atom2 $bondpart2 order $bondorders {
    # puts "DEBUG: vmd_addbond [lindex $addedatoms $atom1] [lindex $addedatoms $atom2] $order"
    vmd_addbond [lindex $addedatoms $atom1] [lindex $addedatoms $atom2] $order
  }

  mol reanalyze $tmpmolid

  # Update molecule characteristic
  update_openvalence

  #Recenter view
  display resetview
  if {[llength $addedatoms] == 0} {
    return -1
  }
  set Molefacture::topocgenff 0
  return $addedatoms

}

proc ::Molefacture::read_fragment_file {dbfile} {
  #Proc to read a fragment DB file
  #This currently clobbers the old fragment DB in memory
  #DB file format will be documented elsewhere

  variable addfrags
  array unset addfrags

  set fragfile [open $dbfile "r"]
  while {[gets $fragfile line] >= 0} {
    set myarr [split $line]
    set name [lindex $myarr 0]
    set fname [lindex $myarr 1]
    set myarr [lreplace $myarr 0 0 [regsub -all _ $name " "]]
    set myarr [lreplace $myarr 1 1 [file join $::env(MOLEFACTUREDIR) lib fragments $fname]]
    array set addfrags $myarr
  }
  close $fragfile

}

proc ::Molefacture::read_basefrag_file {dbfile} {
  #Proc to read a base fragment DB file
  #This currently clobbers the old fragment DB in memory
  #DB file format will be documented elsewhere

  variable basefrags
  array unset basefrags

  set fragfile [open $dbfile "r"]
  while {[gets $fragfile line] >= 0} {
    set myarr [split $line]
    set name [lindex $myarr 0]
    set fname [lindex $myarr 1]
    set myarr [lreplace $myarr 0 0 [regsub -all _ $name " "]]
    set myarr [lreplace $myarr 1 1 [file join $::env(MOLEFACTUREDIR) lib basemol $fname]]
    array set basefrags $myarr
  }
  close $fragfile

}

proc ::Molefacture::add_basefrag_file {dbfile} {
  #Proc to read a base fragment DB file
  #This currently clobbers the old fragment DB in memory
  #DB file format will be documented elsewhere

  variable basefrags

  set fragfile [open $dbfile "r"]
  while {[gets $fragfile line] >= 0} {
    set myarr [split $line]
    set name [lindex $myarr 0]
    set fname [lindex $myarr 1]
    set myarr [lreplace $myarr 0 0 [regsub -all _ $name " "]]
    set myarr [lreplace $myarr 1 1 [file join [file dirname $dbfile] $fname]]
    if {[array get $name] != ""} {
      array unset $name
    }
    array set basefrags $myarr
  }
  close $fragfile

}

proc ::Molefacture::add_frag_file {dbfile} {
  #Proc to read a base fragment DB file
  #DB file format will be documented elsewhere

  variable addfrags

  set fragfile [open $dbfile "r"]
  while {[gets $fragfile line] >= 0} {
    set myarr [split $line]
    set name [lindex $myarr 0]
    set fname [lindex $myarr 1]
    set myarr [lreplace $myarr 0 0 [regsub -all _ $name " "]]
    set myarr [lreplace $myarr 1 1 [file join [file dirname $dbfile] $fname]]
    if {[array get $name] != ""} {
      array unset $name
    }
    array set addfrags $myarr
  }
  close $fragfile
}

#####################################################
# "Grow" a nAtoms template, can be ring or chain    #
# depending on the Molefacture::cyclSelector        #
# equals to cycle or chain                          #
# The rings can be attached to one or two atoms.    #
# The chains can be attached to only one            #
# Useful when designing rings merged to each other  #
# Currently works for 3 or 9-atoms ring             #
#####################################################

proc Molefacture::growTemp { args } {
  variable tmpmolid
  variable picklist
  set index0 ""
  set index1 ""

  set nAtoms $Molefacture::cyclNAtoms
  if {$Molefacture::cyclSelector == "chain"} {
    set nAtoms $Molefacture::alchainNAtoms

    if {[llength $picklist] > 1 && [llength $args] == 0} {
      set picklistaux $picklist
      foreach atmsel $picklistaux {
        set picklist $atmsel
        Molefacture::growTemp $atmsel
      }
      return
    }
  }

  if {$nAtoms == "" || $nAtoms < 1} {
    return
  }

  if {[llength $picklist] == 0} {
    set newatoms [Molefacture::add_atom_gui C]
    set picklist $newatoms
    set index0 $newatoms
  } elseif {[llength $picklist] ==1} {
    set index0 $picklist
  } elseif {[llength $picklist] == 2} {
    set index0 [lindex $picklist 0]
    set index1 [lindex $picklist 1]
  } else {
    return
  }

  if {$Molefacture::cyclSelector == "cycle" && $nAtoms != 9 && $nAtoms != 8 && $nAtoms != 7 && $nAtoms != 6 && $nAtoms != 5 && $nAtoms != 4 && $nAtoms != 3} {
    return
  }

  ### check if the index0 is not empty and get the number
  ### of missing atoms = nAtoms - selected atoms
  ### just when $Molefacture::cyclSelector == cycle
  set remain $nAtoms
  set selIndex ""
  if {$index0 == ""} {
    return
  } else {
    incr remain -1
    lappend selIndex $index0
    if {$index1 != "" && $index1 != $index0} {
      incr remain -1
      lappend selIndex $index1
    }
  }

  Molefacture::saveState

  ### Internal angles for the different cycle groups
  set intangle 120
  if {$Molefacture::cyclSelector == "cycle"} {

    switch $nAtoms {
      9 {
        set intangle 140
      }
      8 {
        set intangle 135
      }
      7 {
        set intangle 129
      }
      6 {
        set intangle 120
      }
      5 {
        set intangle 108
      }
      4 {
        set intangle 90
      }
      3 {
        set intangle 60
      }
    }
  }

  ### grow the template
  set newatom ""
  set listatoms ""
  set prevtarget ""
  set bonds ""
  set mindex ""
  for {set i 0} {$i < $remain} {incr i} {
    set stillsel 1
    set mother ""
    if {$i < 2} {
      set lastnewatom ""
      if {$i == 0 || ($i == 1 && [llength $selIndex] == 2)} {
        set lastnewatom [lindex $selIndex $i]
      } elseif {$i == 1 && [llength selIndex] == 1}  {
        set lastnewatom $newatom
        set stillsel 0
      }
      set mother [atomselect $tmpmolid "index $lastnewatom"]
      set bonds  [join [$mother getbonds]]
      set mindex [$mother get index]
      if {[llength $bonds] > 0} {

        ### Delete hydrogen from the selected atoms
        set hydrogen [atomselect $tmpmolid "index $bonds and occupancy == 0.5 "]
        if {[$hydrogen num] > 0} {
          set picklist [lindex [$hydrogen get index] 0]
          del_atom $picklist
          update_openvalence
          $hydrogen delete
          
        }
      }
      $mother delete
      set picklist $lastnewatom
    } elseif {$newatom != ""} {
      set picklist $newatom

      set stillsel 0
    }

    ## calculate the position of the new atom

    set coors [Molefacture::calc_planar_geo $picklist]

    set redo 0
    if {[llength $coors] == 0} {
      set coors [Molefacture::calc_tetrahedral_geo $picklist]
      set redo 1
    }

    set newatom [Molefacture::add_atom C]
    
    if {$newatom < 0} {
      tk_messageBox -message "Can not add a skeleton to the selected atom." \
      -icon error -type ok -parent $Molefacture::topGui -title "Add Skeleton failed"
      return
    }
    if {[llength $newatom] == 0} {
      continue
    }
    update_openvalence

    set newatomsel [atomselect $tmpmolid "index $newatom"]
    if {$redo == 1} {
      $newatomsel moveto [join $coors]
    }

    ## add the new atom
    if {$i == 0 } {
      set listatoms $selIndex
    }
    lappend listatoms $newatom

    ### get the atoms bonded to the mother
    ### to evaluate angles and dihedral
    set mother [atomselect $tmpmolid "index $picklist"]
    set bonds  [join [$mother getbonds]]
    set mindex [$mother get index]

    set seltext "(index [join $bonds] and occupancy >= 0.8)"
    if {[llength $listatoms] >= 3} {
      append seltext " and index [join $listatoms]"
    }
    set selheavy [atomselect $tmpmolid $seltext]

    if {[$selheavy num] >= 2} {
      set bondatms [$selheavy get index]
      set atm1 [lindex $bondatms 0]
      set atm2 $mindex
      set atm3 [lindex $bondatms 1]
      set coords [$selheavy get {x y z}]
      set mcoor [lindex [$mother get {x y z}] 0]


      if {[llength $coors] == 0} {continue}



      set mcoor [lindex [$mother get {x y z}] 0]
      set coords [$selheavy get {x y z}]

      ### Calculate the current angle
      set curangle [measure angle [list $atm1 $atm2 $atm3] molid $tmpmolid]
      # puts "DEBUG curangle $curangle"

      ### Calculate the delta between the current the
      ### the target angle (intangle)
      set delta [expr ($intangle - $curangle)]
      if {$Molefacture::cyclSelector == "cycle" && $atm3 == $newatom} {
        set delta [expr $delta * -1]
      }



      ### Calculate the transformation to assign the angle
      set mat1 [trans angle [lindex $coords 0] $mcoor [lindex $coords 1]  $delta deg]

      ### move the atoms
      $newatomsel move $mat1


      ### assign the correct dihedral angle
      if {[llength $listatoms] >= 4 } {
        set coords [$selheavy get {x y z}]

        set atmlist ""
        set digbond ""
        set atmsel [atomselect $tmpmolid "index $atm1"]
        set antselbonds [lindex [$atmsel getbonds] 0]

        ### get the atoms bonded to mother that are part of the cycle group
        set antsel [atomselect $tmpmolid "(index [join $antselbonds] and occupancy >= 0.8) and index [join $listatoms]"]

        set antselbonds [$antsel get index]

        ### check if the atom1 has another atom bonded to it
        if {[llength $antselbonds] >= 2} {
          set index [lsearch $antselbonds $atm2]
          set antselbonds [lreplace $antselbonds $index $index]
          set index 0

          set atmlist [list [lindex $antselbonds 0] $atm1 $atm2 $atm3]
          set digbond [list [lindex $coords 0] $mcoor ]
        } else {
          set atmsel [atomselect $tmpmolid "index $atm3"]
          set antselbonds [lindex [$atmsel getbonds] 0]
          set antsel [atomselect $tmpmolid "(index [join $antselbonds] and occupancy >= 0.8) and index [join $listatoms]"]
          set antselbonds [$antsel get index]

          set index [lsearch $antselbonds $atm2]
          set antselbonds [lreplace $antselbonds $index $index]
          set index 0

          set atmlist [list $atm1 $atm2 $atm3 [lindex $antselbonds 0]]
          set digbond [list [lindex $coords 1] $mcoor ]

        }


        set curdih [measure dihed $atmlist molid $tmpmolid]
        set target 0.01
        if {$Molefacture::cyclSelector == "chain"} {
          if {$prevtarget != ""} {
            if {$prevtarget > 0} {
              set target -180.01
            } else {
              set target 180.01
            }

          } else {
            set target 180.01
          }

        }
        set delta [expr $target - $curdih]

        set mat1 [trans bond [lindex $digbond 0 ] [lindex $digbond 1] $delta deg]
        $newatomsel move $mat1

        set curdih [measure dihed $atmlist molid $tmpmolid]



        set prevtarget $target
      }

      ### close the ring
      if {$Molefacture::cyclSelector == "cycle" && $i == [expr $remain -1]} {
        set tobond ""
        set picklistaux ""
        if {[expr $nAtoms - $remain] == 2} {
          if {$nAtoms == 3} {
            set tobond [lindex $listatoms 1]
          } else {
            set tobond [lindex $listatoms 2]
          }
          set picklistaux [lrange $listatoms 0 1]
        } else {
          set tobond [lindex $listatoms 0]
          set picklistaux [lindex $listatoms 0]
        }
        set picklist [list $newatom $tobond]
        Molefacture::addDelBond add 0

        ### delete the last state which is stored when deleting the bond
        set picklist $picklistaux
      }

    }
    $newatomsel delete
    $selheavy delete
    $mother delete



  }
  if {[llength $listatoms] != 0} {
    topo -molid "$tmpmolid" -sel "index [join $listatoms]" guessatom radius element
    topo -molid "$tmpmolid" -sel "index [join $listatoms]" guessatom mass element
    set Molefacture::topocgenff 0
  }

  mol reanalyze $Molefacture::tmpmolid
  Molefacture::update_openvalence
  ::Molefacture::updateAtomTable
  Molefacture::clearTableSelection
  Molefacture::draw_selatoms

}

proc ::Molefacture::updateprotlasth {} {
  variable tmpmolid
  variable protlasthyd
  variable picklist


  ### No protein residue was added until now
  if {$protlasthyd == -2} {
    return
  }
  ### Update the hydrogen HC1 used during the protein building process
  ### this atom is not included in the CHARMM topology, so its' position is
  ### minimize
  if {$protlasthyd != -1} {
    set sel [atomselect $tmpmolid "name HC1"]
    if {[$sel num] != 0 && [$sel get resname] != "DUM"} {
      set mother [lindex [$sel getbonds] 0]

      set picklist $protlasthyd
      ::Molefacture::del_atom_gui
      set protlasthyd -1
      set picklist $mother
      Molefacture::set_prot_cterm
      Molefacture::updateAtomTable
      $sel delete
      return
    }

  }

  set sel [atomselect $tmpmolid "name HC1"]
  if {[$sel num] == 0 || [$sel get resname] == "DUM"} {
    set protlasthyd -1
    tk_messageBox -message "Unsigned growing C-terminal for the new segment. Please reassign the growing C-terminal to continue to add protein residues." -title "Unsigned growing C-terminal" -parent $Molefacture::topGui -icon warning -type ok
    return
  }

}

proc ::Molefacture::set_prot_cterm {} {
# Set the selected atom to be the target of future protein growth
  variable protlasthyd
  variable tmpmolid
  variable picklist
  if {[llength $picklist] != 1} {
    tk_messageBox -icon error -type ok -title "Error" -message "You need to select exactly one atom as the C-terminal end of the protein"
    return
  }

  set newparent [atomselect $tmpmolid "index $picklist"]
  if {[$newparent get element] != "C" && [$newparent get name] != "C" && [$newparent get type] != "C"} {
    tk_messageBox -icon warning -type ok -title "Select a carbon atom" -message "Please select the carbonyl Carbon atom to add the new residue."
    return
  }

  if {[llength $picklist] == 1} {
    set picklistold $picklist

    set bonds [$newparent getbonds]
    set bondSel [atomselect $tmpmolid "index [join $bonds] and element H"]
    set protlasthyd ""
    if {[$bondSel num] == 0} {
      if {[vecsum [join [$newparent getbondorders]] ]  == 4} {
        set protlasthyd -1
        return
      }
      set picklist [add_atom H]
      if {[llength $picklist] > 0} {
        topo -molid "$tmpmolid" -sel "index $picklist" guessatom radius element
        topo -molid "$tmpmolid" -sel "index $picklist" guessatom mass element
      }
      set sel [atomselect $tmpmolid "index $picklist"]
      $sel set name "HC1"
      $sel set type "HC1"
      $sel delete
      Molefacture::updateAtomTable
    } else {
      set picklist [lindex [$bondSel get index] 0]
    }
    $bondSel delete

  } else {
    return
  }
  $newparent delete



  set protlasthyd $picklist
}



proc ::Molefacture::add_aa {aaname} {
  variable aapath
  variable protlasthyd
  variable picklist
  variable tmpmolid

  variable topGui
  Molefacture::saveState

  if {$tmpmolid == -1} {
   Molefacture::molefacture_gui_aux
  }
  set Molefacture::mousemode "pick"
  if {[info exists Molefacture::deftopinfo] == -1 || $Molefacture::deftopinfo == ""} {
    ### pass the name of the amino-acid so the protein topology is added
    Molefacture::loadDefaulTopo $aaname
  }
  set newatoms ""

  ### Initiate the variable to sign that at least one protein residue was added
  set tmpmolidaux $tmpmolid


  if {$protlasthyd == -2} {
    set protlasthyd -1
  }

  if {$protlasthyd ==  -1 } {

    if {$aaname == "PRO"} {
       set newatoms [::Molefacture::new_independent_fragment [file join $aapath PRO.mfrag]]
    } else {
      set newatoms [::Molefacture::new_independent_fragment [file join $aapath GLY-end.mfrag]]
    }

    set newsel [atomselect $tmpmolid "($Molefacture::notDumSel)"]
    $newsel set resname $aaname
    if {$aaname != "PRO"} {
      ::Molefacture::set_newmol_phipsi
    }
  } elseif {$aaname == "PRO"} {
    set newatoms [::Molefacture::add_amino_acid_base [file join $aapath PRO.mfrag] $protlasthyd $aaname 0]
    set newsel [atomselect $tmpmolid "index $newatoms"]
  } else {
    set newatoms [::Molefacture::add_amino_acid_base [file join $aapath GLY.mfrag] $protlasthyd $aaname]
    set newsel [atomselect $tmpmolid "index $newatoms"]
  }
  # Set the new value of protlasthyd
  set atmindex [$newsel get index]
  set atmname [$newsel get name]

  set protlasthyd [lindex $atmindex [lsearch $atmname "HC1"]]

  if {$aaname != "PRO"} {

    set hareplace [lindex $atmindex [lsearch $atmname "HA1"]]
  }


  # Replace the appropriate hydrogen with the side chain
  if {$aaname != "GLY" && $aaname != "PRO"} {
    set haother [atomselect $tmpmolid "index [lindex $atmindex [lsearch $atmname HA2]]"]
    $haother set name "HA"
    $haother delete
    $newsel set resname $aaname

    lappend newatoms [replace_atom_with_fragment [file join $aapath "$aaname.mfrag"] $hareplace]
  }
  set picklist ""
  Molefacture::updateNewAtoms [join $newatoms]
}

proc ::Molefacture::add_nuc {nucname} {
# Add a nucleotide to a growing DNA or RNA chain
  variable nucpath
  variable nuclasthyd
  variable tmpmolid
  variable nuctype
  variable prevo3ind

  set nucbasefrag "dna_backbone"
  if {$nuctype == "RNA"} {
    set nucbasefrag "rna_backbone"
  }

  if {$nuclasthyd <  0 || $nuclasthyd == ""} {
    ::Molefacture::new_independent_fragment [file join $nucpath "$nucbasefrag.mfrag"]
    set newsel [atomselect $tmpmolid all]
    $newsel set resname $nucname
    set prevo3ind -1
  } else {
    # Get the index of the O3' atom, needed for dihedrals
    set lasthydsel [atomselect $tmpmolid "index $nuclasthyd"]
    set prevo3ind [join [join [$lasthydsel getbonds] ] ]
    $lasthydsel delete
    set newselind [::Molefacture::replace_atom_with_fragment [file join $nucpath "${nucbasefrag}_added.mfrag"] $nuclasthyd 1]
    set newsel [atomselect $tmpmolid "index $newselind"]
    $newsel set resname $nucname
  }

  # Set the new value of protlasthyd
  set nuclasthyd [lindex [$newsel get index] [lsearch [$newsel get name] "HM1"]]
  set hareplace [lindex [$newsel get index] [lsearch [$newsel get name] "HX"]]

  $newsel set resname $nucname

  # Replace the appropriate dummy atom with the base
  set baseindices [replace_atom_with_fragment [file join $nucpath "$nucname.mfrag"] $hareplace 1]
  $newsel delete

  mol reanalyze $tmpmolid

  # set dihedral angles
  ::Molefacture::set_nuc_diheds $baseindices
}

proc ::Molefacture::add_amino_acid_base {fragpath oldhydindex fragname {diheds 1}} {
  #This procedure will replace a hydrogen (which must have exactly 1 bond)
  #with a new amino acid. This should ONLY be used by the protein builder
  #This only builds the backbone; you need to add the side chain separately
  #Returns a selection containing the newly added atoms
  # If diheds is set to 0, the phi and psi angles are not adjusted

  variable tmpmolid

  # Set the proper phi/psi values
  # phi is between N and CA, psi between CA and C

  set oldHydSel [atomselect $tmpmolid "index $oldhydindex"]
  $oldHydSel set resid [expr [$oldHydSel get resid] +1]

  set newindex [::Molefacture::replace_atom_with_fragment $fragpath $oldhydindex]

  set newsel [atomselect $tmpmolid "index [join $newindex]"]

  set basesel [atomselect $tmpmolid "index $oldhydindex"]
  set baseSegName [$basesel get segname]
  set baseSegID [$basesel get segid]
  set baseResid [$basesel get resid]

  # $baseRes delete

  $newsel set segname $baseSegName
  $newsel set resid $baseResid
  $newsel set segid $baseSegID
  $newsel set resname $fragname


  set mothers [$basesel getbonds]
  set mindex [lindex [lindex $mothers 0] 0]
  set msel [atomselect $tmpmolid "index $mindex"]
  set mothersbonds [$msel getbonds]

  set mbsel [atomselect $tmpmolid "index [lindex $mothersbonds 0]"]
  set mcoind [lindex [$mbsel get index] [lsearch [$mbsel get type] "O"]]
  $mbsel delete

  variable phi
  variable psi
  set atmindex [$newsel get index]
  set atmname [$newsel get type]

  set nind [lindex $atmindex [lsearch -exact $atmname N]]
  set caind [lindex $atmindex [lsearch -exact $atmname CA]]
  set cind [lindex $atmindex [lsearch -exact $atmname C]]
  set oind [lindex $atmindex [lsearch -exact $atmname HC1]]
  set nhind [lindex $atmindex [lsearch -exact $atmname HN]]

#puts "inds: nind $nind caind $caind cind $cind oind $oind mcoind $mcoind oldhydindex $oldhydindex"
#puts "Source: [$newsel get name]"

  set prevind [lindex [$msel get index] 0]
  set casel [atomselect $tmpmolid "index $caind"]
  set csel [atomselect $tmpmolid "index $cind"]
  set nsel [atomselect $tmpmolid "index $nind"]
  set prevsel [atomselect $tmpmolid "index $prevind"]



  #Properly orient the backbone
  set ormovesel [atomselect $tmpmolid "(index [$newsel get index]) and not type N"]
  #puts "DEBUG  list $mcoind $tmpmolid list $prevind $tmpmolid list $nind $tmpmolid list $caind $tmpmolid"
  set dihedral [measure dihed [list [list $mcoind $tmpmolid] [list $prevind $tmpmolid] [list $nind $tmpmolid] [list $caind $tmpmolid] ]]

  set delta [expr $dihedral + 180]
  set mat [trans bond [lindex [$nsel get {x y z}] 0] [lindex [$prevsel get {x y z}] 0] $delta deg]
  $ormovesel move $mat
  set mat [trans bond [lindex [$nsel get {x y z}] 0] [lindex [$prevsel get {x y z}] 0] 180 deg]
  $ormovesel move $mat
  $ormovesel delete

  #Measure and set phi
  set phimovesel [atomselect $tmpmolid "(index [$newsel get index]) and not (type N HN CA)"]
  set dihedral [measure dihed [list [list $prevind $tmpmolid] [list $nind $tmpmolid] [list $caind $tmpmolid] [list $cind $tmpmolid] ]]

  set delta [expr -1 * ($phi - $dihedral)]
#puts "DEBUG: $dihedral $phi $delta"
#puts "DEBUG: [lindex [$casel get {x y z}] 0] [lindex [$basesel get {x y z}] 0]"
  set mat [trans bond [lindex [$casel get {x y z}] 0] [lindex [$basesel get {x y z}] 0] $delta deg]
  if {$diheds != 0} {$phimovesel move $mat}
#  puts "[lindex [lindex [label list Dihedrals] end] 4]"

  #Make sure we moved it right

  $phimovesel delete

  #Measure and set psi
  set psimovesel [atomselect $tmpmolid "(index [$newsel get index]) and not (type N HN CA HA1 HA2)"]
  set dihedral [measure dihed [list [list $nind $tmpmolid] [list $caind $tmpmolid] [list $cind $tmpmolid] [list $oind $tmpmolid] ]]

  set delta [expr -1 * ($psi - $dihedral)]
#puts "DEBUG: $dihedral $psi $delta"
  set mat [trans bond [lindex [$csel get {x y z}] 0] [lindex [$casel get {x y z}] 0] $delta deg]
  if {$diheds != 0} {$psimovesel move $mat}
#  puts "[lindex [lindex [label list Dihedrals] end] 4]"

  #Make sure we moved it right

  $psimovesel delete

  mol reanalyze $tmpmolid

  # Update molecule characteristics
  update_openvalence

  set newselind [$newsel get index]

  #Delete atom selections
  $nsel delete
  $prevsel delete
  $newsel delete
  $msel delete
  $basesel delete
  $casel delete
  $csel delete

  return $newselind
}

proc ::Molefacture::set_nuc_diheds {baseindices} {
  #Set dihedral angles of a newly added nucleotide
  variable tmpmolid
  variable ralpha
  variable rbeta
  variable rgamma
  variable rdelta
  variable repsilon
  variable rchi
  variable rzeta
  variable nuclasthyd
  variable prevo3ind

  set ressel [atomselect $tmpmolid "same residue as index $nuclasthyd"]

  set atomids [$ressel get index]
  set atomnames [$ressel get name]

  $ressel delete

  set phosind [lindex $atomids [lsearch -exact $atomnames P]]
  set o5ind [lindex $atomids [lsearch -exact $atomnames O5']]
  set c5ind [lindex $atomids [lsearch -exact $atomnames C5']]
  set c4ind [lindex $atomids [lsearch -exact $atomnames C4']]
  set c3ind [lindex $atomids [lsearch -exact $atomnames C3']]
  set o3ind [lindex $atomids [lsearch -exact $atomnames O3']]
  set h3ind [lindex $atomids [lsearch -exact $atomnames H3']]
  set c1ind [lindex $atomids [lsearch -exact $atomnames C1']]
  set c2ind [lindex $atomids [lsearch -exact $atomnames C2']]

  # Measure and move the alpha and zeta angles (which involve the previous residue)
  if {$prevo3ind >= 0} {
    set zetamovesel [atomselect $tmpmolid "(same residue as index $phosind)"]
    set prevPsel [atomselect $tmpmolid "name P and same residue as index $prevo3ind"]
    set prevPind [lindex [$prevPsel get index] 0]
    $prevPsel delete
    

    set alphamovesel [atomselect $tmpmolid "(same residue as index $phosind) and (not index $phosind)"]
    
    $alphamovesel delete
  }

  # Measure and move the beta angle
  set betamovesel [atomselect $tmpmolid "(same residue as index $phosind) and (not index $phosind $o5ind and not name O1P O2P O5T)"]
  #set_dihedral $betamovesel $rbeta $phosind $o5ind $c5ind $c4ind
  $betamovesel delete

  # Measure and move the gamma angle
  set gammamovesel [atomselect $tmpmolid "(same residue as index $phosind) and (not index $phosind $o5ind c5ind $c4ind and not name O1P O2P O5T H5' H5'')"]
  #set_dihedral $gammamovesel $rgamma $o5ind $c5ind $c4ind $c3ind
  $gammamovesel delete

  # Measure and move the delta angle
  # this one is a bit tricky due to the ring closure; we pretend the c4-o1 bond
  #  doesn't exist
  set deltamovesel [atomselect $tmpmolid "(same residue as index $phosind) and (not index $phosind $o5ind c5ind $c4ind $c3ind and not name O1P O2P O5T H5' H5'')"]
  
  $deltamovesel delete

  # measure and move the epsilon angle
  set epsilonmovesel [atomselect $tmpmolid "index $o3ind or (name HM1 and same residue as index $o3ind)"]
  
  $epsilonmovesel delete

  # Measure and move the chi angle, which includes the base
  set chimovesel [atomselect $tmpmolid "index $baseindices and not index [lindex $baseindices 0]"]
  
  $chimovesel delete

}

proc ::Molefacture::set_newmol_phipsi {} {
  #Set phi and psi angles of the first amino acid in a chain
  variable tmpmolid
  variable phi
  variable psi

  set allsel [atomselect $tmpmolid "occupancy >= 0.4"]

  set atmindex [$allsel get index]
  set atmname [$allsel get name]

  set prevind [lindex $atmindex [lsearch -exact $atmname HN]]
  set nind [lindex $atmindex [lsearch -exact $atmname N]]
  set caind [lindex $atmindex [lsearch -exact $atmname CA]]
  set cind [lindex $atmindex [lsearch -exact $atmname C]]
  set oind [lindex $atmindex [lsearch -exact $atmname HC1]]

  set basesel [atomselect $tmpmolid "index $nind"]
  set casel [atomselect $tmpmolid "index $caind"]
  set csel [atomselect $tmpmolid "index $cind"]

  #Measure and set phi
  set phimovesel [atomselect $tmpmolid "not (name N HN CA)"]
  set dihedral [measure dihed [list [list $prevind $tmpmolid] [list $nind $tmpmolid] [list $caind $tmpmolid] [list $cind $tmpmolid] ]]

  set delta [expr -1 * ($dihedral - $phi)]
  set mat [trans bond [lindex [$casel get {x y z}] 0] [lindex [$basesel get {x y z}] 0] $delta deg]
  $phimovesel move $mat
  $phimovesel delete

  #Measure and set psi
  set psimovesel [atomselect $tmpmolid "not (name N HN CA HA1 HA2)"]
  set dihedral [measure dihed [list [list $nind $tmpmolid] [list $caind $tmpmolid] [list $cind $tmpmolid] [list $oind $tmpmolid] ]]

  set delta [expr -1 * ($dihedral - $psi)]
  set mat [trans bond [lindex [$csel get {x y z}] 0] [lindex [$casel get {x y z}] 0] $delta deg]
  $psimovesel move $mat
  $psimovesel delete

  $csel delete
  $allsel delete
  $casel delete
  $basesel delete

}

proc ::Molefacture::add_proline {} {
  variable aapath
  variable protlasthyd
  variable tmpmolid

  set aaname "PRO"
  set newatoms ""
  if {$protlasthyd ==  -1} {
    set newatoms [::Molefacture::new_mol_from_fragment [file join $aapath PRO.mfrag]]
    set newsel [atomselect $tmpmolid all]
    $newsel set resname $aaname

  } else {
    set newatoms [::Molefacture::new_independent_fragment [file join $aapath PRO.mfrag] 0]
    set newsel [atomselect $tmpmolid "index $newatoms"]
  }

  # Set the new value of protlasthyd
  set protlasthyd [lindex [$newsel get index] [lsearch [$newsel get name] "HC1"]]

  $newsel delete

  Molefacture::updateNewAtoms [join $newatoms]

}

proc ::Molefacture::write_topfile {filename selection} {
# Write a topology file containing all the entries necessary for the current molecule
# For now, internal coordinates, DONOR cards, and the like are ignored

  variable improper_centeratoms

  variable tmpmolid

  variable mass_by_element

  unique_atomnames

  set sel [atomselect $tmpmolid "$selection"]

  set topfile [open $filename w+]
  set topcontents [list] ;# List for storing things that should be written at the end

  puts $topfile "*>>>>>> CHARMM topology file generated by Molefacture <<<<<<"
  puts $topfile "36 1"
  puts $topfile ""
  puts $topfile "read rtf card"
  puts $topfile ""

  set resnames [list] ;# Names of all residues encountered so far
  set atomnames [list] ;# All atom types needed
  set atomelems [list] ;# Elements of all atoms in atomnames
  set atommasses [list] ;# Masses of all atoms in atommasses

#Loop through residues and write the necessary entries
    set listTopAdded [list]
    foreach residue [lsort -unique [$sel get residue]] {

      # "($selection) and (occupancy >= 0.5) and (resid $resid)"
      set ressel [atomselect $tmpmolid "($selection) and (residue $residue)"]
      array set namearray {}
      set nbonds 1 ;# Number of bonds written so far

      # Get residue info
      set myresname [lsort -unique -dictionary [$ressel get resname]]
      set segname [lsort -unique [$ressel get segname]]
      set doBreak 0
      if {[lsearch $listTopAdded $myresname] > -1} {
        set doBreak 1
      } else {
        lappend listTopAdded $myresname
      }


      if {$doBreak == 1} {
        $ressel delete
        continue
      }

      set mycharge [vecsum [$ressel get charge]]
      # Make sure it isn't a duplicate
      if {[lsearch -exact $resnames $myresname] != -1} {
        puts "Warning: Found multiple instances of residue $myresname. Only the first will be written."
        continue
      }

      lappend resnames $myresname

      #Write the residue heading
      lappend topcontents [format "RESI %4s    % 5.2f\n" $myresname $mycharge]
      lappend topcontents "GROUP\n"

      # Loop through the atoms and write them all
      # In the process, record impropers that need to be written


      foreach index [$ressel get index] name [$ressel get name] type [$ressel get type] charge [$ressel get charge] element [$ressel get element] mass [$ressel get mass] {

      # If this has a hydrogen-like mass, it could be left over from
      # molecule building
      if {$mass < 2.0} {
#        puts "Looking for mass of element $element"
        set tmpmass [array get mass_by_element $element]
#        puts "got $tmpmass"
        if {$tmpmass != {}} {
          set mass [lindex $tmpmass 1]
        }
      }


      # Add this type to the types database if it isn't already there
        set currloc [lsearch -exact $atomnames $type]
        if {$currloc == -1} {
          # Then list the new atom
          lappend atomnames $type
          lappend atomelems $element
          lappend atommasses $mass
        } else {
          if {[string trim [lindex $atomelems $currloc]] != [string trim $element] || [expr [lindex $atommasses $currloc] - $mass] >= 0.001} {
             puts "Warning: Multiple incompatible definitions of type $type. Using the first"
          }
        }

        lappend topcontents [format "ATOM %4s %4s % .5f\n" $name $type $charge]
        array set namearray [list $index $name]
      }

      # Loop through again and this time write all the bonds we need
       lappend topcontents "BOND "
      foreach name [$ressel get name] index [$ressel get index] bonds [$ressel getbonds] bondorders [$ressel getbondorders] {
#        puts "DEBUG: bond $bonds | bo $bondorders"
        foreach bond $bonds bo $bondorders {
          if {$bo != 0 && $index < $bond} {
           #Make sure neither atom is missing
            if {$name == "" || [lindex [array get namearray $bond] 1] == ""} {
              puts "Warning: Bond found between two residues; it will be ignored"
              continue
            }
            if {[expr $nbonds % 4] == 0} {
       #       puts "DEBUG: nbonds $nbonds"
              lappend topcontents "\n"
              lappend topcontents "BOND "
            }

#            puts "DEBUG: [array get namearray $bond]"
#            puts "DEBUG: name $name | othername: [array get namearray $bond]"
           lappend topcontents "$name [lindex [array get namearray $bond] 1]  "
       #    puts "appended bond $name [lindex [array get namearray $bond] 1]"
           incr nbonds
          }
        }
      }

      # Go through one more time, and add impropers for all atoms that have appropriate types and exactly three bonds
      foreach name [$ressel get name] type [$ressel get type] index [$ressel get index] bonds [$ressel getbonds] {
        if { [lsearch $improper_centeratoms $type] >= 0 && [llength $bonds] == 3} {
          set bondednames [list]
          foreach bind $bonds {
            set tmpind [lsearch [$ressel get index] $bind]
            set othername [lindex [$ressel get name] $tmpind]
            lappend bondednames $othername
          }
          lappend topcontents "\nIMPR $name [lindex $bondednames 0] [lindex $bondednames 1] [lindex $bondednames 2]"
        }
      }


      lappend topcontents "\n\n"

     #Now we're done with the resudie
     $ressel delete
    }

  # Now, write all the MASS tags
  set i 1
  foreach type $atomnames elem $atomelems mass $atommasses {
    set ostring [format "MASS %5i %-4s  %8.5f %2s" $i $type $mass $elem]
    puts $topfile $ostring
  }

  puts $topfile ""
  puts $topfile "AUTO ANGLES DIHE"
  puts $topfile ""

# Finally, output all the residues
  foreach line $topcontents {
    puts -nonewline $topfile $line
  }

  puts $topfile "\nEND\n\n"

  set Molefacture::topocgenff 1
# Clean up
  close $topfile
  foreach res $resnames {
    if {[lsearch -index 6 $Molefacture::topoinfo $res] == -1 } {
      set handler [::Toporead::read_charmm_topology $filename 1]
      lappend Molefacture::topoinfo [::Toporead::topology_from_handler $handler]
      break
    } else {
      ### use the fact that in psfgen only the first definition of a duplicated residue is used
      set handler [::Toporead::read_charmm_topology $filename 1]
      set Molefacture::topoinfo [linsert $Molefacture::topoinfo 0 [::Toporead::topology_from_handler $handler]]
      break
    }
  }


  $sel delete
}

proc ::Molefacture::write_fep_hybrid_topfile {filename constantselout fullselin outgoingsel incomingsel totalcharge} {
# write a hybrid topology file for morphing between outgoingsel and incomingsel

  variable tmpmolid

  variable mass_by_element

  set topfile [open $filename w]
  set topcontents [list] ;# List for storing things that should be written at the end

  puts $topfile "*>>>>>> CHARMM topology file generated by Molefacture <<<<<<"
  puts $topfile "27 1"
  puts $topfile ""

  set resnames [list] ;# Names of all residues encountered so far
  set atomnames [list] ;# All atom types needed
  set atomelems [list] ;# Elements of all atoms in atomnames
  set atommasses [list] ;# Masses of all atoms in atommasses

#Loop through residues and write the necessary entries
  set myresname [lindex [$constantselout get resname] 0]
  array set namearray {}
  set nbonds 1 ;# Number of bonds written so far

  #Write the residue heading
  lappend topcontents [format "RESI %4s    % .5f\n" $myresname $totalcharge]
  lappend topcontents "GROUP\n"

  array set oldcomnamearray {}
  array set newcomnamearray {}

  foreach index [$constantselout get index] name [$constantselout get name] {
    array set oldcomnamearray [list $index $name]
  }

  foreach index [$fullselin get index] name [$fullselin get name] {
    array set newcomnamearray [list $index $name]
  }

  # Loop through the atoms and write them all
  set addedbonds [list]

  foreach ressel [list $constantselout $outgoingsel $incomingsel] selid [list 0 1 2] {
    array unset namearray
    foreach index [$ressel get index] name [$ressel get name] type [$ressel get type] charge [$ressel get charge] element [$ressel get element] mass [$ressel get mass] {

      # If this has a hydrogen-like mass, it could be left over from
      # molecule building
      if {$mass < 2.0} {
        #        puts "Looking for mass of element $element"
        set tmpmass [array get mass_by_element $element]
        #        puts "got $tmpmass"
        if {$tmpmass != {}} {
          set mass [lindex $tmpmass 1]
      }
    }

      # Add this type to the types database if it isn't already there
        set currloc [lsearch -exact $atomnames $type]
        if {$currloc == -1} {
          # Then list the new atom
          lappend atomnames $type
          lappend atomelems $element
          lappend atommasses $mass
        } else {
          #puts "Loc: $currloc"
          #puts [lindex atomelems $currloc]
#          if {[lindex atomelems $currloc] != $element || [lindex $atommasses $currloc] != $mass}
        }

        lappend topcontents [format "ATOM %4s %4s % .5f\n" $name $type $charge]
        array set namearray [list $index $name]
      }

      # Loop through again and this time write all the bonds we need
       #lappend topcontents "BOND "
       set nbonds 1
      foreach name [$ressel get name] index [$ressel get index] bonds [$ressel getbonds] bondorders [$ressel getbondorders] {
#        puts "DEBUG: bond $bonds | bo $bondorders"
        foreach bond $bonds bo $bondorders {
          if {$bo != 0} {
           #Make sure neither atom is missing
            set othername [lindex [array get namearray $bond] 1]
            if {$othername == ""} {
              if {$selid == 1} {
                set othername [lindex [array get oldcomnamearray $bond] 1]
              }
              if {$selid == 2} {
                set othername [lindex [array get newcomnamearray $bond] 1]
              }
            }
            if {$othername == ""} {
              continue
            }

           if {[string compare $name $othername] > 0} {
             set bondstring "$name $othername  "
           } else {
             set bondstring "$othername $name  "
           }
           if {[lsearch $addedbonds $bondstring] < 0} {
              if {[expr $nbonds % 4] == 1} {
                lappend topcontents "BOND "
              }
             lappend topcontents $bondstring
             lappend addedbonds $bondstring
              if {[expr $nbonds % 4] == 0} {
         #       puts "DEBUG: nbonds $nbonds"
                lappend topcontents "\n"
                #lappend topcontents "BOND "
              }
             incr nbonds
           }
          }
        }
      }
      lappend topcontents "\n"
      }
      lappend topcontents "\n"

  # Now, write all the MASS tags
  set i 1
  foreach type $atomnames elem $atomelems mass $atommasses {
    set ostring [format "MASS %5i %-4s  %8.5f %2s" $i $type $mass $elem]
    puts $topfile $ostring
  }

  puts $topfile ""
  puts $topfile "AUTO ANGLES DIHE"
  puts $topfile ""

# Finally, output all the residues
  foreach line $topcontents {
    puts -nonewline $topfile $line
  }

  puts $topfile "\nEND\n\n"

# Clean up
  close $topfile
}

proc ::Molefacture::set_dihedral {movesel newval ind1 ind2 ind3 ind4} {
# Acting on atoms in tmpmolid, move the selection movesel such that
#  the dihedral formed by the input atom indices is $newval
  variable tmpmolid

# Validate input
  if {$ind1 == "" || $ind2 == "" || $ind3 == "" || $ind4 == ""} {
    puts "WARNING: Couldn't find one or more atoms for set_dihedral"
    return
  }

  set dihedral [measure dihed [list [list $ind1 $tmpmolid] [list $ind2 $tmpmolid] [list $ind3 $tmpmolid] [list $ind4 $tmpmolid] ]]

  set bsel1 [atomselect $tmpmolid "index $ind2"]
  set bsel2 [atomselect $tmpmolid "index $ind3"]

  set delta [expr -1 * ($dihedral - $newval)]
  set mat [trans bond [lindex [$bsel1 get {x y z}] 0] [lindex [$bsel2 get {x y z}] 0] $delta deg]
  $movesel move $mat

  $bsel1 delete
  $bsel2 delete
}
