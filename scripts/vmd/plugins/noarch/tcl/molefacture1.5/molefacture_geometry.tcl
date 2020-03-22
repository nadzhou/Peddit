##
## $Id: molefacture_geometry.tcl,v 1.24 2019/06/25 20:49:23 jribeiro Exp $
##

proc ::Molefacture::invert_chir {mindex} {
# Invert the chirality of the selected atom
  variable tmpmolid

  set mother [atomselect $tmpmolid "index $mindex"]
  set numbonds [llength [lindex [$mother getbonds] 0]]
  if {$numbonds != 4} {return}

  set neighborlist [lsort -unique [join [$mother getbonds]]]

# Make atom selections for each neighbor, and their selection trees
  set neighborsels [list]
  set neighbortreesels [list]
  foreach ind $neighborlist {
    set mysel [atomselect $tmpmolid "index $ind"]
    lappend neighborsels $mysel
    set treesel [atomselect $tmpmolid "index [join [::util::bondedsel $tmpmolid $mindex $ind 50]] and not index $mindex"]
    lappend neighbortreesels $treesel
  }

  #puts "Treesels: $neighbortreesels"

# Find the two smallest fragments
  set min1 -1
  set min2 -1
  set minval1 999999
  set minval2 999999

  for {set i 0} {$i < 4} {incr i} {
    set mynum [[lindex $neighbortreesels $i] num]
    if {$mynum < $minval1} {
      set min2 $min1
      set minval2 $minval1
      set minval1 $mynum
      set min1 $i
    } elseif {$mynum < $minval2} {
      set min2 $i
      set minval2 $mynum
    }
  }


# Flip them

  # FIXME: We can have a better way of doing this if there's a ring, but it will be tricky
  if {[havecommonelems [lindex $neighbortreesels $min1] [lindex $neighbortreesels $min2] $mindex] > 0} {
    return
  }

#Calculate the necessary transforms and apply them
  set mcoor [join [$mother get {x y z}]]

  set sel1 [lindex $neighbortreesels $min1]
  set sel2 [lindex $neighbortreesels $min2]

  $sel1 moveby [vecscale -1 $mcoor]
  $sel2 moveby [vecscale -1 $mcoor]

  set coor1 [join [[lindex $neighborsels $min1] get {x y z}]]
  set coor2 [join [[lindex $neighborsels $min2] get {x y z}]]

  set vec1 $coor1
  set vec2 $coor2

  set inv1 [transvec $vec1]
  set inv2 [transvec $vec2]

  set x1 [transvecinv $vec1]
  set x2 [transvecinv $vec2]

  $sel1 move $x1
  $sel1 move $inv2

  $sel2 move $x2
  $sel2 move $inv1

  $sel1 moveby $mcoor
  $sel2 moveby $mcoor

  # Clean up
  $mother delete
  foreach sel $neighborsels {$sel delete}
  foreach sel $neighbortreesels {$sel delete}
}

proc ::Molefacture::set_planar {mindex} {
    variable tmpmolid
  #Accessory proc for set_planar; sets all atoms which are editable at a current
  #center to a planar geometry
    set mother [atomselect $tmpmolid "index $mindex"]
    set mcoor [join [$mother get {x y z}]]
    set numbonds [llength [lindex [$mother getbonds] 0]]
    if {$numbonds > 3} {puts "Error: There are too many bonds for a planar geometry" ; return 1}
    if {$numbonds < 3} {set geo [calc_linear_geo $mindex]}
    if {$numbonds == 3} {set geo [calc_planar_geo $mindex]}

    set oldhydrogens [atomselect $tmpmolid "index [join [$mother getbonds]] and occupancy 0.5"]

    set hindex [$oldhydrogens get index]
    for {set i 0} {$i < [$oldhydrogens num]} {incr i} {
      if {[llength $geo] == 0} {break}
     set moveatom [atomselect $tmpmolid "index [lindex $hindex $i]"]
     set movecoor [::Molefacture::placeAtom $mindex [lindex $hindex $i] [lindex $geo 0]]
     $moveatom moveto $movecoor
      set geo [lreplace $geo 0 0]
      $moveatom delete

    }
    $oldhydrogens delete
}

proc ::Molefacture::set_tetra {mindex} {
  variable tmpmolid
  set mother [atomselect $tmpmolid "index $mindex"]
  set mcoor [join [$mother get {x y z}]]
  set numbonds [llength [lindex [$mother getbonds] 0]]
  if {$numbonds > 4} {puts "Error: There are too many bonds for a tetrahedral geometry" ; return 1}
  set geo [calc_tetrahedral_geo $mindex 0]

  set oldhydrogens [atomselect $tmpmolid "index [join [$mother getbonds]] and occupancy 0.5"]

  set hindex [$oldhydrogens get index]
  for {set i 0} {$i < [$oldhydrogens num]} {incr i} {
     if {[llength $geo] == 0} {break}
     set moveatom [atomselect $tmpmolid "index [lindex $hindex $i]"]
     set movecoor [::Molefacture::placeAtom $mindex [lindex $hindex $i] [lindex $geo 0]]
     $moveatom moveto $movecoor

      $moveatom delete
     set geo [lreplace $geo 0 0]
    }

    $oldhydrogens delete
}

######################################################
## evaluate if has open valence to add another bond ##
######################################################
proc ::Molefacture::newbond {atmind} {
  variable tmpmolid
  variable periodic
  variable valence
  variable octet
  variable oxidation
  

  set sel [atomselect $tmpmolid "index $atmind"]

  set element  [lindex $periodic [$sel get atomicnumber]]
  set bonds [join [$sel getbonds]]
  set numbonds [llength $bonds]
#puts "DEBUG: 1"
#puts "Mother: $mindex"
#puts "Valence: [lindex [lindex [array get valence $element] 1] [lindex $oxidation $mindex]]"
  set curvalency [lindex [lindex [array get valence $element] 1] [lindex $oxidation $atmind]]
  set lonepairs [expr ([lindex [array get octet $element] 1]/2) - $curvalency]

  set numnewadds 0
  if {$numbonds > 0} {
    set newadds [atomselect $tmpmolid "index $bonds and occupancy 0.5"]
    set numnewadds [$newadds num]
    $newadds delete
  } else {
    set numnewadds 0
  }

  # set numnewcon $numadd
  # set newaddnum [$newadds num]
  # puts "numadd=$numnewadds lonepairs=$lonepairs numbonds=$numbonds newadd=$numnewadds openvalencelist=$openvalenceval"
  # set freebonds [expr $numadd + $numoldcon]
  set freebonds [expr $numbonds + $lonepairs + 1]
  $sel delete
  #### The openvalencelist is not used to evaluate the curvalency, so if there is a double or triple bond,
  #### it will no be reflected in this number, so use the openvalence insted. openvalencelist cannot be
  #### always used, because it is not updated in when adding the hydrogen, since update_openvalanece is too slow.
  #### if the proc was triggered during the addition of the atom ($Molefacture::mousemode == "addatom") and the valency is complete, always allow the addition of another atom. Otherwise, follow the open valency number
  set returnval $freebonds
  if {$Molefacture::openvalence($atmind) <= 0} {
    if {$Molefacture::mousemode != "addatom"} {
      set returnval [expr $Molefacture::openvalence($atmind) -1]
    }
  } elseif {[expr $numbonds +1] > $curvalency} {
    if {$Molefacture::mousemode != "addatom"} {
      set returnval [expr $curvalency - [expr $numbonds +1]]
    }
  }
  if {$returnval == 7} {
    set returnval -1
  }

  return $returnval
}

proc ::Molefacture::calc_geo { mindex numbond} {
  variable tmpmolid

  switch -- $numbond {
    0 {return [list]}
    1 {return [calc_linear_geo $mindex]}
    2 {return [calc_linear_geo $mindex]}
    3 {return [calc_planar_geo $mindex]}
    4 {return [calc_tetrahedral_geo $mindex]}
    5 {return [calc_trigbi_geo $mindex]}
    6 {return [calc_rhomb_geo $mindex]}
    default {puts "Unrecgonized geometry"}
  }
}

proc ::Molefacture::calc_e_geo { mindex numadd } {

  #Special function to calculate the geometry for placing
  #lone electrons and lone pairs
  variable tmpmolid
   variable periodic
   variable valence
   variable octet

  #First, figure out how many currently existing bonds there are
  set mother [atomselect $tmpmolid "index $mindex"]
  set mcoor [join [$mother get {x y z}]]
  set numbonds [llength [lindex [$mother getbonds] 0]]
  set element  [lindex $periodic [$mother get atomicnumber]]

  set numnewcon $numadd
  set numoldcon [expr $numbonds]

  switch -- [expr $numnewcon + $numoldcon] {
    0 {return [list]}
    1 {return [calc_linear_geo $mindex 1.0 true]}
    2 {return [calc_linear_geo $mindex 1.0 true]}
    3 {return [calc_planar_geo $mindex 1.0 true]}
    4 {return [calc_tetrahedral_geo $mindex 1.0 true]}
    5 {return [calc_trigbi_geo $mindex 1.0 true]}
    6 {return [calc_rhomb_geo $mindex 1.0 true]}
    default {puts "Unrecgonized geometry"}
  }
}

####################################################################
# Return positions where to put a new atom bonded to the a mother  #
# atom specified by $mindex.                                       #
####################################################################
proc ::Molefacture::calc_linear_geo { mindex {dolp false}} {
#  puts "Entering calc_linear_geo"
   variable tmpmolid

   # The mother atom to which the new atom shall be bonded
   set mother [atomselect $tmpmolid "index $mindex"]
   set mcoor [join [$mother get {x y z}]]

   # The antecedent atoms
   set aindex [join [$mother getbonds]]
   if {[llength [lindex $aindex 0]]==0} {
     set asel [atomselect $tmpmolid none]
   } elseif {$dolp} {
     set asel [atomselect $tmpmolid "index $aindex"]
   } else {
     set asel [atomselect $tmpmolid "index $aindex and occupancy >= 0.8"]
   }

   # If $mother is not bonded to other atoms the ligands are placed along the
   # x-axis
   if {[$asel num] == 0} {
      
      set acoor  [list -1.0 0 0]
      set axis {1 0 0}
      set madir [vecsub $acoor $mcoor]
      set t1 [transabout $axis 120 deg]
      set t2 [transabout $axis -120 deg]
      set hdir1 [coordtrans $t1 $madir]
      set hdir2 [coordtrans $t2 $madir]
      set hcoor1 [vecnorm $hdir1]
      set hcoor2 [vecnorm $hdir2]
      $asel delete
      return [list $hcoor1 $hcoor2]
   }

   # If only one bond exists, the first atom will be opposite that atom
   if {[$asel num]==1} {
      set acoor  [lindex [$asel get {x y z}] 0]
      set madir [vecsub $acoor $mcoor ]

      set hdir1 $madir
      set hcoor1 [vecnorm $hdir1]
      $asel delete
      return [list $hcoor1]
   }
    $asel delete
   return {}
}

proc ::Molefacture::calc_planar_geo { mindex {dolp false}} {
   variable tmpmolid

#   puts "Entering calc planar geo"
   # The mother atom to which the new atom shall be bonded
   set mother [atomselect $tmpmolid "index $mindex"]
   set mcoor [join [$mother get {x y z}]]

   # The antecedent atoms
   set aindex [join [$mother getbonds]]
   # If $mother is not bonded to other atoms the first ligand will be on a
   # position along the x-axis, the other two will be in the x-y plane
   if {[llength [lindex $aindex 0]]==0} {
     set asel [atomselect $tmpmolid none]
   } elseif {$dolp} {
     set asel [atomselect $tmpmolid "index $aindex"]
   } else {
     set asel [atomselect $tmpmolid "index $aindex and occupancy >= 0.8"]
   }
   if {[$asel num]==0} {
      
      set acoor  [list 1.0 0 0]
      set axis {0 0 1}
      set madir [vecsub $acoor $mcoor]
      set t1 [transabout $axis 120 deg]
      set t2 [transabout $axis -120 deg]
      set hdir1 [coordtrans $t1 $madir]
      set hdir2 [coordtrans $t2 $madir]
      set hcoor1 [vecnorm $hdir1]
      set hcoor2 [vecnorm $hdir2]
      $asel delete
      return [list $acoor $hcoor1 $hcoor2]
   } elseif {[$asel num]==1} {
    # If only one bond exists, the first atom will be syn to the first atom that's
   # bonded to the existing ligand. The second one will be in the position
   # that follows from planar geometry VSEPR
      set ante   $asel
      set acoor  [lindex [$ante get {x y z}] 0]
      # The antecedent of the first antecedent atom
      set aante  [atomselect $tmpmolid "index [join [$ante getbonds]] and occupancy 1"]
      set aacoor [lindex [$aante get {x y z}] 0]
      if {$aacoor == ""} {set aacoor [list 1 0 0]}
      # An axis perpendicular to aante ante and mother around which to
      # rotate the new atom to its final position.

      set axis [veccross [vecsub $acoor $aacoor] [vecsub $acoor $mcoor]]
      if {[veclength $axis] == 0.0} { set axis {0 0 1} }
      $aante delete

      set madir [vecsub $acoor $mcoor]
      set t1 [transabout $axis 120 deg]
      set t2 [transabout $axis -120 deg]
      set hdir1 [coordtrans $t1 $madir]
      set hdir2 [coordtrans $t2 $madir]
      set hcoor1 [vecnorm $hdir1]
      set hcoor2 [vecnorm $hdir2]
      $asel delete
      return [list $hcoor1 $hcoor2]
   } elseif {[$asel num]==2} {
      # In case there is only one position left, it can be obtained
      # by adding the existing bond vectors
      set hdir {0 0 0}
      set ante $asel
      foreach acoor [$ante get {x y z}] {
       set dir [vecsub $mcoor $acoor]
       set hdir [vecadd $hdir $dir]
      }

      $asel delete
      return [list [vecnorm $hdir] ]
   }
   $asel delete
   return {}
}

proc ::Molefacture::calc_tetrahedral_geo { mindex {dolp false}} {
   variable tmpmolid

#   puts "Entering calc tetrahedral geo"
   # The mother atom to which the new atom shall be bonded
   set mother [atomselect $tmpmolid "index $mindex"]
   set mcoor [join [$mother get {x y z}]]

   # The antecedent atoms
   set aindex [join [$mother getbonds]]
   if {[llength [lindex $aindex 0]]==0} {
     set asel [atomselect $tmpmolid none]
   } elseif {$dolp} {
     set asel [atomselect $tmpmolid "index $aindex"]
   } else {
     set asel [atomselect $tmpmolid "index $aindex and occupancy >= 0.8"]
   }
   # If $mother is not bonded to other atoms the first ligand will be on a
   # position along the x-axis, the other two will be in the x-y plane
   if {[$asel num]==0} {
      
      set acoor  [list 1.0 0 0]
      set axis {0 0 1}
      set madir [vecsub $acoor $mcoor]
      set t1 [transabout $axis 109 deg]
      set t2 [transabout $madir 120 deg]
      set t3 [transabout $madir -120 deg]
      set hdir1 [coordtrans $t1 $madir]
      set hdir2 [coordtrans $t2 $hdir1]
      set hdir3 [coordtrans $t3 $hdir1]

      if {[veclength $hdir1] == 0.0} {
        set hdir1 [list 0 -1 0]
      }
      if {[veclength $hdir2] == 0.0} {
        set hdir2 [list 0 0 1]
      }
      if {[veclength $hdir3] == 0.0} {
        set hdir3 [list 0 0 -1]
      }
      set hcoor1 [vecnorm $hdir1]
      set hcoor2 [vecnorm $hdir2]
      set hcoor3 [vecnorm $hdir3]
      $asel delete
      return [list $acoor $hcoor1 $hcoor2 $hcoor3]
   } elseif {[$asel num]==1} {
      set ante   $asel
      set acoor  [lindex [$ante get {x y z}] 0]

      # The antecedent of the first antecedent atom
      set aante  [atomselect $tmpmolid "index [join [$ante getbonds]] and occupancy 1"]
      set aacoor [lindex [$aante get {x y z}] 0]
      if {$aacoor == ""} {set aacoor [list 1 0 0]}

      # An axis perpendicular to aante, ante and mother around which to
      # rotate the new atom to its final position.
      set axis [veccross [vecsub $acoor $aacoor] [vecsub $acoor $mcoor]]
      
      if {[veclength $axis] < 0.5} {

        set i 0
        while {[veclength $axis] < 0.5} {
          if {$i == 0} {
            set axis [veccross $acoor $mcoor]

          } elseif {$i == 0} {
            set axis [veccross $acoor $aacoor]
          } else {
            set axis {0 0 1}
          }
         incr i
        }
      }
      set madir [vecsub $acoor $mcoor]
      set t1 [transabout $axis 109 deg]
      set t2 [transabout $madir 120 deg]
      set t3 [transabout $madir -120 deg]
      set hdir1 [coordtrans $t1 $madir]
      set hdir2 [coordtrans $t2 $hdir1]
      set hdir3 [coordtrans $t3 $hdir1]

      if {[expr abs([format %.2f [veclength $hdir1]])] < 0.5} {
        set hdir1 [list 0 -1 0]
      }
      if {[expr abs([format %.2f [veclength $hdir2]])] < 0.5} {
        set hdir2 [list 0 0 1]
      }
      if {[expr abs([format %.2f [veclength $hdir3]])] < 0.5} {
        set hdir3 [list 0 0 -1]
      }

      set hcoor1 [vecnorm $hdir1]
      set hcoor2 [vecnorm $hdir2]
      set hcoor3 [vecnorm $hdir3]
      $asel delete
      return [list $hcoor1 $hcoor2 $hcoor3]
   } elseif {[$asel num]==2} {
      set ante   $asel
      set acoor1  [lindex [$ante get {x y z}] 0]
      set acoor2  [lindex [$ante get {x y z}] 1]

      # An axis perpendicular to ante1 ante2 and mother around which to
      # rotate the new atom to its final position.
      set axis [veccross [vecsub $mcoor $acoor1 ] [vecsub  $mcoor $acoor2]]
      set axis2 [vecadd   [vecsub $mcoor $acoor1 ] [vecsub $mcoor $acoor2 ]]
      if {[veclength $axis] == 0.0} { set axis {0 0 1} }
      
      set madir [vecadd $acoor1 $acoor1]
      set t [transabout $axis 54.5 deg]
      set t1 [transabout $axis2 90 deg]
      set t2 [transabout $axis2 -90 deg]
      set hdir [coordtrans $t $axis2]

      set hdir1 [coordtrans $t1 $hdir]
      set hdir2 [coordtrans $t2 $hdir]
      if {[veclength $hdir1] == 0.0} {
        set hdir1 [list 0 -1 0]
      }
      if {[veclength $hdir2] == 0.0} {
        set hdir2 [list 0 0 1]
      }
      set hcoor1 [vecnorm $hdir1]
      set hcoor2 [vecnorm $hdir2]
      $asel delete
      return [list $hcoor1 $hcoor2]
   } elseif {[$asel num]==3} {
       set ante   $asel
       set bondcoor [$ante get {x y z}]
       set acoor1  [lindex $bondcoor 0]
       set acoor2  [lindex $bondcoor 1]
       set acoor3 [lindex $bondcoor 2]

       set av1 [vecsub $acoor1 $mcoor ]
       set av2 [vecsub $acoor2 $mcoor ]
       set av3 [vecsub $acoor3 $mcoor ]
      set veclist [list $av1 $av2 $av3]
       set hvec [vecadd $av1 $av2 $av3]
       ### if [veclength $hvec] is close to 0 the vecnorm of hvec will be basically
       ### useless.
       for {set i 0} {$i < 3} {incr i} {
        if {[expr abs([format %.2f [veclength $hvec]])] < 0.5} {
          set vecaux [lreplace $veclist $i $i]
          set hvec [eval "vecadd $vecaux"]
        } else {
          break
        }
       }

       set hcoor [vecnorm $hvec]
       $asel delete
       return [list $hcoor]
   }
    $asel delete
   return {}
}

proc ::Molefacture::calc_trigbi_geo { mindex  {dolp false}} {
   variable tmpmolid
   set axatoms [list]
   set eqatoms [list]
   set outcoor [list]
   set nax 0
   set neq 0

#   puts "Entering calc trigbi geo"
   # The mother atom to which the new atom shall be bonded
   set mother [atomselect $tmpmolid "index $mindex"]
   set mcoor [join [$mother get {x y z}]]

   # The antecedent atoms
   set aindex [join [$mother getbonds]]
   set asel ""
   if {[llength [lindex $aindex 0]]==0} {
     set asel [atomselect $tmpmolid none]
   } elseif {$dolp} {
     set asel [atomselect $tmpmolid "index $aindex"]
   } else {
     set asel [atomselect $tmpmolid "index $aindex and occupancy 0.8"]
   }

   #Number of connections already formed
   set numcon 0

   #Number of bonds the mother atom has currently
   set numbonds [$asel num]

   # Now, add the necessary atoms

   #First atom
   if {$numcon >= $numbonds} {
      # If there are no bonded atoms, assume that the first one is along the x axis

      set acoor  [list -1.0 0 0]
      set axis {0 0 1}
      set madir [vecsub $acoor $mcoor]
      set t1 [transabout $axis 180 deg]
      set hdir1 [coordtrans $t1 $madir]
      set acoor1 [vecnorm $hdir1]
      lappend outcoor $acoor1
    } else {
      set acoor1 [lindex [$asel get {x y z}] 0]
    }

    update idletasks
    incr numcon

    #Second atom
    if {$numcon >= $numbonds} {
      # If there is 1 bonded atom, set the second to be opposite to it

      set axis [vecsub $mcoor $acoor1]
      set axis [vecnorm $axis]
      set acoor2 $axis
      lappend outcoor $acoor2
      
    } else {
      set acoor2 [lindex [$asel get {x y z}] 1]
    }
    incr numcon

    update idletasks
    #Third atom
    if {$numcon >= $numbonds} {
      # First, measure the angle between the first 2 atoms
      set angle12 [measure_angle $acoor1 $mcoor $acoor2]
#puts "ANGLE12: $angle12"
      if {$angle12 < 100} {
        # if it's less than 100, assume they're a-e
        set axis [vecsub $mcoor $acoor1]
        set axis [vecnorm $axis]
        set acoor3 $axis
        set neq 1
        set nax 2
        lappend axatoms $acoor1
        lappend eqatoms $acoor2
        lappend axatoms $acoor3
      } elseif {$angle12 < 130} {
        # if it's less than 130, assume they're e-e
        set vec1 [vecsub $mcoor $acoor1]
        set vec2 [vecsub $mcoor $acoor2]
        set acoor3 [vecadd $vec1 $vec2]
        set acoor3 [vecnorm $acoor3]
        set acoor3 $mcoor
        set neq 3
        set nax 0
        lappend eqatoms $acoor1
        lappend eqatoms $acoor2
        lappend eqatoms $acoor3
      } else {
        # else, they're axial-axial
        set vec1 [vecsub $mcoor $acoor1]
        set vec2 [vecsub $mcoor $acoor2]
#puts "Vec1: $vec1 Vec2: $vec2"
        if {[veclength [vecsub $vec1 $vec2]] < 0.1 || [veclength [vecadd $vec1 $vec2]] < 0.1} { set vec2 [vecadd $vec2 {0.01 0.01 0.01}]}
#                puts "diff: [vecadd $vec1 $vec2]"
#                        puts "newv2: $vec2"
        set pvec [veccross $vec1 $vec2]
#puts "PVEC: $pvec"
        set pvec [vecnorm $pvec]
        set pvec $pvec
#puts "PVEC: $pvec"
        set acoor3 $pvec
#puts "Acoor3: $acoor3"
        set nax 2
        set neq 1
        lappend axatoms $acoor1
        lappend axatoms $acoor2
        lappend eqatoms $acoor3
      }
      set angle13 [measure_angle $acoor1 $mcoor $acoor3]
      lappend outcoor $acoor3
    } else {
      #Set up the axial/equatorial information
      set acoor3 [lindex [$asel get {x y z}] 2]
#puts "ac3: $acoor3"
      set angle13 [measure_angle $acoor1 $mcoor $acoor3]
      set angle12 [measure_angle $acoor1 $mcoor $acoor2]
      set angle23 [measure_angle $acoor2 $mcoor $acoor3]
      #puts "Done calculating angles"

      if {$angle13 < 100} {
        #1 and 3 are axial-equatorial
        if {$angle23 < 100} {
          #2 and 3 are axial equatorial
          if {$angle12 < 140} {
            #1 and 2 are equatorial
            set neq 2
            set nax 1
            lappend eqatoms $acoor1
            lappend eqatoms $acoor2
            lappend axatoms $acoor3
          } else {
            #3 is equatorial
            set nax 2
            set neq 1
            lappend axatoms $acoor1
            lappend axatoms $acoor2
            lappend eqatoms $acoor3
          }
        } elseif {$angle23 < 130} {
            #2 and 3 are both equatorial, so 1 is axial
            set nax 1
            set neq 2
            lappend axatoms $acoor1
            lappend eqatoms $acoor2
            lappend eqatoms $acoor3
        }
      } elseif {$angle13 < 130} {
        #1 and 3 are both equatorial
            if {$angle23 < 100} {
              #2 is axial
              set nax 1
              set neq 2
              lappend eqatoms $acoor1
              lappend axatoms $acoor2
              lappend eqatoms $acoor3
            } else {
              #2 is equatorial
              set neq 3
              set nax 0
              lappend eqatoms $acoor1
              lappend eqatoms $acoor2
              lappend eqatoms $acoor3
            }
      } else {
        #1 and 3 are both axial
        set nax 2
        set neq 1
        lappend axatoms $acoor1
        lappend eqatoms $acoor2
        lappend axatoms $acoor3
      }
      #puts "COUNTS: $axatoms $eqatoms"
  }

  incr numcon
   update idletasks
    #Fourth atom
    if {$numcon >= $numbonds} {
        #    puts "Nax: $nax Neq: $neq"
      #add an axial atom if one is missing
      if {$nax == 0} {
        set pvec [veccross $acoor1 $acoor2]
        set pvec [vecnorm $pvec]
        set pvec $pvec
        set acoor4 $pvec
        incr nax
        lappend axatoms $acoor4
      } elseif {$nax == 1} {
        set oppcoor [lindex $axatoms 0]
        set dir [vecsub $oppcoor $mcoor]
        set dir [vecnorm $dir]
        set acoor4 $dir
        incr nax
        lappend axatoms $acoor4
      } else {
        #Add an equatorial atom if necessary
        set ecoor1 [lindex $eqatoms 0]
        set ediff1 [vecsub $ecoor1 $mcoor]
#puts "Ecoor1: $ecoor1"
        set rotax [vecsub $mcoor [lindex $axatoms 0]]
        set t1 [transabout $rotax 120 deg]
        set acoor4 [coordtrans $t1 $ediff1]
#puts "Trans: $acoor4 || $rotax || $t1"
        set acoor4 [vecnorm $acoor4]
        set acoor4 $acoor4
#        set acoor4 [coordtrans $t1 $ecoor1]
        incr neq
        lappend eqatoms $acoor4
      }
      lappend outcoor $acoor4
    } else {
        set acoor4 [lindex [$asel get {x y z}] 3]
        set angle14 [measure_angle $acoor1 $mcoor $acoor4]
        set angle24 [measure_angle $acoor2 $mcoor $acoor4]
      if {$angle13 < 100} {
        if {$angle23 < 100} {
          incr neq
        } else {incr nax}
      } elseif {$angle13 < 130} {
          incr neq
      } else {
          incr nax
      }
    }

    incr numcon

    update idletasks
    #Fifth atom
    #puts "5th: $numcon $numbonds"
    #puts "Coordinates: $outcoor"
    if {$numcon >= $numbonds} {
      if {$nax == 1} {
        set oppcoor [lindex $axatoms 0]
        set dir [vecsub $mcoor $oppcoor]
        set dir [vecnorm $dir]
        set dir $dir
        set acoor5 $dir
      } else {
        #Add an equatorial atom if necessary
        set ecoor1 [lindex $eqatoms 0]
        set ecoor2 [lindex $eqatoms 1]
        set vec1 [vecsub $mcoor $ecoor1]
        set vec2 [vecsub $mcoor $ecoor2]
        set vec3 [vecadd $vec1 $vec2]
#puts "5: $vec3"
        set vec3 [vecnorm $vec3]
#puts "5: $vec3"
        set vec3 $vec3
        set acoor5 $vec3
      }
      #puts "Acoor5: $acoor5"
      lappend outcoor $acoor5
    }
    update idletasks
    #puts "Coordinates: $outcoor"
    return $outcoor

}

proc ::Molefacture::calc_rhomb_geo { mindex {dolp false}} {
   variable tmpmolid
   set axatoms [list]
   set eqatoms [list]
   set outcoor [list]
   set nax 0
   set neq 0

#   puts "Entering calc rhomb geo"
   # The mother atom to which the new atom shall be bonded
   set mother [atomselect $tmpmolid "index $mindex"]
   set mcoor [join [$mother get {x y z}]]

   # The antecedent atoms
   set aindex [join [$mother getbonds]]
   if {[llength [lindex $aindex 0]]==0} {
     set asel [atomselect $tmpmolid none]
   } elseif {$dolp} {
     set asel [atomselect $tmpmolid "index $aindex"]
   } else {
     set asel [atomselect $tmpmolid "index $aindex and occupancy 1"]
   }
   set acoor [$asel get {x y z}]

   #Number of bonds the mother atom has currently
   set numbonds [$asel num]

   #Variables for calculating the geometry
   set axis1 {1 0 0} ;#First axis of the bonded atoms
   set axis2 {0 1 0} ;#Second axis of the bonded atoms
   set axis3 {0 0 1} ;#Third axis of the bonded atoms
   set nax1 0 ;#Number of atoms along axis1
   set nax2 0 ;#Number of atoms along axis2
   set nax3 0 ;#Number of atoms along axis3
   set def_ax1 false ;#True if we have a non-default axis1
   set def_ax2 false ;#True if we have a non-default axis2
   set def_ax3 false ;#True if we have a non-default axis3

   # Now, calculate the axes we're going to use based on the antecedent atoms
   # If there are no current bonds, we use default axes and fall through
   if {$numbonds >= 1} {
      set axis1 [vecsub $mcoor [lindex $acoor 0]]
      set nax1 1
      set def_ax1 true
   }

   if {$numbonds >= 2} {
      set angle [measure_angle [lindex $acoor 0] $mcoor [lindex $acoor 1]]
      if {$angle < 120} {
        # Then they are not opposite each other
        set axis2 [vecsub $mcoor [lindex $acoor 1]]
        set nax2 1
        set def_ax2 true
      } else {
        incr nax1
      }
    }

    if {$numbonds >= 3} {
     if {$def_ax2 == false} {
        set axis2 [vecsub $mcoor [lindex $acoor 2]]
        set nax2 1
        set def_ax2 true
     } else {
        #See if atom 3 lies on axis 1 or axis 2
        set angle13 [measure_angle [vecadd $mcoor $axis1] $mcoor [lindex $acoor 2]]
        set angle23 [measure_angle [vecadd $mcoor $axis2] $mcoor [lindex $acoor 2]]
        if {$angle13 >120} {
          #Then atom 3 is on axis 1
          incr nax1
        } elseif {$angle23 >120} {
          #Then atom 3 is on axis 2
          incr nax2
        } else {
          #Atom 3 defines axis3
          set axis3 [vecsub $mcoor [lindex $acoor 2]]
          set nax3 1
          set def_ax3 true
        }
      }
    }

    if {$numbonds >= 4 && $def_ax3 == false} {
      #See if atom 4 lies on axis 1 or axis 2
      set angle14 [measure_angle [vecadd $mcoor $axis1] $mcoor [lindex $acoor 3]]
      set angle24 [measure_angle [vecadd $mcoor $axis2] $mcoor [lindex $acoor 3]]
      if {$angle14 >120} {
        #Then atom 4 is on axis 1
        incr nax1
      } elseif {$angle23 >120} {
        #Then atom 4 is on axis 2
        incr nax2
      } else {
        #Atom 4 defines axis3
        set axis3 [vecsub $mcoor [lindex $acoor 3]]
        set nax3 1
        set def_ax3 true
      }
    }

    if {$numbonds >= 5 && $def_ax3 == false} {
      #If we get to here, use atom 5 to define axis 3
      set axis3 [vecsub $mcoor [lindex $acoor 4]]
      set nax3 1
      set def_ax3 true
    }

    #Fill in missing axis definitions if needed
    if {$def_ax1 == true} {
      if {$def_ax2 == false} {set axis2 [veccross $mcoor $axis1]}
      if {$def_ax3 == false} {set axis3 [veccross $axis1 $axis2]}
    }

    #Normalize the axes
    set axis1 [vecnorm $axis1]
    set axis2 [vecnorm $axis2]
    set axis3 [vecnorm $axis3]

    #Construct a list of the new atoms needed

    while {$nax1 < 2} {
      set vec $axis1
      set newcoor $vec
      lappend outcoor $newcoor
      incr nax1
      set axis1 [vecscale $vec -1]
    }

    while {$nax2 < 2} {
      set vec $axis2
      set newcoor $vec
      lappend outcoor $newcoor
      incr nax2
      set axis2 [vecscale $vec -1]
    }

    while {$nax3 < 2} {
      set vec $axis3
      set newcoor $vec
      lappend outcoor $newcoor
      incr nax3
      set axis3 [vecscale $vec -1]
    }

    #puts "Rhomb output: $outcoor"
    return $outcoor

}

###########################################################
# Adjust the length of the selected bond.                 #
###########################################################

 proc ::Molefacture::adjust_bondlength { bond } {
    variable tmpmolid
    variable toCapAtoms
    variable picklist



    if {[llength $picklist] != 2} {
      return
    }
    set atom1 [lindex $picklist 0]
    set atom2 [lindex $picklist 1]
    set sel1 [atomselect $::Molefacture::tmpmolid "index $atom1"]
    set sel2 [atomselect $::Molefacture::tmpmolid "index $atom2"]
    set pos1 [join [$sel1 get {x y z}]]
    set pos2 [join [$sel2 get {x y z}]]
    $sel1 delete
    $sel2 delete

    set molsize [$Molefacture::atomTable size]
    set indexes1 [join [::util::bondedsel $tmpmolid $atom2 $atom1 -maxdepth $molsize]]
    set indexes2 [join [::util::bondedsel $tmpmolid $atom1 $atom2 -maxdepth $molsize]]
    #puts $indexes1
    #puts $indexes2
    set seltext1 "(index $indexes1) and not index $atom2 and not segname CAP"
    set seltext2 "(index $indexes2) and not index $atom1 and not segname CAP"


    set sel1 [atomselect $tmpmolid $seltext1]
    set sel2 [atomselect $tmpmolid $seltext2]

    set dir    [vecnorm [vecsub $pos1 $pos2]]
    set curval [veclength [vecsub $pos2 $pos1]]
    $sel1 moveby [vecscale [expr -0.5*($curval-$bond)] $dir]
    $sel2 moveby [vecscale [expr 0.5*($curval-$bond)] $dir]
    $sel1 delete
    $sel2 delete
    draw_selatoms
    draw_openvalence

    set Molefacture::bondguiValentry $bond
 }


###########################################################
# Rotate the selected bond.                               #
###########################################################

 proc ::Molefacture::rotate_bond { newdihed } {
    variable tmpmolid

    variable picklist
    variable atommarktags
    variable anglemarktags
    variable toCapAtoms

    set newdihed [format %.2f $newdihed]
    if {[llength $picklist] != 4 || $Molefacture::dihedrgrp == ""} { return }

    set angledef "$Molefacture::angleind1 $Molefacture::anglemother $Molefacture::angleind2"

    set delta [expr ($newdihed - $Molefacture::curangle)]
    set negdelta [expr $delta * -1]
    if {$Molefacture::dihedrgrp == "group1"} {
      set selText "index $Molefacture::anglemvgrp1list and not (segname CAP)"

      set sel1 [atomselect $tmpmolid $selText]
      set mat1 [trans bond $Molefacture::dihemother1coor $Molefacture::dihemother2coor $delta deg]
      $sel1 move $mat1
      $sel1 delete

    } elseif {$Molefacture::dihedrgrp =="group2"} {
      set selText "index $Molefacture::anglemvgrp2list and not (segname CAP)"

      set sel2 [atomselect $tmpmolid $selText]
      set mat2 [trans bond $Molefacture::dihemother1coor $Molefacture::dihemother2coor $negdelta deg]
      $sel2 move $mat2
      $sel2 delete

    }


    Molefacture::draw_angleatoms
    Molefacture::draw_selatoms
    set Molefacture::dihguiValentry $newdihed
    set Molefacture::curangle $newdihed

}


###########################################################
# Scale the selected angle.                               #
###########################################################

proc ::Molefacture::resize_angle { angle } {
    variable tmpmolid

    variable picklist
    variable atommarktags
    variable anglemarktags
    variable toCapAtoms
    set angle [format %.2f $angle]
    if {[llength $picklist]!=3 || ($Molefacture::anglemvgrp1 == 0 && $Molefacture::anglemvgrp2 == 0)} { return }
    set angledef "$Molefacture::angleind1 $Molefacture::anglemother $Molefacture::angleind2"

    set scale 1
    if {$Molefacture::anglemvgrp1 == 1 && $Molefacture::anglemvgrp2 == 1} {
      set scale 0.5
    }

    set delta [expr ($angle - $Molefacture::curangle) * $scale]
    set negdelta [expr $delta * -1]
    if {$Molefacture::anglemvgrp1 == 1} {
      set selText "index $Molefacture::anglemvgrp1list and not (segname CAP)"

       set sel1 [atomselect $tmpmolid $selText]
      set mat1 [trans angle $Molefacture::angleind1coor $Molefacture::anglemothercoor $Molefacture::angleind2coor $delta deg]
      $sel1 move $mat1
      $sel1 delete

    }
    if {$Molefacture::anglemvgrp2 == 1} {
      set selText "index $Molefacture::anglemvgrp2list and not (segname CAP)"

      set sel2 [atomselect $tmpmolid $selText]
      set mat2 [trans angle $Molefacture::angleind1coor $Molefacture::anglemothercoor $Molefacture::angleind2coor $negdelta deg]
      $sel2 move $mat2
      $sel2 delete

    }


    Molefacture::draw_angleatoms
    Molefacture::draw_selatoms
    set Molefacture::angleguiValentry $angle
    set Molefacture::curangle $angle
}

proc ::Molefacture::draw_openvalence_tr {dummy index op} {
  draw_openvalence
}

proc ::Molefacture::draw_openvalence {} {
   variable tmpmolid
   if {$tmpmolid<0} { return }

   # Delete any previously generated atomlabels
   variable ovalmarktags
   foreach tag $ovalmarktags {
      graphics $tmpmolid delete $tag
   }
   set ovalmarktags {}

   variable labelradius
   variable atomlist
   variable chargelist
   variable tmpmolid

   set cursel [atomselect $tmpmolid "occupancy>0.4"]
   # variable openvalencelist
   variable periodic
   variable octet
   variable valence
   variable showvale
   variable showellp
   variable oxidation

   if {$showvale == 1} {
   foreach index [$cursel get index] coords [$cursel get {x y z}] charge $chargelist {
#      puts "atom=$atom val=$val"
      if {$Molefacture::openvalence($index)== ""} { continue }
      if {$Molefacture::openvalence($index)} {
   lappend ovalmarktags [graphics $tmpmolid color yellow]
   #lappend ovalmarktags [graphics $tmpmolid sphere $coords radius [expr $labelradius*0.9]]
   lappend ovalmarktags [graphics $tmpmolid text $coords [format "  %+2i" $charge]]
      }
   }
  }

  #Draw lone pairs and lone electrons
  if {$showellp == 1} {
  graphics $tmpmolid color green
  foreach atom [$cursel get index] charge $chargelist {
     set mother [atomselect $tmpmolid "index $atom"]
     set element  [lindex $periodic [$mother get atomicnumber]]
#puts "element: $element"
#puts "octet: [array get octet $element]"
#puts "DEBUG: 3"
    set lonepairs [expr ([lindex [array get octet $element] 1]/2) - [lindex [lindex [array get valence $element] 1] [lindex $oxidation $atom] ]]
#puts "lp: $lonepairs"
#puts "DEBUG: 4"
    set nbonds 0
    foreach bo [join [$mother getbondorders]] {
      set bo [expr int($bo)]
      if {$bo < 0} { set bo 1 }
      incr nbonds $bo
    }
#puts "nbonds: $nbonds"
    set numadd [expr $lonepairs - $charge + $nbonds]; #[llength [lindex [$mother getbonds] 0]]
#puts $numadd
#    puts "For atom $atom I want to put on $lonepairs lone pairs, $charge lone electrons, and I think there are $numadd total things to add"
set geos [calc_e_geo [join [$mother get index]] [expr $lonepairs - $charge]]
#puts "geos: $geos"

#Check to see if we're working on a nitrogen near a conjugated system
# If so, make one of the lone pairs a barbell type
# If this is done, also make all the other electrons around it tetrahedral
    if {$lonepairs>0 && [$mother get element] == "N"} {
      #Look to see if there are double bonds on this atom or its immediate neighbors
      set dumbell false

      set bondlist [$mother getbonds]

      set bos [list]

      foreach bltemp $bondlist {
      foreach blindex $bltemp {
        set tempsel [atomselect $tmpmolid "index $blindex"]
        lappend bos [$tempsel getbondorders]
        $tempsel delete
      }
      }

      foreach bosdummy $bos {
      foreach bosdummy2 $bosdummy {
      foreach order $bosdummy2 {
        if {$order > 1.0} {set dumbell true}
      }
      }
      }

      if {$dumbell == true} {
        #Then make a barbell orbital, and decrement the number of lone pairs
        incr lonepairs -1
        set bond1 [lindex [lindex $bondlist 0] 0]
        set bond1sel [atomselect $tmpmolid "index $bond1"]
        set coor1 [lindex [$bond1sel get {x y z}] 0]
        $bond1sel delete
        if {[llength [lindex $bondlist 0]] > 1} {
          set bond2 [lindex [lindex $bondlist 0] 1]
#puts "DEBUG: [lindex $bondlist 0]"
#puts "DEBUG: $bond2"
          set bond2sel [atomselect $tmpmolid "index $bond2"]
          set coor2 [lindex [$bond2sel get {x y z}] 0]
          $bond2sel delete
        } else {
          set coor2 {0 0 0}
        }
        set mcoor [lindex [$mother get {x y z}] 0]
        set vec1 [vecsub $mcoor $coor1]
        set vec2 [vecsub $mcoor $coor2]
        set pvec [veccross $vec1 $vec2]
        if {[veclength $pvec] < 0.0001} {
          set pvec {1.0 0.0 0.0}
        } else {
          set pvec [vecnorm $pvec]
        }
        set coor1 [vecadd $mcoor $pvec]
        set coor2 [vecsub $mcoor $pvec]
        drawlpbarbell $mcoor $coor1 $coor2
        #Make the lone electrons planar
        set geos [calc_planar_geo [join [$mother get index]] 1.0 true]
      }
    }

    for {set i 0} {$i < $lonepairs} {incr i} {
      set mcoor [lindex [$mother get {x y z}] 0]
      set newcoor [lindex $geos end-$i]
      #puts "[llength $mcoor] [llength $newcoor]"
#puts "$mcoor | $newcoor"

#puts "DEBUG H"
      drawlp $mcoor $newcoor
#      puts "DEBUG I"
    }
#    puts "DEBUG J"

    for {set j 0} {$j > $charge} {incr i ; incr j -1} {
      set mcoor [lindex [$mother get {x y z}] 0]
      set newcoor [lindex $geos end-$i]
#puts "$i $newcoor"
      #puts "[llength $mcoor] [llength $newcoor]"

#puts "$mcoor | [$mother get index] | $newcoor"
      drawlone $mcoor $newcoor
    }

   #puts "open3"
  }
  }
  #puts "open4"
  $cursel delete

  # If fep mode is on, delete and update all of its floating spheres
  variable fepmolid
  variable atomarktags
  if {$fepmolid >= 0} {
    foreach graphid [graphics $fepmolid list] {
      graphics $fepmolid delete $graphid
    }

    foreach tag $atomarktags {
      set atomid [lindex $tag 0]
      set color [lindex $tag 1]
      set atomsel [atomselect $tmpmolid "index $atomid"]
      graphics $fepmolid color $color
      graphics $fepmolid sphere [join [$atomsel get {x y z}]] radius [expr $labelradius * 1.5] resolution 16
      $atomsel delete
    }

  }
}
