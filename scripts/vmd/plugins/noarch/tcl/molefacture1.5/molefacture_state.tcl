##
## $Id: molefacture_state.tcl,v 1.22 2019/06/25 20:49:23 jribeiro Exp $
##

proc ::Molefacture::init_oxidation {} {
  #Initialize the oxidation array

  variable oxidation
  variable tmpmolid
  variable periodic 
  variable valence
  set oxidation [list]
  set sel [atomselect $tmpmolid all]

  foreach atom [$sel get index] {
    set mother [atomselect $tmpmolid "index $atom"]
    set element  [lindex $periodic [$mother get atomicnumber]]
#puts "$atom $element [$mother get atomicnumber]"
    set valences [lindex [array get valence $element] 1] 
#puts "DEBUG a"
#    set valences [lindex [array get valence $element] [lindex $oxidation $atom]] 
    set i 0
#    puts "Valence: $valences"
    if {[llength $valences] > 1} {
#      puts "Checking valence for atom $atom"
      set currval [lindex $valences $i]
      set numbonds [llength [lindex [$mother getbonds] 0]]
#puts "Comparing: currval=$currval numbonds=$numbonds"
      while {$currval < $numbonds && $i < [llength $valences]} {
        incr i
        set currval [lindex $valences $i]
      }
#puts "Comparing: currval=$currval numbonds=$numbonds"
      if {[llength $valences] < $i} {incr i -1}

    }
    $mother delete
    lappend oxidation $i
#puts "DEBUG: index: $atom oxlength: [llength $oxidation]"
  }
#  puts $oxidation
#  puts "exit"
}




proc ::Molefacture::drawlone {startcoor endcoor} {

  if {$endcoor == ""} {return}
  variable tmpmolid
  graphics $tmpmolid color purple 
  variable ovalmarktags
  variable tmpmolid

  lappend ovalmarktags [graphics $tmpmolid cone $endcoor $startcoor radius 0.15 resolution 20] 
  lappend ovalmarktags [graphics $tmpmolid sphere $endcoor radius 0.15 resolution 20]
  graphics $tmpmolid color green
}

proc ::Molefacture::drawlp {startcoor endcoor} {

  if {$endcoor == ""} {return}
#puts "Entering drawlp"
  variable tmpmolid
  variable ovalmarktags
#puts "Coors: $startcoor $endcoor"
  set v1 [vecsub $startcoor {0 0 0}]
  set v2 [vecsub $endcoor {0 0 0}]
  set pvec [veccross $v1 $v2]
  if {[veclength $pvec] == 0} {set pvec {1 0 0}}
  set mvec [vecscale 0.3 [vecnorm $pvec]]

  set distvec [vecsub $startcoor $endcoor]
  set distvec [vecscale 0.75 $distvec]

  lappend ovalmarktags [graphics $tmpmolid cone [vecsub $startcoor $distvec] $startcoor radius 0.1 resolution 20]
  lappend ovalmarktags [graphics $tmpmolid sphere [vecadd $endcoor $mvec] radius 0.15 resolution 20]
  lappend ovalmarktags [graphics $tmpmolid sphere [vecsub $endcoor $mvec] radius 0.15 resolution 20] 

}

proc ::Molefacture::drawlpbarbell {startcoor endcoor1 endcoor2} {
  variable tmpmolid
  variable ovalmarktags

  lappend ovalmarktags [graphics $tmpmolid cone $endcoor1 $startcoor radius 0.15 resolution 20] 
  lappend ovalmarktags [graphics $tmpmolid sphere $endcoor1 radius 0.15 resolution 20]
  lappend ovalmarktags [graphics $tmpmolid cone $endcoor2 $startcoor radius 0.15] 
  lappend ovalmarktags [graphics $tmpmolid sphere $endcoor2 radius 0.15 resolution 20 ]
}

proc ::Molefacture::assign_elements {} {
   variable tmpmolid
   variable periodic

   set sel [atomselect $tmpmolid "occupancy>0.4"]
   topo -molid $tmpmolid -sel "occupancy>0.4" guessatom element mass

   $sel delete
   mol reanalyze $tmpmolid
}

proc ::Molefacture::update_openvalence { } {
   variable tmpmolid
   variable periodic
   variable valence
   variable atomlist        [list] 
   
   array unset Molefacture::openvalence *
   array set ::Molefacture::openvalence ""
   variable totalcharge 0
   variable oxidation
   variable atomlistformat
   variable projectsaved
   set projectsaved 0

   set sel [atomselect $tmpmolid "occupancy>0.4"]
   set n 0
   
   set indexlist [$sel list]
   foreach vmdind $indexlist bondorders [$sel getbondorders] elemind [$sel get atomicnumber] {
      if {$vmdind == ""} { continue }
      set nbonds 0.0
      foreach bo $bondorders {
         if {$bo < 0} { set bo 1 } 
         set nbonds [expr $nbonds + $bo]
      }
      set nbonds [expr int($nbonds)]
      set element [lindex $periodic $elemind]
      if {[lindex $oxidation $vmdind] == ""} {set oxidation [lset oxidation $vmdind 0]}
      set val  [lindex [lindex [array get valence $element] 1] [lindex $oxidation $vmdind]]
      set openval [expr $val-$nbonds]

      set Molefacture::openvalence($vmdind) $openval
      set indarray($vmdind) [lindex $indexlist $n]
      incr n
   }
   
   
   foreach sbonds [$sel getbonds] elemind [$sel get atomicnumber] sindex [$sel get index] bondorder [$sel getbondorders] {
      set element [lindex $periodic $elemind]
      ## To fix display of X--(PO4-)--X
      if {$element == "P" && $Molefacture::openvalence($sindex) > 0} {
        set incrbonds 0
        set nearbyopenind [list]
        foreach b $sbonds {
            set bondedatomvalence $Molefacture::openvalence($indarray($b))
            if { $bondedatomvalence > 0 } {
                lappend nearbyopenind $indarray($b)
            }
        }
        foreach n $nearbyopenind {
            if { $Molefacture::openvalence($n) > 0 && $openval > 0 } {
                incr incrbonds
                set Molefacture::openvalence($n) [expr $Molefacture::openvalence($n)-1]
                incr openval -1
            }
        }
        
        set Molefacture::openvalence($sindex)
      } 

      if {$element == "P"} {
        ### Remove the double double bonds in phosphorus atoms
        ### Set by IDATM
        set numOfDoubles [lsearch -all [join $bondorder] "2.0"]

        if {[llength $numOfDoubles] > 1} {
          set pSel [atomselect $tmpmolid "index $sindex"]
          lset bondorder [lindex $numOfDoubles end] "1.0"
          $pSel setbondorders [list $bondorder]
          $pSel delete
        }

      }
   }
 

   # Clean up
   catch {$sel delete}

   # Delete unneeded temporary molecules
   foreach molid [molinfo list] {
     #puts "Checking for deletion: $molid"
     if {$molid == $tmpmolid} {continue}
     if {[string compare -length 12 "Molefacture_" [molinfo $molid get name]] == 0} {mol delete $molid}
   }

   draw_openvalence

}

proc ::Molefacture::fix_changes {} {
  variable tmpmolid
  set sel [atomselect $tmpmolid "occupancy>0.45 and occupancy<1"]
  $sel set occupancy 0.8
  $sel delete

  set sel [atomselect $tmpmolid "occupancy < 0.45"]
  $sel set occupancy 0
  $sel delete

  mol reanalyze $tmpmolid
}

proc ::Molefacture::undo_changes {} {
  # Remove added hydrogens
  variable tmpmolid
  set sel [atomselect $tmpmolid "occupancy 0.5"]
  set indexarray [$sel get index]
  foreach i $indexarray {
    set sel2 [atomselect $tmpmolid "index $i"]
    set bonds2 [$sel2 getbonds]
    foreach j $bonds2 {
      vmd_delbond $i $j 
      vmd_delbond $j $i
    }
    $sel2 delete
  }
  $sel set occupancy 0
    
  $sel delete

  # Undelete deleted atoms
  set sel [atomselect $tmpmolid "occupancy 0.3"]
  set indexarray [$sel get index]
  foreach i $indexarray {
    set sel2 [atomselect $tmpmolid "index $i"]
    set bonds2 [$sel2 getbonds]

    foreach j [join $bonds2] {

      vmd_addbond $i $j 
      vmd_addbond $j $i
    }
    $sel2 delete
  }
  $sel set occupancy 1
  $sel delete
  update_openvalence
}

# proc ::Molefacture::update_openvalence_FEP { } {
#    variable tmpmolid
#    variable periodic
#    variable valence
#    variable FEPatomlist     [list]
#    variable chargelist      [list]
#    variable openvalencelist [list]
#    variable totalcharge 0
#    variable oxidation
#    variable FEPlistformat
#    variable projectsaved
#    set projectsaved 0

# #puts "DEBUG D"
# #puts "$tmpmolid"
#    set sel [atomselect $tmpmolid "occupancy>0.4"]
#    foreach vmdind [$sel list] name [$sel get name] sbonds [$sel getbonds] bondorders [$sel getbondorders] elemind [$sel get atomicnumber] truecharge [$sel get charge] type [$sel get type] beta [$sel get beta] {
#       if {$vmdind == ""} { continue }
#       set nbonds 0
#       foreach bo $bondorders {
#          set bo [expr int($bo)]
#          if {$bo < 0} { set bo 1 } 
#          incr nbonds $bo
#       }
# #      puts "Elemind: $elemind"
#       set element [lindex $periodic $elemind]
# #puts "DEBUG: 5"
#       if {[lindex $oxidation $vmdind] == ""} {set oxidation [lset oxidation $vmdind 0]}
# #    puts "Ox: $vmdind | [lindex $oxidation $vmdind] | [lindex [array get valence $element] 1]"
#       set val  [lindex [lindex [array get valence $element] 1] [lindex $oxidation $vmdind]]
# #set val 4
# #puts "Index: $vmdind Valence: $val ox: [lindex $oxidation $vmdind]"
# #puts "$vmdind | $name | $sbonds | $bondorders | $elemind | $truecharge | $type"
#       set openval [expr $val-$nbonds]
# #puts "DEBUG: 6"
# #puts "For atom of type $element, we have valence of $val, $nbonds bonds, and an open valence of $openval"
#       lappend openvalencelist $openval
#       set charge [expr $nbonds-$val]
#       set totalcharge [expr $charge + $totalcharge]
#       lappend chargelist $charge
#       set myval [lindex [lindex [array get valence $element] 1] [lindex $oxidation $vmdind]]
# #puts "Index: $vmdind oxstate $myval oxnum [lindex $oxidation $vmdind]"
# #      puts "\[format \"$atomlistformat\" $vmdind $name $element $openval $charge $myval\]"
#       if {$type == ""} {
#         set type "X"
#       }
#       set myval [expr int($myval)]
#       if {$element == ""} {set element "X"}
#       if {$type    == ""} {set type   "{}"}
#       set formatcommand "\[format \"$FEPlistformat\" $vmdind \"$name\" \"$type\" \"$element\" $beta\]"
# #      puts "Element: $element"
# #puts "$formatcommand"
#       eval set atomelem "$formatcommand"
#       #lappend atomlist [format "%5i %4s %2s %2i %+2i" $vmdind $name $element $openval $charge]
#       lappend FEPatomlist $atomelem
#    }

#    # Clean up
#    catch {$sel delete}

#    draw_openvalence
#    variable bondlist [bondlist]
#    variable anglelist [anglelist]
# }


##########################################################
# Save all the necessary information to undo processes   #
#                                                        #
##########################################################
proc Molefacture::saveState {} {
  variable tmpmolid
  variable undo
  variable openvalence
  variable firstUndo 
  variable oxidation

  if {$tmpmolid == -1} {
    return
  }
  set newState [list]

  set sel [atomselect $tmpmolid "$Molefacture::notDumSel or segname CAP"]

  # 1st element atoms' index
  lappend newState [$sel get index]

  # 2nd element atoms' name
  lappend newState [$sel get name]

  # 3rd element atoms' type
  lappend newState [$sel get type]

  # 4th element atoms' element
  lappend newState [$sel get element]

  # 5th element atoms' charge
  lappend newState [$sel get charge]

  # 6th element atoms' resname
  lappend newState [$sel get resname]

  # 7th element atoms' resid
  lappend newState [$sel get resid]

  # 8th element atoms' chain
  lappend newState [$sel get chain]

  # 9th element atoms' segname
  lappend newState [$sel get segname]

  # 10th element atoms' coordinates
  lappend newState [$sel get {x y z}]

  # 11th element bonds
  lappend newState [$sel getbonds]

  # 12th element bondorders
  lappend newState [$sel getbondorders]

   # 13th element occupnacy
  lappend newState [$sel get occupancy]

  # 14th element openvalence
  lappend newState [array get openvalence]

  # 15th element oxidation
  lappend newState $oxidation

  $sel delete

  set length [llength $undo]
  if {$length == 25} {
    set undo [lrange $undo 1 end]
  }
  lappend undo $newState
  set firstUndo 1
  return
}


##########################################################
# Undo last action by revert to the last state saved in  #
# undo variable                                          #
##########################################################
proc Molefacture::undoProc {} {
  variable tmpmolid
  variable undo
  variable openvalence
  variable firstUndo
  variable oxidation

  if {[llength $undo] == 0} {
    return
  }

  set lastIndex [expr [llength $undo] -1]
  
  set lastState [lindex $undo $lastIndex]

  set sel [atomselect $tmpmolid "all"]

  $sel set resname "DUM"
  $sel set occupancy 0
  set emptylist [string repeat "{} " [$sel num]]
  $sel setbonds $emptylist
  $sel setbondorders $emptylist
  $sel delete

  set sel [atomselect $tmpmolid "index [lindex $lastState 0]"]
  # 1st element atoms' index
  # $sel set index [lindex $lastIndex 0]

  # 2nd element atoms' name
  $sel set name [lindex $lastState 1]

  # 3rd element atoms' type
  $sel set type [lindex $lastState 2]

  # 4th element atoms' element
  $sel set element [lindex $lastState 3]

  # 5th element atoms' charge
  $sel set charge [lindex $lastState 4]

  # 6th element atoms' resname
  $sel set resname [lindex $lastState 5]

  # 7th element atoms' resid
  $sel set resid [lindex $lastState 6]

  # 8th element atoms' chain
  $sel set chain [lindex $lastState 7]

  # 9th element atoms' segname
  $sel set segname [lindex $lastState 8]

  # 10th element atoms' coordinates
  $sel set {x y z} [lindex $lastState 9]

  # 11th element bonds
  $sel setbonds [lindex $lastState 10]

  # 12th element bondorders
  $sel setbondorders [lindex $lastState 11]

  # 13th element occupnacy
  $sel set occupancy [lindex $lastState 12]
  

  # 14th element openvalence

  array unset Molefacture::openvalence *
  array set ::Molefacture::openvalence [lindex $lastState 13]

  # 15th element oxidation
  set oxidation [lindex $lastState 14]

  $sel delete
  mol reanalyze $tmpmolid
  
  if {[llength $undo] > 1} {
    set undo [lreplace $undo $lastIndex $lastIndex]
  }
  
  ::Molefacture::updateAtomTable
  update_openvalence
  ::Molefacture::clearTableSelection
  return
}