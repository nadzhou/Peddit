##
## $Id: molefacture_edit.tcl,v 1.16 2019/06/25 20:49:23 jribeiro Exp $
##



proc ::Molefacture::set_bondorder { {bondorder 1} } {
   variable tmpmolid

  variable picklist
   if {[Molefacture::isbond] == 0} {
    return
   }
   Molefacture::saveState
   set all [atomselect $tmpmolid "all"]
   set bondorders [$all getbondorders]
   
   set sb0 [lindex [$all getbonds] [lindex $picklist 0]]
   set sb1 [lindex [$all getbonds] [lindex $picklist 1]]
   set vmdind0 [lindex $picklist 0]
   set vmdind1 [lindex $picklist 1]

   set pos0in1 [lsearch $sb1 $vmdind0]

   set pos1in0 [lsearch $sb0 $vmdind1]
   
   lset bondorders $vmdind0 $pos1in0 $bondorder
   lset bondorders $vmdind1 $pos0in1 $bondorder

    set delatom [list]
    foreach atom $picklist {
      set numbonds [Molefacture::newbond $atom]

       if { $numbonds <= -1} {
        set selatom [atomselect $tmpmolid "index $atom"]
        set atombonds [join [$selatom getbonds]]
        if {[llength $atombonds] == 0} {
          continue
        }
        set bondsel [atomselect $tmpmolid "index $atombonds and occupancy 0.5"]
        set bondind [$bondsel get index]
        set bdninc 0
        while {$numbonds  < 0 &&  $bdninc < [llength $bondind]} {
          lappend delatom [lindex $bondind $bdninc]
          incr bdninc
          incr numbonds
        }
        $bondsel delete
        $selatom delete
       }
    }

   $all setbondorders $bondorders
   $all delete
   update_openvalence

    if {[llength $delatom] > 0} {
      set oldpicklist $picklist
      set picklist $delatom
      ::Molefacture::del_atom_gui
      set picklist $oldpicklist
      set delatom [list]
    }
}

proc ::Molefacture::raise_bondorder {} {
   variable tmpmolid

  variable picklist
   
   if {[Molefacture::isbond] == 0} {
    return
   }
   Molefacture::saveState
   set all [atomselect $tmpmolid "all"]
   set bondorders [$all getbondorders]
   # variable atomlist
   set sb0 [lindex [$all getbonds] [lindex $picklist 0]]
   set sb1 [lindex [$all getbonds] [lindex $picklist 1]]
   set vmdind0 [lindex $picklist 0]
   set vmdind1 [lindex $picklist 1]

   set pos0in1 [lsearch $sb1 $vmdind0]

   set pos1in0 [lsearch $sb0 $vmdind1]
   set currorder [lindex $bondorders $vmdind0 $pos1in0]
   if {$currorder < 0} {set currorder 1}
   set currorder [expr $currorder + 1]

   lset bondorders $vmdind0 $pos1in0 $currorder
   lset bondorders $vmdind1 $pos0in1 $currorder

    set delatom [list]
    foreach atom $picklist {
      set numbonds [Molefacture::newbond $atom]

       if { $numbonds <= -1} {
        set selatom [atomselect $tmpmolid "index $atom"]
        set atombonds [join [$selatom getbonds]]
        if {[llength $atombonds] == 0} {
          continue
        }
        set bondsel [atomselect $tmpmolid "index $atombonds and occupancy 0.5"]
        set bondind [$bondsel get index]
        set bdninc 0
        while {$numbonds  < 0 &&  $bdninc < [llength $bondind]} {
          lappend delatom [lindex $bondind $bdninc]
          incr bdninc
          incr numbonds
        }
        $bondsel delete
        $selatom delete
       }
    }

   $all setbondorders $bondorders
   $all delete
   update_openvalence


    if {[llength $delatom] > 0} {
      set oldpicklist $picklist
      set picklist $delatom
      ::Molefacture::del_atom_gui
      set picklist $oldpicklist
      set delatom [list]
    }
}

proc ::Molefacture::lower_bondorder {} {
   variable tmpmolid

    variable picklist
   
   if {[Molefacture::isbond] == 0} {
    return
   }
   Molefacture::saveState
   set all [atomselect $tmpmolid "all"]
   set bondorders [$all getbondorders]
   variable atomlist
   set sb0 [lindex [$all getbonds] [lindex $picklist 0]]
   set sb1 [lindex [$all getbonds] [lindex $picklist 1]]
   set vmdind0 [lindex $picklist 0]
   set vmdind1 [lindex $picklist 1]

   set pos0in1 [lsearch $sb1 $vmdind0]

   set pos1in0 [lsearch $sb0 $vmdind1]
   set currorder [lindex $bondorders $vmdind0 $pos1in0]
   if {$currorder < 0} {set currorder 1}
   if {$currorder == 1} {
     return
   }
   set currorder [expr $currorder - 1]
   lset bondorders $vmdind0 $pos1in0 $currorder
   lset bondorders $vmdind1 $pos0in1 $currorder

   $all setbondorders $bondorders
   $all delete
   
   update_openvalence

}

########################################################
### Add or delete bond between two atoms selected in ###
### the atomTable                                    ###
### opt == add - add bond                            ###
### opt == del - delete bond                         ###
########################################################
proc Molefacture::addDelBond { opt {dosavestate 1}} {
  variable tmpmolid
  variable picklist

  if {[llength $picklist] == 2} {
    if {$dosavestate} {
      Molefacture::saveState
    }

    set index1 [lindex $picklist 0]
    set index2 [lindex $picklist 1]

    ### the command getbonds return the bonds of the atoms order by index
    set selindex [lsort -integer [list $index1 $index2]]
    set sel1 [atomselect $tmpmolid "index [lindex $selindex 0]"]

    set bonds [lindex [$sel1 getbonds] 0]

    #### variable to store the indexes to delete when adding bonds
    set bondtodelete {}
    ### This is the only difference between add and delete bonds

    if {$opt == "add"} {
      set delatom [list]
      if {[lsearch $bonds [lindex $selindex 1]] == -1} {
        
         foreach atom $selindex {

          set numbonds [Molefacture::newbond $atom]

           if { $numbonds <= -1} {
            set selatom [atomselect $tmpmolid "index $atom"]
            set atombonds [join [$selatom getbonds]]
            if {[llength $atombonds] == 0} {
              continue
            }
            set bondsel [atomselect $tmpmolid "index $atombonds and occupancy 0.5"]
            set bondind [$bondsel get index]
            set bdninc 0
            while {$numbonds  < 0 &&  $bdninc < [llength $bondind]} {
              lappend delatom [lindex $bondind $bdninc]
              incr bdninc
              incr numbonds
            }
            $bondsel delete
            $selatom delete
           }
         }
        vmd_addbond $index1 $index2

        set tbinde [$Molefacture::atomTable curselection]
        set segnames [$Molefacture::atomTable getcolumns Segname]
        set segname1 [$Molefacture::atomTable cellcget [lindex $tbinde 0],Segname -text]
        set segname2 [$Molefacture::atomTable cellcget [lindex $tbinde 1],Segname -text]
        if {$segname1 != $segname2} {
          set seg1num [llength [lsearch -all $segnames $segname1]]
          set seg2num [llength [lsearch -all $segnames $segname2]]
          set newsegname $segname1
          if {$seg2num > $seg1num && $seg2num > 1} {
            set newsegname $segname2
          }
          set selall [atomselect $tmpmolid "$Molefacture::notDumSel"]
          $selall set segname $newsegname
          $selall delete
          ::Molefacture::updateAtomTable
        }
      }
      if {[llength $delatom] > 0} {
        set oldpicklist $picklist
        set picklist $delatom
        ::Molefacture::del_atom_gui
        set picklist $oldpicklist
        set delatom [list]
      }
    } else {
      if {[lsearch $bonds [lindex $selindex 1]] != -1} {
        vmd_delbond $index1 $index2
      }
    }

    set Molefacture::topocgenff 0

    update_openvalence
    $sel1 delete
    mol reanalyze $tmpmolid
  }
}



proc ::Molefacture::del_atom {delindex {updateopenv 1}} {
  variable tmpmolid

  variable picklist

    set currdel [atomselect $tmpmolid "index [join $delindex]"]
    $currdel set occupancy 0.0
    $currdel set resname "DUM"
    $currdel set $Molefacture::modFlagPSF 0
    foreach bondedatom [$currdel getbonds] index [$currdel get index] element [$currdel get element] {

      foreach bonds $bondedatom {
        vmd_delbond $bonds $index
        vmd_delbond $index $bonds
        if {$element == "H"} {
          set selaux [atomselect $tmpmolid "index $bonds"]
          $selaux set $Molefacture::modFlagH 0
          $selaux delete
        }
      }

    }

    $currdel delete

   ### reset flags of adding atoms and psf

  if {$updateopenv == 1} {
    update_openvalence
  }

  return
}


###############################################
# New proc to change the element of an atom   #
###############################################
proc ::Molefacture::edit_atom_elem {index atmnum } {
  variable tmpmolid
  variable periodic
  variable oxidation
  variable valence
  Molefacture::saveState
  set sel [atomselect $tmpmolid "index $index"]
  set type [lindex $periodic $atmnum]
  set newelement "[string index $type 0][string tolower [string index $type 1] ]"
  $sel set element $newelement
  set atmnames [$Molefacture::atomTable getcolumns AtmNAME]
  set ind 0
  set name "${type}$ind"
  while {[lsearch $atmnames "$name"] != -1} {
    set name "${type}$ind"
    if {[string length $name] > 4} {
      incr ind -1
      set name "${type}$ind"
      break
    }
    incr ind
  }
  $sel set name $name
  $sel set type $type

  if {$newelement == "H"} {
    $sel set occupancy 0.5
  } else {
    $sel set occupancy 0.8
  }

  topo -molid $tmpmolid -sel "index $index" guessatom radius element
  topo -molid $tmpmolid -sel "index $index" guessatom mass element
  $sel delete

  ### update table with new name, type and element

  set indexlis [$Molefacture::atomTable getcolumns 0]
  set tbindex [lsearch $indexlis $index]
  set oxlist [lindex [array get valence $newelement] 1]

  set oxi [lindex $oxlist [lindex $oxidation $index]]

  if {$oxi == ""} {
    set oxi "0"
  }

  $Molefacture::atomTable cellconfigure $tbindex,AtmNAME -text $name
  $Molefacture::atomTable cellconfigure $tbindex,AtmType -text $type
  $Molefacture::atomTable cellconfigure $tbindex,Element -text $newelement
  $Molefacture::atomTable cellconfigure $tbindex,OxState -text $oxi


}

#Procs to edit an atom's properties
proc ::Molefacture::edit_atom {} {
  variable editatom_name
  variable editatom_type
  variable editatom_element
  variable editatom_charge
  variable editatom_index
  variable tmpmolid
  variable picklist
  variable periodic

  #Make the appropriate changes to the selected atom
  set tmpsel [atomselect $tmpmolid "index [lindex $picklist 0]"]
  puts "|$editatom_name|"
  $tmpsel set name "\"$editatom_name\""
  $tmpsel set type "\"$editatom_type\""
  $tmpsel set charge $editatom_charge
  $tmpsel set element "\"$editatom_element\""
  $tmpsel delete

  #Update the menus and such
  update_openvalence
}


proc ::Molefacture::edit_incoming_atoms_FEP {} {
  variable editatom_beta
  variable tmpmolid
  variable picklist
  if {[llength $picklist] < 1} {
    return
  }
  set tmpsel [atomselect $tmpmolid "index $picklist"]
  $tmpsel set beta 1
  $tmpsel delete
  update_openvalence_FEP
}

proc ::Molefacture::edit_outgoing_atoms_FEP {} {
  variable editatom_beta
  variable tmpmolid
  variable picklist
  if {[llength $picklist] < 1} {
    return
  }
  set tmpsel [atomselect $tmpmolid "index $picklist"]
  $tmpsel set beta -1
  $tmpsel delete
  update_openvalence_FEP
}

proc ::Molefacture::edit_clear_atoms_FEP {} {
  variable editatom_beta
  variable tmpmolid
  variable picklist
  if {[llength $picklist] < 1} {
    return
  }
  set tmpsel [atomselect $tmpmolid "index $picklist"]
  $tmpsel set beta 0
  $tmpsel delete
  update_openvalence_FEP
}


proc ::Molefacture::unique_atomnames {} {
  variable tmpmolid

  set worksel [atomselect $tmpmolid "occupancy >= 0.5"]
  set workelems [$worksel get element]
  set workresnames [$worksel get resname]
  set worknames [$worksel get name]
  set workresid [$worksel get resid]
  array set elemcounts [list]

  foreach elem [lsort -unique $workelems] {
    array set elemcounts [list $elem 0]
  }

  array set namecounts [list]
  set nameslist [list]
  foreach name $worknames resid $workresid {
    set index [lsearch $nameslist "${name}_${resid}"]
    if {$index == -1 } {
      lappend nameslist "${name}_$resid"
    } elseif {$index > -1} {
      array set namecounts [list ${name}_$resid 0]
    }

  }

  set newnames [list]
  set namechanges ""
  foreach ind [$worksel get index] elem $workelems name $worknames resid $workresid {
    set elemcount [array get elemcounts $elem]
    set namecount [array get namecounts ${name}_$resid]
    set nameval [lindex $namecount 1]
    if {[llength $nameval] == 0} {
      continue
    }
    incr nameval
    set val [lindex $elemcount 1]
    incr val
    if { $nameval > 1 } {
        set tmpname "${name}${val}"
        if { [string length $tmpname] > 4 } {
            set tmpname "${elem}${val}"
#            lappend newnames "${elem}${val}"
        }
        lappend newnames $tmpname

        append namechanges "    $ind: $name -> ${tmpname}\n"
    } else {
        lappend newnames "${name}"
    }
    array set elemcounts [list $elem $val]
    array set namecounts [list $name $nameval]
  }

  #puts $newnames
  set proceed "yes"
  if { $namechanges != "" } {
    set proceed [tk_messageBox -type yesno -title "Atom names change" -icon warning \
        -message "Not all the atom names were unique. To fix this the\nfollowing atom names have to be changed in the top file:\n${namechanges}\nDo you wish to proceed?"]
    if {$proceed == "yes"} {
        $worksel set name $newnames
    } else {
      return $proceed
    }
  }
  $worksel delete

  update_openvalence
  return $proceed
}
