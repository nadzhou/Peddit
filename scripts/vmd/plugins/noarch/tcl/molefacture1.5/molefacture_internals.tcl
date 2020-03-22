##
## $Id: molefacture_internals.tcl,v 1.52 2019/07/18 15:16:11 jribeiro Exp $
##

proc ::Molefacture::measure_angle {coor1 coor2 coor3} {
  #Nondestructively measures the angle a1-a2-a3
  #puts "$coor1 $coor2 $coor3"

  set vec1 [vecsub $coor1 $coor2]
  set vec2 [vecsub $coor3 $coor2]
#puts "Vectors: $vec1 $vec2"
  if {[veclength $vec1] != 0} {
    set vec1 [vecnorm $vec1]
  }
  if {[veclength $vec2] != 0} {
    set vec2 [vecnorm $vec2]
  }

#  puts "here"
  set dotprod [vecdot $vec1 $vec2]
#  puts "$dotprod"
  set angle [expr acos($dotprod)]
#  puts "$angle"
  set angle [expr ($angle * 180/3.14159)]
#  puts "Final: $angle"
  return $angle
}


###########################################################
# Returns a list of triples that form angles
###########################################################

 proc ::Molefacture::anglelist {} {
    variable tmpmolid
    variable anglelistformat
    set sel [atomselect $tmpmolid "occupancy >= 0.5"]
    set bondsperatom  [$sel getbonds]
    set atom1arr [$sel get index]

    set atom1 0
    set angles {}
    foreach bonds $bondsperatom atom1 $atom1arr {
       set atom1sel [atomselect $tmpmolid "index $atom1"]
       set atom1coor [lindex [$atom1sel get {x y z}] 0]
       $atom1sel delete
       foreach partner1 $bonds {
          set psel1 [atomselect $tmpmolid "index $partner1"]
          set psel1coor [lindex [$psel1 get {x y z}] 0]
          foreach partner2 $bonds {
            if {$partner1 < $partner2} {

               set psel2 [atomselect $tmpmolid "index $partner2"]
               # puts "DEBUG $atom1coor || $psel1coor || [lindex [$psel2 get {x y z}] 0]]"
               # set thisangle [measure_angle $atom1coor $psel1coor [lindex [$psel2 get {x y z}] 0]]
                # Use this version if we want to display the angle in the menu:
                #    lappend angles [format "$anglelistformat" $partner1 $atom1 $partner2  $thisangle]
               lappend angles [format "$anglelistformat" $partner1 $atom1 $partner2]

              $psel2 delete
            }
          }
          $psel1 delete
       }
    }

    return $angles
 }

##############################################################
# If bond does not yet exist in vmd's bond list then add it. #
##############################################################

proc ::Molefacture::vmd_addbond {atom0 atom1 {order 1}} {

  variable tmpmolid
  if {$atom0 > $atom1} {
    set atomaux $atom1
    set atom1 $atom0
    set atom0 $atomaux

  }
   set sel [atomselect $tmpmolid "index $atom0 $atom1"]
   set bondlist   [$sel getbonds]
   set bondorders [$sel getbondorders]

   set bondatom0  [lindex $bondlist   0]
   set bondorder0 [lindex $bondorders 0]
   set ind1 [lsearch $bondatom0 $atom1]
   if {$ind1<0} {
      set newlist  [join [join [list $bondatom0  $atom1]]]
      set neworder [join [join [list $bondorder0 $order]]]
      #puts "newlist1: $newlist"
      lset bondlist   0 $newlist
      lset bondorders 0 $neworder
   }

   set bondatom1  [lindex $bondlist   1]
   set bondorder1 [lindex $bondorders 1]
   set ind2 [lsearch $bondatom1 $atom0]
   if {$ind2<0} {
      set newlist  [join [join [list $bondatom1  $atom0]]]
      set neworder [join [join [list $bondorder1 $order]]]
      #puts "newlist2: $newlist"
      lset bondlist   1 $newlist
      lset bondorders 1 $neworder
   }

   $sel setbonds      $bondlist
#   puts "$atom0 | $atom1 | $bondorders"
   $sel setbondorders $bondorders
   $sel delete
   if {$ind1>=0 && $ind2>=0} { return 1 }



   return 0
}

proc ::Molefacture::vmd_delbond {atom0 atom1} {
   #GOTCHA WARNING: This only deletes the atom0-atom1 bond, NOT the atom1-atom0
   #bond. This is not a bug but a feature, useful in deleting atoms.
   variable tmpmolid

   if {$atom0 > $atom1} {
    set atomaux $atom1
    set atom1 $atom0
    set atom0 $atomaux

  }
   set sel [atomselect $tmpmolid "index $atom0 $atom1"]
   set bondlist   [$sel getbonds]
   set bondorders [$sel getbondorders]
   set bondatom0  [lindex $bondlist  0]
   set bondorder0 [lindex $bondorders 0]

   set bondatom1  [lindex $bondlist  1]
   set bondorder1 [lindex $bondorders 1]

   set ind1 [lsearch $bondatom0 $atom1]
   set ind0 [lsearch $bondatom1 $atom0]
   #puts "Deleting: $atom0 $atom1"
   #puts "Old: $bondatom0"
   if {$ind1>-1} {
      set newlist  [join [join [lreplace $bondatom0  $ind1 $ind1]]]
      set neworder [join [join [lreplace $bondorder0 $ind1 $ind1]]]
      #puts "New: $newlist"
      lset bondlist   0 $newlist
      lset bondorders 0 $neworder
  }
  if {$ind0>-1} {
      set newlist  [join [join [lreplace $bondatom1  $ind0 $ind0]]]
      set neworder [join [join [lreplace $bondorder1 $ind0 $ind0]]]
      #puts "New: $newlist"
      lset bondlist   1 $newlist
      lset bondorders 1 $neworder
  }

   $sel setbonds      $bondlist
   $sel setbondorders $bondorders
   $sel delete
  return 0
}



proc ::Molefacture::get_element {name resname mass} {
  #Returns the atomic symbol for the element it is called on
  #Tries to use topology file and then mass for this determination
  variable toplist
  variable periodic

  set restopindex [lsearch [::Toporead::topology_get names $toplist] $resname]
  if {$restopindex != -1} {
    set resatomlist [join [lindex [::Toporead::topology_get residues $toplist] $restopindex 3]] ;# This is one list with the atom/type/charge info
    set type [get_atom_type $name $resatomlist]
    if {$type != -1} {
    set element [get_element_type $type [::Toporead::topology_get types $toplist]]
    if {$element != -1} {return $element}
  }
  }
  #If this has failed, try to assign by name

  regexp {([A-Z|a-z]+)\d*} $name -> namehead

  set namehead [string toupper $namehead]
  if {[lsearch $periodic $namehead] >= 0} {
    return $namehead
  }
  #Else, try assignment based on mass
  #WARNING: This is bogus for anything beyond the common organic stuff
  if {$mass != 0} {
    return [lindex $periodic [expr int($mass/2.0+0.5)]]
  }
  return "X"
}

proc ::Molefacture::get_atom_type {name atoms} {
  #This function takes a nested list of atoms, with each atom defined by a
  #list of form {name type charge}, and returns the type of the atom with
  #name $name
  foreach atom $atoms {
    if {[lindex $atom 0] == $name} {return [lindex $atom 1]}
  }
  return -1
}

proc ::Molefacture::get_element_type {type types} {
  #Returns the element name of the type $type

  foreach atom $types {
    if {[lindex $atom 0] == $type} {return [lindex $atom 2]}
  }

  return -1
}

proc ::Molefacture::new_mol {} {
  # Create a new molecule with some dummy atoms, and transition molefacture
  # to work on it
  ## No need to use hydrogens as dummys. Use dummys instead. Less confusion for the user
  # variable availabledummy
  variable protlastdummy
  variable tmpmolid
  variable nuclastdummy
  variable dummyincr

  set protlastdummy -1
  set nuclastdummy -1

  # Create a new  molecule with a pool of dummy atoms
  # set tmpmolid [mol new "Molefacture_newmol.xbgf"]
  set tmpmolid [mol new atoms $dummyincr]
  set atomslist [list]
  for {set index 1} {$index<=$dummyincr} {incr index} {
    # set name "HM$hind"
    ## list format {radius resname resid chain x y z occupancy element segname}
    lappend atomslist [list 0 DUM 9999 X 0.00 0.00 0.00 0 DUM MOLF]
  }

  # initialize the dummy atoms pool
  set sel [atomselect $tmpmolid "all"]
  animate dup [molinfo top]
  $sel set {radius resname resid chain x y z occupancy element segname} $atomslist

  ### Initialize the flags with 0
  $sel set $Molefacture::modFlagPSF 0
  $sel set $Molefacture::modFlagH 0

  mol reanalyze $tmpmolid
  # set availabledummy $dummyincr

  $sel set user -99
  $sel delete
  # Undraw all other molecules in VMD:
  foreach m [molinfo list] {
    if {$m==$tmpmolid} { molinfo $m set drawn 1; continue }
    molinfo $m set drawn 0
  }

  mol rename top "Molefacture_newmol"

  mol addrep top
  mol selection      "not resname DUM"
  mol representation "Bonds 0.1"
  mol color          Element
  mol modrep 0 top
  mol representation "VDW 0.1"
  mol addrep top

   variable oxidation
   set oxidation [list]
   for {set i 0} {$i < $dummyincr} {incr i} {
     lappend oxidation 0
   }

}
###########################################################
## Load the topologies and parameters shipped with VMD   ##
## inside ReadCharmmTop and ReadCharmmPar plugins +      ##
## the stream files added to QwikMD if exists            ##
###########################################################

proc Molefacture::loadDefaulTopo { {reslist ""} } {
   ### Use QwikMD commands to load Topologies to the
    ### variable QWIKMD::topoinfo
    ### if qwikmd package is not loaded already
    variable selTop
    set unknownRes ""
    set knownRes ""
    set pwd [pwd]
    global env
    if {$Molefacture::deftopinfo == -1} {
        if {[lindex [package present qwikmd] 0] == "package"} {
          package require qwikmd
        }

        set i 0

        if {[QWIKMD::checkDeposit 0] == 0} {
          
          QWIKMD::resetBtt 2 0
          
          array set ::Molefacture::selDefTop ""
          set Molefacture::deftopinfo $QWIKMD::topoinfo
        } else {
          ### In case the QwikMD loading topology process fails
          set TopList [file join $env(CHARMMTOPDIR) top_all36_prot.rtf]
          lappend TopList [file join $env(CHARMMTOPDIR) top_all36_lipid.rtf]
          lappend TopList [file join $env(CHARMMTOPDIR) top_all36_na.rtf]
          lappend TopList [file join $env(CHARMMTOPDIR) top_all36_carb.rtf]
          lappend TopList [file join $env(CHARMMTOPDIR) top_all36_cgenff.rtf]
          lappend TopList [file join $env(CHARMMTOPDIR) toppar_all36_carb_glycopeptide.str]
          lappend TopList [file join $env(CHARMMTOPDIR) toppar_water_ions_namd.str]
          set Molefacture::deftopinfo [list]
          foreach topo $TopList {
              if {[file exists $topo]} {
                  set handler [::Toporead::read_charmm_topology $topo 1]
                  lappend Molefacture::deftopinfo [::Toporead::topology_from_handler $handler]
              }
          }
        }

      $Molefacture::defTopMenu delete 0 end
      foreach top $Molefacture::deftopinfo {
        set top [lindex $top 0]
        $Molefacture::defTopMenu add checkbutton -label [file tail $top] -variable Molefacture::selDefTop($i,1) -command Molefacture::updateDefTop
        set Molefacture::selDefTop($i,0) $top
        set Molefacture::selDefTop($i,1) 1
        incr i
      }
      $Molefacture::defTopMenu add command -label "All" -command {Molefacture::updateDefTop "all"}
      $Molefacture::defTopMenu add command -label "None" -command {Molefacture::updateDefTop "none"}

    }

    if {$Molefacture::tmpmolid == "" || $Molefacture::tmpmolid == -1 ||\
     [lsearch [molinfo list] $Molefacture::tmpmolid] == -1} {
      return
    }
    set ind [lsearch -index 0 $Molefacture::deftopinfo "*top_all36_cgenff.rtf"]
    set selTop [list [lindex [lindex $Molefacture::deftopinfo $ind] 0]]
    
    foreach residue $reslist {
      set found 0
      foreach topo $Molefacture::deftopinfo {
        if { [::Toporead::topology_contains_resi $topo $residue] != -1 && $residue != "\*"} {
          lappend selTop [lindex $topo 0]
          set found 1
          break
        }
      }
      
      if {$found == 0 } {
        lappend unknownRes $residue
      } else {
        lappend knownRes $residue
      }
    }
    cd $pwd
    return [list $unknownRes $knownRes]
}

###################################################
# Update the default topology info (deftopinfo)   #
# array based on the options selected on the menu #
# Settings -> Default Topologies                  #
###################################################
proc Molefacture::updateDefTop { {allnone ""} } {
  ### array has two entries
  set asize [expr [array size Molefacture::selDefTop]/2]

  if {$allnone == "none"} {
    set Molefacture::deftopinfo ""
    for {set i 0} {$i < $asize} {incr i} {
      set Molefacture::selDefTop($i,1) 0
    }
    return
  } elseif {$allnone == "all"} {
    for {set i 0} {$i < $asize} {incr i} {
      set Molefacture::selDefTop($i,1) 1
    }
  }

  for {set i 0} {$i < $asize} {incr i} {
    if {$Molefacture::selDefTop($i,1) == 1 && \
      [lsearch -index 0 $Molefacture::deftopinfo "$Molefacture::selDefTop($i,0)"] ==-1} {
        set handler [::Toporead::read_charmm_topology $Molefacture::selDefTop($i,0) 1]
        lappend Molefacture::deftopinfo [::Toporead::topology_from_handler $handler]
    } elseif {$Molefacture::selDefTop($i,1) == 0 && \
      [lsearch -index 0 $Molefacture::deftopinfo "$Molefacture::selDefTop($i,0)"] !=-1} {
        set index [lsearch -index 0 $Molefacture::deftopinfo "$Molefacture::selDefTop($i,0)"]
        set Molefacture::deftopinfo  [lreplace $Molefacture::deftopinfo $index $index]
    }
  }
}

###################################################
# Reload the selection as a temporary molecule    #
# with a bunch of dummy atoms to add them to the  #
# actual molecule. Dummies have occupancy 0 while #
# all other atoms have occupancy 1.               #
###################################################

proc ::Molefacture::reload_selection { {runAutoPsf 1} } {
  
  global env
  variable tmpmolid
  variable origseltext
  variable origmolid
  variable oxidation
  variable valence
  variable unknownRes
  variable selTop
  ### Create a molecule using mol fromsels

  set molList [molinfo list]


  set oldmol $tmpmolid
  set poolsel [atomselect $tmpmolid "all"]
  set origsel [atomselect $origmolid $origseltext]
  set newatomnum [$origsel num]

  set knownRes [list]

  if {[lsearch [join [molinfo $origmolid get filetype]] "psf"] == -1 && $runAutoPsf == 1} {
      set reslist [lsort -unique -dictionary [$origsel get resname]]

     ### Molefacture::loadDefaulTopo returns two lists:
     ### index 0 with the unknown (not defined in the topologies) residues
     ### index 1 wirh the known (defined in the topologies) residues list
      set resultRes [Molefacture::loadDefaulTopo $reslist]
      set unknownRes [lindex $resultRes 0]
      set knownRes [lindex $resultRes 1]
      set knownsel ""
      if {[llength $knownRes] > 0} {

        set restring ""
        foreach resname $knownRes {
          append restring "\'$resname\' "
        }
        set sel [atomselect $origmolid "(resname ${restring} or water) and ($origseltext) "]
        set pwd [pwd]
        cd $env(TMPDIR)
        if {[file exist molefacture_tmp] != 1} {
          file mkdir molefacture_tmp
        }
        cd molefacture_tmp

        set autopsfLog ""
        $sel writepdb knownres.pdb
        set mol [mol new knownres.pdb]
        set atpsfOk [catch {autopsf -mol ${mol} -inlcude "all" -prefix knownres -top $selTop -regen -noterm} autopsfLog ]

        if {$atpsfOk >= 1} {
            puts "Molefacture:  Autopsf error $autopsfLog"
            set residue [lindex $knownRes end]
            set knownRes [lreplace $knownRes end end]
            lappend unknownRes $residue
        }


        $sel delete
       if {$atpsfOk < 1} {

         if {[file exists  knownres_formatted_autopsf.psf] != 1} {
            puts "Molefacture: Please inspect VMD outputs for more information."
         }
         set newmol [mol new knownres_formatted_autopsf.psf]
         mol addfile knownres_formatted_autopsf.pdb waitfor all $newmol

         set knownsel [atomselect $newmol "all"]
         
        }
       cd $pwd
      }
  }

    set selist [list]

   if {[llength $knownRes] > 0}  {
    lappend selist $knownsel

   }

   set unknownsel ""
   if {[llength $unknownRes] > 0} {
    set restring ""
    foreach resname $unknownRes {
      append restring "\'$resname\' "
    }
    set unknownsel [atomselect $origmolid "resname ${restring} and $origseltext"]
    set residListRename ""
    set residListDone ""
    foreach resid [$unknownsel get resid] segname [$unknownsel get segname] {
      if {[lsearch $residListDone $resid] == -1 && $segname == ""} {
        lappend residListRename $resid
      }
    }
    if {[llength $residListRename] > 0} {
      set tempSel [atomselect $origmolid "resname ${restring} and resid [join $residListRename]"]
      $tempSel set segname X
      $tempSel delete
    }
    lappend selist $unknownsel
   }

   if {[llength $knownRes] == 0 && [llength $unknownRes] == 0} {
    set selist [list $origsel $poolsel]
   } else {
    lappend selist $poolsel
   }

  foreach sel $selist {
    $sel global
  }


  Molefacture::mergeSels $molList $selist



  if {[llength $knownRes] > 0}  {
    set sel [atomselect $tmpmolid "resname [join $knownRes]"]
    $sel set $Molefacture::modFlagPSF 1
    $sel delete
    $knownsel delete
  }

  if {[llength $unknownRes] > 0} {
    $unknownsel delete
  }

}

####################################################
# Create the representations of the editing        #
# molecule                                         #
####################################################

proc Molefacture::createRep {} {
  variable tmpmolid
  mol rename $tmpmolid "Molefacture_newmol"
  mol addrep $tmpmolid
  mol selection      "not resname DUM"
  mol representation "Bonds 0.1"
  mol color          Element
  mol modrep 0 $tmpmolid
  mol representation "VDW 0.1"
  mol addrep $tmpmolid
  display resetview

  topo -molid $tmpmolid -sel "$Molefacture::notDumSel" guessatom element mass
  topo -molid $tmpmolid -sel "$Molefacture::notDumSel" guessatom radius element
  topo -molid $tmpmolid -sel "$Molefacture::notDumSel" guessatom mass element

  mol reanalyze $tmpmolid
}

########################################################
# Merge selections into molecules using mol fromsels   #
########################################################

proc Molefacture::mergeSels {molList selist} {

  variable tmpmolid
  variable oxidation


  set tmpmolid [mol fromsels $selist]

  set atmindex 0
  set sel [atomselect $tmpmolid "$Molefacture::notDumSel" ]

  foreach index [$sel get index] element [$sel get element] {
    #### Update the oxidation state of the new atoms
    #### if the multiple valence is allowed
    if {[info exists valence($element)] == 1 && [llength $valence($element)] > 1} {
      set numbondsel [atomselect $origmolid "index $index"]
      set numbonds [llength [lindex [$numbondsel getbonds] 0] ]
      set oxi 0
      for {set i 0} {$i < [llength $valence($element)]} {incr i} {

        if {$numbonds <= [lindex $valence($element) $i]} {
          set oxi $i
          break
        }
      }

      lset oxidation $atmindex $oxi
      $numbondsel delete
    }

    #### append 0 to the end of oxidation to account
    #### the additional atoms introduced
    lappend oxidation 0


    #### update Atoms' name to avoid duplicate names in the same residue

    set atomsel [atomselect $tmpmolid "index $index"]
    set name [::Molefacture::getAtomNewName [$atomsel get name] $element [$atomsel get resid]]
    if {$element == "H"} {
      $atomsel set occupancy 0.5
    } else {
      $atomsel set occupancy 0.8
    }
    $atomsel set name $name
    $atomsel delete

    incr atmindex

  }

  $sel delete

  ##delete temp molecules created by autopsf
  lappend molList $tmpmolid
  foreach mol [molinfo list] {
    if {[lsearch $molList $mol] == -1} { mol delete $mol }
  }

  #### Determine location of the atom repository
  set resid 9999

  #### Load the file including the dummy atom repository
  #### variable tmpmolid [mol new $tmpfile]

   #### Undraw all other molecules in VMD:
   foreach m [molinfo list] {
      if {$m==$tmpmolid} { molinfo $m set drawn 1; continue }
      molinfo $m set drawn 0
   }

  set sel [atomselect $tmpmolid "$Molefacture::notDumSel"]

  if {[lindex [$sel get atomicnumber] 0] == -1} {assign_elements}

  if {[lsort -unique -dictionary [$sel get segname]] == "{}"} {
    $sel set segname X
  }
  $sel delete

  ### Reset the user field

  set sel [atomselect $tmpmolid "$Molefacture::dumSel"]


  $sel set user "-99"
  $sel delete

  Molefacture::createRep

  Molefacture::capConnections

  update_openvalence


  #### run idatm to guess the bond orders. Although this is also an
  #### atom typing procedure, we only care about the bond order, as
  #### ultimately the atom's type will be fetched from the CGenFF
  #### TODO: Fix the double double bonds in case of phosphorus atoms
  #### in phosphates
  Molefacture::run_idatm

  Molefacture::saveState

}


####################################################
# Cap the connections between the selected atoms   #
# and the original molecule to prevent adding      #
# hydrogens or editing the bonds to the parent     #
# molecule                                         #
####################################################

proc Molefacture::capConnections {} {
  variable tmpmolid
  variable origseltext
  variable origmolid
  variable oxidation
  variable valence
  variable capAtoms
  variable toCapAtoms

  if {$origmolid == -1 || [lsearch [molinfo list] $origmolid] == -1} {
    return
  }
  set parentsel [atomselect $origmolid "all"]
  set plresidues [lsort -unique -integer [$parentsel get residue]]
  $parentsel delete

  set cursel [atomselect $origmolid "$origseltext"]
  set clresidues [lsort -unique -integer [$cursel get residue]]
  set clindex [$cursel get index]
  $cursel delete

  ### list storing the indexes to cap and bond orders
  set toCapAtoms [list]
  set capAtoms [list]
  foreach clres $clresidues {
    ### residues of the editing selection
    set sel [atomselect $origmolid "residue $clres"]

    ### indexes and bonds of the editing selection
    set selindex [$sel get index]
    set selbonds [$sel getbonds]
    $sel delete
    if {[llength [join $selbonds] ] > 0} {
      ### Bonds of the editing selection with parent molecule
      set clplconnect [atomselect $origmolid "index [join $selbonds] and not residue [join $clresidues]"]
      set clplbonds [$clplconnect get index]

      ### Find the atoms' index ($index) making bonds with parent
      ### molecule

      foreach clplb $clplbonds {
        set index 0
        foreach bonds $selbonds {

          if {[lsearch $bonds $clplb] > -1} {

            lappend toCapAtoms [lindex $selindex $index]
            lappend capAtoms "$clplb"
            break
          }
          incr index
        }

      }

      $clplconnect delete
    }
  }

  if {[llength $capAtoms] == 0} {
    return
  }

  ### Place the cap atoms and bond to the new molecule
  set dumsel [atomselect $tmpmolid "$Molefacture::dumSel or segname CAP"]
  $dumsel set segname MOLF
  $dumsel set resname DUM
  set dumindexlist [$dumsel get index]
  set dumindex [lindex $dumindexlist end]

  set origatom [atomselect $origmolid "index [join $capAtoms]"]
  set origname [$origatom get name]
  set origelement [$origatom get element]
  set origtype [$origatom get type]
  set origxyz  [$origatom get {x y z}]
  set origresname [$origatom get resname]
  set origresids [$origatom get resid]

  set origToCap [atomselect $origmolid "index [join $toCapAtoms]"]
  set toCapElement [$origToCap get element]
  $origToCap delete

  $origatom delete
  set index 0
  set newToCap ""
  foreach atms $toCapAtoms {

    set capsel [atomselect $tmpmolid "index $dumindex"]

    $capsel set segname "CAP"
    $capsel set occupancy 0.0
    $capsel set element [lindex $origelement $index]
    $capsel set name [lindex $origname $index]
    $capsel set {x y z} [list [lindex $origxyz $index]]
    $capsel set type [lindex $origtype $index]
    set resname "C[lindex $origresname $index]"

    if {[string length $resname] > 4} {
      set resname [string range 0 3]
    }
    $capsel set resname $resname
    $capsel set resid [lindex $origresids $index]

    ### Identify the atoms to cap based on the distance
    ### It would be great if we could mark the atoms and carry this information
    ### when using the mol fromsels command. For now, this is not possible
    set atm [atomselect $tmpmolid "(exwithin 2 of index $dumindex) \
    and element [lindex $toCapElement $index]"]
    set atmindex [$atm get index]
    if {[llength $atmindex] > 0} {


      $capsel set radius 7.5
      lappend newToCap $atmindex
      vmd_addbond $dumindex $atmindex
      vmd_addbond $atmindex $dumindex


      mol reanalyze $tmpmolid
      # mol bondsrecalc $tmpmolid
      $capsel delete
      incr dumindex -1
      incr index
    } else {
       $capsel set resname DUM
       $capsel set segname MOLF
    }
    $atm delete
  }
  set toCapAtoms $newToCap

}

###################################################
# Get color of the elements for the periodic table#
###################################################
### based on the code https://wiki.tcl-lang.org/page/Periodic+Table+of+Chemical+Elements
proc ::Molefacture::addElemet {atmnum} {
   variable periodic
   variable periodicNames
   variable periTab

   set w 30; ## width of elements boxes
   set h 30; ## height of elements

   ## labels
   # am -> Alkali metal
   # aem -> Alakline earth metal
   # tm -> transition metal
   # lnthd -> Lanthanide
   # pstm -> Pos-transition metal
   # mtoid -> metaloids
   # onm -> other non metaloids
   # h -> halogens
   # ng -> noble gases
   set font [font create -family Helvetica -size 7 -weight bold]
   set typelist [list {am #ffcc33} {aem #ffff44} {tm #dbbaba} {pstm #9bdccb} {mtoid #7bdc8d} {onm #3dff4b} {h #3eedcc} {nb #7ccafb} {lnthd #ffbb99}]

   set type ""

   set p 0
   set g 0
   if {$type == ""} {
      set typegp [::Molefacture::getElemtnType $atmnum]
      set type [lindex $typegp 0]
      set g [lindex $typegp 1]
      set p [lindex $typegp 2]
   }
   if {$type == ""} {
      tk_messageBox -title "Element Not Found" -icon error -message "Element Not Found" -type okcancel -parent $Molefacture::topGui
      return
   }
   if {$g > 18} {
      return
   }
   set bgcol [lindex [lindex $typelist [lsearch -index 0 $typelist $type]] 1]
   $periTab create rectangle 0 0 $w $h -tags ${atmnum}_box -fill $bgcol

   $periTab move ${atmnum}_box [expr ($g-1)*$w+10] [expr ($p -1)*$h+10]
   $periTab bind ${atmnum}_box <ButtonPress-1> "::Molefacture::elementPTBind $atmnum"

   set text [lindex $periodic $atmnum]
   $periTab create text [expr [expr ($g-1)*$w+10] + [expr $w/2.0]] [expr [expr ($p -1)*$h+10] + [expr $h/2.0]] -text $text -font $font -anchor center -tags ${atmnum}_main
   $periTab bind ${atmnum}_main <ButtonPress-1> "::Molefacture::elementPTBind $atmnum"

}
##################################################################
# get period and group from atomic number                        #
##################################################################

proc ::Molefacture::getElemtnType {atmnum} {
   set type ""
   set atmaux 1
   for {set p 1} {$p <= 7} {incr p} {
      for {set g 1} {$g <= 18} {incr g} {
         if {$atmnum == $atmaux} {
            if {$atmnum > 71} {
               incr g -1
               if {$g < 1} {
                  set g 18
               }
            }
            if {$p == 1 && $g == 1} {
               set type onm
            } elseif {$g == 1} {
               set type am
            } elseif {$g == 2} {
               set type aem
            } elseif {$g > 2 && $g <= 12} {
               set type tm
            } elseif {($g == 13 && $p >= 3) || ($g == 14 && $p >= 5) || \
               ($g == 15 && $p >= 6) || ($g == 16 && $p == 7)} {
               set type pstm
            } elseif {($g == 13 && $p == 2) || ($g == 14 && ($p == 3 || $p == 4)) ||\
               ($g == 15 && ($p == 4 || $p == 5)) || ($g == 16 && ($p == 5 || $p == 6))} {
               set type mtoid
            } elseif {($p == 2 && ($g == 14 || $g == 15 || $g == 16)) || \
               ($p == 3 && ($g == 15 || $g == 16)) || $p == 4 && $g == 16} {
               set type onm
            } elseif {$g == 17} {
               set type h
            } elseif {$g == 18} {
               set type nb
            }
         }
         if {$atmnum >=57 && $atmnum <=71} {
            set p 8
            set g [expr $atmnum - 57 + 3]
            set type lnthd
         } elseif {$atmnum > 71} {
            set p 6
            set g [expr $atmnum - 71 + 3]
         }
         if {$type != ""} {
            break
         }
         if {$atmaux == 1} {
            incr g 16
         } elseif {$atmaux == 4 || $atmaux == 12} {
            incr g 10
         }
         incr atmaux
      }
      if {$type != ""} {
         return "$type $g $p"
      }
   }
}


##################################################################
# command bond to the selection of element in the periodic table #
##################################################################

proc ::Molefacture::elementPTBind {atmnum} {
   variable periodic
   variable periodicNames
   variable mass_by_element
   variable periTab
   set w 75
   set h 70
   set x 180
   set y 20
   ## labels
   # am -> Alkali metal
   # aem -> Alakline earth metal
   # tm -> transition metal
   # lnthd -> Lanthanide
   # pstm -> Pos-transition metal
   # mtoid -> metaloids
   # onm -> other non metaloids
   # h -> halogens
   # ng -> noble gases

   set labellist "big_box element atmnumtxt mass name"

   foreach label $labellist {
      $periTab delete $label
   }
   set font1 [font create -family Helvetica -size 25 -weight bold]
   set font2 [font create -family Helvetica -size 10 -weight bold]
   set font3 [font create -family Helvetica -size 8 -weight bold]
   set typelist [list {am #ffcc33} {aem #ffff44} {tm #dbbaba} {pstm #9bdccb} {mtoid #7bdc8d} {onm #3dff4b} {h #3eedcc} {nb #7ccafb} {lnthd #ffbb99}]

   set type [lindex [::Molefacture::getElemtnType $atmnum] 0]

   set bgcol [lindex [lindex $typelist [lsearch -index 0 $typelist $type]] 1]
   $periTab create rectangle 0 0 $w $h -tags big_box -fill $bgcol

   $periTab move big_box $x $y

   set element [lindex $periodic $atmnum]
   set name [lindex $periodicNames $atmnum]
   set length 1
   if {[string length $element] > 1} {
      set element "[string index $element 0][string tolower [string index $element 1] ]"
      set length 2
   }
   set mass [lindex [array get mass_by_element $element] 1]
   set mainx [expr $x + [expr $w/2]]
   set mainy [expr [expr $y + [expr $h/2]] -10]

   $periTab create text $mainx $mainy -text $element -font $font1 -anchor center -tags element

   set fspace [font metrics $font2 -linespace]
   set xfspace $fspace
   if {$length == 2} {
      set xfspace [expr $fspace * 2]
   }
   $periTab create text [expr $mainx - $xfspace]  [expr $mainy - $fspace] -text $atmnum -font $font2 -anchor center -tags atmnumtxt

   set fspace [expr 1.6* [font metrics $font3 -linespace] ]

   $periTab create text $mainx  [expr $mainy + $fspace] -text $mass -font $font3 -anchor center -tags mass

   $periTab create text $mainx  [expr $mainy + 1.7*$fspace] -text $name -font $font3 -anchor center -tags name

   set Molefacture::elemcomb "${name}(${element}:${atmnum})"

   set index [$Molefacture::atomTable curselection]
   set atmnum [string trimright [lindex [split $Molefacture::elemcomb ":"] 1] ")"]
   #if {[llength $index] == 1} {
    foreach ind $index {
      ::Molefacture::edit_atom_elem [$Molefacture::atomTable cellcget $ind,Index -text ] $atmnum
    }

   #}

}

proc ::Molefacture::export_molecule {filename} {
   variable tmpmolid
   set sel [atomselect $tmpmolid "occupancy>=0.8"]
   write_xbgf $filename $sel
   $sel delete
}

####################################################
# Writes a selection into a xbgf file and appends  #
# VDW parameters as REMARKs.                       #
####################################################

proc ::Molefacture::write_xbgf {xbgffile sel} {
   $sel writexbgf $xbgffile

   # Find the END tag
   set fid [open $xbgffile r+]
   while {![eof $fid]} {
      set filepos [tell $fid]
      set line [gets $fid]
      if {[lindex $line 0]=="END"} { break }
   }
   seek $fid $filepos

   set VDWdata {}
   foreach remark [molinfo [$sel molid] get remarks] {
      if {[lindex $remark 0]=="VDW"} {
   lappend VDWdata [lrange $remark 1 end]
      }
   }
   set i 1
   foreach index [$sel list] {
      #puts      [format "VDW %6i %s" $i [lindex $VDWdata $index]]
      puts $fid [format "VDW %6i %s" $i [lindex $VDWdata $index]]
      incr i
   }
   variable chargelist
   set i 1
   foreach index [$sel list] {
      #puts      [format "LEWIS %6i %3s" $i [lindex $chargelist $index]]
      puts $fid [format "LEWIS %6i %3s" $i [lindex $chargelist $index]]
      incr i
   }
   puts $fid "END"

   close $fid
}

############################################################
### This function is invoked whenever an atom is picked. ###
############################################################

proc ::Molefacture::atom_picked_fctn { args } {
  global vmd_pick_atom
  global vmd_pick_shift_state
  variable picklist
  variable tmpmolid

  if {$tmpmolid == -1 || [lsearch [molinfo list] $tmpmolid] == -1} {
    return
  }
  ### Check if a CAP atom was selected
  set sel [atomselect $tmpmolid "index $vmd_pick_atom"]
  if {[$sel get segname] == "CAP"} {
    tk_messageBox -message "Cannot Select or Edit the Atoms Use To Mark the Terminals 
    of the Molecules." \
    -parent $Molefacture::topGui -title "Terminal Atom Selected"
     -type ok -icon info
     # return
  }
  set newatoms ""
  ### collect all picked atoms in picklist variable
  if {[lsearch $picklist $vmd_pick_atom] == -1} {
    lappend picklist $vmd_pick_atom
  }

  if {[info exists vmd_pick_shift_state] == 1 && !$vmd_pick_shift_state} {
   set picklist $vmd_pick_atom
   # label delete Atoms all
  }

  if {$Molefacture::mousemode == "addatom" || $Molefacture::mousemode == "delatom"} {
    set picklist [lindex $picklist end]
  }
  ### Check the mouse mode and perform the correspondent
  ### action
  switch $Molefacture::mousemode {
    addatom {
      set newatoms [Molefacture::addAtomMain]
      set picklist $newatoms
    }
    delatom {
      Molefacture::delAtomMain
    }
    addbond {
      Molefacture::addDelBond add
    }
    delbond {
      Molefacture::addDelBond del
    }
    incrorder {
      Molefacture::raise_bondorder
    }
    decrorder {
      Molefacture::lower_bondorder
    }
    default {}
  }
  ### update the Atom Table
  set newpicklist [list]
  set allatomlist [$Molefacture::atomTable getcolumns Index]
  foreach atomind $picklist {
    lappend newpicklist [lsearch $allatomlist $atomind]

  }

  $Molefacture::atomTable selection clear 0 end
  $Molefacture::atomTable selection set $newpicklist
  ::Molefacture::selectAtomTable

  ## highlight atoms in the OpenGL window
  draw_selatoms
}

############################################################
### Get world coordinates from pointer position          ###
############################################################

proc Molefacture::get_coor_from_pointer { x y molid } {
  set view [lindex [molinfo $molid get rotate_matrix] 0]
  set center [lindex [molinfo $molid get center_matrix] 0]
  set scale [lindex [molinfo $molid get scale_matrix] 0]

  set trans [transoffset "$x $y 0"]

  set trans [transmult $trans $scale $view $center]

  return [coordtrans $trans [lindex [molinfo $molid get center] 0] ]

}

##############################################################
### This function is invoked whenever an atom is picked to ###
### update the element combobox value                      ###
##############################################################

proc ::Molefacture::update_combo_sel_atoms { index } {
  set element [$Molefacture::atomTable cellcget $index,Element -text ]
  set atmnum [lsearch $Molefacture::periodic [string toupper $element]]
  set name [lindex $Molefacture::periodicNames $atmnum]
  set Molefacture::elemcomb "${name}(${element}:${atmnum})"
}
proc ::Molefacture::hilight_angle {atom1 atom2 atom3 angleindex} {
# Procedure to select an angle in molefacture, including updating the angle list
    # Blank all item backgrounds

  variable atomlist

      for {set i 0} {$i<[.molefac.angles.list.list index end]} {incr i} {
        .molefac.angles.list.list itemconfigure $i -background {}
      }
      # Get current selection index
      set selangle [.molefac.angles.list.list curselection]

      # Paint the background of the selected bond
      .molefac.angles.list.list itemconfigure $angleindex -background $::Molefacture::selectcolor
      .molefac.angles.list.list activate $angleindex
      .molefac.angles.list.list yview $angleindex

      # Get the selected bond
      set selindex [lrange [lindex $::Molefacture::anglelist $angleindex] 0 2]

      # Get information about this angle
      variable tmpmolid
      variable angle
      variable angleatom
      variable angleaxis
      variable anglesel
      variable anglepicklist
      variable anglemoveatom
      variable anglecoor
      # variable dihedmarktags

      set sel1 [atomselect $tmpmolid "index $angleatom(1)"]
      set sel2 [atomselect $tmpmolid "index $angleatom(2)"]
      set sel3 [atomselect $tmpmolid "index $angleatom(3)"]

      set coor2 [join [$sel2 get {x y z}]]
      set coor1 [join [$sel1 get {x y z}]]
      set coor3 [join [$sel3 get {x y z}]]

      set anglecoor(1) $coor1
      set anglecoor(2) $coor2
      set anglecoor(3) $coor3

#      puts "Subcoors: $coor1 $coor2 $coor3"
      set vec1 [vecsub $coor2 $coor1]
      set vec2 [vecsub $coor2 $coor3]

      set axis [veccross $vec1 $vec2]
#      puts "Doing norm of $axis"
      set angleaxis [vecscale [vecnorm $axis] 1.0]
#      puts "Done"

      # Generate two selections for the two molecule halves
      set indexes1 [join [::util::bondedsel $tmpmolid $angleatom(2) $angleatom(1) -maxdepth [llength $atomlist]]]
      set indexes2 [join [::util::bondedsel $tmpmolid $angleatom(2) $angleatom(3) -maxdepth [llength $atomlist]]]
      if {[havecommonelems $indexes1 $indexes2 [list $angleatom(1) $angleatom(2) $angleatom(3)]] > 0} {
        set indexes1 $angleatom(1)
        set indexes2 $angleatom(3)
      }
      if {[array exists anglesel]} { catch {$anglesel(1) delete}; catch {$anglesel(2) delete }}
      set mysel1 [uplevel 2 "atomselect $tmpmolid \"index $indexes1 and not index $angleatom(2)\""]
      set mysel2 [uplevel 2 "atomselect $tmpmolid \"index $indexes2 and not index $angleatom(2)\""]
      array set anglesel [list 1 $mysel1]
      array set anglesel [list 2 $mysel2]

      #Compute the angle
      set angle [measure angle [list [list $angleatom(1) $tmpmolid] [list $angleatom(2) $tmpmolid] [list $angleatom(3) $tmpmolid] ]]

      .molefac.angles.realangle.scale set $::Molefacture::angle
}

proc ::Molefacture::hilight_bond {atom1 atom2 bondindex} {
# Carry out all the background stuff required when you select a bond in molefacture
# Includes updating the bond list, highlighting dihedral, and so forth
# Must be given the two atom indices and the index of the bond in bondlist
  #puts "DEBUG: Running hilight_bond for $atom1 $atom2 $bondindex"

  variable atomlist

  # Set up the bond and dihedral values

  variable tmpmolid
  if {![llength $tmpmolid]} { return }
  variable bondcoor
  variable dihedatom
  set dihedatom(1) $atom1
  set dihedatom(2) $atom2
  set sel1 [atomselect $tmpmolid "index $dihedatom(1)"]
  set sel2 [atomselect $tmpmolid "index $dihedatom(2)"]
  set bondcoor(1) [join [$sel1 get {x y z}]]
  set bondcoor(2) [join [$sel2 get {x y z}]]
  set ::Molefacture::bondlength [veclength [vecsub $bondcoor(2) $bondcoor(1)]]

  # Choose a dihedral for this bond
  set bonds1 [lsearch -all -inline -not [join [$sel1 getbonds]] $dihedatom(2)]
  set bonds2 [lsearch -all -inline -not [join [$sel2 getbonds]] $dihedatom(1)]
  set dihedatom(0) [lindex $bonds1 0]
  set dihedatom(3) [lindex $bonds2 0]
  # Delete the old marks
  # variable dihedmarktags
  # foreach tag $dihedmarktags {
  #   graphics $tmpmolid delete $tag
  # }

   if {[llength $dihedatom(0)] && [llength $dihedatom(3)]} {
     #puts "dihedatom(0)=$dihedatom(0); [join [$sel1 getbonds]]; $bonds1"
     #puts "dihedatom(1)=$dihedatom(1)"
     #puts "dihedatom(2)=$dihedatom(2)"
     #puts "dihedatom(3)=$dihedatom(3); [join [$sel2 getbonds]]; $bonds2"
     set sel0 [atomselect $tmpmolid "index $dihedatom(0)"]
     set bondcoor(0) [join [$sel0 get {x y z}]]
     set sel3 [atomselect $tmpmolid "index $dihedatom(3)"]
     set bondcoor(3) [join [$sel3 get {x y z}]]

     #puts "DEBUG: Generating molecule half selections"
     # Generate two selections for the two molecule halves
     variable bondsel
     set indexes1 [join [::util::bondedsel $tmpmolid $dihedatom(2) $dihedatom(1) -maxdepth [llength $atomlist]]]
     set indexes2 [join [::util::bondedsel $tmpmolid $dihedatom(1) $dihedatom(2) -maxdepth [llength $atomlist]]]
     if {[havecommonelems $indexes1 $indexes2 [list $dihedatom(1) $dihedatom(2)]] > 0} {
       set indexes1 $dihedatom(1)
       set indexes2 $dihedatom(2)
 }
     #puts [array get bondsel]
 if {[array exists bondsel]} {catch {$bondsel(1) delete}; catch {$bondsel(2) delete}}
 set mysel1 [uplevel 2 "atomselect $tmpmolid \"index $indexes1 and not index $dihedatom(2)\""]
 set mysel2 [uplevel 2 "atomselect $tmpmolid \"index $indexes2 and not index $dihedatom(1)\""]
 array set bondsel [list 1 $mysel1]
 array set bondsel [list 2 $mysel2]

 # Compute the bond dihedral angle
 set ::Molefacture::dihedral [measure dihed [list [list $dihedatom(0) $tmpmolid] [list $dihedatom(1) $tmpmolid] [list $dihedatom(2) $tmpmolid] [list $dihedatom(3) $tmpmolid] ] ]
 }

   variable w
   $w.bonds.f2.angle.scale set $::Molefacture::dihedral

}

proc ::Molefacture::draw_selatoms {} {
   variable tmpmolid

   if {$tmpmolid == -1 || [lsearch [molinfo list] $tmpmolid] == -1} {
    return
   }
   # Delete any previously generated atomlabels
   variable atommarktags
   variable anglemarktags
   foreach tag $atommarktags  {
      graphics $tmpmolid delete $tag
   }
  set atommarktags {}

   variable labelradius
   variable picklist
   set selatoms $picklist;
   graphics $tmpmolid color orange
   graphics $tmpmolid materials on
    graphics $tmpmolid material "Transparent"
   foreach atom $selatoms {
      set sel [atomselect $tmpmolid "index $atom"]
      # lappend atommarktags [graphics $tmpmolid ]
      lappend atommarktags [graphics $tmpmolid sphere [join [$sel get {x y z}]] radius [expr $labelradius*1.5] resolution 15]
      $sel delete
   }

}


##############################################################
## Evaluate if the atom groups were precomputed already     ##
##############################################################

proc ::Molefacture::predraw_angleatoms {} {
  variable picklist
  set checedkangle 0
  if {$Molefacture::anglemvgrp1list == -1} {
    if {[llength $picklist] == 3} {
      set checedkangle [Molefacture::isanglebond]
    } elseif {[llength $picklist] == 4} {
      set checedkangle [Molefacture::isdihbond]
    }
  } elseif {[llength $Molefacture::anglemvgrp1list] > 0 || [llength $Molefacture::anglemvgrp2list] > 0} {
    set checedkangle 1
  }
  if {$checedkangle == 1} {
    Molefacture::draw_angleatoms
    return 1
  } else {
    return 0
  }
}


##############################################################
## Draw the atoms forming the two groups to be moved when   ##
## changing the angle                                       ##
##############################################################

proc ::Molefacture::draw_angleatoms {} {
  variable tmpmolid
  variable anglemarktags
  variable labelradius
  foreach tag $anglemarktags  {
    graphics $tmpmolid delete $tag
  }

  set anglemarktags {}
  graphics $tmpmolid materials on
  graphics $tmpmolid material "Transparent"
  if {($Molefacture::anglemvgrp1 == 1 || $Molefacture::dihedrgrp == "group1") && [llength $Molefacture::anglemvgrp1list] > 0 && $Molefacture::anglemvgrp1list != -1} {
    foreach atom $Molefacture::anglemvgrp1list {
      set sel [atomselect $tmpmolid "index $atom"]
      # lappend anglemarktags [graphics $tmpmolid]
      graphics $tmpmolid color green
      lappend anglemarktags [graphics $tmpmolid sphere [join [$sel get {x y z}]] radius [expr $labelradius*1.6] resolution 15]
      $sel delete
   }
  }

  if {($Molefacture::anglemvgrp2 == 1 || $Molefacture::dihedrgrp == "group2") && [llength $Molefacture::anglemvgrp2list] >0 && $Molefacture::anglemvgrp2list != -1} {
    foreach atom $Molefacture::anglemvgrp2list {
      set sel [atomselect $tmpmolid "index $atom"]
      # lappend anglemarktags [graphics $tmpmolid ]
      graphics $tmpmolid color blue
      lappend anglemarktags [graphics $tmpmolid sphere [join [$sel get {x y z}]] radius [expr $labelradius*1.6] resolution 15]
      $sel delete
    }
  }

}

proc ::Molefacture::flag_incomingatoms {} {
   variable tmpmolid
   variable fepmolid
   # Delete any previously generated atomlabels
   variable atomarktags

   if {![winfo exists .molefac.val.list.list]} { return }
   variable labelradius
   variable picklist
   set selatoms $picklist; #[.molefac.val.list.list curselection]
   foreach atom $selatoms {
     set newmark [list]
     set sel [atomselect $tmpmolid "index $atom"]
     lappend newmark $atom
     lappend newmark "green"

     $sel delete
     lappend atomarktags $newmark
   }
   edit_incoming_atoms_FEP
}

proc ::Molefacture::flag_outgoingatoms {} {
   variable tmpmolid
   variable fepmolid
   # Delete any previously generated atomlabels
   variable atomarktags

   if {![winfo exists .molefac.val.list.list]} { return }
   variable labelradius
   variable picklist
   set selatoms $picklist; #[.molefac.val.list.list curselection]
   foreach atom $selatoms {
     set newmark [list]
     set sel [atomselect $tmpmolid "index $atom"]
     lappend newmark $atom
     lappend newmark "red"
     $sel delete
     lappend atomarktags $newmark
   }
   edit_outgoing_atoms_FEP
}

proc ::Molefacture::flag_commonatoms {} {
   variable tmpmolid
   variable fepmolid
   variable atomarktags
   if {![winfo exists .molefac.val.list.list]} { return }
   variable labelradius
   variable picklist
   set selatoms $picklist; #[.molefac.val.list.list curselection]
   foreach atom $selatoms {
     for {set i 0} {$i < [llength $atomarktags]} {incr i} {
       if {[lindex $atomarktags $i 0] == $atom} {
         set atomarktags [lreplace $atomarktags $i $i]
       }

    }
   }
   edit_clear_atoms_FEP
}


proc ::Molefacture::edit_update_list {taglist} {

   # This function generates a formatted list from the columns in taglist
   # It returns both a formatted header for a table using the given column types
   # and a string and list containing the appropriate formatting codes

   set formatstring {}
   set formatfield {}

   foreach tag $taglist {
      switch $tag {
   Index {
      append formatstring {%5s }
      append formatfield  {%5s }
   }
         Atom1 {
            append formatstring {%5s }
            append formatfield  {%5s }
         }
         Atom2 {
            append formatstring {%5s }
            append formatfield  {%5s }
         }
         Atom3 {
            append formatstring {%5s }
            append formatfield  {%5s }
         }
   Elem {
      append formatstring {%4s }
      append formatfield  {%4s }
   }
   Name    {
      append formatstring {%4s }
      append formatfield  {%4s }
   }
   Flags    {
      append formatstring {%5s }
      append formatfield  {%5s }
   }
   Type    {
      append formatstring {%4s }
      append formatfield  {%4s }
   }
   Rname {
      append formatstring {%5s }
      append formatfield  {%5s }
   }
   Resid   {
      append formatstring {%5s }
      append formatfield  {%5s }
   }
   Segid   {
      append formatstring {%5s }
      append formatfield  {%5s }
   }
   VDWeps  {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   Order  {
      append formatstring {%5s }
      append formatfield  {% 4.1f }
   }
   VDWrmin {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   FormCharge     {
      append formatstring {% 12s }
      append formatfield  {% 12.4f }
   }
   Charge     {
      append formatstring {%7s}
      append formatfield  {% 7.4f }
   }
   Lewis {
      append formatstring {%5s }
      append formatfield  {%5i }
   }
   Mullik     {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   MulliGr {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   NPA    {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   SupraM    {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   ESP     {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   RESP    {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   ForceF    {
      append formatstring {%7s }
      append formatfield  {% 7.4f }
   }
   Open    {
      append formatstring {%5s}
      append formatfield  {% 4i }
   }
   OxState    {
      append formatstring { %7s}
      append formatfield  {% 6i }
   }
   FEPindex  {
      append formatstring {%9s }
      append formatfield  {%7.1f }
   }
      }
   }

   # Build the formatted atomproptags
   set formatcommand "\[format \"$formatstring\" $taglist\]"
   eval set taglist "$formatcommand"
   #puts "Taglist: $taglist Format: $formatstring"
   set returnlist [list]
   lappend returnlist $taglist
   lappend returnlist $formatfield

   return $returnlist

}


proc ::Molefacture::noblanks {mylist} {
  set newlist [list]
  foreach elem $mylist {
    if {$elem != ""} {lappend newlist $elem}
  }

  return $newlist
}

proc ::Molefacture::hrealloc {seltext} {
  ### keep information about the original
  ### molecule and selection in case of
  ### editing an existing molecule
  variable origmolid
  variable origseltext
  variable tmpmolid
  variable protlasthyd
  variable picklist

  set prevmol $origmolid
  set prevsel $origseltext
  set curtmpmolid $tmpmolid
  set prevprotlasthyd $protlasthyd
  set prevpicklist $picklist

  ::Molefacture::molefacture_gui_aux $seltext 0

  set protlasthyd $prevprotlasthyd
  mol delete $curtmpmolid
  set origseltext $prevsel
  set origmolid $prevmol
  set picklist $prevpicklist
}

proc ::Molefacture::set_phipsi {newphi newpsi} {
  variable phi
  variable psi
  set phi $newphi
  set psi $newpsi
}

proc ::Molefacture::build_textseq {sequence} {
  set seqlist [split $sequence ""]
  foreach aa $seqlist {
    switch -nocase $aa {
      A {::Molefacture::add_aa ALA}
      C {::Molefacture::add_aa CYS}
      D {::Molefacture::add_aa ASP}
      E {::Molefacture::add_aa GLU}
      F {::Molefacture::add_aa PHE}
      G {::Molefacture::add_aa GLY}
      H {::Molefacture::add_aa HIS}
      I {::Molefacture::add_aa ILE}
      K {::Molefacture::add_aa LYS}
      L {::Molefacture::add_aa LEU}
      M {::Molefacture::add_aa MET}
      N {::Molefacture::add_aa ASN}
      P {::Molefacture::add_aa PRO}
      Q {::Molefacture::add_aa GLN}
      R {::Molefacture::add_aa ARG}
      S {::Molefacture::add_aa SER}
      T {::Molefacture::add_aa THR}
      V {::Molefacture::add_aa VAL}
      W {::Molefacture::add_aa TRP}
      Y {::Molefacture::add_aa TYR}
      default {error "Warning: I didn't recognize the amino acid $aa. Please use one letter code and the standard 20 amino acids. I'm skipping this entry."}
    }
  }
}

proc ::Molefacture::search_atomlist {mylist myindex} {
  #Search a list of lists for an element where the first field matches myindex
  set index 0
  foreach elem $mylist {
#    puts "Comparing $myindex with [lindex $elem 0]"
    if {[lindex $elem 0] == $myindex} {
      return $index
    }
    incr index
  }

  return -1
}

proc ::Molefacture::havecommonelems {list1 list2 {excludelist {}}} {
# Check to see if list1 and list2 have any common elements
# Return 1 if they do, -1 if they don't
# Common elements that are in excludelist will be ignored
# Yes yes, I know this isn't the most efficient algorithm for this operation

  foreach elem1 $list1 {
    if {[lsearch -exact $list2 $elem1] >= 0 && [lsearch -exact $excludelist $elem1] == -1} {
      return 1
    }
  }

  return -1
}


########################################################################
# Merge the current editable molecule with the original molecule       #
# Using the "mol fromsels" command                                     #
########################################################################
proc ::Molefacture::apply_changes_to_parent_mol {} {


  variable origmolid
  variable tmpmolid
  variable origseltext
  variable capAtoms
  variable toCapAtoms
  variable oricFlagCap
  variable origFlagReg
  variable origFlagSel
  variable newFlagModCap
  variable unknownRes

  if {$origseltext == ""} {
    tk_messageBox -icon warning -type ok -title "No parent molecule exists" -parent $Molefacture::topGui -message "Since molefacture was started with no selection, there's no parent molecule to merge changes into."
    return
  }

  mol reanalyze $tmpmolid

  set caps 0

  set selall [atomselect $origmolid "all"]


  $selall set $oricFlagCap 0
  $selall set $origFlagReg 1
  $selall set $origFlagSel 0
  $selall set $newFlagModCap 0
  $selall set beta 0
  $selall delete

  ### use  $origFlagReg == 1 to mark the modified region
  set notmysel [atomselect $origmolid "not ($origseltext)"]
  $notmysel set $origFlagReg 0

  if { [llength $capAtoms] > 0} {
    set notmyselCAP [atomselect $origmolid "index [join $capAtoms]"]
    #### use the flag0 as cap atoms on the original molecule
    $notmyselCAP set $oricFlagCap 1

    ### use flags when the values get ported in the mol fromsels command
    $notmyselCAP set beta -1
    set caps 1
  }

  set selall [atomselect $tmpmolid "all"]
  set modFlagCap flag3
  set modFlagSel flag1
  $selall set $modFlagCap 0
  $selall set $modFlagSel 0
  $selall delete
  set modifiedsel [atomselect $tmpmolid "$Molefacture::notDumSel and not name HC1"]


  if { [llength $toCapAtoms] > 0} {
    set modifiedseltopCAP [atomselect $tmpmolid "index [join $toCapAtoms]"]
    #### use the flag3 as cap atoms on the modified molecule
    $modifiedseltopCAP set $modFlagCap 1
    ### use flags when the values get ported in the mol fromsels command
    $modifiedseltopCAP set beta 2
    set caps 1
  }


  if {$Molefacture::qwikmd == 1 } {
    ### load the molecule to the "Edit Atom" window of QwikMD

    set topfile -1
    if {$Molefacture::topoinfo == "" || $Molefacture::topocgenff == 0} {
      set resname [lindex $unknownRes 0]
      set topfile [Molefacture::submitCGENFF $resname ]
      if {$topfile == -1} {
        ### Generate the topology with the information available in VMD
        set filename [pwd]/[lindex $unknownRes 0].top
        write_topfile $filename "$Molefacture::notDumSel"
        if {[file exists $filename] == 1} {
          if {[file dirname $filename] == "."} {
            set filename "[pwd]/$filename"
          }
          set topfile $filename

        } else {
          tk_messageBox -title "Topology Error" -message "Error Generating the topology based on the information inserted in VMD." -parente $Molefacture::topGui -icon error -type ok
          return
        }
      }
    }

    global env
    set pwd [pwd]

    cd $env(TMPDIR)
    if {[file exist molefacture_tmp] != 1} {
      file mkdir molefacture_tmp
    }
    cd molefacture_tmp
    set filename [lindex $Molefacture::topoinfo 0 6]
    Molefacture::exportFile $modifiedsel $filename.mol2 mol2
    cd $pwd
    file copy -force ${env(TMPDIR)}/molefacture_tmp/$filename.mol2 ${env(TMPDIR)}/$filename.mol2


    if {$Molefacture::topoinfo != "" && $Molefacture::topocgenff == 1} {

      file copy -force [lindex $Molefacture::topoinfo 0 0] ${env(TMPDIR)}/
      set topfile ${env(TMPDIR)}/[file tail [lindex $Molefacture::topoinfo 0 0]]

    }


    QWIKMD::loadMol2EditAtom ${topfile} ${env(TMPDIR)}/$filename.mol2

    cd $pwd

  } else {

    set selList [list $notmysel $modifiedsel]

    ### Make sure that the mol fromsels get the selection of the modified fragments
    ### in the correct the order to assemble the molecule in the original sequence
    if {$caps == 1} {
      set selList  [list]
      set mysel [atomselect $origmolid "$origFlagReg == 1"]

      set origFirst 0
      set targetFirst 0

      set length [llength $capAtoms]
      if {[llength $toCapAtoms] < [llength $capAtoms]} {
        set length [llength $toCapAtoms]
      }

      set opt -1
      for {set i 0} { $i < $length} {incr i 2} {
        set orgiCap [lindex $capAtoms $i]

        set targetCap [lindex $toCapAtoms [expr $i +1] ]

        ### evaluate if the modified section is in the middle of the chain (opt == 0)
        ### or at the begging (opt == 1)

        if {$i == 0} {
           if {$orgiCap < [lindex [$mysel get index] 0]} {
            set opt 0
           } else {
            set opt 1
           }
        }

        if {$opt == 0} {


          if {$i > [expr [llength $capAtoms] -2]} {
            set orgiCap [lindex [$notmysel get index] end]
          }

          set targetFirst [lindex $toCapAtoms $i]

          if {$i > [expr [llength $toCapAtoms] -2]} {
            set targetCap [lindex [$modifiedsel get index] end]
          }

          ### The origFlagSel is used to keep track which atoms where already selected

          lappend selList [atomselect $origmolid "(same residue as (index $origFirst to $orgiCap)) and ($origFlagReg == 0) and $origFlagSel == 0"]
          [lindex $selList end] set $origFlagSel 1

          if {$targetFirst != "" && $targetCap != ""} {
            lappend selList [atomselect $tmpmolid "same residue as (index $targetFirst to $targetCap)"]
            [lindex $selList end] set $modFlagSel 1
          }

          set origFirst [lindex $capAtoms [expr $i +1]]

          if {$i == [expr [llength $capAtoms] -2] } {
            lappend selList [atomselect $origmolid "(same residue as (index [lindex $capAtoms end] to [lindex [$notmysel get index] end])) and ($origFlagReg == 0) and $origFlagSel == 0"]
            [lindex $selList end] set $origFlagSel 1
          }

          if {$i == [expr [llength $toCapAtoms] -2]} {
            lappend selList [atomselect $tmpmolid "same residue as ($Molefacture::notDumSel and not name HC1 and $modFlagSel == 0)"]
            [lindex $selList end] set $modFlagSel 1
          }


        } else {

          if {$i > [expr [llength $capAtoms] -2]} {
            set orgiCap [lindex [$notmysel get index] end]
          }

          set origFirst [lindex $capAtoms $i]

          if {$i > [expr [llength $toCapAtoms] -2]} {
            set targetCap [lindex [$modifiedsel get index] end]
          }


          lappend selList [atomselect $tmpmolid "same residue as (index $targetFirst to $targetCap) and $modFlagSel == 0"]
          [lindex $selList end] set $modFlagSel 1



          set origFirst [lindex $capAtoms $i]
          set orgiCap [lindex $capAtoms  [expr $i +1]]

          if {$i > [expr [llength $capAtoms] -2]} {
            set orgiCap [lindex [$notmysel get index] end]
          }

          if {$origFirst != "" && $orgiCap != ""} {
            lappend selList [atomselect $origmolid "(same residue as (index $origFirst to $orgiCap)) and ($origFlagReg == 0) and $origFlagSel == 0"]
            [lindex $selList end] set $origFlagSel 1
          }

          set targetFirst [lindex $toCapAtoms [expr $i +1]]

          if {$i == [expr [llength $toCapAtoms] -2]} {
            lappend selList [atomselect $tmpmolid "same residue as ($Molefacture::notDumSel and not name HC1 and $modFlagSel == 0)"]
            [lindex $selList end] set $modFlagSel 1
          }

          if {$i == [expr [llength $capAtoms] -2] } {
            lappend selList [atomselect $origmolid "(same residue as (index [lindex $capAtoms end] to [lindex [$notmysel get index] end])) and ($origFlagReg == 0) and $origFlagSel == 0"]
            [lindex $selList end] set $origFlagSel 1
          }

        }
      }
      $mysel delete
      $modifiedseltopCAP delete
      $notmyselCAP delete
    }

    set newmol [mol fromsels $selList]

    $modifiedsel delete
    $notmysel delete

    set toCapsel [atomselect $newmol "beta == -1"]


    set Capsel [atomselect $newmol "beta == 2"]

    set capBondList ""
    if {[$Capsel num] > 0} {

      foreach tocapIndex [$toCapsel get index] capIndex [$Capsel get index] {
        if {[llength $tocapIndex] == 0 || [llength $capIndex] == 0 } {
          continue
        }


        set atom1sel [atomselect $newmol "index $capIndex"]
        set atom1bonds [lindex [$atom1sel getbonds] 0]

        lappend atom1bonds $tocapIndex
        $atom1sel setbonds [list $atom1bonds]




        set atom2sel [atomselect $newmol "index $tocapIndex"]
        set atom2bonds [lindex [$atom2sel getbonds] 0]
        lappend atom2bonds $capIndex


        $atom2sel setbonds [list $atom2bonds]

        #### change the segname to ensure the bond regions are part of the same fragment
        set segname [lindex [$atom1sel get segname] 0]
        if {[llength $segname] > 0} {

          set fragment [atomselect $newmol "segname $segname"]
          set origsegname [lindex [$atom2sel get segname] 0]
          if {[llength $origsegname ] == 0} {
            set origsegname "\{\}"
          }
          $fragment set segname ${origsegname}

          $fragment delete
        }
        $atom2sel delete
        $atom1sel delete
      }
    }
    $Capsel delete
    $toCapsel delete

    # mol reanalyze $newmol

    set molname [molinfo $origmolid get name]

    set newname [file root $molname]_edited
    if {[file extension $molname] != ""} {
      append newname [file extension $molname]
    }

    set sel [atomselect $newmol "all"]
    if {[lindex [$sel get atomicnumber] 0] == -1} {assign_elements}
    $sel delete

    mol reanalyze $newmol

    mol rename $newmol $newname

    ::CloneRep::clone_reps $origmolid $newmol
    set repindex [molinfo $newmol get numreps]
    mol addrep $newmol
    mol modselect $repindex $newmol "$origseltext"
    mol modstyle $repindex $newmol "Licorice"
    mol modcolor $repindex $newmol Element

  }

    Molefacture::closeMainWindow
    Molefacture::initialize

    
    if {[molinfo num] != 0} {
      mol reanalyze top
      display resetview
    }

}

proc ::Molefacture::do_fep_typing {} {
# Split the current FEP molecule into a before and after molecule, run the
# desired typing scheme on both of them, and then write a corresponding hybrid
# topology
# remember that beta=-1 is outgoing and +1 is incoming

  variable tmpmolid
  variable feptyping
  variable chargelist

  variable fepstartcharge
  variable fependcharge

  variable FEPdelres
  variable FEPreplaceflag

  set myresid 1
  if {$FEPreplaceflag > 0} {
    set myresid $FEPdelres
  }

  unique_atomnames

  set outgoingatoms [atomselect $tmpmolid "occupancy >= 0.5 and beta <= 0"]
  set outbeta [$outgoingatoms get beta]
  set outocc [$outgoingatoms get occupancy]
  $outgoingatoms writexbgf Molefacture_fep_tmp.xbgf
  set outgoingmol [mol new Molefacture_fep_tmp.xbgf]
  $outgoingatoms delete
  set outgoingatoms [atomselect $outgoingmol all]
  $outgoingatoms set beta $outbeta
  $outgoingatoms set occupancy $outocc

  set incomingatoms [atomselect $tmpmolid "occupancy >= 0.5 and beta >= 0"]
  set inbeta [$incomingatoms get beta]
  set inocc [$incomingatoms get occupancy]
  $incomingatoms writexbgf Molefacture_fep_tmp.xbgf
  set incomingmol [mol new Molefacture_fep_tmp.xbgf]
  $incomingatoms delete
  set incomingatoms [atomselect $incomingmol all]
  $incomingatoms set beta $inbeta
  $incomingatoms set occupancy $inocc

  # run typing and load the optimized geometries
  switch -exact $feptyping {
    GAFF {
      ::ANTECHAMBER::ac_type_in_place $outgoingatoms bcc $fepstartcharge
      file copy -force sqm.xyz outgoing.xyz
      load_sqm_coords $outgoingatoms
      ::ANTECHAMBER::ac_type_in_place $incomingatoms bcc $fependcharge
      load_sqm_coords $incomingatoms
      file copy -force sqm.xyz incoming.xyz
    }
    OPLS {
      puts [glob *xyz]
      puts "Running typing of outgoing atoms"
      ::ANTECHAMBER::ac_type_in_place $outgoingatoms cm1 $fepstartcharge "CUSTOM[file join $::env(MOLEFACTUREDIR) lib ATOMTYPE_OPLS.DEF]"
      #::ANTECHAMBER::ac_type_in_place $outgoingatoms bcc $fepstartcharge "CUSTOM[file join $::env(MOLEFACTUREDIR) lib ATOMTYPE_OPLS.DEF]"
      file copy -force sqm.xyz outgoing_opls.xyz
      load_divcon_coords $outgoingatoms
      ::ANTECHAMBER::ac_type_in_place $incomingatoms cm1 $fependcharge "CUSTOM[file join $::env(MOLEFACTUREDIR) lib ATOMTYPE_OPLS.DEF]"
      file copy -force sqm.xyz incoming_opls.xyz
      load_divcon_coords $incomingatoms
    }
  }

  set commonatoms [atomselect $tmpmolid "beta < 0.5 and beta > -0.5 and occupancy >= 0.5"]
  set commonatomnames [$commonatoms get name]
  $commonatoms delete

  # Handle charges as follows:
  # First make all charges in the common portion the average of their
  # charge in the initial and final states
  # then scale the charges of the transformed atoms to give appropriate
  # final charges
  # warn the user if the charge RMSD of any of these changes > 0.05

  set chargermsdthreshold 0.05

  set origcommon [atomselect $tmpmolid "beta < 0.5 and beta > -0.5 and occupancy >= 0.5"]
  set oldcommon [atomselect $outgoingmol "beta < 0.5 and beta > -0.5 and occupancy >= 0.5"]
  set newcommon [atomselect $incomingmol "beta < 0.5 and beta > -0.5 and occupancy >= 0.5"]

  $incomingatoms move [measure fit $newcommon $oldcommon]


  set oldcharges [$oldcommon get charge]
  set newcharges [$newcommon get charge]

  set avgcharges [vecscale 0.5 [vecadd $oldcharges $newcharges] ]
  set olddiff [vecsub $avgcharges $oldcharges]
  set newdiff [vecsub $avgcharges $newcharges]
  set constantcharge [vecsum $avgcharges]

  set oldrmsd [expr sqrt( [vecdot $olddiff $olddiff] )]
  set newrmsd [expr sqrt( [vecdot $newdiff $newdiff] )]

  $oldcommon set charge $avgcharges
  $newcommon set charge $avgcharges

  if {$oldrmsd > $chargermsdthreshold || $newrmsd > $chargermsdthreshold} {
    tk_messageBox -icon warning -type ok -title Warning -message "Warning: RMSD of charges in the constant atom fragment to \
    the average of the incoming and outgoing forms is greater than $chargermsdthreshold ($oldrmsd old, $newrmsd new). You may \
    need to rethink your definition of incoming/outgoing atoms."
  }

  set incomingtotcharge [expr $fependcharge - $constantcharge]
  set outgoingtotcharge [expr $fepstartcharge - $constantcharge]

  set incomingonly [atomselect $incomingmol "beta > 0.5"]
  set outgoingonly [atomselect $outgoingmol "beta < -0.5"]


  set incomingcurrcharge [vecsum [$incomingonly get charge]]
  set outgoingcurrcharge [vecsum [$outgoingonly get charge]]


  set incomingchargedel [expr ($incomingtotcharge - $incomingcurrcharge) / [$incomingonly num]]
  set outgoingchargedel [expr ($outgoingtotcharge - $outgoingcurrcharge) / [$outgoingonly num]]

  set newincomingcharges [list]
  set newoutgoingcharges [list]

  foreach charge [$incomingonly get charge] {
    lappend newincomingcharges [expr $charge + $incomingchargedel]
  }

  foreach charge [$outgoingonly get charge] {
    lappend newoutgoingcharges [expr $charge + $outgoingchargedel]
  }

  set incdiff [vecsub $newincomingcharges [$incomingonly get charge]]
  set outdiff [vecsub $newoutgoingcharges [$outgoingonly get charge]]

  set incrmsd [expr sqrt( [vecdot $incdiff $incdiff])]
  set outrmsd [expr sqrt( [vecdot $outdiff $outdiff])]

  $incomingonly set charge $newincomingcharges
  $outgoingonly set charge $newoutgoingcharges

  if {$incrmsd > $chargermsdthreshold || $outrmsd > $chargermsdthreshold} {
    tk_messageBox -icon warning -type ok -title Warning -message "Warning: RMSD of charges in the incoming or outgoing atom \
    fragments to the exact charges from charge assignment is greater than $chargermsdthreshold ($incrmsd incoming, $outrmsd outgoing). \
    You may need to rethink your definition of incoming/outgoing atoms."
  }

  # Make sure all atoms have unique names
  foreach sel [list $outgoingonly $incomingonly] suffix {A B} {
    set newnames [list]
    foreach name [$sel get name] {
      lappend newnames "${name}${suffix}"
    }
    $sel set name $newnames
  }


  write_fep_hybrid_topfile "Molefacture_fep_temp.top" $oldcommon $incomingatoms $outgoingonly $incomingonly [expr $fepstartcharge + $fependcharge]

  # Make one pdb with the before and after atoms
  $oldcommon set beta 0
  $outgoingonly set beta -1
  $incomingonly set beta 1
  $oldcommon set resid $myresid
  $outgoingonly set resid $myresid
  $incomingonly set resid $myresid
  $oldcommon writepdb "Molefacture_fep_temp_A.pdb"
  $outgoingonly writepdb "Molefacture_fep_temp_B.pdb"
  $incomingonly writepdb "Molefacture_fep_temp_C.pdb"

  puts "writing combined file"
  set ofile [open "Molefacture_fep_temp_combined.pdb" w]

  foreach part {A B C} {
    set infile [open "Molefacture_fep_temp_${part}.pdb" "r"]
    while { [gets $infile line] >= 0 } {
      set match [regexp ATOM $line matchstr]
      if { $match } {
        puts $ofile "$line"
      }
    }
    close $infile
  }

  puts $ofile "END"
  close $ofile

  # Update the coordinates in the original atom to match the optimized geometries
  set outgoingatoms_origmol [atomselect $tmpmolid "occupancy >= 0.5 and beta < 0"]
  set incomingatoms_origmol [atomselect $tmpmolid "occupancy >= 0.5 and beta > 0"]

  $outgoingatoms_origmol set {x y z} [$outgoingonly get {x y z}]
  $incomingatoms_origmol set {x y z} [$incomingonly get {x y z}]
  $origcommon set {x y z} [$oldcommon get {x y z}]

  $oldcommon delete
  $newcommon delete
  $outgoingonly delete
  $incomingonly delete
  $outgoingatoms_origmol delete
  $incomingatoms_origmol delete
  $origcommon delete

}

proc load_mopac_coords {selection} {
# load the coordinates from mopac.pdb into selection
  set mopacmol [mol new mopac.pdb]
  set mopacsel [atomselect top all]
  $selection set {x y z} [$mopacsel get {x y z}]
  $mopacsel delete
  mol delete $mopacmol
}

proc load_divcon_coords {selection} {
# load the coordinates from divcon.pdb into selection
  set mopacmol [mol new divcon.pdb]
  set mopacsel [atomselect top all]
  $selection set {x y z} [$mopacsel get {x y z}]
  $mopacsel delete
  mol delete $mopacmol
}

proc load_sqm_coords {selection} {
# load the coordinates from sqm.xyz into selection
  set sqmmol [mol new sqm.xyz]
  set sqmsel [atomselect top all]
  $selection set {x y z} [$sqmsel get {x y z}]
  $sqmsel delete
  mol delete $sqmmol
}

proc ::Molefacture::lremove {listVariable value} {
    upvar 1 $listVariable var
    set idx [lsearch -exact $var $value]
    set var [lreplace $var $idx $idx]
}

# proc ::Molefacture::parsePSFfile { psffile prmfile newprmfile } {
#     set fp [open $psffile r]
#     set file_data [read $fp]
#     close $fp
#     set data [split $file_data "\n"]
#     set cursection -1
#     set prevsection -1

#     set atoms {}
#     set atomtypes {}
#     set bonds {}
#     set angles {}
#     set dihedrals {}
#     set impropers {}
#     set natoms 0
#     set nbonds 0
#     set nangles 0
#     set ndihedrals 0
#     set nimpropers 0
#     puts "Parsing PSF file $psffile"
#     foreach line $data {
#         if {[string length $line] > 1} {
#         if { [string first "NATOM" $line] != -1 } {
#             set cursection 0
#         } elseif { [string first "NTHETA" $line] != -1 } {
#             set cursection 2
#         } elseif {[string first "NBOND" $line] != -1 } {
#             set cursection 1
#         } elseif {[string first "NPHI" $line] != -1 } {
#             set cursection 3
#         } elseif {[string first "NIMPHI" $line] != -1 } {
#             set cursection 4
#         } elseif {[string first "NDON" $line] != -1} {
#             set cursection 5
#         }

#         if { $cursection == 0 && $prevsection == 0 } {
#             lappend atoms [string trim [string range $line 24 27]]
#             lappend atomtypes [string trim [string range $line 29 32]]
# #            puts "$natoms [lindex $atoms $natoms] [lindex $atomtypes $natoms]"
#             incr natoms
#         } elseif { $cursection == 1 && $prevsection == 1 } {
#             set numbonds [expr [string length $line]/16]
#             for { set i 0 } { $i < $numbonds } {incr i } {
#                 set b0 [string range $line [expr $i*16] [expr $i*16 + 7]]
#                 set b1 [string range $line [expr $i*16+8] [expr $i*16 + 15]]
#                 set b {}
#                 lappend b [expr [string trim $b0] -1 ]
#                 lappend b [expr [string trim $b1] -1 ]
#                 lappend bonds $b
#                 incr nbonds
#             }
#         } elseif { $cursection == 2 && $prevsection == 2 } {
#             set numbonds [expr [string length $line]/24]
#             for { set i 0 } { $i < $numbonds } {incr i } {
#                 set b0 [string range $line [expr $i*24] [expr $i*24 + 7]]
#                 set b1 [string range $line [expr $i*24+8] [expr $i*24 + 15]]
#                 set b2 [string range $line [expr $i*24+16] [expr $i*24 + 23]]
#                 set b {}
#                 lappend b [expr [string trim $b0] -1 ]
#                 lappend b [expr [string trim $b1] -1 ]
#                 lappend b [expr [string trim $b2] -1 ]
#                 lappend angles $b
#                 incr nangles
#             }
#         } elseif { $cursection == 3 && $prevsection == 3 } {
#             set numbonds [expr [string length $line]/32]
#             for { set i 0 } { $i < $numbonds } {incr i } {
#                 set b0 [string range $line [expr $i*32] [expr $i*32 + 7]]
#                 set b1 [string range $line [expr $i*32+8] [expr $i*32 + 15]]
#                 set b2 [string range $line [expr $i*32+16] [expr $i*32 + 23]]
#                 set b3 [string range $line [expr $i*32+24] [expr $i*32 + 31]]
#                 set b {}
#                 lappend b [expr [string trim $b0]  -1 ]
#                 lappend b [expr [string trim $b1] -1 ]
#                 lappend b [expr [string trim $b2] -1 ]
#                 lappend b [expr [string trim $b3] -1 ]
#                 lappend dihedrals $b
#                 incr ndihedrals
#             }
#         } elseif { $cursection == 4 && $prevsection == 4 } {
#             set numbonds [expr [string length $line]/32]
#             for { set i 0 } { $i < $numbonds } {incr i } {
#                 set b0 [string range $line [expr $i*32] [expr $i*32 + 7]]
#                 set b1 [string range $line [expr $i*32+8] [expr $i*32 + 15]]
#                 set b2 [string range $line [expr $i*32+16] [expr $i*32 + 23]]
#                 set b3 [string range $line [expr $i*32+24] [expr $i*32 + 31]]
#                 set b {}
#                 lappend b [expr [string trim $b0] -1 ]
#                 lappend b [expr [string trim $b1] -1 ]
#                 lappend b [expr [string trim $b2] -1 ]
#                 lappend b [expr [string trim $b3] -1 ]
#                 lappend impropers $b
#                 incr nimpropers
#             }
#        }
#         set prevsection $cursection
#     }
#     }
#     #set bonds [$atoms getbonds]
#     set n 0
#     set bondlist {}
#     foreach b $bonds {
#         set n [lindex $b 0]
#         set b2 [lindex $b 1]
#         set cur [format "%4s %4s" [lindex $atomtypes $n] [lindex $atomtypes $b2]]
#         if { [lsearch -exact $bondlist $cur] == -1 } {
#             lappend bondlist $cur
#             #puts $cur
#         }
#     }
#     #puts "ANGLE"
#     set anglelist {}
#     foreach ang $angles {
#         set cur [format "%4s %4s %4s" [lindex $atomtypes [lindex $ang 0 ]] [lindex $atomtypes [lindex $ang 1]] [lindex $atomtypes [lindex $ang 2]]]
#         if { [lsearch -exact $anglelist $cur] == -1 } {
#            lappend anglelist $cur
#            #puts "$cur"
#         }
#     }
#     set dihedrallist {}
#     #puts "DIHEDRAL"
#     foreach dih $dihedrals {
#         set cur [format "%4s %4s %4s %4s" [lindex $atomtypes [lindex $dih 0 ]] [lindex $atomtypes [lindex $dih 1]] [lindex $atomtypes [lindex $dih 2]] [lindex $atomtypes [lindex $dih 3]]]
#         if { [lsearch -exact $dihedrallist $cur] == -1 } {
#            lappend dihedrallist $cur
#            #puts $cur
#         }
#     }
#     # NOTE IMPROPERS FORMATTED SUCH THAT 1st ATOM IS CENTRAL ATOM
#     # NAMD DOESN'T SEARCH FOR ATOM TYPES!!!
#     # set dihedrals [lindex [molinfo $mol get impropers] 0]
#     #puts "IMPROPER"
#     set improperlist {}
#     foreach dih $impropers {
#         set cur [format "%4s %4s %4s %4s" [lindex $atomtypes [lindex $dih 0 ]] [lindex $atomtypes [lindex $dih 1]] [lindex $atomtypes [lindex $dih 2]] [lindex $atomtypes [lindex $dih 3]]]
#         if { [lsearch -exact $improperlist $cur] == -1 } {
#            lappend improperlist $cur
#            #puts $cur
#         }
#     }
#     #puts "ATOMS"
#     set n 0
#     set atomlist {}
#     foreach b $atomtypes {
#         set cur [format "%4s" $b]
#         if { [lsearch -exact $bondlist $cur] == -1 } {
#             lappend atomlist $cur
#             #puts $cur
#         }
#     }

#     set fp [open $prmfile r]
#     set file_data [read $fp]
#     close $fp
#     set data [split $file_data "\n"]
#     set cursection -1
#     set prevsection -1
#     #puts "REQ'D PARAMETERS"
#     set fo [open $newprmfile w]
#     set ncnt -1
#     foreach line $data {
#         incr ncnt
#         if {[string index $line 0] != "!"} {
#         if { [string first "NONBONDED" $line] != -1 } {
#             set line [string trimright $line]
#             #puts "\n$line"
#             puts $fo "\n$line"
#             set m 1
#             while { [string range $line [expr [string length $line]-1] [expr [string length $line]-1]] == "-" } {
#                 #puts "[lindex $data [expr $ncnt+$m]]"
#                 puts $fo "[lindex $data [expr $ncnt+$m]]"
#                 set line "[lindex $data [expr $ncnt+$m]]"
#                 incr m
#             }
#             set cursection 0
#         } elseif { [string first "ANGLE" $line] != -1 } {
#             set cursection 2
#             #puts "ANGLE"
#             puts $fo "\n$line"

#         } elseif {[string first "DIHEDRAL" $line] != -1 } {
#             set cursection 3
#             #puts "DIHEDRAL"
#             puts $fo "\n$line"
#         } elseif {[string first "IMPROPER" $line] != -1 } {
#             set cursection 4
#             #puts "IMPROPER"
#             puts $fo "\n$line"
#         } elseif {[string first "BOND" $line] != -1 } {
#             set cursection 1
#             #puts "BOND"
#             puts $fo "\n$line"
#         }

#         if { $cursection == 0 && $prevsection == 0 } {
#             set a0 [string range $line 0 3]
#             set a0 [string map { "*" "\\*" } $a0]
#             set parms [string range $line 4 [string length $line]]
#             set key [format "%4s" [string trim $a0]]
#             #puts "? $key"
#             if { [lsearch -exact $atomlist $key] != -1 } {
#                     ::Molefacture::lremove atomlist $key
#                     set line "$a0$parms"
#                     puts $fo $line
#                     #puts $line
#            }
#        } elseif { $cursection == 1 && $prevsection == 1 } {
#             set a0 [string range $line 0 3]
#             set a1 [string range $line 5 8]
#             set parms [string range $line 9 [string length $line]]
#             set key [format "%4s %4s" $a0 $a1]
#             set keyrev [format "%4s %4s" $a1 $a0]
#             if { [lsearch -exact $bondlist $key] != -1 || [lsearch -exact $bondlist $keyrev] != -1 } {
#                 set n0 [lsearch -exact $bondlist $key]
#                 if { $n0  > -1 } {
#                     set key [lindex $bondlist $n0]
#                     set bondlist [lreplace $bondlist $n0 $n0]
#                     set line "$key$parms"
#                     puts $fo $line
#                 }
#                 set n1 [lsearch -exact $bondlist $keyrev]
#                 if { $n1 > -1 } {
#                     set key [lindex $bondlist $n1]
#                     set bondlist [lreplace $bondlist $n1 $n1]
#                     set line "$keyrev$parms"
#                     puts $fo $line
#                 }

#             }
#         } elseif { $cursection == 2 && $prevsection == 2 } {
#             set a0 [string range $line 0 3]
#             set a1 [string range $line 5 8]
#             set a2 [string range $line 10 13]
#             set parms [string range $line 14 [string length $line]]
#             set key [format "%4s %4s %4s" $a0 $a1 $a2]
#             set keyrev [format "%4s %4s %4s" $a2 $a1 $a0]
#             if { [lsearch -exact $anglelist $key] != -1 || [lsearch -exact $anglelist $keyrev] != -1 } {
#                 set n0 [lsearch -exact $anglelist $key]
#                 if { $n0 > -1 } {
#                     set anglelist [lreplace $anglelist $n0 $n0]
#                     set line "$key$parms"
#                     puts $fo $line
#                     #puts $line
#                 }
#                 set n0 [lsearch -exact $anglelist $keyrev]
#                 if { $n0 > -1 } {
#                     set anglelist [lreplace $anglelist $n0 $n0]
#                     set line "$keyrev$parms"
#                     puts $fo $line
#                     #puts $line
#                 }
#             }
#         } elseif { $cursection == 3 && $prevsection == 3 } {
#             set a0 [string range $line 0 3]
#             set a0 [string map { "*" "\\*" } $a0]
#             set a1 [string range $line 5 8]
#             set a1 [string map { "*" "\\*" } $a1]
#             set a2 [string range $line 10 13]
#             set a2 [string map { "*" "\\*" } $a2]
#             set a3 [string range $line 15 18]
#             set a3 [string map { "*" "\\*" } $a3]
#             set parms [string range $line 19 [string length $line]]
#             if { [string first " X" $a0] != -1 } { set a0 "...." }
#             if { [string first " X" $a1] != -1 } { set a1 "...." }
#             if { [string first " X" $a2] != -1 } { set a2 "...." }
#             if { [string first " X" $a3] != -1 } { set a3 "...." }
#             set key {}
#             lappend key [format "%4s %4s %4s %4s" $a0 $a1 $a2 $a3]
#             lappend key [format "%4s %4s %4s %4s" $a3 $a2 $a1 $a0]
#             set found -1
#             set n 0
#             set foundlist {}
#             set line ""
#             foreach k $key {
#                 while { [lsearch -regexp $dihedrallist $k] != -1 } {
#                     set found [lsearch -regexp $dihedrallist $k]
#                     if { [lsearch -exact $foundlist [lindex $dihedrallist $found]] == -1 } {
#                         lappend foundlist [lindex $dihedrallist $found]
#                         append line [lindex $dihedrallist $found]
#                         append line $parms
#                         append line "\n"
#                     }
#                     ::Molefacture::lremove dihedrallist [lindex $dihedrallist $found]

#                 }
#                 incr n
#             }
#             if {$found > -1} {
#                #puts [string trimright $line]
#                puts $fo [string trimright $line]
#             }
#         } elseif { $cursection == 4 && $prevsection == 4 } {
#             set a0 [string range $line 0 3]
#             set a0 [string map { "*" "\\*" } $a0]
#             set a1 [string range $line 5 8]
#             set a1 [string map { "*" "\\*" } $a1]
#             set a2 [string range $line 10 13]
#             set a2 [string map { "*" "\\*" } $a2]
#             set a3 [string range $line 15 18]
#             set a3 [string map { "*" "\\*" } $a3]
#             set parms [string range $line 19 [string length $line]]
#             if { [string first "X" $a0] != -1 } { set a0 "...." }
#             if { [string first "X" $a1] != -1 } { set a1 "...." }
#             if { [string first "X" $a2] != -1 } { set a2 "...." }
#             if { [string first "X" $a3] != -1 } { set a3 "...." }
#             set key {}
#             lappend key [format "%4s %4s %4s %4s" $a0 $a1 $a2 $a3]
#             lappend key [format "%4s %4s %4s %4s" $a0 $a1 $a3 $a2]
#             lappend key [format "%4s %4s %4s %4s" $a0 $a2 $a1 $a3]
#             lappend key [format "%4s %4s %4s %4s" $a0 $a2 $a3 $a1]
#             lappend key [format "%4s %4s %4s %4s" $a0 $a3 $a2 $a1]
#             lappend key [format "%4s %4s %4s %4s" $a0 $a3 $a1 $a2]
#             set found -1
#             set n 0
#             set foundlist {}
#             set line ""
#             foreach k $key {
#                 while { [lsearch -regexp $improperlist $k] != -1 } {
#                     set found [lsearch -regexp $improperlist $k]
#                     if { [lsearch -exact $foundlist [lindex $improperlist $found]] == -1 } {
#                         lappend foundlist [lindex $improperlist $found]
#                         append line [lindex $improperlist $found]
#                         append line $parms
#                         append line "\n"
#                     }
#                     ::Molefacture::lremove improperlist [lindex $improperlist $found]

#                 }
#                 incr n
#             }
#             if {$found > -1} {
#                #puts [string trimright $line]
#                puts $fo [string trimright $line]
#             }
#         }
#         set prevsection $cursection
#     }
#     }
#     puts $fo "\n\nEND\n"
#     close $fo
#     set missingparms ""
#     if { [llength $bondlist] > 0 } {
#         append missingparms "BONDS    [llength $bondlist] $bondlist\n"
#     }
#     if { [llength $anglelist] > 0 } {
#         append missingparms "ANGLES   [llength $anglelist] $anglelist\n"
#     }
#     if { [llength $dihedrallist] > 0 } {
#         append missingparms "DIHEDRAL [llength $dihedrallist] $dihedrallist\n"
#     }
#     if { [llength $improperlist] > 0 } {
#         append missingparms "IMPROPER [llength $improperlist] $improperlist\n"
#     }
#     #if { $missingparms != "" } {
#         #puts "Missing Parameters: "
#         #puts $missingparms
#     #}
#     puts "Writing PRM file $newprmfile"
#     return $missingparms
# }

proc ::Molefacture::write_psfpdbfiles { {outpref ""} } {
    variable tmpmolid
    variable unknownRes

    global env
    if {$outpref == ""} {
      set outpref [tk_getSaveFile -title "PSF & PDB File Name" -parent $Molefacture::topGui -filetypes {} -defaultextension ""]
    }
    set outpath [file dirname $outpref]
    set topfile ""
    set notcgenfftopo 0
    if {$outpref == ""} {
      return -1
    }
    set returnlist ""
    #### if the molecule was edited or the topology was never generated,
    #### fetch again from the CGenFF server
    set newtop ""
    if {$Molefacture::topocgenff == 0} {

      set found 0


      set sel [atomselect $tmpmolid "$Molefacture::notDumSel and not water and not ions"]
      set resnames [lsort -unique -dictionary [$sel get resname]]
      $sel delete
      set resultRes [Molefacture::loadDefaulTopo $resnames]
      set unknownRes [lindex $resultRes 0]
      set knownRes [lindex $resultRes 1]

      if {[llength $unknownRes] != 0} {
        
        set topfile [Molefacture::submitCGENFF $outpref]

        ### topfile == -1 in case the user decided not to use the CGenFF server
        if {$topfile == -1} {
          set notcgenfftopo 1
          if {$outpref != ""} {
            set filename $outpref.top
          } else {
            set filename [pwd]/[lindex $unknownRes 0].top
          }

          ### Generate the topology with the information available in VMD
          write_topfile $filename "$Molefacture::notDumSel"
          if {[file exists $filename] == 1} {
            set newtop $filename

            if {[file dirname $filename] == "."} {
              set filename "$filename"
            }
            set topfile $filename

          } else {
            return
          }
        } elseif {$topfile == -2} {
         return -1
        }
      }
      
      
    } else {
      set newtop [lindex $Molefacture::selTop end]
    }

    set topfile ""
    foreach top $Molefacture::topoinfo {
      lappend topfile [lindex ${top} 0]
    }

    if {[info exists Molefacture::deftopinfo] && [lsearch $topfile [lindex $Molefacture::deftopinfo end 0]] == -1} {
      set str ""
      foreach top $Molefacture::deftopinfo {
        lappend topfile [lindex ${top} 0]
      }
    }

    set sel [atomselect $tmpmolid "$Molefacture::notDumSel"]
    if { [ $sel num] > 0 } {

        set segname [lsort -unique -dictionary [$sel get segname]]

        set topolist [list]
        foreach topo $topfile {
            # if {[file exists $outpath/[file tail $topo] ] != 1} {
              set file $topo
              file copy -force -- ${file} ${outpath}/
            # }
            lappend topolist [file tail $topo]
        }


        cd $outpath/
        set outpdbfile "$outpref.pdb"
        set psffile "$outpref.psf"

        set pdbfile ""
      
        
        if {$notcgenfftopo == 0} {
          
          foreach seg $segname {
            set selseg [atomselect $tmpmolid "segname $seg"]
            $selseg writepdb "$seg.pdb"
            lappend pdbfile "$seg.pdb"
            $selseg delete
          }
          psfcontext reset
          resetpsf
          set molefaccontext [psfcontext create]
          psfcontext eval $molefaccontext {
            puts "DEBUG topolist $topolist"
            foreach topo $topolist {
              topology [lindex $topo 0]
            }
            if {[lsearch $topolist "top_all36_cgenff.rtf"] == -1} {
               topology $env(CHARMMTOPDIR)/top_all36_cgenff.rtf
            }
            foreach seg $segname pdb $pdbfile {
                segment $seg {
                  pdb $pdb
                  first none
                  last none
                  auto angles dihedrals
                }
                coordpdb $pdb
            }
            writepsf $psffile
            writepdb $outpdbfile
         }
         foreach seg $segname {
          if {"$seg.pdb" != "$outpdbfile"} {
           file delete -force $seg.pdb
          }
         }
         psfcontext delete $molefaccontext
      } else {
        $sel writepdb $outpdbfile
        $sel writepsf $psffile
      }
       set returnlist [list $psffile $outpdbfile $newtop]
  }
  $sel delete
  return $returnlist
}

###################################################
# Fetch topology file from the CGenFF server      #
###################################################
proc ::Molefacture::submitCGENFF { outpref } {
  set topfile -1
  variable unknownRes
  global env

  if {$env(CGENFF) == "local" && $env(CGENFFPATH) == -1} {
    if {$Molefacture::cgenffcaller == 0} {
       return -1
    } else {
      tk_messageBox -message "Unparameterized residues found. \
      Please define the CGenFF local executable path, export the molecule as mol2 \
      and submit manually to the CGenFF Server, \
      or use UFF and OpenBabel to minimize the structure." \
      -icon error -type ok -parent $Molefacture::topGui
    }
  }

  if {[llength $unknownRes] > 1} {
     set answer [tk_messageBox -message "The molecule presents more than one unparameterized \
     residue - different residue names and/or IDs. This is not optimal to be submitted to\
      CGenFF. Do you wish to continue?" -parent $Molefacture::topGui -icon warning -type yesno -title "Missing Topology"]
      if {$answer == "no"} {
        return -2
      }
  }
  set message ""
  set his 0
  if {[llength $unknownRes] > 0} {
    if {[lsearch $unknownRes "HIS"] != -1} {
      set message "Histidine residues are define in CHARMM force field with different \
      residue names depending on the protonation state. Please change HIS residue(s) \
      name to:\n\n HSD (Neutral HIS, proton on ND1)\n HSE (Neutral His, proton on NE2) \
      or\n HSP (Protonated His).\n\nDo you wish to continue?"
      set his 1
    } else {
      set message "Residues [join $unknownRes]: molecules edited or missing topology"
    }

  } else {
    set message "Molecules edited or missing topology"
  }

  if {$his == 0} {
    set text "server"
    if {$env(CGENFF) == "local"} {
      set text "local executable"
    }
    set message "$message. Do you want to retrieve from the CGenFF ${text}?"
  }
  set answer [tk_messageBox -message $message -parent $Molefacture::topGui -title "Missing Topology" -icon info -type yesnocancel]
  if {$answer == "no"} {
    return -1
  } elseif {$answer == "cancel"} {
    return -2
  }

  set molfiles [Molefacture::run_cgenff export $outpref.mol2]
  if {$molfiles == 1} {
    return -1
  }

  ### env(CGENFF) can have two values "sever" and "local" and defines if 
  ### the CGenFF server or a local executable will be use to generate 
  ### the topology file

  if {$env(CGENFF) == "server" && $Molefacture::cgenffcaller} {

    if {[info exists env(CGENFFUSERNAME)] != 1 || $env(CGENFFUSERNAME) == -1} {
      Molefacture::cgenffUserWindows
      Molefacture::checkCgenffUser
      return -2
    }

    if {[eval "::CGENFFCALLER::callCGENFF -molfiles $molfiles -username $env(CGENFFUSERNAME) \
    -password $env(CGENFFPASSWORD)"] == 0} {

      set totalstr [::CGENFFCALLER::getNumStreamFiles]
      for {set i 0} {$i < $totalstr} {incr i} {
        set stream [::CGENFFCALLER::getStream $i]

        ### If CgenFF can't find any homologous topology, it returns "null" as
        ### message flag in the json object

        if {$stream != "null"} {
          set filename [file root [lindex $molfiles $i]]
          set file $filename.str
          
          set streamfile [open $file w+]
          foreach line [split $stream "\n"] {
             puts $streamfile $line
          }

          close $streamfile

          set continue [Molefacture::run_cgenff import [file root $filename].str]
          if {$continue == -2} {
            return -2
          }
        } else {
          tk_messageBox -message "$message. Please use Molefacture to \
          generate the topology." -title "CGenFF Topology Not Found" \
          -icon error -type ok -parent $Molefacture::topGui
         return -2
        }
        
      }
      set Molefacture::cgenffsub [::CGENFFCALLER::getCount]
    } else {
      tk_messageBox -title "CGenFF Server Not Available" -icon error \
      -message "Molefacture uses many mechanisms to access CGenFF server and \
      non of them worked. Please submit a export the structure to a mol2 file \
      (File->Write MOL2 file) and submit to the CGenFF server manually at \
      https://cgenff.umaryland.edu/." -type ok -parent $Molefacture::topGui
      return -2
    }
  } elseif {$env(CGENFF) == "local" && [file exist $env(CGENFFPATH)]} {
    
    ### Call the CGenFF local executable defined on the $env(CGENFFPATH)
    foreach mol $molfiles {
      set stream ""
      ### As far as I am aware, the CGenFF local executable is only available for
      ### Linux
      catch {exec $env(CGENFFPATH) $mol >& [file root $mol].str} 
      set filestream [open [file root $mol].str r]
      set stream [split [read $filestream] "\n"]
      close $filestream
      if {[lsearch $stream "RESI*"] == -1} {
        tk_messageBox -message "CGenFF local executable failed to \
          generate the topology." -title "CGenFF Topology Not Found" \
          -icon error -type ok -parent $Molefacture::topGui
          return -2
      } else {
       
       ### Add the character '*' to the print of the local CGenFF executable
       ### to make sure that is compliant with psfgen
       set filestream [open [file root $mol].str w+]
       set header 0
       foreach line $stream {
         if {[string index $line 0] == "*" } {
            set header 1
         }
         if {$header == 1} {
            puts $filestream $line
         } else {
            puts $filestream "*$line"
         }
       }
       close $filestream
      }
      set continue [Molefacture::run_cgenff import [file root $mol].str]
      if {$continue == -2} {
        return -2
      }
    }
  }

}
