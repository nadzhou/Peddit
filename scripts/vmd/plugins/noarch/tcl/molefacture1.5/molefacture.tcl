#
# Molefacture -- structure building plugin
#
# $Id: molefacture.tcl,v 1.109 2019/07/09 15:32:24 jribeiro Exp $
#

package require psfgen
package require idatm
package require runante
package require runsqm
package require readcharmmtop
package require utilities
package require exectool
package require topotools
package require clonerep
package require http
package require tktooltip
package require platform

if {$::tcl_version < 8.5} {
    package require dict
}

package provide molefacture 1.5


namespace eval ::Molefacture:: {
  proc format3Dec {val} {
    return [format %.3f [expr double(round(1000*$val))/1000]]
  }

  
   proc initialize {} {

      variable selectcolor lightsteelblue
      global env

      variable slavemode      0;  # Molefacture can be run as a slave for another application
      variable exportfilename "molefacture_save.xbgf"
      variable mastercallback {}; # When we are done editing in slavemode this callback of the master will be called

      variable origsel     {}
      variable origmolid "" ;# molid of parent molecule
      variable origseltext "" ;# Selection text used to start molefacture

      variable picklist     {}
      variable tmpmolid -1
      variable fepmolid -1; # molecule for fep graphics
      variable atommarktags {}
      variable ovalmarktags {}
      variable labelradius  0.2

      variable atmselrep     {}

      variable totalcharge  0
      variable atomlistformat {}

      variable anglelistformat {}

      variable dihedral   0
      variable dihedmarktags {}
      variable projectsaved 1

      variable showvale     0
      variable showellp     0

      variable angle 0
      variable anglemarktags {}

      variable minsteps 500
      variable minff "UFF"
      variable topocgenff 0


      # list of atoms needing markings in fep mode
      variable atomarktags {}

      variable oxidation {}

      variable lonepairs {}
      variable vmdindexes {}

      variable fragmentlist
      variable taglist
      variable templist
      variable atomlistformat
      variable bondtaglist
      variable bondlistformat

      variable protlasthyd -2; # make this value -1 a protein residue is added. Check if a protein is being edited or not
      variable nuclasthyd -1
      ### DEfine phi and psi to originate "straight" backbone conformations
      variable phi -180;# phi value for protein builder
      variable psi 180;# psi value for protein builder
      variable ralpha 180 ;# alpha angle for rna
      variable rbeta  180 ;# beta angle for rna
      variable rgamma 180 ;# gamma angle for rna
      variable rdelta 180 ;# delta angle for rna
      variable repsilon 180 ;# epsilon angle for rna
      variable rchi 180 ;# chi angle for rna
      variable rzeta 180 ;# zeta angle for rna
      variable nuctype "RNA" ;# RNA, DNA, or dsDNA
      variable dihedmoveatom "Atom1"

      variable addfrags ;# Array containing paths to fragments for H replacement, indexed by fragment name
      variable basefrags ;# Array containing paths to fragments used for making new molecules, indexed by basefrag name

     array set openvalence ""
     variable lastAtmSel ""
     variable elemcomb "Carbon(C:6)"
     variable mousemode "pick"
     variable prevmousemode "pick"
     variable cmboxstatvalues [list]
     variable cmboxstat ""

     variable bondguiValscale 0
     variable bondguiValentry 0

     variable angleguiValscale 0
     variable angleguiValentry 0
     ### variables controling the checkbuttons in the gui to
     ### select which group if atoms to move when changing the angle
     variable anglemvgrp1 0
     variable anglemvgrp2 0

     ### list of atoms associated with the group 1 and 2 of the checkbuttons
     ### in the gui
     variable isangle 0; # use this variable to avoid using the return command in isanglecommand command
     variable anglemvgrp1list 0
     variable anglemvgrp2list 0
     variable anglemother -1
     variable angleind1 -1
     variable angleind2 -1
     variable anglemothercoor ""
     variable angleind1coor ""
     variable angleind2coor ""
     variable curangle -1
     #### Dihedral angles
     variable dihedrgrp ""

     variable dihguiValscale 0
     variable dihguiValentry 0

     variable dihemother1 ""
     variable dihemother2 ""
     variable dihemother1coor ""
     variable dihemother2coor ""
     ##### Edit Molecules Variables
     variable resname ""
     variable chain ""
     variable resid ""

     variable segname ""

     ######### load topology variables
     variable topoinfo ""
     variable unknownRes [list]
     variable selTop [list]
     variable capAtoms [list]
     variable toCapAtoms [list]
     ######### Undo list information
     variable undo [list]
     ### tracks if first undo after save state
     variable firstUndo 0

     ### if called from QwikMD
     variable qwikmd 0

     variable cgenffWinVMDRC ""
     variable cgenffWinVMDRCPASS ""


     variable cyclNAtoms 6

     variable alchainNAtoms 6

     ### variable to be used as general info labels
     variable numAtoms 0
     variable charge 0.000

     variable addAtom 1
     variable addAtomSel ""
     
     ### check if the cgenffcaller plugin exists and set the variable to 1, or 0 
     ### if it does not exist.
     ### The command "package present" returns a string "package xxxx is not present"
     ### if the package was not loaded already.
     variable cgenffcaller -1
     set found ""
     ### Check if the CGenFF interface plugin
     ### is available 
      
     if {[catch {package require cgenffcaller}] == 1} {
        set cgenffcaller 0
     } else {
        set cgenffcaller 1
     }
     
     if {[info exist env(CGENFFUSERNAME)] != 1 || [info exist env(CGENFFPASSWORD)] != 1} {
        set env(CGENFFUSERNAME) -1
        set env(CGENFFPASSWORD) -1
     } 
     if {[info exist env(CGENFFPASSWORD)] != 1} {
        set env(CGENFFPASSWORD) -1
     }

     ### env(CGENFF) variable sets which CGenFF to call to get topology and parameter
     ### values can be server or local
     if {[info exist env(CGENFF)] != 1} {
        if {$cgenffcaller == 0} {
          set env(CGENFF) "local"
        } else {
          set env(CGENFF) "server"
        }
        
     } 

     if {[info exist env(CGENFFPATH)] != 1} {
        set env(CGENFFPATH) -1
     } 
     
   }

   ### Variable to use in radiobutton to select the type of template to build - cycle or chain
   variable cyclSelector cycle
   
    variable periodicNames {Unknown Hydrogen Helium Lithium Beryllium Boron Carbon Nitrogen Oxygen Fluorine Neon Sodium Magnesium Aluminum Silicon Phosphorus Sulfur Chlorine Argon Potassium Calcium Scandium Titanium Vanadium Chromium Manganese Iron Cobalt Nickel Copper Zinc Gallium Germanium Arsenic Selenium  Bromine Krypton Rubidium  Strontium Yttrium Zirconium Niobium Molybdenum  Technetium  Ruthenium Rhodium Palladium Silver  Cadmium Indium  Tin Antimony  Tellurium Iodine  Xenon Cesium  Barium  Lanthanum Cerium  Praseodymium  Neodymium Promethium  Samarium  Europium  Gadolinium  Terbium Dysprosium  Holmium Erbium  Thulium Ytterbium Lutetium  Hafnium Tantalum  Tungsten  Rhenium Osmium  Iridium Platinum  Gold  Mercury Thallium  Lead  Bismuth Polonium  Astatine  Radon}
      variable periodic {
        X H HE LI BE B C N O F NE NA MG AL SI P S CL AR K CA SC TI V CR MN FE CO NI CU ZN GA GE AS SE BR KR RB SR Y ZR NB MO TC RU RH PD AG CD IN SN SB TE I XE CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU HF TA W RE OS  IR PT AU HG TL PB BI PO AT RN
       }
      variable valence
      array set valence {
        X {0} H {1} HE {0}  LI {1} BE {2} B {3} C {4} N {3 4} O {2} F {1} NE {0} NA {1} MG {2} AL {3} SI {4}  P {3 5}  S {2 4 6} CL {1} AR {0}  K {1} CA {2} SC {0} TI {0} V {0} CR {0} MN {0} FE {2 3 4 6} CO {0} NI {0} CU {0} ZN {0}  GA {3} GE {4} AS {3} SE {2} BR {1} KR {0} RB {1} SR {2}  Y {0} ZR {0} NB {0} MO {0} TC {0} RU {0} RH {0} PD {0} AG {0} CD {0}  IN {3} SN {4} SB {3} TE {2}  I {1} XE {0}  CS {1} BA {2}  LA {0} CE {0} PR {0} ND {0} PM {0} SM {0} EU {0} GD {0} TB {0} DY {0} HO {0} ER {0} TM {0} YB {0}  LU {0} HF {0} TA {0} W {0} RE {0} OS {0} IR {0} PT {0} AU {0} HG {0}  TL {0} PB {0} BI {0} PO {0} AT {0} RN {0}
      }
      variable octet
      array set octet {X 0 H 2 HE 2 LI 8 BE 8 B 8 C 8 N 8 O 8 F 8 NE 8 NA 8 MG 8 AL 8 SI 8 P 10 S 8 CL 8 AR 8  K 18 CA 18 SC 18 TI 18 V 18 CR 18 MN 18 FE 18 CO 18 NI 18 CU 18 ZN 18 GA 18 GE 18 AS 18 SE 18 BR 18 KR 18 RB 18 SR 18 Y 18 ZR 18  NB 18 MO 18 TC 18 RU 18 RH 18 PD 18 AG 18 CD 18 IN 18 SN 18 SB 18  TE 18 I 18 XE 18 CS 32 BA 0 LA 0 CE 0 PR 0 ND 0 PM 0 SM 0 EU 0 GD 0 TB 0 DY 0 HO 0 ER 0 TM 0 YB 0 LU 0 HF 0 TA 0 W 0 RE 0 OS 0 IR 0 PT 0 AU 0 HG 0 TL 0 PB 0 BI 0 PO 0 AT 0 RN 0 }

      ### Bond distances from CHARMM FF
      ### To define generic bonds, like any Hydrogen bond, the index 0 must be HX. The X comes
      ### always at the end
      variable bondistlist {
        {CC 1.5} {CF 1.4} {CN 1.4} {NO 1.4} {CO 1.4} {OX 1.4} {CS 1.8}
      }
   variable mass_by_element
   array set mass_by_element { H 1.00794 He 4.002602 Li 6.941 Be 9.012182 B 10.811 C 12.011 N 14.00674 O 15.9994 F 18.9984032 Ne 20.1797 Na 22.989768 Mg 24.3050 Al 26.981539 Si 28.0855 P 30.973762 S 32.066 Cl 35.4527 K 39.0983 Ar 39.948 Ca 40.078 Sc 44.955910 Ti 47.88 V 50.9415 Cr 51.9961 Mn 54.93805 Fe 55.847 Ni 58.6934 Co 58.93320 Cu 63.546 Zn 65.39 Ga 69.723 Ge 72.61 As 74.92159 Se 78.96 Br 79.904 Kr 83.80 Rb 85.4678 Sr 87.62 Y 88.90585 Zr 91.224 Nb 92.90638 Mo 95.94 Tc 98 Ru 101.07 Rh 102.90550 Pd 106.42 Ag 107.8682 Cd 112.411 In 114.82 Sn 118.710 Sb 121.757 I 126.90447 Te 127.60 Xe 131.29 Cs 132.90543 Ba 137.327 La 138.9055 Ce 140.115 Pr 140.90765 Nd 144.24 Pm 145 Sm 150.36 Eu 151.965 Gd 157.25 Tb 158.92534 Dy 162.50 Ho 164.93032 Er 167.26 Tm 168.93421 Yb 173.04 Lu 174.967 Hf 178.49 Ta 180.9479 W 183.85 Re 186.207 Os 190.2 Ir 192.22 Pt 195.08 Au 196.96654 Hg 200.59 Tl 204.3833 Pb 207.2 Bi 208.98037 Po 209 At 210 Rn 222 Fr 223 Ra 226.025 Ac 227.028 Pa 231.0359 Th 232.0381 Np 237.048 U 238.0289 Am 243 Pu 244 Bk 247 Cm 247 Cf 251 Es 252 Fm 257 Md 258 No 259 Rf 261 Bh 262 Db 262 Lr 262 Sg 263 Hs 265 Mt 266}

   variable availabledummy 0

   variable notDumSel "not resname DUM and not segname CAP"
   variable dumSel "resname DUM or segname CAP"
   variable maxaasize 26 ;# size of the largest residue that could be added
   variable dummyincr 300 ;# number of dummy hydrogens added at each step

   variable feptyping "None" ;# typing scheme to use in the fep module
   variable fepstartcharge 0
   variable fependcharge 0
   variable filename_FEP
   variable topGui ".molefac"
   variable startGui ".molefacstart"
   variable periTabWindow ".molefacperitab"
   variable protBuildWindow ".protbuilder"
   variable cgenffWindow ".cgenff"
   variable periTab ""
   variable atomTable ""
   variable cmboxstatwidget ""
   variable bondscalewdgt ""

   ### Menu listing the default top
    variable defTopMenu ""

   variable anglescalewdgt ""
   variable angleentrywdgt ""
   variable dihscalewdgt ""
   variable dihentrywdgt ""
   variable rnamewdgt ""
   variable chnwdgt ""
   variable rsidwdgt ""
   variable chrgwdgt ""
   variable sgnmwdgt ""
   variable aapath [file join $::env(MOLEFACTUREDIR) lib amino_acids]
   variable nucpath [file join $::env(MOLEFACTUREDIR) lib nucleic_acids]
   variable anglemvgrp1wdgt ""
   variable anglemvgrp2wdgt ""
   variable bondentrywdgt ""
   variable dihradbtt1 ""
   variable dihradbtt2 ""
   variable bindTop 0
   # atom types that require impropers for OPLS ATOM TYPES only!!
   # ALL 3 bonded C & N types:
   variable improper_centeratoms [list C2 C1 N1 N2 N3 CA1 CA3 CA4 CA2 CM1 CM2 CM4 CM3 CT2 CT7 CT3 CT6 CT4 CT1 CT8  C* CN CB2 CV CW1 CX CR1 N21 NA CS1 CS2]
   variable improper_atomorders [list [list C O] [list N] [list CA] [list CM HC CM HC] [list CM HC CM CT] [list CM CT CM HC] [list CM CT CM CT]]

   ### Array to be used to store the default topologies entries in the menu and their checkbutton value (1 if selected, 0 if not)
  array set selDefTop ""

  ### default Topology list return from read_charmm_topology
  variable deftopinfo -1

  ### definition of the atom flags to mark operations performed on subsets of the structure


  ### use flag0 as cap atoms on the original molecule
  ### use flag3 as cap atoms on the new molecule from modified molecule
  ### use flag1 to mark the modified region
  ### use flag2 to mark already selected atoms to merge
  set oricFlagCap flag0
  set origFlagReg flag1
  set origFlagSel flag2
  set newFlagModCap flag3

  ### modFlagPSF marks the atoms that were submitted to autopsf/psfgen to
  ### add the hydrogens, set atom types, charges and etc. Useful when
  ### during the add hydrogen routine to avoid submit the same structure
  ### to autopsf (slow) all the time
  set modFlagPSF flag4

  ### flag to mark the atoms that H atoms were added
  set modFlagH flag5



  set arrowRight ""
  set arrowDown ""
  variable cyclSelecWidget [list]
  
  
  variable cgenffsub "---"; # number of submissions to the cgenff server
  
  
  variable typecgenff 0; # run cgenff local installation
  variable uernameentry; # path of the cgenff username entry
  variable cgenffvmdrc ""; # path of the frame of vmdrc definition
  variable cgenffpathentry ""; # path of the entry of cgenff path definition
  variable cgenffpathval ""; # value of the entry of cgenff path definition

  variable cgenffbrowser ""; # path of the button of cgenff path browser
  # initialize

  
}


#############################################################
## Get default distances between to atoms of the elements  ##
## ele1 and ele2                                           ##
#############################################################

proc ::Molefacture::getBondDist { ele1 ele2 } {
  variable bondistlist

  set dist 1.0

  if {$ele1 == "H" || $ele2 == "H"} {
    return $dist
  }
  set str1 "${ele1}$ele2"
  set str2 "${ele2}$ele1"

  set ind1 [lsearch -index 0 $bondistlist $str1]
  set ind2 [lsearch -index 0 $bondistlist $str2]
  set ind $ind1
  if {$ind1 == -1} {
    set ind $ind2
    if {$ind2 == -1} {
      set str1 "${ele1}X"
      set str2 "${ele2}X"
      set ind1 [lsearch -index 0 $bondistlist $str1]
      set ind2 [lsearch -index 0 $bondistlist $str2]
      set ind $ind1
      if {$ind1 == -1} {
        set ind $ind2
      }
    }
  }
  if {$ind != -1} {
    set dist [lindex [lindex $bondistlist $ind] 1]
  }


  return $dist
}

proc ::Molefacture::set_slavemode { callback filename } {
   variable slavemode 1
   variable exportfilename $filename
   variable mastercallback $callback
}

## Temporary gui for getting molefacture started
proc ::Molefacture::molefacture_start {} {
  variable topGui
  variable startGui

  if { [winfo exists $startGui] ==1} {
    wm deiconify $startGui
    focus $startGui
    raise $startGui
    return
  }


   set w [toplevel $startGui]
   wm title $w "Molefacture - Molecule Builder"
   wm resizable $w 1 0

   grid columnconfigure $startGui 0 -weight 1
   grid rowconfigure $startGui 0 -weight 1

   variable atomsel

   set atomsel ""
   grid [ttk::frame $startGui.mainframe] -row 0 -column 0 -sticky nwes -padx 4 -pady 4
   grid columnconfigure $startGui.mainframe 0 -weight 1
   grid rowconfigure $startGui.mainframe 0 -weight 1

   grid [ttk::frame $startGui.mainframe.lbl] -row 0 -column 0 -sticky nwes
   grid columnconfigure $startGui.mainframe.lbl 0 -weight 1
   grid rowconfigure $startGui.mainframe.lbl 0 -weight 1

  grid [ttk::label $startGui.mainframe.lbl.warning -text "Enter a selection below and click \"Start\" to start molefacture and edit the atoms of this selection. Please check the documentation (accessible through the Help menu) to learn how to use it." -wraplength 360 -justify center] -row 0 -column 0 -sticky ns -padx 2 -pady 2

   grid [ttk::frame $startGui.mainframe.entry] -row 1 -column 0 -sticky nwes
   grid columnconfigure $startGui.mainframe.entry 0 -weight 0
   grid columnconfigure $startGui.mainframe.entry 1 -weight 2
   grid rowconfigure $startGui.mainframe.entry 0 -weight 1

   grid [ttk::label $startGui.mainframe.entry.sel -text "Selection: "] -row 0 -column 0 -sticky nwes -pady 2 -padx 2
   grid [ttk::entry $startGui.mainframe.entry.entry -textvar [namespace current]::atomsel] -row 0 -column 1 -sticky nwes -pady 2 -padx 2
   grid [ttk::button $startGui.mainframe.entry.go -text "Start Molefacture" -command "[namespace current]::molefacture_gui_aux $[namespace current]::atomsel"] -row 0 -column 2 -sticky nwes -pady 2 -padx 2
}

proc ::Molefacture::molefacture_gui_aux { {seltext ""}  {reset 1}} {
  variable origseltext
  variable origmolid
  variable tmpmolid
  variable startGui
  #### Reset Molefacture variables
  if {$reset == 1} {
    set mouse $Molefacture::mousemode
    Molefacture::initialize
    Molefacture::updateTempSel
    set Molefacture::mousemode $mouse
  }

  set origmolid [molinfo top]

  set mysel ""
  if {[string trim $seltext " "] == ""} {
    set mysel ""
  } elseif {[catch {atomselect $origmolid $seltext} mysel] != 0} {

      tk_messageBox -icon error -type ok -title Error \
      -message "You entered an invalid atom selection." -parent $Molefacture::startGui
      return
  }

  if {$seltext != ""} {
    # if {[$mysel num] > 200} {
    #   tk_messageBox -icon error -type ok -title Warning \
    #      -message "The current version of molefacture is best used on structures of 200 atoms or smaller.\
    #       Future versions will be able to handle larger structures. You may continue, but some features may work slowly.\
    #        See the molefacture documentation for more details."
    # } else

    if {[$mysel num] == 0} {

      tk_messageBox -icon error -type ok -title Warning \
         -message "The atom selection returned 0 atoms. Please select at lease one atom or proceed with an empty molecule." -parent $Molefacture::startGui
           return
    }

    set origseltext $seltext
  }
  ::Molefacture::new_mol

  if {$seltext != ""} {
    ### the reset works to also as a flag to use autopsf in the reloading of the selection
    ::Molefacture::reload_selection $reset
    if {[llength $Molefacture::unknownRes] == 0} {
      set Molefacture::topocgenff 1
    }
  }
  update_openvalence
  ::Molefacture::updateAtomTable
  ::Molefacture::clearTableSelection
  if {$mysel != ""} {$mysel delete}
  if {[winfo exists $startGui] == 1} {
    set ::Molefacture::atomsel ""
    wm withdraw $startGui
  }


}

###################################################
# Update Atom Table                               #
###################################################

proc ::Molefacture::updateAtomTable {} {
  variable tmpmolid
  variable oxidation
  variable valence
  $Molefacture::atomTable delete 0 end

  set sel [atomselect $tmpmolid "$Molefacture::notDumSel"]

  set resnamelist [$sel get resname]
  set residlist [$sel get resid]
  set chainlist [$sel get chain]
  set chargelist [$sel get charge]
  set segnamelist [$sel get segname]
  set chargepenalty [$sel get user]

  ### Add the column for the Charge Penalty values if they exist
  set chrgpenal 0

  if {[lsort -unique $chargepenalty] != "-99" && [$sel num] > 0} {
    if {[$Molefacture::atomTable findcolumnname "CrhgPenal"] == -1} {
      $Molefacture::atomTable insertcolumns 6 0 "Chrg. Penalty" center
      $Molefacture::atomTable columnconfigure 6 -name CrhgPenal -sortmode real -editable false -selectbackground cyan -width 0 -maxwidth 0
    }
    set chrgpenal 1
  } else {
    if {[$Molefacture::atomTable findcolumnname "CrhgPenal"] != -1} {
      $Molefacture::atomTable deletecolumns CrhgPenal
    }

  }

  foreach index [$sel get index] resname $resnamelist resid $residlist chain $chainlist\
  name [$sel get name] type [$sel get type] element [$sel get element] charge $chargelist segname $segnamelist penalty $chargepenalty  {

    set oxlist [lindex [array get valence $element] 1]

    set oxi [lindex $oxlist [lindex $oxidation $index]]

    if {$oxi == ""} {
      set oxi "0"
    }

    if {$chrgpenal == 0} {
      set str "$index $name $type $element $oxi [Molefacture::format3Dec $charge] $resname $resid $chain $segname"
    } else {
      set str "$index $name $type $element $oxi [Molefacture::format3Dec $charge] [Molefacture::format3Dec $penalty] $resname $resid $chain $segname"
    }

    $Molefacture::atomTable insert end $str
  }

  set Molefacture::numAtoms [llength $resnamelist]

  $sel delete
  ::Molefacture::updateEditMolecule
}

###################################################
# Update Edit Molecule Gui values                 #
###################################################

proc ::Molefacture::updateEditMolecule {} {
  variable picklist
  variable tmpmolid
  if {[$Molefacture::atomTable size] > 0} {
    $Molefacture::rnamewdgt configure -state normal
    $Molefacture::chnwdgt configure -state normal
    $Molefacture::rsidwdgt configure -state normal
    $Molefacture::sgnmwdgt configure -state normal
  } else {
    return
  }

  set seltext "$Molefacture::notDumSel"
  if {[llength $picklist] != 0} {
    set seltext "index [join $picklist]"
  }

  set sel [atomselect $tmpmolid "$seltext"]
    set resname [lsort -unique -dictionary [$sel get resname]]
    if {[llength $resname] > 1} {
      set Molefacture::resname "---"
    } else {
      set Molefacture::resname $resname
    }

    set resid [lsort -unique -integer [$sel get resid]]
    if {[llength $resid] > 1} {
      set Molefacture::resid "---"
    } else {
      set Molefacture::resid $resid
    }

    set chain [lsort -unique -dictionary [$sel get chain]]
    if {[llength $chain] > 1} {
      set Molefacture::chain "---"
    } else {
      set Molefacture::chain $chain
    }
    set chargelist [$sel get charge]
    if {[llength $chargelist] > 1} {
      set charge [eval "vecadd $chargelist"]

        set Molefacture::charge [Molefacture::format3Dec $charge]
      
    } else {
      set Molefacture::charge [Molefacture::format3Dec $chargelist]
    }


    set segname [lsort -unique -dictionary [$sel get segname]]
    if {[llength $segname] > 1} {
      set Molefacture::segname "---"
    } else {
      set Molefacture::segname $segname
    }

  $sel delete
}

###################################################
# Update Oxidation State Combobox Values          #
###################################################
proc ::Molefacture::oxiStateCombo {} {
  variable picklist
  variable tmpmolid
  variable valence
  variable oxidation
  if {[$Molefacture::atomTable size] > 0 && [llength $picklist] > 0} {

    set sel [atomselect $tmpmolid "index [join $picklist]"]
    set element [lsort -unique -dictionary [$sel get element]]
    if {[llength $element] > 1} {
      set Molefacture::cmboxstatvalues [list]
      set Molefacture::cmboxstat ""
      $Molefacture::cmboxstatwidget configure -state disabled
      $sel delete
    } else {
      ###
      ### Check if all atoms of the same element
      ### have the same oxidation state. If
      ### not, set cmboxstat oxidation state as ---
      ### else, set cmboxstat the oxidation state.
      ###
      $Molefacture::cmboxstatwidget configure -state normal

      set oxlist [lindex [array get valence $element] 1]

      set curox [lindex $oxidation [lindex $picklist 0] ]
      set newcurox $curox
      if {[llength $picklist] > 1} {
        set pickaux [lrange $picklist 1 end]

        foreach oxiaux $pickaux {
          if {[lindex $oxidation $oxiaux] != $curox} {
            set newcurox -1
            break
          }
        }
      }

      if {$newcurox == -1} {
        set Molefacture::cmboxstat "---"
      } else {
        set Molefacture::cmboxstat [lindex $oxlist $newcurox]
      }
      $sel delete
      $Molefacture::cmboxstatwidget configure -values $oxlist
      set Molefacture::cmboxstatvalues $oxlist
    }
    update

  } else {
    set Molefacture::cmboxstatvalues [list]
    set Molefacture::cmboxstat ""
    $Molefacture::cmboxstatwidget configure -state disabled -values [list]
  }
  return 1
}


###################################################
# change Oxidation State of the selected atoms    #
###################################################
proc ::Molefacture::updateOxiState {} {
  variable picklist
  variable tmpmolid
  variable valence
  variable oxidation

  set tblindex [$Molefacture::atomTable curselection]
  set element [$Molefacture::atomTable cellcget [lindex $tblindex 0],Element -text]
  set oxlist [lindex [array get valence $element] 1]

  Molefacture::saveState

  set oxstat [lsearch $oxlist $Molefacture::cmboxstat]

  if {$oxstat == -1} {
    set list [lappend valence($element) $Molefacture::cmboxstat]
    
    set valence($element) $list

    set oxstat [lsearch $list $Molefacture::cmboxstat]
  }
  if {[$Molefacture::atomTable size] > 0 && [llength $picklist] > 0} {
    set delatom [list]

    foreach ind $picklist tbind $tblindex {

      ### set new oxidation state
      lset oxidation $ind $oxstat
      $Molefacture::atomTable cellconfigure $tbind,OxState -text $Molefacture::cmboxstat
      set selatom [atomselect $tmpmolid "index $ind"]

      set atombonds [join [$selatom getbonds]]
      set numbonds [llength $atombonds]
      if {$numbonds == 0} {
        $selatom delete
        continue
      }

      #### account with bond orders
      set bondorder 1
      if {$numbonds > 1} {
         set bondorder [expr int([vecsum [join [$selatom getbondorders]] ])]
      }

      if {$Molefacture::cmboxstat < $bondorder && [llength $atombonds] > 1 } {

        set bondsel [atomselect $tmpmolid "index $atombonds and occupancy 0.5"]

        set bondind [$bondsel get index]
        set bdninc 0
        set numbonds [expr $Molefacture::cmboxstat - $bondorder]
        while {$numbonds  < 0 &&  $bdninc < [llength $bondind]} {
          lappend delatom [lindex $bondind $bdninc]
          incr bdninc
          incr numbonds
        }
        $bondsel delete

      }
      $selatom delete
    }

    if {[llength $delatom] > 0} {
      set oldpicklist $picklist
      set picklist $delatom
      ::Molefacture::del_atom_gui
      set picklist $oldpicklist
      set delatom [list]
    }
    if {$Molefacture::topocgenff == 1} {
      set Molefacture::topocgenff 0
    }
  }
  return
}

###################################################
# Minimize structure using Open Babel for UFF     #
# or NAMD for CHARMM                              #
###################################################

proc Molefacture::minimize {} {
  global env
  variable tmpmolid
  if {$tmpmolid == -1 || [lsearch [molinfo list] $tmpmolid] == -1} {
    return
  }
  set pwd [pwd]
  cd $env(TMPDIR)
  if {[file exist molefacture_tmp] != 1} {
    file mkdir molefacture_tmp
  }

  Molefacture::saveState

  ### Save the view point before the minimization
  set viewPointsList [list]

  set viewPointsList [list [molinfo $tmpmolid get rotate_matrix]\
      [molinfo $tmpmolid get center_matrix]\
      [molinfo $tmpmolid get scale_matrix]\
      [molinfo $tmpmolid get global_matrix]]


  cd molefacture_tmp

  if {$Molefacture::minff == "UFF"} {

    set openbabel [::ExecTool::find "obminimize"]
    if {[string length $openbabel] == ""} {
      tk_messageBox -message "To minimize the structure using UFF, please install Open Babel software,\
       and include its path in the global Path variable." -icon info -type ok -parent $Molefacture::topGui
       cd $pwd
       return
    }


    set sel [atomselect $tmpmolid "$Molefacture::notDumSel"]
    if {[Molefacture::exportFile $sel minimize_vmd.mol2 mol2] == -1} {
      cd $pwd
      return
    }
    set done 0
    file delete -force minimize_obabel.xyz

    catch {eval ::ExecTool::exec "obminimize -o xyz -n $Molefacture::minsteps \
    -ff UFF minimize_vmd.mol2 > minimize_obabel.xyz"} done
    if {[file exists minimize_obabel.xyz] == 1} {
      set file [open minimize_obabel.xyz r]
      set lines [read $file]
      set lines [split $lines "\n"]
      close $file

      set file [open minimize_obabel.xyz w+]

      foreach line $lines {
        if {[lindex $line 0 0] != "WARNING:"} {
           puts $file $line
        }
      }
      close $file
      set obabelmol [mol new minimize_obabel.xyz]

      set newcoorsel [atomselect $obabelmol "all"]
      $sel set {x y z} [$newcoorsel get {x y z}]

      $newcoorsel delete
      mol delete $obabelmol
    } else {
      puts "Error: $done"
    }
    
    $sel delete

    cd $pwd


  } elseif {$Molefacture::minff == "CHARMM36"} {

    set namd [::ExecTool::find "namd2"]
    if {$namd == ""} {
      tk_messageBox -message "To minimize the structure using CHARMM36, please install NAMD software,\
       and include its path in the global Path variable." -icon info -type ok  -parent $Molefacture::topGui
       cd $pwd
       return
    }

    set sel [atomselect $tmpmolid "$Molefacture::notDumSel or segname CAP"]
    
    # if {[Molefacture::findHalogens $sel] == 1} {
    #   tk_messageBox -message "Halogen elements are currently not supported in the automatic \
    #   generation of topologies. Please use \
    #   UFF to minimize the structure." -title "Minimization failure" -icon error -type ok \
    #   -parent $Molefacture::topGui
    #   $sel delete
    #   return
    # }

    set capSel [atomselect $tmpmolid "segname CAP"]

    $sel set beta 0
    set haveCap 0
    if {[$capSel num] > 0} {
      set haveCap 1
    }
    $sel delete
    set fixAtomSel ""
    foreach bond [$capSel getbonds] index [$capSel get index]  {
        set selaux [atomselect $tmpmolid "index [lindex $bond 0]"]
        lappend fixAtomSel "[$selaux get resid] [$selaux get name]"
        $selaux delete
        set haveCap 1
    }

    set sel [atomselect $tmpmolid "$Molefacture::notDumSel"]
    set resname [lindex [$sel get resname] 0]
    set resultRes [Molefacture::loadDefaulTopo [lsort -unique -dictionary [$sel get resname]]]
    update idletasks
    
  
    if {[llength [lindex $resultRes 0] ] > 0} {
      ### Residue is not defined in the default topologies.
      ### Is it defined in the new defined topologies?
      set Molefacture::unknownRes [lindex $resultRes 0] 
      set found 0 
      foreach top $Molefacture::topoinfo {
        foreach res $Molefacture::unknownRes {
          set index [lsearch [lindex $top 6] $res]
          if {$index > -1} {
              set Molefacture::unknownRes [lreplace $Molefacture::unknownRes $index $index]
          }
        }
      }
      

      if {[llength $Molefacture::unknownRes ] > 0} {
        set Molefacture::topocgenff 0  
        
        ### because of the trace in the varibale Molefacture::topocgenff
        ### the molefacture_tmp fodler is deleted. Check if that is the case
        if {[file exist molefacture_tmp] != 1} {
          file mkdir molefacture_tmp
          cd molefacture_tmp
        }
      }
    }
    
    set continue [Molefacture::write_psfpdbfiles $resname]
    $sel delete

    if {[file exists $resname.psf] != 1 || $continue == -1} {
      cd $pwd
      return
    }

    if {$haveCap == 1} {
      set newmol [mol new $resname.psf]
      mol addfile $resname.pdb waitfor all $newmol
      set selAll [atomselect $newmol "all"]
      $selAll set beta 0

      foreach fix $fixAtomSel  {
        
        set selaux [atomselect $newmol "resid [lindex $fix 0] and backbone"]
        $selaux set beta 1
        $selaux delete

      }

      $selAll writepdb fixatoms.pdb
      $selAll delete
      mol delete $newmol
    }

    set namdfile [open "minimize.conf" w+]

    puts $namdfile "### NAMD configuration file to minimize structures modeled in Molefacture\n\n"

    puts $namdfile "structure $resname.psf"
    puts $namdfile "coordinates $resname.pdb"
    puts $namdfile "temperature 0"
    puts $namdfile "paratypecharmm on"
    set line ""
    if {[info exists QWIKMD::ParameterList]} {
      foreach param $QWIKMD::ParameterList {
        append line "parameters $param\n"
      }
    }
    if {$Molefacture::topocgenff == 1} {
      if {[lsearch -index 0 $QWIKMD::ParameterList "*par_all36_cgenff.prm"] == -1} {
        append line "parameters ${env(CHARMMPARDIR)}/par_all36_cgenff.prm\n"
      }
      foreach top $Molefacture::topoinfo {
        ### Check if the user generated the top using the info in VMD
        ### because he/she loaded the psf and pdb, making $Molefacture::topocgenff == 1 
        ### but if the the residue is not defined in any topology file,
        ### [llength $Molefacture::unknownRes] > 0.
        if {[llength $Molefacture::unknownRes] == 0} {
          append line "parameters [lindex $top 0]\n"
        }
      }

    }
    if {$haveCap == 1} {
      append line "fixedAtoms on\nfixedAtomsFile fixatoms.pdb\nfixedAtomsCol B"
    }

    puts $namdfile $line
    
    puts $namdfile "outputname tmp_namd"
    puts $namdfile "dcdfreq 0"
    puts $namdfile "binaryoutput no"
    puts $namdfile "cutoff 12"
    puts $namdfile "pairlistdist 14.0"
    puts $namdfile "switching on"
    puts $namdfile "switchdist 10"
    puts $namdfile "exclude scaled1-4"
    puts $namdfile "1-4scaling 1.0"
    puts $namdfile "rigidbonds all"

    puts $namdfile "minimize $Molefacture::minsteps"
    close $namdfile

    set done 0
    catch {eval ::ExecTool::exec "namd2 +idlepoll +setcpuaffinity minimize.conf >> minimize.log"} done
    if {[file exists tmp_namd.coor] == 1} {
      set namdmol [mol new $resname.psf]
      mol addfile tmp_namd.coor type pdb
      set sel [atomselect $tmpmolid "$Molefacture::notDumSel"]
      foreach resid [$sel get resid] {
        set newcoorsel [atomselect $namdmol "resid $resid"]
        set names [$newcoorsel get name]
        $newcoorsel delete

        set tmpmSel [atomselect $tmpmolid "(resid $resid and ($Molefacture::notDumSel)) and name [join $names]"]
        set minSel [atomselect $namdmol "(resid $resid and ($Molefacture::notDumSel)) and name [join $names]"]
        set indecesList [list]
        set topSelNames [$tmpmSel get name]
        set minSelNames [$minSel get name]
        foreach name $topSelNames {
          lappend indecesList [lsearch $minSelNames $name]
        }

        ### sort the atoms based on atom's name to ensure that the atoms from both selection match
        ### to avoid possible reordering of atoms from psfgen or addition of hydrogen atoms on
        ### open terminal during protein building operations
        proc sortlist {list indeces} {
          set outlist [list]
          foreach ind $indeces {
            lappend outlist [lindex $list $ind]
          }
          return $outlist
        }

        $tmpmSel set {x y z} [Molefacture::sortlist [$minSel get {x y z}] $indecesList ]
        $tmpmSel set charge [Molefacture::sortlist [$minSel get charge] $indecesList]
        $tmpmSel set type [Molefacture::sortlist  [$minSel get type] $indecesList]
        $tmpmSel delete
        $minSel delete

      }

      mol delete $namdmol
      $sel delete

      Molefacture::updateprotlasth
      Molefacture::updateAtomTable

    } else {
      puts "Error: $done"
    }
  }
  cd $pwd

  ### Restore the view point before the minimization

  molinfo $tmpmolid set rotate_matrix [lindex $viewPointsList 0]
  molinfo $tmpmolid set center_matrix [lindex $viewPointsList 1]
  molinfo $tmpmolid set scale_matrix [lindex $viewPointsList 2]
  molinfo $tmpmolid set global_matrix [lindex $viewPointsList 3]

}

##########################################
# Check if ther is any halogen atom in   #
# the selection                          #
##########################################
proc Molefacture::findHalogens { sel } {
  ### Check if the selection has halogens in it. So far, psfgen does not know 
  ### what to do with them, so bail out. Tell user to export mol2 file and 
  ### submit by hand
  set elements [lsort -unique -dictionary [$sel get element]]
  set hal {"F" "Cl" "Br" "I" "At"}
  set found 0
  foreach ele $hal {
    if {[lsearch $elements $ele] > -1} {
      set found 1
      break
    }
  }
  
  return $found
}

###################################################
# Clean up and quit the program.                  #
###################################################

proc ::Molefacture::done { {force 0} } {
   fix_changes
   variable projectsaved
   variable slavemode

   if {!$projectsaved && !$force && !$slavemode} {
      set reply [tk_dialog .foo "Quit - save file?" "Quit Molefacture - Save molecule?" \
        questhead  0 "Save" "Don't save" "Cancel"]
      switch $reply {
	 0 { fix_changes; export_molecule_gui}
	 1 { }
	 2 { return 0 }
      }
   }

   # Make the master molecule visible again
   variable molidorig
   if {$molidorig != -1} {molinfo $molidorig set drawn 1}

   # Remove_traces
   foreach t [trace info variable ::vmd_pick_event] {
      trace remove variable ::vmd_pick_event write ::Molefacture::atom_picked_fctn
   }

   # Set mouse to rotation mode
   mouse mode 0
   mouse callback off;

   if { [winfo exists .molefac] }    { wm withdraw .molefac }

   if {$slavemode} {
      variable exportfilename
      variable mastercallback
      variable atomlist
      variable bondlist
      fix_changes
      export_molecule $exportfilename
      $mastercallback $exportfilename
   }

   variable tmpmolid
   mol delete $tmpmolid

   # Cleanup
   if {[file exists Molefacture_tmpmol.xbgf]} { file delete Molefacture_tmpmol.xbgf }

   # Close molefacture
}

#######################################################
# Proc to be called from external plugins             #
# This will open Molefacture windows                  #
#######################################################

proc ::Molefacture::extCall { args } {
  set numargs [llength $args]
  if {$numargs == 0} {extCall_usage}

  set argnum 0
  set arglist $args
  set qwikmd 0

  molefacture_tk

  foreach i $args {
    if {$i == "-qwikmd"} {
      set qwikmd 1

      set arglist [lreplace $arglist $argnum $argnum]
      continue
    }
    incr argnum
  }


  set topo ""
  set otherarglist {}
  foreach {i j} $arglist {
    if {$j == ""} {
      continue
    }
    if {$i=="-mol"} {
      set ::Molefacture::molidorig $j
      continue
    }
    if {$i=="-sel"} {
      set ::Molefacture::atomsel $j
      continue
    }
    if {$i=="-topo"} {
      set topo $j
      continue
    }
  }


  menu molefacture on

  if {$::Molefacture::atomsel != ""} {
    ### During molefacture_gui_aux proc, all variables are re-initialized
    ::Molefacture::molefacture_gui_aux $::Molefacture::atomsel

  }

  set Molefacture::qwikmd $qwikmd
  if {[file exists $topo] == 1} {
    set handler [::Toporead::read_charmm_topology $j 1]
    set Molefacture::topoinfo [::Toporead::topology_from_handler $handler]
    set Molefacture::topocgenff 1
  }
}

proc ::Molefacture::extCall_usage { } {
  puts "Usage: ::Molefacture::extCall <option1> <option2>...\n"
  puts "Options:"
  puts "  Switches:"
  puts "    -qwikmd called from QwikMD"
  puts "  Double option arguments:"
  puts "    -mol <molecule> (REQUIRED: Run on molid <molecule>)"
  puts "    -sel <text> (Atom selection text containing the molecule(s) to be edited)"
  puts "    -topo <topology file> (Initial Topology file)"
  error ""
}
proc molefacture_tk {} {

   #### Initialize variables
    Molefacture::initialize
    Molefacture::molefacture_gui
    Molefacture::moleculeChanged
    Molefacture::loadDefaulTopo
    Molefacture::setHotKeys
    return $Molefacture::topGui
}

#Load all procs from other molefacture files
foreach lib { molefacture_builder.tcl \
              molefacture_state.tcl \
              molefacture_geometry.tcl \
              molefacture_gui.tcl \
              molefacture_edit.tcl \
              molefacture_internals.tcl\
              molefacture_balloontext.tcl } {
  if { [catch { source [file join $env(MOLEFACTUREDIR) $lib] } errmsg ] } {
    error "ERROR: Could not load molefacture library file $lib: $errmsg"
  }
}
