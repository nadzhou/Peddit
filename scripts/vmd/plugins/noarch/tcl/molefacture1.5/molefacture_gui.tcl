##
## $Id: molefacture_gui.tcl,v 1.71 2019/07/02 17:52:12 jribeiro Exp $
##

###############################################################################
# NEW Molefacture Version = 1.5 (major gui overhaul)
# 
# The new interface is based on two main areas/panels:
#   The left area where different modes can be implemented, by changing this 
#   panel to a notebook and adding tabs to it. As for the first 
#   implementation for the first version, only a table (and no notebook)
#   with the atoms fields was implemented. One can imagine using this left panel 
#   to extend the functionality of molefacture like add QM input file generator.
#
#   the right panel contains all the buttons and action related widgets, like
#   slide bars and modifiers that affect the molecule.
#
#
#   Many of the functionality previously available in the Molefacture were 
#   disabled to make sure that the modeling part of the tool worked properly
#   slowly these functionalities should be re-implemented.
#
###############################################################################
proc ::Molefacture::molefacture_gui { {sel ""} } {

   variable topGui

   ### define the style of combobox when in readonly state
   ttk::style map TCombobox -fieldbackground [list readonly #ffffff]
   ttk::style map TEntry -fieldbackground [list readonly #ffffff]

   # If already initialized, just turn on
   if {[winfo exists $Molefacture::topGui] == 1} {
      wm deiconify $Molefacture::topGui
      raise $Molefacture::topGui

      ### Re-bind the Atom's table selection lost event
      bind $Molefacture::atomTable <<TablelistSelectionLost>> Molefacture::bindAtomTable

      return $Molefacture::topGui
   }


   toplevel $topGui

   ## Title of the windows
   wm title $topGui "Molefacture - Molecule Builder"
   ###proc triggered when Molefacture is closed
   bind $Molefacture::topGui <Button-1> {
     if {$Molefacture::bindTop == 0} {
         wm protocol $Molefacture::topGui WM_DELETE_WINDOW Molefacture::closeMainWindow
         set Molefacture::bindTop 1
      }
   }

   ### create the images for the arrow right and down to
   ### replace the unicode breaking in DCV

   image create photo arrowRight -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons arrow_right.gif]

   image create photo arrowDown -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons arrow_down.gif]

   set Molefacture::arrowRight "arrowRight"
   set Molefacture::arrowDown "arrowDown"

   grid columnconfigure $topGui 0 -weight 1
   grid columnconfigure $topGui 1 -weight 0
   grid rowconfigure $topGui 0 -weight 0
   grid rowconfigure $topGui 1 -weight 1

   ## main frame
   grid [ttk::frame $topGui.topframe -relief groove] -row 0 -column 0 -sticky nwes -padx 2 -pady 2
   grid columnconfigure $topGui.topframe 0 -weight 1
   grid rowconfigure $topGui.topframe 0 -weight 1


   ## menubar frame
   grid [ttk::frame $topGui.topframe.menubar] -row 0 -column 0 -sticky nswe -pady 2 -padx 2
   grid columnconfigure $topGui.topframe.menubar 3 -weight 1
   grid rowconfigure $topGui.topframe.menubar 0 -weight 1


   grid [menubutton $topGui.topframe.menubar.file -text "File" -menu $topGui.topframe.menubar.file.menu -width 8] -row 0 -column 0 -sticky ew
   grid [menubutton $topGui.topframe.menubar.build -text "Build" -menu $topGui.topframe.menubar.build.menu -width 9] -row 0 -column 1 -sticky ew

   grid [menubutton $topGui.topframe.menubar.settings -text "Settings" -menu $topGui.topframe.menubar.settings.menu -width 9] -row 0 -column 2 -sticky ew
   grid [ttk::frame $topGui.topframe.menubar.strchfrm] -row 0 -column 3 -sticky ew -pady 2
   grid [menubutton $topGui.topframe.menubar.help -text "Help" -menu $topGui.topframe.menubar.help.menu -width 9] -row 0 -column 4 -sticky ew

   menu $topGui.topframe.menubar.file.menu -tearoff no

   menu $topGui.topframe.menubar.file.menu.newmenu -tearoff no
   $topGui.topframe.menubar.file.menu.newmenu add command -label "From selection" -command ::Molefacture::molefacture_start
   $topGui.topframe.menubar.file.menu.newmenu add command -label "From scratch" -command {
      Molefacture::molefacture_gui_aux
      set Molefacture::origseltext ""
      set Molefacture::origmolid -1
   }

   ### New molecule
   $topGui.topframe.menubar.file.menu add cascade -label "New molecule" -menu $topGui.topframe.menubar.file.menu.newmenu

   ### Apply changes
   $topGui.topframe.menubar.file.menu add command -label "Apply changes to parent" -command ::Molefacture::apply_changes_to_parent_mol
   
   ### Export molecule as
   $topGui.topframe.menubar.file.menu add command -label "Write MOL2 file" -command {Molefacture::run_cgenff export ""}
   $topGui.topframe.menubar.file.menu add command -label "Write XYZ file" -command {Molefacture::run_cgenff export "" xyz}
   $topGui.topframe.menubar.file.menu add command -label "Write top/str file" -command ::Molefacture::write_topology_gui
   $topGui.topframe.menubar.file.menu add command -label "Write psf and pdb files" -command ::Molefacture::write_psfpdbfiles

   ## Build menu
   menu $topGui.topframe.menubar.build.menu -tearoff no
   menu $topGui.topframe.menubar.build.menu.abtype -title "Topology & Parameters" -tearoff no


   if {$Molefacture::cgenffcaller == 0} {
      $topGui.topframe.menubar.build.menu.abtype add command -label "Export mol2 to CGenFF" \
      -command {::Molefacture::run_cgenff export}
   }

   $topGui.topframe.menubar.build.menu.abtype add command -label "Import CGenFF Topology" \
    -command {::Molefacture::run_cgenff import}

   $topGui.topframe.menubar.build.menu.abtype add command -label "Submit to CGenFF Server" \
    -command Molefacture::CgenffServerCall

   $topGui.topframe.menubar.build.menu.abtype add command -label "Export to ffTK" \
    -command ::Molefacture::export_ffTK
   $topGui.topframe.menubar.build.menu.abtype add command -label "Guess Bond Order (IDATM)" \
    -command ::Molefacture::run_idatm

   $topGui.topframe.menubar.build.menu add cascade -label "Topology & Parameters" \
    -menu $topGui.topframe.menubar.build.menu.abtype

   menu $topGui.topframe.menubar.build.menu.addfrag -title "Replace atom with fragment" -tearoff no
   
   ### Read the frag.mdb file and add the the global array variable addfrags
   read_fragment_file [file join $::env(MOLEFACTUREDIR) lib fragments frag.mdb]
   ### Populate the fragments menu with the fragments in the lib/fragments folder
   fill_fragment_menu


   menu $topGui.topframe.menubar.build.menu.basefrag -tearoff no
   
   ### Same as the fragments menu
   read_basefrag_file [file join $::env(MOLEFACTUREDIR) lib basemol basefrag.mdb]
   fill_basefrag_menu

   $topGui.topframe.menubar.build.menu add cascade -label "Replace atom with fragment" \
    -menu $topGui.topframe.menubar.build.menu.addfrag
   $topGui.topframe.menubar.build.menu add cascade -label "Add independent fragment" \
    -menu $topGui.topframe.menubar.build.menu.basefrag

   $topGui.topframe.menubar.build.menu add command -label "Duplicate Selected Atoms" -command {
      Molefacture::duplicateSel
   }
   $topGui.topframe.menubar.build.menu add command -label "Protein Builder" -command ::Molefacture::prot_builder_gu

   ## Settings menu

   menu $topGui.topframe.menubar.settings.menu -tearoff no
   $topGui.topframe.menubar.settings.menu add command -label "CGenFF Settings" -command ::Molefacture::cgenffUserWindows

   $topGui.topframe.menubar.settings.menu add cascade -label "Default Topologies" -menu $topGui.topframe.menubar.settings.menu.defTop

   menu $topGui.topframe.menubar.settings.menu.defTop -tearoff no

   set Molefacture::defTopMenu $topGui.topframe.menubar.settings.menu.defTop
   ## Help menu
   menu $topGui.topframe.menubar.help.menu -tearoff no
   $topGui.topframe.menubar.help.menu add command -label "About" \
    -command {tk_messageBox -type ok -title "About Molefacture" \
            -message "Molecule editing tool"}
  $topGui.topframe.menubar.help.menu add command -label "Help..." \
    -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/molefacture"


   ############## Frame for Atom Table List and Editing #################

   set atomframe "$topGui.atomframe"
   grid [ttk::frame $atomframe] -row 1 -column 0 -sticky nswe -padx 2 -pady 2
   grid columnconfigure $atomframe 0 -weight 1
   grid rowconfigure $atomframe 0 -weight 1

   set tbframe "$topGui.atomframe.tbframe"
   grid [ttk::frame $tbframe] -row 0 -column 0 -sticky nswe
   grid columnconfigure $tbframe 0 -weight 1
   grid rowconfigure $tbframe 0 -weight 1

   option add *Tablelist.activeStyle       frame

   option add *Tablelist.movableColumns    no
   option add *Tablelist.labelCommand      tablelist::sortByColumn

   tablelist::tablelist $tbframe.tb \
        -columns { 0 "Index" center
                0 "Atom Name"    center
                0 "Atom Type" center
                0 "Element" center
                0 "OxState" center
                0 "Charge" center
                0 "Resname"  center
                0 "Res ID"   center
                0 "Chain ID"     center
                0 "Segname" center
                } \
                -yscrollcommand [list $tbframe.scr1 set] -xscrollcommand [list $tbframe.scr2 set] \
                -showseparators 0 -labelrelief groove -labelcommand tablelist::sortByColumn  -labelbd 1 -selectforeground black\
                -foreground black -background white -state normal -width 50 -stretch "all" -editselectedonly true -selectmode extended\
                -stripebackgroun white -exportselection true\
                -editstartcommand Molefacture::startEditTable -editendcommand Molefacture::validateEditTable -forceeditendcommand true

   $tbframe.tb columnconfigure 0 -selectbackground cyan -width 0 -maxwidth 0 -name Index -sortmode integer
   $tbframe.tb columnconfigure 1 -selectbackground cyan -width 0 -maxwidth 0 -editable true -editwindow ttk::entry -name AtmNAME -sortmode dictionary
   $tbframe.tb columnconfigure 2 -selectbackground cyan -width 0 -maxwidth 0 -editable true -editwindow ttk::entry -name AtmType -sortmode dictionary
   $tbframe.tb columnconfigure 3 -selectbackground cyan -width 0 -maxwidth 0 -editable false -name Element -sortmode dictionary
   $tbframe.tb columnconfigure 4 -selectbackground cyan -width 0 -maxwidth 0 -editable true -editwindow ttk::combobox -name OxState -sortmode real
   $tbframe.tb columnconfigure 5 -selectbackground cyan -width 0 -maxwidth 0 -editable true -editwindow ttk::entry -name Charge -sortmode real
   $tbframe.tb columnconfigure 6 -selectbackground cyan -width 0 -maxwidth 0 -editable true -editwindow ttk::entry -name Resname -sortmode dictionary
   $tbframe.tb columnconfigure 7 -selectbackground cyan -width 0 -maxwidth 0 -editable true -editwindow ttk::entry -name ResID -sortmode integer
   $tbframe.tb columnconfigure 8 -selectbackground cyan -width 0 -maxwidth 0 -editable true -editwindow ttk::entry -name ChainID -sortmode dictionary
   $tbframe.tb columnconfigure 9 -selectbackground cyan -width 0 -maxwidth 0 -editable true -editwindow ttk::entry -name Segname -sortmode dictionary


   grid $tbframe.tb -row 0 -column 0 -sticky news
   grid columnconfigure $tbframe.tb 0 -weight 1; grid rowconfigure $tbframe.tb 0 -weight 1

   ## Scroll Bar Vertical
   scrollbar $tbframe.scr1 -orient vertical -command [list $tbframe.tb  yview]
   grid $tbframe.scr1 -row 0 -column 1  -sticky ens

   ## Scroll Bar Horizontal
   scrollbar $tbframe.scr2 -orient horizontal -command [list $tbframe.tb xview]
   grid $tbframe.scr2 -row 1 -column 0 -sticky swe

   set Molefacture::atomTable $tbframe.tb


   bind $tbframe.tb <<TablelistSelect>> {
      if {[$Molefacture::atomTable size] > 0} {
         ::Molefacture::selectAtomTable
      }

   }

   bind $tbframe.tb <<TablelistSelectionLost>> Molefacture::bindAtomTable


   ############## buttons frame #################


   set cmdframe "$topGui.atomframe.cmdframe"
   grid [ttk::frame $cmdframe] -row 0 -column 1 -sticky nwe -padx 2 -pady 2
   grid columnconfigure $cmdframe 0 -weight 1
   grid rowconfigure $cmdframe 0 -weight 1

   set mainrow 0

   ############## Mouse Mode Label Frame #################

   grid [ttk::labelframe $cmdframe.mouse -text "Actions"] -row $mainrow -column 0 -sticky ewns -padx 2 -pady 2
   grid columnconfigure $cmdframe.mouse 0 -weight 1
   grid rowconfigure $cmdframe.mouse 0 -weight 1
   grid rowconfigure $cmdframe.mouse 1 -weight 1

   set col 0
   set row 0
   set width 35
   set height 50

   grid [frame $cmdframe.mouse.iconframe -background white] -row 0 -column 0 -sticky ewns -padx 2 -pady 2

   ### add atom

   ### Create image to be applied to the radio/check/tk button to give the look 
   ### of image and the label beneath the image. The second image has the 
   ### background green to make it easier to identify the selected mode.
   ### All the images are stored in the lib/icons

   ### On the most recent MacOS, the radiobuttons are not showing the labels.
   ### Still need to check this.

   set file add_atom.gif
   set selfile add_atom_selected.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file add_atom_mac.gif
      set selfile add_atom_selected_mac.gif
   }

   image create photo imag1 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]
   image create photo imag1sel -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $selfile]

   ### The Molefacture::mousemode keeps track of the active mouse mode
   grid [radiobutton $cmdframe.mouse.iconframe.addatom -borderwidth 1 -anchor center \
    -image imag1 -selectimage imag1sel -compound top \
    -padx 0 -pady 0 \
    -text "Add" -background white -activebackground white \
    -command Molefacture::addAtomMain -indicatoron false\
    -variable Molefacture::mousemode -value addatom -width $width -height $height] \
    -row $row -column $col -padx 1 -sticky nsw
    
   incr col

   ### delete atom
   set file delete.gif
   set selfile delete_selected.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file delete_mac.gif
      set selfile delete_selected_mac.gif
   }

   image create photo img2 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]
   image create photo img2sel -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $selfile]

   grid [radiobutton $cmdframe.mouse.iconframe.delatom -borderwidth 1 -anchor center \
    -image img2 -selectimage img2sel -background white -activebackground white \
    -text "Del" -compound top -padx 0 -pady 0 -command Molefacture::delAtomMain \
    -variable Molefacture::mousemode -value delatom -indicatoron false -width $width \
    -height $height] -row $row -column $col -padx 1 -sticky nsw

   incr col


   ### select atom

   set file selectatom.gif
   set selfile selectatom_sel.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file selectatom_mac.gif
      set selfile selectatom_sel_mac.gif
   }

   image create photo img3 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]
   image create photo img3sel -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $selfile]

   ### Select Atom doesn't need a command as the mousemode variable is tracking
   ### for this variable. The other buttons need a specific command to trigger
   ### a specific command, when the button is pressed.
   grid [radiobutton $cmdframe.mouse.iconframe.selatom -borderwidth 1 -indicatoron false \
    -anchor center -background white -activebackground white -compound top -text "Select" \
    -image img3 -selectimage img3sel -width $width -height $height -padx 0 -pady 0 -variable Molefacture::mousemode \
    -value pick -padx 1] -row $row -column $col -padx 1 -sticky nsw

   incr col

   ### rotate scene

   set file rotate.gif
   set selfile rotate_selected.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file rotate_mac.gif
      set selfile rotate_selected_mac.gif
   }
   image create photo imgrotate -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]
   image create photo imgrotatesel -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $selfile]

   grid [radiobutton $cmdframe.mouse.iconframe.rotate -borderwidth 1 -indicatoron false \
    -anchor center -background white -activebackground white -compound top -text "Rotate" \
    -image imgrotate -selectimage imgrotatesel -variable Molefacture::mousemode -value rotate \
    -padx 1 -width $width -height $height] -row $row -column $col -padx 1 -sticky nsw

   incr col

   ### translate scene
   set file translate.gif
   set selfile translate_selected.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file translate_mac.gif
      set selfile translate_selected_mac.gif
   }

   image create photo imgtranslate -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]
   image create photo imgtranslatesel -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $selfile]

   grid [radiobutton $cmdframe.mouse.iconframe.translate -borderwidth 1 -indicatoron false -anchor center \
    -background white -activebackground white -compound top -text "Trans-\nlate" -image imgtranslate \
    -selectimage imgtranslatesel -variable Molefacture::mousemode -value translate -padx 1 \
    -width $width -height $height] -row $row -column $col -padx 1 -sticky nsw

   incr col



   ### Drop a row
   set row 1
   set col 0

   ### move fragment
   set file move_frag.gif
   set selfile move_frag_selected.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file move_frag_mac.gif
      set selfile move_frag_selected_mac.gif
   }
   image create photo imageaux2 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]
   image create photo imageaux2sel -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $selfile]

   grid [radiobutton $cmdframe.mouse.iconframe.movefrag -borderwidth 1 -indicatoron false -anchor center \
    -background white -activebackground white -compound top -text "Move\nFrag" -image imageaux2 \
    -selectimage imageaux2sel -width $width -height $height -padx 0 \
    -pady 0 -variable Molefacture::mousemode \
    -value movefrag] -row $row -column $col -padx 1 -sticky nsw

   incr col

   ### move atom
   set file move.gif
   set selfile move_selected.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file move_mac.gif
      set selfile move_selected_mac.gif
   }
   image create photo imagemoveatom -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]
   image create photo imagemoveatomsel -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $selfile]

   grid [radiobutton $cmdframe.mouse.iconframe.moveatom -borderwidth 1 -indicatoron false \
    -anchor center -background white -activebackground white -compound top -text "Move\nAtom" \
    -image imagemoveatom -selectimage imagemoveatomsel -width $width -height $height \
    -padx 0 -pady 0 -variable Molefacture::mousemode -value moveatom]\
     -row $row -column $col -padx 1 -sticky nsw

   incr col

   set file add_bond.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file add_bond_mac.gif
   }
   image create photo imageaux5 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]

   grid [button $cmdframe.mouse.iconframe.addbond -borderwidth 1 -anchor center -padx 0 -pady 0\
    -background white -activebackground white -text "Add\nBond"  -compound top -image imageaux5 \
    -width $width -height $height -padx 0 -pady 0 -command {
      set Molefacture::mousemode addbond
      Molefacture::addDelBond add
      set Molefacture::mousemode pick
   }] -row $row -column $col -padx 1 -sticky nsw

   ### Hide the text from the buttons and menubuttons, as these widgets are the only one showing
   ### the text and image on mac
   if {$::tcl_platform(os) == "Darwin"} {
      $cmdframe.mouse.iconframe.addbond configure -text "" 
   }
   

   incr col

   set file delete_bond.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file delete_bond_mac.gif
   }
   image create photo imageaux6 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]

   grid [button $cmdframe.mouse.iconframe.delbond -borderwidth 1 -width $width \
   -height $height -anchor center -background white -activebackground white \
   -text "Dele\nBond" -compound top -image imageaux6 -padx 0 -pady 0 -command {
    set Molefacture::mousemode delbond
    Molefacture::addDelBond del
    set Molefacture::mousemode pick
   } ] -row $row -column $col -padx 1 -sticky nsw

   if {$::tcl_platform(os) == "Darwin"} {
      $cmdframe.mouse.iconframe.delbond configure -text "" 
   }

   incr col

   ### increase or decrease bond order
   set file bond_order.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file bond_order_mac.gif
   }

   image create photo imageaux7 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]

   grid [menubutton $cmdframe.mouse.iconframe.incrorder -menu $cmdframe.mouse.iconframe.incrorder.menu \
    -borderwidth 1 -relief groove -padx 2 -pady 0 -indicatoron false -anchor center -background white \
    -activebackground white -image imageaux7 -compound top -text "Bond\nOrder" \
    -width $width -height $height -padx 0 -pady 0] -row $row -column $col -padx 1 -sticky nsw

   if {$::tcl_platform(os) == "Darwin"} {
      $cmdframe.mouse.iconframe.incrorder configure -text "" 
   }

   ### Add option to set the bond order
   menu $cmdframe.mouse.iconframe.incrorder.menu -tearoff no
   $cmdframe.mouse.iconframe.incrorder.menu add command -background white -label "Single" -command {::Molefacture::set_bondorder 1}
   $cmdframe.mouse.iconframe.incrorder.menu add command -background white -label "Double" -command {::Molefacture::set_bondorder 2}
   $cmdframe.mouse.iconframe.incrorder.menu add command -background white -label "Triple" -command {::Molefacture::set_bondorder 3}

   incr col

   trace remove variable ::Molefacture::mousemode write ::Molefacture::setMouseMode

   trace add variable ::Molefacture::mousemode write ::Molefacture::setMouseMode

   set row 2
   set col 0

   ### add hydrogens

   set file add_hydrogen.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file add_hydrogen_mac.gif
   }
   image create photo imageaux9 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]

   grid [button $cmdframe.mouse.iconframe.addH -borderwidth 1 -background white \
   -activebackground white -compound top -text "Add\nHydr" -width $width -height $height \
   -padx 0 -pady 0 -command Molefacture::addHMain -image imageaux9] \
   -row $row -column $col -padx 1 -sticky nsw

   if {$::tcl_platform(os) == "Darwin"} {
      $cmdframe.mouse.iconframe.addH configure -text "" 
   }

   incr col

   ### duplicate atoms

   set file duplicate.gif
   if {$::tcl_platform(os) == "Darwin"} {
      set file duplicate_mac.gif
   }

   image create photo imageaux10 -format gif -file [file join $::env(MOLEFACTUREDIR) lib icons $file]

   grid [button $cmdframe.mouse.iconframe.duplicate -borderwidth 1 -background white -activebackground white \
    -compound top -text "Dupl" -width $width -height $height\
    -padx 0 -pady 0 -command Molefacture::duplicateSel \
    -image imageaux10] -row $row -column $col -padx 1 -sticky nsw

   if {$::tcl_platform(os) == "Darwin"} {
      $cmdframe.mouse.iconframe.duplicate configure -text "" 
   }
   for {set i 0} {$i < $col} {incr i} {
      grid columnconfigure $cmdframe.mouse.iconframe $i -weight 0
   }


   
   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.selatom [Molefacture::getMouseBText select]
   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.movefrag [Molefacture::getMouseBText move]
   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.moveatom [Molefacture::getMouseBText moveatom]

   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.translate [Molefacture::getMouseBText translate]
   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.rotate [Molefacture::getMouseBText rotate]

   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.addatom [Molefacture::getMouseBText addatom]
   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.delatom [Molefacture::getMouseBText delatom]
   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.addbond [Molefacture::getMouseBText addbond]

   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.delbond [Molefacture::getMouseBText delbond]

   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.incrorder [Molefacture::getMouseBText incrorder]
   # TKTOOLTIP::balloon $cmdframe.mouse.edit.decrorder [Molefacture::getMouseBText decrorder]
   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.addH [Molefacture::getMouseBText addh]
   TKTOOLTIP::balloon $cmdframe.mouse.iconframe.duplicate [Molefacture::getMouseBText duplicate]

   ### Add templates
   grid [ttk::frame $cmdframe.mouse.temp] -row 2 -column 0 -sticky wns -padx 2 -pady 2
   grid columnconfigure $cmdframe.mouse.temp 1 -weight 1
   grid rowconfigure $cmdframe.mouse.temp 1 -weight 1

   grid [ttk::label $cmdframe.mouse.temp.lbl -text "Add\nSkel." -width 5] -row 0 -column 0 -sticky w -padx 2 -pady 2

   grid [ttk::frame $cmdframe.mouse.temp.selector] -row 0 -column 1 -sticky ewns -padx 2 -pady 2
   grid columnconfigure $cmdframe.mouse.temp.selector 1 -weight 1
   grid rowconfigure $cmdframe.mouse.temp.selector 0 -weight 1

   grid [ttk::radiobutton $cmdframe.mouse.temp.selector.selcycle -text "Cycle" -variable Molefacture::cyclSelector -value cycle -command Molefacture::updateTempSel] -row 0 -column 0 -padx 2 -sticky w

   grid [ttk::radiobutton $cmdframe.mouse.temp.selector.selchain -text "Aliph. Chain" -variable Molefacture::cyclSelector -value chain -command Molefacture::updateTempSel] -row 0 -column 1 -padx 2 -sticky w


   grid [ttk::frame $cmdframe.mouse.temp.action] -row 1 -column 1 -sticky ewns -padx 2 -pady 2
   grid columnconfigure $cmdframe.mouse.temp.action 1 -weight 1
   grid rowconfigure $cmdframe.mouse.temp.action 0 -weight 1

   grid [ttk::label $cmdframe.mouse.temp.action.nlbl -text "N. Atoms"] -row 0 -column 0 -padx 2 -sticky w

   grid [ttk::combobox $cmdframe.mouse.temp.action.cycle -state readonly -width 3 -values {3 4 5 6 7 8 9} -textvariable Molefacture::cyclNAtoms] -row 0 -column 1 -padx 2 -sticky w

   grid [spinbox $cmdframe.mouse.temp.action.chain -width 3 -increment 1 -from 1 -to 40 -textvariable Molefacture::alchainNAtoms -background white]  -row 0 -column 1 -padx 2 -sticky w

   grid [ttk::button $cmdframe.mouse.temp.action.add -command Molefacture::growTemp -text "Add" -padding "0 0 0 0"] -row 0 -column 2 -padx 2 -sticky we

   bind $cmdframe.mouse.temp.action.cycle <<ComboboxSelected>> {
      %W selection clear
   }

   lappend Molefacture::cyclSelecWidget $cmdframe.mouse.temp.action.cycle
   lappend Molefacture::cyclSelecWidget $cmdframe.mouse.temp.action.chain

   grid forget $cmdframe.mouse.temp.action.chain
  

   set Molefacture::mousemode "pick"
   #::Molefacture::setMouseMode
   incr mainrow


     ############## General Edits#################

   grid [ttk::frame $cmdframe.bigbuttons] -row $mainrow -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $cmdframe.bigbuttons 0 -weight 1
   grid rowconfigure $cmdframe.bigbuttons 0 -weight 1

   incr mainrow

   set bbrow 0

   grid [ttk::frame $cmdframe.bigbuttons.forcegeo] -row $bbrow -column 0 -pady 2 -sticky nwes
   grid columnconfigure $cmdframe.bigbuttons.forcegeo 0 -weight 1
   grid columnconfigure $cmdframe.bigbuttons.forcegeo 1 -weight 1
   grid rowconfigure $cmdframe.bigbuttons.forcegeo 0 -weight 1

   grid [ttk::button $cmdframe.bigbuttons.forcegeo.frcplanar -text "Force Planar" -padding "2 0 2 0" -command {
      Molefacture::saveState
     foreach mindex $Molefacture::picklist {
       Molefacture::set_planar $mindex
     }
   }] -row 0 -column 0 -sticky we
   # incr bbrow

   grid [ttk::button $cmdframe.bigbuttons.forcegeo.frctetra -text "Force Tetrahedral" -padding "2 0 2 0" -command {
      Molefacture::saveState
     foreach mindex $Molefacture::picklist {
       Molefacture::set_tetra $mindex
     }
   }] -row 0 -column 1 -sticky we -padx 2
   incr bbrow

   grid [ttk::button $cmdframe.bigbuttons.forcegeo.selectall -text "Select All" -padding "2 0 2 0" \
   -command ::Molefacture::selectAllTable] -row 1 -column 0 -sticky we -pady 2

   grid [ttk::button $cmdframe.bigbuttons.forcegeo.clearselection -text "Clear Selection" -padding "2 0 2 0" \
   -command ::Molefacture::clearTableSelection] -row 1 -column 1 -sticky we -pady 2

   TKTOOLTIP::balloon $cmdframe.bigbuttons.forcegeo.clearselection [Molefacture::getMouseBText clearsel]

   incr bbrow

   # grid [ttk::separator $cmdframe.bigbuttons.sep3 -orient horizontal] -row 5 -column 0 -sticky we -pady 4

   grid [ttk::button $cmdframe.bigbuttons.undo -text "Undo" -padding "2 0 2 0" -command Molefacture::undoProc] -row $bbrow -column 0 -sticky we -pady 2

   TKTOOLTIP::balloon $cmdframe.bigbuttons.undo [Molefacture::getMouseBText undo]


   incr bbrow

   grid [ttk::separator $cmdframe.bigbuttons.sep1 -orient horizontal] -row $bbrow -column 0 -sticky we -pady 2
   incr bbrow

   ############## Edit Atoms & Bonds #################


   grid [ttk::frame $cmdframe.edit] -row $mainrow -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $cmdframe.edit 0 -weight 1
   grid rowconfigure $cmdframe.edit 0 -weight 1

   incr mainrow
   grid [ttk::label $cmdframe.edit.prt -text "Atoms & Bonds" -image $Molefacture::arrowDown \
   -compound left] -row 0 -column 0 -sticky w -pady 1

   bind $cmdframe.edit.prt <Button-1> {
      ::Molefacture::hideFrame %W [lindex [grid info %W] 1] "Atoms & Bonds"
   }

   grid [ttk::frame $cmdframe.edit.fcolapse] -row 1 -column 0 -sticky nsew -pady 5
   grid columnconfigure $cmdframe.edit.fcolapse 0 -weight 1


   set framecolapse $cmdframe.edit.fcolapse

   grid [ttk::frame $framecolapse.element] -row 0 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.element 0 -weight 0
   grid columnconfigure $framecolapse.element 1 -weight 1
   grid rowconfigure $framecolapse.element 0 -weight 1

   grid [ttk::label $framecolapse.element.lbl -text "Element: "] -row 0 -column 0 -sticky w -pady 2
   set cmbelements {
      Hydrogen(H:1) Boron(B:5) Carbon(C:6) Nitrogen(N:7) Oxygen(O:8) Flourine(F:9) Phosphorus(P:15)\
      Sulfur(S:16) Clorine(Cl:17) Iron(Fe:26) Bromine(Br:35) Other...
   }
   grid [ttk::combobox $framecolapse.element.ele -values $cmbelements -justify left \
   -state readonly -textvariable Molefacture::elemcomb -style TblEdit.TCombobox -postcommand {
         $Molefacture::atomTable cancelediting
      }] -row 0 -column 1 -sticky w -pady 2

   set Molefacture::elemcomb "Carbon(C:6)"

   ### Bind the selection of the combobox with the changes in the table
   ### and the generation of the periodic table window in case of
   ### "Other ..." is selected
   bind $framecolapse.element.ele <<ComboboxSelected>> {
      if {$Molefacture::mousemode == "pick"} {
         Molefacture::editElement
      }
      %W selection clear
   }

   grid [ttk::frame $framecolapse.oxistate] -row 1 -column 0 -sticky we -padx 2 -pady 4
   grid columnconfigure $framecolapse.oxistate 0 -weight 0
   grid columnconfigure $framecolapse.oxistate 1 -weight 1
   grid rowconfigure $framecolapse.oxistate 0 -weight 1


   grid [ttk::label $framecolapse.oxistate.lbl -text "Oxidation State  " ] -row 0 -column 0 -sticky we -padx 2 -pady 2

   grid [ttk::combobox $framecolapse.oxistate.cmb -state normal -justify left \
   -textvariable Molefacture::cmboxstat -width 4 -style TblEdit.TCombobox -postcommand {
         $Molefacture::atomTable cancelediting
      }] -row 0 -column 1 -sticky w -padx 2


   set Molefacture::cmboxstatwidget $framecolapse.oxistate.cmb
   set Molefacture::cmboxstat ""

   bind $framecolapse.oxistate.cmb <<ComboboxSelected>> {
      Molefacture::oxiStateComboUpdate %W

      %W selection clear
   }

   bind $framecolapse.oxistate.cmb <Return> {
      Molefacture::oxiStateComboUpdate %W
   }

   grid [ttk::frame $framecolapse.editframe] -row 2 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.editframe 0 -weight 1
   grid columnconfigure $framecolapse.editframe 1 -weight 1
   grid rowconfigure $framecolapse.editframe 0 -weight 1


   grid [ttk::frame $framecolapse.bonddist] -row 3 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.bonddist 0 -weight 1
   grid rowconfigure $framecolapse.bonddist 0 -weight 1

   grid [ttk::frame $framecolapse.bonddist.labls] -row 0 -column 0 -sticky we
   grid columnconfigure $framecolapse.bonddist.labls 0 -weight 1

   grid rowconfigure $framecolapse.bonddist.labls 0 -weight 1

   grid [ttk::label $framecolapse.bonddist.labls.lbl -text "Bond Distance (A)"] -row 0 -column 0 -sticky we -pady 2

   grid [ttk::frame $framecolapse.bonddist.scaleframe] -row 1 -column 0 -sticky we
   grid columnconfigure $framecolapse.bonddist.scaleframe 0 -weight 1
   grid columnconfigure $framecolapse.bonddist.scaleframe 1 -weight 0
   grid columnconfigure $framecolapse.bonddist.scaleframe 2 -weight 0
   grid rowconfigure $framecolapse.bonddist.scaleframe 0 -weight 1

   grid [scale $framecolapse.bonddist.scaleframe.scale -from 0.5 -to 5.0 -resolution 0.01 \
   -showvalue false -length 100 -orient horizontal -variable Molefacture::bondguiValscale \
   -state disabled -command ::Molefacture::adjust_bondlength] -row 0 -column 0 -sticky we -pady 2


   grid [ttk::entry $framecolapse.bonddist.scaleframe.entryval -width 6 -state disabled \
   -validate none -textvariable Molefacture::bondguiValentry ] -row 0 -column 1 -sticky we -pady 2 -padx 2

   #-validatecommand Molefacture::validateBondEntry
   bind $framecolapse.bonddist.scaleframe.entryval <Return> {
      Molefacture::validateBondEntry
      Molefacture::saveState
      ::Molefacture::adjust_bondlength $Molefacture::bondguiValentry
   }

   bind $framecolapse.bonddist.scaleframe.scale  <ButtonPress> {
      Molefacture::saveState
   }

   set Molefacture::bondscalewdgt $framecolapse.bonddist.scaleframe.scale
   set Molefacture::bondentrywdgt $framecolapse.bonddist.scaleframe.entryval


   ############## Edit Angles #################

   grid [ttk::frame $cmdframe.angles] -row $mainrow -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $cmdframe.angles 0 -weight 1
   grid rowconfigure $cmdframe.angles 0 -weight 1

   incr mainrow
   grid [ttk::label $cmdframe.angles.prt -text "Angles" -image $Molefacture::arrowRight -compound left] -row 0 -column 0 -sticky w -pady 1

   bind $cmdframe.angles.prt <Button-1> {
      ::Molefacture::hideFrame %W [lindex [grid info %W] 1] "Angles"
   }

   grid [ttk::frame $cmdframe.angles.fcolapse] -row 1 -column 0 -sticky nsew -pady 5
   grid columnconfigure $cmdframe.angles.fcolapse 0 -weight 1


   set framecolapse $cmdframe.angles.fcolapse

   grid [ttk::frame $framecolapse.editang] -row 0 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.editang 0 -weight 1
   grid rowconfigure $framecolapse.editang 0 -weight 1


   grid [ttk::frame $framecolapse.editang.grps] -row 0 -column 0 -sticky we -padx 2
   grid columnconfigure $framecolapse.editang.grps 1 -weight 1
   grid columnconfigure $framecolapse.editang.grps 2 -weight 1
   grid rowconfigure $framecolapse.editang.grps 0 -weight 1

   grid [ttk::label $framecolapse.editang.grps.lbl -text "Move: "] -row 0 -column 0 -sticky w -pady 2

   ttk::style map group1.TCheckbutton -indicatorcolor [list selected green]
   ttk::style map group2.TCheckbutton -indicatorcolor [list selected blue]

   grid [ttk::checkbutton $framecolapse.editang.grps.grp1 -text "Group 1" -style group1.TCheckbutton -state disabled -variable ::Molefacture::anglemvgrp1 -command ::Molefacture::predraw_angleatoms] -row 0 -column 1 -sticky we -pady 2 -padx 2
   grid [ttk::checkbutton $framecolapse.editang.grps.grp2 -text "Group 2" -style group2.TCheckbutton -state disabled -variable ::Molefacture::anglemvgrp2 -command ::Molefacture::predraw_angleatoms] -row 0 -column 2 -sticky we -pady 2 -padx 2

   set Molefacture::anglemvgrp1wdgt $framecolapse.editang.grps.grp1
   set Molefacture::anglemvgrp2wdgt $framecolapse.editang.grps.grp2
   grid [ttk::frame $framecolapse.editang.scllbframe] -row 1 -column 0 -sticky we -padx 2
   grid columnconfigure $framecolapse.editang.scllbframe 0 -weight 1
   grid rowconfigure $framecolapse.editang.scllbframe 0 -weight 1

   grid [ttk::label $framecolapse.editang.scllbframe.lbl -text "Angle Value (deg)"] -row 0 -column 0 -sticky we -pady 2

   grid [ttk::frame $framecolapse.editang.sclframe] -row 2 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.editang.sclframe 0 -weight 1
   grid rowconfigure $framecolapse.editang.sclframe 0 -weight 1

   grid [scale $framecolapse.editang.sclframe.scale -from 0.1 -to 360.0 -resolution 0.1 -showvalue false -length 100 -orient horizontal -variable Molefacture::angleguiValscale -state disabled -command ::Molefacture::validateAngleScale ] -row 0 -column 0 -sticky we -pady 2
   grid [ttk::entry $framecolapse.editang.sclframe.entryval -width 6 -textvariable ::Molefacture::angleguiValentry -state disabled] -row 0 -column 1 -sticky e -pady 2 -padx 2

   set Molefacture::anglescalewdgt $framecolapse.editang.sclframe.scale
   set Molefacture::angleentrywdgt $framecolapse.editang.sclframe.entryval

   bind $framecolapse.editang.sclframe.entryval <Return> {
      Molefacture::validateAngleEntry
      Molefacture::saveState
      ::Molefacture::resize_angle $Molefacture::angleguiValentry
   }

   bind $framecolapse.editang.sclframe.scale  <ButtonPress> {
      Molefacture::saveState
   }

   ############## Edit Dihedral Angles#################

   grid [ttk::frame $cmdframe.dihedral] -row $mainrow -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $cmdframe.dihedral 0 -weight 1
   grid rowconfigure $cmdframe.dihedral 0 -weight 1

   incr mainrow

   grid [ttk::label $cmdframe.dihedral.prt -text "Dihedral Angles" -image $Molefacture::arrowRight -compound left] -row 0 -column 0 -sticky w -pady 1

   bind $cmdframe.dihedral.prt <Button-1> {
      ::Molefacture::hideFrame %W [lindex [grid info %W] 1] "Dihedral Angles"
   }

   grid [ttk::frame $cmdframe.dihedral.fcolapse] -row 1 -column 0 -sticky nsew -pady 5
   grid columnconfigure $cmdframe.dihedral.fcolapse 0 -weight 1


   set framecolapse $cmdframe.dihedral.fcolapse

   grid [ttk::frame $framecolapse.editdih] -row 0 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.editdih 0 -weight 1
   grid rowconfigure $framecolapse.editdih 0 -weight 1

   grid [ttk::frame $framecolapse.editdih.grps] -row 0 -column 0 -sticky we -padx 2
   grid columnconfigure $framecolapse.editdih.grps 0 -weight 1
   grid rowconfigure $framecolapse.editdih.grps 0 -weight 1

   grid [ttk::label $framecolapse.editdih.grps.lbl -text "Move: "] -row 0 -column 0 -sticky w -pady 2

   grid [ttk::radiobutton $framecolapse.editdih.grps.grp1 -text "Group 1" -value "group1" -variable Molefacture::dihedrgrp -state disabled -style group1.TCheckbutton -command ::Molefacture::predraw_angleatoms] -row 0 -column 1 -sticky w -pady 2 -padx 2
   grid [ttk::radiobutton $framecolapse.editdih.grps.grp2 -text "Group 2" -value "group2" -variable Molefacture::dihedrgrp -state disabled -style group2.TCheckbutton -command ::Molefacture::predraw_angleatoms] -row 0 -column 2 -sticky w -pady 2 -padx 2

   set Molefacture::dihradbtt1 $framecolapse.editdih.grps.grp1
   set Molefacture::dihradbtt2 $framecolapse.editdih.grps.grp2

   grid [ttk::frame $framecolapse.editdih.scllbframe] -row 1 -column 0 -sticky we -padx 2
   grid columnconfigure $framecolapse.editdih.scllbframe 0 -weight 1
   grid rowconfigure $framecolapse.editdih.scllbframe 0 -weight 1

   grid [ttk::label $framecolapse.editdih.scllbframe.lbl -text "Angle Value (deg)"] -row 0 -column 0 -sticky we -pady 2

   grid [ttk::frame $framecolapse.editdih.sclframe] -row 2 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.editdih.sclframe 0 -weight 1
   grid rowconfigure $framecolapse.editdih.sclframe 0 -weight 1

   grid [scale $framecolapse.editdih.sclframe.scale -from 0.1 -to 360.0 -resolution 0.1 -showvalue false -length 100 -orient horizontal -variable Molefacture::dihguiValscale -state disabled -command ::Molefacture::validateAngleScale] -row 0 -column 0 -sticky we -pady 2
   grid [ttk::entry $framecolapse.editdih.sclframe.entryval -width 6 -textvariable ::Molefacture::dihguiValentry -state disabled] -row 0 -column 1 -sticky e -pady 2 -padx 2

   set Molefacture::dihscalewdgt $framecolapse.editdih.sclframe.scale
   set Molefacture::dihentrywdgt $framecolapse.editdih.sclframe.entryval

   bind $framecolapse.editdih.sclframe.entryval <Return> {
      Molefacture::validateAngleEntry
      Molefacture::saveState
      ::Molefacture::rotate_bond $Molefacture::dihguiValentry
   }

   bind $framecolapse.editdih.sclframe.scale <ButtonPress> {
      Molefacture::saveState
   }

   ############## Edit Molecule#################

   grid [ttk::frame $cmdframe.molecule] -row $mainrow -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $cmdframe.molecule 0 -weight 1
   grid rowconfigure $cmdframe.molecule 0 -weight 1

   incr mainrow
   grid [ttk::label $cmdframe.molecule.prt -text "Edit Molecule" -image $Molefacture::arrowRight -compound left] -row 0 -column 0 -sticky w -pady 1

   bind $cmdframe.molecule.prt <Button-1> {
      ::Molefacture::hideFrame %W [lindex [grid info %W] 1] "Edit Molecule"
   }

   grid [ttk::frame $cmdframe.molecule.fcolapse] -row 1 -column 0 -sticky nsew -pady 5
   grid columnconfigure $cmdframe.molecule.fcolapse 0 -weight 1


   set framecolapse $cmdframe.molecule.fcolapse

   grid [ttk::frame $framecolapse.editmol] -row 0 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.editmol 0 -weight 1
   grid rowconfigure $framecolapse.editmol 0 -weight 1

   grid [ttk::frame $framecolapse.editmol.reschain] -row 0 -column 0 -sticky we -padx 2 -pady 2
   # grid columnconfigure $framecolapse.editmol.resseg 1 -weight 1
   grid columnconfigure $framecolapse.editmol.reschain 2 -weight 1
   grid rowconfigure $framecolapse.editmol.reschain 0 -weight 1

   grid [ttk::label $framecolapse.editmol.reschain.reslbl -text "Resname:" -width 9] -row 0 -column 0 -sticky w -pady 2
   grid [ttk::entry $framecolapse.editmol.reschain.resname -width 6 -state disabled -textvariable Molefacture::resname] -row 0 -column 1 -sticky we -pady 2

   grid [ttk::frame $framecolapse.editmol.reschain.tmp] -row 0 -column 2 -sticky we -padx 2 -pady 2

   grid [ttk::label $framecolapse.editmol.reschain.chainlbl -text "Chain:" -width 8] -row 0 -column 3 -sticky w -pady 2
   grid [ttk::entry $framecolapse.editmol.reschain.chain -width 6 -state disabled -textvariable Molefacture::chain] -row 0 -column 4 -sticky we -pady 2


   grid [ttk::frame $framecolapse.editmol.residchrg] -row 1 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.editmol.residchrg 2 -weight 1
   grid rowconfigure $framecolapse.editmol.residchrg 0 -weight 1


   grid [ttk::label $framecolapse.editmol.residchrg.reslbl -text "Resid:" -width 9 ] -row 0 -column 0 -sticky w -pady 2
   grid [ttk::entry $framecolapse.editmol.residchrg.residentry -width 6 -state disabled -textvariable Molefacture::resid] -row 0 -column 1 -sticky we -pady 2


   grid [ttk::frame $framecolapse.editmol.segc] -row 2 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.editmol.segc 2 -weight 1
   grid rowconfigure $framecolapse.editmol.segc 0 -weight 1

   grid [ttk::label $framecolapse.editmol.segc.seglbl -text "Segname:" -width 9] -row 0 -column 0 -sticky w -pady 2
   grid [ttk::entry $framecolapse.editmol.segc.segname -width 6 -state disabled -textvariable Molefacture::segname] -row 0 -column 1 -sticky we -pady 2

   grid [ttk::frame $framecolapse.editmol.segc.tmp] -row 0 -column 2 -sticky we -padx 2 -pady 2


   grid [ttk::frame $framecolapse.apply] -row 1 -column 0 -sticky we -padx 2
   grid columnconfigure $framecolapse.apply 0 -weight 1
   grid rowconfigure $framecolapse.apply 0 -weight 1

   grid [ttk::button $framecolapse.apply.applybtt -text "Apply" -padding "2 0 2 0" -command Molefacture::apply_moledits] -row 0 -column 0 -sticky e


   set Molefacture::rnamewdgt $framecolapse.editmol.reschain.resname
   set Molefacture::chnwdgt $framecolapse.editmol.reschain.chain
   set Molefacture::rsidwdgt $framecolapse.editmol.residchrg.residentry
   set Molefacture::sgnmwdgt $framecolapse.editmol.segc.segname

   ### Minimize structure

   grid [ttk::frame $cmdframe.minimize] -row $mainrow -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $cmdframe.minimize 0 -weight 1
   grid rowconfigure $cmdframe.minimize 0 -weight 1

   incr mainrow

   grid [ttk::label $cmdframe.minimize.prt -text "Minimize Structure" -image $Molefacture::arrowRight \
   -compound left] -row 0 -column 0 -sticky w -pady 1

   bind $cmdframe.minimize.prt <Button-1> {
      ::Molefacture::hideFrame %W [lindex [grid info %W] 1] "Minimize Structure"
   }

   grid [ttk::frame $cmdframe.minimize.fcolapse] -row 1 -column 0 -sticky nsew -pady 5
   grid columnconfigure $cmdframe.minimize.fcolapse 0 -weight 1


   set framecolapse $cmdframe.minimize.fcolapse

   grid [ttk::frame $framecolapse.minstrct] -row 0 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.minstrct 0 -weight 1
   grid rowconfigure $framecolapse.minstrct 0 -weight 1

   grid [ttk::frame $framecolapse.minstrct.selff] -row 0 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.minstrct.selff 1 -weight 1
   grid rowconfigure $framecolapse.minstrct.selff 0 -weight 1

   grid [ttk::label $framecolapse.minstrct.selff.lbl -text "Force Field" -width 10] -row 0 -column 0 -sticky w -padx 2 -pady 2

   grid [ttk::combobox $framecolapse.minstrct.selff.cmb -values {UFF CHARMM36} -textvariable Molefacture::minff \
   -width 12 -state readonly -style TblEdit.TCombobox] -row 0 -column 1 -sticky w -padx 2 -pady 2

   bind $framecolapse.minstrct.selff.cmb <<ComboboxSelected>> {
      %W selection clear
   }

    set Molefacture::minff "UFF"

   grid [ttk::frame $framecolapse.minstrct.labelbtt] -row 1 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.minstrct.labelbtt 0 -weight 1
   grid rowconfigure $framecolapse.minstrct.labelbtt 0 -weight 1

   set font [$framecolapse.minstrct.selff.lbl cget -font]
   set fontname [lindex $font 0]
   set fontsize [expr int([lindex $font 1] *.9)]


   grid [ttk::label $framecolapse.minstrct.labelbtt.label -text "CHARMM36 - NAMD\nUFF - OpenBabel" \
   -font "$fontname $fontsize" -wraplength 150 -padding "4 2 0 2"] -row 0 -column 0 -sticky new -pady 0 -padx 2


   grid [ttk::frame $framecolapse.minstrct.step] -row 2 -column 0 -sticky we -padx 2 -pady 2
   grid columnconfigure $framecolapse.minstrct.step 1 -weight 1
   grid rowconfigure $framecolapse.minstrct.step 0 -weight 1

   grid [ttk::label $framecolapse.minstrct.step.lbl -text "Num Steps" -width 10] -row 0 -column 0 -sticky w -padx 2 -pady 2
   grid [ttk::entry $framecolapse.minstrct.step.cmb -textvariable Molefacture::minsteps -width 9 \
   -validate focusout -validatecommand Molefacture::validateMinimizeSteps -justify right] -row 0 -column 1 -sticky w -padx 1 -pady 2

   bind $framecolapse.minstrct.step.cmb <Return> {
      Molefacture::validateMinimizeSteps
   }

   grid [ttk::button $framecolapse.minstrct.step.min -text "Minimize" -padding "2 0 2 0" -command Molefacture::minimize] -row 0 -column 2 -sticky e
   
   
   grid [ttk::frame $framecolapse.minstrct.cgenffsub] -row 3 -column 0 -sticky w -padx 2 -pady 2
   grid columnconfigure $framecolapse.minstrct.cgenffsub 0 -weight 0
    grid columnconfigure $framecolapse.minstrct.cgenffsub 1 -weight 0
     grid columnconfigure $framecolapse.minstrct.cgenffsub 2 -weight 0
   grid rowconfigure $framecolapse.minstrct.cgenffsub 0 -weight 1
   
   grid [ttk::label $framecolapse.minstrct.cgenffsub.lbl -text "CGenFF\nSubmissions" -width 12] -row 0 -column 0 -sticky w -padx 2 -pady 2
   grid [ttk::label $framecolapse.minstrct.cgenffsub.val -textvariable Molefacture::cgenffsub -width 3 -anchor e] -row 0 -column 2 -sticky e -padx 0 -pady 2
   grid [ttk::label $framecolapse.minstrct.cgenffsub.total -text "/100" -width 4 -anchor w] -row 0 -column 3 -sticky w -padx 0 -pady 2

   ### Collapse frames

   grid forget $cmdframe.angles.fcolapse
   grid forget $cmdframe.dihedral.fcolapse
   grid forget $cmdframe.molecule.fcolapse
   grid forget $cmdframe.minimize.fcolapse



   ### General molecule information

   grid [ttk::frame $cmdframe.geninfo] -row $mainrow -column 0 -sticky w -padx 2 -pady 2
   # grid columnconfigure $cmdframe.geninfo 0 -weight 1
   # grid columnconfigure $cmdframe.geninfo 0 -weight 1
   grid rowconfigure $cmdframe.geninfo 0 -weight 1

   incr mainrow

   grid [ttk::frame $cmdframe.geninfo.numatom] -row 0 -column 0 -sticky w -padx 2 -pady 2

   grid [ttk::label $cmdframe.geninfo.numatom.natomlbl -text "Num\nAtoms:" -width 6] -row 0 -column 0 -sticky w -pady 2
   grid [ttk::label $cmdframe.geninfo.numatom.natomval -width 6 -textvariable Molefacture::numAtoms -justify left] -row 0 -column 4 -sticky w -pady 2



   grid [ttk::frame $cmdframe.geninfo.charge] -row 0 -column 1 -sticky w -padx 0 -pady 2

   grid [ttk::label $cmdframe.geninfo.charge.chrglbl -text "Total\nCharge:" -width 8] -row 0 -column 0 -sticky w -pady 2
   grid [ttk::label $cmdframe.geninfo.charge.chargeval -width 6 -textvariable Molefacture::charge -justify left] -row 0 -column 4 -sticky w -pady 2

   trace add variable ::Molefacture::topocgenff write ::Molefacture::moleculeChanged

   raise $topGui
   return $topGui
}
######################################################
# Proc triggered by the "CGenFF Server" Menu item in #
# the Build Menu                                     #
######################################################
proc Molefacture::CgenffServerCall {} {
   global env
  if {$Molefacture::topocgenff == 1} {
    tk_messageBox -message "The most updated topology is already assigned to the molecule." -icon info -type ok -parent $Molefacture::topGui
    return
  }
  if {[$Molefacture::atomTable size] ==0} {
    return
  }

  set sel [atomselect $Molefacture::tmpmolid "$Molefacture::notDumSel"]
  if {[Molefacture::findHalogens $sel] == 1} {
    tk_messageBox -message "It is not currently possible to generate psf and\
    pdb files from a CHARMM topology containing halogen elements." -title "Atom typing failure" -icon error -type ok \
    -parent $Molefacture::topGui
    $sel delete
    return
  }
  set unknownRes ""

  set resname [lsort -unique -dictionary [$Molefacture::atomTable getcolumns Resname]]

  if {$Molefacture::topocgenff != 1} {
    set resultRes [Molefacture::loadDefaulTopo $resname]
    set unknownRes [lindex $resultRes 0]
    if {[llength $unknownRes] == 0} {
      Molefacture::add_all_hydrogen
      # Molefacture::updateprotlasth
      return
    }
  }

  set pwd [pwd]

  cd $env(TMPDIR)
  if {[file exist molefacture_tmp] != 1} {
    file mkdir molefacture_tmp
  }

  cd molefacture_tmp

  ::Molefacture::submitCGENFF $unknownRes

  cd $pwd
}
#######################################################
# Proc triggered when the oxidation state combo value #
# is changed                                          #
#######################################################
proc Molefacture::oxiStateComboUpdate { W } {
  $Molefacture::cmboxstatwidget set $Molefacture::cmboxstat

  set oxlist [$W cget -values]
  set oxstat [lsearch $oxlist $Molefacture::cmboxstat]
  if {$oxstat == -1} {
     set list [lappend oxlist $Molefacture::cmboxstat]
     # set list [lsort -unique -integer -increasing $list]

     $W configure -values $list
     set oxlist [%W cget -values]

     set oxstat [lsearch $oxlist $Molefacture::cmboxstat]
  }
  $W current $oxstat
  ::Molefacture::updateOxiState

}



######################################################
# Set keyborad hot keys to be used on the OpenGL     #
# window                                             #
######################################################
proc Molefacture::setHotKeys {} {

   ### link the hot keys with the mouse modes in the action frame
   ### Alt key doesn't work on Mac

   user add key Control-s {set Molefacture::mousemode pick}
   user add key Control-m {set Molefacture::mousemode movefrag}
   user add key Control-a {set Molefacture::mousemode addatom}
   user add key A {set Molefacture::mousemode delatom}
   user add key Control-b {set Molefacture::mousemode addbond}
   user add key B {set Molefacture::mousemode delbond}
   user add key Control-i {set Molefacture::mousemode incrorder}
   user add key I {set Molefacture::mousemode decrorder}

   ### link the hot keys with the buttons in the action frame

   user add key Control-h {Molefacture::addHMain}
   user add key Control-d {Molefacture::duplicateSel}

   user add key O {Molefacture::setElemComb "Oxygen(O:8)"}
   user add key C {Molefacture::setElemComb "Carbon(C:6)"}
   user add key H {Molefacture::setElemComb "Hydrogen(H:1)"}
   user add key S {Molefacture::setElemComb "Sulfur(S:16)"}
   user add key N {Molefacture::setElemComb "Nitrogen(N:7)"}
   user add key P {Molefacture::setElemComb "Phosphorus(P:15)"}

   user add key Control-c {Molefacture::clearTableSelection}
   user add key Control-z {Molefacture::undoProc}
}
######################################################
# Set element combobox text for the hot keys         #
# definition                                         #
######################################################
proc Molefacture::setElemComb { text } {
   set Molefacture::elemcomb $text
   Molefacture::editElement
}
######################################################
# Event triggered by the element combobox            #
######################################################
proc Molefacture::editElement {} {
      set cmbtext $Molefacture::elemcomb
      ### get table current selection

      set index [::Molefacture::keepTableSelection]
     
      if {$cmbtext == "Other..."} {

         ::Molefacture::periodicTableWindow


         if {$index != ""} {
            set element [$Molefacture::atomTable cellcget [lindex $index 0],Element -text ]
            set atmnum [lsearch $Molefacture::periodic $element]
            ::Molefacture::elementPTBind $atmnum
         }

      } elseif {$index != ""} {
         Molefacture::saveState

         set atmnum [string trimright [lindex [split $cmbtext ":"] 1] ")"]
         foreach tbind $index {
            ::Molefacture::edit_atom_elem [$Molefacture::atomTable cellcget $tbind,Index -text ] $atmnum
         }

         ::Molefacture::clearTableSelection
      }
       if {$Molefacture::topocgenff == 1} {
         set Molefacture::topocgenff 0
      }
      update_openvalence
}
######################################################
# Update the "Add Template" frame                    #
######################################################
proc Molefacture::updateTempSel {} {

   if {$Molefacture::cyclSelector == "cycle"} {
      grid forget [lindex $Molefacture::cyclSelecWidget 1]
      grid conf [lindex $Molefacture::cyclSelecWidget 0]  -row 0 -column 1 -padx 2 -sticky w
   } else {
      grid forget [lindex $Molefacture::cyclSelecWidget 0]
      grid conf [lindex $Molefacture::cyclSelecWidget 1]  -row 0 -column 1 -padx 2 -sticky w
   }
}

######################################################
# Event triggered from the add atom button           #
######################################################
proc Molefacture::addAtomMain {} {


   set atmnum [string trimright [lindex [split $Molefacture::elemcomb ":"] 1] ")"]
   set element [lindex $Molefacture::periodic $atmnum]
   set newatoms [Molefacture::add_atom_gui $element]
   Molefacture::keepTableSelection
   mol reanalyze $Molefacture::tmpmolid
   Molefacture::update_openvalence

   return $newatoms
}

######################################################
# Event triggered from the delete atom button        #
######################################################
proc Molefacture::delAtomMain {} {
   Molefacture::del_atom_gui
   Molefacture::clearTableSelection
   mol reanalyze $Molefacture::tmpmolid
}

######################################################
# Event triggered from the +H  button                #
######################################################
proc Molefacture::addHMain {} {
  variable picklist
   if {[molinfo top] == -1} {
      return
    }
   
   ### keep initial selection
   set auxpicklist $picklist
   trace remove variable ::Molefacture::mousemode write ::Molefacture::setMouseMode
   set prevmousemode $Molefacture::mousemode

   set Molefacture::mousemode "pick"
   ::Molefacture::add_all_hydrogen

   ::Molefacture::keepTableSelection
   mol reanalyze $Molefacture::tmpmolid
   set Molefacture::topocgenff 0

   set Molefacture::mousemode $prevmousemode
   trace add variable ::Molefacture::mousemode write ::Molefacture::setMouseMode
   
   ### Reassign atom selection
   if {[llength $auxpicklist] > 0} {
     set picklist $auxpicklist 
   } else {
     ::Molefacture::clearTableSelection
   }
}


######################################################
# Duplicate selected atoms                           #
######################################################

proc Molefacture::duplicateSel {} {
   set newatoms [::Molefacture::new_independent_fragment "" 1]
   ::Molefacture::updateNewAtoms [join $newatoms]
   Molefacture::clearTableSelection
}


######################################################
# Event invoked when Molefacture window is           #
# closed                                             #
######################################################
proc ::Molefacture::closeMainWindow {} {

   if {[lsearch [molinfo list] $Molefacture::tmpmolid] != -1} {
      set answer [tk_messageBox -message "Closing Molefacture will delete current editing molecule. Do you want to close?" -type yesno -title "Closing Molefacture" -icon info -parent $Molefacture::topGui]
      if {$answer == "no"} {
         return
      }
   }
  wm withdraw $Molefacture::topGui

  ### Removing trace events
  trace remove variable ::Molefacture::topocgenff write ::Molefacture::moleculeChanged

  ### Remove the binding of selection lost on the main table
  bind $Molefacture::atomTable <<TablelistSelectionLost>> {break}

  mouse mode rotate
  if {[$Molefacture::atomTable size] > 0} {
   Molefacture::clearTableSelection
   $Molefacture::atomTable delete 0 end
  }

  Molefacture::moleculeChanged
  if {[lsearch [molinfo list] $Molefacture::tmpmolid] != -1} {
   mol delete $Molefacture::tmpmolid
  }
  set Molefacture::tmpmolid -1
  foreach m [molinfo list] {
    if {$Molefacture::qwikmd == 1 && $Molefacture::topocgenff == 1 && $m == $QWIKMD::topMol} {
      molinfo $m set drawn 0
      continue
    }
    molinfo $m set drawn 1
  }
  if {$Molefacture::qwikmd == 0 } {
   set QWIKMD::topoinfo ""
  }
  display resetview
}
######################################################
# Binding Proc invoked when the Atom's table loose   #
# selection                                          #
######################################################
proc ::Molefacture::bindAtomTable {} {
    if {[lsearch [molinfo list] $Molefacture::tmpmolid] == -1} {
      return
    }
   if {[$Molefacture::atomTable size] > 0} {
      set index [::Molefacture::keepTableSelection]
      if {$index == ""} {
         ::Molefacture::clearTableSelection
      }
   }
}

######################################################
# Window to set the username to access the CGenff    #
# websever                                           #
######################################################
proc ::Molefacture::cgenffUserWindows {} {
  global env
  variable topGui
  variable cgenffWindow

  if { [winfo exists $cgenffWindow] ==1} {
    wm deiconify $cgenffWindow
    focus $cgenffWindow
    raise $cgenffWindow
    return
  }


   set w [toplevel $cgenffWindow]
   wm title $w "CGenFF Settings"
   wm resizable $w 1 1

   grid columnconfigure $cgenffWindow 0 -weight 1
   grid rowconfigure $cgenffWindow 1 -weight 1

   variable atomsel

   set atomsel ""
   
   grid [ttk::frame $cgenffWindow.mainframe] -row 0 -column 0 -sticky nwes -padx 4 -pady 4
   grid columnconfigure $cgenffWindow.mainframe 0 -weight 1
   grid rowconfigure $cgenffWindow.mainframe 0 -weight 0
   grid rowconfigure $cgenffWindow.mainframe 1 -weight 1

   set mainrow 0
   
   grid [ttk::frame $cgenffWindow.mainframe.disc1] -row $mainrow -column 0 -sticky we
   grid columnconfigure $cgenffWindow.mainframe.disc1 0 -weight 1
   grid rowconfigure $cgenffWindow.mainframe.disc1 0 -weight 0
   incr mainrow
   
   set text1 "Privacy information: Your activities on \
   cgenff.umaryland.edu may be monitored for the purpose of gathering usage statistics, and \
   the molecules you upload may be used for improving the quality of the CGenFF program and \
   force field."
   
   set text2 "Non-Commercial Use: The access to the CGenFF server should be used for non-commercial purpose only. \
   For more information, please contact: support@silcsbio.com."
   
   grid [text $cgenffWindow.mainframe.disc1.txt -width 60 -wrap word -relief flat -bg white -exportselection yes -pady 2 -padx 2 -height 10] -row 0 -column 0 
   $cgenffWindow.mainframe.disc1.txt insert end $text1
    $cgenffWindow.mainframe.disc1.txt tag add title 1.0 1.19;# "Privacy information - 19 characters"
   $cgenffWindow.mainframe.disc1.txt tag configure title -font "helvetica 12 bold"
   
   
   $cgenffWindow.mainframe.disc1.txt insert end "\n\n"

   $cgenffWindow.mainframe.disc1.txt insert end $text2
   
   set ini [$cgenffWindow.mainframe.disc1.txt search -exact $text2 1.0 end]
   set line [split $ini "."]

   $cgenffWindow.mainframe.disc1.txt tag add title2 $ini [lindex $line 0].18;# "Non-Commercial Use - 18 characters"
   $cgenffWindow.mainframe.disc1.txt tag configure title2 -font "helvetica 12 bold"
   
   $cgenffWindow.mainframe.disc1.txt configure -state disabled
   
   ### Use CGenFF locally
   grid [ttk::frame $cgenffWindow.mainframe.runlocal] -row $mainrow -column 0 -sticky nwes
   grid columnconfigure $cgenffWindow.mainframe.runlocal 0 -weight 1
   grid rowconfigure $cgenffWindow.mainframe.runlocal 0 -weight 1
   incr mainrow  
    
   grid [ttk::radiobutton $cgenffWindow.mainframe.runlocal.run -text "Run CGenFF locally" \
   -variable env(CGENFF) -value "local" -command { 
     Molefacture::checkCgenffWindow
     } -state normal] -row 0 -column 0 -sticky nwes  
   
   grid [ttk::frame $cgenffWindow.mainframe.cgenffpath] -row $mainrow -column 0 -sticky nwes -padx 2 -pady 3
   grid columnconfigure $cgenffWindow.mainframe.cgenffpath 1 -weight 1
   grid rowconfigure $cgenffWindow.mainframe.cgenffpath 0 -weight 1
   incr mainrow  
  
   
   
   grid [ttk::label $cgenffWindow.mainframe.cgenffpath.lbl -text "CGenFF Path"] -row 0 -column 0 -padx 2 -pady 2 -sticky news
   
   grid [ttk::entry $cgenffWindow.mainframe.cgenffpath.entry -textvariable env(CGENFFPATH) -validate all -validatecommand {
     if {[string trim $env(CGENFFPATH)] != "" && \
     [file exists $env(CGENFFPATH)] == 1} {
       $Molefacture::cgenffbrowser configure -text "Apply"
     } else {
       $Molefacture::cgenffbrowser configure -text "Browser"
     }
     return 1
     } ] -row 0 -column 1 -padx 2 -pady 2 -sticky ew
   
  
   grid [ttk::button $cgenffWindow.mainframe.cgenffpath.btt -text "Browser" \
   -command Molefacture::cgenffBrwApply -padding "2 0 2 0"] -row 0 -column 2 -padx 2 -pady 2 -sticky news
   
   set Molefacture::cgenffpathval ""

   
   grid [ttk::frame $cgenffWindow.mainframe.lbl] -row $mainrow -column 0 -sticky nwes -pady 2 -padx 2
   grid columnconfigure $cgenffWindow.mainframe.lbl 0 -weight 1
   grid rowconfigure $cgenffWindow.mainframe.lbl 0 -weight 1
   incr mainrow  
   
   grid [ttk::radiobutton $cgenffWindow.mainframe.lbl.run -text "Run CGenFF Server" \
   -variable env(CGENFF) -value "server" -command Molefacture::checkCgenffWindow ] -row 0 -column 0 -sticky nwes  
   
   grid [ttk::label $cgenffWindow.mainframe.lbl.txt1 -text "Please enter the username of \
   your account on the the CGenFF webserver." -justify left] -row 1 -column 0 -sticky we -padx 0 -pady 2

  
   grid [ttk::frame $cgenffWindow.mainframe.lbl2] -row $mainrow -column 0 -sticky nwes
   grid columnconfigure $cgenffWindow.mainframe.lbl2 0 -weight 1
   grid rowconfigure $cgenffWindow.mainframe.lbl2 0 -weight 1
   incr mainrow  
  grid [ttk::label $cgenffWindow.mainframe.lbl2.text1 -text "For more information please visit "  -wraplength 250 -justify left] -row 1 -column 0 -sticky nws -padx 1 -pady 2

  grid [ttk::label $cgenffWindow.mainframe.lbl2.text2 -foreground blue -text "https://cgenff.umaryland.edu/" -anchor w -cursor hand1] -row 1 -column 1 -sticky nws -padx 0 -pady 2

   bind $cgenffWindow.mainframe.lbl2.text2 <Button-1> { vmd_open_url "https://cgenff.umaryland.edu/" }


   grid [ttk::frame $cgenffWindow.mainframe.entry] -row $mainrow -column 0 -sticky nwes
   grid columnconfigure $cgenffWindow.mainframe.entry 0 -weight 0
   grid columnconfigure $cgenffWindow.mainframe.entry 1 -weight 2
   grid rowconfigure $cgenffWindow.mainframe.entry 0 -weight 1
   incr mainrow
   
   grid [ttk::label $cgenffWindow.mainframe.entry.sel -text "Username: "] -row 0 -column 0 -sticky nwes -pady 2 -padx 2
   grid [ttk::entry $cgenffWindow.mainframe.entry.entry -textvariable env(CGENFFUSERNAME) -width 15] -row 0 -column 1 -sticky nswe -pady 2 -padx 2

   grid [ttk::label $cgenffWindow.mainframe.entry.selpss -text "Password: "] -row 0 -column 2 -sticky nwes -pady 2 -padx 2
   grid [ttk::entry $cgenffWindow.mainframe.entry.entrypass -textvariable env(CGENFFPASSWORD) -width 15] -row 0 -column 3 -sticky nswe -pady 2 -padx 2

   grid [ttk::button $cgenffWindow.mainframe.entry.go -text "Apply" -command Molefacture::checkCgenffUser] -row 0 -column 4 -sticky nwes -pady 2 -padx 2

   
   grid [ttk::frame $cgenffWindow.mainframe.vmdrc] -row $mainrow -column 0 -sticky nwes
   grid columnconfigure $cgenffWindow.mainframe.vmdrc 0 -weight 0
   grid columnconfigure $cgenffWindow.mainframe.vmdrc 1 -weight 2
   grid rowconfigure $cgenffWindow.mainframe.vmdrc 0 -weight 1
   incr mainrow
   
   grid [ttk::label $cgenffWindow.mainframe.vmdrc.sel -text "You can add the following line to your .vmdrc (or vmd.rc on Windows) file: " -wraplength 150 -justify left -padding "0 2 2 0"] -row 0 -column 0 -sticky nwes -pady 2 -padx 2
   
   grid [ttk::frame $cgenffWindow.mainframe.vmdrc.frmentry] -row 0 -column 1 -sticky we

   grid [ttk::entry $cgenffWindow.mainframe.vmdrc.frmentry.entry -textvariable Molefacture::cgenffWinVMDRC -state readonly ] -row 0 -column 0 -sticky we -pady 2 -padx 0
   grid [ttk::entry $cgenffWindow.mainframe.vmdrc.frmentry.entrypass -textvariable Molefacture::cgenffWinVMDRCPASS -state readonly ] -row 1 -column 0 -sticky we -pady 0 -padx 2
   grid forget $cgenffWindow.mainframe.vmdrc
   
   set Molefacture::cgenffvmdrc "$cgenffWindow.mainframe.vmdrc"
   set Molefacture::uernameentry "$cgenffWindow.mainframe.entry.entry"
   set Molefacture::cgenffpathentry "$cgenffWindow.mainframe.cgenffpath.entry"
   set Molefacture::cgenffbrowser "$cgenffWindow.mainframe.cgenffpath.btt"
   
   ### Initialize the env variables
   if {[info exist env(CGENFFPATH)] != 1} {
      set command [::ExecTool::find "cgenff"]
      if {[string length $command] != 0} {
         set env(CGENFFPATH) $command
         $Molefacture::cgenffbrowser configure -text "Apply"
      } else {
      	set env(CGENFFPATH) -1
      	$Molefacture::cgenffbrowser configure -text "Browser"
      }
   }
   
   if {$Molefacture::cgenffcaller == 0} {
      $cgenffWindow.mainframe.lbl.run configure -state disabled -text "Run CGenFF Server (Coming Soon)"
      set env(CGENFF) "local"
   }

   Molefacture::checkCgenffWindow
}
###############################
# Check the CGenFF username   #
###############################
proc Molefacture::checkCgenffUser {} {
  global env

  if {[winfo exists $Molefacture::cgenffWindow.mainframe.entry] == 1 && \
    $env(CGENFFUSERNAME) != -1 && \
     $env(CGENFFUSERNAME) != ""}  {
    set Molefacture::cgenffWinVMDRC "set env(CGENFFUSERNAME) \"$env(CGENFFUSERNAME)\""
    set Molefacture::cgenffWinVMDRCPASS "set env(CGENFFPASSWORD) \"$env(CGENFFPASSWORD)\""
    $Molefacture::cgenffWindow.mainframe.vmdrc.frmentry.entry configure -width [string length $Molefacture::cgenffWinVMDRC]
    $Molefacture::cgenffWindow.mainframe.vmdrc.frmentry.entrypass configure -width [string length $Molefacture::cgenffWinVMDRCPASS]  
    set info [grid info $Molefacture::cgenffWindow.mainframe.entry]
    grid conf $Molefacture::cgenffWindow.mainframe.vmdrc -row [expr [lindex $info 5] +2] -column 0 -sticky nwes
    return
  } elseif {[winfo exists $Molefacture::cgenffWindow.mainframe.entry] == 1 && \
   $env(CGENFFUSERNAME) == ""} {
    tk_messageBox -message "Something went wrong. Please try again." -title "CGenFF Username & Password" -icon error \
    -type ok -parent $Molefacture::cgenffWindow
  }
  set env(CGENFFUSERNAME) -1
}
#######################################################
# Disable and hide section of the CGenFF window       #
# of the depending the value of localcgenff           #
#######################################################
proc Molefacture::checkCgenffWindow {} {
   global env
  if {$env(CGENFF) == "local"} {
    grid forget $Molefacture::cgenffvmdrc 
    $Molefacture::uernameentry configure -state disabled
    $Molefacture::cgenffpathentry configure -state normal
    $Molefacture::cgenffbrowser configure -state normal
    
  } else {
    $Molefacture::uernameentry configure -state normal
    $Molefacture::cgenffpathentry configure -state disabled
    $Molefacture::cgenffbrowser configure -state disabled
    
    
  }
}

#######################################################
# Command triggered by the Browser/Apply command      #
#######################################################
proc Molefacture::cgenffBrwApply {} {
    global env
    set btttext [$Molefacture::cgenffbrowser cget -text]
    if {$btttext == "Apply"} {
      if {[string trim $env(CGENFFPATH)] == "-1" || [string trim $env(CGENFFPATH)] == "" || \
       [file exists $env(CGENFFPATH)] != 1} {
         tk_messageBox -message "Please specify the path to the CGenFF Executable." \
         -title "CGenFF Executable Path" -icon error \
            -type ok -parent $Molefacture::cgenffWindow
        set env(CGENFFPATH) -1
      } else {
         tk_messageBox -message "CGenFF Path set." -title "CGenFF Executable Path" -icon info \
            -type ok -parent $Molefacture::cgenffWindow
         $Molefacture::cgenffbrowser configure -text "Browser"
      }
    } else {
      set file [tk_getOpenFile -title "Select CGenFF Executable"]
      if {$file != ""} {
        set env(CGENFFPATH) $file
        $Molefacture::cgenffbrowser configure -text "Apply"
      }
    }
}
#######################################################
# Delete the temp folder whenever the molecule        #
# is edited and save a state to enable revert to that #
# state                                               #
#######################################################
proc ::Molefacture::moleculeChanged { args } {
   global env
   variable tmpmolid
   if {$Molefacture::topocgenff == 0 && [file exists $env(TMPDIR)/molefacture_tmp] == 1} {
      file delete -force $env(TMPDIR)/molefacture_tmp
   }
   if {[$Molefacture::atomTable size] > 0} {
      # Molefacture::saveState
   } else {
      set Molefacture::undo [list]
   }

   if {$tmpmolid == -1 || [lsearch [molinfo list] $tmpmolid] == -1} {
      return
   }
   ### remove the topology of the residues being edited from the list of topologies
   ### only do this when editing only one residue (most common scenario when editing ligands)
   set sel [atomselect $tmpmolid "($Molefacture::notDumSel) and not water and not ions and not protein and not lipid and not glycan"]
   set resnames [lsort -unique -dictionary [$sel get resname]]

   if {[llength $resnames] > 0 && $Molefacture::topocgenff == 0} {
      foreach residue $resnames {
        set index 0
        foreach topo $Molefacture::deftopinfo {
          if { [::Toporead::topology_contains_resi $topo $residue] != -1} {
            set Molefacture::deftopinfo [lreplace $Molefacture::deftopinfo $index $index]
            break
          }
          incr index
        }
        
      }
      
      foreach residue $resnames {
        set index 0
        
        foreach topo $Molefacture::topoinfo {
          if { [::Toporead::topology_contains_resi $topo $residue] != -1} {
            if {[file exists [lindex $Molefacture::topoinfo $index]]} {
                file delete -force [lindex $Molefacture::topoinfo $index]
            }
            set Molefacture::topoinfo [lreplace $Molefacture::topoinfo $index $index]
            break
          }
          incr index
        }
        
      }
   }
   $sel delete
}


#####################################################
# Proc set Mouse mode. The args is a requirement    #
# from the trace command                            #
#####################################################
proc ::Molefacture::setMouseMode { args } {
   global vmd_pick_atom
   variable picklist

   if {$Molefacture::mousemode  != $Molefacture::prevmousemode && \
      $Molefacture::mousemode != "rotate" && $Molefacture::mousemode != "translate" &&\
      $Molefacture::mousemode != "addbond" && $Molefacture::mousemode != "delbond"} {
      set picklist ""

      Molefacture::clearTableSelection
   }

   trace remove variable ::vmd_pick_event write ::Molefacture::atom_picked_fctn
   if {$Molefacture::mousemode == "pick"} {
      mouse mode pick
   } elseif {$Molefacture::mousemode == "movefrag"} {
      mouse mode movefrag
      # return
   } elseif {$Molefacture::mousemode == "moveatom"} {
      mouse mode moveatom
      # return
   } elseif {$Molefacture::mousemode == "rotate"} {
      mouse mode rotate
      # return
   } elseif {$Molefacture::mousemode == "translate"} {
      mouse mode translate
      # return
   } else {
      mouse mode pick
   }
   mouse callback on

    trace add variable ::vmd_pick_event write ::Molefacture::atom_picked_fctn
   

   set Molefacture::prevmousemode $Molefacture::mousemode
   

}

#####################################################
# Select atom on the Table                          #
#####################################################
proc ::Molefacture::selectAtomTable {} {
   variable tmpmolid
   if {[lsearch [molinfo list] $tmpmolid] == -1} {
     Molefacture::clearTableSelection
     return
   }
   set index [$Molefacture::atomTable curselection]
   if {$index != "" && [llength $index] == 1} {
      ::Molefacture::update_combo_sel_atoms $index
      if {[winfo exists $Molefacture::periTabWindow] == 1} {
         set element [$Molefacture::atomTable cellcget $index,Element -text ]
         set atmnum [lsearch $Molefacture::periodic [string toupper $element]]
         ::Molefacture::elementPTBind $atmnum
      }
   }
   set Molefacture::lastAtmSel $index
   set Molefacture::picklist {}
   set pickelement [list]
   foreach pickatm [$Molefacture::atomTable curselection] {
      lappend Molefacture::picklist [$Molefacture::atomTable cellcget $pickatm,Index -text]
      lappend pickelement [$Molefacture::atomTable cellcget $pickatm,Element -text]
   }
   

   ### Cancel oxidation state edit combobox if different elements
   ### are selected
   set pickelement [lsort -unique -dictionary $pickelement]
   set editcolumnind [lindex [$Molefacture::atomTable editinfo] end]

   if {[llength $pickelement] > 1 && $editcolumnind !=  -1\
      && [$Molefacture::atomTable columncget $editcolumnind -name] == "OxState"} {

      $Molefacture::atomTable cancelediting
   }

   ### update the graphics of the selected atoms
   Molefacture::draw_selatoms
   set resetangles 0
   set resetdih 0

   if {[llength $index] != 2} {
      ::Molefacture::disablebond
   }
   if {[llength $index] != 3} {
      ::Molefacture::disableangle
   }
   if {[llength $index] != 4} {
      ::Molefacture::disabledih
   }

   ### Enable the bond length changes Gui Events
   if {[llength $index] == 2} {

       if {[Molefacture::isbond] == 1} {
         $Molefacture::bondentrywdgt configure -state normal
         $Molefacture::bondscalewdgt configure -state normal
      }

   } elseif {[llength $index] == 3} {
      if {[Molefacture::predraw_angleatoms] == 1} {

         $Molefacture::anglemvgrp1wdgt configure -state normal
         $Molefacture::anglemvgrp2wdgt configure -state normal
         $Molefacture::angleentrywdgt configure -state normal
         $Molefacture::anglescalewdgt configure -state normal
      } else {
         set resetangles 1
      }

   } elseif {[llength $index] == 4} {
      if {[Molefacture::predraw_angleatoms] == 1} {

         $Molefacture::dihradbtt1 configure -state normal
         $Molefacture::dihradbtt2 configure -state normal
         $Molefacture::dihscalewdgt configure -state normal
         $Molefacture::dihentrywdgt configure -state normal

      } else {
         set resetdih 1
      }
   }


   if {$resetangles == 1} {
      Molefacture::disableangle
   } elseif {$resetdih == 1} {
      ::Molefacture::disabledih
   }

   if {[llength $Molefacture::picklist] > 0} {
    set Molefacture::numAtoms [llength $Molefacture::picklist]
  } else {
    set Molefacture::numAtoms [$Molefacture::atomTable size]
  }

   ::Molefacture::updateEditMolecule
   ::Molefacture::oxiStateCombo
}
#########################################################
# Evaluate if the selected atoms (two) are bonded       #
# and if any of the groups is selected, draw the atoms  #
#########################################################
proc ::Molefacture::isbond {} {
   variable tmpmolid
   variable picklist
   if {[llength $picklist] != 2} {
      return 0
   }

   set index1 [lindex $picklist 0]
   set index2 [lindex $picklist 1]

   ### the command getbonds return the bonds of the atoms order by index
   set selindex [lsort -integer [list $index1 $index2]]
   set sel1 [atomselect $tmpmolid "index [lindex $selindex 0]"]

   set bonds [lindex [$sel1 getbonds] 0]
   $sel1 delete
   if {[lsearch $bonds [lindex $selindex 1]] != -1} {
      set Molefacture::bondguiValscale [format %.2f [measure bond $selindex]]
      set Molefacture::bondguiValentry $Molefacture::bondguiValscale
      return 1
   } else {
      return 0
   }

}


#########################################################
# Evaluate if the selected atoms make an angle (bonded) #
# and if any of the groups is selected, draw the atoms  #
#########################################################
proc ::Molefacture::isanglebond {} {
   variable tmpmolid
   variable picklist
   
   ### Check if they are bonded and select the two different groups and represent the
   ### the angle atoms
   set index1 [lindex $picklist 0]
   set index2 [lindex $picklist 1]
   set index3 [lindex $picklist 2]

   set idx1sel [atomselect $tmpmolid "index $index1"]
   set idx2sel [atomselect $tmpmolid "index $index2"]
   set idx3sel [atomselect $tmpmolid "index $index3"]

   set bondlist [list [lindex [$idx1sel getbonds] 0] [lindex [$idx2sel getbonds] 0] [lindex [$idx3sel getbonds] 0]]
   set Molefacture::anglemother ""

   set ind 0
   set indlist "$index1 $index2 $index3"
   set motherind -1
   foreach bonds $bondlist {
      set indlistaux [lreplace $indlist $ind $ind]
      # puts "DEBUG is mother [lindex $indlist $ind]? ||  $bonds [lindex $indlistaux 0] || $bonds [lindex $indlistaux 1]"
      if {[lsearch $bonds [lindex $indlistaux 0]] != -1 && [lsearch $bonds [lindex $indlistaux 1]] != -1} {
         set Molefacture::anglemother [lindex $indlist $ind]
         set motherind $ind
         break
      }
      incr ind
   }
   $idx1sel delete
   $idx2sel delete
   $idx3sel delete
   if {$Molefacture::anglemother == ""} {
      return 0
   }

   ### sorting assuming the smallest index were added before when building the
   ### the molecule from scratch. No practical effect on the editing the molecule based
   ### on a selection from a VMD molecule
   set indlist [lsort -integer [lreplace $indlist $motherind $motherind]]


   ### When dividing the molecule in the mother atom,
   ### Set the group 2 the group of bonded atoms making the smallest bundle
   set molsize [$Molefacture::atomTable size]
   set indexes1 [join [::util::bondedsel $tmpmolid $Molefacture::anglemother [lindex $indlist 0] -maxdepth $molsize]]
   set indexes2 [join [::util::bondedsel $tmpmolid $Molefacture::anglemother [lindex $indlist 1] -maxdepth $molsize]]

   set Molefacture::anglemvgrp1list $indexes1
   set Molefacture::anglemvgrp2list $indexes2
   set ind1 [lindex $indlist 0]
   set ind2 [lindex $indlist 1]
   set indaux -1
   if {[llength $indexes2] > [llength $indexes1]} {
      set Molefacture::anglemvgrp1list $indexes2
      set Molefacture::anglemvgrp2list $indexes1
      set indaux $ind1
      set ind1 $ind2
      set ind2 $indaux
   }

   ### Removing the mother from the group of atoms
   set grplist [list $Molefacture::anglemvgrp1list $Molefacture::anglemvgrp2list]
   set incrm 0
   foreach grp $grplist {
      set index [lsearch $grp $Molefacture::anglemother]
      if {$incrm == 0 && $index != -1} {
         set Molefacture::anglemvgrp1list [lreplace $Molefacture::anglemvgrp1list $index $index]
      } elseif {$incrm == 1 && $index != -1} {
         set Molefacture::anglemvgrp2list [lreplace $Molefacture::anglemvgrp2list $index $index]
      }
      incr incrm
   }

   ### Store indexes and selections to avoid doing all
   ### these operations per resize_angle
   set Molefacture::angleind1 $ind1
   set Molefacture::angleind2 $ind2

   set mothersel [atomselect $tmpmolid "index $Molefacture::anglemother"]
   set ind1sel [atomselect $tmpmolid "index $Molefacture::angleind1"]
   set ind2sel [atomselect $tmpmolid "index $Molefacture::angleind2"]

   set Molefacture::anglemothercoor [lindex [$mothersel get {x y z}] 0]
   set Molefacture::angleind1coor [lindex [$ind1sel get {x y z}] 0]
   set Molefacture::angleind2coor [lindex [$ind2sel get {x y z}] 0]

   $mothersel delete
   $ind1sel delete
   $ind2sel delete

   ### Store the angle in memory and not compute every time
   ### because negative angles might be difficult to differentiate
   ### due to the symmetry of the molecule
   set Molefacture::curangle [format %.2f [measure angle [list $Molefacture::angleind1 $Molefacture::anglemother $Molefacture::angleind2] molid $tmpmolid]]

   set Molefacture::angleguiValscale $Molefacture::curangle
   set Molefacture::angleguiValentry $Molefacture::angleguiValscale
   return 1
}

#########################################################
# Evaluate if the selected atoms make an dihedral      ##
# angle (bonded) and if any of the groups is selected, ##
#draw the atoms                                        ##
#########################################################
proc ::Molefacture::isdihbond {} {
   variable tmpmolid
   variable picklist
   
   ### Check if they are bonded and select the two different groups and represent the
   ### the angle atoms
   set index1 [lindex $picklist 0]
   set index2 [lindex $picklist 1]
   set index3 [lindex $picklist 2]
   set index4 [lindex $picklist 3]

   # set idx1sel [atomselect $tmpmolid "index $index1"]
   # set idx2sel [atomselect $tmpmolid "index $index2"]
   # set idx3sel [atomselect $tmpmolid "index $index3"]
   # set idx4sel [atomselect $tmpmolid "index $index4"]

   set allsel [atomselect $tmpmolid "index $index1 $index2 $index3 $index4" ]

   set bondlist [$allsel getbonds]
   set dihlist [::util::dihedlist -bonds $bondlist -sel $allsel -all]
   set dihdef ""
   foreach dih $dihlist {
      if {[lsearch $dih $index1] != -1 && \
         [lsearch $dih $index2] != -1 && \
         [lsearch $dih $index3] != -1 && \
         [lsearch $dih $index4] != -1} {
            set dihdef $dih
           break
      }
   }
   $allsel delete
   if {$dihdef == ""} {
      return 0
   }

   set mother ""
   set indexlist "$index1 $index2 $index3 $index4"
   set countlist [list]
   foreach ind $indexlist {
      set count [llength [lsearch -all [join $bondlist] $ind ]]
      lappend countlist $count
   }

   set orderindex [lsort -integer -decreasing -indices $countlist]
   set orderlist [lsort -integer -decreasing $countlist]
   set Molefacture::dihemother1 [lindex $indexlist [lindex $orderindex 0]]
   set Molefacture::dihemother2 [lindex $indexlist [lindex $orderindex 1]]

   ### set ind1 the atom bonded to the mother1 and ind2 the atom bonded to mother2
   ### to define group1 and group2 to be moved

   set ind1 [lindex $indexlist [lindex $orderindex 2]]
   set ind2 [lindex $indexlist [lindex $orderindex 3]]
   set mother1bond [lindex $bondlist [lindex $orderindex 0]]
   set mother2bond [lindex $bondlist [lindex $orderindex 0]]
   set indaux -1
   if {[lsearch $mother1bond $ind1] == -1} {
      set indaux $ind1
      set ind1 $ind2
      set ind2 $indaux
   }

   ### When dividing the molecule in the mother atom,
   ### Set the group 2 the group of bonded atoms making the smallest bundle
   set molsize [$Molefacture::atomTable size]
   set indexes1 [join [::util::bondedsel $tmpmolid $Molefacture::dihemother1 $Molefacture::dihemother2 -maxdepth $molsize]]
   set indexes2 [join [::util::bondedsel $tmpmolid $Molefacture::dihemother2 $Molefacture::dihemother1 -maxdepth $molsize]]

   set Molefacture::anglemvgrp1list $indexes1
   set Molefacture::anglemvgrp2list $indexes2

   set indaux -1
   if {[llength $indexes2] > [llength $indexes1]} {
      set Molefacture::anglemvgrp1list $indexes2
      set Molefacture::anglemvgrp2list $indexes1
      set indaux $ind1
      set ind1 $ind2
      set ind2 $indaux
   }

   ### Removing the mother from the group of atoms
   set grplist [list $Molefacture::anglemvgrp1list $Molefacture::anglemvgrp2list]

   ### Removing the mothers from the group of atoms
   set index1 [lsearch $Molefacture::anglemvgrp1list $Molefacture::dihemother1]
   if {$index1 != -1} {
      set Molefacture::anglemvgrp1list [lreplace $Molefacture::anglemvgrp1list $index1 $index1]
   }
   set index2 [lsearch $Molefacture::anglemvgrp1list $Molefacture::dihemother2]
   set Molefacture::anglemvgrp1list [lreplace $Molefacture::anglemvgrp1list $index2 $index2]

   set index1 [lsearch $Molefacture::anglemvgrp2list $Molefacture::dihemother1]
   if {$index1 != -1} {
      set Molefacture::anglemvgrp2list [lreplace $Molefacture::anglemvgrp2list $index1 $index1]
   }
   set index2 [lsearch $Molefacture::anglemvgrp2list $Molefacture::dihemother2]
   set Molefacture::anglemvgrp2list [lreplace $Molefacture::anglemvgrp2list $index2 $index2]



   ### Store indexes and selections to avoid doing all
   ### these operations per resize_angle
   set Molefacture::angleind1 $ind1
   set Molefacture::angleind2 $ind2

   set mothersel1 [atomselect $tmpmolid "index $Molefacture::dihemother1"]
   set mothersel2 [atomselect $tmpmolid "index $Molefacture::dihemother2"]
   

   set Molefacture::dihemother1coor [lindex [$mothersel1 get {x y z}] 0]
   set Molefacture::dihemother2coor [lindex [$mothersel2 get {x y z}] 0]


   $mothersel1 delete
   $mothersel2 delete

   ### Store the angle in memory and not compute every time
   ### because negative angles might be difficult to differentiate
   ### due to the symmetry of the molecule
   set Molefacture::curangle [format %.2f [measure dihed [list $Molefacture::angleind1 $Molefacture::dihemother1 $Molefacture::dihemother2 $Molefacture::angleind2] molid $tmpmolid]]

   set Molefacture::dihguiValscale $Molefacture::curangle
   set Molefacture::dihguiValentry $Molefacture::dihguiValscale
   return 1
}

#########################################################
# Validate bond length entry and assign the same        #
# value to the bond scale (this value cannot be empty)  #
#########################################################
proc ::Molefacture::validateBondEntry {} {
   if {$Molefacture::bondguiValentry == "" || $Molefacture::bondguiValentry > 5.0} {
      set Molefacture::bondguiValentry 1.00
   }
   set Molefacture::bondguiValscale $Molefacture::bondguiValentry
   return 1
}


#########################################################
# Validate angle entry and assign the same              #
# value to the angle scale (this value cannot be empty) #
#########################################################
proc ::Molefacture::validateAngleEntry {} {
   variable picklist
   set val ""
   if {[llength $picklist] == 3} {
      set val $Molefacture::angleguiValentry
   } elseif {[llength $picklist] == 4} {
      set val $Molefacture::dihguiValentry
   }

   if {$val == "" || $val <= 0.0} {
      set val 0.01
   } elseif {$val > 360 } {
      set val 360.00
   } elseif {$val == 180.00} {
      set val 180.01
   }

   if {[llength $picklist] == 3} {
      set Molefacture::angleguiValentry $val
      set Molefacture::angleguiValscale $Molefacture::angleguiValentry
      # set Molefacture::curangle 0.1
   } elseif {[llength $picklist] == 4} {
      set Molefacture::dihguiValentry $val
      set Molefacture::dihguiValscale $Molefacture::dihguiValentry
   }


   return 1
}

#########################################################
# Validate angle entry and assign the same              #
# value to the angle scale (this value cannot be empty) #
#########################################################
proc ::Molefacture::validateAngleScale { angle } {
   variable picklist
   set val ""
   if {[llength $picklist] == 3} {
      set val $Molefacture::angleguiValscale
   } else {
       set val $Molefacture::dihguiValscale
   }

   if {$val == 180.00} {
      set val 180.01
   } elseif {$val > 360.00} {
      set val 360.00
   }

   if {[llength $picklist] == 3} {
      set Molefacture::angleguiValscale $val
      set Molefacture::angleguiValentry $val
      ::Molefacture::resize_angle $angle
   } elseif {[llength $picklist] == 4} {
      set Molefacture::dihguiValscale $val
      set Molefacture::dihguiValentry $val
      ::Molefacture::rotate_bond $angle
   }

}

#########################################################
# Validate the Charge Value                             #
#########################################################
proc Molefacture::validateCharge {} {
   if {$Molefacture::charge == "" || [string is double [string trim $Molefacture::charge " "]] != 1} {
      tk_messageBox -message "Please choose a valid value for Charge" -icon error -type ok -parent $Molefacture::topGui -title "Charge Error"
      set Molefacture::charge "0.000"

      return 0
   } else {
      set Molefacture::charge [Molefacture::format3Dec $Molefacture::charge]
      return 1
   }

}


#########################################################
# Validate the Residue Name                             #
#########################################################
proc Molefacture::validateResname {} {
   set Molefacture::resname [string trim $Molefacture::resname]
   if {$Molefacture::resname == "" || [string length $Molefacture::resname] > 4} {
      tk_messageBox -message "Please choose a valid value for Resname. Residue's name cannot be more than 4 characters long." -icon error -type ok -parent $Molefacture::topGui -title "Resname Error"
      set Molefacture::resname "TMP"
      return  0
   }

   return 1
}

#########################################################
# Validate the Chain ID                                 #
#########################################################
proc Molefacture::validateChain {} {
   set Molefacture::chain [string trim $Molefacture::chain]
   if {$Molefacture::chain == "" || [string length $Molefacture::chain] > 1} {
      tk_messageBox -message "Please choose a valid value for Chain ID. Chain's ID cannot be more than 1 character long." -icon error -type ok -parent $Molefacture::topGui -title "Chain ID Error"
      set Molefacture::chain "A"
      return 0
   }
   return 1
}

#########################################################
# Validate the Segment Name                             #
#########################################################
proc Molefacture::validateSegname {} {
   set Molefacture::segname [string trim $Molefacture::segname]
   if {[string length $Molefacture::segname] > 4 || [string length $Molefacture::segname] == 0 || $Molefacture::segname == "{}"} {
      tk_messageBox -message "Please choose a valid value for Segment Name. Segment's Name cannot be more than 4 characters long or empty." -icon error -type ok -parent $Molefacture::topGui -title "Segname Error"
      set Molefacture::segname "A"
      return 0
   }
   return 1
}

#########################################################
# Validate the Residue ID                               #
#########################################################
proc Molefacture::validateResid {} {
   set Molefacture::resid [string trim $Molefacture::resid]
   if {$Molefacture::resid == "---"} {
      return 1
   }
   if {$Molefacture::resid == "" || [string is integer $Molefacture::resid] != 1} {
      tk_messageBox -message "Please choose a valid value for Residue ID. Residue's ID can only accept integer (positive and negative) numbers." -icon error -type ok -parent $Molefacture::topGui -title "Chain ID Error"
      set Molefacture::resid "1"
      return 0
   }
   return 1
}
#######################################################
# Start of the edit process of the Atom's Table       #
#######################################################
proc Molefacture::startEditTable {tbl row col txt} {
   set txt [string trim $txt]
   set tblsel [$tbl curselection]
   set indexlist ""
   set w [$tbl editwinpath]

   switch [$tbl columncget $col -name] {
      OxState {
         $w configure -state normal -style TblEdit.TCombobox -values $Molefacture::cmboxstatvalues
         bind $w <<ComboboxSelected>> {
            if {[winfo exists %W]} {
               $Molefacture::atomTable finishediting
            }
         }

      }
   }
   return $txt
}
#######################################################
# Validation of the edit process of the Atom's Table  #
#######################################################
proc Molefacture::validateEditTable {tbl row col txt} {
   # set prevtext [$tbl cellcget $row,$col -text]
   set txt [string trim $txt]
   set tblsel [$tbl curselection]
   set indexlist ""
   set tableindexlist ""
   foreach row $tblsel {
      lappend indexlist [$tbl cellcget $row,Index -text]
   }
   if {[llength $indexlist] == 0} {
      return
   }
   if {[$tbl canceledediting] == 1 } {
      return
   }
   ### Column 1 - Atom's name
   switch [$tbl columncget $col -name] {

      AtmNAME {
         set atmnameslist [lreplace [$tbl getcolumns AtmNAME] $row $row]
         set equal [lsearch $atmnameslist $txt]
         if {$txt == "" || [string length $txt] > 4 || $equal > -1} {
            tk_messageBox -message "Please choose a unique valid value for Atom's Name. Atom's Name cannot be more than 4 characters long." -icon error -type ok -parent $Molefacture::topGui -title "Atom's Name Error"

            $tbl cancelediting
         }  else {
            Molefacture::saveState
            set sel [atomselect $Molefacture::tmpmolid "index [$tbl cellcget $row,Index -text]"]
            $sel set name $txt
            $sel delete
         }

      }
      AtmType {
         if {$txt == "" || [string length $txt] > 6} {
            tk_messageBox -message "Please choose a unique valid value for Atom's Type. Atom's Name should not be more than 6 characters long." -icon error -type ok -parent $Molefacture::topGui -title "Atom's Type Error"
            $tbl rejectinput
            $tbl cancelediting
         } else {
            Molefacture::saveState
            set sel [atomselect $Molefacture::tmpmolid "index [join $indexlist]"]
            foreach ind $tblsel {
               $tbl cellconfigure $ind,$col -text "$txt"
            }
            $sel set type $txt
            $sel delete
         }
      }
      OxState {
         if {$txt != ""} {
            Molefacture::saveState
            set Molefacture::cmboxstat $txt
            Molefacture::updateOxiState
         }
      }
      Charge {
         set Molefacture::charge $txt
         if {![Molefacture::validateCharge]} {
            set Molefacture::charge [$tbl cellcget $row,$col -text]
            $tbl rejectinput
            $tbl cancelediting
         } else {
            Molefacture::saveState
            foreach ind $tblsel {
               $tbl cellconfigure $ind,$col -text [Molefacture::format3Dec "$txt"]
            }
            set sel [atomselect $Molefacture::tmpmolid "index [join $indexlist]"]
            $sel set charge [Molefacture::format3Dec "$txt"]
            $sel delete
         }

      }
      Resname {
         set Molefacture::resname $txt
         if {![Molefacture::validateResname]} {
            $tbl rejectinput
            $tbl cancelediting
         } else {
            set Molefacture::resname $txt
            Molefacture::saveState
            foreach ind $tblsel {
               $tbl cellconfigure $ind,$col -text "$txt"
            }
            set sel [atomselect $Molefacture::tmpmolid "index [join $indexlist]"]
            $sel set resname $txt
            $sel delete
            Molefacture::loadDefaulTopo
         }

      }
      ResID {
         set Molefacture::resid $txt

         if {![Molefacture::validateResid]} {
            set Molefacture::resid [$tbl cellcget $row,$col -text]
            $tbl rejectinput
            $tbl cancelediting
         } else {
            Molefacture::saveState
            foreach ind $tblsel {
               $tbl cellconfigure $ind,$col -text "$txt"
            }
            set sel [atomselect $Molefacture::tmpmolid "index [join $indexlist]"]
            $sel set resid $txt
            $sel delete
         }

      }
      ChainID {
         set Molefacture::chain $txt
         if {![Molefacture::validateChain]} {
            set Molefacture::chain [$tbl cellcget $row,$col -text]
            $tbl rejectinput
            $tbl cancelediting
         } else {
            Molefacture::saveState
            foreach ind $tblsel {
               $tbl cellconfigure $ind,$col -text "$txt"
            }
            set sel [atomselect $Molefacture::tmpmolid "index [join $indexlist]"]
            $sel set chain $txt
            $sel delete
         }

      }
      Segname {
         set Molefacture::segname $txt
         if {![Molefacture::validateSegname]} {
            set Molefacture::segname [$tbl cellcget $row,$col -text]
            $tbl rejectinput
            $tbl cancelediting
         } else {
            Molefacture::saveState
            foreach ind $tblsel {
               $tbl cellconfigure $ind,$col -text "$txt"
            }
            set sel [atomselect $Molefacture::tmpmolid "index [join $indexlist]"]
            $sel set segname $txt
            $sel delete
         }

      }
   }
   mol reanalyze $Molefacture::tmpmolid
   if {$Molefacture::topocgenff == 1} {
      set Molefacture::topocgenff 0
   }

   return $txt
}

#########################################################
# Validate the number of steps to minimize structure    #
#########################################################
proc ::Molefacture::validateMinimizeSteps {  } {
   set Molefacture::minsteps [expr abs(int($Molefacture::minsteps))]
   return 1
}
#####################################################
# Reselect table row(s) after other Gui calls       #
#####################################################
proc ::Molefacture::keepTableSelection {} {
   set index $Molefacture::lastAtmSel
   $Molefacture::atomTable selection set $index
   ::Molefacture::updateEditMolecule
   return $index
}

#####################################################
# Select all atoms                                  #
#####################################################
proc ::Molefacture::selectAllTable {} {
  set ::Molefacture::mousemode pick
   if {[$Molefacture::atomTable size] > 0} {
      $Molefacture::atomTable selectio set 0 end
      ::Molefacture::selectAtomTable
   }

}

#####################################################
# Clear Table selection                             #
#####################################################
proc ::Molefacture::clearTableSelection {} {
   global vmd_pick_atom
   $Molefacture::atomTable selectio clear 0 end
   set Molefacture::picklist {}
   set vmd_pick_atom {}
   Molefacture::draw_selatoms
   set Molefacture::lastAtmSel ""
   ### Reset and disable the bond modifiers

   set Molefacture::numAtoms [$Molefacture::atomTable size]

   if {$Molefacture::tmpmolid == -1 || [lsearch [molinfo list] $Molefacture::tmpmolid] == -1} {
    return
   }

   Molefacture::disablebond
   Molefacture::disableangle
   Molefacture::disabledih
   Molefacture::updateEditMolecule
   Molefacture::oxiStateCombo
}

#########################################################
# Disable and reset variables related to the bond gui  #
#########################################################
proc ::Molefacture::disablebond {} {
   ### Reset and disable the bond modifiers
   set Molefacture::bondguiValscale 0
   set Molefacture::bondguiValentry 0
   $Molefacture::bondentrywdgt configure -state disabled
   $Molefacture::bondscalewdgt configure -state disabled
}

#########################################################
# Disable and reset variables related to the angle gui  #
#########################################################
proc ::Molefacture::disableangle {} {
   ### Reset and disable the angle modifiers
   set Molefacture::angleguiValscale 0
   set Molefacture::angleguiValentry 0
   $Molefacture::angleentrywdgt configure -state disabled
   $Molefacture::anglescalewdgt configure -state disabled
   $Molefacture::anglemvgrp1wdgt configure -state disabled
   $Molefacture::anglemvgrp2wdgt configure -state disabled
   if {$Molefacture::anglemothercoor != ""} {
      set Molefacture::anglemothercoor ""
      set Molefacture::angleind1coor ""
      set Molefacture::angleind2coor ""
   }
   set Molefacture::anglemvgrp1 0
   set Molefacture::anglemvgrp2 0
   set Molefacture::anglemvgrp1list -1
   set Molefacture::anglemvgrp2list -1
   set Molefacture::curangle -1
   Molefacture::draw_angleatoms
}

#########################################################
# Disable and reset variables related to the angle gui  #
#########################################################
proc ::Molefacture::disabledih {} {
   ### Reset and disable the angle modifiers
   set Molefacture::dihguiValscale 0
   set Molefacture::dihguiValentry 0
   $Molefacture::dihscalewdgt configure -state disabled
   $Molefacture::dihentrywdgt configure -state disabled
   $Molefacture::dihradbtt1 configure -state disabled
   $Molefacture::dihradbtt2 configure -state disabled

   set Molefacture::dihemother1coor ""
   set Molefacture::dihemother2coor ""

   set Molefacture::dihedrgrp ""
   set Molefacture::anglemvgrp1list -1
   set Molefacture::anglemvgrp2list -1
   set Molefacture::curangle -1
   Molefacture::draw_angleatoms
}

proc ::Molefacture::hideFrame {w frame txt} {
    set frameaux "$frame.fcolapse"
    set arrow [lindex [$w cget -image] 0]
    set pad 0
    if {$arrow != $Molefacture::arrowRight} {
        $w configure -text "$txt" -compound left -image $Molefacture::arrowRight
        grid forget $frameaux
        set info [grid info [lindex [grid info $w] 1] ]
        grid rowconfigure [lindex $info 1] [lindex $info 5] -weight 0
    } else {
        $w configure -text "$txt" -compound left -image $Molefacture::arrowDown
        set info [grid info $w]
        grid conf $frameaux -row [expr [lindex $info 5] +1] -column [lindex $info 3] -pady 1 -padx 2 -sticky ewns
        set info [grid info [lindex [grid info $w] 1] ]
        grid rowconfigure [lindex $info 1] [lindex $info 5] -weight 1

        set pad 150
    }
   regexp {([0-9]+)x[0-9]+[\+\-]+[0-9]+[\+\-]+[0-9]+} [wm geometry $Molefacture::topGui] all dimW
   set dimH [winfo reqheight $Molefacture::topGui]
   set dimH [expr {$dimH + $pad}]

   wm geometry $Molefacture::topGui [format "%ix%i" $dimW $dimH]
}
###################################################
# Periodic table window                           #
###################################################

proc ::Molefacture::periodicTableWindow {} {
   variable periTabWindow
   variable periodic
   variable periTab
   if {[winfo exists $periTabWindow] != 1} {
     toplevel $periTabWindow
   } else {
      wm deiconify $periTabWindow
      focus $periTabWindow
      raise $periTabWindow
      return
   }

   ## periodic table canvas option
   set canh   290;   # canvas height
   set canw   560;  # canvas width
   set linecolor black;    # color of lines connecting data points
   set bkgcolor white;    # color of the canvas background

   grid columnconfigure $periTabWindow 0 -weight 1
   grid rowconfigure $periTabWindow 0 -weight 1

   ## Title of the windows

   wm title $periTabWindow "Molefacture - Element Selection"
   wm resizable $periTabWindow 0 0
   ## main frame
   grid [ttk::frame $periTabWindow.topframe] -row 0 -column 0 -sticky nwes -padx 2 -pady 2
   grid columnconfigure $periTabWindow.topframe 0 -weight 1
   grid rowconfigure $periTabWindow.topframe 0 -weight 1

   grid [canvas $periTabWindow.topframe.cv -relief flat -borderwidth 0 -width $canw -height $canh -bg $bkgcolor] -row 0 -column 0 -sticky nsew -padx 2 -pady 2

   set periTab $periTabWindow.topframe.cv
   ### Fill the canvas


   set atmnum 1

   for {set atmnum 1} {$atmnum <= [llength $periodic]} { incr atmnum} {
      ::Molefacture::addElemet $atmnum
   }
}



proc ::Molefacture::select_atoms_byvmdindex {indexlist} {
  #Helper proc that translates vmd indices into indices from .molefac.val.list.list
  variable atomlist
  set outputlist [list]
  set translist [list]
  foreach molefind $atomlist {
    lappend translist [lindex $molefind 0]
  }
#  puts "Translation indices $translist"

  foreach vmdind $indexlist {
#    puts "Looking for $vmdind in $translist"
    set i [lsearch -exact -integer $translist $vmdind]
    if {$i != -1} {lappend outputlist $i}
#    puts "Found $i"
  }

#  puts "found indices $outputlist"
  select_atoms $outputlist
}

#######################################
# proc to get representation index    #
#######################################
proc ::Molefacture::getrepnum {repname mol} {
    return [mol repindex $mol $repname]
}

####################################################################
# This function is used to select atoms.                           #
# It paints the background of the selected list items accordingly  #
# and appends the atomindexes to variable picklist in order to be  #
# independent of the window focus.                                 #
####################################################################

proc ::Molefacture::select_atoms { indexlist } {

   variable picklist {}

   draw_selatoms
   foreach index $indexlist {
      lappend picklist $indexatomind
   }
}

proc ::Molefacture::add_atoms_to_selection { indexlist } {
   variable picklist
   variable atomlist
   variable selectcolor
   if {![llength $indexlist] || ![winfo exists .molefac.val.list.list]} { return }

#   puts "DEBUG: Indexlist: $indexlist"
   foreach index $indexlist {
#      puts "DEBUG: found $i"
      .molefac.val.list.list selection set $index
      .molefac.val.list.list itemconfigure $index -background $selectcolor
      set indexatomind [lindex [.molefac.val.list.list get $index] 0]
      lappend picklist $indexatomind
   }

   draw_selatoms
   .molefac.val.list.list see $index
   .molefac.val.list.list activate $index
}


####################################################################
# This function can be used by external programs to select atoms.  #
####################################################################

proc ::Molefacture::user_select_atoms { atomdeflist } {
   variable tmpmolid

   # Select the corresponding atoms
   set indexlist {}
   foreach atomdef $atomdeflist {
      set sel [atomselect $tmpmolid "segid [lindex $atomdef 0] and resid [lindex $atomdef 1] and name [lindex $atomdef 2] "]
      lappend indexlist [join [$sel get index]]
      $sel delete
   }

   select_atoms $indexlist
}

proc ::Molefacture::raise_bondorder_gui { } {
   variable tmpmolid

   variable picklist; #[.molefac.val.list.list curselection]
   if {[llength $picklist]!=2} {
      tk_messageBox -icon error -type ok -parent $Molefacture::topGui -title Message \
	 -message "To modify the bond order, you should select exactly two atoms!"
      return
   }

   set Molefacture::topocgenff 0

   raise_bondorder $picklist
}

proc ::Molefacture::lower_bondorder_gui { } {
   variable tmpmolid

   variable picklist; # [.molefac.val.list.list curselection]
   if {[llength $picklist]!=2} {
      tk_messageBox -icon error -type ok -title Message \
	 -message "To modify the bond order, you should select exactly two atoms!"
      return
   }
   lower_bondorder $picklist
   set Molefacture::topocgenff 0
}

proc ::Molefacture::invert_gui {} {
  variable tmpmolid

  variable picklist
  foreach mindex $picklist {
    invert_chir $mindex
  }
  draw_openvalence
}

proc ::Molefacture::set_planar_gui { } {
  variable tmpmolid

  variable picklist;# [.molefac.val.list.list curselection]
  foreach mindex $picklist {
    set_planar $mindex
  }
  draw_openvalence
}


proc ::Molefacture::set_tetra_gui { } {
  variable tmpmolid

  variable picklist;
  foreach mindex $picklist {
    set_tetra $mindex
  }
}

proc ::Molefacture::del_atom_gui {} {
  variable picklist
  variable tmpmolid
  # variable openvalencelist
  variable picklist
  # set mother 0
  # set curcoor {0 0 0}

  if {[llength $picklist] == 0} {
   return
  }
   Molefacture::saveState
   set retlist [del_atom $picklist 0]
   set tbindex [list]
   set tbcolumns [$Molefacture::atomTable getcolumns Index]
   foreach delindex $picklist {
    lappend tbindex [lsearch $tbcolumns $delindex]
   }
   $Molefacture::atomTable delete $tbindex

   set Molefacture::numAtoms [$Molefacture::atomTable size]
  # update_bondlist
  update_openvalence
  variable bondlist [bondlist]
  variable anglelist [anglelist]

   if {[$Molefacture::atomTable size] == 0} {
      set Molefacture::resname ""
      set Molefacture::resid ""
      set Molefacture::chain ""
      set Molefacture::charge 0.000
      set Molefacture::segname ""
      $Molefacture::rnamewdgt configure -state disabled
      $Molefacture::chnwdgt configure -state disabled
      $Molefacture::rsidwdgt configure -state disabled
      $Molefacture::sgnmwdgt configure -state disabled
   }
   mol reanalyze $tmpmolid
    if {$Molefacture::topocgenff == 1} {
      set Molefacture::topocgenff 0
   }
   
}

#######################
## Add new atom Atom ##
#######################
proc ::Molefacture::add_atom_gui { element {updateTable 1}} {
  variable tmpmolid
  variable picklist
  variable topGui
  Molefacture::saveState

  if {$tmpmolid == -1} {
   Molefacture::molefacture_gui_aux
  }
   set newatoms {}
   ### select mode of adding the atoms
   ### do == 0; don't add any atom
   ### do == 1; add atom
   set do 0
  if {[llength $picklist] == 0} {
    if {[$Molefacture::atomTable size] > 0} {
      if {[$Molefacture::atomTable size] == 1} {
        set do 1
      } else {
        return
      }
    } else {
      set do 1
    }
  } else {
     set do 1
  }
  if {$do == 0} {
    tk_messageBox -message "Please select an atom to attach the new atom" -icon info -type ok -parent $topGui
    return ""
  }
  if {$do == 1} {
   set index [add_atom $element]
   if {$index == -1} {
      tk_messageBox -message "No more atoms are allowed to this atom." -icon warning -type ok -parent $topGui
      return
   }
   lappend newatoms $index
  }

   set newatoms [join $newatoms]

   if {[llength $newatoms] != 0} {
      topo -molid "$tmpmolid" -sel "index $newatoms" guessatom radius element
      topo -molid "$tmpmolid" -sel "index $newatoms" guessatom mass element
   } else {
      return
   }

   mol reanalyze $tmpmolid
   if {$updateTable == 1} {
      ::Molefacture::updateAtomTable
   }

   set Molefacture::topocgenff 0
   return $newatoms
}


proc ::Molefacture::export_molecule_gui {} {

   variable tmpmolid
   set types {
           {{XBGF Files} {.xbgf} }
           {{PDB Files} {.pdb} }
           {{MOL2 Files} {.mol2} }
           }
   set filename [tk_getSaveFile -parent .molefac -filetypes $types -defaultextension ".pdb"]

   if {[llength $filename] == 0} {
      return
   }

   fix_changes
   set sel [atomselect $tmpmolid "occupancy>=0.8"]
   if {[regexp {xbgf$} $filename] > 0 } {
     write_xbgf "$filename" $sel
   } elseif {[regexp {mol2$} $filename] > 0} {
     $sel writemol2 "$filename"
   } else {
     $sel writepdb "$filename"
   }
   $sel delete

   variable projectsaved
   set projectsaved 1
}

proc ::Molefacture::export_molecule_gui_FEP {} {
   fix_changes
   variable tmpmolid
   variable filename_FEP
   set filename_FEP [tk_getSaveFile -parent .molefac ]
   write_topfile $filename_FEP.top "occupancy >= 0.8"
   set sel [atomselect $tmpmolid "occupancy>=0.8"]

   $sel writepdb "$filename_FEP.pdb"

   $sel writepdb "$filename_FEP.fep"
   variable projectsaved
   set projectsaved 1
}

proc ::Molefacture::apply_moledits { } {
   variable picklist
   variable tmpmolid
   
   if {![Molefacture::validateResid] || ![Molefacture::validateResname] || ![Molefacture::validateSegname] || ![Molefacture::validateChain]} {
      return 0
   }

   if {$Molefacture::resname == "---" && $Molefacture::segname == "---" && $Molefacture::chain == "---" && $Molefacture::resid == "---"} {
      return 0
   }

   #### true if molecule was edited
   set edited 0
   Molefacture::saveState
   set seltext "$Molefacture::notDumSel"
   set indexlist [$Molefacture::atomTable curselection]
  
   if {[llength $picklist] != 0} {
      set seltext "index [join $picklist]"
   } else {
      set indexlist ""
      for {set i 0} {$i < [$Molefacture::atomTable size]} { incr i } {
         lappend indexlist $i
      }
      set seltex "index [join $indexlist]"

   }
   
   set tmpsel [atomselect $tmpmolid "$seltext"]
   
   if {$Molefacture::resname != "---"} {
     set Molefacture::topocgenff 0
      $tmpsel set resname "$::Molefacture::resname"
      foreach ind $indexlist {
         $Molefacture::atomTable cellconfigure $ind,Resname -text $::Molefacture::resname
      }

      set edited 1
      foreach topo $Molefacture::topoinfo {
         if {[::Toporead::topology_contains_resi $topo $Molefacture::resname] > -1} {
            set edited 0
            break
         }
      }
      set curres [lsort -unique [$Molefacture::atomTable getcolumns Resname]]
      set Molefacture::unknownRes [lindex [Molefacture::loadDefaulTopo $curres] 0]
   }
   if {$Molefacture::segname != "---"} {
      $tmpsel set segname "$::Molefacture::segname"
      foreach ind $indexlist {
         $Molefacture::atomTable cellconfigure $ind,Segname -text $::Molefacture::segname
      }
      set Molefacture::topocgenff 0
      set edited 2
   }
   if {$Molefacture::chain != "---"} {
      $tmpsel set chain "$::Molefacture::chain"
      foreach ind $indexlist {
         $Molefacture::atomTable cellconfigure $ind,ChainID -text $::Molefacture::chain
      }
      set Molefacture::topocgenff 0
      set edited 2
   }
   if {$Molefacture::resid != "---"} {
      $tmpsel set resid "$::Molefacture::resid"
      foreach ind $indexlist {
         $Molefacture::atomTable cellconfigure $ind,ResID -text $::Molefacture::resid
      }
      set Molefacture::topocgenff 0
      set edited 2
   }

   #### Update Charges
   set size [llength $picklist]
   if {$size == 0} {
      set size [$Molefacture::atomTable size]
   }
   if {$Molefacture::topocgenff == 1 && $Molefacture::resname != "---"} {
      set topocharge [::Toporead::topology_get_resid [lindex $Molefacture::topoinfo 0] $Molefacture::resname charge]
      if {$Molefacture::charge != "---" && $topocharge != $Molefacture::charge} {
         set Molefacture::topocgenff 0
      }
   }

   ### Block this functionality until evaluate if it is needed. This functions 
   ### devided the total charge by all atoms. Might not be a good idea
   if {$Molefacture::charge != "---" && $Molefacture::topocgenff == 0 && 0} {
     set charge [Molefacture::format3Dec [expr $Molefacture::charge / [expr $size * 1.0]]]

     if {$size > 0} {

         set chrglist [list]
         set totaux ""

         if {$size > 1} {
            set totaux [eval "vecadd [$tmpsel get charge]"]
         } else {
            set totaux [$tmpsel get charge]
         }
         set diff [Molefacture::format3Dec [expr $totaux - $Molefacture::charge]]

         set row 0
         set sign "+"
         if {$diff < 0} {
             set sign "-"
         }
         set rstart 0
         # set indexcolumn [$Molefacture::atomTable getcolumns Index]
         for {set i 0} {$i < $size} {incr i} {
            # set ind $i
            # if {$alltbl == 0 } {
            #    set ind [lsearch $indexcolumn [lindex $indexlist $row]]
            # }
            set ind [lindex $indexlist $i]
            $Molefacture::atomTable cellconfigure $ind,Charge -text 0.000
            lappend chrglist "0.000"
         }


         if {[Molefacture::format3Dec $Molefacture::charge] != 0.000} {
            set chrglist [list]
            while {$diff != 0.000} {
                set ind [lindex $indexlist $row]

                set current [$Molefacture::atomTable cellcget $ind,Charge -text]
                set chrg [Molefacture::format3Dec [expr $current - ${sign}0.001]]
                $Molefacture::atomTable cellconfigure $ind,Charge -text $chrg
                set diff [Molefacture::format3Dec [expr $diff - ${sign}0.001]]
                # puts "DEBUG row $row || size $size ||rstart $rstart ||  chrglist size [llength $chrglist]"
                if {$rstart == 0} {
                  lappend chrglist $chrg
                } else {
                  lset chrglist $row $chrg
                }
                incr row
                if {$row == $size } {
                    set row 0
                    set rstart 1
                }
            }

         }

         $tmpsel set charge $chrglist
     }

      set edited 1

  }

  if {$edited == 1} {
      set Molefacture::topocgenff 0
  }
  mol reanalyze $tmpmolid
  $tmpsel delete
}


proc ::Molefacture::check_at_gui { } {
  variable topGui
  if { [winfo exists .checkatgui] } {
    wm deiconify .checkatgui
    raise .checkatgui
    return
  }


  set w [toplevel ".checkatgui"]
  wm title $topGui "Check Atom Typing"
  wm resizable $topGui 0 0

  set rownum 0
  frame $topGui.directories
  set refdir ""
  grid [label $topGui.directories.refdirlabel -text "Directory containing reference structures:"] -row $rownum -column 0 -sticky e
  grid [entry $topGui.directories.refdir -width 20 -textvar ::Molefacture::refdir] \
    -row $rownum -column 1
  grid [button $topGui.directories.refbrowse -text "Browse" -command ::Molefacture::find_refdir] -row $rownum -column 2 -columnspan 1 -sticky w
  incr rownum
  set testdir ""
  grid [label $topGui.directories.testdirlabel -text "Directory to run tests:"] -row $rownum -column 0 -sticky e
  grid [entry $topGui.directories.testdir -width 20 -textvar ::Molefacture::testdir] \
    -row $rownum -column 1
  grid [button $topGui.directories.testbrowse -text "Browse" -command ::Molefacture::find_testdir] -row $rownum -column 2 -columnspan 1 -sticky w
  incr rownum
  grid [label $topGui.directories.at_test -text "Atom type to check:"] -row $rownum -column 0 -sticky e
  grid [entry $topGui.directories.at_var -width 4 -textvar ::Molefacture::at_test] \
    -row $rownum -column 1 -sticky w
  incr rownum
  grid [button $topGui.directories.done -text "Run Comparisons" -command ::Molefacture::run_at_comps] -row $rownum -column 2 -columnspan 1 -sticky e
  pack $topGui.directories

}

proc ::Molefacture::find_refdir {} {
  set filename [tk_chooseDirectory -title "Choose a directory" -initialdir $::Molefacture::refdir]
  set ::Molefacture::refdir "$filename"
}

proc ::Molefacture::find_testdir {} {
  set filename [tk_chooseDirectory -title "Choose a directory" -initialdir $::Molefacture::testdir]
  set ::Molefacture::testdir "$filename"
}


#Procs to raise and lower oxidation state
proc ::Molefacture::raise_ox_gui {} {
  variable tmpmolid

  variable picklist;# [.molefac.val.list.list curselection]
  foreach mindex $picklist {
    raise_ox $mindex
  }
  update_openvalence
}

proc ::Molefacture::lower_ox_gui {} {
  variable tmpmolid

  variable picklist;# [.molefac.val.list.list curselection]
  foreach mindex $picklist {
    lower_ox $mindex
  }
  update_openvalence
}

proc ::Molefacture::fill_fragment_menu {} {
  # Looks through the current fragment database, and fills out the fragment menu
  # Each entry in the menu runs the replace hydrogen with fragment proc, using
  # the appropriate fragment
  # Currently clobbers all entries in the old menu, if any

  variable topGui
  variable addfrags

  foreach fragname [lsort [array names addfrags]] {
    set fragfile $addfrags($fragname)
    $topGui.topframe.menubar.build.menu.addfrag add command -label $fragname -command "::Molefacture::replace_atom_with_fragment_gui {$fragfile}"
  }

}


proc ::Molefacture::replace_atom_with_fragment_gui {fragpath} {
  # GUI dummy proc to find the atoms for replacement, and then replace them
  # with the appropriate fragment
  variable tmpmolid
  variable picklist
  foreach mindex $picklist {
    set returncode [replace_atom_with_fragment $fragpath $mindex]
    if {$returncode == -1} {
      tk_messageBox -type ok -title "Error" -icon error -message "You can only replace singly bonded hydrogen atoms! The atom you picked doesn't meet one or both of these criteria" -parent $Molefacture::topGui
    }
   }
   ::Molefacture::updateNewAtoms [join $returncode]
}

################################################################
## Update table and atoms' elements
################################################################
proc ::Molefacture::updateNewAtoms {newAtomsList} {

     # Select a newly added atom (usually a hydrogen)
    global vmd_pick_atom
    variable tmpmolid


   if {[llength $newAtomsList] != 0} {
      ##### Guess the elements raidus and mass for proper display
      topo -molid "$tmpmolid" -sel "index $newAtomsList" guessatom radius element
      topo -molid "$tmpmolid" -sel "index $newAtomsList" guessatom mass element
   } else {
      return
   }
    mol reanalyze $tmpmolid
   ::Molefacture::updateAtomTable

    set sel [atomselect $tmpmolid "occupancy >= 0.5"]
    set newind [lindex [$sel get index] end]
    set vmd_pick_atom $newind
    atom_picked_fctn
    $sel delete


}

proc ::Molefacture::add_custom_frags {} {
#First let them navigate to the file of interest
  variable topGui
  set fragfile [tk_getOpenFile]
  if {$fragfile == ""} {return}

  add_frag_file $fragfile
  $topGui.menubar.build.addfrag delete 3 end
  fill_fragment_menu
}

proc ::Molefacture::add_custom_basefrags {} {
#First let them navigate to the file of interest
  variable topGui
  set fragfile [tk_getOpenFile]
  if {$fragfile == ""} {return}

  add_basefrag_file $fragfile
  $topGui.menubar.build.basefrag delete 3 end
  fill_basefrag_menu
}

proc ::Molefacture::reset_basefrags_menu {} {
  variable topGui

  $topGui.menubar.build.basefrag delete 3 end
  read_basefrag_file [file join
   lib basemol basefrag.mdb]
  fill_basefrag_menu
}

proc ::Molefacture::fill_basefrag_menu {} {
  # Looks through the current base fragment database, and fills out the fragment menu
  # Each entry in the menu creates a new molecule from
  # the appropriate fragment
  # Currently clobbers all entries in the old menu, if any

  variable topGui
  variable basefrags

  foreach fragname [lsort [array names basefrags]] {
    set fragfile $basefrags($fragname)
    $topGui.topframe.menubar.build.menu.basefrag add command -label $fragname -command "::Molefacture::new_independent_fragment_gui {$fragfile}"
  }

}

proc ::Molefacture::new_independent_fragment_gui_gui {fragpath} {
  # GUI dummy proc to find the atoms for replacement, and then replace them
  # with the appropriate fragment
  variable tmpmolid
  variable picklist

    set returncode [new_independent_fragment $fragpath]
    if {$returncode == -1} {
      tk_messageBox -type ok -title "Error" -icon error -message "Error when adding an independent fragment" -parent $Molefacture::topGui
    }
   set newatoms [join $returncode]

   if {[llength $newatoms] != 0} {
      topo -molid "$tmpmolid" -sel "index $newatoms" guessatom radius element
      topo -molid "$tmpmolid" -sel "index $newatoms" guessatom mass element
   } else {
      return
   }

    mol reanalyze $tmpmolid
    display update
    display resetview
   ::Molefacture::updateAtomTable

}

proc ::Molefacture::new_mol_gui {} {
  # Use whenever you want to start working on a blank molecule
  # Prompts the user to make sure they don't need to save anything, and then
  # opens up a blank molecule with some preallocated hydrogens
  # Return 0 if they say no at the prompt, 1 otherwise

  set answer [tk_messageBox -message "This will abandon all editing on the current molecule. Are you sure?" -type yesno -icon question -parent $Molefacture::topGui]
  switch -- $answer {
    no {return 0}
    yes {}
  }

  new_mol
  return 1
}

proc ::Molefacture::nuc_builder_gui {} {

  variable ralpha
  variable rbeta
  variable rgamma
  variable rdelta
  variable repsilon
  variable rchi
  variable rzeta
  variable nucseq
  variable nucpath

  if {[winfo exists .nucbuilder]} {
    wm deiconify .nucbuilder
    raise .nucbuilder
    return
  }

  set w [toplevel ".nucbuilder"]
  wm title $w "Molefacture Nucleic Acid Builder"
  wm resizable $w no no

  #Frame for buttons for individual nucleotides
  labelframe $topGui.nucs -bd 2 -relief ridge -text "Add nucleotide" -padx 1m -pady 1m
  frame $topGui.nucs.buttons
  grid [button $topGui.nucs.buttons.ade -text "A" -command {::Molefacture::add_nuc A ; ::Molefacture::update_openvalence}] -row 0 -column 0 -columnspan 1 -sticky nsew
  grid [button $topGui.nucs.buttons.cyt -text "C" -command {::Molefacture::add_nuc C ; ::Molefacture::update_openvalence}] -row 0 -column 1 -columnspan 1 -sticky nsew
  grid [button $topGui.nucs.buttons.gua -text "G" -command {::Molefacture::add_nuc G ; ::Molefacture::update_openvalence}] -row 0 -column 2 -columnspan 1 -sticky nsew
  grid [button $topGui.nucs.buttons.thy -text "T" -command {::Molefacture::add_nuc T ; ::Molefacture::update_openvalence}] -row 0 -column 3 -columnspan 1 -sticky nsew
  grid [button $topGui.nucs.buttons.ura -text "U" -command {::Molefacture::add_nuc U ; ::Molefacture::update_openvalence}] -row 0 -column 4 -columnspan 1 -sticky nsew

  pack $topGui.nucs.buttons -side top
  pack $topGui.nucs -side top

  # Alternative: build from a sequence
  frame $topGui.buildseq
  button $topGui.buildseq.set_parent_hyd -text "Set parent hydrogen" -command ::Molefacture::set_nuc_parent
  label $topGui.buildseq.label -text "Add a sequence: "
  entry $topGui.buildseq.seq -textvar [namespace current]::nucseq
  button $topGui.buildseq.go -text "Build" -command {::Molefacture::build_nuctextseq $::Molefacture::nucseq; set nucseq ""}

  pack $topGui.buildseq.set_parent_hyd -side left
  pack $topGui.buildseq.go -side right
  pack $topGui.buildseq.seq -side right
  pack $topGui.buildseq.label -side right -fill x
  pack $topGui.buildseq -side top

  labelframe $topGui.type -bd 2 -relief ridge -text "Structure Type" -padx 1m -pady 1m
  frame $topGui.type.buttons
  radiobutton $topGui.type.buttons.rna -text "RNA" -width 6 -value "RNA" -variable [namespace current]::nuctype
  radiobutton $topGui.type.buttons.ssdna -text "ssDNA" -width 6 -value "DNA" -variable [namespace current]::nuctype
  radiobutton $topGui.type.buttons.adna -text "ADNA" -width 6 -value "ADNA" -variable [namespace current]::nuctype
  radiobutton $topGui.type.buttons.bdna -text "BDNA" -width 6 -value "BDNA" -variable [namespace current]::nuctype
  radiobutton $topGui.type.buttons.zdna -text "ZDNA" -width 6 -value "ZDNA" -variable [namespace current]::nuctype

  pack $topGui.type.buttons.rna $topGui.type.buttons.ssdna $topGui.type.buttons.adna $topGui.type.buttons.bdna $topGui.type.buttons.zdna -side left -fill x
  pack $topGui.type.buttons
  pack $topGui.type -side top

  labelframe $topGui.struct -bd 2 -relief ridge -text "Dihedral Angles" -padx 1m -pady 1m
  frame $topGui.struct.ss
  frame $topGui.struct.ss.row1
  frame $topGui.struct.ss.row2
  frame $topGui.struct.ss.row3
  frame $topGui.struct.ss.row4
  frame $topGui.struct.ss.row5
  frame $topGui.struct.ss.row6
  frame $topGui.struct.ss.row7
  label $topGui.struct.ss.row1.alphalabel -text "Alpha angle: "
  entry $topGui.struct.ss.row1.alpha -textvar [namespace current]::ralpha
  label $topGui.struct.ss.row2.betalabel -text "Beta angle: "
  entry $topGui.struct.ss.row2.beta -textvar [namespace current]::rbeta
  label $topGui.struct.ss.row3.gammalabel -text "Gamma angle: "
  entry $topGui.struct.ss.row3.gamma -textvar [namespace current]::rgamma
  label $topGui.struct.ss.row4.deltalabel -text "Delta angle: "
  entry $topGui.struct.ss.row4.delta -textvar [namespace current]::rdelta
  label $topGui.struct.ss.row5.epsilonlabel -text "Epsilon angle: "
  entry $topGui.struct.ss.row5.epsilon -textvar [namespace current]::repsilon
  label $topGui.struct.ss.row6.chilabel -text "Chi angle: "
  entry $topGui.struct.ss.row6.chi -textvar [namespace current]::rchi
  label $topGui.struct.ss.row7.zetalabel -text "Zeta angle: "
  entry $topGui.struct.ss.row7.zeta -textvar [namespace current]::rzeta

  frame $topGui.struct.im
  canvas $topGui.struct.im.pic -height 250 -width 152
  image create photo $topGui.struct.im.dihedpic -format gif -file [file join $nucpath dna_dihedrals.gif]
  $topGui.struct.im.pic create image 1 1 -anchor nw -image $topGui.struct.im.dihedpic

  pack $topGui.struct.ss.row1.alphalabel $topGui.struct.ss.row1.alpha -side left
  pack $topGui.struct.ss.row2.betalabel $topGui.struct.ss.row2.beta -side left
  pack $topGui.struct.ss.row3.gammalabel $topGui.struct.ss.row3.gamma -side left
  pack $topGui.struct.ss.row4.deltalabel $topGui.struct.ss.row4.delta -side left
  pack $topGui.struct.ss.row5.epsilonlabel $topGui.struct.ss.row5.epsilon -side left
  pack $topGui.struct.ss.row6.chilabel $topGui.struct.ss.row6.chi -side left
  pack $topGui.struct.ss.row7.zetalabel $topGui.struct.ss.row7.zeta -side left

  pack $topGui.struct.ss.row1 $topGui.struct.ss.row2 $topGui.struct.ss.row3 $topGui.struct.ss.row4 $topGui.struct.ss.row5 $topGui.struct.ss.row6 $topGui.struct.ss.row7 -side top

  pack $topGui.struct.im.pic -side left

  pack $topGui.struct.ss -side left
  pack $topGui.struct.im -side right

  pack $topGui.struct -side top

}

proc ::Molefacture::prot_builder_gui {} {

  variable phi
  variable psi
  variable aaseq
  variable protBuildWindow

  if {[winfo exists $protBuildWindow]} {
    wm deiconify $protBuildWindow
    raise $protBuildWindow
    return
  }

  set w [toplevel $protBuildWindow]
  wm title $w "Molefacture Protein Builder"
  wm resizable $w no no

  #Frame for buttons for individual amino acids
  labelframe $w.aas -bd 2 -relief ridge -text "Add amino acids" -padx 1m -pady 1m
  frame $w.aas.buttons
  grid [button $w.aas.buttons.ala -text "ALA" -command {::Molefacture::add_aa ALA}] -row 0 -column 0 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.arg -text "ARG" -command {::Molefacture::add_aa ARG}] -row 0 -column 1 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.asn -text "ASN" -command {::Molefacture::add_aa ASN}] -row 0 -column 2 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.asp -text "ASP" -command {::Molefacture::add_aa ASP}] -row 0 -column 3 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.cys -text "CYS" -command {::Molefacture::add_aa CYS}] -row 0 -column 4 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.gln -text "GLN" -command {::Molefacture::add_aa GLN}] -row 0 -column 5 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.glu -text "GLU" -command {::Molefacture::add_aa GLU}] -row 0 -column 6 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.gly -text "GLY" -command {::Molefacture::add_aa GLY}] -row 0 -column 7 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.his -text "HIS" -command {::Molefacture::add_aa HIS}] -row 0 -column 8 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.ile -text "ILE" -command {::Molefacture::add_aa ILE}] -row 0 -column 9 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.leu -text "LEU" -command {::Molefacture::add_aa LEU}] -row 1 -column 0 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.lys -text "LYS" -command {::Molefacture::add_aa LYS}] -row 1 -column 1 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.met -text "MET" -command {::Molefacture::add_aa MET}] -row 1 -column 2 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.phe -text "PHE" -command {::Molefacture::add_aa PHE}] -row 1 -column 3 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.pro -text "PRO" -command {::Molefacture::add_aa PRO}] -row 1 -column 4 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.ser -text "SER" -command {::Molefacture::add_aa SER}] -row 1 -column 5 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.thr -text "THR" -command {::Molefacture::add_aa THR}] -row 1 -column 6 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.trp -text "TRP" -command {::Molefacture::add_aa TRP}] -row 1 -column 7 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.tyr -text "TYR" -command {::Molefacture::add_aa TYR}] -row 1 -column 8 -columnspan 1 -sticky nsew
  grid [button $w.aas.buttons.val -text "VAL" -command {::Molefacture::add_aa VAL}] -row 1 -column 9 -columnspan 1 -sticky nsew

  pack $w.aas.buttons -side top
  pack $w.aas -side top

  # Alternative: build from a sequence
  frame $w.buildseq
  button $w.buildseq.set_parent_hyd -text "Set growing C-Term" -command ::Molefacture::set_prot_cterm
  label $w.buildseq.label -text "Add a sequence: "
  entry $w.buildseq.seq -textvar [namespace current]::aaseq
  button $w.buildseq.go -text "Build" -command {::Molefacture::build_textseq $::Molefacture::aaseq; set aaseq ""}

  pack $w.buildseq.set_parent_hyd -side left
  pack $w.buildseq.go -side right
  pack $w.buildseq.seq -side right
  pack $w.buildseq.label -side right -fill x
  pack $w.buildseq -side top




  labelframe $w.ss -bd 2 -relief ridge -text "Phi/Psi Angles" -padx 1m -pady 1m
  frame $w.ss.row1
  frame $w.ss.row2
  label $w.ss.row1.philabel -text "Phi angle: "
  entry $w.ss.row1.phi -textvar [namespace current]::phi
  button $w.ss.row1.ahel -text "Alpha helix" -width 15 -command {::Molefacture::set_phipsi -57 -47}
  button $w.ss.row1.bs -text "Beta sheet" -width 15 -command {::Molefacture::set_phipsi -120 113}
  label $w.ss.row2.psilabel -text "Psi angle: "
  entry $w.ss.row2.psi -textvar [namespace current]::psi
  button $w.ss.row2.turn -text "Turn" -width 15 -command {::Molefacture::set_phipsi -60 30}
  button $w.ss.row2.straight -text "Straight" -width 15 -command {::Molefacture::set_phipsi -180 180}

  pack $w.ss.row1.philabel $w.ss.row1.phi $w.ss.row1.ahel $w.ss.row1.bs -side left
  pack $w.ss.row2.psilabel $w.ss.row2.psi $w.ss.row2.turn $w.ss.row2.straight -side left

  pack $w.ss.row1 $w.ss.row2 -side top
  pack $w.ss -side top

}

##################################################################################
# Export topology File(s). If the topology retrieved from the CGenFF is still    #
# current (no changes to the molecule after getting the topology                 #
# - Molefacture::topocgenff == 1 if no changes, 0 otherwise), copy the files.    #
# If there is no topology or opt != "nocgenff", try to fetch from CGenFF.        #
# if opt == "nocgenff", use information stored in VMD to produce the topology.   #
##################################################################################
proc ::Molefacture::write_topology_gui { {opt cgenff} } {
   variable tmpmolid
   global env
   # fix_changes

   set types {
           {{CHARMm topology file} {.top} }
           {{CHARMm Stream file} {.str} }
           }


   set defaultext ".top"
   if {$Molefacture::topocgenff == 1} {
      set defaultext ".str"
      set types {
        {{CHARMm Stream file} {.str} }
        {{CHARMm topology file} {.top} }
     }

   }
   set filename [tk_getSaveFile -parent $Molefacture::topGui -title "Save Topology As" -filetypes $types -defaultextension $defaultext]

   if {$filename == ""} {
      return
   }
   set outputtop ""
   if {$Molefacture::topocgenff == 0 && $opt != "nocgenff"} {
      set outpref "tmp"
      set pwd [pwd]
      cd $env(TMPDIR)
      if {[file exist molefacture_tmp] != 1} {
        file mkdir molefacture_tmp
      }
      cd molefacture_tmp
      set topfile [Molefacture::submitCGENFF $outpref]
      cd $pwd

      if {$topfile == -1} {
         set proceed [::Molefacture::unique_atomnames]
         if {$proceed == "yes"} {

            if {$filename != ""} {
               write_topfile $filename "$Molefacture::notDumSel"
               set outputtop [pwd]/$filename
            }

         }
      } elseif {$topfile == -2} {
         return
      }
   }

   ## If the topology was already generated by the CGenFF and
   ## the structure was not changed.
   if {$Molefacture::topocgenff == 1 && $opt != "nocgenff"} {
      set filename_orig $filename

      if {$Molefacture::topoinfo != ""} {
         foreach top $Molefacture::topoinfo {
            set topfile [lindex $top 0]

            if {[llength $Molefacture::topoinfo] > 1} {
               set filename "[file root $filename_orig]_[lindex $top 6][file extension $filename_orig]"
            }
            set outputtop ${filename}
            file copy -force ${topfile} ${filename}
         }
      } else {
         if {[llength $Molefacture::selTop] > 1 } {
            set topfile [lindex $Molefacture::selTop end]

            set outputtop ${filename}
            file copy -force ${topfile} ${filename}
         }
      }


   }

   return $outputtop
}
################################################################
### Proc to prepare the structure to to the cgenff server    ###
### to produce the topology and parameters to be simulated.  ###
### Currently the export and import has to be done separated ###
### and explicitly since we don't have connection to cgenff  ###
### server.                                                  ###
### opt == export - export to mol2                           ###
### opt == import - import str+mol2 file                     ###
################################################################
proc ::Molefacture::run_cgenff { opt {filename ""} {format mol2}} {

   variable tmpmolid
   variable topGui
   variable unknownRes
   if {$opt == "export"} {
      set text "($Molefacture::notDumSel)"
      if {[llength $unknownRes] > 0} {
         append text " and resname [join $unknownRes]"
      }
      set tmpsel [atomselect $tmpmolid $text]
      set resnames [lsort -unique -dictionary [$tmpsel get resname]]
      set residues  [lsort -unique -integer [$tmpsel get residue]]
      set filewntitle "Save Molecule"
      if {[llength $resnames] > 1} {
         set answer [tk_messageBox -message "The molecule contains more than one Residue Name. Do you\
         want to treat them separately?" -icon info -type yesnocancel -parent $topGui]

         if {$answer == "no"} {
            tk_messageBox -message "Please assign an unique Residue Name to the molecule." -icon info -type ok -parent $Molefacture::topGui
            $tmpsel delete -parent $topGui
            return 1
         } elseif {$answer == "cancel"} {
            return 1
         }
         set residues [list]
         set listaux $resnames
         foreach name $resnames {
            set resnamesel [atomselect $tmpmolid "resname $name"]
            lappend residues [lindex [$resnamesel get residue] 0]
            $resnamesel delete
         }

         set filewntitle "Topo Files Prefix"
      } elseif {[llength $resnames] == 1} {
        set residues [lindex $residues 0]
      }

      if {$filename == ""} {
         set filetypes {{{MOL2 Files} {.mol2}}}
         if {$format == "xyz"} {
            set filetypes {{{XYZ Files} {.xyz}}}
         }
         set filename [tk_getSaveFile -title $filewntitle -title "Save Structure As" -parent $topGui -filetypes $filetypes -defaultextension ".$format"]
      }
      if {[llength $filename] == 0} {
         $tmpsel delete
         return 1
      }
      # fix_changes
      set proceed [::Molefacture::unique_atomnames]
      if {$proceed != "yes"} {
         $tmpsel delete
         return 1
      }
      
      set seltext "($Molefacture::notDumSel) and resname"
      set increlement $resnames
      $tmpsel delete
      set filename_orig $filename
      set molfiles [list]
      foreach ele $increlement residue $residues {
         set tmpsel [atomselect $tmpmolid "$seltext $ele and residue $residue"]
         if { [$tmpsel num] > 0 } {
            if {[llength $increlement] > 1} {
               set filename "${ele}[file extension $filename_orig]"
            } else {
               set filename $filename_orig
            }
            Molefacture::exportFile $tmpsel $filename $format
            lappend molfiles $filename
         }
         $tmpsel delete
      }
      update_openvalence
      return $molfiles
   } elseif {$opt == "import"} {
      # set Molefacture::topoinfo [list]
      set filenamelist [list]
      if {$filename == ""} {
          set filenamelist [tk_getOpenFile -title "Import Topology File" -filetypes {{{Stream Files} {.str}} {{Topology Files} {.rtf}} {{All Files} {.*}}} -defaultextension ".str" -multiple 1]
      } else {
         set filenamelist $filename
      }
      foreach filename $filenamelist {
         if {$filename != ""} {
            if {[file dirname $filename] == "."} {
               set filename "[pwd]/$filename"
            }

            set handler [::Toporead::read_charmm_topology $filename 1]
            set topo [::Toporead::topology_from_handler $handler]

            set resnamelist [Toporead::topology_get resnames $topo]
            
            set lp [::Toporead::topology_resi_contains_type $topo $resnamelist "LP*"]
            
            if {$lp != -1} {
              tk_messageBox -message "Molefacture is not yet compatible with Lone \
              Pairs." -title "Import Topology Error failure" -icon error -type ok \
              -parent $Molefacture::topGui
              return -2
            }

            lappend Molefacture::topoinfo $topo


            set row 0
            
            set resindex [lsearch $unknownRes $resnamelist]
            if {$resindex != -1} {
               set unknownRes [lreplace $unknownRes $resindex $resindex]
            }
            set atmsnamelist [list]
            set list [Toporead::topology_get_resid $topo $resnamelist atomlist]
            foreach atom $list {
               lappend atmsnamelist [lindex $atom 0]
            }

            set atmtype [list]
            set atmcharge [list]
            set chargepenalty [list]
            foreach name $atmsnamelist {
               set atminfo [::Toporead::topology_resi_contains_atom $topo $resnamelist $name]
               lappend atmtype [lindex $atminfo 1]
               lappend atmcharge [lindex $atminfo 2]
               lappend chargepenalty [lindex $atminfo 3]
            }

            set tmpsel [atomselect $tmpmolid "resname [join $resnamelist]"]

            $tmpsel set name $atmsnamelist
            $tmpsel set type $atmtype
            $tmpsel set charge $atmcharge
            $tmpsel set $Molefacture::modFlagPSF 1
            ### Store the charge penalties in the user field
            if {[llength $chargepenalty] > 0} {
               $tmpsel set user $chargepenalty
            }
            $tmpsel delete
            ::Molefacture::updateAtomTable
            set Molefacture::topocgenff 1
         }
      }
      ::Molefacture::clearTableSelection
   }
}

################################################################
### Proc to export the structure and str (top) file to ffTK  ###
################################################################
proc ::Molefacture::export_ffTK {} {
   package require forcefieldtoolkit

   global env

   ### Stream file generated by CGenFF
   set files [Molefacture::write_psfpdbfiles]

   set psf [lindex $files 0]
   set pdb [lindex $files 1]
   set top [lindex $files 2]
   
   if {[llength $Molefacture::topoinfo] > 0} {
     set top [lindex $Molefacture::topoinfo 0 0]
   } 
   

   if {$pdb == "" || $files == ""} {
      return
   }

   fftk
   set ::ForceFieldToolKit::BuildPar::idMissingPSF $psf
   set ::ForceFieldToolKit::BuildPar::idMissingPDB $pdb

   set ForceFieldToolKit::BuildPar::cgenffMol $pdb
   set ForceFieldToolKit::BuildPar::cgenffStr $top

}


################################################################
### Proc to export the structure to a mol2 or xyz file       ###
################################################################
proc ::Molefacture::exportFile { tmpsel filename format} {
   variable topGui


   ### The package on the "else" part of the statment is also 
   ### implemented in the cgenffcaller plugin as one need to ensure
   ### that the files are compatible with CgenFF. Once the cgenffcaller
   ### is distributed with VMD, the "else" statement can be deleted

   if {$Molefacture::cgenffcaller == 1} {
      #puts "tmpsel $tmpsel || filename $filename || format $format "
      if {[::CGENFFCALLER::saveMol $tmpsel $filename $format] == -1} {
         tk_messageBox -message "Residue's Name Not Found. Please Attribute a Name to the Residue." \
          -title "Residue's Name Not Found" -icon error -type ok -parent $topGui
          $tmpsel delete
          return -1
      }
   } else {
      set resname [lsort -unique [$tmpsel get resname]]
      if {[llength $resname] == 0} {
          tk_messageBox -message "Residue's Name Not Found. Please Attribute a Name to the Residue." \
          -title "Residue's Name Not Found" -icon error -type ok -parent $topGui
          $tmpsel delete
          return -1
      }
      ### set the element as the atom's type. This tends to produce
      ### better results with the cgenff, as it eliminates confusion
      ### with atom types
      set typelist [$tmpsel get type]
      set elementlist [$tmpsel get element]
      $tmpsel set type $elementlist
      if {$format == "mol2"} {
         $tmpsel writemol2 $filename
      } elseif {$format == "xyz"} {
         $tmpsel writexyz $filename
      }

      $tmpsel set type $typelist

      if {[file exists $filename] == 1 && $format == "mol2"} {
         file copy -force ${filename} ${filename}_bkup
         set filein [open ${filename}_bkup r]
         set fileout [open ${filename} w+]
         while {[eof $filein] != 1} {
            set line [gets $filein]
            ### Send the "generated by VMD" line to the end of the file
            ### as a comment to ensure compatibility with CGenFF
            if {$line == "generated by VMD"} {
               puts $fileout "$resname"
               seek $filein [tell $filein]
               set lines [read $filein]
               puts $fileout $lines
               puts $fileout "\#$line"
               break
            } else {
               puts $fileout $line
            }
         }
         close $fileout
         close $filein
         file delete -force ${filename}_bkup
      }
   }
   

   
}

##############################################################
# Deprecated function previously used to assign atoms' types #
##############################################################

proc ::Molefacture::run_idatm {} {
  # fix_changes
  variable tmpmolid

  set tmpsel [atomselect $tmpmolid "occupancy > 0.5 and not resname DUM"]

  set type [$tmpsel get type]
  if { [$tmpsel num] > 0 } {
   set output ""
   if {[catch {::IDATM::runtyping $tmpsel} output] == 1} {
      tk_messageBox -message "The automatic assignment of bond order failed. Please revise the bond order." -title "Bond Order Failed" -parent $Molefacture::topGui -icon error -type ok
      puts "IDATM error $output"
   }

  }
  $tmpsel set type $type
  $tmpsel delete
  update_openvalence
}

proc ::Molefacture::run_ante {} {
  fix_changes
  variable tmpmolid
  variable totalcharge

  set tmpsel [atomselect $tmpmolid "occupancy >= 0.8"]
  
  if { [$tmpsel num] > 0 } {
    set result [::ANTECHAMBER::ac_type_in_place $tmpsel rc $totalcharge]
    if { $result != "0" } {
        tk_messageBox -message $result -type ok -title "Antechamber Error - Check Structure" -icon error
    }
  }
  $tmpsel delete
  update_openvalence
  draw_openvalence
  draw_selatoms
  clean_dihed_tags
}

proc ::Molefacture::am1_geo_opt {} {
  fix_changes
  variable tmpmolid
  variable totalcharge
  set tmpsel [atomselect $tmpmolid "occupancy >= 0.8"]
  if { [$tmpsel num] > 1 } {
    ::runsqm::SQMopt $tmpsel $totalcharge
  }
  
  update_openvalence
  draw_openvalence
  draw_selatoms
  clean_dihed_tags
}



proc ::Molefacture::run_ante_opls {} {
  fix_changes
  variable tmpmolid
  variable totalcharge

  set tmpsel [atomselect $tmpmolid "occupancy >= 0.8"]
  
  if { [$tmpsel num] > 1 } {
    set result [::ANTECHAMBER::ac_type_in_place $tmpsel rc $totalcharge "CUSTOM$::Molefacture::OPLSatomdef"]
    if { $result != "0" } {
        tk_messageBox -message $result -type ok -title "Antechamber Error - Check Structure" -icon error
    }
  }
  $tmpsel delete
  update_openvalence
  draw_openvalence
  draw_selatoms
  clean_dihed_tags
}

proc ::Molefacture::ante_gui {} {
  fix_changes
  ::ANTECHAMBER::antechamber_gui 1
  update_openvalence
}

proc ::Molefacture::fep_simulations_gui {} {
  variable selectcolor
  variable fepmolid
  variable tmpmolid
  set w [toplevel ".fepsimul"]
  wm title $w "FEP Simulations options"
  wm resizable $w no no

  labelframe $topGui.buttons -bd 2 -relief ridge -text "Actions" -padx 2 -pady 2
  button $topGui.buttons.incomming -text "Define incoming atoms (g)" -command ::Molefacture::flag_incomingatoms
  button $topGui.buttons.outgoing -text "Define outgoing atoms (v)" -command ::Molefacture::flag_outgoingatoms
  button $topGui.buttons.alchemify -text "Run Alchemify" -command ::Molefacture::run_alchemify
  button $topGui.buttons.common -text "Clear (x)" -command ::Molefacture::flag_commonatoms
  pack $topGui.buttons.incomming $topGui.buttons.outgoing $topGui.buttons.common $topGui.buttons.alchemify -side left  -fill x -expand 1
  pack $topGui.buttons  -side top -fill x

  labelframe $topGui.options -bd 2 -relief ridge -text "Typing and charges"
  frame $topGui.options.line1
  label $topGui.options.line1.label -text "Typing:"
  radiobutton $topGui.options.line1.notypingopt -text "None" -value "None" -variable [namespace current]::feptyping
  radiobutton $topGui.options.line1.gafftypingopt -text "GAFF" -value "GAFF" -variable [namespace current]::feptyping
  radiobutton $topGui.options.line1.oplstypingopt -text "OPLS" -value "OPLS" -variable [namespace current]::feptyping
  pack $topGui.options.line1.label $topGui.options.line1.notypingopt $topGui.options.line1.gafftypingopt $topGui.options.line1.oplstypingopt -side left -fill x -expand 1
  pack $topGui.options.line1 -side left -fill x -expand 1

  frame $topGui.options.line2
  label $topGui.options.line2.sclabel -text "Starting charge:"
  entry $topGui.options.line2.sc -textvar [namespace current]::fepstartcharge -width 5
  label $topGui.options.line2.eclabel -text "Ending charge:"
  entry $topGui.options.line2.ec -textvar [namespace current]::fependcharge -width 5

  pack $topGui.options.line2.sclabel $topGui.options.line2.sc $topGui.options.line2.eclabel $topGui.options.line2.ec -side left -fill x -expand 1
  pack $topGui.options.line2 -side left -fill x -expand 1

  pack $topGui.options -side top -fill x

  labelframe $topGui.saveoptions -bd 2 -relief ridge -text "Output options" -padx 2 -pady 2
  frame $topGui.saveoptions.line1
  label $topGui.saveoptions.line1.outpreflabel -text "Output prefix:"
  entry $topGui.saveoptions.line1.outpref -textvar [namespace current]::FEPoutprefix
  checkbutton $topGui.saveoptions.line1.mpflag -text "Merge into parent molecule" -variable [namespace current]::FEPreplaceflag
  pack $topGui.saveoptions.line1.outpreflabel $topGui.saveoptions.line1.outpref $topGui.saveoptions.line1.mpflag -side left -fill x -expand 1
  pack $topGui.saveoptions.line1 -side top -fill x
  frame $topGui.saveoptions.line2
  label $topGui.saveoptions.line2.parlabel -text "Parent prefix:"
  entry $topGui.saveoptions.line2.parmolid -textvar [namespace current]::FEPparentmol
  label $topGui.saveoptions.line2.seglabel -text "Segname:"
  entry $topGui.saveoptions.line2.segname -textvar [namespace current]::FEPdelseg
  label $topGui.saveoptions.line2.reslabel -text "Resid:"
  entry $topGui.saveoptions.line2.resid -textvar [namespace current]::FEPdelres
  pack $topGui.saveoptions.line2.parlabel $topGui.saveoptions.line2.parmolid $topGui.saveoptions.line2.seglabel $topGui.saveoptions.line2.segname $topGui.saveoptions.line2.reslabel $topGui.saveoptions.line2.resid -side left -fill x -expand 1
  pack $topGui.saveoptions.line2 -side top -fill x
  pack $topGui.saveoptions -side top -fill x



  labelframe $topGui.betalist -bd 2 -relief ridge -text "FEP selections" -padx 2 -pady 2

  label $topGui.betalist.format -font {tkFixed 9} -textvariable ::Molefacture::FEPlist -relief flat -bd 2 -justify left;
  scrollbar $topGui.betalist.scroll -command "$topGui.betalist.list yview" -takefocus 0
  listbox $topGui.betalist.list -activestyle dotbox -yscroll "$topGui.betalist.scroll set" -font {tkFixed 9} \
  -width 60 -height 5 -setgrid 1 -selectmode extended -selectbackground $selectcolor -listvariable ::Molefacture::FEPatomlist
  pack $topGui.betalist.format -side top -anchor w
  pack $topGui.betalist.list -side left -fill both -expand 1
  pack $topGui.betalist.scroll -side left -fill y
  pack $topGui.betalist -padx 1m -pady 1m -fill both -expand 1

  set fepmolid [mol new]
  graphics $fepmolid material Transparent
  mol top $tmpmolid
  display resetview

#  Hotkeys forme $topGui.val.list

  user add key g {
    ::Molefacture::flag_incomingatoms
  }
  user add key v {
    ::Molefacture::flag_outgoingatoms
  }
  user add key x {
	 ::Molefacture::flag_commonatoms
  }
}

proc ::Molefacture::run_psfgen_FEP {} {

  variable FEPreplaceflag
  variable FEPdelseg

  if {$FEPreplaceflag > 0} {
    set segn $FEPdelseg
  } else {
    set segn U
  }

   package require psfgen
   psfcontext new delete
   topology "Molefacture_fep_temp.top"
   segment $segn {pdb Molefacture_fep_temp_combined.pdb }
   coordpdb Molefacture_fep_temp_combined.pdb $segn
   writepsf Molefacture_fep_temp_combined.psf

   set tmpmol [mol new Molefacture_fep_temp_combined.pdb]
   set tmpsel [atomselect $tmpmol all]
   $tmpsel set segname $segn
   $tmpsel writepdb Molefacture_fep_temp_combined.pdb
   $tmpsel delete
   mol delete $tmpmol
}

proc ::Molefacture::run_alchemify {} {
   variable tmpmolid
   variable filename_FEP
   variable FEPreplaceflag
   variable FEPoutprefix
   variable FEPdelseg
   variable FEPdelres
   variable FEPparentmol

   mol reanalyze $tmpmolid
   do_fep_typing
   run_psfgen_FEP
   alchemify Molefacture_fep_temp_combined.psf Molefacture_fep_ready.psf Molefacture_fep_temp_combined.pdb

   if {$FEPreplaceflag != 0} {
     mol new Molefacture_fep_ready.psf
     mol addfile Molefacture_fep_temp_combined.pdb

     package require psfgen
     psfcontext new delete

     readpsf ${FEPparentmol}.psf
     coordpdb ${FEPparentmol}.pdb

     delatom ${FEPdelseg} ${FEPdelres}

     readpsf Molefacture_fep_ready.psf
     coordpdb Molefacture_fep_temp_combined.pdb

     writepsf ${FEPoutprefix}_temp.psf
     writepdb ${FEPoutprefix}.pdb

     set fepatoms [atomselect top all]
     mol new ${FEPoutprefix}.pdb
     set sel [atomselect top all]
     $sel set beta 0
     set fepnewatoms [atomselect top "segname ${FEPdelseg} and resid ${FEPdelres}"]
     $fepnewatoms set beta [$fepatoms get beta]
     $fepnewatoms set occupancy [$fepatoms get occupancy]

     $sel writepdb ${FEPoutprefix}.pdb

     file rename -force Molefacture_fep_temp.top ${FEPoutprefix}.top

     alchemify ${FEPoutprefix}_temp.psf ${FEPoutprefix}.psf ${FEPoutprefix}.pdb

     $sel delete
     $fepatoms delete
     $fepnewatoms delete

     mol delete top
     mol delete top

     file delete -force ${FEPoutprefix}_temp.psf

   } else {
     file rename -force Molefacture_fep_ready.psf ${FEPoutprefix}.psf
     file rename -force Molefacture_fep_temp_combined.pdb ${FEPoutprefix}.pdb
     file rename -force Molefacture_fep_temp.top ${FEPoutprefix}.top

     atom_picked_fctn
     update_openvalence
     draw_openvalence
   }
}


proc ::Molefacture::runante_am1bcc {} {
  fix_changes
  variable tmpmolid
  variable totalcharge

  set tmpsel [atomselect $tmpmolid "occupancy >= 0.8"]
  #::IDATM::runtyping $tmpsel
  if { [$tmpsel num] > 0 } {
    ::Molefacture::run_sqm_mullik AM1
    set result [::ANTECHAMBER::getAM1BCC $tmpsel bcc $totalcharge none ]
    if { $result != "0" } {
        tk_messageBox -message $result -type ok -title "Antechamber Error - Check Structure" -icon error
    }
  }
  $tmpsel delete
  update_openvalence
  draw_openvalence
  draw_selatoms
  clean_dihed_tags
}


proc ::Molefacture::run_sqm_mullik { {qmlvl ""} } {
  fix_changes
  variable tmpmolid
  variable totalcharge
  set curqmtheory $::runsqm::qmtheory
  if { $qmlvl !=  "" } {
    variable ::runsqm::qmtheory $qmlvl
  }
  set tmpsel [atomselect $tmpmolid "occupancy >= 0.5"]
  if { [$tmpsel num] > 0 } {
      ::runsqm::get_charges $tmpsel $totalcharge "Mulliken"
  }
  variable ::runsqm::qmtheory $curqmtheory
  $tmpsel delete
  update_openvalence
  return
}

proc ::Molefacture::run_sqm_cm1a {} {
  fix_changes
  variable tmpmolid
  variable totalcharge
  set curqmtheory $::runsqm::qmtheory
  variable ::runsqm::qmtheory "AM1"
  set tmpsel [atomselect $tmpmolid "occupancy >= 0.5"]
  if { [$tmpsel num] > 0 } {
    ::runsqm::get_charges $tmpsel $totalcharge "cm1a"
  }
  variable ::runsqm::qmtheory $curqmtheory
  $tmpsel delete
  update_openvalence
  return
}

proc ::Molefacture::antechamber_setting_gui { } {
  variable topGui

  if { [winfo exists .antechambersettingsgui] } {
    wm deiconify .antechambersettingsgui
    raise .antechambersettingsgui
    return
  }


  set w [toplevel ".antechambersettingsgui"]
  wm title $w "Antechamber Settings"
  wm resizable $w 0 0

  set rownum 0
  frame $topGui.settings

  grid [label $topGui.settings.amberhomelabel -text "AMBERHOME:"] -row $rownum -column 0 -sticky w
  grid [entry $topGui.settings.maxcyc -width 20 -textvar env(AMBERHOME)] \
    -row $rownum -column 1
  grid [button $topGui.settings.browse -text "Browse" -command ::Molefacture::find_amberhome] -row $rownum -column 2 -columnspan 1 -sticky w
  incr rownum

  grid [label $topGui.settings.oplsdefloclabel -text "OPLS Atom Typing Rules:"] -row $rownum -column 0 -sticky w
  grid [entry $topGui.settings.oplsdefloc -width 20 -textvar ::Molefacture::OPLSatomdef] \
    -row $rownum -column 1
  grid [button $topGui.settings.oplsbrowse -text "Browse" -command ::Molefacture::find_OPLSfile] -row $rownum -column 2 -columnspan 1 -sticky w
  incr rownum

  grid [button $topGui.settings.done -text "Done" -command ::Molefacture::set_antechamber_settings] -row $rownum -column 2 -columnspan 1 -sticky e
  pack $topGui.settings

}

proc ::Molefacture::find_amberhome {} {
  set filename [tk_chooseDirectory -title "Choose a directory" -initialdir "$::env(HOME)"]
#  puts $filename
  set ::env(AMBERHOME) "$filename"
}

proc ::Molefacture::find_sqm {} {
  set filename [tk_getOpenFile -title "Please locate sqm executable" -initialdir "$::env(HOME)"]
#  puts $filename
  set ::runsqm::sqmpath "$filename"
}

proc ::Molefacture::find_OPLSfile {} {
  set filename [tk_getOpenFile -title "Select OPLS atom typing .DEF file" -initialdir [file dirname $::Molefacture::OPLSatomdef]]
#  puts $filename
  set ::Molefacture::OPLSatomdef "$filename"
}

proc ::Molefacture::set_antechamber_settings {} {
    if { [info exist ::env(AMBERHOME)] == 1 && [string length $::env(AMBERHOME)] > 0 } {
        set  ::ANTECHAMBER::acpath \
            [::ExecTool::find -interactive \
            -description "Antechamber" \
            -path [file join $::env(AMBERHOME) bin $::ANTECHAMBER::acbin ] $::ANTECHAMBER::acbin]
    } else {
       set  ::ANTECHAMBER::acpath \
            [::ExecTool::find -interactive \
            -description "Antechamber" $::ANTECHAMBER::acbin]
    }
   ::ANTECHAMBER::set_oplsdef_loc $::Molefacture::OPLSatomdef
   ::Molefacture::checkAntechamber
   destroy .antechambersettingsgui
}

proc ::Molefacture::sqm_settings_gui { } {
  variable topGui

  if { [winfo exists .sqmsettingsgui] } {
    wm deiconify .sqmsettingsgui
    raise .sqmsettingsgui
    return
  }

  set w [toplevel ".sqmsettingsgui"]
  wm title $w "SQM Settings"
  wm resizable $w 0 0

  set rownum 0
  frame $topGui.settings
  grid [label $topGui.settings.sqmloclalbel -text "Location of sqm executable:"] -row $rownum -column 0 -sticky w
  incr rownum
  grid [entry $topGui.settings.sqmloc -width 25 -textvar ::runsqm::sqmpath] \
    -row $rownum -column 0 -columnspan 1
  grid [button $topGui.settings.browse -text "Browse" -command ::Molefacture::find_sqm] -row $rownum -column 1 -columnspan 3 -sticky w
  incr rownum

  grid [label $topGui.settings.qmtheorylabel -text "QM Theory:"] -row $rownum -column 0 -sticky w
  grid [menubutton $topGui.settings.qmtheory -menu $topGui.settings.qmtheory.menu -textvar ::runsqm::qmtheory -relief raised -width 10] \
    -row $rownum -column 1 -columnspan 3 -sticky w
  menu $topGui.settings.qmtheory.menu -tearoff no
  $topGui.settings.qmtheory.menu add radiobutton -label "AM1" -variable ::runsqm::qmtheory -value "AM1"
  $topGui.settings.qmtheory.menu add radiobutton -label "PM3" -variable ::runsqm::qmtheory -value "PM3"
  $topGui.settings.qmtheory.menu add radiobutton -label "MNDO" -variable ::runsqm::qmtheory -value "MNDO"
  $topGui.settings.qmtheory.menu add radiobutton -label "RM1" -variable ::runsqm::qmtheory -value "RM1"
  $topGui.settings.qmtheory.menu add radiobutton -label "PM3-PDDG" -variable ::runsqm::qmtheory -value "PM3-PDDG"
  $topGui.settings.qmtheory.menu add radiobutton -label "MNDO-PDDG" -variable ::runsqm::qmtheory -value "MNDO-PDDG"
  $topGui.settings.qmtheory.menu add radiobutton -label "PM3-CARB1" -variable ::runsqm::qmtheory -value "PM3-CARB1"
  $topGui.settings.qmtheory.menu add radiobutton -label "PM6" -variable ::runsqm::qmtheory -value "PM6"
  incr rownum
  grid [label $topGui.settings.maxcyclabel -text "Max minimisation cycles:"] -row $rownum -column 0 -sticky w
  grid [entry $topGui.settings.maxcyc -width 4 -textvar ::runsqm::maxcyc] \
    -row $rownum -column 1 -sticky ew
  incr rownum
  grid [label $topGui.settings.printchargeslabel -text "Print charges during minimisation:"] -row $rownum -column 0 -sticky ew
  grid [checkbutton $topGui.settings.printcharges -variable ::runsqm::printcharges] \
    -row $rownum -column 1 -sticky ew
  incr rownum
  grid [label $topGui.settings.ntprlabel -text "Number of steps between printing\nminimisation output:" -justify left] -row $rownum -column 0 -sticky w
  grid [entry $topGui.settings.ntpr -width 3 -textvar ::runsqm::ntpr] \
    -row $rownum -column 1 -sticky ew
  incr rownum
  grid [label $topGui.settings.peptidecorrlabel -text "Apply an MM correction to peptide\nlinkages:" -justify left] -row $rownum -column 0 -sticky w
  grid [checkbutton $topGui.settings.peptidecorr -variable ::runsqm::peptidecorr] \
    -row $rownum -column 1 -sticky ew
  incr rownum
  grid [label $topGui.settings.verbositylabel -text "Verbosity Level"] -row $rownum -column 0 -sticky w
  grid [menubutton $topGui.settings.verbosity -menu $topGui.settings.verbosity.menu -textvar ::runsqm::verbosity -relief raised -width 2] \
    -row $rownum -column 1 -columnspan 1 -sticky ew
  menu $topGui.settings.verbosity.menu -tearoff no
  $topGui.settings.verbosity.menu add radiobutton -label "0" -variable ::runsqm::verbosity -value "0"
  $topGui.settings.verbosity.menu add radiobutton -label "1" -variable ::runsqm::verbosity -value "1"
  $topGui.settings.verbosity.menu add radiobutton -label "2" -variable ::runsqm::verbosity -value "2"
  $topGui.settings.verbosity.menu add radiobutton -label "3" -variable ::runsqm::verbosity -value "3"
  $topGui.settings.verbosity.menu add radiobutton -label "4" -variable ::runsqm::verbosity -value "4"
  $topGui.settings.verbosity.menu add radiobutton -label "5" -variable ::runsqm::verbosity -value "5"
  incr rownum
  grid [label $topGui.settings.itrmaxlabel -text "Maximin number of scf iterations:" -justify left] -row $rownum -column 0 -sticky w
  grid [entry $topGui.settings.itrmax -width 3 -textvar ::runsqm::itrmax] \
    -row $rownum -column 1 -sticky ew
  incr rownum


  grid [label $topGui.settings.scfconvlabel -justify left -text "SCF convergence criteria:"] -row $rownum -column 0 -sticky w
  grid [entry $topGui.settings.scfconvcoeff -width 2 -textvar ::runsqm::scfconvcoeff] \
    -row $rownum -column 1 -sticky w
  grid [label $topGui.settings.scfconvexp -justify left -text "x10^"] -row $rownum -column 1 -sticky e
  grid [scale $topGui.settings.scfconvscale -orient horizontal -length 100 -from -2 -to -14 -variable runsqm::scfconvexp  -tickinterval -4] -row $rownum -column 2 -columnspan 1 -sticky ew
  incr rownum

  grid [label $topGui.settings.scferrabel -justify left -text "Maximum absolute SCF error:"] -row $rownum -column 0 -sticky w
  grid [entry $topGui.settings.scferrcoeff -width 2 -textvar ::runsqm::errconvcoeff] \
    -row $rownum -column 1 -sticky w
  grid [label $topGui.settings.scfcerrexp -justify left -text "x10^"] -row $rownum -column 1 -sticky e
  grid [scale $topGui.settings.scferrscale -orient horizontal -length 100 -from -1 -to -16 -variable runsqm::errconvexp  -tickinterval -5] -row $rownum -column 2 -columnspan 1 -sticky ew
  incr rownum

  grid [button $topGui.settings.defaults -text "Restore Defaults" -command ::runsqm::setdefaults] -row $rownum -column 0 -columnspan 1 -sticky e
  grid [button $topGui.settings.apply -text "Done" -command ::Molefacture::checkSQMpars] -row $rownum -column 1 -columnspan 1
  pack $topGui.settings

}

proc ::Molefacture::checkAntechamber {} {

    set m .molefac.menubar.build
    ::ANTECHAMBER::acinit
    if { $::ANTECHAMBER::achere == 0 } {
        $m.abtype entryconfigure 1 -state disable
        $m.abtype entryconfigure 2 -state disable
        $m.abtype entryconfigure 3 -state disable
        $m.charges entryconfigure 2 -state disable
    } else {
        $m.abtype entryconfigure 1 -state normal
        $m.abtype entryconfigure 2 -state normal
        $m.abtype entryconfigure 3 -state normal
        $m.charges entryconfigure 2 -state normal
    }

}

proc ::Molefacture::checkSQMpars {} {

  set display_error 0
  set msg ""
  set originalpath $::runsqm::sqmpath
  set ::runsqm::sqmhere 0
  ::Molefacture::checkSQM
   if {$::runsqm::sqmhere == 0} {
      set display_error 1
      append msg "Warning: $originalpath not executable - sqm calculations not possible."
   }
   if {$::runsqm::maxcyc == 0} {
      set display_error 1
      append msg "Warning: With maxcyc = 0 only a single point calculation will be run.\n"
   }
   if {$::runsqm::itrmax > 1000} {
      set display_error 1
      append msg "Warning: If SCF convergence hasn't been reached by 1000 steps it is unlikely to converge. It is recommended that the maximum number of SCF steps is set to 1000 or less.\n"
   }
   if {$::runsqm::scfconvcoeff < 1 || $::runsqm::scfconvcoeff > 9} {
     set display_error 2
     append msg "Error: SCF convergence criteria must have a coefficient between 1 and 9.\n"
   }
   if {$::runsqm::errconvcoeff < 1 || $::runsqm::errconvcoeff > 9} {
     set display_error 2
     append msg "Error: Maximum SCF error criteria must have a coefficient between 1 and 9.\n"
   }
   if { $display_error == 2 } {
    tk_messageBox -message $msg -type ok -title "Check SQM parameters" -icon error
   } elseif { $display_error == 1} {
     tk_messageBox -message $msg -type ok -title "Check SQM parameters" -icon warning
     destroy .sqmsettingsgui
   } else {
     destroy .sqmsettingsgui
   }
   return

}

proc ::Molefacture::checkSQM {} {

    set ::runsqm::smqhere [::runsqm::sqminit]
    set m .molefac.menubar.build
    if { $::runsqm::sqmhere == 0 } {
        $m.charges entryconfigure 0 -state disable
        $m.charges entryconfigure 1 -state disable
        $m.menu entryconfigure 2 -state disable
    } else {
        $m.charges entryconfigure 0 -state normal
        $m.charges entryconfigure 1 -state normal
        $m.menu entryconfigure 2 -state normal
    }

}

# proc ::Molefacture::write_namdfiles_gui { } {
#   variable topGui

#   if { [winfo exists .writenamdfilesgui] } {
#     wm deiconify .writenamdfilesgui
#     raise .writenamdfilesgui
#     return
#   }


#   set w [toplevel ".writenamdfilesgui"]
#   wm title $w "Write NAMD input files"
#   wm resizable $w 0 0

#   set rownum 0
#   frame $topGui.writefiles
#   variable prmloc
#   variable prmprefix
#   variable outloc
#   if {![info exists prmloc]} {
#     set prmloc ""
#   }
#   if {![info exists namdfilesprefix]} {
#     set prmprefix "molecule"
#   }
#   if {![info exists prefix]} {
#     set outloc "./"
#   }
#   grid [label $topGui.writefiles.info -text "This tools writes out the topology (top), structure (psf) and coordinates (pdb) files"] -row $rownum -column 0 -columnspan 4 -sticky w
#   #grid [label $topGui.writefiles.info -text "This tools writes out the structure (psf), coordinates (pdb)\nand parameter (prm) files that you need to run NAMD\n"] -row $rownum -column 0 -columnspan 4 -sticky w
#   #grid [button $topGui.settings.browse -text "Browse" -command ::Molefacture::find_amberhome] -row $rownum -column 2 -columnspan 1 -sticky w
#   incr rownum

#   #grid [label $topGui.writefiles.prmloclabel -text "Location of force-field\nparameter file:"] -row $rownum -column 0 -columnspan 2 -sticky w
#   #grid [entry $topGui.writefiles.prmloc -width 20 -textvar ::Molefacture::prmloc ] -row $rownum -column 2
#   #grid [button $topGui.writefiles.prmbrowse -text "Browse" -command ::Molefacture::findprmfile] -row $rownum -column 3 -columnspan 1 -sticky w
#   #incr rownum

#   grid [label $topGui.writefiles.prefixlabel -text "Prefix for output files:"] -row $rownum -column 0 -columnspan 2 -sticky w
#   grid [entry $topGui.writefiles.prefix -width 20 -textvar ::Molefacture::prmprefix ] -row $rownum -column 2
#   #grid [button $topGui.writefiles.prmbrowse -text "Browse" -command ::Molefacture::findprmfile] -row $rownum -column 3 -columnspan 1 -sticky w
#   incr rownum

#   grid [label $topGui.writefiles.dirlabel -text "Directory to write output:"] -row $rownum -column 0 -columnspan 2 -sticky w
#   grid [entry $topGui.writefiles.dirloc -width 20 -textvar ::Molefacture::dirloc ] -row $rownum -column 2
#   grid [button $topGui.writefiles.dirrowse -text "Browse" -command ::Molefacture::findoutputdir] -row $rownum -column 3 -columnspan 1 -sticky w
#   incr rownum

#   grid [button $topGui.writefiles.writefiles -text "Write output files" -command {::Molefacture::write_namdfiles [file join $::Molefacture::dirloc $::Molefacture::prmprefix] $::Molefacture::prmloc }] -row $rownum -column 2 -columnspan 1 -sticky e
#   grid [button $topGui.writefiles.cancel -text "Cancel" -command {destroy .writenamdfilesgui}] -row $rownum -column 3 -columnspan 1 -sticky e
#   pack $topGui.writefiles
# }

proc ::Molefacture::findprmfile {} {
  set filename [tk_getOpenFile -title "Select force-field parameter file"]
#  puts $filename
  set ::Molefacture::prmloc "$filename"
}

proc ::Molefacture::findoutputdir {} {
  set filename [tk_chooseDirectory -title "Choose an output directory" -initialdir "$::env(PWD)"]
#  puts $filename
  set ::Molefacture::dirloc "$filename"
}
