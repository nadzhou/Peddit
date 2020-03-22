package require mdff
package require ssrestraints
package require namdrun
package require tktooltip

package provide mdff_gui 0.2
namespace eval MDFFGUI:: {
  namespace export mdffgui
  namespace eval gui { 
    variable w 
    variable wad 
    variable xwad 
    variable currentMol
    variable nullMolString
    variable nullVolString
    variable CurrentPDBFile
    variable CurrentPSFFile
    variable DensityID
    variable CButtonChi
    variable MapID
    variable xMDFFMapIDs
    variable CButtonCis
    variable CCResInput
    variable CCSelInput
    variable CurrentDXFile
    variable Extrabonds 
    variable CurrentRefs
    variable MolMenuText   
    variable CurrentCrystPDB
    variable xMDFFDensityViewAdd
    variable CurrentDirLabel
    variable env
    variable DensityView
    variable xMDFFDensityView
    variable xMDFFSelInput
        
    variable CCStep
    variable CCThreshInput
    variable UseCCThresh
    variable CCMolInput

    variable IMDSubmit
    variable IMDConnect
    variable IMDPause
    variable IMDFinish
    variable IMDStatusLabel
    variable IMDProcInput   
    variable IMDCCPlot
    variable Notebook
    variable justChanged
    variable GridpdbGen
    variable FixedpdbGen
    variable MDFFSetup
    variable GenerateXBonds
    
    variable CrystOption
    variable CurrentRefsEntry
    variable RefStepsInput
    variable CurrentCrystEntry
    variable GetRefs
    variable GetCryst
    variable CButtonMask
    variable CButtonBFS
    variable MaskResInput ""
    variable MaskCutoffInput ""
   
    variable REMDFFNumReplicasInput
    variable CButtonREMDFFAutoBlur
    
    variable IMDInpOrList
    variable IMDHost
    variable namdbin
    variable CurrentServer 
    variable ServerList
    variable IMDHostInput
    variable IMDServer 
    variable CButtonIMDWait
    variable CButtonIMDIgnore
    variable IMDPortInput
    variable IMDFreqInput
    variable IMDKeepInput
    variable HLF
    variable CCx
    variable CCy
    variable ShowBTNPadX

    variable MapToolsMolMenuText
    variable MapToolsVolMenuText
    variable MapToolsVol2MenuText
    variable MapToolsStructMolMenuText
    variable CCOutput ""
    variable MapToolsMol2MenuText
    #variable HistPlot
  }
  namespace eval settings:: {
    namespace export settings
    variable MolID ""
    variable CurrentPDBPath ""
    variable CurrentPSFPath ""
    variable GridPDBSelList ""
    variable xMDFFGridPDBSelList ""
    variable xMDFFMaskList ""
    variable xMDFFBSharpList ""
    variable xMDFFBSharpFactorList ""
    variable xMDFFAverageMapList ""
    variable xMDFFMaskResList ""
    variable xMDFFMaskCutoffList ""
    variable DensityList ""
    variable RefsList ""
    variable xMDFFSel "noh and not water"
    variable GScaleList ""
    variable ThresholdList ""
    variable xMDFFGScaleList ""
    variable CurrentDXPath ""
    variable CurrentRefsPath ""
    variable CurrentDir [pwd]
    variable CurrentCrystPath ""
    variable Temperature 300
    variable FTemperature 300
    variable Minsteps 200
    variable Numsteps 50000
    variable GScale 0.3
    variable GriddxThreshold 0
    variable GridOff 0
    variable GridPDBSelText "protein and noh"
    variable IMDSelText "all"
    variable FixedPDBSelText "none"
    variable FixedColumn "occupancy"
    variable SimulationName "mdff"
    variable SimulationStep "1"
    variable SSRestraints 0
    variable ChiralityRestraints 0
    variable CispeptideRestraints 0
    variable REMDFF 0
    variable REMDFFNumReplicas 6
    variable REMDFFAutoBlur 1
    variable xMDFF 0
    variable BFS 0
    variable BSharp "none"
    variable BSharpFactor "0"
    variable Mask 0
    variable AverageMap 0
    variable MaskRes 5
    variable MaskCutoff 5 
    variable RefSteps 20000
    variable PBCorGBIS ""
    variable IMD 0
    variable IMDPort 2000
    variable IMDFreq 1
    variable IMDWait 0
    variable IMDIgnore 0
    variable IMDProcs ""
    variable IMDKeep 0
    variable CrossCorrelation 0
    variable CCSel "protein and noh"
    variable CCRes ""
    variable UseCCThresh 0
    variable CCThresh ""
    variable GridforceLite 0
    variable QwikMDLogFile ""

         
    variable ParameterList [list [file join $env(CHARMMPARDIR) par_all36_prot.prm]\
    [file join $env(CHARMMPARDIR) par_all36_lipid.prm] \
  [file join $env(CHARMMPARDIR) par_all36_na.prm] [file join $env(CHARMMPARDIR) par_all36_carb.prm] \
  [file join $env(CHARMMPARDIR) par_all36_cgenff.prm] [file join $env(CHARMMPARDIR) toppar_all36_carb_glycopeptide.str]\
  [file join $env(CHARMMPARDIR) toppar_water_ions_namd.str]]
  }
  
  #MapTools Settings. separate because we don't really need to save these in our settings file.
  namespace eval maptoolssettings {
    variable MapToolsMolID ""  
    variable MapToolsVolID 0  
    variable MapToolsVol2ID 0  
    variable MapToolsPlotX ""  
    variable HistPlotBins 30
    variable HistLog 1 
    variable MapToolsXsize ""
    variable MapToolsYsize ""
    variable MapToolsZsize ""
    variable MapToolsMin ""
    variable MapToolsMax ""
    variable MapToolsCOM ""
    variable MapToolsOrigin ""
    variable TrimX1 0
    variable TrimX2 0
    variable TrimY1 0
    variable TrimY2 0
    variable TrimZ1 0
    variable TrimZ2 0
    variable CropX1 0
    variable CropX2 0
    variable CropY1 0
    variable CropY2 0
    variable CropZ1 0
    variable CropZ2 0
    variable MapToolsStructMolID ""
    variable FitSel "all"
    variable FitRes ""
    variable MaskSel "all"
    variable MaskRes ""
    variable MaskCutoff ""
    variable SimSel "noh"
    variable SimRes ""
    variable CCSel "noh"
    variable CCRes ""
    variable ClampMin 0
    variable ClampMax 0
    variable SmultAmt 1
    variable SAddAmt 0
    variable RangeMin 0
    variable RangeMax 0
    variable BinmaskAmt 0
    variable SmoothAmt 0
    variable PotAmt 0
    variable SegAmt 5
    variable MapToolsMol2ID ""
    variable UnionOrIntersection ""
    variable InterpolationOrNo ""

  }
  
  namespace eval servers {} ;# create the "servers" namespace
  
}

proc MDFFGUI::gui::mdffgui {} {
  variable w 
  variable wad 
  variable xwad 
  variable nullMolString "none"
  variable nullVolString 0
  variable CurrentPDBFile
  variable CurrentPSFFile
  variable DensityID ""
  variable CButtonChi
  variable MapID
  variable xMDFFMapIDs ""
  variable CButtonCis
  variable CCResInput
  variable CCSelInput
  variable MolMenuText  
  variable CurrentDXFile
  variable Extrabonds 
  variable CurrentRefs
  variable NumFrames 0
  variable DensityView 
  variable xMDFFDensityView
  variable IMDProcInput   
  variable CurrentCrystPDB
  variable CCMolInput
  variable xMDFFSelInput
  
  variable CurrentDirLabel
  
  variable env
  
  variable xMDFFDensityViewAdd
  variable CCStep
  variable CCThreshInput
  variable UseCCThresh

  variable IMDSubmit
  variable IMDConnect
  variable IMDPause
  variable IMDFinish
  variable IMDStatusLabel
  variable IMDCCPlot
  variable Notebook
  variable justChanged
  variable GridpdbGen
  variable FixedpdbGen
  variable MDFFSetup
  variable GenerateXBonds
  
  variable CrystOption
  variable CurrentRefsEntry
  variable RefStepsInput
  variable CurrentCrystEntry
  variable GetRefs
  variable GetCryst
  variable CButtonMask
  variable CButtonBFS
  variable MaskResInput
  variable MaskCutoffInput
  
  variable REMDFFNumReplicasInput
  variable CButtonREMDFFAutoBlur
  
  variable IMDInpOrList
  variable IMDHost
  variable namdbin
  variable CurrentServer 
  variable ServerList
  variable IMDHostInput
  variable IMDServer 
  variable CButtonIMDWait
  variable CButtonIMDIgnore
  variable IMDPortInput
  variable IMDFreqInput
  variable IMDKeepInput
  variable HLF
  variable ShowBTNPadX
  

  #variable HistPlot
#  # If already initialized, just deiconify and return
  if { [winfo exists .mdffgui] } {
    wm deiconify $w
    raise $w
    return
  }
  

  set downPoint \u25BC
  set rightPoint \u25B6
  
  variable CCx
  variable CCy
  
  set Extrabonds {}


  set CrystOption ""

  set IMDHost "localhost"
  set namdbin "namd2"

  #set MolID ""
  set MapID ""  

  set CCx ""
  set CCy ""
  set CCStep ""
  set justChanged 0
  
  set ShowBTNPadX 0
 
  
  # setup the theme depending on what is available
  set themeList [ttk::style theme names]
  if { [lsearch -exact $themeList "aqua"] != -1 } {
      ttk::style theme use aqua
      set ShowBTNPadX 18
  } elseif { [lsearch -exact $themeList "clam"] != -1 } {
      ttk::style theme use clam
  } elseif { [lsearch -exact $themeList "classic"] != -1 } {
      ttk::style theme use classic
  } else {
      ttk::style theme use default
  }
  MDFFGUI::gui::add_server "default localhost" {
    jobtype  local
    namdbin  {namd2 +p%d}
    maxprocs 1
    numprocs 1
    timeout  20

  }
  trace add variable MDFFGUI::settings::CrossCorrelation write MDFFGUI::gui::cc_toggle
  trace add variable MDFFGUI::settings::xMDFF write MDFFGUI::gui::xmdff_toggle
  trace add variable MDFFGUI::settings::Mask write MDFFGUI::gui::mask_toggle
  trace add variable MDFFGUI::settings::IMD write MDFFGUI::gui::imd_toggle
  trace add variable MDFFGUI::settings::UseCCThresh write MDFFGUI::gui::thresh_toggle



 # set CurrentServer ""
  #source [file join $env(MDFFDIR) settings.tcl]
  set ServerList [list]
  foreach s [lsort [info vars ::MDFFGUI::servers::*]] {
    lappend ServerList [namespace tail $s]
  }
  set CurrentServer [lindex $ServerList 0]

  set w [toplevel ".mdffgui"]
  wm title $w "MDFF GUI"
  #wm resizable $w 1 1
  grid columnconfigure $w 0 -weight 1
  grid rowconfigure $w 0 -weight 1
  
  set HLF [ttk::frame $w.hlf]
  grid $w.hlf -column 0 -row 0 -sticky nsew
  # allow hlf to resize with window
  grid columnconfigure $w.hlf 0 -weight 1
  grid rowconfigure $w.hlf 1 -weight 1
 
  set FileMenu [ttk::menubutton $w.hlf.filemenu -text "File" -menu $w.hlf.filemenu.save -width 5]
  set FileMenuSave [menu $w.hlf.filemenu.save -tearoff no] 
  $FileMenuSave add command -label "Save Settings..." -command "MDFFGUI::gui::save_settings $w"
  $FileMenuSave add command -label "Load Settings..." -command "MDFFGUI::gui::load_settings $w"
   
  
  set Notebook [ttk::notebook $w.hlf.n]
  set NBTab1 [ttk::frame $w.hlf.n.f1]; # first page, which would get widgets gridded into it 
  $Notebook add $NBTab1 -text "MDFF Setup"
  ttk::notebook::enableTraversal $Notebook
  
  grid columnconfigure $w.hlf.n 0 -weight 1
  grid rowconfigure $w.hlf.n 0 -weight 1
  
  grid columnconfigure $w.hlf.n.f1 0 -weight 1
  grid rowconfigure $w.hlf.n.f1 0 -weight 1
  #MDFF Files
  set MFrame [ttk::frame $w.hlf.n.f1.main]
  grid columnconfigure $MFrame 0 -weight 1
  grid rowconfigure $MFrame 6 -weight 1

  set FileFrame [ttk::labelframe $w.hlf.n.f1.main.fileframe -labelanchor nw]
  grid columnconfigure $FileFrame 1 -weight 1
  #grid rowconfigure $FileFrame 0 -weight 1
  #set BFrame [frame $w.hlf.n.f1.main.bframe -relief raised -bd 2]
  set CurrentPDBFile [ttk::label $w.hlf.n.f1.main.fileframe.pdblabel -text "PDB File:"]
  set CurrentPSFFile [ttk::label $w.hlf.n.f1.main.fileframe.psflabel -text "PSF File:"]
  set CurrentPDBEntry [ttk::entry $w.hlf.n.f1.main.fileframe.pdbentry -textvariable MDFFGUI::settings::CurrentPDBPath -width 40]
  set CurrentPSFEntry [ttk::entry $w.hlf.n.f1.main.fileframe.psfentry -textvariable MDFFGUI::settings::CurrentPSFPath -width 40]
  set SelectPDB [ttk::button $w.hlf.n.f1.main.fileframe.button5 -text "Browse" -command {MDFFGUI::gui::get_pdb} ]
  set SelectPSF [ttk::button $w.hlf.n.f1.main.fileframe.button6 -text "Browse" -command {MDFFGUI::gui::get_psf} ]
  set LoadFileSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.seperator  -orient horizontal]   
  set LoadFiles [ttk::button $w.hlf.n.f1.main.fileframe.button7 -text "Load PSF/PDB" -command {MDFFGUI::gui::load_struct}]

  set CurrentMol [ttk::label $w.hlf.n.f1.main.fileframe.mollabel -text "Mol ID (PDB and PSF):"]
::TKTOOLTIP::balloon $CurrentMol [MDFFGUI::gui::currentmol]

  #set CurrentMolEntry [ttk::entry $w.hlf.n.f1.main.fileframe.molentry -textvariable MDFFGUI::settings::MolID -width 10]

  set CurrentMolMenuButton [ttk::menubutton $w.hlf.n.f1.main.fileframe.molmenubutton -textvar MDFFGUI::gui::MolMenuText -menu $w.hlf.n.f1.main.fileframe.molmenubutton.molmenu]
  set CurrentMolMenu [menu $w.hlf.n.f1.main.fileframe.molmenubutton.molmenu -tearoff no]
  
  #set_mol_text 
  fill_mol_menu $CurrentMolMenu
  trace add variable ::vmd_initialize_structure write "MDFFGUI::gui::fill_mol_menu $CurrentMolMenu"
  trace add variable MDFFGUI::settings::MolID write MDFFGUI::gui::set_mol_text
  set MDFFGUI::settings::MolID $nullMolString

  
  set CurrentDirLabel [ttk::label $w.hlf.n.f1.main.fileframe.dirlabel -text "Working Directory:"]
  ::TKTOOLTIP::balloon $CurrentDirLabel [MDFFGUI::gui::workingdir]

  set CurrentDirEntry [ttk::entry $w.hlf.n.f1.main.fileframe.direntry -textvariable MDFFGUI::settings::CurrentDir -width 40]
  set SelectDir [ttk::button $w.hlf.n.f1.main.fileframe.selectdir -text "Browse" -command "MDFFGUI::gui::get_dir $CurrentDirEntry" ] 
  set CurrentDirSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.cdseperator  -orient horizontal]   
  
  set RestraintFrame [ttk::frame $w.hlf.n.f1.main.fileframe.restraints]
  set RestraintLabel [ttk::label $w.hlf.n.f1.main.fileframe.restraints.label -text "Restraints: "]
  ::TKTOOLTIP::balloon $RestraintLabel [MDFFGUI::gui::restraints]
  #grid columnconfigure $FileFrame 1 -weight 1
  set CButtonSecondary [ttk::checkbutton $w.hlf.n.f1.main.fileframe.restraints.ss -text "Secondary Structure Restraints" -variable MDFFGUI::settings::SSRestraints]
  ::TKTOOLTIP::balloon $CButtonSecondary [MDFFGUI::gui::ssrestraintstip]
  set CButtonChirality [ttk::checkbutton $w.hlf.n.f1.main.fileframe.restraints.chi -text "Chirality Restraints" -variable MDFFGUI::settings::ChiralityRestraints]
  ::TKTOOLTIP::balloon $CButtonChirality [MDFFGUI::gui::chiralityrestraints]
  set CButtonCispeptide [ttk::checkbutton $w.hlf.n.f1.main.fileframe.restraints.cis -text "Cispeptide Restraints" -variable MDFFGUI::settings::CispeptideRestraints]
  ::TKTOOLTIP::balloon $CButtonCispeptide [MDFFGUI::gui::cispeptiderestraints]
  set GenerateXBonds [ttk::button $w.hlf.n.f1.main.fileframe.xbonds -text "Generate Restraints" -command {MDFFGUI::gui::generate_xbonds} -state disabled]
  set XBondsSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.xbondsseperator  -orient horizontal]   
  
#  set GridpdbLabel [ttk::label $w.hlf.n.f1.main.fileframe.label2 -text "Gridpdb selection text:"]
#  set GridpdbInput [ttk::entry $w.hlf.n.f1.main.fileframe.entry2 -textvariable MDFFGUI::settings::GridPDBSelText -width 40]
  #set GridpdbGen [ttk::button $w.hlf.n.f1.main.fileframe.gridpdb -text "Generate GridPDB" -command {MDFFGUI::gui::make_gridpdb} -state disabled]
  set GridpdbSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.gridpdbseperator  -orient horizontal]   

  set FixedpdbLabel [ttk::label $w.hlf.n.f1.main.fileframe.fixedlabel -text "Fixed PDB selection text:"]
  ::TKTOOLTIP::balloon $FixedpdbLabel [MDFFGUI::gui::fixedpdbtext]  
  set FixedpdbInput [ttk::entry $w.hlf.n.f1.main.fileframe.fixedentry -textvariable MDFFGUI::settings::FixedPDBSelText -width 40]
  set FixedColumnLabel [ttk::label $w.hlf.n.f1.main.fileframe.fixedclabel -text "Fixed PDB Column:"]
  ::TKTOOLTIP::balloon $FixedColumnLabel [MDFFGUI::gui::fixedpdbcolumn]  
  set FixedColumnInput [ttk::entry $w.hlf.n.f1.main.fileframe.fixedcentry -textvariable MDFFGUI::settings::FixedColumn -width 40]
  set FixedpdbGen [ttk::button $w.hlf.n.f1.main.fileframe.fixedpdb -text "Generate Fixed PDB" -command {MDFFGUI::gui::make_fixedpdb} -state disabled]
  set FixedpdbSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.fixedpdbseperator  -orient horizontal]   
 
  set ParamListLabel [ttk::label $w.hlf.n.f1.main.fileframe.paramlistlabel -text "Parameter Files: "]
  ::TKTOOLTIP::balloon $ParamListLabel [MDFFGUI::gui::parameters]  
  set ParamListScroll [ttk::scrollbar $w.hlf.n.f1.main.fileframe.paramlistscroll -command "$w.hlf.n.f1.main.fileframe.paramlistbox yview" -orient vertical]
  set ParamListBox [listbox $w.hlf.n.f1.main.fileframe.paramlistbox -height 5 -yscrollcommand "$ParamListScroll set" -listvariable MDFFGUI::settings::ParameterList]
  set ParamListAdd [ttk::button $w.hlf.n.f1.main.fileframe.paramlistadd -text "Add" -command "MDFFGUI::gui::add_paramfile $ParamListBox"]
  set ParamListRemove [ttk::button $w.hlf.n.f1.main.fileframe.paramlistremove -text "Remove" -command {
   foreach ind [.mdffgui.hlf.n.f1.main.fileframe.paramlistbox curselection] {set MDFFGUI::settings::ParameterList [lreplace $MDFFGUI::settings::ParameterList $ind $ind]} 
  }]
  
  set CurrentDXFile [ttk::label $w.hlf.n.f1.main.fileframe.dxlabel -text "Density Map:"]
  ::TKTOOLTIP::balloon $CurrentDXFile [MDFFGUI::gui::currentdensitymap]
  
  set CurrentDXEntry [ttk::entry $w.hlf.n.f1.main.fileframe.dxentry -textvariable MDFFGUI::settings::CurrentDXPath -width 40]
  set SelectDX [ttk::button $w.hlf.n.f1.main.fileframe.dx -text "Browse" -command "MDFFGUI::gui::get_density $CurrentDXEntry" ]
  #set MakeDX [ttk::button $w.hlf.n.f1.main.fileframe.griddx -text "Generate MDFF Potential" -command {MDFFGUI::gui::make_griddx} ]
  set DXSeparator [ttk::separator $w.hlf.n.f1.main.fileframe.dxseperator  -orient horizontal]   

  set DensityViewScroll [ttk::scrollbar $w.hlf.n.f1.main.fileframe.denviewscroll -command "$w.hlf.n.f1.main.fileframe.denview yview" -orient vertical]
  set DensityViewScrollX [ttk::scrollbar $w.hlf.n.f1.main.fileframe.denviewscrollx -command "$w.hlf.n.f1.main.fileframe.denview xview" -orient horizontal]
  set DensityView [ttk::treeview $w.hlf.n.f1.main.fileframe.denview -selectmode browse -yscrollcommand "$DensityViewScroll set" -xscrollcommand "$DensityViewScrollX set" ]

  $DensityView configure -columns {file selection gscale threshold} -show {headings} -height 3
  $DensityView heading file -text "Density Filename"
  $DensityView heading selection -text "gridpdb selection"
  $DensityView heading gscale -text "gscale"
  $DensityView heading threshold -text "threshold"
  $DensityView column file -width 200 -stretch 0 -anchor w
  $DensityView column selection -width 75 -stretch 1 -anchor w
  $DensityView column gscale -width 50 -stretch 0 -anchor e
  $DensityView column threshold -width 50 -stretch 0 -anchor e
  set DensityViewAdd [ttk::button $w.hlf.n.f1.main.fileframe.denviewadd -text "Add" -command {MDFFGUI::gui::add_density}]
  
  set DensityViewRemove [ttk::button $w.hlf.n.f1.main.fileframe.denviewremove -text "Remove" -command {MDFFGUI::gui::remove_densityview  $MDFFGUI::gui::DensityView [$MDFFGUI::gui::DensityView focus] }]

  
  set NameLabel [ttk::label $w.hlf.n.f1.main.fileframe.simnamelabel -text "Simulation Output Name:"]
  ::TKTOOLTIP::balloon $NameLabel [MDFFGUI::gui::simoutputname]  
  set NameInput [ttk::entry $w.hlf.n.f1.main.fileframe.entry3 -textvariable MDFFGUI::settings::SimulationName -width 40]
  
  set StepLabel [ttk::label $w.hlf.n.f1.main.fileframe.simsteplabel -text "Simulation Output Step:"]
  ::TKTOOLTIP::balloon $StepLabel [MDFFGUI::gui::simstep]  
  set StepInput [ttk::entry $w.hlf.n.f1.main.fileframe.stepentry -textvariable MDFFGUI::settings::SimulationStep -width 40]
  
  set SetupSeparator [ttk::separator $w.hlf.n.f1.main.setupseperator  -orient horizontal]
#there is code in place to allow this to start disabled and only open once a psf/pdb and map are
  #specified. however, this doesn't account for xMDFF where the user just puts in a name and
  #doesn't click browse, so this button is normal now by default
  set MDFFSetup [ttk::button $w.hlf.n.f1.main.setup -text "Generate NAMD files" -command {MDFFGUI::gui::mdff_setup} -state normal]
 
  set ShowMDFF [ttk::label $w.hlf.n.f1.main.showmdff -text "$rightPoint MDFF Files..." -anchor w]
  ::TKTOOLTIP::balloon $ShowMDFF [MDFFGUI::gui::mdfffiles]  
  set HideMDFF [ttk::label $w.hlf.n.f1.main.fileframe.hidemdff -text "$downPoint MDFF Files" -anchor w]
  ::TKTOOLTIP::balloon $HideMDFF [MDFFGUI::gui::mdfffiles]  
  $FileFrame configure -labelwidget $HideMDFF
  bind $HideMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.fileframe
      grid .mdffgui.hlf.n.f1.main.showmdff
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showmdff 
      grid .mdffgui.hlf.n.f1.main.fileframe -row 0 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.fileframe 1 -weight 1
      MDFFGUI::gui::resizeToActiveTab
  }
  grid $HLF -row 0 -column 0 -sticky nswe
  grid $FileMenu -row 0 -column 0 -sticky nsw
  grid $Notebook -row 1 -column 0 -sticky nsew
  grid $MFrame -row 0 -column 0 -sticky nswe
  grid $ShowMDFF -row 0 -column 0 -sticky nswe -pady 5 -padx $ShowBTNPadX
 # grid $HideMDFF -row 0 -column 0 -sticky nsew
  #grid $BFrame -row 1 -column 0 -columnspan 5 -sticky nsew
 
  grid $CurrentDirLabel -row 1 -column 0 -columnspan 2 -sticky nsew
  grid $CurrentDirEntry -row 1 -column 1 -columnspan 2 -sticky nsew
  grid $SelectDir -row 1 -column 3 -sticky nsew
  #grid $CurrentDirSeparator -row 2 -column 0 -columnspan 5 -sticky ew
 
  grid $CurrentMol -row 3 -column 0 -sticky nswe
 # grid $CurrentMolEntry -row 3 -column 1 -sticky nsew
  grid $CurrentMolMenuButton -row 3 -column 1 -sticky nsw
#  grid $CurrentPSFFile -row 3 -column 0 -sticky nswe
#  grid $CurrentPSFEntry -row 3 -column 1 -columnspan 2 -sticky nsew
#  grid $CurrentPDBFile -row 4 -column 0 -sticky nswe
#  grid $CurrentPDBEntry -row 4 -column 1 -columnspan 2 -sticky nsew
#  grid $SelectPSF -row 3 -column 3 -sticky nswe
#  grid $SelectPDB -row 4 -column 3 -sticky nswe
#  #grid $LoadFileSeparator -row 6 -column 0 -columnspan 5 -sticky ew
#  grid $LoadFiles -row 5 -column 0 -columnspan 5 -sticky nswe -padx 8 -pady 8

  
  grid $CurrentDXFile -row 7 -column 0 -sticky nsew
  #grid $CurrentDXEntry -row 7 -column 1 -columnspan 2 -sticky nsew 
  #grid $SelectDX -row 7 -column 3 -sticky nsew

  grid $DensityView -row 7 -column 1 -sticky nsew -rowspan 2
  grid $DensityViewScrollX -row 9 -column 1 -sticky nsew
  grid $DensityViewAdd -row 7 -column 3 -sticky nsew
  grid $DensityViewRemove -row 8 -column 3 -sticky nsew
  grid $DensityViewScroll -row 7 -column 2 -sticky nsew -rowspan 2

  #grid $MakeDX -row 8 -column 0 -sticky nsew -columnspan 5 -padx 8 -pady 8
  #grid $DXSeparator -row 9 -column 0 -columnspan 5 -sticky ew

#  grid $GridpdbLabel -row 10 -column 0 -sticky nsew
#  grid $GridpdbInput -row 10 -column 1 -columnspan 2 -sticky nsew
 # grid $GridpdbGen -row 11 -column 0 -sticky nsew -columnspan 5 -padx 8 -pady 8
 # grid $GridpdbSeparator -row 12 -column 0 -columnspan 5 -sticky ew
  

  grid $FixedpdbLabel -row 13 -column 0 -sticky nsew
  grid $FixedpdbInput -row 13 -column 1 -columnspan 2 -sticky nsew
  grid $FixedColumnLabel -row 14 -column 0 -sticky nsew
  grid $FixedColumnInput -row 14 -column 1 -columnspan 2 -sticky nsew
  #grid $FixedpdbGen -row 15 -column 0 -sticky nsew -columnspan 5 -padx 8 -pady 8
 # grid $FixedpdbSeparator -row 16 -column 0 -columnspan 5 -sticky ew


  grid $NameLabel -row 17 -column 0 -sticky nsew
  grid $NameInput -row 17 -column 1 -columnspan 2 -sticky nsew
  
  grid $StepLabel -row 18 -column 0 -sticky nsew
  grid $StepInput -row 18 -column 1 -columnspan 2 -sticky nsew
  #grid $XBondsSeparator -row 18 -column 0 -columnspan 5 -sticky ew
 
  grid $RestraintFrame -row 19 -column 0 -sticky nsw -columnspan 4
  grid $RestraintLabel -row 0 -column 0 -sticky nsew
  grid $CButtonSecondary -row 0 -column 1 -sticky nsew
  grid $CButtonChirality -row 0 -column 2 -sticky nsew
  grid $CButtonCispeptide -row 0 -column 3 -sticky nsew
  #grid $GenerateXBonds -row 18 -column 0 -sticky nsew -columnspan 5 -padx 8 -pady 8
  
  grid $ParamListLabel -row 20 -column 0 -sticky nsew
  grid $ParamListBox -row 21 -column 0 -sticky nsew -columnspan 2 -rowspan 2
  grid $ParamListAdd -row 21 -column 3 -sticky nsew
  grid $ParamListRemove -row 22 -column 3 -sticky nsew
  grid $ParamListScroll -row 21 -column 2 -sticky nsew -rowspan 2
  
  grid $SetupSeparator -row 4 -column 0 -columnspan 5 -sticky ew
  grid $MDFFSetup -row 5 -column 0 -sticky nsew -padx 8 -pady 8
 
 
  #Parameters

 # set NBTab2 [ttk::frame $w.hlf.n.f2]; # second page
 # $Notebook add $NBTab2 -text "Parameters" 
  
  set ParamFrame [ttk::labelframe $w.hlf.n.f1.main.paramframe -labelanchor nw]
  grid columnconfigure $ParamFrame 1 -weight 1

  #grid columnconfigure $w.n.f2.main2 0 -weight 1
  set TempLabel [ttk::label $w.hlf.n.f1.main.paramframe.templabel -text "Temperature (K):"]
  set TempInput [ttk::entry $w.hlf.n.f1.main.paramframe.tempentry -textvariable MDFFGUI::settings::Temperature]
  ::TKTOOLTIP::balloon $TempLabel [MDFFGUI::gui::itemp]  
  
  set FTempLabel [ttk::label $w.hlf.n.f1.main.paramframe.ftemplabel -text "Final Temperature (K):"]
  set FTempInput [ttk::entry $w.hlf.n.f1.main.paramframe.ftempentry -textvariable MDFFGUI::settings::FTemperature]
  ::TKTOOLTIP::balloon $FTempLabel [MDFFGUI::gui::ftemp]  
  
  set MinStepLabel [ttk::label $w.hlf.n.f1.main.paramframe.minsteplabel -text "Minimization Steps:"]
  set MinStepInput [ttk::entry $w.hlf.n.f1.main.paramframe.minstepentry -textvariable MDFFGUI::settings::Minsteps]
  ::TKTOOLTIP::balloon $MinStepLabel [MDFFGUI::gui::minsteps]  
  
  set NumStepLabel [ttk::label $w.hlf.n.f1.main.paramframe.numsteplabel -text "Time steps:"]
  set NumStepInput [ttk::entry $w.hlf.n.f1.main.paramframe.numstepentry -textvariable MDFFGUI::settings::Numsteps]
  ::TKTOOLTIP::balloon $NumStepLabel [MDFFGUI::gui::numsteps]  
 
  #set GScaleLabel [ttk::label $w.hlf.n.f1.main.paramframe.gscalelabel -text "Grid Scaling Factor:"]
  #set GScaleInput [ttk::entry $w.hlf.n.f1.main.paramframe.gscaleentry -textvariable MDFFGUI::settings::GScale]

  set CButtonLite [ttk::checkbutton $w.hlf.n.f1.main.paramframe.gridforcelite -text "gridforcelite" -variable MDFFGUI::settings::GridforceLite]
  ::TKTOOLTIP::balloon $CButtonLite [MDFFGUI::gui::gridforcelite]  

  set CButtonGridOff [ttk::checkbutton $w.hlf.n.f1.main.paramframe.gridoff -text "Turn off grid forces" -variable MDFFGUI::settings::GridOff]
  ::TKTOOLTIP::balloon $CButtonGridOff [MDFFGUI::gui::gridforceoff]  
  
  set EnvironFrame [ttk::frame $w.hlf.n.f1.main.paramframe.envframe]
  set EnvironLabel [ttk::label $w.hlf.n.f1.main.paramframe.envframe.label -text "System Environment: "]
  set PBC [ttk::radiobutton $w.hlf.n.f1.main.paramframe.envframe.pbc -variable MDFFGUI::settings::PBCorGBIS -value "-pbc" -text "Periodic Boundary Conditions"]
  ::TKTOOLTIP::balloon $PBC [MDFFGUI::gui::pbc]  
  
  set GBIS [ttk::radiobutton $w.hlf.n.f1.main.paramframe.envframe.gbis -variable MDFFGUI::settings::PBCorGBIS -value "-gbis" -text "Implicit Solvent"]
  ::TKTOOLTIP::balloon $GBIS [MDFFGUI::gui::gbis]  
  
  set VAC [ttk::radiobutton $w.hlf.n.f1.main.paramframe.envframe.vac -variable MDFFGUI::settings::PBCorGBIS -value "" -text "Vacuum"]
  ::TKTOOLTIP::balloon $VAC [MDFFGUI::gui::vac]  
  
  set ShowParams [ttk::label $w.hlf.n.f1.main.showparams -text "$rightPoint Simulation Parameters..." -anchor w]
  ::TKTOOLTIP::balloon $ShowParams [MDFFGUI::gui::simparams]  
  set HideParams [ttk::label $w.hlf.n.f1.main.paramframe.hideparams -text "$downPoint Simulation Parameters" -anchor w]
  ::TKTOOLTIP::balloon $HideParams [MDFFGUI::gui::simparams]  
  $ParamFrame configure -labelwidget $HideParams
  bind $HideParams <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.paramframe
      grid .mdffgui.hlf.n.f1.main.showparams
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowParams <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showparams
      grid .mdffgui.hlf.n.f1.main.paramframe -row 1 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.paramframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }
   
  #pack $Notebook
  #grid $MFrameP -row 0 -column 0

  grid $ShowParams -row 1 -column 0 -sticky w -pady 5 -padx $ShowBTNPadX
  grid $TempLabel -row 0 -column 0 -sticky nw
  grid $TempInput -row 0 -column 1 -sticky nsew
  
  grid $FTempLabel -row 1 -column 0 -sticky nw
  grid $FTempInput -row 1 -column 1 -sticky nsew

  grid $MinStepLabel -row 2 -column 0 -sticky nw
  grid $MinStepInput -row 2 -column 1 -sticky nsew
  
  grid $NumStepLabel -row 3 -column 0 -sticky nw
  grid $NumStepInput -row 3 -column 1 -sticky nsew

  #grid $GScaleLabel -row 4 -column 0 -sticky nw
  #grid $GScaleInput -row 4 -column 1 -sticky nsew
  
  grid $CButtonLite -row 5 -column 0 -sticky nw
  grid $CButtonGridOff -row 6 -column 0 -sticky nw

  grid $EnvironFrame -row 7 -column 0 -sticky nsew -columnspan 4
  grid $EnvironLabel -row 0 -column 0 -sticky nsew
  grid $VAC -row 0 -column 1 -sticky nsew
  grid $PBC -row 0 -column 2 -sticky nsew
  grid $GBIS -row 0 -column 3 -sticky nsew

  #xMDFF 
  #set NBTab3 [ttk::frame $w.hlf.n.f3];
  #$Notebook add $NBTab3 -text "xMDFF" 


  set xMDFFFrame [ttk::labelframe $w.hlf.n.f1.main.xmdffframe ]
  grid columnconfigure $xMDFFFrame 1 -weight 1
  
  set CButtonxMDFF [ttk::checkbutton $w.hlf.n.f1.main.xmdffframe.xmdff -text "xMDFF" -variable MDFFGUI::settings::xMDFF -command {MDFFGUI::gui::xmdff_toggle}]
  ::TKTOOLTIP::balloon $CButtonxMDFF [MDFFGUI::gui::doxmdff]  
  
  set xMDFFDensityLabel [ttk::label $w.hlf.n.f1.main.xmdffframe.denlabel -text "xMDFF Density File:"] 
  ::TKTOOLTIP::balloon $xMDFFDensityLabel [MDFFGUI::gui::xmdffdensity]  
  set xMDFFDensityViewScroll [ttk::scrollbar $w.hlf.n.f1.main.xmdffframe.denviewscroll -command "$w.hlf.n.f1.main.xmdffframe.denview yview" -orient vertical]
  set xMDFFDensityViewScrollX [ttk::scrollbar $w.hlf.n.f1.main.xmdffframe.denviewscrollx -command "$w.hlf.n.f1.main.xmdffframe.denview xview" -orient horizontal]
  set xMDFFDensityView [ttk::treeview $w.hlf.n.f1.main.xmdffframe.denview -selectmode browse -yscrollcommand "$xMDFFDensityViewScroll set" -xscrollcommand "$xMDFFDensityViewScrollX set"]

  $xMDFFDensityView configure -columns {file selection gscale mask maskres maskcutoff average bsharp bsharpfactor} -show {headings} -height 3
  $xMDFFDensityView heading file -text "Reflections Filename"
  $xMDFFDensityView heading selection -text "gridpdb selection"
  $xMDFFDensityView heading gscale -text "gscale"
  $xMDFFDensityView heading mask -text "mask"
  $xMDFFDensityView heading maskres -text "mask res (A)"
  $xMDFFDensityView heading maskcutoff -text "mask cutoff (A)"
  $xMDFFDensityView heading average -text "Average Map"
  $xMDFFDensityView heading bsharp -text "Sharpen/Smooth"
  $xMDFFDensityView heading bsharpfactor -text "Sharpen/Smooth Factor"
  $xMDFFDensityView column file -width 200 -stretch 0 -anchor w
  $xMDFFDensityView column selection -width 150 -stretch 0 -anchor w
  $xMDFFDensityView column gscale -width 50 -stretch 0 -anchor w
  $xMDFFDensityView column mask -width 50 -stretch 0 -anchor w
  $xMDFFDensityView column maskres -width 90 -stretch 0 -anchor w
  $xMDFFDensityView column maskcutoff -width 110 -stretch 0 -anchor w
  $xMDFFDensityView column average -width 85 -stretch 0 -anchor w
  $xMDFFDensityView column bsharp -width 110 -stretch 0 -anchor w
  $xMDFFDensityView column bsharpfactor -width 150 -stretch 0 -anchor w
  set xMDFFDensityViewAdd [ttk::button $w.hlf.n.f1.main.xmdffframe.denviewadd -text "Add" -command {MDFFGUI::gui::add_density_xmdff} -state disabled]
  
  set xMDFFDensityViewRemove [ttk::button $w.hlf.n.f1.main.xmdffframe.denviewremove -text "Remove" -command {MDFFGUI::gui::remove_densityview $MDFFGUI::gui::xMDFFDensityView [$MDFFGUI::gui::xMDFFDensityView focus]}]
  
  set CurrentRefs [ttk::label $w.hlf.n.f1.main.xmdffframe.refs -text "Reflection Data File:"] 
  set CurrentRefsEntry [ttk::entry $w.hlf.n.f1.main.xmdffframe.refsentry -textvariable MDFFGUI::settings::CurrentRefsPath -width 40 -state disabled]
  set GetRefs [ttk::button $w.hlf.n.f1.main.xmdffframe.getrefs -text "Browse" -command {MDFFGUI::gui::get_refs} -state disabled ]
  
  set xMDFFSelLabel [ttk::label $w.hlf.n.f1.main.xmdffframe.xmdffsellabel -text "Refinement Selection:"]
  ::TKTOOLTIP::balloon $xMDFFSelLabel [MDFFGUI::gui::xmdffsel]  
  set xMDFFSelInput [ttk::entry $w.hlf.n.f1.main.xmdffframe.xmdffselentry -textvariable MDFFGUI::settings::xMDFFSel -width 40 -state disabled]
  
  set RefStepsLabel [ttk::label $w.hlf.n.f1.main.xmdffframe.refstepslabel -text "Refinement Steps:"]
  ::TKTOOLTIP::balloon $RefStepsLabel [MDFFGUI::gui::xmdffsteps]  
  set RefStepsInput [ttk::entry $w.hlf.n.f1.main.xmdffframe.refstepsentry -textvariable MDFFGUI::settings::RefSteps -width 40 -state disabled]
  
  set CurrentCrystPDB [ttk::label $w.hlf.n.f1.main.xmdffframe.crystpdb -text "Symmetry PDB (optional):"] 
  ::TKTOOLTIP::balloon $CurrentCrystPDB [MDFFGUI::gui::xmdffsymmetry]  
  set CurrentCrystEntry [ttk::entry $w.hlf.n.f1.main.xmdffframe.crystentry -textvariable MDFFGUI::settings::CurrentCrystPath -width 40 -state disabled]
  set GetCryst [ttk::button $w.hlf.n.f1.main.xmdffframe.getcrystpdb -text "Browse" -command "MDFFGUI::gui::get_cryst $CurrentCrystEntry" -state disabled ]
  
#  set CButtonMask [ttk::checkbutton $w.hlf.n.f1.main.xmdffframe.mask -text "Mask Density" -variable MDFFGUI::settings::Mask -state disabled -command {MDFFGUI::gui::mask_toggle}]
  
  set CButtonBFS [ttk::checkbutton $w.hlf.n.f1.main.xmdffframe.bfs -text "Calculate Beta Factors" -variable MDFFGUI::settings::BFS -state disabled]
  ::TKTOOLTIP::balloon $CButtonBFS [MDFFGUI::gui::xmdffbf]  

#  set MaskResLabel [ttk::label $w.hlf.n.f1.main.xmdffframe.maskreslabel -text "Mask Resolution (A):"]
#  set MaskResInput [ttk::entry $w.hlf.n.f1.main.xmdffframe.maskresentry -textvariable MDFFGUI::settings::MaskRes -width 40 -state disabled]
  
#  set MaskCutoffLabel [ttk::label $w.hlf.n.f1.main.xmdffframe.maskcutofflabel -text "Mask Cutoff (A):"]
#  set MaskCutoffInput [ttk::entry $w.hlf.n.f1.main.xmdffframe.maskcutoffentry -textvariable MDFFGUI::settings::MaskCutoff -width 40 -state disabled]
  
  set ShowxMDFF [ttk::label $w.hlf.n.f1.main.showxmdff -text "$rightPoint xMDFF..." -anchor w]
  ::TKTOOLTIP::balloon $ShowxMDFF [MDFFGUI::gui::xmdff]  
  set HidexMDFF [ttk::label $w.hlf.n.f1.main.xmdffframe.hidexmdff -text "$downPoint xMDFF" -anchor w]
  ::TKTOOLTIP::balloon $HidexMDFF [MDFFGUI::gui::xmdff]  
  $xMDFFFrame configure -labelwidget $HidexMDFF
  bind $HidexMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.xmdffframe
      grid .mdffgui.hlf.n.f1.main.showxmdff
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowxMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showxmdff
      grid .mdffgui.hlf.n.f1.main.xmdffframe -row 2 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.xmdffframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }
  
  #grid $xMDFFFrame -row 0 -column 0
  grid $ShowxMDFF -row 2 -column 0 -sticky nsew -pady 5 -padx $ShowBTNPadX
  grid $CButtonxMDFF -row 0 -column 0 -sticky nw
  
  grid $xMDFFDensityLabel -row 1 -column 0 -sticky nw
  grid $xMDFFDensityView -row 1 -column 1 -sticky nsew -rowspan 2
  grid $xMDFFDensityViewScrollX -row 3 -column 1 -sticky nsew
  grid $xMDFFDensityViewAdd -row 1 -column 3 -sticky nsew
  grid $xMDFFDensityViewRemove -row 2 -column 3 -sticky nsew
  grid $xMDFFDensityViewScroll -row 1 -column 2 -sticky nsew -rowspan 2
  
#  grid $CurrentRefs -row 1 -column 0 -sticky nw
#  grid $CurrentRefsEntry -row 1 -column 1 -sticky nsew
#  grid $GetRefs -row 1 -column 2 -sticky nsew

  grid $CurrentCrystPDB -row 4 -column 0 -sticky nw -columnspan 2
  grid $CurrentCrystEntry -row 4 -column 1 -sticky nsew -columnspan 2
  grid $GetCryst -row 4 -column 3 -sticky nsew
  
  grid $RefStepsLabel -row 5 -column 0 -sticky nw
  grid $RefStepsInput -row 5 -column 1 -sticky nsew
  
  grid $xMDFFSelLabel -row 6 -column 0 -sticky nw
  grid $xMDFFSelInput -row 6 -column 1 -sticky nsew

 # grid $CButtonMask -row 5 -column 0 -sticky nsew
 # grid $MaskResLabel -row 6 -column 0 -sticky nw
#  grid $MaskResInput -row 6 -column 1 -sticky nsew
 
#  grid $MaskCutoffLabel -row 7 -column 0 -sticky nw
#  grid $MaskCutoffInput -row 7 -column 1 -sticky nsew
  
  grid $CButtonBFS -row 8 -column 0 -sticky nsew

  
  #IMDFF
  
  set IMDFrame [ttk::labelframe $w.hlf.n.f1.main.imdframe]
  grid columnconfigure $IMDFrame 1 -weight 1
  
  set CButtonIMD [ttk::checkbutton $w.hlf.n.f1.main.imdframe.imdff -text "IMD" -variable MDFFGUI::settings::IMD -command {MDFFGUI::gui::imd_toggle}]
  ::TKTOOLTIP::balloon $CButtonIMD [MDFFGUI::gui::doimd]  
  
  #set IMDSelectionLabel [ttk::label $w.hlf.n.f1.main.fileframe.imdsellabel -text "IMD selection text:"]
  #set IMDSelectionInput [ttk::entry $w.hlf.n.f1.main.fileframe.imdselentry -textvariable MDFFGUI::settings::IMDSelText -width 40]

  set IMDPortLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdportlabel -text "IMD Port:"]
  ::TKTOOLTIP::balloon $IMDPortLabel [MDFFGUI::gui::imdport]  
  set IMDPortInput [ttk::entry $w.hlf.n.f1.main.imdframe.imdportentry -textvariable MDFFGUI::settings::IMDPort -state disabled]
  
  set IMDFreqLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdfreqlabel -text "IMD Frequency:"]
  ::TKTOOLTIP::balloon $IMDFreqLabel [MDFFGUI::gui::imdfreq]  
  set IMDFreqInput [ttk::entry $w.hlf.n.f1.main.imdframe.imdfreqentry -textvariable MDFFGUI::settings::IMDFreq -state disabled]
  
  set IMDKeepLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdkeeplabel -text "IMD Keep Frames:"]
  ::TKTOOLTIP::balloon $IMDKeepLabel [MDFFGUI::gui::imdkeep]  
  set IMDKeepInput [ttk::entry $w.hlf.n.f1.main.imdframe.imdkeepentry -textvariable MDFFGUI::settings::IMDKeep -state disabled]
  
  set CButtonIMDWait [ttk::checkbutton $w.hlf.n.f1.main.imdframe.imdffwait -text "Wait for IMD Connection" -variable MDFFGUI::settings::IMDWait -state disabled]
  ::TKTOOLTIP::balloon $CButtonIMDWait [MDFFGUI::gui::imdwait]  
  
  set CButtonIMDIgnore [ttk::checkbutton $w.hlf.n.f1.main.imdframe.imdffignore -text "Ignore IMD forces" -variable MDFFGUI::settings::IMDIgnore -state disabled]
  ::TKTOOLTIP::balloon $CButtonIMDIgnore [MDFFGUI::gui::imdignore]  


  #set IMDInputSwitch [radiobutton $w.hlf.n.f4.main.imdradioinput -variable IMDInpOrList -value "input" -text "Input IMD Hostname" -command {MDFFGUI:select_imd_host}]
  #set IMDHostLabel [label $w.hlf.n.f4.main.imdhostlabel -text "IMD Hostname"]
  #set IMDHostInput [entry $w.hlf.n.f4.main.imdhostentry -textvariable IMDHost -state disabled]
  #set IMDHostLabel [label $w.hlf.n.f4.main.imdhostlabel -text "OR"]
 # set IMDListSwitch [radiobutton $w.hlf.n.f4.main.imdradiolist -variable IMDInpOrList -value "list" -text "IMD Host from List" -command {MDFFGUI:select_imd_host}]
  set IMDHostLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdhostlabel -text "IMD Server:"]
  ::TKTOOLTIP::balloon $IMDHostLabel [MDFFGUI::gui::imdhost]  
  set IMDServer [ttk::combobox $w.hlf.n.f1.main.imdframe.imdserver -textvariable MDFFGUI::gui::CurrentServer -values $MDFFGUI::gui::ServerList -state disabled]

  set IMDProcLabel [ttk::label $w.hlf.n.f1.main.imdframe.imdproclabel -text "Processors:"]
  ::TKTOOLTIP::balloon $IMDProcLabel [MDFFGUI::gui::imdproc]  
  set IMDProcInput [ttk::entry $w.hlf.n.f1.main.imdframe.imdprocentry -textvariable MDFFGUI::settings::IMDProcs -state disabled]
 
   
  bind $IMDServer <<ComboboxSelected>> {MDFFGUI::gui::imd_server}
  
  set ShowIMD [ttk::label $w.hlf.n.f1.main.showimd -text "$rightPoint IMD Parameters..." -anchor w]
  ::TKTOOLTIP::balloon $ShowIMD [MDFFGUI::gui::imd]  
  set HideIMD [ttk::label $w.hlf.n.f1.main.imdframe.hideimd -text "$downPoint IMD Parameters" -anchor w]
  ::TKTOOLTIP::balloon $HideIMD [MDFFGUI::gui::imd]  
  $IMDFrame configure -labelwidget $HideIMD
  bind $HideIMD <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.imdframe
      grid .mdffgui.hlf.n.f1.main.showimd
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowIMD <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showimd
      grid .mdffgui.hlf.n.f1.main.imdframe -row 3 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.imdframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }

  #grid $MFrameI -row 0 -column 0

  grid $ShowIMD -row 3 -column 0 -sticky nsew -pady 5 -padx $ShowBTNPadX
    
  grid $CButtonIMD -row 0 -column 0 -sticky nw
  
  grid $CButtonIMDWait -row 1 -column 0 -sticky nw
  
  grid $CButtonIMDIgnore -row 2 -column 0 -sticky nw
  
 # grid $IMDSelectionLabel -row 3 -column 0 -sticky nsew
 # grid $IMDSelectionInput -row 3 -column 1 -columnspan 2 -sticky nsew

  grid $IMDPortLabel -row 4 -column 0 -sticky nw
  grid $IMDPortInput -row 4 -column 1 -sticky nsew

  grid $IMDFreqLabel -row 5 -column 0 -sticky nw
  grid $IMDFreqInput -row 5 -column 1 -sticky nsew
  
  grid $IMDKeepLabel -row 6 -column 0 -sticky nw
  grid $IMDKeepInput -row 6 -column 1 -sticky nsew

  grid $IMDHostLabel -row 7 -column 0 -sticky nw
#  grid $IMDListSwitch -row 5 -column 3
  grid $IMDServer -row 7 -column 1 -sticky nsew
 
  grid $IMDProcLabel -row 8 -column 0 -sticky nw
  grid $IMDProcInput -row 8 -column 1 -sticky nsew
  
  #REMDFF
  
  set REMDFFFrame [ttk::labelframe $w.hlf.n.f1.main.remdffframe]
  grid columnconfigure $REMDFFFrame 1 -weight 1
 
  set REMDFFLabel [ttk::label $w.hlf.n.f1.main.remdffframe.remdfflabel -text "Test Label"]
  set CButtonREMDFF [ttk::checkbutton $w.hlf.n.f1.main.remdffframe.remdff -text "Resolution Exchange MDFF (BETA)" -variable MDFFGUI::settings::REMDFF -command {MDFFGUI::gui::remdff_toggle}]

  set CButtonREMDFFAutoBlur [ttk::checkbutton $w.hlf.n.f1.main.remdffframe.remdffautoblur -text "Automatically Generate Replica Potentials" -variable MDFFGUI::settings::REMDFFAutoBlur -state disabled]

  set REMDFFReplicasLabel [ttk::label $w.hlf.n.f1.main.remdffframe.remdffreplicaslabel -text "Number of Replicas:"]
  set REMDFFNumReplicasInput [ttk::entry $w.hlf.n.f1.main.remdffframe.remdffreplicasentry -textvariable MDFFGUI::settings::REMDFFNumReplicas -state disabled]

  set ShowREMDFF [ttk::label $w.hlf.n.f1.main.showremdff -text "$rightPoint Resolution Exchange Parameters (BETA) ..." -anchor w]
  set HideREMDFF [ttk::label $w.hlf.n.f1.main.remdffframe.hideremdff -text "$downPoint Resolution Exchange Parameters (BETA)" -anchor w]
  
  $REMDFFFrame configure -labelwidget $HideREMDFF
  bind $HideREMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.remdffframe
      grid .mdffgui.hlf.n.f1.main.showremdff
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowREMDFF <Button-1> {
      grid remove .mdffgui.hlf.n.f1.main.showremdff
      grid .mdffgui.hlf.n.f1.main.remdffframe -row 4 -column 0 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f1.main.remdffframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }
  
  grid $ShowREMDFF -row 4 -column 0 -sticky nsew -pady 5 -padx $ShowBTNPadX
 
  grid $CButtonREMDFF -row 0 -column 0 -sticky nw 
  grid $CButtonREMDFFAutoBlur -row 1 -column 0 -sticky nw 
  
  grid $REMDFFReplicasLabel -row 2 -column 0 -sticky nw
  grid $REMDFFNumReplicasInput -row 2 -column 1 -sticky nw
  
  #IMD Connect
  set NBTab4 [ttk::frame $w.hlf.n.f4];
  $Notebook add $NBTab4 -text "IMDFF Connect" 
  grid columnconfigure $NBTab4 0 -weight 1
  grid rowconfigure $NBTab4 0 -weight 1
  set IMDConnectFrame [ttk::frame $w.hlf.n.f4.main]
  grid columnconfigure $IMDConnectFrame 0 -weight 1
  set IMDButtonFrame [ttk::frame $w.hlf.n.f4.main.buttons]
  grid columnconfigure $IMDButtonFrame {0 1 2 3} -weight 1
  
  set IMDSubmit [ttk::button $w.hlf.n.f4.main.buttons.imdsubmit -text "Submit" -command {MDFFGUI::gui::submit_job} -state disabled ]
  
  set IMDConnect [ttk::button $w.hlf.n.f4.main.buttons.imdconnect -text "Connect" -command {MDFFGUI::gui::imd_connect} -state disabled ]
  
  set IMDPause [ttk::button $w.hlf.n.f4.main.buttons.imdpause -text "Pause" -command {MDFFGUI::gui::imd_pause} -state disabled ]
  
  set IMDFinish [ttk::button $w.hlf.n.f4.main.buttons.imdfinish -text "Finish" -command {MDFFGUI::gui::imd_kill} -state disabled ]
    
  set IMDStatusLabel [ttk::label $w.hlf.n.f4.main.imdstatuslabel -text "IMD Status:"]
  
  set CCFrame [ttk::labelframe $w.hlf.n.f4.main.ccframe -labelanchor nw]
  grid columnconfigure $CCFrame 1 -weight 1
  grid rowconfigure $CCFrame 3 -weight 1
  set CButtonCC [ttk::checkbutton $w.hlf.n.f4.main.ccframe.cc -text "Calculate Real-Time Cross Correlation" -variable MDFFGUI::settings::CrossCorrelation -command {MDFFGUI::gui::cc_toggle}]
  
  set CCMolLabel [ttk::label $w.hlf.n.f4.main.ccframe.ccmollabel -text "Experimental Density (Mol ID):"]
  set CCMolInput [ttk::entry $w.hlf.n.f4.main.ccframe.ccmolentry -textvariable MDFFGUI::gui::DensityID -width 40 -state disabled]
  
  set CCSelLabel [ttk::label $w.hlf.n.f4.main.ccframe.ccsellabel -text "Selection:"]
  set CCSelInput [ttk::entry $w.hlf.n.f4.main.ccframe.ccselentry -textvariable MDFFGUI::settings::CCSel -width 40 -state disabled]
  
  set CCResLabel [ttk::label $w.hlf.n.f4.main.ccframe.ccreslabel -text "Map Resolution:"]
  set CCResInput [ttk::entry $w.hlf.n.f4.main.ccframe.ccresentry -textvariable MDFFGUI::settings::CCRes -width 40 -state disabled]
  
  set CButtonCCThresh [ttk::checkbutton $w.hlf.n.f4.main.ccframe.ccthresh -text "Use Threshold" -variable MDFFGUI::settings::UseCCThresh -command {MDFFGUI::gui::thresh_toggle}]
  set CCThreshLabel [ttk::label $w.hlf.n.f4.main.ccframe.ccthreshlabel -text "Threshold:"]
  set CCThreshInput [ttk::entry $w.hlf.n.f4.main.ccframe.ccthreshentry -textvariable MDFFGUI::settings::CCThresh -width 40 -state disabled]
  
  set IMDCCPlotFrame [ttk::frame $w.hlf.n.f4.main.ccframe.plot]
  grid columnconfigure $IMDCCPlotFrame 0 -weight 1
  set IMDCCPlot [multiplot embed $w.hlf.n.f4.main.ccframe.plot -xsize 600 -ysize 400 -title "Real-Time Cross Correlation" -xlabel "Timestep" -ylabel "Cross Correlation" -lines -linewidth 2 -marker point -radius 2 -autoscale ]
  set imdccplot [$IMDCCPlot getpath]
  set ShowCC [ttk::label $w.hlf.n.f4.main.showcc -text "$rightPoint Cross Correlation Analysis..." -anchor w]
  set HideCC [ttk::label $w.hlf.n.f4.main.ccframe.hidecc -text "$downPoint Cross Correlation Analysis" -anchor w]
  $CCFrame configure -labelwidget $HideCC
  bind $HideCC <Button-1> {
      grid remove .mdffgui.hlf.n.f4.main.ccframe
      grid .mdffgui.hlf.n.f4.main.showcc
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowCC <Button-1> {
      grid remove .mdffgui.hlf.n.f4.main.showcc
      grid .mdffgui.hlf.n.f4.main.ccframe -row 2 -column 0 -columnspan 4 -sticky nsew -pady 5
      grid columnconfigure .mdffgui.hlf.n.f4.main.ccframe 1 -weight 1 
      #grid rowconfigure .mdffgui.hlf.n.f1.main 1 -weight 1 
      MDFFGUI::gui::resizeToActiveTab
  }
 # grid $IMDInputSwitch -row 5 -column 0 
 # grid $IMDHostInput -row 5 -column 1

  grid $IMDConnectFrame -row 0 -column 0 -sticky nsew
  grid $IMDButtonFrame -row 0 -column 0 -sticky nsew
  grid $IMDSubmit -row 0 -column 0 -sticky nsew
  grid $IMDConnect -row 0 -column 1 -sticky nsew
  grid $IMDPause -row 0 -column 2 -sticky nsew
  grid $IMDFinish -row 0 -column 3 -sticky nsew

  grid $IMDStatusLabel -row 1 -column 0 -columnspan 4 -sticky nsw
  grid $ShowCC -row 2 -column 0 -sticky nsew -pady 5 -columnspan 4 -padx $ShowBTNPadX
 
  grid $CButtonCC -row 0 -column 0 -sticky nsw -columnspan 4
   
  grid $CCMolLabel -row 1 -column 0 -sticky nsw
  grid $CCMolInput -row 1 -column 1 -sticky nsew
  
  grid $CCSelLabel -row 2 -column 0 -sticky nsw
  grid $CCSelInput -row 2 -column 1 -sticky nsew
  
  grid $CCResLabel -row 3 -column 0 -sticky nsw
  grid $CCResInput -row 3 -column 1 -sticky nsw
  
  grid $CButtonCCThresh -row 4 -column 0 -sticky nsw
  grid $CCThreshLabel -row 5 -column 0 -sticky nsw
  grid $CCThreshInput -row 5 -column 1 -sticky nsw

  grid $IMDCCPlotFrame -row 6 -column 0 -sticky nsew -columnspan 4
  grid $imdccplot -row 0 -column 0 -sticky nsew
  
  #Density Map Tools
  set NBTab5 [ttk::frame $w.hlf.n.f5];
  $Notebook add $NBTab5 -text "Map Tools" 
  grid columnconfigure $NBTab5 0 -weight 1
  grid rowconfigure $NBTab5 0 -weight 1
  
  set MapToolsFrame [ttk::frame $w.hlf.n.f5.main]
  grid columnconfigure $MapToolsFrame 1 -weight 1
  
  grid $MapToolsFrame -row 0 -column 0 -sticky nsew
  
  set MapToolsFileFrame [ttk::frame $w.hlf.n.f5.main.fileframe]
  #grid columnconfigure $MapToolsFileFrame 1 -weight 1 
 #Following is for using already loaded volumes
  set MapToolsMol [ttk::label $w.hlf.n.f5.main.fileframe.mollabel -text "Mol ID (Density Map):"]
  set MapToolsMolMenuButton [ttk::menubutton $w.hlf.n.f5.main.fileframe.molmenubutton -textvar MDFFGUI::gui::MapToolsMolMenuText -menu $w.hlf.n.f5.main.fileframe.molmenubutton.molmenu]
  set CurrentMapToolsMolMenu [menu $w.hlf.n.f5.main.fileframe.molmenubutton.molmenu -tearoff no]
  
  #set_mol_text 
  fill_mapmol_menu $CurrentMapToolsMolMenu
  trace add variable ::vmd_molecule write "MDFFGUI::gui::fill_mapmol_menu $CurrentMapToolsMolMenu"
  trace add variable MDFFGUI::maptoolssettings::MapToolsMolID write MDFFGUI::gui::set_mapmol_text
  set MDFFGUI::maptoolssettings::MapToolsMolID $nullMolString
  
  set MapToolsVol [ttk::label $w.hlf.n.f5.main.fileframe.vollabel -text "Vol ID (Density Map):"]
  set MapToolsVolMenuButton [ttk::menubutton $w.hlf.n.f5.main.fileframe.volmenubutton -textvar MDFFGUI::gui::MapToolsVolMenuText -menu $w.hlf.n.f5.main.fileframe.volmenubutton.volmenu]
  set CurrentMapToolsVolMenu [menu $w.hlf.n.f5.main.fileframe.volmenubutton.volmenu -tearoff no]
  
  trace add variable MDFFGUI::maptoolssettings::MapToolsMolID write MDFFGUI::gui::set_map_stats
  trace add variable MDFFGUI::maptoolssettings::MapToolsVolID write MDFFGUI::gui::set_map_stats
  
  #set_vol_text 
  fill_mapvol_menu $CurrentMapToolsVolMenu
  #nuclear option trace to capture volume data sets added to an existing molecule.
  trace add variable ::vmd_logfile write "MDFFGUI::gui::fill_mapvol_menu $CurrentMapToolsVolMenu"
  trace add variable MDFFGUI::maptoolssettings::MapToolsVolID write MDFFGUI::gui::set_mapvol_text
  trace add variable MDFFGUI::maptoolssettings::MapToolsMolID write "MDFFGUI::gui::fill_mapvol_menu $CurrentMapToolsVolMenu"
  trace add variable MDFFGUI::maptoolssettings::MapToolsMolID write MDFFGUI::gui::revert_volid

  set MDFFGUI::maptoolssettings::MapToolsVolID $nullVolString
  
  set MapToolsSelectMap [ttk::button $w.hlf.n.f5.main.fileframe.mapbutton -text "Browse" -command {MDFFGUI::gui::get_map} ]
  set MapToolsSaveMap [ttk::button $w.hlf.n.f5.main.fileframe.savemapbutton -text "Save" -command {MDFFGUI::gui::save_map} ]
 
 
  grid $MapToolsFileFrame -row 0 -column 0 -sticky nsw -columnspan 10
  grid $MapToolsMol -row 0 -column 0 -sticky nsw
  grid $MapToolsMolMenuButton   -row 0 -column 1 -sticky nsw
  grid $MapToolsVol -row 0 -column 2 -sticky nsw
  grid $MapToolsVolMenuButton   -row 0 -column 3 -sticky nsw
  grid $MapToolsSelectMap   -row 0 -column 4 -sticky nsw
  grid $MapToolsSaveMap   -row 0 -column 5 -sticky nsw
  
  set MapToolsInfoFrame [ttk::labelframe $w.hlf.n.f5.main.infoframe -labelanchor nw]
  #grid columnconfigure $MapToolsInfoFrame 1 -weight 1
  set MapToolsMapOriginLabel [ttk::label $w.hlf.n.f5.main.infoframe.originlabel -text "Origin: "]
  set MapToolsMapOriginLabelDisp [ttk::label $w.hlf.n.f5.main.infoframe.originlabeldisp -textvar MDFFGUI::maptoolssettings::MapToolsOrigin]
  set MapToolsMapCOMLabel [ttk::label $w.hlf.n.f5.main.infoframe.comlabel -text "Center of Mass: "]
  set MapToolsMapCOMLabelDisp [ttk::label $w.hlf.n.f5.main.infoframe.comlabeldisp -textvar MDFFGUI::maptoolssettings::MapToolsCOM]
  set MapToolsMapXsizeLabel [ttk::label $w.hlf.n.f5.main.infoframe.xsizelabel -text "X Size: "]
  set MapToolsMapXsizeLabelDisp [ttk::label $w.hlf.n.f5.main.infoframe.xsizelabeldisp -textvar MDFFGUI::maptoolssettings::MapToolsXsize]
  set MapToolsMapYsizeLabel [ttk::label $w.hlf.n.f5.main.infoframe.ysizelabel -text "Y Size: "]
  set MapToolsMapYsizeLabelDisp [ttk::label $w.hlf.n.f5.main.infoframe.ysizelabeldisp -textvar MDFFGUI::maptoolssettings::MapToolsYsize]
  set MapToolsMapZsizeLabel [ttk::label $w.hlf.n.f5.main.infoframe.zsizelabel -text "Z Size: "]
  set MapToolsMapZsizeLabelDisp [ttk::label $w.hlf.n.f5.main.infoframe.zsizelabeldisp -textvar MDFFGUI::maptoolssettings::MapToolsZsize]
  set MapToolsMapMinLabel [ttk::label $w.hlf.n.f5.main.infoframe.minlabel -text "Min Value: "]
  set MapToolsMapMinLabelDisp [ttk::label $w.hlf.n.f5.main.infoframe.minlabeldisp -textvar MDFFGUI::maptoolssettings::MapToolsMin]
  set MapToolsMapMaxLabel [ttk::label $w.hlf.n.f5.main.infoframe.maxlabel -text "Max Value: "]
  set MapToolsMapMaxLabelDisp [ttk::label $w.hlf.n.f5.main.infoframe.maxlabeldisp -textvar MDFFGUI::maptoolssettings::MapToolsMax]
  
  set ShowMapToolsInfo [ttk::label $w.hlf.n.f5.main.showinfo -text "$rightPoint Map Info..." -anchor w]
  set HideMapToolsInfo [ttk::label $w.hlf.n.f5.main.infoframe.hideinfo -text "$downPoint Map Info" -anchor w]

  $MapToolsInfoFrame configure -labelwidget $HideMapToolsInfo
  bind $HideMapToolsInfo <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.infoframe
      grid .mdffgui.hlf.n.f5.main.showinfo
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowMapToolsInfo <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.showinfo
      grid .mdffgui.hlf.n.f5.main.infoframe -row 1 -column 0 -sticky nsw -pady 5 -columnspan 10
#      grid columnconfigure .mdffgui.hlf.n.f5.main.infoframe 1 -weight 1
      MDFFGUI::gui::resizeToActiveTab
  }
  
  
  grid $ShowMapToolsInfo -row 1 -column 0 -sticky nswe -pady 5 -padx $ShowBTNPadX
  #grid $MapToolsInfoFrame -row 1 -column 0 -sticky nsew
  grid $MapToolsMapOriginLabel -row 0 -column 0 -sticky nsw 
  grid $MapToolsMapOriginLabelDisp -row 0 -column 1 -sticky nsw
  grid $MapToolsMapCOMLabel -row 1 -column 0 -sticky nsw 
  grid $MapToolsMapCOMLabelDisp -row 1 -column 1 -sticky nsw -columnspan 5
  grid $MapToolsMapXsizeLabel -row 2 -column 0 -sticky nsw 
  grid $MapToolsMapXsizeLabelDisp -row 2 -column 1 -sticky nsw
  grid $MapToolsMapYsizeLabel -row 2 -column 2 -sticky nsw 
  grid $MapToolsMapYsizeLabelDisp -row 2 -column 3 -sticky nsw
  grid $MapToolsMapZsizeLabel -row 2 -column 4 -sticky nsw 
  grid $MapToolsMapZsizeLabelDisp -row 2 -column 5 -sticky nsw
  grid $MapToolsMapMinLabel -row 3 -column 0 -sticky nsw
  grid $MapToolsMapMinLabelDisp -row 3 -column 1 -sticky nsw
  grid $MapToolsMapMaxLabel -row 3 -column 2 -sticky nsw
  grid $MapToolsMapMaxLabelDisp -row 3 -column 3 -sticky nsw
  
  set HistPlotFrame [ttk::labelframe $w.hlf.n.f5.main.plot -labelanchor nw]
 # grid columnconfigure $HistPlotFrame 1 -weight 1  
  set GenerateHistPlot [ttk::button $w.hlf.n.f5.main.plot.histbutton -text "Generate Histogram" -command {MDFFGUI::gui::generate_histogram} -state enabled]
  set HistPlotBinLabel [ttk::label $w.hlf.n.f5.main.plot.histbinlabel -text "Number of Bins:"]
  set HistPlotBinEntry [ttk::entry $w.hlf.n.f5.main.plot.histbinentry -textvariable MDFFGUI::maptoolssettings::HistPlotBins -width 5]
  set HistPlotLog [ttk::radiobutton $w.hlf.n.f5.main.plot.log -variable MDFFGUI::maptoolssettings::HistLog -value 1 -text "Log10"]
  set HistPlotNot [ttk::radiobutton $w.hlf.n.f5.main.plot.not -variable MDFFGUI::maptoolssettings::HistLog -value 0 -text "Unchanged"]

  global MapToolsPlotX
  set MapToolsPlotX 0
  set MapToolsPlotXLabel [ttk::label $w.hlf.n.f5.main.plot.plotxlabel -text "Selected Voxel Value: "]
  set MapToolsPlotXLabelX [ttk::label $w.hlf.n.f5.main.plot.plotxlabelx -textvar MapToolsPlotX -width 5]
  
  set ShowMapToolsHist [ttk::label $w.hlf.n.f5.main.showhist -text "$rightPoint Histogram..." -anchor w]
  set HideMapToolsHist [ttk::label $w.hlf.n.f5.main.plot.hidehist -text "$downPoint Histogram" -anchor w]

  $HistPlotFrame configure -labelwidget $HideMapToolsHist
  bind $HideMapToolsHist <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.plot
      grid .mdffgui.hlf.n.f5.main.showhist
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowMapToolsHist <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.showhist
      grid .mdffgui.hlf.n.f5.main.plot -row 2 -column 0 -sticky nsw -pady 5 -columnspan 10
   #   grid columnconfigure .mdffgui.hlf.n.f5.main.plot 1 -weight 1
      MDFFGUI::gui::resizeToActiveTab
  }
   
  grid $ShowMapToolsHist -row 2 -column 0 -sticky nswe -pady 5 -padx $ShowBTNPadX
  grid $GenerateHistPlot -row 0 -column 0 -sticky nsew
  grid $HistPlotBinLabel -row 0 -column 1 -sticky nsew
  grid $HistPlotBinEntry -row 0 -column 2 -sticky nsew
  grid $HistPlotLog -row 0 -column 3 -sticky nsew
  grid $HistPlotNot -row 0 -column 4 -sticky nsew
  grid $MapToolsPlotXLabel -row 1 -column 0 -sticky nsew
  grid $MapToolsPlotXLabelX -row 1 -column 1 -sticky nsew

  #structure ops
  set MapToolsStructFrame [ttk::labelframe $w.hlf.n.f5.main.structframe -labelanchor nw]

  set MapToolsStructMol [ttk::label $w.hlf.n.f5.main.structframe.mollabel -text "Mol ID (Structure):"]
  set MapToolsStructMolMenuButton [ttk::menubutton $w.hlf.n.f5.main.structframe.molmenubutton -textvar MDFFGUI::gui::MapToolsStructMolMenuText -menu $w.hlf.n.f5.main.structframe.molmenubutton.molmenu]
  set CurrentMapToolsStructMolMenu [menu $w.hlf.n.f5.main.structframe.molmenubutton.molmenu -tearoff no]
  
  #set_mol_text 
  fill_mapstructmol_menu $CurrentMapToolsStructMolMenu
  trace add variable ::vmd_molecule write "MDFFGUI::gui::fill_mapstructmol_menu $CurrentMapToolsStructMolMenu"
  trace add variable MDFFGUI::maptoolssettings::MapToolsStructMolID write MDFFGUI::gui::set_mapstructmol_text
  set MDFFGUI::maptoolssettings::MapToolsStructMolID $nullMolString
  
  set MapToolsSelectStruct [ttk::button $w.hlf.n.f5.main.structframe.structbutton -text "Browse" -command {MDFFGUI::gui::get_maptools_struct} ]

  #Fit
  set MapToolsFitLabel [ttk::label $w.hlf.n.f5.main.structframe.fitlabel -text "Rigid Body Fitting:"]
  set MapToolsFitLabelSel [ttk::label $w.hlf.n.f5.main.structframe.fitlabelsel -text "atom selection:"]
  set MapToolsFitEntrySel [ttk::entry $w.hlf.n.f5.main.structframe.fitentrysel -textvariable MDFFGUI::maptoolssettings::FitSel -width 20]
  set MapToolsFitLabelRes [ttk::label $w.hlf.n.f5.main.structframe.fitlabelres -text "map resolution (A):"]
  set MapToolsFitEntryRes [ttk::entry $w.hlf.n.f5.main.structframe.fitentryres -textvariable MDFFGUI::maptoolssettings::FitRes -width 5]
  set MapToolsFitButton [ttk::button $w.hlf.n.f5.main.structframe.fitbutton -text "Fit" -command { voltool fit [atomselect $MDFFGUI::maptoolssettings::MapToolsStructMolID $MDFFGUI::maptoolssettings::FitSel] -res $MDFFGUI::maptoolssettings::FitRes -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]

  #Mask
  set MapToolsMaskLabel [ttk::label $w.hlf.n.f5.main.structframe.masklabel -text "Map Masking:"]
  set MapToolsMaskLabelSel [ttk::label $w.hlf.n.f5.main.structframe.masklabelsel -text "atom selection:"]
  set MapToolsMaskEntrySel [ttk::entry $w.hlf.n.f5.main.structframe.maskentrysel -textvariable MDFFGUI::maptoolssettings::MaskSel -width 20]
  set MapToolsMaskLabelRes [ttk::label $w.hlf.n.f5.main.structframe.masklabelres -text "map resolution (A):"]
  set MapToolsMaskEntryRes [ttk::entry $w.hlf.n.f5.main.structframe.maskentryres -textvariable MDFFGUI::maptoolssettings::MaskRes -width 5]
 # set MapToolsMaskLabelCutoff [ttk::label $w.hlf.n.f5.main.structframe.masklabelcutoff -text "cutoff (A):"]
 # set MapToolsMaskEntryCutoff [ttk::entry $w.hlf.n.f5.main.structframe.maskentrycutoff -textvariable MDFFGUI::settings::MaskCutoff -width 5]
  set MapToolsMaskButton [ttk::button $w.hlf.n.f5.main.structframe.maskbutton -text "Mask" -command { voltool mask [atomselect $MDFFGUI::maptoolssettings::MapToolsStructMolID $MDFFGUI::maptoolssettings::MaskSel] -res $MDFFGUI::maptoolssettings::MaskRes -cutoff $MDFFGUI::maptoolssettings::MaskRes -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
 
  #Sim
  set MapToolsSimLabel [ttk::label $w.hlf.n.f5.main.structframe.simlabel -text "Simulated Density:"]
  set MapToolsSimLabelSel [ttk::label $w.hlf.n.f5.main.structframe.simlabelsel -text "atom selection:"]
  set MapToolsSimEntrySel [ttk::entry $w.hlf.n.f5.main.structframe.simentrysel -textvariable MDFFGUI::maptoolssettings::SimSel -width 20]
  set MapToolsSimLabelRes [ttk::label $w.hlf.n.f5.main.structframe.simlabelres -text "map resolution (A):"]
  set MapToolsSimEntryRes [ttk::entry $w.hlf.n.f5.main.structframe.simentryres -textvariable MDFFGUI::maptoolssettings::SimRes -width 5]
  set MapToolsSimButton [ttk::button $w.hlf.n.f5.main.structframe.simbutton -text "Generate" -command { voltool sim [atomselect $MDFFGUI::maptoolssettings::MapToolsStructMolID $MDFFGUI::maptoolssettings::SimSel] -res $MDFFGUI::maptoolssettings::SimRes } ]
  
  #CC
  set MapToolsCCLabel [ttk::label $w.hlf.n.f5.main.structframe.cclabel -text "Cross Correlation:"]
  set MapToolsCCLabelSel [ttk::label $w.hlf.n.f5.main.structframe.cclabelsel -text "atom selection:"]
  set MapToolsCCEntrySel [ttk::entry $w.hlf.n.f5.main.structframe.ccentrysel -textvariable MDFFGUI::maptoolssettings::CCSel -width 20]
  set MapToolsCCLabelRes [ttk::label $w.hlf.n.f5.main.structframe.cclabelres -text "map resolution (A):"]
  set MapToolsCCEntryRes [ttk::entry $w.hlf.n.f5.main.structframe.ccentryres -textvariable MDFFGUI::maptoolssettings::CCRes -width 5]
  set MapToolsCCButton [ttk::button $w.hlf.n.f5.main.structframe.ccbutton -text "Calculate" -command { set MDFFGUI::gui::CCOutput [voltool cc [atomselect $MDFFGUI::maptoolssettings::MapToolsStructMolID $MDFFGUI::maptoolssettings::CCSel] -res $MDFFGUI::maptoolssettings::CCRes -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID] } ]
  set MapToolsCCLabelOutput [ttk::label $w.hlf.n.f5.main.structframe.cclabelout -text "CC:"]
  set MapToolsCCOutput [ttk::label $w.hlf.n.f5.main.structframe.ccout -textvariable MDFFGUI::gui::CCOutput -width 10]
  
  set ShowMapToolsStruct [ttk::label $w.hlf.n.f5.main.showstruct -text "$rightPoint Structure Ops..." -anchor w]
  set HideMapToolsStruct [ttk::label $w.hlf.n.f5.main.structframe.hidestruct -text "$downPoint Structure Ops" -anchor w]

  $MapToolsStructFrame configure -labelwidget $HideMapToolsStruct
  bind $HideMapToolsStruct <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.structframe
      grid .mdffgui.hlf.n.f5.main.showstruct
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowMapToolsStruct <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.showstruct
      grid .mdffgui.hlf.n.f5.main.structframe -row 3 -column 0 -sticky nsw -pady 5 -columnspan 10
 #     grid columnconfigure .mdffgui.hlf.n.f5.main.structframe 1 -weight 1
      MDFFGUI::gui::resizeToActiveTab
  }
  
  grid $ShowMapToolsStruct -row 3 -column 0 -sticky nswe
 
  grid $MapToolsStructMol -row 0 -column 0 -sticky nswe
  grid $MapToolsStructMolMenuButton -row 0 -column 1 -sticky nswe
  grid $MapToolsSelectStruct -row 0 -column 2 -sticky nswe

  grid $MapToolsFitLabel -row 1 -column 0 -sticky nswe
  grid $MapToolsFitLabelSel -row 1 -column 1 -sticky nswe
  grid $MapToolsFitEntrySel -row 1 -column 2 -sticky nswe
  grid $MapToolsFitLabelRes -row 1 -column 3 -sticky nswe
  grid $MapToolsFitEntryRes -row 1 -column 4 -sticky nswe
  grid $MapToolsFitButton -row 1 -column 5 -sticky nswe
  
  grid $MapToolsMaskLabel -row 2 -column 0 -sticky nswe
  grid $MapToolsMaskLabelSel -row 2 -column 1 -sticky nswe
  grid $MapToolsMaskEntrySel -row 2 -column 2 -sticky nswe
  grid $MapToolsMaskLabelRes -row 2 -column 3 -sticky nswe
  grid $MapToolsMaskEntryRes -row 2 -column 4 -sticky nswe
  grid $MapToolsMaskButton -row 2 -column 5 -sticky nswe
  
  grid $MapToolsSimLabel -row 3 -column 0 -sticky nswe
  grid $MapToolsSimLabelSel -row 3 -column 1 -sticky nswe
  grid $MapToolsSimEntrySel -row 3 -column 2 -sticky nswe
  grid $MapToolsSimLabelRes -row 3 -column 3 -sticky nswe
  grid $MapToolsSimEntryRes -row 3 -column 4 -sticky nswe
  grid $MapToolsSimButton -row 3 -column 5 -sticky nswe
  
  grid $MapToolsCCLabel -row 4 -column 0 -sticky nswe
  grid $MapToolsCCLabelSel -row 4 -column 1 -sticky nswe
  grid $MapToolsCCEntrySel -row 4 -column 2 -sticky nswe
  grid $MapToolsCCLabelRes -row 4 -column 3 -sticky nswe
  grid $MapToolsCCEntryRes -row 4 -column 4 -sticky nswe
  grid $MapToolsCCButton -row 4 -column 5 -sticky nswe
  grid $MapToolsCCLabelOutput -row 4 -column 6 -sticky nswe
  grid $MapToolsCCOutput -row 4 -column 7 -sticky nsw
  
  #unary ops
  set MapToolsUnaryFrame [ttk::labelframe $w.hlf.n.f5.main.unaryframe -labelanchor nw]
  
  #segmentation
  set MapToolsSegLabel [ttk::label $w.hlf.n.f5.main.unaryframe.seglabel -text "Segmentation:"]
  set MapToolsSegLabelAmt [ttk::label $w.hlf.n.f5.main.unaryframe.seglabelamt -text "Segments:"]
  set MapToolsSegEntryAmt [ttk::entry $w.hlf.n.f5.main.unaryframe.segentryamt -textvariable MDFFGUI::maptoolssettings::SegAmt -width 5]
  set MapToolsSegButton [ttk::button $w.hlf.n.f5.main.unaryframe.segbutton -text "Segment" -command { segmentation segment -groups $MDFFGUI::maptoolssettings::SegAmt -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID -separate_groups } ]
  
  #Trim
  set MapToolsTrimLabel [ttk::label $w.hlf.n.f5.main.unaryframe.trimlabel -text "Trim map in each direction:"]
  set MapToolsTrimLabelX1 [ttk::label $w.hlf.n.f5.main.unaryframe.trimlabelx1 -text "-x:"]
  set MapToolsTrimEntryX1 [ttk::entry $w.hlf.n.f5.main.unaryframe.trimentryx1 -textvariable MDFFGUI::maptoolssettings::TrimX1 -width 5]
  set MapToolsTrimLabelX2 [ttk::label $w.hlf.n.f5.main.unaryframe.trimlabelx2 -text "+x:"]
  set MapToolsTrimEntryX2 [ttk::entry $w.hlf.n.f5.main.unaryframe.trimentryx2 -textvariable MDFFGUI::maptoolssettings::TrimX2 -width 5]
  set MapToolsTrimLabelY1 [ttk::label $w.hlf.n.f5.main.unaryframe.trimlabely1 -text "-y:"]
  set MapToolsTrimEntryY1 [ttk::entry $w.hlf.n.f5.main.unaryframe.trimentryy1 -textvariable MDFFGUI::maptoolssettings::TrimY1 -width 5]
  set MapToolsTrimLabelY2 [ttk::label $w.hlf.n.f5.main.unaryframe.trimlabely2 -text "+y:"]
  set MapToolsTrimEntryY2 [ttk::entry $w.hlf.n.f5.main.unaryframe.trimentryy2 -textvariable MDFFGUI::maptoolssettings::TrimY2 -width 5]
  set MapToolsTrimLabelZ1 [ttk::label $w.hlf.n.f5.main.unaryframe.trimlabelz1 -text "-z:"]
  set MapToolsTrimEntryZ1 [ttk::entry $w.hlf.n.f5.main.unaryframe.trimentryz1 -textvariable MDFFGUI::maptoolssettings::TrimZ1 -width 5]
  set MapToolsTrimLabelZ2 [ttk::label $w.hlf.n.f5.main.unaryframe.trimlabelz2 -text "+z:"]
  set MapToolsTrimEntryZ2 [ttk::entry $w.hlf.n.f5.main.unaryframe.trimentryz2 -textvariable MDFFGUI::maptoolssettings::TrimZ2 -width 5]
  set MapToolsTrimButton [ttk::button $w.hlf.n.f5.main.unaryframe.trimbutton -text "Trim" -command { voltool trim -amt "$MDFFGUI::maptoolssettings::TrimX1 $MDFFGUI::maptoolssettings::TrimX2 $MDFFGUI::maptoolssettings::TrimY1 $MDFFGUI::maptoolssettings::TrimY2 $MDFFGUI::maptoolssettings::TrimZ1 $MDFFGUI::maptoolssettings::TrimZ2" -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  set MapToolsCropLabel [ttk::label $w.hlf.n.f5.main.unaryframe.croplabel -text "Crop map in coordinate space:"]
  set MapToolsCropLabelX1 [ttk::label $w.hlf.n.f5.main.unaryframe.croplabelx1 -text "minx:"]
  set MapToolsCropEntryX1 [ttk::entry $w.hlf.n.f5.main.unaryframe.cropentryx1 -textvariable MDFFGUI::maptoolssettings::CropX1 -width 5]
  set MapToolsCropLabelX2 [ttk::label $w.hlf.n.f5.main.unaryframe.croplabelx2 -text "maxx:"]
  set MapToolsCropEntryX2 [ttk::entry $w.hlf.n.f5.main.unaryframe.cropentryx2 -textvariable MDFFGUI::maptoolssettings::CropX2 -width 5]
  set MapToolsCropLabelY1 [ttk::label $w.hlf.n.f5.main.unaryframe.croplabely1 -text "miny:"]
  set MapToolsCropEntryY1 [ttk::entry $w.hlf.n.f5.main.unaryframe.cropentryy1 -textvariable MDFFGUI::maptoolssettings::CropY1 -width 5]
  set MapToolsCropLabelY2 [ttk::label $w.hlf.n.f5.main.unaryframe.croplabely2 -text "maxy:"]
  set MapToolsCropEntryY2 [ttk::entry $w.hlf.n.f5.main.unaryframe.cropentryy2 -textvariable MDFFGUI::maptoolssettings::CropY2 -width 5]
  set MapToolsCropLabelZ1 [ttk::label $w.hlf.n.f5.main.unaryframe.croplabelz1 -text "minz:"]
  set MapToolsCropEntryZ1 [ttk::entry $w.hlf.n.f5.main.unaryframe.cropentryz1 -textvariable MDFFGUI::maptoolssettings::CropZ1 -width 5]
  set MapToolsCropLabelZ2 [ttk::label $w.hlf.n.f5.main.unaryframe.croplabelz2 -text "maxz:"]
  set MapToolsCropEntryZ2 [ttk::entry $w.hlf.n.f5.main.unaryframe.cropentryz2 -textvariable MDFFGUI::maptoolssettings::CropZ2 -width 5]
  set MapToolsCropButton [ttk::button $w.hlf.n.f5.main.unaryframe.cropbutton -text "Crop" -command { voltool crop -amt "$MDFFGUI::maptoolssettings::CropX1 $MDFFGUI::maptoolssettings::CropY1 $MDFFGUI::maptoolssettings::CropZ1 $MDFFGUI::maptoolssettings::CropX2 $MDFFGUI::maptoolssettings::CropY2 $MDFFGUI::maptoolssettings::CropZ2" -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #Clamp
  set MapToolsClampLabel [ttk::label $w.hlf.n.f5.main.unaryframe.clamplabel -text "Clamp map to range:"]
  set MapToolsClampLabelMin [ttk::label $w.hlf.n.f5.main.unaryframe.clamplabelmin -text "min:"]
  set MapToolsClampEntryMin [ttk::entry $w.hlf.n.f5.main.unaryframe.clampentrymin -textvariable MDFFGUI::maptoolssettings::ClampMin -width 5]
  set MapToolsClampLabelMax [ttk::label $w.hlf.n.f5.main.unaryframe.clamplabelmax -text "max:"]
  set MapToolsClampEntryMax [ttk::entry $w.hlf.n.f5.main.unaryframe.clampentrymax -textvariable MDFFGUI::maptoolssettings::ClampMax -width 5]
  set MapToolsClampButton [ttk::button $w.hlf.n.f5.main.unaryframe.clampbutton -text "Clamp" -command { voltool clamp -min $MDFFGUI::maptoolssettings::ClampMin -max $MDFFGUI::maptoolssettings::ClampMax -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #Smult
  set MapToolsSmultLabel [ttk::label $w.hlf.n.f5.main.unaryframe.smultlabel -text "Multiply by scalar:"]
  set MapToolsSmultLabelAmt [ttk::label $w.hlf.n.f5.main.unaryframe.smultlabelamt -text "Amount:"]
  set MapToolsSmultEntryAmt [ttk::entry $w.hlf.n.f5.main.unaryframe.smultentryamt -textvariable MDFFGUI::maptoolssettings::SmultAmt -width 5]
  set MapToolsSmultButton [ttk::button $w.hlf.n.f5.main.unaryframe.smultbutton -text "Mult" -command { voltool smult -amt $MDFFGUI::maptoolssettings::SmultAmt -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #SAdd
  set MapToolsSAddLabel [ttk::label $w.hlf.n.f5.main.unaryframe.saddlabel -text "Add scalar:"]
  set MapToolsSAddLabelAmt [ttk::label $w.hlf.n.f5.main.unaryframe.saddlabelamt -text "Amount:"]
  set MapToolsSAddEntryAmt [ttk::entry $w.hlf.n.f5.main.unaryframe.saddentryamt -textvariable MDFFGUI::maptoolssettings::SAddAmt -width 5]
  set MapToolsSAddButton [ttk::button $w.hlf.n.f5.main.unaryframe.saddbutton -text "Add" -command { voltool sadd -amt $MDFFGUI::maptoolssettings::SAddAmt -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #Range
  set MapToolsRangeLabel [ttk::label $w.hlf.n.f5.main.unaryframe.rangelabel -text "Fit values to range:"]
  set MapToolsRangeLabelMin [ttk::label $w.hlf.n.f5.main.unaryframe.rangelabelmin -text "min:"]
  set MapToolsRangeEntryMin [ttk::entry $w.hlf.n.f5.main.unaryframe.rangeentrymin -textvariable MDFFGUI::maptoolssettings::RangeMin -width 5]
  set MapToolsRangeLabelMax [ttk::label $w.hlf.n.f5.main.unaryframe.rangelabelmax -text "max:"]
  set MapToolsRangeEntryMax [ttk::entry $w.hlf.n.f5.main.unaryframe.rangeentrymax -textvariable MDFFGUI::maptoolssettings::RangeMax -width 5]
  set MapToolsRangeButton [ttk::button $w.hlf.n.f5.main.unaryframe.rangebutton -text "Fit" -command { voltool range -minmax "$MDFFGUI::maptoolssettings::RangeMin $MDFFGUI::maptoolssettings::RangeMax" -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]

  #Binmask
  set MapToolsBinmaskLabel [ttk::label $w.hlf.n.f5.main.unaryframe.binmasklabel -text "Binary Mask:"]
  set MapToolsBinmaskLabelAmt [ttk::label $w.hlf.n.f5.main.unaryframe.binmasklabelamt -text "Threshold:"]
  set MapToolsBinmaskEntryAmt [ttk::entry $w.hlf.n.f5.main.unaryframe.binmaskentryamt -textvariable MDFFGUI::maptoolssettings::BinmaskAmt -width 5]
  set MapToolsBinmaskButton [ttk::button $w.hlf.n.f5.main.unaryframe.binmaskbutton -text "Mask" -command { voltool binmask -threshold $MDFFGUI::maptoolssettings::BinmaskAmt -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #Smooth
  set MapToolsSmoothLabel [ttk::label $w.hlf.n.f5.main.unaryframe.smoothlabel -text "Guassian Blur:"]
  set MapToolsSmoothLabelAmt [ttk::label $w.hlf.n.f5.main.unaryframe.smoothlabelamt -text "Sigma:"]
  set MapToolsSmoothEntryAmt [ttk::entry $w.hlf.n.f5.main.unaryframe.smoothentryamt -textvariable MDFFGUI::maptoolssettings::SmoothAmt -width 5]
  set MapToolsSmoothButton [ttk::button $w.hlf.n.f5.main.unaryframe.smoothbutton -text "Smooth" -command { voltool smooth -sigma $MDFFGUI::maptoolssettings::SmoothAmt -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #Pot
  set MapToolsPotLabel [ttk::label $w.hlf.n.f5.main.unaryframe.potlabel -text "MDFF Potential:"]
  set MapToolsPotLabelAmt [ttk::label $w.hlf.n.f5.main.unaryframe.potlabelamt -text "Threshold:"]
  set MapToolsPotEntryAmt [ttk::entry $w.hlf.n.f5.main.unaryframe.potentryamt -textvariable MDFFGUI::maptoolssettings::PotAmt -width 5]
  set MapToolsPotButton [ttk::button $w.hlf.n.f5.main.unaryframe.potbutton -text "Convert" -command { voltool pot -threshold $MDFFGUI::maptoolssettings::PotAmt -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #Downsample
  set MapToolsDownsampleLabel [ttk::label $w.hlf.n.f5.main.unaryframe.downsamplelabel -text "Downsample by x2:"]
  set MapToolsDownsampleButton [ttk::button $w.hlf.n.f5.main.unaryframe.downsamplebutton -text "Downsample" -command { voltool downsample -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #Supersample
  set MapToolsSupersampleLabel [ttk::label $w.hlf.n.f5.main.unaryframe.supersamplelabel -text "Supersample by x2:"]
  set MapToolsSupersampleButton [ttk::button $w.hlf.n.f5.main.unaryframe.supersamplebutton -text "Supersample" -command { voltool supersample -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  #Sigma
  set MapToolsSigmaLabel [ttk::label $w.hlf.n.f5.main.unaryframe.sigmalabel -text "Convert to Sigma Scale:"]
  set MapToolsSigmaButton [ttk::button $w.hlf.n.f5.main.unaryframe.sigmabutton -text "Convert" -command { voltool sigma -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID } ]
  
  set ShowMapToolsUnary [ttk::label $w.hlf.n.f5.main.showunary -text "$rightPoint Single Map Ops..." -anchor w]
  set HideMapToolsUnary [ttk::label $w.hlf.n.f5.main.unaryframe.hideunary -text "$downPoint Single Map Ops" -anchor w]

  $MapToolsUnaryFrame configure -labelwidget $HideMapToolsUnary
  bind $HideMapToolsUnary <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.unaryframe
      grid .mdffgui.hlf.n.f5.main.showunary
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowMapToolsUnary <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.showunary
      grid .mdffgui.hlf.n.f5.main.unaryframe -row 4 -column 0 -sticky nsw -pady 5 -columnspan 10
  #    grid columnconfigure .mdffgui.hlf.n.f5.main.unaryframe 1 -weight 1
      MDFFGUI::gui::resizeToActiveTab
  }
  
  grid $ShowMapToolsUnary -row 4 -column 0 -sticky nswe
  
  grid $MapToolsTrimLabel -row 0 -column 0 -sticky nsw 
  grid $MapToolsTrimLabelX1 -row 0 -column 1 -sticky nsw 
  grid $MapToolsTrimEntryX1 -row 0 -column 2 -sticky nsw 
  grid $MapToolsTrimLabelX2 -row 0 -column 3 -sticky nsw 
  grid $MapToolsTrimEntryX2 -row 0 -column 4 -sticky nsw 
  grid $MapToolsTrimLabelY1 -row 0 -column 5 -sticky nsw 
  grid $MapToolsTrimEntryY1 -row 0 -column 6 -sticky nsw 
  grid $MapToolsTrimLabelY2 -row 0 -column 7 -sticky nsw 
  grid $MapToolsTrimEntryY2 -row 0 -column 8 -sticky nsw 
  grid $MapToolsTrimLabelZ1 -row 0 -column 9 -sticky nsw 
  grid $MapToolsTrimEntryZ1 -row 0 -column 10 -sticky nsw 
  grid $MapToolsTrimLabelZ2 -row 0 -column 11 -sticky nsw 
  grid $MapToolsTrimEntryZ2 -row 0 -column 12 -sticky nsw 
  grid $MapToolsTrimButton -row 0 -column 13 -sticky nsw -columnspan 2
  
  grid $MapToolsCropLabel -row 1 -column 0 -sticky nswe 
  grid $MapToolsCropLabelX1 -row 1 -column 1 -sticky nswe 
  grid $MapToolsCropEntryX1 -row 1 -column 2 -sticky nswe 
  grid $MapToolsCropLabelX2 -row 1 -column 3 -sticky nswe 
  grid $MapToolsCropEntryX2 -row 1 -column 4 -sticky nswe 
  grid $MapToolsCropLabelY1 -row 1 -column 5 -sticky nswe 
  grid $MapToolsCropEntryY1 -row 1 -column 6 -sticky nswe 
  grid $MapToolsCropLabelY2 -row 1 -column 7 -sticky nswe 
  grid $MapToolsCropEntryY2 -row 1 -column 8 -sticky nswe 
  grid $MapToolsCropLabelZ1 -row 1 -column 9 -sticky nswe 
  grid $MapToolsCropEntryZ1 -row 1 -column 10 -sticky nswe 
  grid $MapToolsCropLabelZ2 -row 1 -column 11 -sticky nswe 
  grid $MapToolsCropEntryZ2 -row 1 -column 12 -sticky nswe 
  grid $MapToolsCropButton -row 1 -column 13 -sticky nsw -columnspan 2
  
  grid $MapToolsClampLabel -row 2 -column 0 -sticky nswe 
  grid $MapToolsClampLabelMin -row 2 -column 1 -sticky nswe 
  grid $MapToolsClampEntryMin -row 2 -column 2 -sticky nswe 
  grid $MapToolsClampLabelMax -row 2 -column 3 -sticky nswe 
  grid $MapToolsClampEntryMax -row 2 -column 4 -sticky nswe 
  grid $MapToolsClampButton -row 2 -column 5 -sticky nsw -columnspan 2
  
  grid $MapToolsSmultLabel -row 3 -column 0 -sticky nswe 
  grid $MapToolsSmultLabelAmt -row 3 -column 1 -sticky nswe 
  grid $MapToolsSmultEntryAmt -row 3 -column 2 -sticky nswe 
  grid $MapToolsSmultButton -row 3 -column 3 -sticky nsw -columnspan 2
  
  grid $MapToolsSAddLabel -row 4 -column 0 -sticky nswe 
  grid $MapToolsSAddLabelAmt -row 4 -column 1 -sticky nswe 
  grid $MapToolsSAddEntryAmt -row 4 -column 2 -sticky nswe 
  grid $MapToolsSAddButton -row 4 -column 3 -sticky nsw -columnspan 2
  
  grid $MapToolsRangeLabel -row 5 -column 0 -sticky nswe 
  grid $MapToolsRangeLabelMin -row 5 -column 1 -sticky nswe 
  grid $MapToolsRangeEntryMin -row 5 -column 2 -sticky nswe 
  grid $MapToolsRangeLabelMax -row 5 -column 3 -sticky nswe 
  grid $MapToolsRangeEntryMax -row 5 -column 4 -sticky nswe 
  grid $MapToolsRangeButton -row 5 -column 5 -sticky nsw -columnspan 2
  
  grid $MapToolsBinmaskLabel -row 6 -column 0 -sticky nswe 
  grid $MapToolsBinmaskLabelAmt -row 6 -column 1 -sticky nswe 
  grid $MapToolsBinmaskEntryAmt -row 6 -column 2 -sticky nswe 
  grid $MapToolsBinmaskButton -row 6 -column 3 -sticky nsw -columnspan 2
  
  grid $MapToolsSmoothLabel -row 7 -column 0 -sticky nswe 
  grid $MapToolsSmoothLabelAmt -row 7 -column 1 -sticky nswe 
  grid $MapToolsSmoothEntryAmt -row 7 -column 2 -sticky nswe 
  grid $MapToolsSmoothButton -row 7 -column 3 -sticky nsw -columnspan 2
  
  grid $MapToolsPotLabel -row 8 -column 0 -sticky nswe 
  grid $MapToolsPotLabelAmt -row 8 -column 1 -sticky nswe 
  grid $MapToolsPotEntryAmt -row 8 -column 2 -sticky nswe 
  grid $MapToolsPotButton -row 8 -column 3 -sticky nsw -columnspan 2
  
  grid $MapToolsDownsampleLabel -row 9 -column 0 -sticky nswe 
  grid $MapToolsDownsampleButton -row 9 -column 1 -sticky nsw -columnspan 2

  grid $MapToolsSupersampleLabel -row 10 -column 0 -sticky nswe 
  grid $MapToolsSupersampleButton -row 10 -column 1 -sticky nsw -columnspan 2
 
  grid $MapToolsSigmaLabel -row 11 -column 0 -sticky nswe  
  grid $MapToolsSigmaButton -row 11 -column 1 -sticky nsw -columnspan 2
  
  grid $MapToolsSegLabel -row 12 -column 0 -sticky nswe 
  grid $MapToolsSegLabelAmt -row 12 -column 1 -sticky nswe 
  grid $MapToolsSegEntryAmt -row 12 -column 2 -sticky nswe 
  grid $MapToolsSegButton -row 12 -column 3 -sticky nsw -columnspan 2
  
  #Binary Ops
  set MapToolsBinaryFrame [ttk::labelframe $w.hlf.n.f5.main.binaryframe -labelanchor nw]
  
  set MapToolsMol2 [ttk::label $w.hlf.n.f5.main.binaryframe.mol2label -text "Mol2 ID (Density Map):"]
  set MapToolsMol2MenuButton [ttk::menubutton $w.hlf.n.f5.main.binaryframe.mol2menubutton -textvar MDFFGUI::gui::MapToolsMol2MenuText -menu $w.hlf.n.f5.main.binaryframe.mol2menubutton.mol2menu]
  set CurrentMapToolsMol2Menu [menu $w.hlf.n.f5.main.binaryframe.mol2menubutton.mol2menu -tearoff no]
  
  #set_mol_text 
  fill_mapmol2_menu $CurrentMapToolsMol2Menu
  trace add variable ::vmd_molecule write "MDFFGUI::gui::fill_mapmol2_menu $CurrentMapToolsMol2Menu"
  trace add variable MDFFGUI::maptoolssettings::MapToolsMol2ID write MDFFGUI::gui::set_mapmol2_text
  set MDFFGUI::maptoolssettings::MapToolsMol2ID $nullMolString
  
  set MapToolsVol2 [ttk::label $w.hlf.n.f5.main.binaryframe.vollabel -text "Vol ID (Density Map):"]
  set MapToolsVol2MenuButton [ttk::menubutton $w.hlf.n.f5.main.binaryframe.volmenubutton -textvar MDFFGUI::gui::MapToolsVol2MenuText -menu $w.hlf.n.f5.main.binaryframe.volmenubutton.volmenu]
  set CurrentMapToolsVol2Menu [menu $w.hlf.n.f5.main.binaryframe.volmenubutton.volmenu -tearoff no]
  
  #set_vol2_text 
  fill_mapvol2_menu $CurrentMapToolsVol2Menu
  #nuclear option trace to capture volume data sets added to an existing molecule.
  #trace add variable ::vmd_logfile write "MDFFGUI::gui::fill_mapvol2_menu $CurrentMapToolsVol2Menu"
  trace add variable MDFFGUI::maptoolssettings::MapToolsVol2ID write MDFFGUI::gui::set_mapvol2_text
  trace add variable MDFFGUI::maptoolssettings::MapToolsMol2ID write "MDFFGUI::gui::fill_mapvol2_menu $CurrentMapToolsVol2Menu"
  trace add variable MDFFGUI::maptoolssettings::MapToolsMol2ID write MDFFGUI::gui::revert_vol2id
  
  set MDFFGUI::maptoolssettings::MapToolsVol2ID $nullVolString
  
  set MapToolsSelectMap2 [ttk::button $w.hlf.n.f5.main.binaryframe.mapbutton -text "Browse" -command {MDFFGUI::gui::get_map2} ]
  
 #Union or Intersection
  set MapToolsIntersection [ttk::radiobutton $w.hlf.n.f5.main.binaryframe.intersection -variable MDFFGUI::maptoolssettings::UnionOrIntersection -value "" -text "Intersection"]
  set MapToolsUnion [ttk::radiobutton $w.hlf.n.f5.main.binaryframe.union -variable MDFFGUI::maptoolssettings::UnionOrIntersection -value "-union" -text "Union"]
 
 #Interpolation or No
  set MapToolsInterpolation [ttk::radiobutton $w.hlf.n.f5.main.binaryframe.interpolation -variable MDFFGUI::maptoolssettings::InterpolationOrNo -value "" -text "Interpolate"]
  set MapToolsNoInterp [ttk::radiobutton $w.hlf.n.f5.main.binaryframe.nointerp -variable MDFFGUI::maptoolssettings::InterpolationOrNo -value "-nointerp" -text "No Interpolation"]
  
  #Add
  set MapToolsAddButton [ttk::button $w.hlf.n.f5.main.binaryframe.addbutton -text "Add Maps" -command { voltool add -mol1 $MDFFGUI::maptoolssettings::MapToolsMolID -vol1 $MDFFGUI::maptoolssettings::MapToolsVolID -mol2 $MDFFGUI::maptoolssettings::MapToolsMol2ID -vol2 $MDFFGUI::maptoolssettings::MapToolsVol2ID $MDFFGUI::maptoolssettings::UnionOrIntersection $MDFFGUI::maptoolssettings::InterpolationOrNo } ]
  
  #Diff
  set MapToolsDiffButton [ttk::button $w.hlf.n.f5.main.binaryframe.diffbutton -text "Subtract Maps" -command { voltool diff -mol1 $MDFFGUI::maptoolssettings::MapToolsMolID -vol1 $MDFFGUI::maptoolssettings::MapToolsVolID -mol2 $MDFFGUI::maptoolssettings::MapToolsMol2ID -vol2 $MDFFGUI::maptoolssettings::MapToolsVol2ID $MDFFGUI::maptoolssettings::UnionOrIntersection $MDFFGUI::maptoolssettings::InterpolationOrNo } ]
  
  #Mult
  set MapToolsMultButton [ttk::button $w.hlf.n.f5.main.binaryframe.multbutton -text "Multiply Maps" -command { voltool mult -mol1 $MDFFGUI::maptoolssettings::MapToolsMolID -vol1 $MDFFGUI::maptoolssettings::MapToolsVolID -mol2 $MDFFGUI::maptoolssettings::MapToolsMol2ID -vol2 $MDFFGUI::maptoolssettings::MapToolsVol2ID $MDFFGUI::maptoolssettings::UnionOrIntersection $MDFFGUI::maptoolssettings::InterpolationOrNo } ]
  
  #Avg
  set MapToolsAvgButton [ttk::button $w.hlf.n.f5.main.binaryframe.avgbutton -text "Average Maps" -command { voltool avg -mol1 $MDFFGUI::maptoolssettings::MapToolsMolID -vol1 $MDFFGUI::maptoolssettings::MapToolsVolID -mol2 $MDFFGUI::maptoolssettings::MapToolsMol2ID -vol2 $MDFFGUI::maptoolssettings::MapToolsVol2ID $MDFFGUI::maptoolssettings::UnionOrIntersection $MDFFGUI::maptoolssettings::InterpolationOrNo } ]
  
  set ShowMapToolsBinary [ttk::label $w.hlf.n.f5.main.showbinary -text "$rightPoint Multiple Maps Ops..." -anchor w]
  set HideMapToolsBinary [ttk::label $w.hlf.n.f5.main.binaryframe.hidebinary -text "$downPoint Multiple Maps Ops" -anchor w]

  $MapToolsBinaryFrame configure -labelwidget $HideMapToolsBinary
  bind $HideMapToolsBinary <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.binaryframe
      grid .mdffgui.hlf.n.f5.main.showbinary
      MDFFGUI::gui::resizeToActiveTab
  }
  bind $ShowMapToolsBinary <Button-1> {
      grid remove .mdffgui.hlf.n.f5.main.showbinary
      grid .mdffgui.hlf.n.f5.main.binaryframe -row 5 -column 0 -sticky nsw -pady 5 -columnspan 10
 #     grid columnconfigure .mdffgui.hlf.n.f5.main.binaryframe 1 -weight 1
      MDFFGUI::gui::resizeToActiveTab
  }
  
  grid $ShowMapToolsBinary -row 5 -column 0 -sticky nswe
 
   
  grid $MapToolsMol2 -row 0 -column 0 -sticky nswe
  grid $MapToolsMol2MenuButton -row 0 -column 1 -sticky nswe
  grid $MapToolsVol2 -row 0 -column 2 -sticky nswe
  grid $MapToolsVol2MenuButton -row 0 -column 3 -sticky nswe
  grid $MapToolsSelectMap2 -row 0 -column 4 -sticky nswe
  grid $MapToolsIntersection -row 1 -column 0 -sticky nswe 
  grid $MapToolsUnion -row 1 -column 1 -sticky nswe 
  grid $MapToolsInterpolation -row 2 -column 0 -sticky nswe 
  grid $MapToolsNoInterp -row 2 -column 1 -sticky nswe 
  grid $MapToolsAddButton -row 3 -column 0 -sticky nswe
  grid $MapToolsDiffButton -row 3 -column 1 -sticky nswe
  grid $MapToolsMultButton -row 3 -column 2 -sticky nswe
  grid $MapToolsAvgButton -row 3 -column 3 -sticky nswe
  #Basic MDFF analysis
  #set NBTab6 [ttk::frame $w.hlf.n.f6];
  #$Notebook add $NBTab6 -text "Analysis" 
  
  
  bind $Notebook <<NotebookTabChanged>> {MDFFGUI::gui::resizeToActiveTab}
}

proc MDFFGUI::gui::mask_toggle {args} {
  variable MaskResInput
  variable MaskCutoffInput

  if {$MDFFGUI::settings::Mask && $MaskResInput != "" && $MaskCutoffInput != "" && [winfo exists ".add_density_xmdff"] } {
    $MaskResInput configure -state normal  
    $MaskCutoffInput configure -state normal  
  } elseif {$MaskResInput != "" && $MaskCutoffInput != "" && [winfo exists ".add_density_xmdff"]} {
    $MaskResInput configure -state disabled 
    $MaskCutoffInput configure -state disabled  
  }


}

proc MDFFGUI::gui::add_density {} {
  variable wad

  if { [winfo exists ".add_density"] } {
    wm deiconify ".add_density"
    return
  }

  set wad [toplevel ".add_density"]
  wm title $wad "Add density map"
  wm resizable $wad no no

  set AddDensityHLF [ttk::frame $wad.hlf]
  grid $wad.hlf -column 0 -row 0 -sticky nsew
  # allow hlf to resize with window
  grid columnconfigure $wad.hlf 0 -weight 1
  grid rowconfigure $wad.hlf 1 -weight 1

  
  set ADCurrentDXFile [ttk::label $wad.hlf.dxlabel -text "Density Map:"]
  set ADCurrentDXEntry [ttk::entry $wad.hlf.dxentry -textvariable MDFFGUI::settings::CurrentDXPath -width 20]
  set ADSelectDX [ttk::button $wad.hlf.dx -text "Browse" -command "MDFFGUI::gui::get_density $ADCurrentDXEntry" ]

  set ADGridpdbLabel [ttk::label $wad.hlf.seltextlabel -text "Gridpdb selection text:"]
  set ADGridpdbInput [ttk::entry $wad.hlf.seltextentry -textvariable MDFFGUI::settings::GridPDBSelText -width 20]

  set ADGScaleLabel [ttk::label $wad.hlf.gscalelabel -text "Grid Scaling Factor:"]
  set ADGScaleInput [ttk::entry $wad.hlf.gscaleentry -textvariable MDFFGUI::settings::GScale]
  
  set ADThresholdLabel [ttk::label $wad.hlf.thresholdlabel -text "Threshold:"]
  set ADThresholdInput [ttk::entry $wad.hlf.thresholdentry -textvariable MDFFGUI::settings::GriddxThreshold]

  set ADLoadDensity [ttk::button $wad.hlf.loaddensity -text "Add Density" -command {MDFFGUI::gui::load_density}]

  grid $ADCurrentDXFile -row 0 -column 0 -sticky nsew
  grid $ADCurrentDXEntry -row 0 -column 1 -sticky nsew
  grid $ADSelectDX -row 0 -column 2 -sticky nsew

  grid $ADGridpdbLabel -row 1 -column 0 -sticky nsew
  grid $ADGridpdbInput -row 1 -column 1 -sticky nsew

  grid $ADGScaleLabel -row 2 -column 0 -sticky nsew
  grid $ADGScaleInput -row 2 -column 1 -sticky nsew
  
  grid $ADThresholdLabel -row 3 -column 0 -sticky nsew
  grid $ADThresholdInput -row 3 -column 1 -sticky nsew

  grid $ADLoadDensity -row 4 -column 0 -sticky nsew -columnspan 3
}

proc MDFFGUI::gui::add_density_xmdff {} {
  variable xwad
  variable CButtonMask
  variable MaskResLabel
  variable MaskResInput
  variable MaskCutoffLabel
  variable MaskCutoffInput

  if { [winfo exists ".add_density_xmdff"] } {
    wm deiconify ".add_density_xmdff"
    return
  }

  set xwad [toplevel ".add_density_xmdff"]
  wm title $xwad "Add density map"
  wm resizable $xwad no no

  set xAddDensityHLF [ttk::frame $xwad.hlf]
  grid $xwad.hlf -column 0 -row 0 -sticky nsew
  # allow hlf to resize with window
  grid columnconfigure $xwad.hlf 0 -weight 1
  grid rowconfigure $xwad.hlf 1 -weight 1

  
  set xADCurrentDXFile [ttk::label $xwad.hlf.dxlabel -text "Reflections file:"]
  set xADCurrentDXEntry [ttk::entry $xwad.hlf.dxentry -textvariable MDFFGUI::settings::CurrentRefsPath -width 20]
  set xADSelectDX [ttk::button $xwad.hlf.dx -text "Browse" -command "MDFFGUI::gui::get_refs $xADCurrentDXEntry" ]

  set xADGridpdbLabel [ttk::label $xwad.hlf.seltextlabel -text "Gridpdb selection text:"]
  set xADGridpdbInput [ttk::entry $xwad.hlf.seltextentry -textvariable MDFFGUI::settings::GridPDBSelText -width 20]

  set xADGScaleLabel [ttk::label $xwad.hlf.gscalelabel -text "Grid Scaling Factor:"]
  set xADGScaleInput [ttk::entry $xwad.hlf.gscaleentry -textvariable MDFFGUI::settings::GScale]

  set CButtonMask [ttk::checkbutton $xwad.hlf.mask -text "Mask Density" -variable MDFFGUI::settings::Mask -command {MDFFGUI::gui::mask_toggle}]
  
  set MaskResLabel [ttk::label $xwad.hlf.maskreslabel -text "Mask Resolution (A):"]
  set MaskResInput [ttk::entry $xwad.hlf.maskresentry -textvariable MDFFGUI::settings::MaskRes -width 40 -state disabled]
  
  set MaskCutoffLabel [ttk::label $xwad.hlf.maskcutofflabel -text "Mask Cutoff (A):"]
  set MaskCutoffInput [ttk::entry $xwad.hlf.maskcutoffentry -textvariable MDFFGUI::settings::MaskCutoff -width 40 -state disabled]
  
  set CButtonAverage [ttk::checkbutton $xwad.hlf.average -text "Average Map" -variable MDFFGUI::settings::AverageMap]
  
  set BSharpFrame [ttk::frame $xwad.hlf.bsharpframe]
  set BSharpNone [ttk::radiobutton  $xwad.hlf.bsharpframe.bsharpnone -variable MDFFGUI::settings::BSharp -value "none" -text "None"]
  set BSharpSharp [ttk::radiobutton  $xwad.hlf.bsharpframe.bsharpsharp -variable MDFFGUI::settings::BSharp -value "sharpen" -text "Sharpen"]
  set BSharpSmooth [ttk::radiobutton  $xwad.hlf.bsharpframe.bsharpsmooth -variable MDFFGUI::settings::BSharp -value "smooth" -text "Smooth"]
  set BSharpLabel [ttk::label $xwad.hlf.bsharpframe.bsharplabel -text "Sharpening/Smoothing Factor \n(0 for automatic, Sharpening only):"]
  set BSharpInput [ttk::entry $xwad.hlf.bsharpframe.bsharpentry -textvariable MDFFGUI::settings::BSharpFactor]
  
  set xADLoadDensity [ttk::button $xwad.hlf.loaddensity -text "Add Density" -command {MDFFGUI::gui::load_density_xmdff}]

  grid $xADCurrentDXFile -row 0 -column 0 -sticky nsew
  grid $xADCurrentDXEntry -row 0 -column 1 -sticky nsew
  grid $xADSelectDX -row 0 -column 2 -sticky nsew

  grid $xADGridpdbLabel -row 1 -column 0 -sticky nsew
  grid $xADGridpdbInput -row 1 -column 1 -sticky nsew

  grid $xADGScaleLabel -row 2 -column 0 -sticky nsew
  grid $xADGScaleInput -row 2 -column 1 -sticky nsew

  grid $CButtonMask -row 3 -column 0 -sticky nsew
  grid $MaskResLabel -row 4 -column 0 -sticky nw
  grid $MaskResInput -row 4 -column 1 -sticky nsew
 
  grid $MaskCutoffLabel -row 5 -column 0 -sticky nw
  grid $MaskCutoffInput -row 5 -column 1 -sticky nsew

  grid $CButtonAverage -row 6 -column 0 -sticky nsew
  
  grid $BSharpFrame -row 7 -column 0 -sticky nsew -columnspan 3
  grid $BSharpNone -row 0 -column 0 -sticky nsew
  grid $BSharpSharp -row 0 -column 1 -sticky nsew
  grid $BSharpSmooth -row 0 -column 2 -sticky nsew
  grid $BSharpLabel -row 0 -column 3 -sticky nsew
  grid $BSharpInput -row 0 -column 4 -sticky nsew
  
  grid $xADLoadDensity -row 8 -column 0 -sticky nsew -columnspan 3
}

proc MDFFGUI::gui::add_paramfile {parent} {
  set pathname [tk_getOpenFile -defaultextension "CHARMM Parameter Files" -filetypes {{inp {.inp}} {str {.str}} {prm {.prm}} {all {*}}} -parent $parent -initialdir $::env(CHARMMPARDIR) -multiple 1]
  
  if {$pathname != ""} {
    foreach path $pathname {
      set MDFFGUI::settings::ParameterList [lappend MDFFGUI::settings::ParameterList $path]
    }
  }
}

proc MDFFGUI::gui::imd_server {} {  
  variable CurrentServer
  upvar MDFFGUI::servers::$CurrentServer currentserver
  set MDFFGUI::settings::IMDProcs "$currentserver(numprocs)"
}


proc MDFFGUI::gui::load_density_list {} {
  variable DensityView
  variable MapID

  if {[llength $MDFFGUI::settings::DensityList] !=  [llength $MDFFGUI::settings::GridPDBSelList] && \
     [llength $MDFFGUI::settings::DensityList] !=  [llength $MDFFGUI::settings::GScaleList]} {
        tk_messageBox -message "Warning: Uneven number of density parameters!\nPlease fix the settings file."
        return 0
     } else {
       foreach file $MDFFGUI::settings::DensityList gscale $MDFFGUI::settings::GScaleList seltext $MDFFGUI::settings::GridPDBSelList threshold $MDFFGUI::settings::ThresholdList {
          set MDFFGUI::settings::CurrentDXPath $file 
          set MDFFGUI::settings::GridPDBSelText $seltext
          set MDFFGUI::settings::GScale $gscale
          set MDFFGUI::settings::GriddxThreshold $threshold
          load_density
          #set MapID [mol new $file]
          #$DensityView insert {} 0 -id $MapID -values [list $file $seltext $gscale]
        }      
     }
  if {[llength $MDFFGUI::settings::RefsList] !=  [llength $MDFFGUI::settings::xMDFFGridPDBSelList] && \
     [llength $MDFFGUI::settings::RefsList] !=  [llength $MDFFGUI::settings::xMDFFGScaleList]} {
        tk_messageBox -message "Warning: Uneven number of density parameters!\nPlease fix the settings file."
        return 0
     } else {
       foreach file $MDFFGUI::settings::RefsList gscale $MDFFGUI::settings::xMDFFGScaleList seltext $MDFFGUI::settings::xMDFFGridPDBSelList mask $MDFFGUI::settings::xMDFFMaskList maskres $MDFFGUI::settings::xMDFFMaskResList maskcutoff $MDFFGUI::settings::xMDFFMaskCutoffList average $MDFFGUI::settings::xMDFFAverageMapList bsharp $MDFFGUI::settings::xMDFFBSharpList bsharpfactor $MDFFGUI::settings::xMDFFBSharpFactorList { 
          set MDFFGUI::settings::CurrentRefsPath $file 
          set MDFFGUI::settings::GridPDBSelText $seltext
          set MDFFGUI::settings::GScale $gscale
          set MDFFGUI::settings::Mask $mask
          set MDFFGUI::settings::MaskRes $maskres
          set MDFFGUI::settings::MaskCutoff $maskcutoff
          set MDFFGUI::settings::AverageMap $average
          set MDFFGUI::settings::BSharp $bsharp
          set MDFFGUI::settings::BSharpFactor $bsharpfactor
          load_density_xmdff
        }      
     }

     
  return 1
}

proc MDFFGUI::gui::save_density_list {} {
  variable DensityView
  variable xMDFFDensityView
   
  set MDFFGUI::settings::RefsList ""
  set MDFFGUI::settings::DensityList ""
  set MDFFGUI::settings::GridPDBSelList ""
  set MDFFGUI::settings::GScaleList ""
  set MDFFGUI::settings::xMDFFGridPDBSelList ""
  set MDFFGUI::settings::xMDFFGScaleList ""
  set MDFFGUI::settings::xMDFFMaskList ""
  set MDFFGUI::settings::xMDFFMaskResList ""
  set MDFFGUI::settings::xMDFFMaskCutoffList ""
  set MDFFGUI::settings::xMDFFAverageMapList ""
  set MDFFGUI::settings::xMDFFBSharpList ""
  set MDFFGUI::settings::xMDFFBSharpFactorList ""
 
  foreach density [$DensityView children {}] {
    set vals [$DensityView item $density -values]
    lappend MDFFGUI::settings::DensityList [lindex $vals 0] 
    lappend MDFFGUI::settings::GridPDBSelList [lindex $vals 1]
    lappend MDFFGUI::settings::GScaleList [lindex $vals 2]
    lappend MDFFGUI::settings::ThresholdList [lindex $vals 3]
  }
  
  foreach density [$xMDFFDensityView children {}] {
    set vals [$xMDFFDensityView item $density -values]
    lappend MDFFGUI::settings::RefsList [lindex $vals 0] 
    lappend MDFFGUI::settings::xMDFFGridPDBSelList [lindex $vals 1]
    lappend MDFFGUI::settings::xMDFFGScaleList [lindex $vals 2]
    lappend MDFFGUI::settings::xMDFFMaskList [lindex $vals 3]
    lappend MDFFGUI::settings::xMDFFMaskResList [lindex $vals 4]
    lappend MDFFGUI::settings::xMDFFMaskCutoffList [lindex $vals 5]
    lappend MDFFGUI::settings::xMDFFAverageMapList [lindex $vals 6]
    lappend MDFFGUI::settings::xMDFFBSharpList [lindex $vals 7]
    lappend MDFFGUI::settings::xMDFFBSharpFactorList [lindex $vals 8]
  }
  

}

proc MDFFGUI::gui::save_map {} {
  set pathname [tk_getSaveFile -defaultextension ".dx" -filetypes {{dx {.dx}} } -parent $MDFFGUI::gui::w -initialdir $MDFFGUI::settings::CurrentDir]
  voltool write -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID -o $pathname
}

proc MDFFGUI::gui::save_settings {parent} {
  
  
  set pathname [tk_getSaveFile -defaultextension ".tcl" -filetypes {{tcl {.tcl}} {all {*}}} -parent $parent -initialdir $MDFFGUI::settings::CurrentDir]

  if {$pathname != "" } {
    save_struct
    save_density_list
    set f [open $pathname "w"] 
    foreach setting [lsort [info vars ::MDFFGUI::settings::*]] {
      puts $f "set $setting [list [set $setting]]" 
    }
    foreach server [lsort [info vars ::MDFFGUI::servers::*]] {
      upvar $server currentserver
      puts $f "MDFFGUI::gui::add_server \"[namespace tail $server]\" \{"
      foreach e [array names $server] {
        puts $f "  $e [list $currentserver($e)]"  
      }
      puts $f "\}"
    }
    close $f
  }
}

proc MDFFGUI::gui::load_settings {parent} {  
  set pathname [tk_getOpenFile -defaultextension ".tcl" -filetypes {{tcl {.tcl}} {all {*}}} -parent $parent -initialdir $MDFFGUI::settings::CurrentDir]

  if {$pathname != "" } {

    set MDFFGUI::settings::RefsList ""
    set MDFFGUI::settings::DensityList ""
    set MDFFGUI::settings::GridPDBSelList ""
    set MDFFGUI::settings::GScaleList ""
    set MDFFGUI::settings::xMDFFGridPDBSelList ""
    set MDFFGUI::settings::xMDFFGScaleList ""
    set MDFFGUI::settings::xMDFFMaskList ""
    set MDFFGUI::settings::xMDFFMaskResList ""
    set MDFFGUI::settings::xMDFFMaskCutoffList ""
    set MDFFGUI::settings::xMDFFAverageMapList ""
    set MDFFGUI::settings::xMDFFBSharpList ""
    set MDFFGUI::settings::xMDFFBSharpFactorList ""
    mol delete $MDFFGUI::settings::MolID

    source $pathname
    if {$::MDFFGUI::settings::CrossCorrelation} {MDFFGUI::gui::cc_toggle}
    if {$::MDFFGUI::settings::xMDFF} {MDFFGUI::gui::xmdff_toggle}
   # if {$::MDFFGUI::settings::Mask} {MDFFGUI::gui::mask_toggle}
    if {$::MDFFGUI::settings::IMD} {MDFFGUI::gui::imd_toggle}
    if {$::MDFFGUI::settings::UseCCThresh} {MDFFGUI::gui::thresh_toggle}
    foreach density [$MDFFGUI::gui::DensityView children {}] {
      remove_densityview $MDFFGUI::gui::DensityView $density
    }
    foreach density [$MDFFGUI::gui::xMDFFDensityView children {}] {
      remove_densityview $MDFFGUI::gui::xMDFFDensityView $density
    }

    load_density_list
    #loads structure and density...nice for me testing, but does a user want this?
#    if {$::MDFFGUI::settings::CurrentDXPath != ""} {
#      set MDFFGUI::gui::MapID [mol new $MDFFGUI::settings::CurrentDXPath]
#      set MDFFGUI::gui::DensityID $MDFFGUI::gui::MapID  
#    }
    if {$::MDFFGUI::settings::CurrentPSFPath != "" && $::MDFFGUI::settings::CurrentPDBPath !=""} {MDFFGUI::gui::load_struct}
  }
}

proc MDFFGUI::gui::resizeToActiveTab {} {
  variable Notebook
  variable HLF

  update idletasks

#  set dimW [winfo reqwidth [$Notebook select]]
#  set dimMinW [winfo reqwidth $Notebook]

  set dimH [winfo reqheight [$Notebook select]]
  set dimH [expr $dimH + 75]
  
#  if {$dimW > $dimMinW} {
#    set useW $dimW 
#  } else {
#    set useW $dimMinW
#  }

  #other W code for variable widths, but that doesn't look as good on MACOS so go with
  #static width for now
  set useW 800
  
  wm geometry $MDFFGUI::gui::w [format "%ix%i" $useW $dimH] 

}

proc MDFFGUI::gui::xmdff_toggle {args} {
  variable CurrentRefsEntry
  variable RefStepsInput
  variable CurrentCrystEntry
  variable GetRefs
  variable GetCryst
  variable CButtonMask
  variable CButtonBFS
  variable xMDFFDensityViewAdd
  variable xMDFFSelInput
    
  if {$MDFFGUI::settings::xMDFF} {
    $CurrentRefsEntry configure -state normal
    $RefStepsInput configure -state normal
    $CurrentCrystEntry configure -state normal
    $GetRefs configure -state normal
    $GetCryst configure -state normal
#    $CButtonMask configure -state normal
    $CButtonBFS configure -state normal
    $xMDFFDensityViewAdd configure -state normal
    $xMDFFSelInput configure -state normal
  } else {
    $CurrentRefsEntry configure -state disabled
    $RefStepsInput configure -state disabled
    $CurrentCrystEntry configure -state disabled
    $GetRefs configure -state disabled
    $GetCryst configure -state disabled
#    $CButtonMask configure -state disabled
    $CButtonBFS configure -state disabled
    $xMDFFDensityViewAdd configure -state disabled
    $xMDFFSelInput configure -state disabled
  }

}

proc MDFFGUI::gui::remdff_toggle {args} {
  variable REMDFFNumReplicasInput
  variable CButtonREMDFFAutoBlur
    
  if {$MDFFGUI::settings::REMDFF} {
    $REMDFFNumReplicasInput configure -state normal
    $CButtonREMDFFAutoBlur configure -state normal
  } else {
    $REMDFFNumReplicasInput configure -state disabled
    $CButtonREMDFFAutoBlur configure -state disabled
  }

}

proc MDFFGUI::gui::thresh_toggle {args} {
  variable CCThreshInput

  if {$MDFFGUI::settings::UseCCThresh} {
    $CCThreshInput configure -state normal
  } else {
    $CCThreshInput configure -state disabled
  }

}


proc MDFFGUI::gui::cc_toggle {args} {
  variable CCResInput
  variable CCSelInput
  variable CCMolInput

  if {$MDFFGUI::settings::CrossCorrelation} {
    $CCSelInput configure -state normal
    $CCResInput configure -state normal
    $CCMolInput configure -state normal
  } else {
    $CCSelInput configure -state disabled
    $CCResInput configure -state disabled
    $CCMolInput configure -state disabled
  }

}

proc MDFFGUI::gui::imd_kill {} {
  imd kill
  $MDFFGUI::gui::IMDStatusLabel configure -text "IMD Status: NAMD job terminated."
}

proc MDFFGUI::gui::imd_pause {} {

  imd pause toggle
  #$IMDStatusLabel configure -text "IMD Status: Paused."
}

proc MDFFGUI::gui::imd_toggle {args} {
  variable CButtonIMDWait
  variable CButtonIMDIgnore
  variable IMDPortInput
  variable IMDFreqInput
  variable IMDKeepInput
  variable IMDServer

  if {$MDFFGUI::settings::IMD} {
    $CButtonIMDWait configure -state normal
    $CButtonIMDIgnore configure -state normal
    $IMDPortInput configure -state normal
    $IMDFreqInput configure -state normal
    $IMDKeepInput configure -state normal
    $IMDServer configure -state normal
    $MDFFGUI::gui::IMDSubmit configure -state normal
    $MDFFGUI::gui::IMDConnect configure -state normal
    $MDFFGUI::gui::IMDPause configure -state normal
    $MDFFGUI::gui::IMDFinish configure -state normal
    $MDFFGUI::gui::IMDProcInput configure -state normal
  } else {
    $CButtonIMDWait configure -state disabled
    $CButtonIMDIgnore configure -state disabled
    $IMDPortInput configure -state disabled
    $IMDFreqInput configure -state disabled
    $IMDKeepInput configure -state disabled
    $IMDServer configure -state disabled
    $MDFFGUI::gui::IMDSubmit configure -state disabled
    $MDFFGUI::gui::IMDConnect configure -state disabled
    $MDFFGUI::gui::IMDPause configure -state disabled
    $MDFFGUI::gui::IMDFinish configure -state disabled
    $MDFFGUI::gui::IMDProcInput configure -state disabled
  }
}

proc MDFFGUI::gui::select_imd_host {} {
  variable IMDInpOrList
  variable IMDHostInput
  variable IMDServer 
  
  if {$IMDInpOrList == "input"} {
    $IMDServer configure -state disabled
    $IMDHostInput configure -state normal
  } elseif {$IMDInpOrList == "list"} {
    $IMDHostInput configure -state disabled
    $IMDServer configure -state normal
  }

}

proc MDFFGUI::gui::save_struct {} {
#Note: This section would be closer to how autoimd works if FixedPDBSelText is set
#to be: within 8 of ($IMDSelText) and (not ($IMDSelText)). Instead, for now, I am not truncating the
#structure for IMD purposes because this way the output dcd's give a complete history of the changes.
 # if $MDFFGUI::settings::IMD {
 #   set sel_mol [atomselect $MDFFGUI::settings::MolID "($MDFFGUI::settings::IMDSelText) or ($MDFFGUI::settings::FixedPDBSelText)"]
 # } else {}

  file mkdir $MDFFGUI::settings::CurrentDir
  if {[lsearch [molinfo list] $MDFFGUI::settings::MolID] != -1} {
    set sel_mol [atomselect $MDFFGUI::settings::MolID all]
    $sel_mol frame now
    set struct_name "$MDFFGUI::settings::SimulationName-step$MDFFGUI::settings::SimulationStep"
    
    set MDFFGUI::settings::CurrentPDBPath [file join $MDFFGUI::settings::CurrentDir "$struct_name.pdb"]
    set MDFFGUI::settings::CurrentPSFPath [file join $MDFFGUI::settings::CurrentDir "$struct_name.psf"]
    
    $sel_mol writepdb $MDFFGUI::settings::CurrentPDBPath
    $sel_mol writepsf $MDFFGUI::settings::CurrentPSFPath
  }
}

proc MDFFGUI::gui::mdff_setup {} {
  
  variable DensityView  
  variable xMDFFDensityView
  variable Extrabonds
  
  #variable MDFFGUI::settings::CurrentDir 
  set CM [molinfo $MDFFGUI::settings::MolID get center_matrix]
  set RM [molinfo $MDFFGUI::settings::MolID get rotate_matrix]
  set SM [molinfo $MDFFGUI::settings::MolID get scale_matrix]
  set GM [molinfo $MDFFGUI::settings::MolID get global_matrix]
  
  file mkdir $MDFFGUI::settings::CurrentDir
  save_struct 
  variable CrystOption
  
 
  #variable MDFFGUI::settings::CurrentDir 
   
  #file rename -force $MDFFGUI::settings::CurrentDXPath $MDFFGUI::settings::CurrentDir 
#  if {$Extrabonds != 0} {
#    foreach extrabond $Extrabonds {
#      file rename -force $extrabond $MDFFGUI::settings::CurrentDir 
#    }
#  }
#  file rename -force "gridpdb.pdb" $MDFFGUI::settings::CurrentDir 
#  file rename -force "griddx.dx" $MDFFGUI::settings::CurrentDir 
  
  if {$MDFFGUI::settings::CurrentPSFPath == ""} {
    tk_messageBox -message "Please provide a PSF file."
  } elseif {$MDFFGUI::settings::CurrentPDBPath == ""} {
    tk_messageBox -message "Please provide a PDB file."
  } elseif {!$MDFFGUI::settings::REMDFF && !$MDFFGUI::settings::xMDFF && [llength [$MDFFGUI::gui::DensityView children {}]] == 0 && !$MDFFGUI::settings::GridOff} {
    tk_messageBox -message "Please provide a density map file or turn grid forces off in simulation parameters."
  } elseif {$MDFFGUI::settings::MolID == ""} {
    tk_messageBox -message "Please select a molecule containing a PSF and PDB."
  } elseif {$MDFFGUI::settings::xMDFF && [llength [$MDFFGUI::gui::xMDFFDensityView children {}]] == 0} {
    tk_messageBox -message "Please provide a reflections data file or turn off xMDFF."
  } elseif {($MDFFGUI::settings::REMDFF && $MDFFGUI::settings::xMDFF) || ($MDFFGUI::settings::REMDFF && $MDFFGUI::settings::xMDFF)} {
    tk_messageBox -message "Replica exchange not currently compatible with IMD or xMDFF."
  } elseif {$MDFFGUI::settings::REMDFF && [llength [$DensityView children {}]] > 1} {
    tk_messageBox -message "Replica exchange not currently compatible with multiple input density maps."
  } else {
    set refs ""
    set refsteps ""
    set xmdff ""
    set bfs ""
    set mask ""
    set maskres ""
    set maskcutoff ""
    set averagemap ""
    set bsharp ""
    set remdff ""
    set autosmooth ""
    set replicasnum ""

  #  file mkdir $MDFFGUI::settings::CurrentDir 
    
    if {[llength $Extrabonds] == 0} {set Extrabonds 0}
    if {$MDFFGUI::settings::xMDFF} {
      set xmdff "--xmdff"
      set refsteps "-refsteps $MDFFGUI::settings::RefSteps"
      set refs "-refs"
      #set refs "-refs [file tail $MDFFGUI::settings::CurrentRefsPath]"
      #set maskres "-mask_res $MDFFGUI::settings::MaskRes"
      #set maskcutoff "-mask_cutoff $MDFFGUI::settings::MaskCutoff"
      #file copy -force $MDFFGUI::settings::CurrentRefsPath $MDFFGUI::settings::CurrentDir 
      if {$CrystOption != ""} {
        file copy -force $MDFFGUI::settings::CurrentCrystPath $MDFFGUI::settings::CurrentDir 
      }
      if {$MDFFGUI::settings::BFS} {
        set bfs "--bfs"
      }
    #  if {$MDFFGUI::settings::Mask} {
    #   set mask "--mask"
    #  }
    }

    set fixedpdb ""
    set fixedcol ""

    if {$MDFFGUI::settings::FixedPDBSelText != "none"} {
      set fixedpdb "-fixpdb fixed.pdb"
      set fixedcol "-fixcol O"
#    file rename -force "fixed.pdb" MDFFGUI::$MDFFGUI::settings::CurrentDir 
    }
     
    set parlist ""
    foreach par $MDFFGUI::settings::ParameterList {
      if {[file dirname $par] != $MDFFGUI::settings::CurrentDir} {
        file copy -force $par $MDFFGUI::settings::CurrentDir
        set parlist [lappend parlist [file tail $par]]
      }
    }
    set imd ""
    set imdport ""
    set imdfreq ""
    set imdwait ""
    set imdignore ""
    if {$MDFFGUI::settings::IMD} {
      set imd "--imd"
      set imdport "-imdport $MDFFGUI::settings::IMDPort"
      set imdfreq "-imdfreq $MDFFGUI::settings::IMDFreq"
      if {$MDFFGUI::settings::IMDWait} {
        set imdwait "--imdwait"
      }
      if {$MDFFGUI::settings::IMDIgnore} {
        set imdignore "--imdignore"
      }
    }
    
    set lite ""
    if {$MDFFGUI::settings::GridforceLite} {
      set lite "--lite"
    }
    
    set gridoff ""
    if {$MDFFGUI::settings::GridOff} {
      set gridoff "--gridoff"
    }
    

    
    MDFFGUI::gui::make_fixedpdb
    
    set maps 0
    set gscales ""
    set potentials ""
    set gridpdbs ""
    foreach density [$DensityView children {}] {
      set vals [$DensityView item $density -values]
      set file [lindex $vals 0] 
      set seltext [lindex $vals 1]
      set gscale [lindex $vals 2]
      lappend gscales $gscale  
      MDFFGUI::gui::make_gridpdb $maps $seltext
      lappend gridpdbs "gridpdb$maps.pdb"
      if {$MDFFGUI::settings::REMDFF && $MDFFGUI::settings::REMDFFAutoBlur} {
        # 0.dx is the original map with no smoothing, so do this outside the loop first
        set replicadxfile [file join $MDFFGUI::settings::CurrentDir "0.dx"]
        file delete -force $replicadxfile
        #don't use make_griddx since we aren't following the usual naming convention in the case of replicas
        mdff griddx -i $file -o $replicadxfile
        lappend potentials "0.dx"
        #for now each subsequent blur is only incremented by 1 sigma, but could
        #introduce other stepping later.
        for {set i 1} {$i < $MDFFGUI::settings::REMDFFNumReplicas} {incr i} {
          set replicadxfile [file join $MDFFGUI::settings::CurrentDir "$i.dx"]
          file delete -force $replicadxfile
          voltool smooth -sigma $i -i $file -o $replicadxfile
          #don't use make_griddx since we aren't following the usual naming convention in the case of replicas
          mdff griddx -i $replicadxfile -o $replicadxfile
          lappend potentials "$i.dx"
        }
      } else {
        MDFFGUI::gui::make_griddx $maps $file 
        lappend potentials "griddx$maps.dx"
      }
      incr maps
    }
    set reffiles ""
    set masks ""
    set masksres ""
    set maskscutoff ""
    set xmdffsel ""
    set averagemaps ""
    set bsharps ""
    if {$MDFFGUI::settings::xMDFF} {
      foreach density [$xMDFFDensityView children {}] {
        set vals [$xMDFFDensityView item $density -values]
        set file [lindex $vals 0] 
        set seltext [lindex $vals 1]
        set gscale [lindex $vals 2]
        set maski [lindex $vals 3]
        set maskresi [lindex $vals 4]
        set maskcutoffi [lindex $vals 5]
        set averagemapi [lindex $vals 6]
        set bsharpi [lindex $vals 7]
        set bsharpfactori [lindex $vals 8]

        file copy -force $file $MDFFGUI::settings::CurrentDir 
        lappend gscales $gscale  
        lappend reffiles [file tail $file]
        lappend masks $maski
        lappend masksres $maskresi
        lappend maskscutoff $maskcutoffi
        lappend averagemaps $averagemapi
        if {$bsharpi == "none"} {
          lappend bsharps 0
        } elseif {$bsharpi == "sharpen"} {
          if {$bsharpfactori == 0} {
            lappend bsharps "a"
          } else {  
           lappend bsharps [expr abs($bsharpfactori)]
          }
        } elseif {$bsharpi == "smooth"} {
           lappend bsharps [expr -abs($bsharpfactori)]
        }
        
        MDFFGUI::gui::make_gridpdb $maps $seltext
        lappend gridpdbs "gridpdb$maps.pdb"
        lappend potentials "griddx$maps.dx"
        incr maps
      }
      if {$masks != ""} {
        set mask "-mask"      
        set maskres "-mask_res"
        set maskcutoff "-mask_cutoff"
      }
      if {$averagemaps != ""} {
        set averagemap "-averagemap"
      }
      if {$bsharps != ""} {
        set bsharp "-bsharp"
      }
      set xmdffsel "-xmdffsel \"$MDFFGUI::settings::xMDFFSel\""  
    }

    if {$MDFFGUI::settings::REMDFF} {
      set remdff "--remdff"
      set replicasnum "-replicas $MDFFGUI::settings::REMDFFNumReplicas"
      set initialMapDir [file join $MDFFGUI::settings::CurrentDir "initialmaps/"]
      file mkdir $initialMapDir
      if {$MDFFGUI::settings::REMDFFAutoBlur} {
        #set this whe you make mdff_setup handle this instead. only makes
        #sense if mdff_setup can also auto-make potential files, which we might want to do.
        #set autosmooth "--autosmooth"
        for {set i 0} {$i < $MDFFGUI::settings::REMDFFNumReplicas} {incr i} {
          file mkdir [file join $MDFFGUI::settings::CurrentDir "output/$i"]
          file rename -force [file join $MDFFGUI::settings::CurrentDir [lindex $potentials $i]] [file join $initialMapDir "$i.dx"] 
          file delete -force [file join $MDFFGUI::settings::CurrentDir [lindex $potentials $i]]
          file link "[file join $MDFFGUI::settings::CurrentDir [lindex $potentials $i]]" [file join $initialMapDir "$i.dx"] 
        }
      }
    }

    MDFFGUI::gui::generate_xbonds

    mdff setup -o $MDFFGUI::settings::SimulationName -dir $MDFFGUI::settings::CurrentDir -psf [file tail $MDFFGUI::settings::CurrentPSFPath] \
    -pdb [file tail $MDFFGUI::settings::CurrentPDBPath] -griddx $potentials \
    -gridpdb $gridpdbs -extrab $Extrabonds -gscale $gscales -minsteps $MDFFGUI::settings::Minsteps -numsteps $MDFFGUI::settings::Numsteps \
    -minsteps $MDFFGUI::settings::Minsteps -numsteps $MDFFGUI::settings::Numsteps -temp $MDFFGUI::settings::Temperature -ftemp $MDFFGUI::settings::FTemperature $MDFFGUI::settings::PBCorGBIS \
    $xmdff $refs $reffiles {*}$refsteps {*}$CrystOption {*}$xmdffsel $bfs $mask $masks $maskres $masksres $bsharp $bsharps \
    $maskcutoff $maskscutoff $averagemap $averagemaps $imd {*}$imdport {*}$imdfreq $imdwait $imdignore \
    {*}$fixedpdb {*}$fixedcol $lite $gridoff -parfiles $parlist -step $MDFFGUI::settings::SimulationStep $remdff {*}$replicasnum 
    #put this option back in when mdff_setup handles this instead.
    #$autosmooth 
  
    #write QwikMD log file if GUI was run from QwikMD
#    if {$MDFFGUI::settings::QwikMDLogFile != ""} {
#      set qwikmdlogfile [open $MDFFGUI::settings::QwikMDLogFile "a"]

#    }
  }
  foreach mol [molinfo list] {
    molinfo $mol set center_matrix $CM 
    molinfo $mol set rotate_matrix $RM
    molinfo $mol set scale_matrix $SM
    molinfo $mol set global_matrix $GM
  }
}

proc MDFFGUI::gui::add_server {args} {
  variable ServerList
  variable CurrentServer
  variable IMDServer

  set servername "[lindex $args 0]"
  set opts [lindex $args 1]
  
  array set $servername $opts

#  lappend ServerList $servername
#   lappend ServerList [array get $servername]
 
  if {[info exists IMDServer] && [lsearch $ServerList $servername] == -1} {
    array set "::MDFFGUI::servers::$servername" $opts
    set ServerList [lappend ServerList $servername]
    $IMDServer configure -values $ServerList
    set CurrentServer [lindex $ServerList 0]
    upvar MDFFGUI::servers::$CurrentServer currentserver
    set MDFFGUI::settings::IMDProcs "$currentserver(numprocs)"
  } elseif {![info exists IMDServer]} {
    array set "::MDFFGUI::servers::$servername" $opts
    set ServerList [lappend ServerList $servername]
    set CurrentServer [lindex $ServerList 0]
    upvar MDFFGUI::servers::$CurrentServer currentserver
    set MDFFGUI::settings::IMDProcs "$currentserver(numprocs)"
  }
 #puts "adding server: $servername into $ServerList" 
 # set ServerList [list]
 # foreach s [lsort [info vars ::MDFFGUI::servers::*]] {
 #   lappend ServerList [namespace tail $s]
 # }
 # set CurrentServer [lindex $ServerList 0]
}

proc MDFFGUI::gui::get_dir {parent} {
  #variable MDFFGUI::settings::CurrentDir 
  variable CurrentDirLabel

  set pathname [tk_chooseDirectory -parent $parent -initialdir $MDFFGUI::settings::CurrentDir]
  if {$pathname != ""} {
 #   set fname [lindex [split $pathname "/"] end]
 #   $CurrentDirLabel configure -text "Current Working Directory: $fname/"
    set MDFFGUI::settings::CurrentDir $pathname
  }
  
}

proc MDFFGUI::gui::submit_job {} {
  variable namdbin
  variable CurrentServer
  variable IMDHost
  variable CCStep
  variable CCx
  variable CCy
  #variable MDFFGUI::settings::CurrentDir 
  variable IMDInpOrList
  variable xMDFFMapIDs
   
   set xMDFFMapIDs ""  
   mol top $MDFFGUI::settings::MolID
   $MDFFGUI::gui::IMDCCPlot clear
   set CCx ""
   set CCy ""
  # $IMDCCPlot add $CCx $CCy
   $MDFFGUI::gui::IMDCCPlot replot
#    if {$IMDInpOrList == "list"}{   
    upvar MDFFGUI::servers::$CurrentServer currentserver
    set namdbin "$currentserver(namdbin)"
    set exec_command "$namdbin $MDFFGUI::settings::SimulationName-step$MDFFGUI::settings::SimulationStep.namd"
#    } elseif {$IMDInpOrList == "input"}{
#      array set currentserver {
#        jobtype  local
#        namdbin  {namd2 +p%d}
#        maxprocs 12
#        numprocs 12
#        timeout  20
#      }
#    }
 # set currdir [exec pwd]
  set imdjobspec [NAMDRun::submitJob currentserver $exec_command $MDFFGUI::settings::CurrentDir $MDFFGUI::settings::IMDProcs automdff.log] 

  set IMDHost [lindex $imdjobspec 1]
  $MDFFGUI::gui::IMDStatusLabel configure -text "IMD Status: NAMD job submitted."
}

proc MDFFGUI::gui::imd_update {args} {
  
  variable IMDFreq
  variable justChanged
  variable CCx
  variable CCy
  variable CC
  variable CCStep
  variable NumFrames
  variable xMDFFDensityView

  variable DensityID
  variable MapID
  variable xMDFFMapIDs
  #puts $xMDFF
  set ts [molinfo $MDFFGUI::settings::MolID get timesteps]
  $MDFFGUI::gui::IMDStatusLabel configure -text "IMD Status: Step $ts"
  set CMMap ""
  set RMMap ""
  set SMMap ""
  set GMMap ""
  set CM ""
  set RM ""
  set SM ""
  set GM ""
  set newmapids ""
  set DensityIDName ""
  
  #xMDFF map loading
  if {$MDFFGUI::settings::xMDFF} {
    if {$xMDFFMapIDs == ""} {
      for {set i 0} {$i < [llength [$xMDFFDensityView children {}]]} {incr i} {
        lappend xMDFFMapIDs [mol new [file join $MDFFGUI::settings::CurrentDir "xmdff_density$i.dx"]] 
      }
    }
    if {[expr ($ts-$MDFFGUI::settings::Minsteps-100) % $MDFFGUI::settings::RefSteps] == 0 && $justChanged != $ts || $ts == 1 && $justChanged != $ts} {
      for {set id 0} {$id < [llength $xMDFFMapIDs]} {incr id} {
        set mapid [lindex $xMDFFMapIDs $id]
        if {$mapid != ""} {
          set MapReps ""
          set MapColors ""
          set MapMats ""
          #set CMMap [molinfo $MapID get center_matrix]
          #set RMMap [molinfo $MapID get rotate_matrix]
          #set SMMap [molinfo $MapID get scale_matrix]
          #set GMMap [molinfo $MapID get global_matrix]
          for {set i 0} {$i < [molinfo $mapid get numreps] } {incr i} {
            set MapReps [lappend MapReps [lindex [molinfo $mapid get "{rep $i}"] 0]]
            set MapColors [lappend MapColors [lindex [molinfo $mapid get "{color $i}"] 0]]
            set MapMats [lappend MapMats [lindex [molinfo $mapid get "{material $i}"] 0]]
          }
          if {$DensityID != "" && [lsearch [molinfo list] $DensityID] != -1} {
            set DensityIDName [molinfo $DensityID get name]
          }
          mol delete $mapid 
        }
        set CM [molinfo $MDFFGUI::settings::MolID get center_matrix]
        set RM [molinfo $MDFFGUI::settings::MolID get rotate_matrix]
        set SM [molinfo $MDFFGUI::settings::MolID get scale_matrix]
        set GM [molinfo $MDFFGUI::settings::MolID get global_matrix]
        set newmap [mol new [file join $MDFFGUI::settings::CurrentDir "xmdff_density$id.dx"]]
     #   puts "$DensityIDName [molinfo $newmap get name]"
        if {[molinfo $newmap get name] == $DensityIDName || $DensityIDName == ""  } { 
          set DensityID $newmap 
          set DensityIDName [molinfo $newmap get name]
         }
        lappend newmapids $newmap
        mol top $MDFFGUI::settings::MolID
        #molinfo $MDFFGUI::settings::MolID set center_matrix $CM
        #molinfo $MDFFGUI::settings::MolID set rotate_matrix $RM
        #molinfo $MDFFGUI::settings::MolID set scale_matrix $SM
        #molinfo $MDFFGUI::settings::MolID set global_matrix $GM
        #molinfo $MapID set center_matrix $CMMap
        #molinfo $MapID set rotate_matrix $RMMap
        #molinfo $MapID set scale_matrix $SMMap
        #molinfo $MapID set global_matrix $GMMap
        foreach mol [molinfo list] {
          molinfo $mol set center_matrix $CM 
          molinfo $mol set rotate_matrix $RM
          molinfo $mol set scale_matrix $SM
          molinfo $mol set global_matrix $GM
        }


        for {set i 0} {$i < [llength $MapReps]} {incr i} {
          if {$i > 0} { mol addrep $newmap }
          mol modstyle $i $newmap [lindex $MapReps $i]
          mol modcolor $i $newmap [lindex $MapColors $i]
          mol modmaterial $i $newmap [lindex $MapMats $i]
        }

      }
      set justChanged $ts
      set xMDFFMapIDs $newmapids
    }
  }
  #CC plot
  if {$MDFFGUI::settings::CrossCorrelation} {
    MDFFGUI::gui::resizeToActiveTab
    if {$MDFFGUI::settings::CCRes != "" && $CCStep != $ts} {
      set sel [atomselect $MDFFGUI::settings::MolID $MDFFGUI::settings::CCSel]
      #if {$MDFFGUI::settings::xMDFF} { 
       # set DensityMol $MapID 
     # } else {
        set DensityMol $DensityID
     # }
      if {$DensityMol != "" && [lsearch [molinfo list] $DensityMol ] != -1} {
        if {!$MDFFGUI::settings::UseCCThresh} {
          set CCy [lappend CCy [mdffi cc $sel -res $MDFFGUI::settings::CCRes -mol $DensityMol]]
          set CCx [lappend CCx $ts]
        } elseif {[string is double $MDFFGUI::settings::CCThresh] || [string is integer $MDFFGUI::settings::CCThresh] } {
          set CCy [lappend CCy [mdffi cc $sel -res $MDFFGUI::settings::CCRes -mol $DensityMol -thresholddensity $MDFFGUI::settings::CCThresh]]
          set CCx [lappend CCx $ts]
        }
          $MDFFGUI::gui::IMDCCPlot clear
          $MDFFGUI::gui::IMDCCPlot add $CCx $CCy
          $MDFFGUI::gui::IMDCCPlot replot
        set CCStep $ts
      }
    } else {
    }
  }

#real time timeline analysis
#  set currentframes [molinfo $MDFFGUI::settings::MolID get numframes]
#  if {$currentframes > $NumFrames} {
#    ::timeline::timeLineMainCC
#    set ::timeline::lastAnalysisFrame $currentframes
#    set ::timeline::firstAnalysisFrame $NumFrames
#    ::timeline::calcGlobalCcss 22
#    set NumFrames $currentframes 
#    set ::timeline::prevLastFrame [expr $::timeline::dataOrigin + $NumFrames]
#  }

}

proc MDFFGUI::gui::imd_connect {} {
  
  variable IMDCCPlot
  #variable vmd_timestep

  set attempt_delay 1500   ;# delay between attepmts to connect in ms
  set attempt_timeout 10000 ;# timeout in ms
   
   
  mol top $MDFFGUI::settings::MolID
  #$IMDCCPlot clear
  $MDFFGUI::gui::IMDStatusLabel configure -text "IMD Status: Trying to connect (waiting for NAMD)..."
  
  update idletasks
  for { set timecounter 0 } { $timecounter <= $attempt_timeout } {incr timecounter $attempt_delay} {
    if ![catch { imd connect $MDFFGUI::gui::IMDHost $MDFFGUI::settings::IMDPort }] {
      imd keep 0
      #setstate "connected"
      #trace add variable ::vmd_timestep($imdmol) write AutoMDFF::tracetimestep
      trace add variable ::vmd_timestep($MDFFGUI::settings::MolID) write ::MDFFGUI::gui::imd_update   
      imd keep $MDFFGUI::settings::IMDKeep
      return
    }
    # else give NAMD more time
    after $attempt_delay
  }

}


proc MDFFGUI::gui::get_cryst {parent} {
  variable CurrentCrystPDB
  
  variable CrystOption

  set pathname [tk_getOpenFile -defaultextension ".pdb" -filetypes {{pdb {.pdb}} {all {*}}} -parent $parent -initialdir $MDFFGUI::settings::CurrentDir]
  if {$pathname != ""} {
#    $CurrentCrystPDB configure -text "Symmetry PDB: $fname"
    set MDFFGUI::settings::CurrentCrystPath $pathname
    set CrystOption "-crystpdb [file tail $pathname]"
  }
}

proc MDFFGUI::gui::get_refs {parent} {
  variable CurrentRefs
  

  set pathname [tk_getOpenFile -defaultextension ".mtz, .cif" -filetypes {{cif {.cif}} {mtz {.mtz}} {all {*}}} -parent $parent -initialdir $MDFFGUI::settings::CurrentDir]
  if {$pathname != ""} {
    #set fname [lindex [split $pathname "/"] end]
   # $CurrentRefs configure -text "Reflection Data File: $fname"
    set MDFFGUI::settings::CurrentRefsPath $pathname
  }
}

proc MDFFGUI::gui::make_gridpdb {num seltext} {
  file mkdir $MDFFGUI::settings::CurrentDir 
  mdff gridpdb -psf $MDFFGUI::settings::CurrentPSFPath -pdb $MDFFGUI::settings::CurrentPDBPath -o [file join $MDFFGUI::settings::CurrentDir "gridpdb$num.pdb"] -seltext $seltext
  mol delete top
}

proc MDFFGUI::gui::make_fixedpdb {} {
  #variable MDFFGUI::settings::CurrentDir 
  
  
  file mkdir $MDFFGUI::settings::CurrentDir
  set selall [atomselect $MDFFGUI::settings::MolID all] 
  set sel [atomselect $MDFFGUI::settings::MolID $MDFFGUI::settings::FixedPDBSelText]
  $selall set $MDFFGUI::settings::FixedColumn 0
  $sel set $MDFFGUI::settings::FixedColumn 1
  $selall writepdb [file join $MDFFGUI::settings::CurrentDir "fixed.pdb"]
}

proc MDFFGUI::gui::make_griddx {num file} {
  
  #variable MDFFGUI::settings::CurrentDir 
  variable MDFFSetup
  
  variable DensityID

  if {$MDFFGUI::settings::CurrentDXPath == ""} {
    tk_messageBox -message "Please provide a density map file."
  } else {
    file mkdir $MDFFGUI::settings::CurrentDir 
    mdff griddx -i $file -o [file join $MDFFGUI::settings::CurrentDir "griddx$num.dx"] -threshold $MDFFGUI::settings::GriddxThreshold
    $MDFFSetup configure -state normal
   # set DensityID [mol new $MDFFGUI::settings::CurrentDXPath]
   # mol inactive $DensityID
   # mol off $DensityID
#    mol top $MDFFGUI::settings::MolID
  }
}

proc MDFFGUI::gui::get_density {parent} {
  variable CurrentDXFile
  variable MDFFSetup 
  variable MapID

  set pathname [tk_getOpenFile -defaultextension ".dx, .situs" -filetypes {{dx {.dx}} {situs {.situs}} {mrc {.mrc}} {ccp4 {.ccp4}} {all {*}}} -parent $parent -initialdir $MDFFGUI::settings::CurrentDir]
  if {$pathname != ""} {
    set MDFFGUI::settings::CurrentDXPath $pathname
#    set MapID [mol new $MDFFGUI::settings::CurrentDXPath]
    if {$MDFFGUI::settings::CurrentPSFPath != "" && $MDFFGUI::settings::CurrentPDBPath != ""} {
      $MDFFSetup configure -state normal
    }
  }
}


proc MDFFGUI::gui::load_density {} {
  variable MapID
  variable DensityView 

  set MapID [mol new $MDFFGUI::settings::CurrentDXPath]
 
  $DensityView insert {} end -id $MapID -values [list $MDFFGUI::settings::CurrentDXPath $MDFFGUI::settings::GridPDBSelText $MDFFGUI::settings::GScale $MDFFGUI::settings::GriddxThreshold]
 
   
  if { [winfo exists ".add_density"] } {
    wm withdraw ".add_density"
  }
}

proc MDFFGUI::gui::load_density_xmdff {} {
 # variable MapID
  variable xMDFFDensityView 

  #set MapID [mol new $MDFFGUI::settings::CurrentDXPath]
 
  set id [$xMDFFDensityView insert {} end -values [list $MDFFGUI::settings::CurrentRefsPath $MDFFGUI::settings::GridPDBSelText $MDFFGUI::settings::GScale $MDFFGUI::settings::Mask $MDFFGUI::settings::MaskRes $MDFFGUI::settings::MaskCutoff $MDFFGUI::settings::AverageMap $MDFFGUI::settings::BSharp $MDFFGUI::settings::BSharpFactor]]
 
 set MDFFGUI::settings::Mask 0
 set MDFFGUI::settings::AverageMap 0
 set MDFFGUI::settings::BSharp "none"
 set MDFFGUI::settings::BSharpFactor 0

  if { [winfo exists ".add_density_xmdff"] } {
    wm withdraw ".add_density_xmdff"
  }
}

proc MDFFGUI::gui::remove_densityview {densityview id} {
  mol delete $id 
  $densityview delete $id
}

proc MDFFGUI::gui::load_struct {} {
  
  
  variable GridpdbGen
  variable FixedpdbGen
  
  variable MDFFSetup
  variable GenerateXBonds
  
  if {$MDFFGUI::settings::CurrentPSFPath == ""} {
    tk_messageBox -message "Please provide a PSF file."
  } elseif {$MDFFGUI::settings::CurrentPDBPath == ""} {
    tk_messageBox -message "Please provide a PDB file."
  } else {  
    set MDFFGUI::settings::MolID [mol new $MDFFGUI::settings::CurrentPSFPath]
    mol addfile $MDFFGUI::settings::CurrentPDBPath
    
    if {$MDFFGUI::settings::CurrentDXPath != ""} {
      $MDFFSetup configure -state normal
    }
  }
#  $GridpdbGen configure -state normal
#  $FixedpdbGen configure -state normal
#  $GenerateXBonds configure -state normal
}

proc MDFFGUI::gui::generate_xbonds {} {
  
  
  variable Extrabonds
  #variable MDFFGUI::settings::CurrentDir 
#  variable CButtonSS
#  variable CButtonCis
#  variable CButtonChi


#  variable CButtonSS
 # variable CButtonChi
 # variable CButtonCis
  
  set Extrabonds ""

  file mkdir $MDFFGUI::settings::CurrentDir 
  if {$MDFFGUI::settings::SSRestraints} { 
    ssrestraints -psf $MDFFGUI::settings::CurrentPSFPath -pdb $MDFFGUI::settings::CurrentPDBPath -o [file join $MDFFGUI::settings::CurrentDir "ssrestraints.txt"] -hbonds
    lappend Extrabonds "ssrestraints.txt" 
    mol delete top 
  }
  if {$MDFFGUI::settings::CispeptideRestraints} { 
    cispeptide restrain -o [file join $MDFFGUI::settings::CurrentDir "cispeptide.txt"] -mol $MDFFGUI::settings::MolID
    lappend Extrabonds "cispeptide.txt"
  }
  if {$MDFFGUI::settings::ChiralityRestraints} {
    chirality restrain -o [file join $MDFFGUI::settings::CurrentDir "chirality.txt"] -mol $MDFFGUI::settings::MolID
    lappend Extrabonds "chirality.txt"
  }
  if {[llength $Extrabonds] == 0} {
    set Extrabonds 0
  }

}

proc MDFFGUI::gui::get_pdb {} {
  variable CurrentPDBFile
 
  set pathname [tk_getOpenFile -defaultextension ".pdb" -filetypes {{pdb {.pdb}} {all {*}}}]
  if {$pathname != ""} {
#    $CurrentPDBFile configure -text "PDB File: $fname"
    set MDFFGUI::settings::CurrentPDBPath $pathname
  }
}

proc MDFFGUI::gui::get_psf {} {
  variable CurrentPSFFile
  

  set pathname [tk_getOpenFile -defaultextension ".psf" -filetypes {{psf {.psf}} {all {*}}}]
  if {$pathname != ""} {
 #   $CurrentPSFFile configure -text "PSF File: $fname"
    set MDFFGUI::settings::CurrentPSFPath $pathname
  }
}

proc MDFFGUI::gui::get_map {} {
  set pathname [tk_getOpenFile -defaultextension ".mrc" -filetypes {{mrc {.mrc}} {ccp4 {.ccp4}} {situs {.sit}} {dx {.dx}} {all {*}}}]
  if {$pathname != ""} {
    if {$MDFFGUI::maptoolssettings::MapToolsMolID != ""} {mol off $MDFFGUI::maptoolssettings::MapToolsMolID}  
    set mapmol [mol new $pathname] 
    set MDFFGUI::maptoolssettings::MapToolsMolID $mapmol
    set MDFFGUI::maptoolssettings::MapToolsVolID 0
  }
}

proc MDFFGUI::gui::get_map2 {} {
  set pathname [tk_getOpenFile -defaultextension ".mrc" -filetypes {{mrc {.mrc}} {ccp4 {.ccp4}} {situs {.sit}} {dx {.dx}} {all {*}}}]
  if {$pathname != ""} {
    if {$MDFFGUI::maptoolssettings::MapToolsMol2ID != ""} {mol off $MDFFGUI::maptoolssettings::MapToolsMol2ID}  
    set mapmol [mol new $pathname] 
    set MDFFGUI::maptoolssettings::MapToolsMol2ID $mapmol
    set MDFFGUI::maptoolssettings::MapToolsVol2ID 0
  }
}

proc MDFFGUI::gui::get_maptools_struct {} {
  set pathname [tk_getOpenFile -defaultextension ".pdb" -filetypes {{pdb {.pdb}} {all {*}}}]
  if {$pathname != ""} {
    if {$MDFFGUI::maptoolssettings::MapToolsStructMolID != ""} {mol off $MDFFGUI::maptoolssettings::MapToolsStructMolID}  
    set mapstructmol [mol new $pathname] 
    set MDFFGUI::maptoolssettings::MapToolsStructMolID $mapstructmol
  }
}

proc MDFFGUI::gui::set_mol_text {args} {
  variable MolMenuText 
  if { ! [catch { molinfo $MDFFGUI::settings::MolID get name } name ] } {
    set MolMenuText "$MDFFGUI::settings::MolID: $name"
  } elseif {$MDFFGUI::settings::MolID == -1} {
    set MolMenuText "none"    
  } else {
    set MolMenuText $MDFFGUI::settings::MolID    
  }
}

proc MDFFGUI::gui::set_mapmol_text {args} {
  variable MapToolsMolMenuText
  if { ! [catch { molinfo $MDFFGUI::maptoolssettings::MapToolsMolID get name } name ] } {
    set MapToolsMolMenuText "$MDFFGUI::maptoolssettings::MapToolsMolID: $name"
  } elseif {$MDFFGUI::maptoolssettings::MapToolsMolID == -1} {
    set MapToolsMolMenuText "none"    
  } else {
    set MapToolsMolMenuText $MDFFGUI::maptoolssettings::MapToolsMolID    
  }
}

proc MDFFGUI::gui::set_mapvol_text {args} {
  variable MapToolsVolMenuText
  if {$MDFFGUI::maptoolssettings::MapToolsVolID == ""} {
    set MapToolsVolMenuText 0    
  } else {
    set MapToolsVolMenuText $MDFFGUI::maptoolssettings::MapToolsVolID    
  }
}

proc MDFFGUI::gui::set_mapvol2_text {args} {
  variable MapToolsVol2MenuText
  if {$MDFFGUI::maptoolssettings::MapToolsVol2ID == ""} {
    set MapToolsVol2MenuText 0    
  } else {
    set MapToolsVol2MenuText $MDFFGUI::maptoolssettings::MapToolsVol2ID    
  }
}

proc MDFFGUI::gui::set_mapmol2_text {args} {
  variable MapToolsMol2MenuText
  if { ! [catch { molinfo $MDFFGUI::maptoolssettings::MapToolsMol2ID get name } name ] } {
    set MapToolsMol2MenuText "$MDFFGUI::maptoolssettings::MapToolsMol2ID: $name"
  } elseif {$MDFFGUI::maptoolssettings::MapToolsMol2ID == -1} {
    set MapToolsMol2MenuText "none"    
  } else {
    set MapToolsMol2MenuText $MDFFGUI::maptoolssettings::MapToolsMol2ID    
  }
}

proc MDFFGUI::gui::set_mapstructmol_text {args} {
  variable MapToolsStructMolMenuText 
  if { ! [catch { molinfo $MDFFGUI::maptoolssettings::MapToolsStructMolID get name } name ] } {
    set MapToolsStructMolMenuText "$MDFFGUI::maptoolssettings::MapToolsStructMolID: $name"
  } elseif {$MDFFGUI::maptoolssettings::MapToolsStructMolID == -1} {
    set MapToolsStructMolMenuText "none"    
  } else {
    set MapToolsStructMolMenuText $MDFFGUI::maptoolssettings::MapToolsStructMolID    
  }
}

proc MDFFGUI::gui::set_map_stats {args} { 
  if {$MDFFGUI::maptoolssettings::MapToolsMolID != "none"} {
    set MDFFGUI::maptoolssettings::MapToolsOrigin [voltool info origin -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID]    
    set MDFFGUI::maptoolssettings::MapToolsCOM [voltool com -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID]    
    set MDFFGUI::maptoolssettings::MapToolsXsize [voltool info xsize -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID]    
    set MDFFGUI::maptoolssettings::MapToolsYsize [voltool info ysize -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID]    
    set MDFFGUI::maptoolssettings::MapToolsZsize [voltool info zsize -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID]    
    
    set minmax [voltool info minmax -mol $MDFFGUI::maptoolssettings::MapToolsMolID -vol $MDFFGUI::maptoolssettings::MapToolsVolID]    
    set MDFFGUI::maptoolssettings::MapToolsMin [lindex $minmax 0]    
    set MDFFGUI::maptoolssettings::MapToolsMax [lindex $minmax 1]    
    #For clamping command, set defaults to actual range
    set MDFFGUI::maptoolssettings::ClampMin [lindex $minmax 0]    
    set MDFFGUI::maptoolssettings::ClampMax [lindex $minmax 1]    
    #For range command, set defaults to actual range
    set MDFFGUI::maptoolssettings::RangeMin [lindex $minmax 0]    
    set MDFFGUI::maptoolssettings::RangeMax [lindex $minmax 1]    
  }
}

proc MDFFGUI::gui::fill_mol_menu {args} {
  variable nullMolString
  set name [lindex $args 0]
  $name delete 0 end

  set molList {}
  foreach mm [array names ::vmd_initialize_structure] {
    if { $::vmd_initialize_structure($mm) != 0} {
      lappend molList $mm
      $name add radiobutton -variable MDFFGUI::settings::MolID \
      -value $mm -label "$mm [molinfo $mm get name]" \
    }
  }

  #set if any non-Graphics molecule is loaded
  if {[lsearch -exact $molList $MDFFGUI::settings::MolID] == -1} {
    if {[lsearch -exact $molList [molinfo top]] != -1} {
      set MDFFGUI::settings::MolID [molinfo top]
    } else { set MDFFGUI::settings::MolID $nullMolString }
  }
}

proc MDFFGUI::gui::fill_mapmol_menu {args} {
  variable nullMolString
  set name [lindex $args 0]
  $name delete 0 end

  set molList {}
  foreach mm [array names ::vmd_molecule] {
    if { $::vmd_molecule($mm) != 0 && [lsearch [molinfo list] $mm] != -1 } {
      lappend molList $mm
      $name add radiobutton -variable MDFFGUI::maptoolssettings::MapToolsMolID \
      -value $mm -label "$mm [molinfo $mm get name]" \
    }
  }

  #set if any non-Graphics molecule is loaded
  if {[lsearch -exact $molList $MDFFGUI::maptoolssettings::MapToolsMolID] == -1} {
    if {[lsearch -exact $molList [molinfo top]] != -1} {
      set MDFFGUI::maptoolssettings::MapToolsMolID [molinfo top]
    } else { set MDFFGUI::maptoolssettings::MapToolsMolID $nullMolString }
  }
}

proc MDFFGUI::gui::fill_mapvol_menu {args} {
  variable nullVolString
  set name [lindex $args 0]
  $name delete 0 end
  if { $MDFFGUI::maptoolssettings::MapToolsMolID != "none" && [lsearch -exact [molinfo list] $MDFFGUI::maptoolssettings::MapToolsMolID] != -1 } {
    for {set i 0} {$i < [molinfo $MDFFGUI::maptoolssettings::MapToolsMolID get numvolumedata]} {incr i} {
      $name add radiobutton -variable MDFFGUI::maptoolssettings::MapToolsVolID \
      -value $i -label $i \
    }
  }

}

proc MDFFGUI::gui::fill_mapvol2_menu {args} {
  variable nullVolString
  set name [lindex $args 0]
  $name delete 0 end
  
  if { $MDFFGUI::maptoolssettings::MapToolsMol2ID != "none" && [lsearch -exact [molinfo list] $MDFFGUI::maptoolssettings::MapToolsMol2ID] != -1 } {
    for {set i 0} {$i < [molinfo $MDFFGUI::maptoolssettings::MapToolsMol2ID get numvolumedata]} {incr i} {
      $name add radiobutton -variable MDFFGUI::maptoolssettings::MapToolsVol2ID \
      -value $i -label $i \
    }
  }
}

proc MDFFGUI::gui::fill_mapmol2_menu {args} {
  variable nullMolString
  set name [lindex $args 0]
  $name delete 0 end

  set molList {}
  foreach mm [array names ::vmd_molecule] {
    if { $::vmd_molecule($mm) != 0 && [lsearch [molinfo list] $mm] != -1 } {
      lappend molList $mm
      $name add radiobutton -variable MDFFGUI::maptoolssettings::MapToolsMol2ID \
      -value $mm -label "$mm [molinfo $mm get name]" \
    }
  }

  #set if any non-Graphics molecule is loaded
  if {[lsearch -exact $molList $MDFFGUI::maptoolssettings::MapToolsMol2ID] == -1} {
    if {[lsearch -exact $molList [molinfo top]] != -1} {
      set MDFFGUI::maptoolssettings::MapToolsMol2ID [molinfo top]
    } else { set MDFFGUI::maptoolssettings::MapToolsMol2ID $nullMolString }
  }
}

proc MDFFGUI::gui::revert_volid {args} {
  set MDFFGUI::maptoolssettings::MapToolsVolID 0
}

proc MDFFGUI::gui::revert_vol2id {args} {
  set MDFFGUI::maptoolssettings::MapToolsVol2ID 0
}

proc MDFFGUI::gui::fill_mapstructmol_menu {args} {
  variable nullMolString
  set name [lindex $args 0]
  $name delete 0 end

  set molList {}
  foreach mm [array names ::vmd_initialize_structure] {
    if { $::vmd_initialize_structure($mm) != 0} {
      lappend molList $mm
      $name add radiobutton -variable MDFFGUI::maptoolssettings::MapToolsStructMolID \
      -value $mm -label "$mm [molinfo $mm get name]" \
    }
  }

  #set if any non-Graphics molecule is loaded
  if {[lsearch -exact $molList $MDFFGUI::maptoolssettings::MapToolsStructMolID] == -1} {
    if {[lsearch -exact $molList [molinfo top]] != -1} {
      set MDFFGUI::maptoolssettings::MapToolsStructMolID [molinfo top]
    } else { set MDFFGUI::maptoolssettings::MapToolsStructMolID $nullMolString }
  }
}

proc MDFFGUI::gui::generate_histogram {} {
  global HistPlot 
  set nbins $MDFFGUI::maptoolssettings::HistPlotBins
  if {[info exists HistPlot]} { catch {$HistPlot quit}}
  
  global MAPMOL
  global MAPVOL
  set MAPMOL $MDFFGUI::maptoolssettings::MapToolsMolID 
  set MAPVOL $MDFFGUI::maptoolssettings::MapToolsVolID
    
  set histreturn [voltool hist -nbins $nbins -mol $MAPMOL -vol $MAPVOL]
  set minmax [voltool info minmax -mol $MAPMOL -vol $MAPVOL]
  set min [lindex $minmax 0]
  set max [lindex $minmax 1]
  set xlist [list]
  set delta [expr {($max - $min) / $nbins}]
 
  foreach {midpt hist} $histreturn {
    lappend xlist $midpt
    lappend histogram $hist
  }
  
  #Make histogram values on log scale to reduce drastic difference between highest and lowest bins.
  if {$MDFFGUI::maptoolssettings::HistLog} {
    for {set i 0} {$i < $nbins} {incr i} {
      lset histogram $i [expr log10([lindex $histogram $i])]
    }
    set sorted [lsort $histogram]   
    set ylabel "Number of voxels (log10)"   
  } else {
    set sorted [lsort -integer $histogram]
    set ylabel "Number of voxels"
  }
  global ymin
  set ymin [lindex $sorted 0]
  global ymax
  set ymax [lindex $sorted end]

  set highhistval 0.5
  set highhist 0
  for {set l 0} {$l < $nbins} {incr l} {
    if {[lindex $histogram $l] > $highhist} {
     set highhist [lindex $histogram $l]
     set highhistval [lindex $xlist $l] 
    }
  }

  mol modstyle 0 $MAPMOL Isosurface $highhistval $MAPVOL 0 0 1 1
  
  set HistPlot [multiplot -x $xlist -y $histogram -title "Density histogram" -xlabel "Density" -ylabel $ylabel -nolines -marker square -fill black -xmin [expr [lindex $xlist 0] - (0.5*$delta)] -xmax [expr [lindex $xlist end] + (0.5*$delta)]]
 
  for {set j 0} {$j < $nbins} {incr j} {
    set left [expr [lindex $xlist $j] - (0.5 * $delta)]
    set right [expr [lindex $xlist $j] + (0.5 * $delta)]
    $HistPlot draw rectangle $left $ymin $right [lindex $histogram $j] -fill "#0000ff" -tags rect$j
  #  $HistPlot add [lindex $xlist $j] [lindex $histogram $j] -marker square -fillcolor black
  
  }

  $HistPlot replot
  
  global bpress
  set bpress 0
    
  global xmaxg
  global xplotming
  global xplotmaxg
  global scalexg
  global xming
  global MapToolsPlotX
  variable [$HistPlot namespace]::xplotmin
  variable [$HistPlot namespace]::xplotmax
  variable [$HistPlot namespace]::scalex
  variable [$HistPlot namespace]::xmin
  variable [$HistPlot namespace]::xmax
  #variable [$HistPlot namespace]::ymin
  set xplotming $xplotmin
  set xplotmaxg $xplotmax
  set scalexg $scalex
  set xming $xmin
  set xmaxg $xmax
  
  bind [$HistPlot getpath].f.cf <ButtonPress> {
    set bpress 1    
    #variable [$HistPlot namespace]::xplotmin
    set x [expr (%x - $xplotming)/$scalexg + $xming]
    if {$x >= $xming && $x <= $xmaxg} {
      set MapToolsPlotX $x
      $HistPlot undraw "line"
      $HistPlot draw line $x $ymin $x $ymax -tag "line"
      mol modstyle 0 $MAPMOL Isosurface $x $MAPVOL 0 0 1 1
    }
  }
  
  bind [$HistPlot getpath].f.cf <ButtonRelease> {
    set bpress 0 
  }
  
  bind [$HistPlot getpath].f.cf <Motion> {
    #variable [$HistPlot namespace]::xplotmin
    set x [expr (%x - $xplotming)/$scalexg + $xming]
    if {$bpress && $x >= $xming && $x <= $xmaxg} {
      set MapToolsPlotX $x
      $HistPlot undraw "line"
      $HistPlot draw line $x $ymin $x $ymax -tag "line"
      mol modstyle 0 $MAPMOL Isosurface $x $MAPVOL 0 0 1 1
    }
  }
}

proc mdffgui_tk {} {
  MDFFGUI::gui::mdffgui
  return $MDFFGUI::gui::w
}
