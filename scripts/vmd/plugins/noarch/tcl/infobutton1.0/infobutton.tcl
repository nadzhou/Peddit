#
# $Id: infobutton.tcl,v 1.1 2019/06/07 03:57:15 johns Exp $
#
#==============================================================================
# InfoButton plugin to create an info icon and the window to display the text.
# the bind of the button should be done in the code using these buttons to keep
# some flexibility in terms of the content and where to use these info windows.

#     Example of binding code
#     Creating the icon and placing it inside  the frame. Returns the name of the 
#     icon. The command return the info button name to be use to bind the mouse
#     click.
#     set info [INFOBUTTON::createInfoButton $frame.fcontrol 0 0]
#     
#     Bind the mouse click with the widgent name and send the title of the window ([lindex val 0]),
#     the link to display on the button of the page ([lindex val 1]) and the nested 
#     list with the "text" and font to be use.
#     The text arg is a nested list and each element inside the list is composed by
#     a the text and the font of the text. Their is only three fonts possible, the title
#     = "helvetica 15 bold"; the subtitle "helvetica 12 bold" and the body "helvetica 12".
#
#     bind $info <Button-1> {
#         set val [QWIKMD::MDControlsinfo]
#         #set QWIKMD::link [lindex $val 1]
#         INFOBUTTON::infoWindow mdControlsinfo [lindex $val 0] [lindex $val 1] [lindex $val 2]
#     }
#
#     Example of a proc to format the info to go into the window:
#     proc QWIKMD::nscaPlotInfo {} {
#     
#         set text1 "Contacting Surface Area \n\n"
#         set font1 "title"
#     
#         set link {http://www.ks.uiuc.edu/Research/vmd/}
#         set text [list [list $text1 $font1]]
#         set title "Analysis with the Computational Microscope"
#         return [list ${text} $link $title]
#     
#     }
#
#
# Authors:
#   João V. Ribeiro
#   Beckman Institute for Advanced Science and Technology
#   University of Illinois, Urbana-Champaign
#   jribeiro@ks.uiuc.edu
#   http://www.ks.uiuc.edu/~jribeiro
#=============================================================================

package provide infobutton 1.0

namespace eval ::INFOBUTTON:: {
  ## The link has to be defined at the namespace level so be visable
  ## in the text widget
  set link ""
}

proc INFOBUTTON::createInfoButton {frame row column} {
    image create photo imageaux -format gif -file [file join $::env(INFOBUTTONDIR) info.gif]

    grid [ttk::label $frame.info -image imageaux -anchor center -background [ttk::style lookup TFrame -background] ] -row $row -column $column -sticky e -padx 0 -pady 0

    $frame.info configure -cursor hand1
    
    return $frame.info
}

################################################################################
# Creation of the window with the the text. The systax is very picky, so keep
# the format.
################################################################################
 
proc INFOBUTTON::infoWindow {name text link title} {
    
    set wname ".$name"
    if {[winfo exists $wname] != 1} {
        toplevel $wname
    } else {
        wm deiconify $wname
        return
    }
    
    set INFOBUTTON::link $link
    
    wm geometry $wname 600x400
    grid columnconfigure $wname 0 -weight 2
    grid rowconfigure $wname 0 -weight 2
    ## Title of the windows
    wm title $wname $title ;# titulo da pagina

    grid [ttk::frame $wname.txtframe] -row 0 -column 0 -sticky nsew
    grid columnconfigure  $wname.txtframe 0 -weight 1
    grid rowconfigure $wname.txtframe 0 -weight 1

    grid [text $wname.txtframe.info -wrap word -width 420 -bg white -yscrollcommand [list $wname.txtframe.scr1 set] -xscrollcommand [list $wname.txtframe.scr2 set] -exportselection true] -row 0 -column 0 -sticky nsew -padx 2 -pady 2
    
    
    for {set i 0} {$i <= [llength $text]} {incr i} {
        set txt [lindex [lindex $text $i] 0]
        set font [lindex [lindex $text $i] 1]
        $wname.txtframe.info insert end $txt
        set ini [$wname.txtframe.info search -exact $txt 1.0 end]
        
        set line [split $ini "."]
        set fini [expr [lindex $line 1] + [string length $txt] ]
         
        $wname.txtframe.info tag add $wname$i $ini [lindex $line 0].$fini
        if {$font == "title"} {
            set fontarg "helvetica 15 bold"
        } elseif {$font == "subtitle"} {
            set fontarg "helvetica 12 bold"
        } else {
            set fontarg "helvetica 12"
        } 
        $wname.txtframe.info tag configure $wname$i -font $fontarg
    }


        ##Scroll_BAr V
    scrollbar $wname.txtframe.scr1  -orient vertical -command [list $wname.txtframe.info yview]
    grid $wname.txtframe.scr1  -row 0 -column 1  -sticky ens

    ## Scroll_Bar H
    scrollbar $wname.txtframe.scr2  -orient horizontal -command [list $wname.txtframe.info xview]
    grid $wname.txtframe.scr2 -row 1 -column 0 -sticky swe

    grid [ttk::frame $wname.linkframe] -row 1 -column 0 -sticky ew -pady 2 -padx 2
    grid columnconfigure $wname.linkframe 0 -weight 2
    grid rowconfigure $wname.linkframe 0 -weight 2

    grid [tk::text $wname.linkframe.text -bg [ttk::style lookup $wname.linkframe -background ] -width 100 -height 1 -relief flat -exportselection yes -foreground blue] -row 1 -column 0 -sticky w
    $wname.linkframe.text configure -cursor hand1
    $wname.linkframe.text see [expr [string length $link] * 1.0 -1]
    $wname.linkframe.text tag add link 1.0 [expr [string length $link] * 1.0 -1]
    $wname.linkframe.text insert 1.0 $INFOBUTTON::link link
    $wname.linkframe.text tag bind link <Button-1> {vmd_open_url $INFOBUTTON::link}
      # bind link <Button-1> <Enter>
      $wname.linkframe.text tag configure link -foreground blue -underline true
      $wname.linkframe.text configure -state disabled

     
     $wname.txtframe.info configure -state disabled
}