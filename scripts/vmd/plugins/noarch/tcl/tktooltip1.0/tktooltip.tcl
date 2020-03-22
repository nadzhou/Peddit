#
# $Id: tktooltip.tcl,v 1.1 2019/06/07 03:59:49 johns Exp $
#
#==============================================================================
# TKTOOLTIP auxiliary package to create help balloons (tooltip balloons)
# to be used in the tk GUIs
# 
# Authors:
#   João V. Ribeiro
#   Beckman Institute for Advanced Science and Technology
#   University of Illinois, Urbana-Champaign
#   jribeiro@ks.uiuc.edu
#   http://www.ks.uiuc.edu/~jribeiro
#=============================================================================

package provide tktooltip 1.0

namespace eval ::TKTOOLTIP:: {

    proc balloon {w text} {
        bind $w <Any-Enter> "after 2500 [list [namespace current]::balloonShow %W [list $text]]"
        bind $w <Any-Leave> "destroy %W.balloon"
    }
      
    proc balloonShow {w arg} {
        if {[eval winfo containing  [winfo pointerxy .]]!=$w} {return}
        set top $w.balloon
        catch {destroy $top}
        toplevel $top -bd 1 -bg black
        wm overrideredirect $top 1
        if {[string equal [tk windowingsystem] aqua]}  {
            ::tk::unsupported::MacWindowStyle style $top help none
        }   
        pack [message $top.txt -aspect 10000 -bg lightyellow \
                -font fixed -text $arg]
        set wmx [winfo rootx $w]
        set wmy [expr [winfo rooty $w]+[winfo height $w]]
        wm geometry $top \
          [winfo reqwidth $top.txt]x[winfo reqheight $top.txt]+$wmx+$wmy
        raise $top
    }
}