#
# $Id: molefacture_balloontext.tcl,v 1.1 2019/06/25 20:49:23 jribeiro Exp $
#
#==============================================================================

######################################################
# retrieve the texto display in the balloon          #
# for the actions buttons and radiobuttons           #
######################################################
proc Molefacture::getMouseBText { widget } {
    set text ""
    switch $widget {
        select {
            set text "Select Atom (Control + s)"
        }
        move {
            set text "Move Fragment (Control + m)"
        }
        moveatom {
            set text "Move Atom (5)"
        }
        rotate {
            set text "Rotate Scene (r)"
        }
        translate {
            set text "Translate Scene (t)"
        }
        addatom {
            set text "Add Atom (Control + a)"
        }
        delatom {
            set text "Delete Atom (Shift + a)"
        }
        addbond {
            set text "Create Bond (Control + b)"
        }
        delbond {
            set text "Delete Bond (Shift + b)"
        }
        incrorder {
            set text "Raise Bond Order (Control + i)\n\
            Lower Bond Order (Shift + i)"
        }
        addh {
            set text "Add Hydrogens (Control + h)"
        }
        duplicate {
            set text "Duplicate Selected Atoms (Control + d)"
        }
        clearsel {
            set text "Clear Table Selection (Control + c)"
        }
        undo {
            set text "Undo Action (Control + z)"
        }
        default {}
    }
    return $text
}