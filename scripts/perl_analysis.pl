set filename "/home/nadzhou/Desktop/1yu5.pdb"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   set filename "C:/Users/Hamdan/Desktop/bonds.txt"
set fbonds [open $filename "w"]

set all [atomselect top "resname LIG"]
set at [$all list]
set mt [molinfo top]
foreach i $at { 
  set p [atomselect top "within 5 of index $i"] 
  set pt [$p list]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
  foreach tl $pt {
    set c 0
    foreach tt $at {                                                                                              
      if { $tt == $tl } { 
	set c 1                                                                                                       
      }	                                                                        
    }
    if { $c == 0} {
      lappend tn $tl
    }
  }
  foreach j $tn { 
    set rj [atomselect top "index $j"]
    set mm [measure bond "$j $i"] 
    if { $mm < 4} { 
      label add Bonds $mt/$j $mt/$i 
      label add Atoms $mt/$j
      set ri [atomselect top "index $i"]
      set pt1a [ $ri get name ]
      set pt1r [ $ri get resid ]
      set pt1rn [ $ri get resname ]
      set pt2a [ $rj get name ]
      set pt2r [ $rj get resid ]
      set pt2rn [ $rj get resname ]
      puts -nonewline $fbonds "$pt1r $pt1rn  $pt1a  : $pt2r $pt2rn $pt2a \t"
      puts -nonewline $fbonds $mm
      puts -nonewline $fbonds "\n"
      unset -nocomplain rj
      unset -nocomplain ri
    } 
  }
  set tn [lreplace $tn 0 0]
}
label hide Bonds all
label hide Atoms all

close $fbonds
