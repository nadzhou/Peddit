#
# Read charmm topology file
#
# $Id: readcharmmtop.tcl,v 1.19 2019/06/24 21:28:02 jribeiro Exp $
#

package provide readcharmmtop 1.2

namespace eval ::Toporead:: {
  variable topologylist {}
  variable handlerc 0; ## counter to be used to originate the handler key

}

#############################################################
# Read a topology file and return its contents in form of a #
# structured list.                                          #
# The gui variable sets if a GUI (plugin calles it)         #
#############################################################

proc ::Toporead::read_charmm_topology { file {gui 0}} {
   variable topologylist
   set fd [open $file]
   set atomtypes    {}
   set defaultfirst {}
   set defaultlast  {}
   set declarations {}
   set autogenerate {}
   set residuedata  {}
   set residuenames {}
   set residue      {}
   set grouplist    {}
   set group        {}
   set singlebonds  {}
   set doublebonds  {}
   set triplebonds  {}
   set angles       {}
   set dihedrals    {}
   set impropers    {}
   set donors       {}
   set acceptors    {}
   set residefaultfirst {}
   set residefaultlast  {}
   set cgenff  0 ; # topology file coming from CGenFF server
   set nresid 0
   set ngroupatoms 0
   set unsavedresi 0

   ### print warning on the first time that the readcharmmtop plugin is called to warn the user
   ### of the current changes
   if {[llength $topologylist] == 0 && $gui == 0} {
      puts "Toporead: WARNING The ::Toporead::read_charmm_topology now returns a handler (key) value that can be used to access\
      the digested list of the elements present in the topology files loaded."
      puts "Toporead: WARNING To obtain the digested list previously return by the ::Toporead::read_charmm_topology command,\
      please run the command ::Toporead::topology_from_handler handler."
   }
   foreach line [split [read -nonewline $fd] \n] {
      set line [string trim $line]

      # Skip comments
      if {[string index $line 0]=="!"} { continue }
      if {[string index $line 0]=="*"} { continue }

      set comment [string trimleft [string range $line [expr [string first "!" $line] + 1] end]]
      if {$cgenff == 0 && [string first "RESI" $line] > -1 && [string first "param penalty" $line] > -1} {
         set cgenff 1
      }

      if { [string first "!" $line] != -1 && $cgenff == 0} {
         set line [string range $line 0 [expr [string first "!" $line] - 1]]
      }

      set keyword [string toupper [lindex $line 0]]

      # Look for atomtypes
      if { [string equal $keyword "MASS"] } {
	 set num  [lindex $line 1]
	 set type [string toupper [lindex $line 2]]
	 set mass [lindex $line 3]
	 set elem [lindex $line 4]
	 if {![string is alpha $elem]} { puts "WARNING: Invalid element element symbol!" }
##	 set comment [string trimleft [string trimleft [lrange $line 5 end] "!"]]
	 lappend atomtypes [list $type $mass $elem "$comment"]
	 continue
      }

      # Look for default patches
      if { [string match "DEFA*" $keyword ] } {
	 foreach {key value} [lrange $line 1 end] {
	    if {[string match "FIRS*" $key]} {
	       set defaultfirst $value
	    }
	    if {[string match "LAST" $key]} {
	       set defaultlast $value
	    }
	 }
	 continue
      }

      # Look for upstrea/downstream atom declarations
      if { [string match "DECL*" $keyword ] } {
	 lappend declarations [lrange $line 1 end]
	 continue
      }

      # Look for autogeneration preferences
      if { [string match "AUTO*" $keyword ] } {
	 set autogenerate [lrange $line 1 end]
	 continue
      }

      # Residue definitions
      if { [regexp {(RESI|PRES).*} $keyword] } {
	 if {$nresid} {
	    # Append the data of the previous residue
	    lappend grouplist $group
	    set residue [list $reshead $resname $rescharge $grouplist $singlebonds \
				$doublebonds $triplebonds $angles $dihedrals $impropers $donors $acceptors \
				$residefaultfirst $residefaultlast]
	    lappend residuedata  $residue
	    lappend residuenames $resname
	    set resname   [lindex $line 1]
	    set rescharge [lindex $line 2]
            set reshead [string range $keyword 0 3]
	    set unsavedresi 0
	 } else {
	    set resname   [lindex $line 1]
	    set rescharge [lindex $line 2]
            set reshead [string range $keyword 0 3]
	    set unsavedresi 1
         }

	 incr nresid

	 set residue      {}
	 set grouplist    {}
	 set group        {}
	 set singlebonds  {}
	 set doublebonds  {}
	 set triplebonds  {}
	 set angles       {}
	 set dihedrals    {}
	 set impropers    {}
	 set donors       {}
	 set acceptors    {}
	 set residefaultfirst {}
	 set residefaultlast  {}
	 set ngroupatoms 0
	 continue
      }

      # Group
      if { [string equal {GROUP} $keyword ] } {
	 if {$ngroupatoms} { lappend grouplist $group }
	 set group {}
	 incr ngroupatoms
	 continue
      }

      # Atom
      if { [string equal "ATOM" $keyword] } { 
	 set name   [string toupper [lindex $line 1]]
	 set type   [string toupper [lindex $line 2]]
	 set charge [lindex $line 3]
    set list [list $name $type $charge]
    if {$cgenff == 1} {
      lappend list [string trimright [lindex $line 5] "!"]
    }
	 lappend group $list
	 incr ngroupatoms
	 continue
      }

      # Single bonds
      if { [string equal "BOND" $keyword] } { 
	 foreach {atom1 atom2} [string toupper [lrange $line 1 end]] {
	    lappend singlebonds [list $atom1 $atom2] 
	 }
	 continue
      }

      # Double bonds
      if { [string match "DOUB*" $keyword] } { 
	 foreach {atom1 atom2} [string toupper [lrange $line 1 end]] {
	    lappend doublebonds [list $atom1 $atom2] 
	 }
	 continue
      }

      # Triple bonds
      if { [string match "TRIP*" $keyword] } { 
	 foreach {atom1 atom2} [string toupper [lrange $line 1 end]] {
	    lappend triplebonds [list $atom1 $atom2] 
	 }
	 continue
      }

      # Angles
      if { [string match "ANGL*" $keyword] } { 
	 foreach {atom1 atom2 atom3} [string toupper [lrange $line 1 end]] {
	    lappend angles [list $atom1 $atom2 $atom3] 
	 }
	 continue
      }

      # Dihedrals
      if { [string match "DIHE*" $keyword] } { 
	 foreach {atom1 atom2 atom3 atom4} [string toupper [lrange $line 1 end]] {
	    lappend dihedrals [list $atom1 $atom2 $atom3 $atom4] 
	 }
	 continue
      }

      # Impropers
      if { [string match "IMPR*" $keyword] } { 
	 foreach {atom1 atom2 atom3 atom4} [string toupper [lrange $line 1 end]] {
	    lappend impropers [list $atom1 $atom2 $atom3 $atom4] 
	 }
	 continue
      }

      # Donors
      if { [string match "DONO*" $keyword] } { 
	 foreach {hyd don} [string toupper [lrange $line 1 end]] {
	    lappend donors [list $hyd $don] 
	 }
	 continue
      }

      # Donors
      if { [string match "ACCE*" $keyword] } { 
	 foreach {acc ant} [string toupper [lrange $line 1 end]] {
	    lappend acceptors [list $acc $ant] 
	 }
	 continue
      }

      # Look for residue patches
      if { [string match "PATC*" $keyword ] } {
	 foreach {key value} [lrange $line 1 end] {
	    if {[string match "FIRS*" $key]} {
	       set residefaultfirst $value
	    }
	    if {[string match "LAST" $key]} {
	       set residefaultlast $value
	    }
	 }
	 continue
      }

      # Look for end of file
      if { [string equal "END" $keyword] } {   
	 lappend grouplist $group
	 set residue [list $reshead $resname $rescharge $grouplist $singlebonds \
			 $doublebonds $triplebonds $angles $dihedrals $impropers $donors $acceptors \
			 $residefaultfirst $residefaultlast]
	 lappend residuedata  $residue
	 lappend residuenames $resname
	 set unsavedresi 0
	 break
    }
   }

   close $fd

   if {$unsavedresi} {
      lappend grouplist $group
      set residue [list $reshead $resname $rescharge $grouplist $singlebonds \
		      $doublebonds $triplebonds $angles $dihedrals $impropers $donors $acceptors \
		      $residefaultfirst $residefaultlast]
      lappend residuedata  $residue
      lappend residuenames $resname
   }
   # return [list $file $atomtypes $defaultfirst $defaultlast $declarations $autogenerate \
		    $residuenames $residuedata]

   lappend topologylist [list $file $atomtypes $defaultfirst $defaultlast $declarations $autogenerate \
          $residuenames $residuedata]

   return [::Toporead::originate_handler]
}

proc ::Toporead::topology_get { key {top {}} } {
   if {![llength $top]} {}
   switch $key {
      file   { return [lindex $top 0] }
      types  { return [lindex $top 1] }
      first  { return [lindex $top 2] }
      last   { return [lindex $top 3] }
      decl   { return [lindex $top 4] }
      auto   { return [lindex $top 5] }
      names  { return [lindex $top 6] }
      residues { return [lindex $top 7] }  
      resnames { 
	 set resnamelist {}
	 foreach resi [lindex $top 7] {
	    lappend resnamelist [lindex $resi 1]
	 }
	 return $resnamelist
      }  
   }
   error "topology_get: Unknown key $key"
}

proc ::Toporead::topology_get_resid { top resi {key {}} } {
   # Search for the residue name
   set ind [lsearch [lindex $top 6] $resi]

   if {$ind<0} { return {} }

   # Retrieve residue data
   set residue [lindex [lindex $top 7] $ind]
   if {![llength $key]} {
      return $residue
   }
   
   switch $key {
      head        { return [lindex $residue 0] }
      name        { return [lindex $residue 1] }
      charge      { return [lindex $residue 2] }
      grouplist   { return [lindex $residue 3] }
      singlebonds { return [lindex $residue 4] }
      doublebonds { return [lindex $residue 5] }
      triplebonds { return [lindex $residue 6] }
      angles      { return [lindex $residue 7] }
      dihedrals   { return [lindex $residue 8] }
      impropers   { return [lindex $residue 9] }
      donors      { return [lindex $residue 10] }
      acceptors   { return [lindex $residue 11] }
      patchfirst  { return [lindex $residue 12] }
      patchlast   { return [lindex $residue 13] }
      atomlist    { return [join [lindex $residue 3]] }
      bonds       { return [join [lrange $residue 4 6]] }
   }
   error "topology_get_resid: Unknown key $key"
}

proc ::Toporead::topology_contains_type { top type } {
   return [lsearch [lindex $top 1] "$type *"]
}

proc ::Toporead::topology_contains_resi { top resi } {
   return [lsearch [lindex $top 6] $resi]
}

proc ::Toporead::topology_contains_pres { top pres } {
   set ind [lsearch [lindex $top 6] $pres]
   if {$ind>=0 && [lindex [lindex [lindex $top 7] $ind] 0]=="PRES"} { return $ind }
   return -1
}

proc ::Toporead::topology_resi_contains_type { top resi type } {
   return [lsearch [topology_get_resid $top $resi atomlist] "*? $type ?*"]
}

proc ::Toporead::topology_resi_contains_atom { top resi name } {
   return [lsearch -inline [topology_get_resid $top $resi atomlist] "$name ?*"]
}

proc ::Toporead::topology_resi_contains_atom_with_type { top resi name type} {
   return [lsearch -regexp [topology_get_resid $top $resi atomlist] "$name\\s+$type\\s+.*"]
}

###################################################################
## proc to export the full list of topologies read and stored in ##
## the handler variable.                                         ##
## This proc is to be used while the refactoring of this plugin  ##
## is not complete, and to keep the same output format to be     ##
## used in the current plugins.                                  ##
###################################################################
proc ::Toporead::topology_from_handler { handler } {
   variable topologylist
   if {$handler == ""} {
      puts "Toporead: Please provide a valid topology list handler."
      return
   }

   set listindex [string trimleft $handler "toporead"]
   if {$listindex == -1} {
      puts "Toporead: Invalid topology list handler."
      return
   }
   return [lindex $topologylist $listindex]
}

###################################################################
## Originate the key to identify the handler.
## key = toporead + $handlerc
###################################################################
proc ::Toporead::originate_handler {} {
   variable handlerc
   set key "toporead$handlerc"
   incr handlerc
   return $key
}
