#
# cmap.tcl                                                  
#
# v1.0 Heiko Lammert 09/2009
#
# draws contacts in VMD into the given molecule, between all pairs of atoms
# given in the contact list, if both atoms are also part of the given
# selection, or if one each is in both selections if two are specified.
#
# commands:
#
# cmap create <molecule> <contacts> <selection> [selection2]
#
# cmap load <molecule> <path> <selection> [selection2]
#
# - initialize a new contact map
#   return the name cmap<N>, with N a unique integer
#   contacts are specified as a list { { i1 j1 [rmax1] [rmin1] } ... {iN jN [rmaxN] [rminN]} }
#   or loaded from an ASCII file with one contact per line.
#   Indices i,j correspond to the counter "serial" in VMD.
#   if rmax or rmin are given, they limit for which distances the contact is shown.
#   rmax can also be used alone, rmin alone requires an empty rmax.
#
# cmap list
#
# - print the names of all defined contact maps cmap<N>
#
# cmap<N> show / hide
#
# - draw / erase the contacts from map cmap<N>
#   repeated calls to show erase before drawing
#
# cmap<N> list / listR / mol / sel / sel2
#
# - print the contact pairs (& distances) / associated molecule / selections for map<N>
#
# cmap<N> style [newstyle] / color [newcolor]
#
# - print / change the drawing style / color.
#   the defaults are {cylinder $ipos $jpos radius 0.1 resolution 6} / 1 (red).
#   an alternative {line $ipos $jpos width 2 style dashed/solid}.
#   the previous value is returned, and without [newvalue] it is retained. 
#   the globally set material is used, e.g. %graphics top material Opaque.
#
# cmap<N> delete
#
# - destroys the map cmap<N>
#
# common arguments:
#
# <molecule> can be "top", a name, or a number(*).
#
# <selection> can be a number or a name(*).
#
# options marked by (*) are unique and permanent in VMD,
# the others may change their meaning when molecules are renamed
# or representations are replaced, but assignment to a cmap<N>
# is permanently established at the time of its creation.
# 

proc {#} args {} ; # multi line comment hack


namespace eval ::cmap_impl {
namespace export cmap loadmap

variable id -1

variable _Shown
variable _Selection
variable _Selection2
variable _Molecule
variable _GFX
variable _Style
variable _Color
variable _Conts
variable _Ncont
variable _Bounds
}

# global control function for user
proc ::cmap_impl::cmap {key args} {

switch $key {

"create" {
eval "::cmap_impl::create $args"
}
"load" {
eval "::cmap_impl::load $args"
}
"list" {
eval "::cmap_impl::list_maps $args"
}
default {
puts "invalid argument $key\n usage: cmap \[ create / load / list \]"
}
} ; # switch key
} ; # proc cmap


# constructor
proc ::cmap_impl::create { molecule pairs selection {selection2 ""} } {
variable id

variable _Shown
variable _Selection
variable _Selection2
variable _Molecule
variable _GFX
variable _Style
variable _Color
variable _Conts
variable _Ncont
variable _Bounds

set tag "[incr id]"

# translate molecule and selections to permanent values, i.e. molid & repname

if { $molecule=="top"} {
set molecule [molinfo top]
} elseif { [ regexp {[^0-9]} {$molecule} ] } {
for {set i 0} {$i < [molinfo num]} {incr i} {
  set molid [molinfo index $i]
  if { [molinfo $molid get name] == "$molecule" } {
    set molecule $molid
    break
  }
}
}

if { $selection2==$selection } {
set selection2 ""
} ;
# if the same set is selected by different criteria that is not detected

if { [ regexp {^[0-9]+$} $selection ] } {
set selection [ mol repname $molecule $selection ]
}

if { [ regexp {^[0-9]+$} $selection2 ] } {
set selection2 [ mol repname $molecule $selection2 ]
}

set _Shown($tag) 0
set _Selection($tag) "$selection"
set _Selection2($tag) "$selection2"
set _Molecule($tag) "$molecule"
set _Color($tag) "1"
set _GFX($tag) [list]
set _Style($tag) {cylinder $ipos $jpos radius 0.1 resolution 6}

# separate atom indices from contact distances
set n 0
foreach p $pairs {
set _Conts(${tag}_$n) [lsort -integer "[lindex $p 0] [lindex $p 1]"]
set min [lindex $p 3]
set max [lindex $p 2]
if { $min != "" } { set min [ expr "$min*$min" ] }
if { $max != "" } { set max [ expr "$max*$max" ] }
set _Bounds(${tag}_$n) "$max $min"
incr n
}
set _Ncont($tag) $n

# autogenerate a proc cmap<N> forwarding user input to ::cmap_impl::redirect
set cmd      "proc ::cmap$tag {key args}"
set body     ""
set body [ concat $body "variable tag \{ $tag \} ;"         ]
set body [ concat $body "eval \"::cmap_impl::redirect \$key \$tag \$args\" " ]
lappend cmd $body
eval $cmd

# autogenerate a handler to redraw cmap<N>
set cmd "proc ::cmap_impl::cmap${tag}_frame_handler {name mol event}"
set body ""
set body [ concat $body "variable _Shown ;" ]
set body [ concat $body "if {\$mol==$_Molecule($tag)} {" ]
set body [ concat $body "if {\$_Shown($tag)} {" ]
set body [ concat $body "::cmap_impl::show $tag } }" ]
lappend cmd $body
eval $cmd
# register this handler
trace add variable ::vmd_frame write "::cmap_impl::cmap${tag}_frame_handler"

return "cmap$tag"
}


# constructor loading contacts from a file
proc ::cmap_impl::load { molecule path selection {selection2 ""} } {

#set pairs [list]
#set mapfile [open "$path" r]
#while  { [ gets $mapfile cont ] > 0 } {
#lappend pairs "$cont"
#}
# close $mapfile
set pairs [ ::cmap_impl::loadmap "$path" ]

::cmap_impl::create $molecule $pairs $selection $selection2
}

# function to load a contact map, also exported
proc ::cmap_impl::loadmap { path } {

set pairs [list]
set mapfile [open "$path" r]
while  { [ gets $mapfile cont ] > 0 } {
lappend pairs "$cont"
}
close $mapfile
return $pairs
}


proc ::cmap_impl::list_maps {} {
variable _Shown
foreach tag [ lsort  [ array names _Shown ] ] {
puts "cmap$tag"
}
}


# common implementation called from the autogenerated proc cmap<N>
proc ::cmap_impl::redirect {key tag args} {
switch $key {
"list" {
eval "::cmap_impl::list_pairs $tag $args"  }
"listR" {
eval "::cmap_impl::list_pairs_dist $tag $args" }
"show" {
eval "::cmap_impl::show $tag $args"  }
"hide" {
eval "::cmap_impl::hide $tag $args"  }
"mol" {
eval "::cmap_impl::getmol $tag $args" }
"sel" {
eval "::cmap_impl::getsel $tag $args" }
"sel2" {
eval "::cmap_impl::getsel2 $tag $args" }
"style" {
eval "::cmap_impl::getset_style $tag $args" }
"color" {
eval "::cmap_impl::getset_color $tag $args" }
"delete" {
eval "::cmap_impl::delete $tag $args" }
default {
puts "invalid argument \"$key\"\nusage: cmap$tag \[ show / hide / list / listR / mol / sel / sel2 / style / color / delete \]" }

} ; # switch
}


# individual methods for cmap<N>, called by users via ::cmap_impl::redirect

proc ::cmap_impl::getset_style {tag {cmd ""}} {
variable _Style
if {$cmd == ""} {
  return $_Style($tag)
} else {
set prev $_Style($tag)
set _Style($tag) "$cmd"
return $prev } 
}


proc ::cmap_impl::getset_color {tag {cmd ""}} {
variable _Color
if {$cmd == ""} {
  return $_Color($tag)
} else {
set prev $_Color($tag)
set _Color($tag) "$cmd"
return $prev }
}


proc ::cmap_impl::delete {tag} {
variable _Shown
variable _Molecule
variable _Selection
variable _Selection2
variable _GFX
variable _Style
variable _Color
variable _Conts
variable _Ncont
variable _Bounds

trace remove variable ::vmd_frame write "::cmap_impl::cmap${tag}_frame_handler"

# 2nd condition detects if the molecule was deleted
if { $_Shown($tag) != 0 && [molinfo index $_Molecule($tag)] >= 0 } {
  eval "::cmap_impl::hide $tag"
}
rename "::cmap$tag" ""
unset _Shown($tag)
unset _Molecule($tag)
unset _Selection($tag)
unset _Selection2($tag)
unset _GFX($tag)
unset _Style($tag)
unset _Color($tag)
for {set n 0} {$n < $_Ncont($tag)} {incr n} {
  unset _Conts(${tag}_$n)
  unset _Bounds(${tag}_$n)
}
unset _Ncont($tag)
}


proc ::cmap_impl::getmol {tag} {
variable _Molecule
return $_Molecule($tag)
}


proc ::cmap_impl::getsel {tag} {
variable _Selection
variable _Molecule
set mol $_Molecule($tag)
set sel [ mol repindex $mol $_Selection($tag)] 
set cmd [ molinfo $mol get "{selection $sel}" ]
return "$_Selection($tag) $sel {$cmd}"
}


proc ::cmap_impl::getsel2 {tag} {
variable _Molecule
variable _Selection
variable _Selection2
set mol $_Molecule($tag)
if { $_Selection2($tag) == "" } {
  set name $_Selection($tag)
} else {
  set name $_Selection2($tag)
}
set sel [ mol repindex $mol $name ]
set cmd [ molinfo $mol get "{selection $sel}" ]
return "$name $sel {$cmd}"
}


proc ::cmap_impl::list_pairs {tag} {
variable _Conts
variable _Ncont
set lst [list]
for {set n 0} {$n < $_Ncont($tag)} {incr n} {
  lappend lst $_Conts(${tag}_$n) }
return $lst
}


proc ::cmap_impl::list_pairs_dist {tag} {
variable _Conts
variable _Ncont
variable _Bounds
set lst [list]
for {set n 0} {$n < $_Ncont($tag)} {incr n} {
  lappend lst "$_Conts(${tag}_$n) $_Bounds(${tag}_$n)" }
return $lst
}

# tha main drawing / update function
proc ::cmap_impl::show {tag} {
variable _Shown
variable _Selection
variable _Selection2
variable _Molecule
variable _GFX
variable _Style
variable _Color
variable _Conts
variable _Ncont
variable _Bounds

#puts "show cmap$tag"

set mol  "$_Molecule($tag)"
set sel  "[ mol repindex $mol $_Selection($tag) ]"
if { $_Selection2($tag) != "" } {
set sel2 "[ mol repindex $mol $_Selection2($tag) ]"
} else {
set sel2 ""
}

set asl [ atomselect $mol [ molinfo $mol get "{ selection $sel }" ] ]
set atoms [ $asl get serial ]
set pos   [ $asl get { x y z } ]

if { $sel2 != "" } {
set asl2 [ atomselect $mol [ molinfo $mol get "{ selection $sel2 }" ] ]
set atoms2 [ $asl2 get serial ]
set pos2   [ $asl2 get { x y z } ]
} else {
};

# puts "atoms: $atoms"
# puts "atoms2: $atoms2"

set postodraw [list]
if { $sel2 == "" } {
  for {set n 0} {$n < $_Ncont($tag)} {incr n} {
    set c $_Conts(${tag}_$n)
    set i [ lsearch $atoms [lindex $c 0] ]
    set j [ lsearch $atoms [lindex $c 1] ]
    if  { $i >= 0 && $j >= 0 }   {
      lappend postodraw $n [ lindex $pos  $i ] [ lindex $pos $j ]
    } ;
  } ; 
} else {
  for {set n 0} {$n < $_Ncont($tag)} {incr n} {
    set c $_Conts(${tag}_$n)
    set i [ lsearch $atoms  [lindex $c 0] ]
    set j [ lsearch $atoms2 [lindex $c 1] ]
    if  { $i >= 0 && $j >= 0 }   {
      lappend postodraw $n [ lindex $pos $i ] [ lindex $pos2 $j ]
    } else {
    set i [ lsearch $atoms2 [lindex $c 0] ]
    set j [ lsearch $atoms  [lindex $c 1] ]
      if  { $i >= 0 && $j >= 0 } {
        lappend postodraw $n [ lindex $pos2 $i ] [ lindex $pos $j ]
      } ;
    } ;
  } ;
} ;

# filter by distance
set tmp $postodraw
set postodraw [list]
foreach { n ipos jpos } $tmp {
  set max [lindex $_Bounds(${tag}_$n) 0]
  set min [lindex $_Bounds(${tag}_$n) 1]

  if { $min != "" || $max != "" } {
    set d2 0.0
    foreach d "0 1 2" {
      set d [ expr [lindex $ipos $d] - [lindex $jpos $d] ]
      set d2 [ expr $d2 + $d * $d ]
    }
  }

# puts "$n  $_Conts(${tag}_$n)    \"$min\"  <=  $d2  <=  \"$max\""

if { [ expr { $min == "" || $d2 >= $min } ] \
  && [ expr { $max == "" || $d2 <= $max } ]  } {
  lappend postodraw $n $ipos $jpos  
}
} ; # foreach n
unset tmp

# puts "todraw: $todraw"

if { $_Shown($tag) != 0 } {
  eval "::cmap_impl::hide $tag"
}

eval "graphics $mol color $_Color($tag)"
foreach { n ipos jpos } $postodraw {
set id [ eval "graphics $mol $_Style($tag)" ]
lappend _GFX($tag) $id
}

$asl delete
if { $sel2 != ""} {
$asl2 delete }

set _Shown($tag) 1

return [ llength $_GFX($tag) ]
}


proc ::cmap_impl::hide {tag} {
variable _Molecule
variable _GFX
variable _Shown

# puts "hide cmap$tag"

foreach id $_GFX($tag) {
graphics $_Molecule($tag) delete $id
}
lreplace _GFX($tag) 0 end ; # empty it

set _Shown($tag) 0
}

# enable reloading during testing
#rename cmap ""
#rename loadmap ""

# brings "cmap" and "loadmap" into global namespace
# leave that to the user ?
namespace import ::cmap_impl::*