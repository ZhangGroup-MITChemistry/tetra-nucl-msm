
set sel [atomselect top all]

set reverse [list]
for {set i 7872} {$i>0} {incr i -1} {
    lappend reverse $i
}
#puts $reverse

set numframes [molinfo top get numframes]
set sel [atomselect top "all"]
for {set i 0} {$i<$numframes} {incr i} {
    animate goto $i

    $sel frame $i
    puts "Setting User data for frame [$sel frame] ..."
    $sel set user $reverse
}

$sel delete
