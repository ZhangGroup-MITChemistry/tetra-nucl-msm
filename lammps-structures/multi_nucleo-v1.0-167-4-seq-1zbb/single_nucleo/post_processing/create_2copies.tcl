set cv_value [lindex $argv 0]
#puts "The first argument passed was $cv_value"
#set starting_file_name "cv-$cv_value.dcd"
set starting_file_name "./starting_structure_for2copies/mono.dcd"

mol new {./start_mono.pdb} type {pdb} first 0 last -1 step 1 waitfor -1
animate delete  beg 0 end 0 skip 0 0

puts "File name: $starting_file_name"

mol addfile {./starting_structure_for2copies/mono.dcd} type {dcd} first 0 last -1 step 1 waitfor -1 0

set sel [atomselect top "all"]
$sel writepdb mono1.pdb

# Move by 250 A
$sel moveby {250 0 0}
$sel writepdb mono2.pdb
exit
