mol new 1kx5.pdb
mol new fiber-172-4_chain_id_changed.pdb

set length 172
set all0 [atomselect 0 "all"]
set sel0 [atomselect 0 "chain I and name P and resid >= -72 and resid <= 72"]
set start [expr 172 - 144]
set end [expr 172]

set sel1 [atomselect 1 "chain A and name P and resid >= ${start} and resid <= ${end}"]
set m1 [measure fit $sel0 $sel1]
$all0 move $m1
set histone1 [atomselect 0 "not chain I and not chain J"]
$histone1 writepdb histone-1.pdb
measure rmsd $sel0 $sel1

set sel2 [atomselect 1 "chain C and name P and resid >= ${start} and resid <= ${end}"]
set m2 [measure fit $sel0 $sel2]
$all0 move $m2
set histone2 [atomselect 0 "not chain I and not chain J"]
$histone2 writepdb histone-2.pdb
measure rmsd $sel0 $sel2

set sel3 [atomselect 1 "chain E and name P and resid >= ${start} and resid <= ${end}"]
set m3 [measure fit $sel0 $sel3]
$all0 move $m3
set histone3 [atomselect 0 "not chain I and not chain J"]
$histone3 writepdb histone-3.pdb
measure rmsd $sel0 $sel3

set sel4 [atomselect 1 "chain G and name P and resid >= ${start} and resid <= ${end}"]
set m4 [measure fit $sel0 $sel4]
$all0 move $m4
set histone4 [atomselect 0 "not chain I and not chain J"]
$histone4 writepdb histone-4.pdb
measure rmsd $sel0 $sel4

exit