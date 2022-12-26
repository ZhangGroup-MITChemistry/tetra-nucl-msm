topo readlammpsdata data.prot_dna2_ions.4vmd

set sel [atomselect top "index 0 to 15563"]

animate delete  beg 0 end 0 skip 0 0

mol addfile {./starting_structure_for2copies/150mM_NaCl_5mM_MgCl2/start_ions_7.5_250.0.dcd} type {dcd} first 0 last -1 step 1 waitfor 1 0

$sel writepdb tmp.pdb

mol delete 0

mol new {./tmp.pdb} type {pdb} first 0 last -1 step 1 waitfor 1

topo writelammpsdata data.tmp full

exit
