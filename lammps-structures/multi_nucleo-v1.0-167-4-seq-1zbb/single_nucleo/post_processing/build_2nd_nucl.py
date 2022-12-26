#!/usr/bin/env python

from merge_mols import *

# Load in the data of the first nucleosome
nucl1 = lmp.Data()
nucl1.read_from_file('./data.prot_dna', peptideFlag=1)

print 'start replicating the data file ...'
# Replicate the nucleosome for another one, the output is the combined two nucleosomes
nucl2 = replicate(nucl1, shiftz=False)

nucl2.write_to_file('data.prot_dna2', peptideFlag=1)
nucl2.write_to_file('data.prot_dna2.4vmd', peptideFlag=0)
