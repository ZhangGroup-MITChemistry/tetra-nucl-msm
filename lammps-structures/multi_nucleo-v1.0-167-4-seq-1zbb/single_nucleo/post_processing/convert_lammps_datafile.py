#!/usr/bin/env python

import sys
sys.path.append('/Users/xl23/bin/lmp_tools/lib')
import lammps_tools as lmp
import fileinput
import numpy as np
import time
from copy import *



if __name__ == '__main__':
    
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]

    mol = lmp.Data()
    mol.read_from_file(inputfile, peptideFlag=1)

    mol.write_to_file(outputfile, peptideFlag=0)

