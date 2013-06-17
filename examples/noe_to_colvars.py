# Usage: from within VMD "gopython noe_to_colvars.py"
# Reads an X-PLOR style list of assign commands for NOE restraints
# Currently uses a simple harmonic potential, using the geometric average 
# of d_minus and d_plus to define the "allowed" fluctuation amplitude (width keyword).

# Tcl would be more portable: but this was written in a hurry.  Please feel free to improve it.

# Used in Wang et al: http://www.pnas.org/content/110/4/1315


import os, sys, re, math

import VMD
from VMD import atomsel

input = file ('distance.tbl')
output = file ('noe.colvars.in', 'w')

pattern = "\((.*?)\)"

# As often, the chemical notation used in NMR for a new molecule
# is not the same as the one used in MD
drug_conversion_table = {
    "H1" : "HZ",
    "H2" : "HE1",
    "H3" : "HD1",
    "H4" : "HI2",
    "H6#" : "H21 H22 H81 H82 H101 H102",
    "H7#" : "H3 H5 H7",
    "H8#" : "H41 H42 H61 H62 H91 H92"
}


colvar_names = []
colvar_centers = []

# loop all the assign statements
for line in input:

    # verify that this is a NOE "assign" statement
    if (line[:1] == '!'): continue
    words = line.split()
    if (len (words) == 0): continue
    if (words[0] != 'assign'): continue
    selection_strings = re.findall (pattern, line)
    if (len (selection_strings) != 2):
        sys.stderr.write ('Error: found "assign" line with more than two selections; likely not an NOE?\n')
        sys.stderr.write ('Error: string follows:\n')
        sys.stderr.write (line)
        sys.exit (1)


    # Extract selection data
    words_sel_1 = selection_strings[0].split()
    words_sel_2 = selection_strings[1].split()
    segid_1 = words_sel_1[1]
    segid_2 = words_sel_2[1]
    resid_1 = words_sel_1[4]
    resid_2 = words_sel_2[4]
    name_1  = words_sel_1[7]
    name_2  = words_sel_2[7]

    noe_d       = float (words[-3])
    noe_d_minus = float (words[-2])
    noe_d_plus  = float (words[-1])

    # The protein is a tetramer, with chains labeled A through D in the NMR table, which I translated to M2A ... M2D
    # Segid E is the drug molecule in the NMR table, which I translated to "DRUG"
    segid_1 = 'M2'+segid_1.upper()
    segid_2 = 'M2'+segid_2.upper()
    if (segid_1 == 'M2E'): 
        segid_1 = 'DRUG'
        resid_1 = str (1)
    if (segid_2 == 'M2E'):
        segid_2 = 'DRUG'
        resid_2 = str (1)

    colvar_name = segid_1+resid_1+":"+name_1.replace ('#', '*')+"-"+segid_2+resid_2+":"+name_2.replace ('#', '*')
    colvar_names += [colvar_name]
    colvar_centers += [("%8.3f" % noe_d).rjust (len (colvar_name))]

    # Convert drug atoms
    if (segid_1 == 'DRUG'): 
        name_1 = drug_conversion_table[name_1]
    if (segid_2 == 'DRUG'):
        name_2 = drug_conversion_table[name_2]

    # Convert wildcards to regexp syntax for VMD
    if (name_1[-1] == '#'):
        name_1 = '"'+name_1.replace ('#', '.*')+'"'
    if (name_2[-1] == '#'):
        name_2 = '"'+name_2.replace ('#', '.*')+'"'

    group_1 = atomsel.atomsel ("(segid %s and resid %s and name %s)" % (segid_1, resid_1, name_1), 0)
    group_2 = atomsel.atomsel ("(segid %s and resid %s and name %s)" % (segid_2, resid_2, name_2), 0)
    if (len (group_1) == 0):
        sys.stderr.write ("Error: cannot find atoms corresponding to \""+selection_strings[0]+"\"\n")
        sys.exit (1)
    if (len (group_2) == 0):
        sys.stderr.write ("Error: cannot find atoms corresponding to \""+selection_strings[1]+"\"\n")
        sys.exit (1)


    group_1_numbers = (' ').join ([str (anum) for anum in group_1.get ('serial')])
    group_2_numbers = (' ').join ([str (anum) for anum in group_2.get ('serial')])

    output.write ('''
colvar {
    name %s
    width %8.3f
    distanceInv {
        group1 { atomNumbers %s }
        group2 { atomNumbers %s }
    }
}
''' % (colvar_name, math.sqrt (noe_d_minus * noe_d_plus), group_1_numbers, group_2_numbers))


# Write the multi-dimensional harmonic restraint
output.write ('''
harmonic {
    colvars %s
    centers %s
    outputEnergy yes
    forceConstant 30.0
    targetForceConstant 2.0
    targetNumSteps 12000000 # 24 ns
}
''' % ((' ').join (colvar_names), (' ').join (colvar_centers)) )

output.close()
sys.exit()
