#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 09:28:26 2018

@author: sarahmccomas
"""

import numpy as np
import tables as tb
import sys
import MDAnalysis as md


input_tab = sys.argv[1]
#input_tab = 'POPA.6'
lipid_choice = input_tab.split('.')[0]  #if lipid name is always first this will work!
chain_n = input_tab.split('.')[1]
##########################################
### get system info ready to split ##############
#########################################

system_p = md.Universe('../protein_selection.pdb')
system_l = md.Universe('../lipid_selection.pdb')

protein = system_p.select_atoms('name BB and segid %s' %chain_n)
protein = protein[np.where(protein.atoms.occupancies == 1)]

lipid = system_l.select_atoms('resname %s' %lipid_choice)
lipid = lipid[np.where(lipid.atoms.occupancies == 1)]   #take only the values that are 1 (ie not within 1 of prot)

protein_names = []
for resname, resid in zip(list(protein.residues.resnames), list(protein.residues.resids)):
    protein_names.append(resname + '_' + str(resid))

lipid_names = []
for lipid_id in list(lipid.residues.resids):
    lipid_names.append(lipid_choice + '_' + str(lipid_id))


##########################################
#### open gmx output files##############
#########################################

converters = {0: lambda s: float(s.strip('"'))}

if lipid_choice != "POPC":
    prot_vs_lip = np.loadtxt('%s.xvg' %input_tab, dtype = bool, converters=converters)
elif lipid_choice == "POPC":
    prot_vs_lip1 = np.loadtxt('%s.0-400.xvg' %input_tab, dtype = bool, converters=converters)
    prot_vs_lip2 = np.loadtxt('%s.400-800.xvg' %input_tab, dtype = bool, converters=converters)
    prot_vs_lip = np.vstack((prot_vs_lip1[:-1], prot_vs_lip2)) #they will have one part matching that we need to nix


prot_vs_lip = prot_vs_lip[:,1:]   #cut off ts column

frame_chunk = np.reshape(prot_vs_lip, (len(prot_vs_lip), len(lipid.residues), len(protein.residues)))
### ^^ this will split it into chunks per frame, where each row is a lipid index
## and each column is a residue index
intxn_array = np.transpose((frame_chunk), (2, 1, 0))
### ^^ now we have moved this to chunk per residue, and row is lipid number
### and column is frame ##





#intxn_array = np.random.choice(a=[False, True], size = (437, 115, 201))
#############################################
#### put into tables ##############
#########################################

h5 = tb.open_file('%s.table.h5' %input_tab, 'a')
tab_filter = tb.Filters(complevel=1, complib='blosc', fletcher32=True)
atom = tb.BoolAtom()
group = h5.create_group('/', 'intxn_tab', 'individual group')
tab_array = h5.create_earray(group, 'intxn', atom, (0, np.shape(intxn_array)[1], np.shape(intxn_array)[2]), filters = tab_filter)
for n, resid in enumerate(intxn_array):  
    tab_array.append([resid])

h5.close()



'''

## if you want to read table:
    
import tables as tb
import numpy as np

h5 = tb.open_file('POPA.6.h5', 'r')

a = [] 
for i in xrange(0, 457):
    a.append(h5.root.intxn_tab.intxn[i])

'''

