#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import numpy as np
import MDAnalysis as MD
import matplotlib.pyplot as plt
import MDAnalysis.analysis.distances
import readline
readline.parse_and_bind("tab: complete")

"""
#coord = raw_input('Coordinate file? ')
#trajxtc = raw_input('Trajectory file? ')
#lipid_name = raw_input('Which lipid? ')

## CDL min 1.3 => 13A
## POPE min 1.014 => 10.14A
## POPC min 1.012 => 10.12A
## POPA min 0.97 => 9.7A

coord = "PC_PE_PS_fc0.1.gro"
trajxtc = "PC_PE_PS_fc0.1.100ns.xtc"
lipid_name = 'POPE'
#chain = 'A' # 'A' or 'B' for each dimer
"""


coord = sys.argv[1]
trajxtc = sys.argv[2]
lipid_name = sys.argv[3]
chain = sys.argv[4]

lipid_head_dict= {'POPE':'name NH3 PO4 GL1 GL2', 
                  'POPC':'name NC3 PO4 GL1 GL2',
                  'POPA':'name PO4 GL1 GL2',
                  'CDL2':'name GL0 PO41 GL11 GL21 PO42 GL21 GL22'}

system = MD.Universe(coord, trajxtc)

protein = system.select_atoms('name BB SC1 SC2 SC3 SC4 and segid %s' %chain)
lipid = system.select_atoms('resname %s' %lipid_name) #want only head groups of PE, need to change for each lipid type


def get_com(protein, lipid):
    res_COM = []
    for line in protein.residues:
        res_COM.append(line.atoms.center_of_mass())
    liphead_COM = []
    for line in lipid.residues:
        headgroup = system.select_atoms('resid %s and %s' %(line.resid, lipid_head_dict[lipid_name]))  #doing it this way because we want to only group residues one time, but want all 4 head groups in same atomselection
        liphead_COM.append(headgroup.center_of_mass())

    
    res_COM = np.array(res_COM, dtype=np.float32)
    liphead_COM= np.array(liphead_COM, dtype=np.float32)
    return res_COM, liphead_COM

dist_list = []
for ts in system.trajectory:
    
    print str(system.trajectory.frame) + ' /' + str(len(system.trajectory)) #update user on progress
    res_COM, liphead_COM = get_com(protein, lipid)
    dist = MDAnalysis.analysis.distances.distance_array(res_COM, liphead_COM, backend = "OpenMP")
    #intxn = dist < 10
    dist_list.append(dist)  #do for every frame 
    print dist
    print
    print
    
dist_list = np.array(dist_list)
#intxn = dist_list < 10 
np.savetxt('distance_array_test.txt', dist_list)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



'''



#intxn_matrix = np.zeros((len(protein.residues),len(lipid.residues), len(system.trajectory))) # #residues, #lipids, #frames in sim
#resid_list = []

#for ts in system.trajectory:
#	pope_resids.append((system.trajectory.time, system.select_atoms('name BB SC1 SC2 SC3 SC4 and around 6 resname %s' %lipid).resids))

#system.trajectory.frame = frame # starting @ 0, last step is 5000 => total # frames 5001

#for ts in system.trajectory:
    #for residue_item in protein.residues:



#RDF VALUE PLOTTING FROM GMX
#command line = gmx rdf -f PC_PE_PS_fc0.1.100ns.xtc -s PC_PE_PS_fc0.1.tpr -tu ns -xvg none -n membreg.ndx  -selrpos res_com
#where membreg.ndx is the headgroup of POPE as well as resid 327 beads

rdfGMX = open("rdf.xvg").read().splitlines()

rdf = []
count = 0
for line in rdfGMX:
    line = line.split()
    tup = (float(line[0]), float(line[1]))
    count = count + 2
    rdf.append(tup)  #updating every other frame => count, line[1] is density val for that frame

rdf = np.array(rdf)
plt.plot(rdf[:,0], rdf[:,1])
plt.xlabel('Distance in nm')
plt.ylabel('Density value')
plt.show()



'''
