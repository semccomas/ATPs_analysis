#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import numpy as np
import MDAnalysis as MD
import matplotlib.pyplot as plt
import MDAnalysis.analysis.distances
import readline
readline.parse_and_bind("tab: complete")

coord = sys.argv[1]
trajxtc = sys.argv[2]
replica_num = sys.argv[3]

lipid_head_dict= {'POPE':'name NH3 PO4 GL1 GL2', 
                  'POPC':'name NC3 PO4 GL1 GL2',
                  'POPA':'name PO4 GL1 GL2',
                  'CDL2':'name GL0 PO41 GL11 GL21 PO42 GL21 GL22'}

system = MD.Universe(coord, trajxtc)

def make_dist_list(lipid_name, dimer_name, outname_raw, outname_intxn):   
    #print lipid_name + ' ' + dimer_name
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
        dist = MDAnalysis.analysis.distances.distance_array(res_COM, liphead_COM)
        print dist
        dist_list.append(dist)  #do for every frame 

        
    dist_list = np.array(dist_list)
    print 'saving ' + outname_raw
    np.save(outname_raw, dist_list)
    
    lipid_cutoffs={'POPE':10.14, 
                  'POPC':10.12,
                  'POPA':9.7,
                  'CDL2':13}  ##should be dist in A 
    
    intxn_cutoff = lipid_cutoffs[lipid_name]
    intxn = dist_list < intxn_cutoff
    np.save(outname_intxn, intxn)


for lip_head in lipid_head_dict:
    chain = 'A'
    make_dist_list(lip_head, chain, '%s.%s.%s.raw_distances.npy'%(lip_head, chain, replica_num), '%s.%s.%s.intxns.npy'%(lip_head, chain, replica_num))
    
    chain = 'B'
    make_dist_list(lip_head, chain, '%s.%s.%s.raw_distances.npy'%(lip_head, chain, replica_num), '%s.%s.%s.intxns.npy'%(lip_head, chain, replica_num))

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



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
