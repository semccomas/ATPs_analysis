#!/usr/bin/env python2
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
import readline
import sys
readline.parse_and_bind("tab: complete")


resname1 = sys.argv[1]
resname2 = sys.argv[2]
resname3 = sys.argv[3]
dimerval = sys.argv[4]

rep_1 = np.loadtxt('AN_replica_1/RDF/' + dimerval + '/' + resname1 + '.xvg')
rep_2 = np.loadtxt('AN_replica_2/RDF/' + dimerval + '/' + resname2 + '.xvg')
rep_3 = np.loadtxt('AN_replica_3/RDF/' + dimerval + '/' + resname3 + '.xvg')


def find_min(rdf, col, lab):
    density_rough = rdf[:,1]
    distance = rdf[:,0]
    
    density = gaussian_filter(density_rough,sigma=20,mode='mirror')
        
    maxpt = int(np.argwhere(density == max(density)))  #where does first RDF peak start
    step = maxpt                                        
    upcount = 0                                         # same for smooth below, use this to make sure you found a 'global' min
    ## and not just a little bump. sigma is 10 for this reason too
    smooth = 5
    for line in density[maxpt:]:
        if density[step] >= density[step + 1]:   # if current density is larger than the following, we are still on peak
           #print density[step] - density[step + 1], step, 'going down'
           pass
        else:
            #print density[step] - density[step + 1], step, 'going up'
            upcount = upcount + 1               #as above, for smoothing
            final = step - smooth
            if upcount >= smooth:
                break                           #when reached cutoff, stop and report initial base
        step = step + 1
    
    plt.plot(distance, density, color = col, label = lab)
    
    return final, distance, density


fin_1, dist_1, den_1 = find_min(rep_1, 'blue', resname1)
fin_2, dist_2, den_2 = find_min(rep_2, 'orange', resname2)
fin_3, dist_3, den_3 = find_min(rep_3, 'red', resname3)


print dist_1[fin_1], resname1
print dist_2[fin_2], resname2
print dist_3[fin_3], resname3


plt.legend()

np.savetxt('rdf_plotting_dir/'+ resname1 + '.' + dimerval + '.txt', den_1)
np.savetxt('rdf_plotting_dir/'+ resname2 + '.' + dimerval + '.txt', den_2)
np.savetxt('rdf_plotting_dir/'+ resname3 + '.' + dimerval + '.txt', den_3)