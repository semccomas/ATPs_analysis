dirname=/data2/segmented_ATP_synthase/replica_2/production
tprfile=ATPs.2.3us.tpr
traj=ATPs2.0-4us.skip10.xtc

for i in *.*.ndx; do

j=${i//.ndx}


gmx rdf -f $dirname/$traj -n $i -s $dirname/$tprfile -xvg none -tu ns -seltype part_res_com -selrpos part_res_com -ref 'com of group 20' -sel 'group 15 and name GL0 PO41 GL11 GL21 PO42 GL21 GL22' -o CDL.1.$j.xvg


gmx rdf -f $dirname/$traj -n $i -s $dirname/$tprfile -xvg none -tu ns -seltype part_res_com -selrpos part_res_com -ref 'com of group 20' -sel 'group 14 and name NH3 PO4 GL1 GL2' -o POPE.1.$j.xvg


gmx rdf -f $dirname/$traj -n $i -s $dirname/$tprfile -xvg none -tu ns -seltype part_res_com -selrpos part_res_com -ref 'com of group 20' -sel 'group 13 and name NC3 PO4 GL1 GL2' -o POPC.1.$j.xvg


gmx rdf -f $dirname/$traj -n $i -s $dirname/$tprfile -xvg none -tu ns -seltype part_res_com -selrpos part_res_com -ref 'com of group 20' -sel 'group 16 and name PO4 GL1 GL2' -o POPA.1.$j.xvg

done
