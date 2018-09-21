rep=1
traj=/data2/segmented_ATP_synthase/replica_$rep/production/ATPs$rep.0-4us.skip10.xtc

python distance_matrix.py AN_replica_$rep/ATPs$rep.0ns.chainref.pdb $traj POPC dimer_1 
python distance_matrix.py AN_replica_$rep/ATPs$rep.0ns.chainref.pdb $traj POPC dimer_2



