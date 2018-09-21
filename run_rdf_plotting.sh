lip=CDL
dimer=dimer_1

echo $lip
echo DIMER 1 
for i in AN_replica_1/RDF/$dimer/$lip.1.*xvg; do  #this is just for looping purposes, python will access all 3 directories of relicas
rep1=`basename $i .xvg`
rep2="${rep1/.1./.2.}"  #replace .1. with .2. in name
rep3="${rep1/.1./.3.}" 

python rdf_plotting_SM.py $rep1 $rep2 $rep3 $dimer
echo

done

echo
echo ____________________________________________
echo
echo DIMER 2 

dimer=dimer_2
for i in AN_replica_1/RDF/$dimer/$lip.1.*xvg; do  #this is just for looping purposes, python will access all 3 directories of relicas
rep1=`basename $i .xvg`
rep2="${rep1/.1./.2.}"  #replace .1. with .2. in name
rep3="${rep1/.1./.3.}"

python rdf_plotting_SM.py $rep1 $rep2 $rep3 $dimer
echo

done

echo 
echo 
echo
echo ______________________________________________
