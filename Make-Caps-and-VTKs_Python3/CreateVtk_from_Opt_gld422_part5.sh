#!/bin/bash
ModelName=gld422

root="/Volumes/Accelsior4M2/gld422/"

timeFile="/Volumes/Accelsior4M2/gld422/TimeFile/gld422.timese"

for TimeStep in $(awk '{print $1}' ${timeFile});  do  
#for TimeStep in 46351; do 
echo $TimeStep
for capNumber in 08 09; do

python OptToVtk_Omer_gld422.py ${root}/Opt/${ModelName}.opt${capNumber}.${TimeStep} pidBak.cfg  

done
done
#print "Hadi Bakalim..."

echo "All Done! Good luck!"

# My default modules to run Models on Gadi are:

#1) pbs   3) intel-mkl/2019.3.199        5) python3/3.7.4   7) openmpi/4.0.2  
# 2) dot   4) intel-compiler/2019.3.199   6) gmt/4.5.18    
