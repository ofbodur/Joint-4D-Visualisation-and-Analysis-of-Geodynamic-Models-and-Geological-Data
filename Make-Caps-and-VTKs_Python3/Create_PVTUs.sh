#!/bin/bash

# Omer  Bodur, 30/5/2022
# Sydney, Australia

root="/Volumes/Accelsior4M2/gld422/PVTUs"

timeFile="/Volumes/Accelsior4M2/gld422/TimeFile/gld422.timese"

for TimeStep in $(awk '{print $1}' ${timeFile});  do  
echo $TimeStep

cp gld422.example.pvtu gld422.${TimeStep}.pvtu

sed -i '' 's/.0.vtu/.'${TimeStep}'.vtu/g' gld422.${TimeStep}.pvtu


done
