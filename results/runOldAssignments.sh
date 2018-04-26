#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
for i in 1 2 3
  do 
     echo "Welcome $i times"
     python basemapper.py "new_results/old_assignments/${i}__clust5_beta100.0.out" "new_results/raw_data/${i}/latlong.csv" "new_results/pics/${i}__clust5_beta100"
 done