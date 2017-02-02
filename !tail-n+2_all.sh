#!/bin/bash
echo "cat and sort all strings (except the 1st ones) from all files"
counter=0
for file in `find ./ -name "points_process*"`
do
   tail -n +2 $file > out_$counter;
   ((counter++))
done
echo "cb rhob R tau cw0 cw1 cw2 cw3 cw4 residual" > ./tmp1
cat ./out_* | sort > ./tmp2
rm ./out_*
cat ./tmp1 ./tmp2 > out_all_points
rm ./tmp*