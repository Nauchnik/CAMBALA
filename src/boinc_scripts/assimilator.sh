#!/bin/bash

d=$(date +%Y-%m-%d-%T)
res_file=result_$d

for file in `find ./ -name "*-wu*"`
do
    (cat "${file}"; echo) >> $res_file;
    rm $file;
done
mv ./$res_file /boinc-data/results/$res_file
