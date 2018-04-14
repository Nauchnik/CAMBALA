#!/bin/bash
echo "cat and sort all strings from all files"
counter=0
for file in `find ./ -name "*-wu*"`
do
    sed -i -e '$a\' $file
done
echo "residual cb rhob R tau cws_array depths_array" > ./!results
find ./ -name "*-wu*" | xargs cat | sort > ./!results
