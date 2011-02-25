#!/bin/bash

# Retrieve only the proper columns from the files and store them in another file

for file in *.txt
	do
	filecorr=${file/.txt/_reduce.txt}
	echo $file
	gawk 'BEGIN {FS = "\t"} {print $6 "\t" $9 "\t" $10 "\t" $19}' $file > $filecorr
done
