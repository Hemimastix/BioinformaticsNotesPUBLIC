#!/bin/bash
# v1.2 Yana Eglit 17 Feb 2021
# (Updated for perun; "TMP" removed after sed -i)
# Usage: Recode.sh [binsfile.tsv] [alignment]
# for recoding AAs into 4 categories (to use with MinMax-ChiSq ; or SR4, etc; 0:3 converted to ATCG)
# binsfile: first column: what to recode into; subsequent columns: the AA to change; tab separated

#clear log
>Recode.log

#loading into array
let "i=0"
#check if 4 bins in file:
if [ `wc -l <$1` != 4 ] ; then echo "Are you sure your input bins file contains FOUR categories? Check blank spaces. Exiting." ; exit ; fi

#parsing bins file and loading into array
while read line ; do
	categories+=( "$line" )
	echo "${categories[$i]}"
	let "i+=1"
done<"$1"

#duplicate alignment file
cp "$2" "$2.tmp"

let "changes=0" #counting number of changed AA types
# change AAs to numeric bin representations (to be later changed to ATCG)
for ((j=0 ; j<4 ; j++)) do
	row=( ${categories[$j]} )
	echo "Category ${row[0]} contains: ${row[@]:1}"
	echo "Category ${row[0]} contains: ${row[@]:1}" >>Recode.log
	for AA in ${row[@]:1} ; do
		`sed -i "/^>/! s/$AA/${row[0]}/g" "$2.tmp"`
		let "changes+=1"
	done
done

# testing the correct number of AAs (20) are accounted for.
if [ $changes == 20 ] ; then
	echo "All $changes AAs accounted for"
	# Now convert 0-3 to fake-ATCG for GTR model
	sed "/^>/! s/0/A/g" "$2.tmp" | sed "/^>/! s/1/T/g" | sed "/^>/! s/2/C/g" | sed "/^>/! s/3/G/g" >${2%.*}.recoded.fasta
	echo "Recoding complete"
else
	echo "Incorrect number of AAs! Something wrong. Exiting."
	exit
fi

# Now convert 0-3 to fake-ATCG for GTR model
# sed "/^>/! s/0/A/g" "$2.tmp" | sed "/^>/! s/1/T/g" | sed "/^>/! s/2/C/g" | sed "/^>/! s/3/G/g" >${2%.*}.recoded.fasta
# echo "Recoding complete"
