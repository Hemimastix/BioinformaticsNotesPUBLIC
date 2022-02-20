#! /bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -o Kludge.out
#$ -cwd
#$ -pe threaded 1

#### USAGE ####
# v2.0, for qsub
#
# This script is a workaround for euk transcriptomes expected to be contaminated by other euks.
# Bordor, as configured, spits out one gene per input transcriptome. This can probably be changed, but
# the pipeline is soon to be deprecated... this is a kludge, written by a hack.
# Reference genes (eg. from Bordor.350.refdat.txt, see below) are searched for via blastp and $maxhits
# are split up into individual header files, which are then used to extract contigs from the original
# transcriptome. In other words, top blast hit for each reference gene forms Transcriptome1.fas, second
# hits go into Transcriptome2.fas, etc. These "transcriptomes" will be chimaeric and have no real meaning!
# Bordor will then do the rest with those "transcriptomes" treating them as separate (numbered) input,
# which must then be manually screened by the user for unwanted eukaryotic sequences. No elegance here.
#
# WARNING: change blastp to blastn for nucleotide input!!!
# First, make refdat (Bordor.350.refdat.txt) directory by converting each fasta entry into single line
# fasta file. Call this ./Bordor350refdat/
# Check that input transcriptomes are in InputTranscriptomes/ -- in the directory from which script is
# launched. Output will be in Done/ in same directory level as script called, log files are in same
# level above Done/, blast outputs will be in InputTranscriptomes/.
# Then finish the BordorPrescript.txt as appropriate and incorporate in Bordor batch submission script,
# using the transcriptomes in the Done/ directory.
# Bonus: if filenames linger from previous run, it backs them up as [filename]OLD -- ONCE. You get one undo.
#
# qsub thisscript.sh

#### VARIABLES ####

#number of hits (recommended 2-3)
maxhits=3

quote='"'
path=$PWD

##### checking files #######

if test -e BordorPrescript.txt ; then
	mv BordorPrescript.txt BordorPrescriptOLD.txt
	echo "BordorPrescript.txt already exists; backed up as BordorPrescriptOLD.txt"
fi

if test -e Kludge.log ; then
	mv Kludge.log KludgeOLD.log
	echo "Kludge.log already exists; backed up as KludgeOLD.log"
fi

if test -e headers.1.txt ; then
	mkdir headersOLD
	mv headers.*.txt headersOLD/
	echo "Old headers.*.txt moved to headersOLD/"
fi

####check if single line; make single line version if not
for transcheck in InputTranscriptomes/*.fas ; do
	if test `head -n 50 $transcheck | grep -c '>'` -lt 24  ; then
        mkdir OriginalMultilineBckup
        cp $transcheck ${transcheck%*.fas}.ml.fas
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${transcheck%*.fas}.ml.fas >$transcheck
        mv ${transcheck%*.fas}.ml.fas OriginalMultilineBckup/
        echo "Multiline fasta detected! Backed up to "${transcheck%*.fas}.ml.fas" in OriginalMultilineBackup/ and single line version used" >>Kludge.log
else
        echo "Singleline fasta okay!"
fi
done

###### MAIN ######

mkdir Done/
for transin in InputTranscriptomes/*.fas; do
	transcriptome=${transin##.*/}
#### var $transin contains the path from shell script's CWD to transcriptome; $transcriptome is the name only
	makeblastdb -in $transin -dbtype prot -out ${transcriptome%*.fas*}
	mkdir ${transcriptome%*.fas*}BlastOut
	for query in ./Bordor350refdat/*.fasta ; do
		nquery=${query##.*/}
		echo $nquery >>Kludge.log
		echo ${transcriptome%*.fas*}.${nquery%*.fasta}.blast6 >>Kludge.log
		OUTFMT=${quote}"6 ""qseqid ""sseqid ""pident"${quote}
		COMMAND="blastp -query "$query" -db "$path"/"${transcriptome%*.fas*}" -max_target_seqs "$maxhits" -outfmt "$OUTFMT" >"${transcriptome%*.fas*}"."${nquery%*.fasta}".blast6"
		echo $COMMAND >blastp.sh
		bash blastp.sh
		mv ${transcriptome%*.fas*}.${nquery%*.fasta}.blast6 ${transcriptome%*.fas*}BlastOut/
	done
#### this part extracts headers from blast output (sseqid)
	for blastout in ${transcriptome%*.fas*}BlastOut/*.blast6 ; do
		i=1 ; while [[ $i -le ${maxhits} ]] ; do
			cut -f 2 $blastout | uniq | sed -n ${i}p >>headers.${i}.txt
			i=$((${i}+1))
		done
	done
#### parses headers textfile and extracts sequences from original transcriptome, for each line of blast output ($i)
	i=1 ; while [[ $i -le ${maxhits} ]] ; do
		>${transcriptome%*.fas}${i}.fas
		while read -r line; do
			line="'"$line
			line=$line"'"
			echo $line >>Kludge.log
			GREPC="grep -A 1 "$line" "$transin" >>"${transcriptome%*.fas}${i}".fas"
			echo $GREPC >grep.sh
			bash grep.sh
		done<headers.${i}.txt
#### prevents script from outputting into input folder... (recursion bad!)
		mv $path/${transcriptome%*.fas}${i}.fas $path/Done/
#### this part prepares part of the Bordor script, to be manually checked and edited
        	echo "cd "${transcriptome%*.fas}${i}".fas" >>BordorPrescript.txt
        	echo "python AddPipeline3.0a.py "${transcriptome%*.fas}${i}" ***FULL_NAME*** 1 AA Bordor.350.refdat.txt ./ ***NEXT*** no" >>BordorPrescript.txt
        	echo "cd ../" >>BordorPrescript.txt
        	echo "tar -czf "${transcriptome%*.fas}${i}".tgz "${transcriptome%*.fas}${i} >>BordorPrescript.txt
        	echo "rm -r "${transcriptome%*.fas}${i} >>BordorPrescript.txt
        	echo "mv "${transcriptome%*.fas}${i}".* DONE_DATASET/" >>BordorPrescript.txt
        	echo "--" >>BordorPrescript.txt
#### clears headers temp file (or else they linger until next transcriptome)
		>headers.${i}.txt
		i=$((${i}+1))
	done

done
