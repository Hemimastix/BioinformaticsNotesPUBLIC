#!/bin/bash
# v1.1 Yana Eglit 24 Feb 2021
# (Updated for Perun: "TMP" removed after sed -i)
# Usage: RecodeMixtureModel.sh [binsfile.tsv] [esmodel.nex]
# For recoding mixture model categories into 4 bins
# Use RecodeAAsV1.sh on alignment files with the same binsfile.tsv
# binsfile: first column: what to recode into; subsequent columns: the AA to change; tab separated
# Check that amino acid alphabet order is conventional in model file! (should be)

usage="Usage: RecodeMixtureModel.sh [binsfile.tsv] [esmodel.nex]"
logfile="RecodeMixtureModel.log"
AAalphabet=( A R N D C Q E G H I L K M F P S T W Y V )
>|"$logfile" #clear logfile

#function to round float values for checksum; use round [NUM] [decimals]
round() {
  printf "%.${2}f" "${1}"
}

#checking file type is correct
if [`wc -l <"$1"` != 4 ] ; then echo "Are you sure your input bins file contains FOUR categories? Check blank spaces. Exiting." ; exit ; fi
[[ ! "`head -n 1 <"${2:?"ERROR: Missing an argument. ${usage}"}"`" =~ "#nexus" ]] && echo "Are you sure the model file is in the expected Nexus format?" && echo "${usage}" && exit

#loading into array
let "i=0"
while read line ; do
	categories+=( "$line" )
	let "i+=1"
done<"$1"

#loading arrays into dynamically generated variables (can be combined w above)
for ((j=0 ; j<4 ; j++)) do
	row=( ${categories[$j]} )
  declare -a "Category_${j}" #generating dynamic variable; can later expand to more bins
  for i in "${row[@]:1}"; do
    eval "Category_${j}"+=\("$i"\) #kludge for dynamically generated arrays
  done
done

#Output of what was loaded into variables
#This section and the above (j<4) can be edited to allow more bins
echo "Category_0 ('A') is ${Category_0[@]}"
echo "Category_1 ('T') is ${Category_1[@]}"
echo "Category_2 ('C') is ${Category_2[@]}"
echo "Category_3 ('G') is ${Category_3[@]}"

#writing above to a log file
printf "Category_0 ('A') is %s\nCategory_1 ('T') is %s\nCategory_2 ('C') is %s\nCategory_3 ('G') is %s\n\n" "${Category_0[*]}" "${Category_1[*]}" "${Category_2[*]}" "${Category_3[*]}" >>"$logfile"

#making a copy of the Nexus file
cp "${2}" "${2%.*}.recoded.nex"
outfile="${2%.*}.recoded.nex"
printf "Parsed %s as binfile.\nAA Alphabet is: %s\nNow reading input Nexus file: %s\n" "${1}" "${AAalphabet[*]}" "${2}" >>"$logfile"

#Now parsing Nexus file
while read line ; do
  #only read lines with frequencies:
  if [[ "${line}" =~ "frequency ESclass" ]]; then
    ESclass=`printf "%s" "${line}" | cut -f 2 -d " "` #extract EsClass
    echo "$ESclass" #keep in to see progress in std ouput
    freq_array=( `printf "%s" "${line}" | cut -f 4-23 -d " "` ) #extract numbers
    tosub="${freq_array[@]}" #to be used for sed substitution later
    unset NumCat_0 #clears temporary variables to be sure
    unset NumCat_1
    unset NumCat_2
    unset NumCat_3
    #Now sort AAs by position in number array and classify into bins for summing
    for i in "${!freq_array[@]}" ; do
      case "${AAalphabet[$i]}" in
        `printf "[%s]" "${Category_0[*]}" | sed "s/ //g"`)
          NumCat_0+=( "${freq_array[$i]}" )
          ;;
        `printf "[%s]" "${Category_1[*]}" | sed "s/ //g"`)
          NumCat_1+=( "${freq_array[$i]}" )
          ;;
        `printf "[%s]" "${Category_2[*]}" | sed "s/ //g"`)
          NumCat_2+=( "${freq_array[$i]}" )
          ;;
        `printf "[%s]" "${Category_3[*]}" | sed "s/ //g"`)
          NumCat_3+=( "${freq_array[$i]}" )
          ;;
        *)
          echo "Invalid character in bins! Exiting"
          exit
          ;;
        esac
    done
    #Adding up values in each bin
    #Corresponds to ATCG in recoded alignment! Keep track of order.
    SUM0_A=`echo "${NumCat_0[@]/%/ +}0" | bc`
    SUM1_T=`echo "${NumCat_1[@]/%/ +}0" | bc`
    SUM2_C=`echo "${NumCat_2[@]/%/ +}0" | bc`
    SUM3_G=`echo "${NumCat_3[@]/%/ +}0" | bc`
    #Final substitution:
    sub="0${SUM0_A} 0${SUM1_T} 0${SUM2_C} 0${SUM3_G}"
    #Checksum: all frequencies must add up to 1 in each line
    checksum=`echo "${sub// / + }0" | bc`
    #Checksum has to be rounded... note: decimals must match exactly!
    [[ `round "${checksum}" 2` != 1.00 ]] && echo "CHECKSUM ERROR! Check bins + math" && exit
    printf "Checksum across all bins for %s is: %s\n\n" "${ESclass}" "${checksum}" >>"$logfile"
    `sed -i "/${ESclass}/ s/${tosub}/${sub}/" "${outfile}"` #The Substitution!
    printf "For %s:\n%s\nchanged to\n%s\n\n" "${ESclass}" "${tosub}" "${sub}" >>"${logfile}"
  fi
done<"$2"

printf "File recoded successfully. Output saved as: %s\nLog file is: %s\n" "${outfile}" "${logfile}"
