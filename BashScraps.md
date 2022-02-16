[Back](https://github.com/Hemimastix/PrivateNotes#readme)



# Bits and Pieces

### Profiles

https://medium.com/@abhinavkorpal/bash-profile-vs-bashrc-c52534a787d3

`~/.profile` – login shell (eg Bourne shell)

`~/.bashrc` – interactive shell profile, ie bash reads this file; this shell is a child process of the login shell

`~/.bash_profile` – login shell, but bash – instead of `~/.profile`

(this is why on my local MacOS machine, the `$PATH` moficiation above was done on `~/.bash_profile` and not on `~/.bashrc` as on Perun.)

### Bash Arrays

Example 1: (OLD)

```bash
#!/bin/bash
#reads list of taxa and extracts each taxon from all .fas files in current working directory
#(next version: can give it gene list too)
#usage: Trim2Taxonlist.sh [taxon list.txt]
#ONLY WORKS WITH SINGLE LINE FASTA FILES!!!

#load taxonlist into array
let "i=0"
while read line ; do
	taxon[$i]="${line}"
	let "i+=1"
done<$1

echo "Final taxa in selection are: ${taxon[*]}"

#process each of the fasta files:
for file in *.fas ; do
	>${file%.*}.subspl.faa
	for tax in ${taxon[@]}; do
		grep -A 1 ">"$tax ${file} >>${file%.*}.subspl.faa
	done
done
```

`${var[*]}`: expands all as one string with spaces between values

`${var[@]}`: expands all as multiple strings separated by spaces

i.e. `for tax in ${taxon[@]}` because `for tax in ${taxon[*]}` would attempt to iterate over a single string of concatenated variables (and not work: would grep for entire array later on)

Another example with arrays of directories:

Consense in batch is complicated because the outfile is just ‘outfile’ and would overwrite. It does not play well with bash (too old?) Solution: copy files into directories named after them:

    for file in *.treefile ; do mkdir ${file%%.*} ; cp $file ${file%%.*}/ ; done

Then run Consense from within each directory.

Then to rename outfile s in batch (provided the ONLY directories are ones with Consense files): (from the directory on top of consense directories)

```bash
dirlist=`find . -mindepth 1 -type d`
for dir in ${dirlist[@]} ; do cp $dir/outfile ${dir:2}.outfile ; done
for dir in ${dirlist[@]} ; do cp $dir/outtree ${dir:2}.outtree ; done
```

`${dir:2}` cuts off the first two characters, which are `./` and would confuse bash.



#### To insert extension before existing extension:
  
    "${1%.*}.UX.${1##*.}"

https://stackoverflow.com/questions/12426659/how-to-extract-last-part-of-string-in-bash

>Note that the `'*'` needs to swap places with the `' '` depending on whether you use `#` or `%`. (The `*` is just a wildcard, so you may need to take off your "regex hat" while reading.)
>
>*	`${A% *}` - remove shortest trailing  `*` (strip the last word)
>*	`${A%% *}` - remove longest trailing  `*` (strip the last words)
>*	`${A#* }` - remove shortest leading `*` (strip the first word)
>*	`${A##* }` - remove longest leading `*` (strip the first words)

##### 
