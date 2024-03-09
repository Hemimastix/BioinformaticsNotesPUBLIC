### To extract sequences from multiline fasta:

To match between two patterns, keeping only the first pattern:

~~`sed -n '/>$match$/,/^>/{/0/!p}'`~~ doesn't work

(for extracting sequences from multiline fasta)

`//,//` -- match everything between two patterns

`{//p}` -- print pattern inclusive of matches

`{//!p}` -- print no pattern match but only between them

` sed -n '/>AFGERG_foo.0$/,/^>.*/{x;p}' Foo.fasta | sed /^$/d `

`-n` don't print

`//,//` -- match everything between two patterns

`{x;p}` -- exchange hold buffer and pattern buffer, then print

shifts frame up for some reason; delete first blank line. Creepy that this works, I don't understand it yet.

Further references:

[Unix Sed Tutorial : 7 Examples for Sed Hold and Pattern Buffer Operations](https://www.thegeekstuff.com/2009/12/unix-sed-tutorial-7-examples-for-sed-hold-and-pattern-buffer-operations/)

[Sed One-Liners Explained, Part II: Selective Printing of Certain Lines](https://catonmat.net/sed-one-liners-explained-part-two)

Match pattern and collapse all new lines after the header:

```shell
sed -n '/>AFGERG_foo.0$/,/^>.*/{x;p}' Foo.fasta | sed '/^$/d' | \
sed ':a;N;ba; s/\n//2g'
```

`sed ':a;N;ba; s/\n//2g` -- loop pattern across end of file; think rolling circle replication ; then  delete everything starting from second instance of newline (to keep newline after header)

It would be better to remove newlines everywhere not after />/... maybe /^>.*\n/! not sure. The delete first line is important... how to make one-liner?

* later: alias this function! Use xargs?

```shell
#not tested
function ExtractSeq{
    echo "$1" | xargs -I % sed -n '/>%$/,/^>.*/{x;p}' "$2" | sed '/^$/d' | \
    sed ':a;N;ba; s/\n//2g'
}

ExtractSeq(){
    echo "$1" | xargs -I % sed -n '/>%$/,/^>.*/{x;p}' "$2" | sed '/^$/d' | \
    sed ':a;N;ba; s/\n//2g'
}
```

### Return sequence lengths in single-line alignment file

Returns sequence names with lengths excluding spaces `-`

```bash
while read line; do if [ "${line:0:1}" == ">" ] ; then echo "$line"; else echo -n "$line" | sed '/>/! s/-//g' | wc -c; fi; done<SingleLineMSA.fastaa
```

As script:

`GetSeqLengths_fromAln.sh`:

```bash
#!/bin/bash
printf "Expects single-line fasta input. All characters except - counted.\n"
while read line; do if [ "${line:0:1}" == ">" ] ; then echo "$line"; else echo -n "$line" | sed '/>/! s/-//g' | wc -c; fi; done<$1
```

`echo -n` -- echoes variable without adding a newline `\n` (this is important for `wc -c`)



### Find words ending in .sh in a file:

Eg. I want a list of every shell script (conventionally named) mentioned in a makefile:

```shell
sed -n 's/^.*\( .*\.sh\).*$/\1/p' Snakefile
```

(that space is important!)

* `-n` -- only print matching patterns

* `s/^.* ... .*$/` -- match line containing `...`

* `\( ... \)` -- contents inside form a named pattern block

* ` .*\.sh` -- match everything from space to word ending in `.sh` followed by space (could also use word boundaries `\b` instead)
  
  * do not follow by space, or else it missed .sh at the end of a line

* `/\1/` -- replace matched pattern (line) with the contents of pattern block 1

* `p` -- print only matched patterns (note: this method will include whitespaces around matched filenames)

Could in principle define another pattern block `\( ... \)` downstream and probably refer to that as `\2`? Note: this is GNU sed.

```shell
function getSnakeShellScripts{
    sed -n '/^\#/! s/^.*\( .*\.sh\)[ $].*$/\1/p' Snakefile | sort | uniq
}
# adds extra spaces around output file matches
```

`/^\#/!` ignores commented lines starting with `#`

remove `sort | uniq` for order of appearance

Match script name preceded by space OR `"`:

```shell
 function genSnakeShellScripts{
     sed -n '/^\#/! s/^.*[ "]\(.*\.sh\)[ $].*$/\1/p' Snakefile | uniq
 }
```

`[ "]` matches space or double quote, then the `\(` pattern capture begins.

the other `[ $]` only matches .sh followed by space or end of line; avoids capturing .short.fasta or whatever

Note absense of `sort`: this will remove immediate duplicates but maintain order

## Useful patterns

### To match GenBank accession numbers

including those in the older single alphacharacter 5 digit format:

```regex
[A-Z]{1,2}[0-9]{5,6}[.]{0,1}[1-9]{0,1}
```

(this assumes the final .1, .2 etc are optional)

Eg. to fix accession numbers with . replaced by - in the last column:

```regex
sed -r 's/([A-Z]{1,2}[0-9]{5,6})-([0-2]$)/\1.\2/g' 18S28Staxatoexpand_wRCs.tsv >| 18S28Staxatoexpand_wRCs_fixedaccessionid.tsv
```


