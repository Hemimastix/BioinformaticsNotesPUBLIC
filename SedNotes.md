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
























