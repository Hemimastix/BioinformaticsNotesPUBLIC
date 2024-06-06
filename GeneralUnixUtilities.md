[Back](https://github.com/Hemimastix/PrivateNotes/blob/main/README.md)

# Unix utilities

In no particular order.

### $PATH and how unix finds commands

https://linuxize.com/post/how-to-add-directory-to-path-in-linux/

    echo $PATH

brings up directories unix searches for executables, in order of search. For example:

    /usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/Users/yana/genesis/bin

If two executable files have the same name, the one in the path that shows up first will be called.

To temporarily add a path, do:

    export PATH=” [your_path]:$PATH”

Then type the executable script name. This will reset after terminal session. For permanent modifications, for all users edit `/etc/profile`; for your own shell, edit `~/.bashrc` ; eg

    nano ~/.bashrc

(on our lab Mac OS machine: `~/.bash_profile`) ; and add:

 ~~export PATH=\"\$PATH:/Users/yana/Library/Python/2.7/bin\"~~ **DO NOT USE QUOTES**

    export PATH=$PATH:/Users/yana/Library/Python/2.7/bin

To load this path into the current shell session, use source:

    source ~/.bash_profile

### Miscellany

#### checking qstat or any other verbose command without it clogging your terminal:

    qstat | less

#### Max line length in file

For alignment files in FASTA (singleline format!!! if length = 60, suspicious!)

    wc -L

if wc -L option does not exist (eg. on Mac):

```shell
cat [filename] | awk '{print length}' | sort -un | tail -n 1
```

batch:

```shell
for file in *.fas ; do echo $file ; cat $file | awk '{print length}' | sort -un | tail -n 1; done
```

easier to read:

```shell
for file in *.fas ; do echo $file ':' >>Lengths.txt ; cat $file | awk '{print length}' | sort -un | tail -n 1 >> Lengths.txt ; echo "" >> Lengths.txt ; done
```

#### unzip

```shell
gunzip [file.tgz]
tar -zxvf yourfile.tar.gz
xvf for file.tar
```

#### To Check phred33/64:

For older raw reads (check SRA source)

```shell
head -n 40 file.fastq | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding\!";}'
```

from https://www.biostars.org/p/63225/



#### Extract taxa with corresponding numbers from Nexus file

```bash
#!/bin/bash
grep '\[[0-9]*\]' "$1" | sed 's/ /\t/' | sed "s/'//g" > "${1%.*}.taxa.txv"
# keeps [num], inserts tab, removes quotes
echo "Taxa extracted to ${1%.*}.taxa.txv"
```



#### Move empty files to a directory before deletion:

```shell
    mkdir EmptyForDeletion
    find \*.[fileextensions] -type f -size 0 -exec mv {} EmptyForDeletion/ \;
    cd EmptyForDeletion/
```

\[check!!!\]

    rm -i etc...

\* no space between `{}` !!!

#### Running background processes

To be able to put computer to sleep and then resume process:

`nohup` -- No HangUp

```shell
nohup bash script.sh & #redirects standard output to Nohup.out
nohup bash script.sh > Log.out & #to redirect standard output
```

The `&` to indicate end of command appears important -- not sure though

#### Copying delicate stuff

`cp -av` : `-a` is equal to `-dR`, ie preserve attributes (links, ownership, permissions, timestamp, etc) and copy recursively; `-v` is verbose (redirect output to log)

To transfer a very large database on Perun, submitted the following job:

```shell
#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -o ~/CopyDBtoScratch4.18Jan2022.out

cp -av /home/yana/db/ /scratch4/yana/
```

Here, output was redirected to logfile specified to scheduler via `-o`. Note: `cp` does quote strings in `'`, weird characters should not be an issue, can see this in verbose output.

`cp -l` creates symlinks instead of copying them directly.

### Modifying BASH prompt: \$PS1

Bash prompt (the blahblah/path/ \$ bit) lives in the `$PS1` variable. By default, on Mac OS X terminal this is `\h:\W \u\$ ` (mind the spaces!), where `\h` is the host; `\W` is the abbreviated location; and `\u` is user. It's not ideal for our type of work. Before doing anything, copy current `$PS1` contents to another variable as backup.

First, `\u` can be demoted to the start of the line, and by convention, is followed by an `@` sign. Move `\u@` to the start of the sequence. Next, `\W` is rather useless and only shows the current directory, not the full path. Replace with `${PWD}` for full path, and follow by a `/` which is not there by default.

Lastly, colours can be modified with the following construct: `\[\033[x;ym\]FEATURE_TO_COLOUR\[\033[0m\]` , where `x` is the weight (`0` for regular, `1` for bold), and `y` (followed by `m`!) is a numerical representation for the colour -- that maps on colour settings in terminal preferences. For example:

```shell
\[\033[0;30m\] # Black
\[\033[0;31m\] # Red
\[\033[0;32m\] # Green
\[\033[0;33m\] # Yellow
\[\033[0;34m\] # Blue
\[\033[0;35m\] # Purple
\[\033[0;36m\] # Cyan
\[\033[0;37m\] # White

\[\033[0;11m\] # main text colour
```

After a bit of fiddling, the following prompt:

```shell
export PS1="\[\033[1;31m\]\u@\h:\[\033[0;36m\]\${PWD}/ \\[\033[0m\]$ "
```

Looks like `yana@Excavates:/Users/yana/FormattingBS-TMP/ $ ` with `user@host:` in bold dark red and the full path in cyan.

Add to bash profile to make permanent.

**The rest of this section is for nerds!** 

#### Escape characters, `$'string'` expansion, and additional details

For [Ubuntu colours](https://askubuntu.com/questions/941734/whats-the-rgb-values-for-ubuntus-default-terminal-unity): `P0`-`PF` (hex) are the colours in order; followed by the new colour in hex (note: I have not tried this yet)

```shell
#change default colors
echo -en "\e]P0000000" #black
echo -en "\e]P1FF0000" #lightred
echo -en "\e]P200FF00" #lightgreen
echo -en "\e]P3FFFF00" #yellow
echo -en "\e]P40000FF" #lightblue
echo -en "\e]P5FF00FF" #lightmagenta
echo -en "\e]P600FFFF" #lightcyan
echo -en "\e]P7FFFFFF" #highwhite
echo -en "\e]P8808080" #grey
echo -en "\e]P9800000" #red
echo -en "\e]PA008000" #green
echo -en "\e]PB808000" #brown
echo -en "\e]PC000080" #blue
echo -en "\e]PD800080" #magenta
echo -en "\e]PE008080" #cyan
echo -en "\e]PFC0C0C0" #white
```

`echo` `-n` suppresses newline and `-e` prints escape characters as is. Eg `$ echo "\n"` will print `\n` whereas `$ echo -e "\n"` prints a newline.

Note: `\033` and `\x1B` are octal and hex values, respectively (27 in decimal), for the escape character `\e` (apparently also `^[`) ; note extra `[` (or `]` above -- different shell?). `\e` is not part of POSIX standard and went on to be included in GNU Unix (eg. Ubuntu) but not BSD Unix (Mac OS). It was apparently added in Bash 2.0, but not in Mac OS' versions, even Mac OS' Bash 3.x. (bug alert!)

But it gets weirder! Expansion of `$'string'` does accept `\e` and a broader range of escape sequences, even in BSD Unix. For example:

```shell
$ echo $'\e[34m''COLOURS'
```

Prints COLOURS in `[34m` colour.  (`echo $'\a\a\a\a'` plays the bell sound four times). 

To add background but not have it bleed over to next line, command to reset cursor, etc back to default after:

```shell
echo $'\e[44;1;31m''STUFF'$'\e[0;0;0m'
```

From BASH man page:

```
 Words  of  the  form  $'string' are treated specially.  The word expands to string, with
   backslash-escaped characters replaced as specified by the ANSI  C  standard.   Backslash
   escape sequences, if present, are decoded as follows:
          \a     alert (bell)
          \b     backspace
          \e     an escape character
          \f     form feed
          \n     new line
          \r     carriage return
          \t     horizontal tab
          \v     vertical tab
          \\     backslash
          \'     single quote
          \nnn   the  eight-bit  character whose value is the octal value nnn (one to three
                 digits)
          \xHH   the eight-bit character whose value is the hexadecimal value  HH  (one  or
                 two hex digits)
          \cx    a control-x character

   The expanded result is single-quoted, as if the dollar sign had not been present.

   A  double-quoted string preceded by a dollar sign ($) will cause the string to be trans‐
   lated according to the current locale.  If the current locale is C or POSIX, the  dollar
   sign  is  ignored.  If the string is translated and replaced, the replacement is double-
   quoted.
```

Note: `printf` is more standardised, and therefore portable, than `echo`! There may be differences, even within one machine, between `/bin/echo` and built-in `bash` `echo`.

Now, going back to `$PS1` prompt:

```shell
$ echo $'\033[34m''STUFF'
or
'[\033[x;ym\]FEATURE_TO_COLOUR\[\033[0m\]'
```

The `m` indicates text colour.

```shell
Black       0;30     Dark Gray     1;30  
Blue        0;34     Light Blue    1;34  
Green       0;32     Light Green   1;32  
Cyan        0;36     Light Cyan    1;36  
Red         0;31     Light Red     1;31  
Purple      0;35     Light Purple  1;35  
Brown       0;33     Yellow        1;33  
Light Gray  0;37     White         1;37 
```

Note that the closing bracket in `[\033[x;ym]` closes the first `[`, I believe. Text colour is actually indicated by `\033[x;ym` alone.

Other things than can be controlled: cursor position `\033[<1>;<5>f` (or `H` instead of `f`) will place cursor at line 1 and column 5. Cursor can be moved up 2 lines by `\033[<2>A` ; down with `B` instead of `A`; `C` is forward, `D` is back.

Clear screen and move cursor to 0,0: `\033[2J`

Erase line: `\033[K`

Save and restore cursor position: `\033[s` and `\033[u`, respectively (in xterm and nxterm, not all terminal emulators)

Difference between `\033[` and `\033]`: 'control sequence introducer' (CSI) vs 'operating system command' (OSC). With CSI you refer to what the operating system calls a colour, eg "red", whereas with OSC you specify red in hex. ([source](https://askubuntu.com/questions/831971/what-type-of-sequences-are-escape-sequences-starting-with-033)) See Ubuntu example above!

[More reliable reference](https://tldp.org/HOWTO/Bash-Prompt-HOWTO/x329.html) for BASH prompt.

Can combine background and foreground colours`\[\033[44;1;31m\]`  (here, 44 is blue). Other codes instead of 0 and 1 for reg vs. bold: 4: Underscore, 5: Blink, 7: Inverse, and 8: Concealed
