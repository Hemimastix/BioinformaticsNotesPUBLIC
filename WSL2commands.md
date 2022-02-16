[Back](https://github.com/Hemimastix/PrivateNotes#readme)

# Useful WSL2 commands and DOS/UNIX stuff

### Clipboard

Copy standard output in Linux terminal to Windows clipboard (example until pipe):

    cat blah.txt | clip.exe

### Working with Unix *vs.* Dos paths using `wslpath`

    wslpath 'c:\users\rest\of\path'

returns `'/mnt/c/users/rest/of/path'` on Ubuntu. Works with drag-and-drop from Windows Explorer (which drops the Dos path).

To get Windows path to current Unix directory: **WARNING: DO NOT ALTER UNIX FILES IN WINDOWS!!!**

    pwd | xargs -I % wslpath -w '%'

> Usage:
>
>*    -a    force result to absolute path format
>*    -u    translate from a Windows path to a WSL path (default)
>*    -w    translate from a WSL path to a Windows path
>*    -m    translate from a WSL path to a Windows path, with '/' instead of '\'

### Dos to Unix line endings

https://blog.codinghorror.com/the-great-newline-schism/

```bash
#!/bin/bash
# removes ^M end of line characters (carriage return/line feed thing) from Dos files for Unix compatability
tr -d '\r' < "$1" > "${1%.*}.UX.${1##*.}"
```


#### Running pandoc in Powershell

Simply execute `pandoc.exe` --help to test if works; seems to work on 64bit version. Note: ~ doesn't work like in unix, use tab to autocomplete and write full path.

`pandoc.exe C:\Users\yegli\Documents\MakeNotes.md -f gfm -t latex -s -o C:\Users\yegli\Documents\test3.tex`

* gfm = github flavoured markdown
* latex probably needs a lot more options to work properly; does not do well with codeblocks
* html export works okay 
