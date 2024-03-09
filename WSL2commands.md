[Back](https://github.com/Hemimastix/PrivateNotes#readme)

# Useful WSL2 commands and the wonderful world of the DOS/UNIX frontier

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
> * -a    force result to absolute path format
> * -u    translate from a Windows path to a WSL path (default)
> * -w    translate from a WSL path to a Windows path
> * -m    translate from a WSL path to a Windows path, with '/' instead of '\'

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

## BUG: WSL2 clock stops running while laptop asleep

This causes system date to diverge from reality and may cause problems. To fix:

install `ntpdate` with `sudo apt install ntpdate`

then run:

```
sudo ntpdate pool.ntp.org
```

This connects to a cluster of timeservers: [pool.ntp.org: the internet cluster of ntp servers](https://www.ntppool.org/en/) and retrieves correct time. Apparently there are better ways to do this but should not be necessary for standard bioinformatics applications, probably. (from: https://github.com/microsoft/WSL/issues/4245 )



Some networks can block ntp update, a workaround involves syncing time over http:

```
sudo apt-get install htpdate
sudo htpdate -a google.com
```

[time - ntpdate: no server suitable for synchronization found - Ask Ubuntu](https://askubuntu.com/questions/429306/ntpdate-no-server-suitable-for-synchronization-found)



### Consequences of bug: `apt` fails to update

The WSL2 clock lag issue breaks `apt` updates. If `sudo apt update` produces errors like:

```
bionic-security/InRelease is not valid yet (invalid for another 4h 32min 36s)
```

This indicates WSL2 time does not match Windows system time, nor server time. To get the system time from the Windows machine, run:

```
sudo hwclock --hctosys 
```

Then re-run `sudo apt update` and proceed from there.

The above failure will lead to inability to install any packages, appearing as server connection errors.

[sudo apt update error: &quot;Release file is not yet valid&quot; - Ask Ubuntu](https://askubuntu.com/questions/1096930/sudo-apt-update-error-release-file-is-not-yet-valid)

If in the process of switching from `archive.ubuntu.com` to `old-releases.ubuntu.com` in `/etc/apt/sources.list`  -- from assuming the package issue came from broken repository paths and an outdated OS -- you break `apt`, it can be reconfigured by replacing the above file with the contents here: [How to configure apt-get after you had ubuntu installed. · GitHub](https://gist.github.com/abruzzi/8432986)

```
deb http://archive.ubuntu.com/ubuntu/ raring main restricted universe multiverse
deb http://archive.ubuntu.com/ubuntu/ raring-security main restricted universe multiverse
deb http://archive.ubuntu.com/ubuntu/ raring-updates main restricted universe multiverse
deb http://archive.ubuntu.com/ubuntu/ raring-proposed main restricted universe multiverse
deb http://archive.ubuntu.com/ubuntu/ raring-backports main restricted universe multiverse
deb-src http://archive.ubuntu.com/ubuntu/ raring main restricted universe multiverse
deb-src http://archive.ubuntu.com/ubuntu/ raring-security main restricted universe multiverse
deb-src http://archive.ubuntu.com/ubuntu/ raring-updates main restricted universe multiverse
deb-src http://archive.ubuntu.com/ubuntu/ raring-proposed main restricted universe multiverse
deb-src http://archive.ubuntu.com/ubuntu/ raring-backports main restricted universe multiverse
```

tl;dr check the clock first if on WSL2. Don't go down the rabbitholes I did. There's nothing good there.



**Note:** update `apt-get` if you haven't in a while, if the usual repo sites return 404'd

## `matplotlib findfont()` doesn't see Windows system fonts from WSL2

Sometimes you need the absolute path for a font, here's a way to get it:

From https://stackoverflow.com/questions/40290004/how-can-i-configure-matplotlib-to-be-able-to-read-fonts-from-a-local-path :

```python
#!/usr/bin/env python3
# Imports
import os
import re
import shutil
from glob import glob
from matplotlib import matplotlib_fname
from matplotlib import get_cachedir

# Copy files over
dir_source = '<your-font-directory-here>'  # on WSL: /mnt/c/Windows/Fonts/
dir_data = os.path.dirname(matplotlib_fname())
dir_dest = os.path.join(dir_data, 'fonts', 'ttf')
print(f'Transfering .ttf and .otf files from {dir_source} to {dir_dest}.')
for file in glob(os.path.join(dir_source, '*.[ot]tf')):
    if not os.path.exists(os.path.join(dir_dest, os.path.basename(file))):
        print(f'Adding font "{os.path.basename(file)}".')
        shutil.copy(file, dir_dest)

# Delete cache
dir_cache = get_cachedir()
for file in glob(os.path.join(dir_cache, '*.cache')) + glob(os.path.join(dir_cache, 'font*')):
    if not os.path.isdir(file): # don't dump the tex.cache folder... because dunno why
        os.remove(file)
        print(f'Deleted font cache {file}.')
```

To get font path for a PyQt text instance:

```python
from PyQt5 import QtGui
from matplotlib import font_manager

x = QtGui.QFont().family()  # will be a string
font_manager.find_font(x) 
```

## BUG: running docker in WSL2 causes insane memory leak

Something similar happened after running dpkg I think, either that or mamba. Anyway, to clear up memory if `vmmem` is eating up the bulk of it in task manager: 

 in PowerShell (as admin), run `Restart-Service LxssManager`

**Note:** this will kill all WSL2 processes

See: https://github.com/microsoft/WSL/issues/8725
