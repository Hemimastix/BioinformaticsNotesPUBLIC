# Image/graphical utilities in unix

Useful tips and utilities for working with image and vector files in commandline.

## Timestamps and other metadata from TIFF, JPEG

JPEG metadata can be inspected using the `exif` unix utility. Standard exif tags can be looked up in a [table](https://exiv2.org/tags.html) -- they are standardised with the file format. A sample JPEG file can be looked up using your OS (`properties` in Windows or `info` on MacOS) to determine what values are stored. For example, the original date and time the image data were created is stored here:

`0x9003     36867     Image     Exif.Image.DateTimeOriginal     Ascii     The date and time when the original image data was generated.`

So if the metadata are present, `exif -t 0x9003 Test.jpg` will output the following:

```textile
Tag: 0x9003 ('DateTimeOriginal')
  Format: 2 ('ASCII')
  Components: 20
  Size: 20
  Value: 2018:09:02 03:10:55
```

Of this you need the `Value`, which you can get by grep: `exif -t 0x9003 Test.jpg | grep 'Value:'`. This trick can be used to obtain frame rates from time lapse microscopy images saved individually. Enter the directory with the time lapse JPEGs and run:

```bash
for f in *.jpg ; do exif -t 0x9003 $f | grep 'Value:' ; done
```

This will generate a list like this:

```textile
  Value: 2016:07:14 20:03:05
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:06
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:07
  Value: 2016:07:14 20:03:08
  Value: 2016:07:14 20:03:08
  Value: 2016:07:14 20:03:09
```

Pick a second that is between either end of the sequence (note that this is at the end of the time lapse, and the final image was recorded with a significant delay), use grep to obtain all the values matching it and simply count the lines to get the frame rate for that second. Repeat for a couple other timepoints to ensure the value is consistent. With the example above,

```bash
for f in *.jpg ; do exif -t 0x9003 $f | grep 'Value:' ; done | grep '20:03:07' | wc -l
```

will yield a count of 9 (as will `grep 20:03:06`). This matches the max framerate of the camera used to record this particular image sequence. A better approach, of course, is to record the stated frame rate shortly after image collection.



A similar approach can be used for TIFF metadata. `tiffinfo` obtains the metadata (note: without the filename), and `DateTime`, if present, can be grepped out of it. To obtain a table with the filename as well as the `DateTime`, the following command can be used:

```bash
for f in *.TIF ; do printf "$f\t" ; tiffinfo "$f" | grep DateTime | sed 's/^.* //'; done
```

In addition to time information, these tools can be used to extract other properties from JPEG and TIFF metadata in batch. 

## Converting SVG to PNG and JPEG

`cairosvg` [(source)](https://cairosvg.org/) converts SVG to PNG, PDF, PS (and is also a Python 3.6+ library). It seems to do a better job at SVG --> PNG conversion than `ImageMagick` (below). `-s` denotes scaling factor -- in this example, 5x. Setting `width` without scaling will not increase (or decrease) image resolution. `width` (in pixels) might be unnecessary here, I haven't tested without it. `-d` is DPI, which AFAIK is meaningless in this context but sometimes required a certain way by publishers. However, `cairosvg` lacks JPEG conversion. There also does not seem to be a function to set background colour or transparency in the PNG export, although apparently there are additional functions available via Python script, [e.g. for transparent backgrounds](https://stackoverflow.com/questions/48323809/cairo-library-produce-a-png-file-with-white-background).

The following command worked for converting SVGs (generated by `Illustrator`) to decent quality PNGs:

```bash
$ cairosvg --width 1000 -s 5 -d 300 -o Test.png 'Test.svg'
```

`ImageMagick` [(source)](https://www.imagemagick.org/script/download.php) converts between image types, but also has a variety of other [tools](https://www.imagemagick.org/script/command-line-tools.php). The conversion tool is called as `convert` (without `ImageMagick`).

**warning:** the command `mogrify` modifies the original file!

**warning #2:** excess memory use crashes can lead to the original file being corrupted -- make copies before working with files.

I used the following script on the PNG generated by `cairosvg` to convert to JPEG:

```bash
$ convert 'Test.png' -background White -alpha background 'Test.jpg'
```

`-alpha background` modified the transparent background of the PNG to `-background White` (default is black since JPEG does not support transparency). Note that this method does not save any `exif` metadata with the JPEG file. 

This was converted to a batch script run inside the directory with SVG files to be converted:

```bash
for file in *.svg; do
    cairosvg --width 1000 -s 5 -d 300 -o "${file%.*}.png" "${file}"
    convert "${file%.*}.png" -background White -alpha background "${file%.*}.jpg" 
done
```

This seemed to not have issues with whitespace in filenames even though the individual commands above would fail. Inspect generated JPEGs manually as SVGs of different dimensions might need different `--width` and `--scale`. Note that fonts may be changed in the process.

Other options (untested by author):

- [rasterizer](https://manpages.ubuntu.com/manpages/jammy/man1/rasterizer.1.html)

- [InkScape](https://inkscape.org/doc/inkscape-man.html) CLI (running the Windows version via command line kept crashing on my machine, probably works fine on unix systems). [More CLI documentation on wiki](https://wiki.inkscape.org/wiki/Using_the_Command_Line)