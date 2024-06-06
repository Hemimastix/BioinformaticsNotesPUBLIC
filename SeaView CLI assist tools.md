# SeaView CLI assist tools and tricks

Short scripts (and tricks) for fixing some SeaView (https://doua.prabi.fr/software/seaview) bugs and/or help with routine edits. Yes, SeaView is getting somewhat advanced in age, but I have not yet found any MSA editing GUI that allows profile alignment on the spot as well as (some) manual editing capacity. (also, AliView (v1.27) has a bug when saving NEXUS files, beware when switching between the two...) Besides, the turn-of-the-millennium graphical interface comforts some of us by reminding us of simpler times :-)

**Wait, it's 2023 and you're still manually curating alignments?** Depends on what for; key reference alignments, especially of ribosomal SSU and LSU rRNA genes, require intimiate touch and correction to truly shine (and work). Algorithms don't have eyes yet.

## Editing a sequence

Opening `Edit sequence` leads to a window where you can (somewhat dangerously) edit a sequence, eg. when aligning partial sequences causes bits of them to end up in non-sensical places -- especially when the sequence end falls in a variable region -- which looks something like this (truncated; this is from a base euk-wide SSU rRNA gene alignment with lots of divergent sequences and introns, hence the spacey-ness):

```
----------------------------------------------------------------------------------------------------   100
----------------------------------------------------------------------------------------------------   200
------------------------------------------------------------------------------------------G---------   300
---------------------A-----------------------------------------------------------------------GGTGGCG   400
TTGTG-----------------------------------------TG----------------------------------------------------   500
-------------------------------------GGCACGCCTG-------------------------------C---------------------   600
----------------------------------------------------------------------------------------------------   700
----------------------------------------------------------------------------------------------------   800
----------------------------------------------------------------------------------------------------   900
----------------------------------------------------------------------------------------------------  1000
----------------------------------------------------------------------------------------------------  1100
----------------------------------------------------------------------------------------------------  1200
----------------------------------------------------------------------------------------------------  1300
----------------------------------------------------------------------------------------------------  1400
---------------------------------------------ACTG-CATCAT---------CTCCTCTGTCG------------------------  1500
----------------------------------------------------------------------------------------------------  1600
----------------------------------------------------------------------------------------------------  1700
----------------------------------------------------------------------------------------------------  1800
----------------------------------------------------------------------------------------------------  1900
----------------------------------------------------------------------------------------------------  2000
------------------------------------------------------------------------AAGGACGG-TAGC---------------  2100
----------C-----------------------------------------------------------------------------------------  2200
----------------------------------------------------------------------------------------------------  2300
----------------------------------------------------------------------------------------------------  2400
----------------------------------------------------------------------------------------------------  2500
----------------------------------------------------------------------------------------------------  2600
----------------------------------------------------------------------------------------------------  2700
---------------------------GTTTGC-------------------------------------------------------------------  2800
----------------------------------------------------------------------------------------------------  2900
----------------------------------------------------------------------------------------------------  3000
----------------------------------------------------------------------------------------------------  3100
----------------------------------------------------------------------------------------------------  3200
------------------------GCAAG-CGG-------------------------------------------------------------------  3300
----------------------------------------------------------------------------------------------------  3400
----------------------------------------------------------------------------------------------------  3500
----------------------------------------------------------------------------CGGTTATCGGAGTT----------  3600
CG--------------------------GC-----------------------------CG---------------------------------------  3700
----------------------------------------------------------------------------------------------------  3800
----------------------------------------------------------------------------------------------------  3900
----------------------------------------------------------------------------------------------------  4000
----------------------------------------------------------------------------------------------------  4100
----------------------------------------------------------------------------------------------------  4200
----------------------------------------------------------------------------------------------------  4300
----------------------------------------------------------------------------------------------------  4400
----------------------------------------------------------------------------------------------------  4500
-------------------------------------------------TTTA-----------------------------------------------  4600
----------------------------------------------------------------------------------------------------  4700
-------------------------------------CTTTGAAAAAATTAGAGTGTTCAAAGCAGG-CGAT----------------------------  4800
--------------------------------------------------------------TGCAATTGAATA-A-CCTAGCATGGGA-----------  4900
--------------------------------------------------------------------------TAATGGA---ATAGGAC-TGTGG---  5000
----------------------------------------------------------------------------------------------------  5100
```

Ie, `-` where inserions fall with numbering on the right-hand side (starting with 100, after the first 100 sites). This is risky to edit by hand, but usually you just need to contain all the incorrectly scattered nucleotides and manually align them or place them outside the trimming mask. In the above example, I know from looking at the alignment that the region begining with `TTTA` and the following `CTTTGAAAAAATT...` block is aligned correctly. The preceding nucleotides can be selected together with `-` and line numbers, copied, and pasted into some temporary text file in a \*nix terminal, eg. `foo.foo`, then:

```bash
 sed 's/[0-9\-]//g' foo.foo | awk NF=NF RS= OFS=
```

`sed 's/[0-9\-]//g'` : removes digits and `-` throughout

`awk NF=NF` : deletes blank lines (`awk NF` alone works as a stand-alone step)

`awk NF=NF RS= OFS=` : concatenates lines, spaces them with space (`=` seems to be equivalent to `=""`);  

This produces `GAGGTGGCGTTGTGTGGGCACGCCTGCACTGCATCATCTCCTCTGTCGAAGGACGGTAGCCGTTTGCGCAAGCGGCGGTTATCGGAGTTCGGCCG`, which can be inserted in the insertions-only line preceding the correctly-aligned block (or over several consecutive lines if need be). Then edit the insertions by deleting or adding `-` *until the line numbers align* -- the fixed width font does the work here. Now go through the incorrectly placed nucleotides and replace with `-` by simply highlight over them and hitting the key until the line numbers align (this is important!). Luckily, this process is fast. Hit `Apply` and edit the rest in alignment viewer as the authors intended.

**Is there a more elegant way to do this?** Probably, but I waste too much time writing single-use scripts and convoluted explanatory files already.

**Can you just delete that part of the sequence instead?** Yes, but lately I've taken to keeping the sequences as whole as I can and sticking the phylogenetically-uninformative-anyway bits in the masked regions to differentiate from sequencing errors or other funny things (and to not inevitably forget about the deleted bits later on).

## Clearing species group names from profile alignments

1. Open `.nexus` file in a text editor

2. Go to `BEGIN SETS;` , it will look something like this:
   
   ```
   BEGIN SETS;
     TAXSET 'aa' = 1-192;
     TAXSET 'b' = 1-193;
     TAXSET 'c' = 1-194;
     TAXSET 'd' = 1-195;
     TAXSET 'e' = 1-196;
     TAXSET 'f' = 1-197;
     TAXSET 'g' = 1-198;
     TAXSET 'h' = 1-199;
     TAXSET 'i' = 1-200;
     TAXSET 'j' = 1-201;
     TAXSET 'k' = 1-202;
     TAXSET 'l' = 1-203;
     TAXSET 'm' = 1-204;
     TAXSET 'n' = 1-205; 
     CHARSET 'Gblocks' = 260-270 272-275 277-280 282-287 312-316 326-330 382-386 395-405 407-442 539 541-547 552-575 859-865 884-886 914-920 1240-1247 1249-1253 1277-1285 1454-1456 1458-1476 1478-1479 1481-1482 1484-1494 1499-1508 1510-1532 1534-1550 1552-1566 2053-2080 2082-2085 2087-2111 2113-2115 2449-2471 2473-2474 2619-2638 2643-2648 2695-2700 2713-2726 2728-2733 3225-3229 3231-3236 3577-3590 3592-3598 3601-3630 3632-3638 3660-3661 4738-4767 4870-4874 4878-4889 4975-4981 4985-4991 4993-4995 5543-5546 5549 5552-5556 5558-5560 5562-5563 5671-5720 5722-5740 5744-5748 5752-5773 5920-5925 5927-5935 5937-5965 5967-5974 5984-5996 6159-6168 6603-6608 6628-6638 6640-6641 6648-6658 6789-6792 7038-7046 7335-7383 7668-7675 7677-7678 7680-7686 8056-8075 8077-8098 8106-8121 8197-8215 8217-8225 8227-8269 8272-8283 8390-8404 8412-8415 8425-8430 8979-8982 9001-9002 9009-9010 9012 9018-9024 9063-9068 9070-9103 9105-9108 9110-9117 9123-9126 9128-9161 9293-9297 9299-9302 9368-9374 9497-9507 9524-9527 9530-9535 9537-9545 9557-9568 9574-9575 9577-9598 9605-9609 9611-9627 9629-9670 9674-9676 9679-9684 9686-9700 9790-9803 9805-9825 9833-9843;
   END;
   ```

3. Delete lines containing TAXSET (careful not to delete `CHARSET 'Gblocks'` and other site groups)

(one could write a script for creating a specific group for all currently included taxa after saving the file...)

## Creating taxon sets from custom list of taxa

The taxa sidebar can be a bit hard to read, and anyway it is often more convenient to select taxa from a calculated tree. I haven't found a way to copy taxon names from the tree viewer in SeaView, so I use FigTree's `taxa` selection mode to generate lists of taxa. SeaView's taxon sets can be a versatile tool for generating multiple subtrees from one main alignment, and I've been toying with using Git to control alignment versions that can quickly grow out of control (bonus: Git lets you annotate changes using commit messages, e.g. "expanded clade X" or "removed poor quality sequences Y and Z"). And, of course, if you break your alignment and don't notice until after hitting 'save', you can always go back to a previous version. Ideally, I'd have one master aligment for each of the small handful of groups I work with, and then subsample as needed to 'spin off' datasets. 

Conveniently, taxon sets look like this in Nexus:

```
BEGIN SETS;
  TAXSET 'a' = 1-96 103-119;
  TAXSET 'b' = 1-97 103-119;
  TAXSET 'c' = 1-98 103-119;
  TAXSET 'sdf' = 1-99 103-119;
  CHARSET 'Gblocks' = 105-141 176-207;
END;
```

(ignore `CHARSET` here, just pointing out where it goes)

The numbers correspond to the taxon number, which is formatted at the start of the file like this:

```
[10] ArachnulaimpatiensEU567294
--------------------------GCCAGTAGTCATATGCTTGTCTCAAA-GATTAAG
CCAT-GCAAGTCTAAGTATAAGC---ATTTATACTGTGA--AACTGCGGAAAGCTCATTA
```

I wrote a quick 'n' dirty Python script that extracts the corresponding numbers to desired taxa, and converts them to a `TAXSET` (and an inverse set). The taxon names must match exactly, but one can include taxa that are not in the Nexus file. This means one can apply common lists to Nexus files with taxa in a different order.

[For the script and its usage click here](https://github.com/Hemimastix/BioinformaticsNotesPUBLIC/blob/main/Scripts/taxon_list_to_Nexus_group.py). 

**TODO**: for profile alignments, one must select all taxa and create a new group, then align against it... which gets repetitive. One could make a script that creates a `TAXSET` called "all" that simply includes `TAXSET 'all' = 1-[numtaxa];` where `numtaxa` is derived from this block at the start of Nexus (and mandatory AFAIK):

```
BEGIN DATA;
  DIMENSIONS NTAX=140 NCHAR=4359;
  FORMAT DATATYPE=DNA
  GAP=-
  ;
```

Here `NTAX=140` can be used to generate `  TAXSET 'all' = 1-140;`; this should overwrite the group called  "all".



## Nexus mask fixer

Sometimes SeaView glitches during profile alignments and parts of the mask ended up shifted by one or more sites. Seems to have something to do with gap-only sites, perhaps the aligner itself quietly ignores or deletes them and the mask ends up decoupled from the sequences -- I recommend regularly running `Edit > Del. gap-only sites`, especially when working with divergent sequences. Finding the starting site of the shift is fairly easy but fixing this issue manually is painful, so I wrote a script: [NexusMaskFixer.py](https://github.com/Hemimastix/BioinformaticsNotesPUBLIC/blob/main/Scripts/NexusMaskFixer.py) -- usage is simply, if a bit annoying: first coordinate of first mask block to begin shifting from followed by how much to shift it by. There can be several shifts, so inspect sequence and run it again. The source file is not modified in the process, a new file is saved with original filename + `.shifted`, so mistakes cost only time, not data.
