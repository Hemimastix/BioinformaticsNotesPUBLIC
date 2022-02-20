#!/usr/bin/env python3
# Yana Eglit 19.02.22 v1.0
# For shifting SeaView site selection mask (bug in software sometimes unsyncs it from alignment)
# By default, assumes Nexus file containing [CHARSET 'GBLOCKS' = ] or [CHARSET 'GBLOCKS'= ] (should be case insensitive)

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Script shifts blocks at [LOCATION] by [SHIFT]
    Use coordinate of first base of block to shift! (can't shift in between)
    Note: assumes mask is named "charset 'GBLOCKS'", change charset variable if named otherwise.
    Check Nexus file with tail to see region name + format
    Warning: single-site blocks converted to 'fake' pairs; works for SeaView but undo manually if
    problems with other downstream software.
    """)
    parser.add_argument("nexus", type=str, help="Path to Nexus file")
    parser.add_argument("location", type=int, help="First coordinate of block to shift from\n\
                                                   (Note: can only shift whole blocks!)")
    parser.add_argument("shift", type=int, help="Number of bases to shift; + to the right, - to the left")
    parser.add_argument("-c", "--charset", type=str,
                        help="Alternate name of mask/regions if not gblocks (case-insensitive).",
                        nargs="?", default="gblocks")

args = parser.parse_args()

# charset = "gblocks"  # change this if different mask name!
charset = args.charset

searchstr = ("CHARSET \'%s\' =" % charset).upper()  # search function converts target to uppercase and strips whitespace
searchlen = len(searchstr)  # to avoid searching entire very long lines


def generate_lines_that_startwith(string, f):  # function to extract line with mask
    for line in f:
        if (line[0:(searchlen+10)].strip("\t").lstrip().upper().startswith(string) or
                line[0:(searchlen+10)].strip("\t").lstrip().upper().replace("=", " =").startswith(string)):
            # trims possible tabs and whitespaces, converts to uppercase
            return line
    else:
        print("Pattern [%s] not found! Quitting..." % string)
        exit()


# Open file and search for mask, extract list of coordinate pairs
with open(args.nexus, "r") as infile:
    mask = str(generate_lines_that_startwith(searchstr, infile))
    coordinate_pairs = mask.split("= ", 1)[1].rstrip(";\n").split(" ")  # \n important after ; for strip
    for i, v in enumerate(coordinate_pairs):
        if "-" not in str(v):
            coordinate_pairs[i] = str(v+"-"+v)
            print("Warning: single-site length block found at %s; this may cause problems. Fake pair made as: %s" %
                  (v, coordinate_pairs[i]))
    if (list(map(lambda x: x.isdigit(), coordinate_pairs[0].split("-"))) == [True, True] and
            list(map(lambda x: x.isdigit(), coordinate_pairs[(len(coordinate_pairs)-1)].split("-"))) == [True, True]):
        pass
    else:
        print("Incorrect mask format. Expect: number-number\\n\n\nFirst, last pair %s, %s" %
              (coordinate_pairs[0], coordinate_pairs[(len(coordinate_pairs)-1)]))
        exit()

# sort coordinate pairs into those to keep and those to shift
blockstoshift = []
blockstokeep = []
for i, n in enumerate(coordinate_pairs):
    if int(n.split("-")[0]) == args.location:  # int is important!
        blockstoshift = coordinate_pairs[i:]  # write the rest of the coordinates
        break
    else:
        blockstokeep.append(n)  # coordinates prior to LOCATION
        if i == (len(coordinate_pairs)-1):
            print("Location %s not found in mask coordinates. Quitting." % args.location)
            exit()
print("\nBlocks until %s kept as is." % blockstokeep[len(blockstokeep)-1])
print("Blocks from %s to %s to be shifted by %s\n" %
      (blockstoshift[0].split("-")[0],
       blockstoshift[len(blockstoshift)-1].split("-")[1],
       args.shift))

# the math magic, generating new list of shifted blocks
newblocks = []
for c in blockstoshift:
    pair = c.strip().split("-")
    # if len(pair) == 1:  # in case of 1-length blocks in mask, make a fake "pair" THIS IS PROBABLY NOW OBSOLETE
    #     pair.append(pair[0])
    #     print("Single-site block found, %s converted to %s.\n" % (pair[0], "-".join(pair)))
    for i, n in enumerate(pair):
        pair[i] = (int(n) + args.shift)  # convert to int var ; check what happens to int(num;)
    newstr = "%s-%s" % (pair[0], pair[1])
    newblocks.append(newstr)

finalblocks = blockstokeep + newblocks  # now combine with unshifted blocks

unchanged = " ".join([str(i) for i in blockstokeep])
changed = " ".join([str(i) for i in newblocks])
print("Final string of blocks: %s \033[1;91m%s\033[0m" % (unchanged, changed))  # colours only work in Unix environments

# finalise formatting for output
finalblocks = [str(i) for i in finalblocks]  # convert to strings
final = " ".join(finalblocks)+";\n"  # include final ;!
# print("Final string of blocks:\n%s\n" % final)  # kept for debugging

# write output!
newfile = str(args.nexus + ".shifted")
with open(newfile, "w") as outfile:
    with open(args.nexus, "r") as infile:
        for ln in infile:
            if ln[0:(searchlen+10)].strip("\t").lstrip().upper().startswith(searchstr):
                final = ln.split("=")[0] + "= " + final
                outfile.write(final)
            else:
                outfile.write(ln)

print("\nMask named '%s' in %s shifted by %s starting at %s and saved as %s\n" %
      (charset, args.nexus, args.shift, args.location, newfile))
