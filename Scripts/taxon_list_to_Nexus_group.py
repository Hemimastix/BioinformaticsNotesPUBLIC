# Converts list of taxa into a NEXUS taxonomic group recognised by SeaView.
# Example usage: Load tree corresponding to alignment in FigTree, select taxa as desired and paste them into text file, then run this script.
# By Yana Eglit, 05 June 2024

import pandas as pd
import sys
import os
import re

nexfile = sys.argv[1]
tax_list = sys.argv[2]
taxset_name = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(tax_list).split(".")[0]

usage = """
USAGE:
Converts list of taxa into a NEXUS taxonomic group recognised by SeaView. 
Example use: Load tree corresponding to alignment in FigTree, select taxa as desired and paste them into text file, 
then run this script.

Run this script with the following arguments: [nexus file] [list of taxa] [(optional) name for taxon set]
nexus file -- alignment file to apply taxon groups to
list of taxa -- plaintext file with taxon names exactly as they appear in the Nexus alignment (without single quotes),
separated by newlines.
name for taxon set -- defaults to base portion of list of taxa file, but can be manually specified. I'd avoid spaces.

Should save output modified Nexus file with .newgroup.nxs as extension.
Note: additionally creates an inverse taxon set called [name]_inverse
"""

print(usage)

try:
    #load nexus file
    with open(nexfile, "r") as f:
        nexus = f.readlines()
except FileNotFoundError:
    print(f"Nexus file {nexfile} not found")

try:
    #load taxon list
    with open(tax_list, "r") as f:
        taxa = f.readlines()
except FileNotFoundError:
    print(f"Taxon list file {tax_list} not found")

print(f"Loaded {nexfile} as Nexus file and {tax_list} as list of taxa\n")

nexus_taxblock = []
for line in nexus:
    if re.match('\[[0-9]*\] ', line):
        nexus_taxblock.append(line.replace("'", "").strip("\n").replace(" ", "\t"))
for i, l in enumerate(nexus_taxblock):
    nexus_taxblock[i] = l.split("\t")

nexus_taxa = pd.DataFrame(nexus_taxblock)
print(f"Detected {len(nexus_taxa)} taxa in Nexus file")

# parse taxa list; remove trailing newlines and single quotation marks
numerical_list = []
for taxon in taxa:
    t = taxon.strip("\n").replace("'","")
    if t in nexus_taxa[1].values:
        match = nexus_taxa[nexus_taxa[1] == t].values[0]
        # print(match)  # for debugging TODO remove
        match_num = match[0].strip("[]")
        if match_num in numerical_list:
            print(f"{match} number {match_num} an apparent duplicate, ignoring")
        else:
            numerical_list.append(int(match_num))  # convert to integer from string, needed for sort
    else:
        print(f"Warning: {t} not found among taxa!")

if len(numerical_list) == 0:
    print("Numerical list is empty, no matches found?")
    exit()
else:
    print(f"Detected {len(numerical_list)} taxa in list")

# making inverse taxon list
inverse_numerical_list=[]
for num in nexus_taxa[0].values:
    if int(num.strip("[]")) not in numerical_list:
        inverse_numerical_list.append(int(num.strip("[]")))
print(f"Also generating inverse TAXSET with {len(inverse_numerical_list)} taxa\n")

numerical_list.sort()
for i, n in enumerate(numerical_list):
    numerical_list[i] = str(n)

inverse_numerical_list.sort()
for i, n in enumerate(inverse_numerical_list):
    inverse_numerical_list[i] = str(n)

# for n in numerical_list:
#     " ".join()

taxset = f"  TAXSET '{taxset_name}' = "+" ".join(numerical_list)+";"
print(f"Taxset line to add: {taxset}")

inv_taxset = f"  TAXSET '{taxset_name}_inverse' = "+" ".join(inverse_numerical_list)+";"
print(f"Taxset line to add: {inv_taxset}")


# insert taxset:
for i, line in enumerate(nexus):
    if line.startswith("BEGIN SETS;"):
        nexus[i] = f"{line}{taxset}\n{inv_taxset}\n"
print(f"{nexfile}")
with open(f"{os.path.basename(nexfile).split('.')[0]}.newgroup.nxs", "w") as f:  # TODO: fix this, truncates after first .
    f.writelines(nexus)
