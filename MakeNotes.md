# *make* and *Snakemake* notes

Yana Eglit; 30 Nov 2021

### Installation

Document to cite for Snakemake: [Sustainable data analysis with Snakemake | F1000Research](https://f1000research.com/articles/10-33/v1)

Attempting snakemake. Installed on local Ubuntu via conda.

`make` itself relies on prespecified targets and dependencies.

        target : dependency

There are shortcuts (`$^` and `$@`) and variables `$(VAR)` but apparently Snakemake is much more intuitive for modern scientific computing. 

To install Snakemake, first use conda to install (apparently much faster, written in C++) [Mamba](https://github.com/mamba-org/mamba):

        `conda install mamba -n base -c conda-forge`

Then  run the following from the **base** conda environment:

        `mamba create -c conda-forge -c bioconda -n snakemake snakemake`

(this is dizzyingly fast compared to regular conda, savour the parallel installation!)

Then, of course, `source activate snakemake`

### Getting started

Preexisting [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/) and [The Snakemake Wrappers repository](https://snakemake-wrappers.readthedocs.io/en/stable/).

Snakemake uses bash 'strict mode' ([Bash Strict Mode](http://redsymbol.net/articles/unofficial-bash-strict-mode/)) by default!

##### Working through example [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) on [Gitpod](https://violet-turtle-kveq7nva.ws-us20.gitpod.io/)

###### Step 1:

[Basics: An example workflow](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) -- example using readmapping with bwa:

```yaml
#./Snakefile
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

`{input}` and `{output}` are lists of strings that are by default concatenated with a whitespace. I.e. here `{input}` becomes `data/genome.fa data/samples/A.fastq`.  A dry run `$ snakemake -np mapped_reads/A.bam`  (`-n` is dry run, `-p` is print shell commands that will be executed) produces:

    `bwa mem data/genome.fa data/samples/A.fastq | samtools view -Sb - > mapped_reads/A.bam`

Note that substituting `mapped_reads/B.bam` or dependency `data/genome.fa` leads to error as snakemake expects the **target**.

(note to self: use shell script to generate Snakefiles in batch?)

Snakemake defines 'job' as an application of a **rule** to generate a set of **output files** or targets.

Now execute the workflow -- the number of cores must always be specified:

    `$ snakemake --cores 1 mapped_reads/A.bam`

If run again, snakemake will not recreate A.bam, unless one of the unput files is newer than the latest output file. Rerunning the command anyway returns the following:

```text
(snakemake-tutorial) /workspace/snakemake-tutorial-data$ snakemake --cores 1 mapped_reads/A.bam
Building DAG of jobs...
Nothing to be done (all requested files are present and up to date).
Complete log: /workspace/snakemake-tutorial-data/.snakemake/log/2021-12-01T000335.641810.snakemake.log
```

###### Step 2: Generalising the rule

Rules can be generalised by applying the wildcard `{sample}`:

```yaml
#./Snakefile
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

You can have multiple wildcards, but if a rule has multiple output files, it must have the same wildcard to prevent two parallel jobs overwriting the same file.

Now executing `$ snakemake -np mapped_reads/B.bam` returns:

    `bwa mem data/genome.fa data/samples/B.fastq | samtools view -Sb - > mapped_reads/B.bam` 

Snakemake matches the B in `mapped_reads/B.bam` to replace `{sample}` in the shell command. Multiple targets can be specified as well: `$ snakemake -np mapped_reads/A.bam mapped_reads/B.bam`; or, even combined with Bash brace expansion (or other [Shell expansions](https://tldp.org/LDP/Bash-Beginners-Guide/html/sect_03_04.html)), for example:

`$ snakemake -np mapped_reads/{A,B}.bam` -- which will produce both options. Note: executing the latter only produces the shell command for B.bam as A.bam has already been created. You can `touch` A.bam to update timestamp and re-run the rule:

    `$ touch data/samples/A.fastq`

###### Step 3: Sorting read alignments

To sort alignments using `samtools`, add this block below `bwa_map` rule:

```yaml
# ./Snakefile (continued)
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```

Files from input directory taken, sorted, and saved as files in output directory.

Note: Snakemake creates new directories where any are missing!

Note `{wildcards.sample}`: this runs the wildcard **within** the shell command; this is probably because `{}` already means something in Bash itself, so Snakemake looks for the wildcards instruction? 

Running `snakemake -np sorted_reads/B.bam` now skips over the `bwa_map` rule as their targets are already up to date. Snakemake automatically determines dependencies via filenames.

###### Step 4: Indexing read alignments, and making graphs!

Now using samtools to index read alignments by genomic position:

```yaml
#./Snakefile (continued)
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
```

To look at interdependencies between jobs, a **directed acyclic graph** can be generated with:

    `$ snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg`

The `dot` command by Graphviz creates the visualisation (Snakemake's output is in .dot format, in standard output). Jobs that currently need not be run are represented by dashed frames. Example before and after running `samtools_sort` on B: 

<img src="file:///C:/Users/yegli/AppData/Roaming/marktext/images/2021-11-30-20-50-39-image.png" title="" alt="" width="219">    <img title="" src="file:///C:/Users/yegli/AppData/Roaming/marktext/images/2021-11-30-21-05-08-image.png" alt="" width="205">

###### Step 5: Calling genomic variants, and `expand()`

Expand() with samples takes inputs and formats into an array of strings with a specified name:

eg ```expand("sorted_reads/{sample}.bam", sample=SAMPLES)``` returns `["sorted_reads/A.bam", "sorted_reads/B.bam"]` under the name `SAMPLES`.

while `expand("sorted_reads/{sample}.{replicate}.bam", sample=SAMPLES, replicate=[0, 1])` iterates the substitutions for replicate over those of sample to produce: `["sorted_reads/A.0.bam", "sorted_reads/A.1.bam", "sorted_reads/B.0.bam", "sorted_reads/B.1.bam"]`

This is exceptionally useful especially since Snakemake works backward from requested output, not available input files like a for loop over \*.fas for example.

To the top of the snakefile add, in straight Python:

```python
SAMPLES = ["A", "B"]
```

And the following rule:

```yaml
# ./Snakefile (continued)
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```

Note that order of files is not preserved, and positional input files precede named ones (don't fully understand this yet...).

In the example above, multiple inputs can be specifically referred to via `{input.VARNAME}`, eg `{input.fa}`. Note the shell command is split into multiple indented lines, to be merged together by Python. (I think this requires that pipe... otherwise use `\`?) Input/output file lists can have any Python statements *as long as they produce strings*. Here, all three samples are iterated over each input file list.

    `$ snakemake --dag calls/all.vcf | dot -Tsvg > dag.svg` shows our graph now looks like this:

<img src="file:///C:/Users/yegli/AppData/Roaming/marktext/images/2021-12-01-00-39-53-image.png" title="" alt="" width="210">

###### Step 6: custom scripts

Snakemake has a `script:` directive. Added to end of Snakefile, to plot histograms:

```yaml
# ./Snakefile (continued)
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```

(wow, this does inputs and outputs for you instead of executing `python script.blah input output`) Script path is relative to Snakefile (pay attention to location!) . In this case, under `scripts/` create `plot-quals.py`:

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)

plt.savefig(snakemake.output[0])
```

This way, there is no need to include command line argument parsing code in the script!!! Note `snakemake.input[0]` (Suppose this also means the python script can be tailored to the Snakemake flow, eg positional arguments pre-defined in `input:` files list. Could be risky though...)

R likewise contains an S4 object `snakemake` , accessing the first element of input using `snakemake@input[[1]]` (note first element starts with 1 in R, not zero -- **this can be a source of error!**)  Named input/output can also be accessed via `snakemake@input[["myfile"]]` (not sure if this means if the Snakefile has` input: "myfile"` in the script rule...?) [More info here.](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-external-scripts)

###### Step 7: Adding a target rule

If the requested rule does not have wildcards, Snakemake accepts rulenames as targets. By default, if no target is specified, Snakemake defines the 1st rule as target. It is good practice to have a `rule all:` at the top of the Snakefile with the typically desired [final?] targets as input files. To top of Snakefile, add:

```yaml
# ./Snakefile (continued)
rule all:
    input:
        "plots/quals.svg"
```

Simply running `$ snakemake -n` will run the entire path of rules to generate the target plot from Step 6. Otherwise, *the order of rules is arbitrary in Snakemake*. One can also generate multiple target rules if a set of reasonable endpoints exists, which can be called via ``snakemake -n mytarget``. 

A dry run `-n` of Snakemake generates a planned job stats table; with `rule all:` it is summarised. (indeed, commenting out `rule all:` breaks this!)

<img title="" src="file:///C:/Users/yegli/AppData/Roaming/marktext/images/2021-12-01-01-11-53-image.png" alt="" width="342">

`snakemake --dag all | dot -Tsvg > dag.svg` leads to the full map:

<img title="" src="file:///C:/Users/yegli/AppData/Roaming/marktext/images/2021-12-01-01-15-07-image.png" alt="" width="186">

`snakemake --help` generates a fairly long helpfile. Eg:

```text
  --forcerun [TARGET ...], -R [TARGET ...]
                        Force the re-execution or creation of the given rules or files. Use
                        this option if you changed a rule and want to have all its output
                        in your workflow updated. (default: None)
```

`snakemake --cores 1 --forcerun samtools_sort` executes up to and including the rule `samtools_sort:`.  adding `--reason` adds a field explaining why this rule was run, eg:

```text
rule plot_quals:
    input: calls/all.vcf
    output: plots/quals.svg
    jobid: 1
    reason: Input files updated by another job: calls/all.vcf
    resources: tmpdir=/tmp
```

###### Summary: final workflow

```yaml
SAMPLES = ["A", "B"]


rule all:
    input:
        "plots/quals.svg"


rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}" #this is a continuation of above line


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```

#### Next Tutorial: [Advanced: Decorating the example workflow](https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html)

## Personal scraps

Shell block must only contain one line within quotes (accepts only 1 argument!)... or a quoteblock like so:

```yaml
shell:    
    """
    shell command 1
    shell command 2
    ...
    """
```

Commas after inputs and outputs important:

```yaml
input:
    file1,
    file2
output:
    file3
shell:
    "shell command 1"
```

The input line actually reads:

`input: file1, file2` 

When running Snakemake:

`-p` prints commands!

`-n` is a dry run

`-f` forces a run

Best way to run the whole pipeline from scratch is to `touch` the first input files 

Use the following function to extract shell scripts mentioned in file:

```shell
 function genSnakeShellScripts{
     sed -n '/^\#/! s/^.*[ "]\(.*\.sh\)[ $].*$/\1/p' Snakefile | uniq
 }
```


