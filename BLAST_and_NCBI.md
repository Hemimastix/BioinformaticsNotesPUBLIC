




### SRA download

    wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR179/SRR1793207/*

    fastq-dump -I --split-files  ./SRR1793207.sra

use sratoolkit (./opt/perun/sratoolkit/)'s fastq-dump to convert .sra to fastq. For paired end reads (generates two files):

\***note: not all SRA data is paired =) (in that case, omit --split-files)


