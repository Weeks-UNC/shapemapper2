<!---
NOTE:
If you're reading this, instead try opening README.html in a web browser 
or view this file from within the github repository website.

This is a github-flavored markdown file not meant to be easily readable.
-->

ShapeMapper file format descriptions
====================================

Input files
-----------

### FASTA file(s)
This file specifies reference (target) sequences. 
File extension should be .fa or .fasta. File should contain one or more
target DNA sequences ('`T`' not '`U`'). Sequences must not contain spaces or tabs, but
may be broken down into multiple lines. Each sequence must be preceded by a line with 
`>RNA_name`, where `RNA_name` is replaced with the name of the RNA of interest. 
Lowercase positions will be excluded from reactivity profiles, and should be used to indicate
primer-binding sites if using amplicon primers on either end of the sequence. 
If multiple primer pairs were used, or primer sites are not on the ends of the sequence, 
provide the primer sequences in a separate file with <kbd>--primers</kbd> (see below).

__Example:__

    >TPP_riboswitch
    ggccttcgggccaaggaCTCGGGGTGCCCTTCTGCGTGAAGGCTGAGAAATACCCGT
    ATCACCTGATCTGGATAATGCCAGCGTAGGGAAGTTCTCGATCCGGTTCGCCGGATC
    CAAATcgggcttcggtccggttc

### Primers file
Pairs of primers should be listed following the
name of each RNA preceded by '>'. RNA names in this file 
must match RNA names in provided target sequence .fa files.

    >RNA_name
    forward_primer_sequence reverse_primer_sequence

__Example:__

    >U1_snRNA
    ATACTTACCTGGCA CAGGGGAAAGCGCGAA

Each primer pair should be on its own line.
Multiple pairs can be provided, as can multiple RNAs.

### FASTQ file(s)
These files must have the extension .fastq or .fq, or .fastq.gz if they are compressed,
and must be FASTQ formatted.

If <kbd>--folder</kbd> is used to
pass a folder of FASTQ files to ShapeMapper, ShapeMapper will attempt to
identify and match up files corresponding to paired reads by 
finding '`R1`'|'`r1`' and '`R2`'|'`r2`' in
the filenames, separated by '`.`' or '`_`' characters.

Reads from separate instrument barcode indices must
be in separate files, and should not contain index sequences.



&nbsp;&nbsp;&nbsp;&nbsp;


Output files
------------

### `<name>_shapemapper_log.txt`

Run progress and summary outputs. Includes mate pair merging stats,
read alignment stats, reactivity profile quality control checks, and 
amplicon primer pair read depths.

&nbsp;&nbsp;&nbsp;&nbsp;

### `<name>_<RNA>_profile.txt`

Tab-delimited text columns. First line is column names.                     

| Column name               |  Content                                          |
|---------------------------|---------------------------------------------------|
|`Nucleotide`               |  Nucleotide number (1-based)                      |
|`Sequence`                 |  Nucleotide (AUGCaugc)                            |
|`<Sample>_mutations`       |  Mutation counts                                  |
|`<Sample>_read_depth`      |  Read depth                                       |
|`<Sample>_effective_depth` |  see [Effective read depth](analysis_steps.md#effective-read-depth) |
|`<Sample>_rate`            |  Effective mutation rate calculated as <br> <tt>(mutation count / effective read depth)</tt> |
|`<Sample>_off_target_mapped_depth` | Simple mapped read depths for reads <br> not meeting <kbd>--amplicon</kbd>/<kbd>--primer</kbd> location <br> requirements |
|`<Sample>_low_mapq_mapped_depth` | Simple mapped read depths for reads <br> not meeting <kbd>--min-mapq</kbd>  |
|`<Sample>_mapped_depth` <br> or `<Sample>_primer_pair_<n>_mapped_depth` | Simple mapped read depths for included <br> reads, broken down by primer pair if <br> applicable |
|`Reactivity_profile`       |  Calculated reactivity profile                    |
|`Std_err`                  |  Standard error                                   |
|`HQ_profile`               |  Reactivity profile with high-background and <br> low-depth positions excluded (set to nan)     |
|`HQ_stderr`                |  Standard error with high-background and <br> low-depth positions excluded |
|`Norm_profile`             |  Reactivity profile after normalization <br> (scaling) | 
|`Norm_stderr`              |  Standard error after normalization               |

&nbsp;&nbsp;&nbsp;&nbsp;

### `<name>_<RNA>.shape`

Reactivity profile in format expected by RNAstructure software. 
Tab-delimited text. Two columns. First column is 1-based nucleotide
position. Second is normalized reactivity, with excluded positions set to `-999`.

&nbsp;&nbsp;&nbsp;&nbsp;

### `<name>_<RNA>.map`

Same as .shape file, but with an additional two columns. Third column is
stderr, fourth is nucleotide sequence.

&nbsp;&nbsp;&nbsp;&nbsp;

### `<name>_<RNA>_varna_colors.txt` and `<name>_<RNA>_ribosketch_colors.txt`

Simplified reactivity profile suitable for import into VARNA or 
Ribosketch. Single column of normalized reactivity values. 
Reactivities above 0.85 are set to 0.85, reactivities below
0 are set to 0, and missing data positions are set to 0.

&nbsp;&nbsp;&nbsp;&nbsp;

### `<name>_<RNA>_profiles.pdf`

Figures showing read depths, mutation rates, and reactivity profile.

&nbsp;&nbsp;&nbsp;&nbsp;

### `<name>_<RNA>_histograms.pdf`

Figures showing read depth, mutation rate and reactivity histograms.

&nbsp;&nbsp;&nbsp;&nbsp;

### `<name>_<RNA>_mapped_depths.pdf`

Figures showing simple mapped read depths. Shows reads excluded due to low 
aligner-reported MAPQ (mapping quality score), and shows off-target reads excluded
due to not aligning near expected amplicon primer pair locations. 
Reads included in analysis are further broken down by primer pair. 
See [example plots](analysis_steps.md#mapped-depth-plots).

&nbsp;&nbsp;&nbsp;&nbsp;

### `<name>_<RNA>_per-amplicon_abundance.txt`

Nucleotide locations and maximum mapped read depths associated with each 
amplicon primer pair.



&nbsp;&nbsp;&nbsp;&nbsp;




Optional intermediate output files
----------------------------------

### Processed reads

Commandline option: <kbd>--output-processed-reads</kbd>
These files contain reads after initial quality trimming and paired
read merging steps are performed.

&nbsp;&nbsp;&nbsp;&nbsp;

### Aligned reads

Commandline option: <kbd>--output-aligned-reads</kbd>
Filename: `*_aligned.sam` if bowtie2 used or `*_aligned_paired.sam` and `*_aligned_unpaired.sam` if STAR used.
Format: SAM

&nbsp;&nbsp;&nbsp;&nbsp;

### Parsed mutations

Commandline option: <kbd>--output-parsed-mutations</kbd>

Filename: `<name>_<sample>_<RNA>_parsed.mut`

#### Format:

Text, one line per mapped read. Major fields are tab-delimited. Final
field contains internal space-delimited fields.

| Field                     |  Content                                          |
|---------------------------|---------------------------------------------------|
| 1 | read type (see below) |
| 2 | read name |
| 3 | 0-based leftmost mapping position (inclusive) |
| 4 | 0-based rightmost mapping position (inclusive) |
| 5 | read mapping category (see below)|
| 6 | primer pair index (0-based), <br> or -999 if no associated primers |
| 7 | mapped depth array (see below) |
| 8 | effective depth array (see below) |
| 9 | mutation count array (see below) |
| 10 | mutations (see below) |


#### _Read type:_

Read type is one of
- `PAIRED_R1`
- `PAIRED_R2`
- `UNPAIRED_R1`
- `UNPAIRED_R2`
- `UNPAIRED`
- `MERGED`
- `PAIRED`

#### _Read mapping category:_

Read mapping category is one of 
- `INCLUDED`
- `LOW_MAPQ`
- `OFF_TARGET`

Off-target and low-MAPQ reads are excluded from
analysis (with the exception of showing up as dashed lines
in mapped depth plots).

#### _Arrays:_

Arrays are <tt>right - left + 1</tt> characters long (that is,
the same length as the reference sequence over the aligned region), 
and contain only '`0`' and '`1`' characters.

_Mapped depths:_
Simple read coverage over the aligned region. For unpaired reads or
paired reads without a mapping mate pair, this array is full of '`1`'s.
For paired reads that do not overlap, this array will contain a region
of interior '`0`'s that indicate non-covered sequence between the mate
pairs.

_Effective read depths:_
Read coverage excluding primer regions, low-quality basecalls, and
the covered region of multinucleotide mutations excepting the inferred
adduct site. These arrays are summed to give the denominator used in calculating
mutation rate.

_Mutation counts:_
'`1`' indicates inferred adduct sites for a single read. These arrays are summed
to give the numerator used in calculating mutation rate.


#### _Mutations:_

Mutation fields are space-delimited in groups of five. Each group represents
a mutation from the target (reference) sequence.

| Field                     |  Content                                          |
|---------------------------|---------------------------------------------------|
| 1 | 0-based nearest unchanged target sequence position on the left |
| 2 | 0-based nearest unchanged target sequence position on the right |
| 3 | Double-quoted read sequence that replaces the target sequence <br> between the two positions. |
| 4 | Double-quoted basecall quality scores for each read position <br> in the previous field. |
| 5 | Double-quoted mutation classification (see below) |


#### _Mutation classifications:_

<tt>"A-","T-","G-","C-"</tt> (single-nucleotide deletions)

<tt>"-A","-T","-G","-C","-N"</tt> (single-nucleotide insertions)

<tt>"AT", "AG", "AC",
"TA", "TG", "TC",
"GA", "GT", "GC",
"CA", "CT", "CG" </tt>(single-nucleotide mismatches)

<tt>"multinuc_deletion"</tt>

<tt>"multinuc_insertion"</tt>

<tt>"multinuc_mismatch"</tt>

<tt>"complex_deletion"</tt>

<tt>"complex_insertion"</tt>

<tt>"N_match"</tt> (not a real mutation, just an ambiguous basecall)
    
All classifications may also have `_ambig` appended if a mutation involves 
any ambiguously aligned nucleotides.

&nbsp;&nbsp;&nbsp;&nbsp;



### Mutation counts

Commandline option: <kbd>--output-counted-mutations</kbd>

Filename: `<name>_<sample>_<RNA>_mutation_counts.txt`

#### Format: 

Tab-delimited text columns. First line is column headers. One column for
each mutation classification listed in the previous section, with the
exception of `N_match`. These columns contain mutation counts listed
5′ to 3′. 

Columns `read_depth` and `effective_depth` contain sequencing depths
(see [Effective read depth](analysis_steps.md#effective-read-depth)).

Columns `off_target_mapped_depth` and `low_mapq_mapped_depth` indicate
simple mapped read depths for reads excluded due to low aligner-reported 
MAPQ or due to failing to align near expected amplicon primer sites.

Columns `mapped_depth` or `primer_pair_<n>_mapped_depth` indicate
simple mapped read depths for all reads included in analysis, broken
down by associated amplicon primer pair if applicable.


&nbsp;&nbsp;&nbsp;&nbsp;

[← back to README](../README.md)
