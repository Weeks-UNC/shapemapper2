<!---
NOTE:
If you're reading this, instead try opening README.html in a web browser 
or view this file from within the github repository website.

This is a github-flavored markdown file not meant to be easily readable.
-->

Frequently asked questions
==========================

- [Why am I getting a warning about possible low-quality reactivity profiles?](#warning)
- [How should I handle primer trimming?](#trimming)
- [How do I model RNA secondary structures?](#how-do-i-model)
- [How do I render SHAPE reactivity-colored secondary structures?](#how-do-i-render)
- [Can I run ShapeMapper on multiple RNAs at the same time?](#multiple-rnas)
- [Can I run ShapeMapper on an existing alignment?](#existing-alignment)
- [What should I do if I am studying an RNA transcribed from multiple loci with slightly different sequences?](#multiple-loci)
- [Can I run ShapeMapper using reads from PacBio or other long-read sequencing platforms?](#long-reads)
- [How should I cite ShapeMapper?](#citation)

---

<a name="warning"></a>

### Why am I getting a warning about possible low-quality reactivity profiles?
The run has failed to meet one or more run quality checks.
Check the log file to see which check(s) failed, and refer to 
[Quality control checks](analysis_steps.md#quality-control-checks)
for explanations and possible solutions.

---

<a name="trimming"></a> 

### How should I handle primer trimming?

#### Random primers
If a Nextera or similar prep is used, no explicit primer
trimming should be necessary, since the ends of input cDNA
molecules will largely be digested away before sequencing.

Otherwise, the length of random primers used should be
provided to ShapeMapper using <kbd>--random-primer-len</kbd>.

#### Amplicon primers
For most users, simply set primer binding sites to lowercase sequence in the input
.fa file, and run ShapeMapper with the <kbd>--amplicon</kbd> option. This will handle
cases with a single pair of primers on either end of the target sequence.

For more complex cases with internal PCR primer locations or multiple pairs of primers, run 
ShapeMapper with a [primers file](file_formats.md#primers-file) 
input with 
<kbd>--primers <primers_file></kbd>
. ShapeMapper will automatically determine
the expected locations of input PCR primers by comparing their sequences with those of
reference targets.

Multiple amplicon primer pairs within a single dataset are also supported. Some groups
use a tiled amplicon strategy to cover low-abundance target RNAs.

For additional details on primer trimming and read location requirements, 
see [Primer trimming](analysis_steps.md#primer-trimming-and-enforcement-of-read-location-requirements).

---

<a name="how-do-i-model"></a>

### How do I model RNA secondary structures?
see [Other software](other_software.md#rna-secondary-structure-modeling)

---

<a name="how-do-i-render"></a>

### How do I render SHAPE reactivity-colored secondary structures?
see [Other software](other_software.md#traditional-secondary-structure-diagrams)

---

<a name="multiple-rnas"></a>

### Can I run ShapeMapper on multiple RNAs at the same time?
Yes. The <kbd>--target</kbd> parameter will accept multiple 
FASTA files, or multiple sequences can be included within
a single FASTA file.

However, inputting large numbers of sequences (such as those 
from whole-transcriptome studies) will cause ShapeMapper to 
run out of memory or crash. The [shapemapper-txome](https://github.com/Weeks-UNC/shapemapper-txome)
wrapper can be a useful workaround in this case, or see
[Modular workflow](modular_workflow.md) for help running alignment parsing,
mutation counting, and reactivity profile calculation steps
in isolation.

---

<a name="existing-alignment"></a>

### Can I run ShapeMapper on an existing alignment?

The standard shapemapper executable does not currently accept 
alignment files as input. However, individual shapemapper modules can be run
"manually" to process alignments in SAM format 
(see [Modular workflow](modular_workflow.md) for guidance).

---

<a name="multiple-loci"></a>

### What should I do if I am studying an RNA transcribed from multiple loci with slightly different sequences?
The presence of a subpopulation of sequence variants or transcribed pseudogenes 
is often initially indicated by an unusual number of nucleotides with high
apparent mutation rates in the untreated control sample.

If the sequences of contaminating/alternate transcripts are available, we recommend
providing all sequences to ShapeMapper (the <kbd>--target</kbd> FASTA file can have multiple 
RNAs within it, or multiple FASTA files can be provided). 
There is currently no good fully automated tool for determining the
sequences of contaminating transcripts - this would require a sequence assembly,
which remains an art, especially with noisy reads from reverse transcription. 

Care should be taken that identical sequences are not provided to ShapeMapper,
as this will result in reads aligning to multiple targets with similar alignment
scores and therefore failing to meet <kbd>--min-mapq</kbd>.

If sequences are not available, the <kbd>--render-mutations</kbd>
ShapeMapper option can be useful in some cases. This will render individual reads 
and their mutations (and debugging info that can be ignored) up to <kbd>--max-pages</kbd>. 
For short transcripts, this
allows quickly spotting frequent patterns of mutations that result from the presence
of contaminating transcripts. By default, the page size is scaled to include reads
up to 800 nucleotides long. This is often too zoomed out for convenient visualization,
so the <kbd>--max-paired-fragment-length</kbd> parameter may need to be lowered.

Warning: running ShapeMapper with more than 15 highly similar target sequences 
can exhaust bowtie2's search limits and result in mis-mapping across sequences. 
In this case, we recommend increasing the aligner's effort options
with <kbd>--max-search-depth</kbd> and perhaps <kbd>--max-reseed</kbd>, corresponding 
to bowtie2's -D and -R options. See
[bowtie2 documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#effort-options).

---

<a name="long-reads"></a>

### Can I run ShapeMapper using reads from PacBio or other long-read sequencing platforms?
This is currently unsupported (see [Read length](analysis_steps.md#read-length)).

---

<a name="citation"></a>

### How should I cite ShapeMapper?
#### For ShapeMapper2 software, please cite: 

Busan S, Weeks KM. Accurate detection of chemical modifications in RNA by mutational profiling (MaP) with ShapeMapper 2. _RNA_. 2018, 24(2):143-148.
[link](http://rnajournal.cshlp.org/content/early/2017/11/07/rna.061945.117)

#### For the MaP (mutational profiling) RNA adduct readout strategy, please cite either: 

Siegfried NA, Busan S, Rice GM, Nelson JA, Weeks KM. RNA motif discovery by SHAPE and mutational profiling (SHAPE-MaP). _Nat Methods_. 2014, 11(9):959-65.
[link](http://www.ncbi.nlm.nih.gov/pubmed/25028896)

Smola MJ, Rice GM, Busan S, Siegfried NA, Weeks KM. Selective 2'-hydroxyl acylation analyzed by primer extension and mutational profiling (SHAPE-MaP) for direct, versatile and accurate RNA structure analysis. _Nat Protoc_. 2015, 10(11):1643-69.
[link](http://www.ncbi.nlm.nih.gov/pubmed/21979276)

#### For DMS-specific analyses, please cite:

David Mitchell, III et al, Mutation signature filtering enables high-fidelity RNA structure probing at all four nucleobases with DMS, Nucleic Acids Research, 2023;, gkad522,
[link](https://doi.org/10.1093/nar/gkad522)

#### For msDMS_MaP (N7-G) analyses, please cite:

Irfana Saleem, Thomas Miller, Lucas Kearns, David Mitchell, Ritwika Bose, Chase Weidman, Anthony Mustoe. Title to be determined. Journal to be determined. 202X.

---

&nbsp;&nbsp;&nbsp;&nbsp;

[← back to README](../README.md)
