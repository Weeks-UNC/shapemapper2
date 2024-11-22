<!---
NOTE:
If you're reading this, instead try opening README.html in a web browser 
or view this file from within the github repository website.

This is a github-flavored markdown file not meant to be easily readable.
-->

![](docs/images/header_profile.png)
**ShapeMapper2**
===============
*Copyright 2019 Steven Busan; 2024 Anthony Mustoe*. This project is licensed under the terms of the 
MIT license.

ShapeMapper automates the calculation of RNA chemical probing reactivities 
from mutational profiling (MaP) experiments, in which chemical adducts on RNA
are detected as internal mutations in cDNA through reverse transcription and 
read out by massively parallel sequencing. While originally built for analyzing SHAPE structure 
probing data, ShapeMapper is broadly useful for other types probing experiments, and includes
a DMS mode for analyzing DMS experiments. ShapeMapper performs 
- [Reference sequence correction](docs/analysis_steps.md#optional-reference-sequence-correction)
- [Read basecall quality trimming](docs/analysis_steps.md#initial-basecall-quality-trimming)
- [Paired read merging](docs/analysis_steps.md#paired-read-merging) (using [`BBmerge`](https://sourceforge.net/projects/bbmap/))
- [Alignment to reference sequences](docs/analysis_steps.md#alignment-to-reference-sequences) (using [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [`STAR`](https://github.com/alexdobin/STAR))
- [Enforcement of read location requirements and primer trimming](docs/analysis_steps.md#primer-trimming-and-enforcement-of-read-location-requirements)
- [Multinucleotide](docs/analysis_steps.md#multinucleotide-mutation-handling) and [ambiguously aligned](docs/analysis_steps.md#ambiguously-aligned-mutation-handling) mutation handling
- [Post-alignment basecall quality filtering](docs/analysis_steps.md#post-alignment-basecall-quality-filtering)
- [Chemical adduct location inference from detected mutations](docs/analysis_steps.md#chemical-adduct-location-inference)
- [Mutation rates calculation from mutation counts and effective read depths](docs/analysis_steps.md#calculation-of-mutation-rates)
- [Reactivity profile calculation and normalization](docs/analysis_steps.md#reactivity-profile-calculation-and-normalization)
- [Heuristic quality control checks](docs/analysis_steps.md#quality-control-checks)

<!---
FIXME?: these links don't seem to jump to the right place in the page for any headers after an image (works fine within page)
Update - seems to only be an issue in some browsers (like older Firefox on Linux). Looks like when opening
a link from a separate page, the browser can scroll to a location calculated before images are fully loaded.
-->

Installation
------------
ShapeMapper will only run on 64-bit Linux systems (Mac and Windows are not currently supported).

- Download latest [release](https://github.com/Weeks-UNC/shapemapper2/releases/download/2.2/shapemapper2-2.3.tar.gz)
    - On most systems, typing `wget https://github.com/Weeks-UNC/shapemapper2/releases/download/2.2/shapemapper2-2.3.tar.gz`
      will download the file on the commandline.
    - Be sure to download from the `shapemapper2-2.3.tar.gz` link, _not_ the source code-only links, which
      do not include executables.

- Extract release tarball using
    
    `tar -xvf shapemapper2-2.3.tar.gz`

- Add shapemapper executable to PATH (optional - google this if you don't know how)

- Run the script `run_example.sh` to check if shapemapper successfully 
  runs on a small subset of example data. (optional) 
    - This should produce two folders: `shapemapper_out` 
    and `shapemapper_temp`

- To run all unit and end-to-end tests, run `internals/test/run_all_tests.sh`. 
  This should take about 25-30 minutes. (optional) 

- Occassionally a single module failure detection test may fail. We attribute this error message to idiosyncracies in computational environment as opposed to an issue with shapemapper. This error message may be safely ignored.

- Alternatively, you can build ShapeMapper if the provided binaries do not run on your platform. 
  Building is relatively straightforward using conda. See <a href="docs/building.md">building</a>.

<!-- #### -->

Usage
-----
```
shapemapper <parameters> <inputs> | --version | --help
```

### Inputs
```
--<sample>   [--folder <fastq_folder> | --R1 <file_R1.fastq> --R2 <file_R2.fastq> |
              --unpaired-folder <fastq_folder> | --U <file.fastq> ]
             Samples must be from the following:
               "--modified", "--untreated/--unmodified", "--denatured"
             Folder or files may be specified, but not both.

             Note: reads from separate instrument barcode indices are expected to
             be in separate files, and should not contain index sequences.

--correct-seq [--folder <fastq_folder> | --R1 <file_R1.fastq> --R2 <file_R2.fastq> |
               --unpaired-folder <fastq_folder> | --U <file.fastq> ]
              Files to be used to identify sequence variants prior to SHAPE analysis.
              If a dedicated sequencing experiment is available, use those samples.
              In a typical MaP experiment, it is recommended to use the least-mutated 
              sample (untreated).
```

### Parameters

#### Required:
```
--target     FASTA file or list of files (.fa or .fasta) containing one or more
             target DNA sequences ('T' not 'U'). Lowercase positions will be
             excluded from reactivity profile, and should be used to indicate
             primer-binding sites if using directed primers. If multiple primer
             pairs were used, provide the primer sequences in a separate file with 
             '--primers' (see below).
```
#### Optional:
```
--dms        Run using DMS mode. Data will be normalized on a per-nucleotide basis.
             Note that this option is optimized for DMS-MaP data collected using 
             Marathon and/or TGIRT enzmyes. See docs/dmsmode.md for more details.

--name       Unique name to prefix all output filenames. Highly recommended if 
             running multiple ShapeMapper instances from the same folder.

--out        Output folder. Default="shapemapper_out"

--temp       Temporary file folder. Default="shapemapper_temp"

--overwrite  Overwrite existing files in output and temporary file folders
             without warning. Default=False

--log        Location of output log file. Default="<name>_shapemapper_log.txt"

--verbose    Display full commands for each executed process, and display more
             process output messages in the event of an error. Default=False

--random-primer-len <n>
             Length of random primers used (if any). Mutations within (length+1)
             of the 3-prime end of a read will be excluded, as will read depths
             over this region. Unused if '--amplicon' and/or '--primers' are 
             provided. Default=0

--amplicon   Require reads to align near expected primer pair locations, and
             intelligently trim primer sites. If a single pair of primers on 
             the ends of the RNA sequence is used, simply set primer sequences
             to lowercase in the '--target' fasta file. If multiple pairs or 
             internal locations are needed, specify primers with a '--primers'
             file.

     --primers <primers_file>
             Amplicon primer pair sequences. Each line should contain a pair of 
             primer sequences: the forward primer first followed by the reverse 
             primer, separated by whitespace. To specify primers for multiple RNAs,
             add a line with each RNA name preceded by '>' before each group of
             primer pairs. RNA names must match those in any provided .fa files.

     --max-primer-offset <n>
             If '--amplicon' and/or '--primers' used, require read ends to align to 
             within +/- this many nucleotides of expected amplicon primer pairs. 
             Default=10

--star-aligner
             Use STAR instead of Bowtie2 for sequence alignment. Recommended for
             sequences longer than several thousand nucleotides. Default=False
             Note: STAR slows down considerably in the presence of non-mapping
             sequences (i.e. if the target fasta files don't contain all the
             sequences present in the input reads). With current parameters, STAR 
             may also be slightly less sensitive than Bowtie2 (fewer aligned reads).

     --genomeSAindexNbase <n>
             Manually set STAR index building parameter. Default=0, indicating that
             ShapeMapper should recompute this parameter using the formula 
             recommended by the STAR manual.

     --rerun-on-star-segfault
             Automatically rerun ShapeMapper analyses that fail due to STAR segfault.
             The value of '--genomeSAindexNbase' will be replaced with the value of
             '--rerun-genomeSAindexNbase'. Default=False

     --rerun-genomeSAindexNbase <n>
             Default=3

     --star-shared-index
             Enable shared memory index. Default=False

--preserve-order
             Preserve the order of input reads through all analysis stages. May
             slow down execution, but can be useful for debugging. Default=False

--nproc <n>  Number of processors to use for sequence alignment (corresponding
             to bowtie2's '-p' parameter and STAR's '--runThreadN' parameter). 
             Default=4

--max-paired-fragment-length <n>
             Maximum distance between aligned ends of non-overlapping mate pairs 
             to be merged into a single read (analogous to bowtie2 '--maxins').
             Default=800

--max-search-depth <n>
             Set bowtie2 '-D' parameter. If negative, shapemapper calls bowtie2 
             with a default -D 15. Unused with --star-aligner.
             Default=-1
             
--max-reseed <n>
             Set bowtie2 '-R' parameter. If negative, shapemapper calls bowtie2
             with a default -R 2. Unused with --star-aligner.
             Default=-1

--min-depth  <n>
             Minimum effective sequencing depth for including data (threshold must 
             be met for all provided samples). Default=5000

--max-bg <n> 
             Maximum allowed mutation frequency in untreated sample. Default=0.05

--min-mapq <n>
             Minimum aligner-reported mapping quality for included reads. Default=10
             Note: When using Bowtie2, mutations contribute to lower mapping
             quality. Therefore, raising this threshold will have the side effect
             of excluding highly mutated reads.
             Note: This option does not apply to sequence correction, which uses
             a threshold of 10 regardless of this option

--min-qual-to-trim <n>
             Minimum phred score in initial basecall quality trimming. 
             Default=20

--window-to-trim <n>
             Window size in initial basecall quality trimming. Default=5

--min-length-to-trim
             Minimum trimmed read length in initial basecall quality trimming.
             Default=25

--min-qual-to-count <n>
             Only count mutations with all basecall quality scores meeting this
             minimum score (including the immediate upstream and downstream 
             basecalls). This threshold is also used when computing the 
             effective read depth. Default=30

--indiv-norm Normalize multiple reactivity profiles individually, instead of as a
             group. Default=False

--min-seq-depth <n>
             Minimum sequencing depth for making a sequence correction (with
             '--correct-seq'). Default=50

--min-freq <n>
             Minimum mutation frequency for making a sequence correction (with
             '--correct-seq'). Default=0.6

--disable-soft-clipping
             Disable soft-clipping (i.e. perform end-to-end rather than local 
             alignment). Default=False
             Note: this does not apply to sequence correction, which uses 
             soft-clipping regardless.

--right-align-ambig
             Realign ambiguous deletions/insertions to their rightmost valid position
             instead of leftmost. Not recommended, since left-realignment produces
             empirically better reactivity profiles than right-realignment.
             Default=False

--min-mutation-separation <n>
             For two mutations to be treated as distinct, they must be separated by at
             least this many unchanged reference sequence nucleotides. Otherwise, they
             will be merged and treated as a single mutation. Does not apply to 
             sequence correction. Default=6

--output-processed-reads
--output-aligned-reads
--output-parsed-mutations
--output-counted-mutations
             Produce output files for selected intermediate components. Default=False

--render-flowchart
             Render a flowchart (SVG format) in the output folder. This will depict 
             all data processing components and input and output files for the 
             current analysis pipeline. Default=False

--render-mutations
             Render pdf files showing detailed read and mutation processing steps 
             for each sample and RNA target, up to '--max-pages'. Primarily a debugging
             tool, but can be useful to visually inspect individual reads for the
             presence of pseudogenes.

     --max-pages <n>
             Maximum pages to render for '--render-mutations'. Default=100
             
     --render-must-span <n>-<n>
            Only render reads that cover a given nucleotide range. Disabled by default

--per-read-histograms
             Output read length and per-read mutation frequency histogram tables in 
             log file.

--serial     Run pipeline components one at a time and write all intermediate files
             to disk. Useful for debugging, but not generally recommended, as this will 
             use large amounts of disk space. Default=False

--N7         Add N7 information to data visualization. Adds a graph of mutation rates
             and reactivities specific to N7 data in profiles.pdf. Prior to usage, 
             ensure proper protocol was followed to generate valid N7 data.

--output-temp     
             Preserves temp files. Default=False.
             
--pernt-norm-factor-threshold
             Set the number of NTs needed for effective per-nt normalization factor
             calculation. May need to change in the case of short RNAs. Default=20

--ignore_low_N7
             Bypass N7 quality control filters.

--bypass_filters
             Bypass N7 quality control filters and set threshold for NTs needed for effective
             per-nt normalization factor calculation to 1.
             (Equivalent to "--ignore_low_N7 --theshold 1")
```

&nbsp;&nbsp;&nbsp;&nbsp;

Usage examples
--------------
(Note: commandline argument examples only; will not produce output. 
 For a runnable example, execute `run_example.sh`)

&nbsp;&nbsp;&nbsp;&nbsp;

Three-sample experiment, input FASTQ files:

  ``shapemapper --name example --target TPP.fa --out TPP_shapemap --amplicon --modified --R1 TPPplus_R1.fastq.gz --R2 TPPplus_R2.fastq.gz --untreated --R1 TPPminus_R1.fastq.gz --R2 TPPminus_R2.fastq.gz --denatured --R1 TPPdenat_R1.fastq.gz --R2 TPPdenat_R2.fastq.gz``

&nbsp;&nbsp;&nbsp;&nbsp;

Two-sample experiment, input from folders:

  ``shapemapper --name example2 --target TPP.fa --out TPP_shapemap --amplicon --modified --folder TPPplus --untreated --folder TPPminus``

&nbsp;&nbsp;&nbsp;&nbsp;

Only generate corrected sequence:

  ``shapemapper --name example3 --target TPP.fa --out TPP_mutant --amplicon --correct-seq --folder sequence_variant/A100``

&nbsp;&nbsp;&nbsp;&nbsp;

Generate corrected sequence using untreated sample,
then perform SHAPE-MaP analysis:

  ``shapemapper --name example4 --target TPP.fa --out --amplicon TPP_mutant --correct-seq --folder TPPminus --modified --folder TPPplus --untreated --folder TPPminus --denatured --folder TPPdenat``

&nbsp;&nbsp;&nbsp;&nbsp;

Multiple RNAs, randomly-primed experiment, STAR aligner:

  ``shapemapper --name example5 --target 16S.fa 23S.fa --out ribosome --random-primer-len 9 --star-aligner --modified --folder ribosome_plus --untreated --folder ribosome_minus --denatured --folder ribosome_denat``

&nbsp;&nbsp;&nbsp;&nbsp;

Process single DMS modified sample using DMS mode:

``shapemapper --name example6 --target add.fa --out add_dms --dms --amplicon --modified --folder ribosome_plus``

<!-- #### -->&nbsp;&nbsp;&nbsp;&nbsp;


Additional documentation
------------------------

### Frequently asked questions
see [FAQ](docs/FAQ.md)

### DMS mode
see [DMSmode](docs/dmsmode.md)

### N7-G related functionality
see [N7-G](docs/N7-G.md)

### Low-quality profile warning message
If ShapeMapper gives a red warning message about possible low-quality
reactivity profiles, read the log file to see which quality control
checks failed, and refer to 
[Quality control checks](docs/analysis_steps.md#quality-control-checks)
 for possible remedies.

### Modeling RNA structure
ShapeMapper does not perform RNA structure modeling. 
See [Other software](docs/other_software.md).

### Analysis steps
see [Analysis steps](docs/analysis_steps.md)

### File formats
see [File formats](docs/file_formats.md)

### Dependencies and build requirements
All third-party executables and compiled executables are included in 
the release. These should be compatible with most 64-bit Linux platforms,
even fairly old ones.

In the rare case that a rebuild is necessary, see [Building](docs/building.md)

### Modular workflow
For guidance running components of ShapeMapper in isolation,
see [Modular workflow](docs/modular_workflow.md)

### Version history
see [Version history](docs/changelog.md)


&nbsp;&nbsp;&nbsp;&nbsp;


Citation
--------

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

&nbsp;&nbsp;&nbsp;&nbsp;
