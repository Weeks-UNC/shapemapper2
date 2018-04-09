
**ShapeMapper**
===============
*Copyright 2017 Steven Busan*. This project is licensed under the terms of the 
MIT license.

ShapeMapper automates the calculation of RNA structure probing reactivities 
from mutational profiling (MaP) experiments, in which chemical adducts on RNA
are detected as internal mutations in cDNA through reverse transcription and 
read out using massively parallel sequencing. Using custom and third-party 
components, ShapeMapper performs 
- Read basecall quality trimming
- Paired read merging (with `BBmerge`)
- Alignment to reference sequences (with `bowtie2` or `STAR`) 
- Optional reference sequence correction
- Multinucleotide and ambiguously aligned mutation handling
- Chemical adduct location inference from detected mutations
- Post-alignment basecall quality filtering
- Mutation rates calculation from mutation counts and effective read depths
- Reactivity profile calculation and normalization
- Heuristic quality control checks

Publications:

[Accurate detection of chemical modifications in RNA by mutational
 profiling (MaP) with ShapeMapper 2 _RNA_. 2017](http://rnajournal.cshlp.org/content/early/2017/11/07/rna.061945.117)

[RNA motif discovery by SHAPE and mutational profiling (SHAPE-MaP). _Nat 
 Methods_. 2014](http://www.ncbi.nlm.nih.gov/pubmed/25028896)

[Selective 2'-hydroxyl acylation analyzed by primer extension and 
 mutational profiling (SHAPE-MaP) for direct, versatile and accurate RNA 
 structure analysis. _Nat Protoc_. 2015](http://www.ncbi.nlm.nih.gov/pubmed/21979276)


Installation
------------
- Extract release tarball using
    
    `tar -xvf shapemapper-2.1.3.tar.gz`

- Add shapemapper executable to PATH (optional)

- Run the script `run_example.sh` to check if shapemapper successfully 
  runs on a small subset of example data. (optional) 
    - This should produce two folders: `example_data/shapemapper_out` 
    and `example_data/shapemapper_temp`

- To run all unit and end-to-end tests, run `test/run_all_tests.sh`. 
  This should take about 5-10 minutes. (optional) 

- See "Dependencies" section for how to install third-party 
  components if not included in release.

- See "C++ build requirements" and "Building C++ modules" sections 
  for how to rebuild shapemapper components if the provided binaries 
  do not run on your platform. 


Usage
-----
```
shapemapper <global params> <sample-specific params> | --version | --help

Sample-specific params
----------------------
--<sample>   [--folder <fastq_folder> | --R1 <file_R1.fastq> --R2 <file_R2.fastq> |
              --unpaired-folder <fastq_folder> | --U <file.fastq> ]
             Sample must be one of the following:
               "modified", "untreated/unmodified", "denatured"
             Folder or files may be specified, but not both.

             Note: reads from separate instrument barcode indices are expected to
             be in separate files, and should not contain index sequences.

--correct-seq [--folder <fastq_folder> | --R1 <file_R1.fastq> --R2 <file_R2.fastq> |
               --unpaired-folder <fastq_fodler> | --U <file.fastq> ]
              Files to be used to identify sequence variants prior to SHAPE analysis.
              If a dedicated sequencing experiment is available, use those samples.
              In a typical MaP experiment, it is recommended to use the least-mutated 
              sample (untreated).

Global params
-------------
Required:
--target     FASTA file or list of files (.fa or .fasta) containing one or more
             target DNA sequences ('T' not 'U'). Lowercase positions will be
             excluded from reactivity profile, and should be used to exclude
             primer-binding sites if using non-random primers.

Optional:
--name       Unique name to prefix all output filenames. Highly recommended if 
             running multiple ShapeMapper instances from the same folder.
--out        Output folder. Default="shapemapper_out"
--temp       Temporary file folder. Default="shapemapper_temp"
--overwrite  Overwrite existing files in output and temporary file folders
             without warning. Default=False
--log        Location of output log file. Default="<name>_shapemapper_log.txt"
--verbose    Display full commands for each executed process. Default=False

--random-primer-len
             Length of random primers used (if any). Mutations within (length+1)
             of the 3-prime end of a read will be excluded. This region will also 
             be excluded from the effective read depth. Default=0

--star-aligner
             Use STAR instead of Bowtie2 for sequence alignment. Recommended for
             sequences longer than several thousand nucleotides. Default=False
             Note: STAR slows down considerably in the presence of non-mapping
             sequences (i.e. if the target fasta files don't contain all the
             sequences present in the input reads).

--nproc      Number of processors to use for sequence alignment (corresponds
             to Bowtie2's "-p" parameter and STAR's "runThreadN" parameter). 
             Default=4

--min-depth  Minimum effective sequencing depth for including data (threshold must 
             be met for all provided samples). Default=5000
--max-bg     Maximum allowed mutation frequency in untreated sample. Default=0.05
--min-mapq   Minimum reported mapping quality for including reads. Default=10
             Note: When using Bowtie2, mutations contribute to lower mapping
             quality. Therefore, raising this threshold will have the side effect
             of excluding highly mutated reads.
             Note: This option does not apply to sequence correction, which uses
             a threshold of 10 regardless of this option
--min-qual-to-count
             Only count mutations with all basecall quality scores meeting this
             minimum score (including the immediate upstream and downstream 
             basecalls). This threshold is also used when computing the 
             effective read depth. Default=30
--indiv-norm Normalize multiple reactivity profiles individually, instead of as a
             group. Default=False
--min-seq-depth
             Minimum sequencing depth for making a sequence correction (with
             --correct-seq). Default=50
--min-freq   Minimum mutation frequency for making a sequence correction (with
             --correct-seq). Default=0.6
--disable-soft-clipping
             Disable soft-clipping (i.e. perform end-to-end rather than local 
             alignment). Default=False
             Note: this option does not apply to sequence correction, which uses
             soft-clipping regardless of this option.
--right-align-ambig
             Realign ambiguous deletions/insertions to their rightmost valid position
             instead of leftmost. Not recommended, since left-realignment produces
             empirically better reactivity profiles than right-realignment.
             Default=False
--min-mutation-separation
             For two mutations to be treated as distinct, they must be separated by at
             least this many unchanged reference sequence nucleotides. Otherwise, they
             will be merged and treated as a single mutation. Default=6

--output-processed | --output-processed-reads
--output-aligned
--output-parsed | --output-parsed-mutations
--output-counted | --output-counted-mutations
--output-classified | --output-classified-mutations
             Produce output files for selected intermediate components. Default=False

--render-flowchart
             Render a flowchart (SVG format) in the temp file output folder. This will
             depict all data processing components and input and output files for the
             current analysis pipeline. Default=False

--structured-output
             Output all files in a folder hierarchy matching the pipeline organization.
             If this option is provided, "--temp" folder will not be used.
             Default behavior is to generate files and folders in a hierarchy in
             the temp folder, except for important output files which are
             generated in the main output folder.

--serial     Run pipeline components one at a time and write all intermediate files
             to disk. Not recommended, as this will use large amounts of
             disk space. Default=False

Examples
--------
Three-sample experiment, input FASTQ files:

  shapemapper --name example --target TPP.fa --out TPP_shapemap --modified --R1 TPPplus_R1.fastq.gz --R2 TPPplus_R2.fastq.gz --untreated --R1 TPPminus_R1.fastq.gz --R2 TPPminus_R2.fastq.gz --denatured --R1 TPPdenat_R1.fastq.gz --R2 TPPdenat_R2.fastq.gz

Two-sample experiment, input from folders:

  shapemapper --name example2 --target TPP.fa --out TPP_shapemap --modified --folder TPPplus --untreated --folder TPPminus

Only generate corrected sequence:

  shapemapper --name example3 --target TPP.fa --out TPP_mutant --correct-seq --folder sequence_variant/A100

Generate corrected sequence using untreated sample,
then perform SHAPE-MaP analysis:

  shapemapper --name example4 --target TPP.fa --out TPP_mutant --correct-seq --folder TPPminus --modified --folder TPPplus --untreated --folder TPPminus --denatured --folder TPPdenat

Multiple RNAs, randomly-primed experiment, STAR aligner:

  shapemapper --name example5 --target 16S.fa 23S.fa --out ribosome --modified --folder ribosome_plus --untreated --folder ribosome_minus --denatured --folder ribosome_denat --random-primer-len 9 --star-aligner



```

Quality control checks
----------------------
- Documentation at `docs/quality_control.md`


File formats
------------
- Documentation at `docs/file_formats.md`


Dependencies
------------
All third-party executables are included in the release, within
the subdirectory `thirdparty`. Deleting this folder and running
`install/build_thirdparty.sh` will regenerate these files (requires 
an active Internet connection).

- Bowtie2>=2.2.7, STAR>=2.5.2a
- Python>=3.5, matplotlib>=1.5.1, numpy
- BBmerge>=35.85 (requires java)
- Pipe Viewer
- Graphviz


Source code documentation
-------------------------
- C++ documentation at `cpp-src/doc/html/index.html`
- Python documentation in source code comments (folder `python`)


C++ build requirements
----------------------
- recent `cmake` (>=3.3)
- `gcc` compiler with C++11 support
- `git` with https support
- `build-essential` (should include `dpkg-dev`, `g++`, `libc-dev`, `libstdc++`)
- `libboost-dev` (>=1.60.0)
- `libboost-filesystem-dev`
- `libboost-iostreams-dev`
- `libboost-program-options-dev`
- `doxygen` (optional)
- active internet connection for getting bamtools source


Building C++ modules
--------------------

To rebuild ShapeMapper executables
- Delete the folders `build` and `build_deps` if present
- Run `install/build_binaries.sh --use-system-libs` 
- This assumes the tools and libraries listed above are already 
  present at the system level.

Alternatively, run `install/build_all.sh`
- This will take substantial disk space and time
- Compiler and library dependencies listed above (with the exception
  of Doxygen) will be downloaded and built locally within the 
  `build_deps` folder


C++ code organization notes
---------------------------
- Generally using namespaces to modularize instead of public static classes.
 
- Related functions are in a single namespace in a similarly named `.cpp` file, 
with utility functions in a `detail` nested namespace in a `.h` file. Commandline 
executable versions of modules currently located in e.g. `ReadTrimmerExe.cpp`.


History
-------

### 2.1.3
- Fix for crash with large `--random-primer-len`

### 2.1.2
- Fix for execution in Slurm cluster environments
- Added error message for all-lowercase input sequence
- Third-party conda package fixes
- Bamtools CMake fixes to accommodate recent repo changes

### 2.1.1
- Added simple read length and mutations per read histogram outputs
- Support passing target sequences directly on commandline
- Various fixes to ease remote builds
- Exit with error if tests fail
- Error message if FASTQ files present in a provided `--folder` without
  R1 or R2 in filename
- Bugfix for sequence names with spaces
- Softened data quality warning message
- Added more detailed output for intermediate single read classified 
  mutation files.

### 2.0-rc3
- Added grip-rendered README.

### 2.0-rc2
- Added license
- Updated README
- Allow pipeline to run to completion even if some RNAs have no mapped
  reads
- Updated STAR aligner suggestion message for long RNAs
- Updated thirdparty package management scripts
- Moved some utility scripts to separate repo
- Bugfix for partial argument parsing
- Fix for bowtie2 component failure detection
- Clarified error message for filename collisions
- Do not render pipeline flowchart by default (controlled with 
  `--render-flowchart`)
- Removed some large unused test files

### 0.1.5
- Quality control checks more readable
- Detailed quality control descriptions in `docs/quality_control.md`
- File format descriptions in `docs/file_formats.md`
- Softened quality control warning text
- Better `--version` and `--help` handling
- Files under active development conditionally excluded from tarball
- Changed default `--min-mutation-separation` to 6
- Added post-alignment basecall quality filter for mutation counting.
  Controlled by the `--min-qual-to-count` parameter (default=30). 
  Effective read depths are now calculated using only positions with 
  high-quality basecalls.
- Fixed speed problem with STAR end-to-end tests
- Fixed issue with `--correct-seq` crashing with regions of zero 
  coverage
- Fixed issue with `--correct-seq` causing downstream crash if sequence
  length changed.

### 0.1.4
- Bugfix release. Updated run_example.sh, log file location, STAR 
  aligner warning, and quality control warning for RNAs with long names.

### 0.1.3
- Multiple RNA support
    - Provide one or more FASTA files with one or multiple target 
      sequences in each file.
    - By default, normalize reactivity profiles as a group (disable
      with `--indiv-norm`)
- Unpaired read support
- Masked region support 
    - Lowercase nucleotides in sequence will be excluded from reactivity
      profile calculation.
    - Useful for primer binding regions in targeted primer experiments
- Quality control checks and warnings
    - Good read depths
    - Mutation rates higher in modified sample than in untreated
    - Expected number of highly reactive nucleotides
    - Not too many high background positions
- Exclude mutations and depths over 3-prime random primer binding
  portion of reads with the `--random-primer-len` parameter
- STAR aligner support (`--star-aligner`). Recommended for large RNAs,
  as it can be much faster than Bowtie2.
- Intermediate/debug file output options
    - Aligned reads (`--output-aligned`)
    - Parsed mutations (`--output-parsed`)
    - Classified mutations (`--output-classified`)
    - Counted mutations (`--output-counted`)
    - Rendered mutations (`--render-mutations`). This will generate a 
      postscript image showing a subset of reads with parsed mutations
      and the adjusted mutations that ultimately contribute to profile
      calculation shown above and below the read, respectively.
- Realign ambiguously-located deletions and insertions to their leftmost
  valid positions and include in reactivity profile calculation. This
  empirically produces more accurate profiles than excluding ambiguous
  mutations.
- Combine mutations separated by up to 11 unchanged reference
  nucleotides. This empirically produces more accurate profiles than
  only combining immediately adjacent mutations. This threshold can be
  changed with `--min-mutation-separation`.
- Occluded depth correction support. Under this scheme, multinuc 
  mutations do not contribute to the read depth calculation, since they
  effectively prevent the detection of any modification within that 
  region. Disable with `--no-occluded-depth-correction`
- Options to exclude inserts, deletions, or ambiguously-aligned inserts
  or deletions from reactivity profile calculation
- Sped up some tests
- Added `--overwrite` option, otherwise give an error if existing files 
  conflict with output files
- Verbose option to show each subprocess command
- Default `--min-depth` raised to 5000
- Default `--min-mapq` lowered to 10
- Flowchart legend
- Simplified pipeline-building framework
- Clearer output filenames
- Minor fixes
    - Counted mutations and read depth files are now guaranteed to have 
      the same lengths, even if the 3-prime end of the RNA is covered
      by no reads
    - Error message when denatured and modified samples provided but
      no untreated sample
    - Profile normalization creates a new file rather than overwriting
      input file
    - Error for FASTA file with whitespace in sequence
    - Explicit filename argument handling
    - Sequence variants reported to user in 1-based coordinates
    - Various test runner fixes
    - Fix for build using local libs

### 0.1.2
- Histogram rendering (mutation rates, sequencing depths, reactivities)
- Sequence variant correction integrated into command-line interface
    - Generate updated FASTA file
    - Report sequence changes to user
    - Warn user about high-frequency but sub-threshold mutations
- Alignment stats reported to user
    - Spurious reads making it through BBmerge are now filtered out, so 
      alignment stats better reflect actual mapping percentages
- Option to generate output without excluding ambiguously aligned
  mutations
    - Histograms generated this way should better reflect actual 
      mutation rates
    - Nucleotide-resolution profiles generated this way may be 
      misleading in regions of homopolymeric or repeated sequence
- Log file generation
- Testing
    - Unit tests built again
    - End-to-end pipeline success tests
    - End-to-end module failure detection tests
    - Sequence variant correction tests
    - Reduced size of test dataset to speed up execution
    - Single script to run all tests and summarize results
- Serial execution working again
- Mapping quality filter
- FASTA format checks
- Segfaults reported to user
- Bugfix for 2-sample run
- Correct handling of ambiguous mutations for poor alignments (corner case)
- Flowchart rendering working on UNC cluster
- Quote filenames in shell wrappers
- CMake rebuild much faster (don't rebuild BamTools every time)
- Misc. fixes to build scripts


TODO
----
- Primer site filters
- Clearer syntax for accessing split-to-file wrapper outputs
- Stranded support
- Move sequence variant messages to within overall QC message
- Dump all logs to single file on failure in verbose mode?

&nbsp;
