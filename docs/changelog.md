<!---
NOTE:
If you're reading this, instead try opening README.html in a web browser 
or view this file from within the github repository website.

This is a github-flavored markdown file not meant to be easily readable.
-->

Version history
---------------

### 2.1.5 (August 2019)
No future updates are planned beyond this final release.

- Fixed bug in reported median mutation rates in histogram figures 
  (previous ShapeMapper versions actually displayed the 5th percentile)
- Added guidance for running ShapeMapper components piecemeal (see
  [Modular workflow](modular_workflow.md))
- Added output histogram plot of un-normalized ln(mut_rate_modified/mut_rate_untreated)
- Reorganized list of ShapeMapper dependencies (see [Building](building.md))
- Adjusted <kbd>--min-mutation-separation 0</kbd> behavior to match expectations
- Added missing <kbd>--output-aligned-reads</kbd> option
- Allow <kbd>--input_is_unpaired</kbd> argument to shapemapper_mutation_parser
  to override flags present in SAM alignment
- Fixed issue with unit test paths on machines outside the build environment
- Fixed issue with undocumented --separate-ambig-counts param
- Exposed bowtie2 effort params <kbd>-R</kbd> and <kbd>-D</kbd> as <kbd>--max-reseed</kbd> 
  and <kbd>--max-search-depth</kbd>
- Added a more helpful error message in certain cases of missing input files
- Added error messages when attempting to run on a Mac or run when executables are not present
- Added a top-level CMakeLists to simplify most common build situation

### 2.1.4 (Oct. 2018)
- Documentation reorganized and expanded
- Fix for crash with folder names shorter than 3 characters
- [Amplicon primer pair mapping filters](analysis_steps.md#directed-primer-trimming)
- Now using paired-end alignment mode for paired reads that fail to merge
- Added test for overall pipeline accuracy on a small example dataset
- Added simplified reactivity profile outputs suitable for direct import into VARNA or Ribosketch (see
  [Coloring by SHAPE reactivity](other_software.md#coloring-by-shape-reactivity)).
- Refactored much of MutationCounter and MutationParser
   - Refactored Read class and used throughout
   - Moved all mutation processing and filtering functions from MutationCounter
     to MutationParser
   - Incorporated previously debug outputs into primary output of 
     MutationParser
- Debug mutation rendering (<kbd>--render-mutations</kbd>) reworked
   - Provides more detailed information about each mutation 
     processing step and quality filters
   - Outputs a multi-page pdf file scaled to fit the width specified
     by <kbd>--max-paired-fragment-length</kbd>
- Added end-to-end tests for unpaired inputs
- Bugfix for log file path when <kbd>--name</kbd> provided
- Exposed STAR <kbd>--genomeSAindexNbase</kbd> parameter
  - Added option to automatically rerun with defined <kbd>--genomeSAindexNbase</kbd> in the case of 
    STAR segfault (see [STAR parameters](analysis_steps.md#star-parameters))
- Excluded lowercase sequence from mutation rate histogram plots
- Added <kbd>--per-read-histograms</kbd> option and disabled by default
- Print all subprocess stdout/stderrs to main log file if run failed
  and <kbd>--verbose</kbd>

### 2.1.3
- Fix for crash with large <kbd>--random-primer-len</kbd>

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
- Error message if FASTQ files present in a provided <kbd>--folder</kbd> without
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
  <kbd>--render-flowchart</kbd>)
- Removed some large unused test files

### 0.1.5
- Quality control checks more readable
- Detailed quality control descriptions in <kbd>docs/quality_control.md</kbd>
- File format descriptions in <kbd>docs/file_formats.md</kbd>
- Softened quality control warning text
- Better <kbd>--version</kbd> and <kbd>--help</kbd> handling
- Files under active development conditionally excluded from tarball
- Changed default <kbd>--min-mutation-separation</kbd> to 6
- Added post-alignment basecall quality filter for mutation counting.
  Controlled by the <kbd>--min-qual-to-count</kbd> parameter (default=30). 
  Effective read depths are now calculated using only positions with 
  high-quality basecalls.
- Fixed speed problem with STAR end-to-end tests
- Fixed issue with <kbd>--correct-seq</kbd> crashing with regions of zero 
  coverage
- Fixed issue with <kbd>--correct-seq</kbd> causing downstream crash if sequence
  length changed.

### 0.1.4
- Bugfix release. Updated run_example.sh, log file location, STAR 
  aligner warning, and quality control warning for RNAs with long names.

### 0.1.3
- Multiple RNA support
    - Provide one or more FASTA files with one or multiple target 
      sequences in each file.
    - By default, normalize reactivity profiles as a group (disable
      with <kbd>--indiv-norm</kbd>)
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
  portion of reads with the <kbd>--random-primer-len</kbd> parameter
- STAR aligner support (<kbd>--star-aligner</kbd>). Recommended for large RNAs,
  as it can be much faster than Bowtie2.
- Intermediate/debug file output options
    - Aligned reads (<kbd>--output-aligned</kbd>)
    - Parsed mutations (<kbd>--output-parsed</kbd>)
    - Classified mutations (<kbd>--output-classified</kbd>)
    - Counted mutations (<kbd>--output-counted</kbd>)
    - Rendered mutations (<kbd>--render-mutations</kbd>). This will generate a 
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
  changed with <kbd>--min-mutation-separation</kbd>.
- Occluded depth correction support. Under this scheme, multinuc 
  mutations do not contribute to the read depth calculation, since they
  effectively prevent the detection of any modification within that 
  region. Disable with <kbd>--no-occluded-depth-correction</kbd>
- Options to exclude inserts, deletions, or ambiguously-aligned inserts
  or deletions from reactivity profile calculation
- Sped up some tests
- Added <kbd>--overwrite</kbd> option, otherwise give an error if existing files 
  conflict with output files
- Verbose option to show each subprocess command
- Default <kbd>--min-depth</kbd> raised to 5000
- Default <kbd>--min-mapq</kbd> lowered to 10
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


&nbsp;&nbsp;&nbsp;&nbsp;

[← back to README](../README.md)

