<!---
NOTE:
If you're reading this, instead try opening README.html in a web browser 
or view this file from within the github repository website.

This is a github-flavored markdown file not meant to be easily readable.
-->

Modular workflow
================

General info
------------

Normal ShapeMapper execution automatically performs a series of 
[analysis steps](analysis_steps.md) internally. This is usually the
most convenient option for general users. However, some users may wish
to perform individual analyses semi-manually in isolation, 
allowing sequence alignment with a different aligner, 
or the introduction of specialized filtering stages, or the inclusion of mutation 
processing steps within a larger workflow.

Under the hood, ShapeMapper is fairly modular. This page attempts
to summarize the main steps power users will be most interested in.
Core shapemapper binary executables and scripts are located in `internals/bin`.
Python scripts are designed for python>=3.5. Thirdparty executables are in various
locations within `internals/thirdparty`.

To temporarily put bundled binaries and thirdparty executables 
into the active shell PATH, run `source internals/paths/bin_paths.sh`. This
will enable, for example, simply typing `make_reactivity_profiles.py` instead of
`internals/thirdparty/miniconda/bin/python3 internals/bin/make_reactivity_profiles.py`.


Reverse engineering
-------------------

To aid in understanding the various files, executables, and commandline
parameters used in a run, we recommend executing `run_example_modular.sh`, or
running shapemapper on a small dataset with the following additional parameters:
<kbd>--serial</kbd>
<kbd>--verbose</kbd>
<kbd>--render-flowchart</kbd>
<kbd>--output-processed-reads</kbd>
<kbd>--output-aligned-reads</kbd>
<kbd>--output-parsed-mutations</kbd>
<kbd>--output-counted-mutations</kbd>

Examine the log file to see the exact commands used by shapemapper to run
each module, and examine primary intermediate and output files in the `shapemapper_out`
folder. Additional intermediate files will be present within the
`shapemapper_temp` folder. A workflow graphic will be written to `shapemapper_out/*flowchart.svg`;
this can be inspected in software such as [Inkscape](https://inkscape.org/).

Main steps
----------

### 1. Alignment

ShapeMapper performs sequence alignment using 
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
by default, or optionally
[STAR](https://github.com/alexdobin/STAR). 

Running an alignment typically consists of two separate commands: building an index,
then aligning reads using the index. There are many tutorials online on how to
perform these steps.

See [Alignment to reference sequences](analysis_steps.md#alignment-to-reference-sequences)
for an overview geared toward SHAPE-MaP and see [Aligner parameters](analysis_steps.md#aligner-parameters) 
for commandline details.

#### Important notes:

- Usage of a gapped aligner is highly recommended, since sequence deletions are
  a large component of the MaP signal (at least for backbone adducts read out with
  the SuperScript II reverse transcriptase under relaxed-fidelity conditions). See 
  Fig. 1C in [Busan and Weeks, 2018](http://rnajournal.cshlp.org/content/early/2017/11/07/rna.061945.117).
- Alignments must be in [SAM format](https://github.com/samtools/hts-specs/raw/master/SAMv1.pdf) 
  for processing with ShapeMapper.
- ShapeMapper requires the `MD` field to be present in each SAM read
- Paired reads (if present) must be located on adjacent lines for correct handling. 
- Sorting reads by mapped location is not required.

### 2. Alignment parsing and mutation processing

This stage parses a single `.sam` alignment file as input and produces a single `.mut` file as output containing
read mapping information and processed mutations for each read passing quality filters. 

To run this module in isolation, use `internals/bin/shapemapper_mutation_parser`. Run with
no arguments for commandline help.

See [Parsed mutations](file_formats.md#parsed-mutations) for output file format.

Analyses performed by `shapemapper_mutation_parser` are documented in 
[Analysis steps](analysis_steps.md#analysis-steps):
- [Primer trimming and enforcement of read location requirements](analysis_steps.md#primer-trimming-and-enforcement-of-read-location-requirements)
- [Ambiguously aligned mutation handling](analysis_steps.md#ambiguously-aligned-mutation-handling)
- [Multinucleotide mutation handling](analysis_steps.md#multinucleotide-mutation-handling)
- [Post-alignment basecall quality filtering](analysis_steps.md#post-alignment-basecall-quality-filtering)
- [Chemical adduct location inference](analysis_steps.md#chemical-adduct-location-inference)

### 3. Mutation counting

This stage takes a single `.mut` parsed mutations file as input and outputs a table of
counted mutations and read depths.

Use `shapemapper_mutation_counter` to run this module. Run with no arguments
for commandline help.

See [Mutation counts](file_formats.md#mutation-counts) for output file format.

### 4. Reactivity profile calculation

This stage takes mutation counts and read depth tables from one, two, or three samples 
(modified, untreated, denatured control), computes an overall reactivity profile, normalizes
that profile, generates output figures, and performs quality control checks.

#### Compute profile
To calculate a reactivity profile from mutation counts tables, run
`make_reactivity_profiles.py`. Run with <kbd>--help</kbd> for usage.
This script produces a single reactivity profile table as output.

See [column descriptions](file_formats.md#profile-format).
Also see [Calculation of mutation rates](analysis_steps.md#calculation-of-mutation-rates) and
[Reactivity profile calculation](analysis_steps.md#reactivity-profile-calculation-and-normalization).

#### Normalize profile
To normalize a reactivity profile, run `normalize_profiles.py`.
For most situations, provide a single reactivity profile table as input with 
<kbd>--tonorm</kbd> (the other arguments can be safely ignored). 
The file will be overwritten with added columns for normalized reactivity and reactivity stderr.

See [Reactivity profile normalization](analysis_steps.md#normalization).

#### Convert file formats
To convert a reactivity profile into various output formats convenient for
downstream software, run `tab_to_shape.py`. 
See [File formats](file_formats.md#shape-format).

#### Render figures
To perform quality control checks and render summary figures, run
`render_figures.py`. 

See [Quality control checks](analysis_steps.md#quality-control-checks).

&nbsp;&nbsp;&nbsp;&nbsp;

[← back to README](../README.md)
