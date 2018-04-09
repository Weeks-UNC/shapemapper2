# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import string
import sys
import argparse

from pyshapemap.pipeline_builder import build_pipeline
from pyshapemap.util import version

# TODO: move some of these functions to util
# TODO: add strand-preserving mode. This will require maintaining unmerged reads
#       in separate files and possibly multiple aligner invokations with explicit
#       strand arguments (--norc/--nofw for Bowtie2)

def check_fasta(filename):
    """
    Check extension for fasta. Raise ValueError if
    extension not recognized.

    """
    fa_exts = [".fa", ".fasta"]
    p, ext = os.path.splitext(filename)
    if not ext.lower() in fa_exts:  # TODO: check if bowtie2, STAR handle gzipped fa files
        raise ValueError("Error: \"" + filename + "\" does not match expected extensions: " + str(fa_exts))


def check_fastq(filename):
    """
    Check extension for fastq or gzipped fastq. Raise
    ValueError if extension not recognized.

    """
    fq_exts = [".fq", ".fastq"]
    p, ext = os.path.splitext(filename)
    val = False
    if ext.lower() in [".gz"]:
        p2, ext2 = os.path.splitext(p)
        val = ext2.lower() in fq_exts
    else:
        val = ext.lower() in fq_exts
    if not val:
        raise ValueError("Error: \"" + filename + "\" does not match expected extensions: " + str(fq_exts))


def check_folder(path):
    """
    Check path contains no extension-like pattern, raise
    ValueError if it does.

    """
    p, ext = os.path.splitext(path)
    if not len(ext) == 0:
        raise ValueError("Error: \"" + path + "\" does not look like a folder (no extensions allowed)")


def check_folder_exists(path):
    """
    Check folder exists, raise ValueError if not.

    """
    if not os.path.isdir(path):
        raise ValueError("Error: \"" + path + "\" is not a folder")


def string_distance(s1, s2):
    """
    Calculate the number of characters that differ
    between two strings of identical length. Returns
    1 if lengths do not match.

    """
    if len(s1) != len(s2):
        return 1
    diff_count = 0
    for c1, c2, in zip(s1, s2):
        if c1 != c2:
            diff_count += 1
    return diff_count


def split_sample_args(args):
    """
    argparse chokes on repeated subparser options,
    so need to "manually" group args by sample name. There's
    probably a better way to do this.
    """
    samples = ["modified", "untreated", "unmodified", "denatured", "correct-seq"]
    groups = {}
    rest = []
    current_sample = ""
    for i in range(len(args)):
        # match against sample names
        matches = []
        v = args[i][2:]
        for sample in samples:
            if v == sample[:len(v)]:
                matches.append(sample)
        if len(matches) > 1:
            msg = 'Error: ambiguous sample argument "{}". Use "modified", "untreated",'
            msg += ' "unmodified", "denatured", or "correct-seq".'
            msg = msg.format(args[i])
            raise RuntimeError(msg)
        elif len(matches) == 1:
            current_sample = matches[0].replace("correct-seq","correct_seq")
            groups[current_sample] = []
        elif current_sample != "":
            groups[current_sample].append(args[i])
        else:
            rest.append(args[i])
    if "untreated" in groups and "unmodified" in groups:
        msg = 'Error: specify either "--untreated" or "--unmodified", not both '
        msg += '(these are two names for the same sample)'
        raise RuntimeError(msg)

    if "unmodified" in groups:
        groups["untreated"] = groups.pop("unmodified")

    return groups, rest


usage = """
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



    """

def parse_args(args):
    """

    Args:
        args: commandline arguments, usually sys.argv[1:]

    Returns:
        Namespace, dict: namespace containing global parameters,
            dictionary containing fastq files or lists, grouped
            by sample
    """


    # FIXME: remove/consolidate some redundant defaults

    # Multiple parsers screw up usage string, so just using my own.
    parser = argparse.ArgumentParser(prog="shapemapper",
                                     usage=usage)
    parser.add_argument('--target', type=str, nargs='+', required=False)
    # FIXME: document this argument (pass fasta sequence directly on commandline)
    parser.add_argument('--target-raw', type=str, default="")

    parser.add_argument('--out', type=str, default="shapemapper_out")
    parser.add_argument('--temp', type=str, default="shapemapper_temp")
    parser.add_argument('--log', type=str, default="shapemapper_log.txt")
    parser.add_argument('--structured-output', action="store_true", default=False)
    parser.add_argument('--render-flowchart', action="store_true", default=False)

    parser.add_argument('--nproc', type=int, default=4)
    parser.add_argument('--serial', action="store_true", default=False)
    parser.add_argument('--verbose', action="store_true", default=False)
    parser.add_argument('--overwrite', action="store_true", default=False)

    parser.add_argument('--name', type=str, default='')
    parser.add_argument('--random-primer-len', type=int, default=0)
    parser.add_argument('--min-depth', type=int, default=5000)
    parser.add_argument('--max-bg', type=float, default=0.05)
    parser.add_argument('--min-mapq', type=int, default=10)
    parser.add_argument('--preserve-order', action="store_true", default=False)

    parser.add_argument('--star-aligner', action="store_true", default=False)

    parser.add_argument('--disable-soft-clipping', action="store_true", default=False)

    parser.add_argument('--right-align-ambig',action="store_true", default=False)
    parser.add_argument('--right-align-ambig-dels',action="store_true", default=False)
    parser.add_argument('--right-align-ambig-ins',action="store_true", default=False)

    parser.add_argument('--min-mutation-separation', type=int, default=6)

    parser.add_argument('--mutation-type-to-count', type=str, default="")

    parser.add_argument('--min-qual-to-trim', type=int, default=20)
    parser.add_argument('--window-to-trim', type=int, default=5)
    parser.add_argument('--min-length-to-trim', type=int, default=25)

    parser.add_argument('--min-qual-to-count', type=int, default=30)

    parser.add_argument('--indiv-norm', action="store_true", default=False)

    parser.add_argument('--min-seq-depth', type=int, default=50)
    parser.add_argument('--min-freq', type=float, default=0.6)

    parser.add_argument('--output-processed', '--output-processed-reads',
                        dest="output_processed_reads", action="store_true", default=False)
    parser.add_argument('--output-aligned', action="store_true", default=False)
    parser.add_argument('--output-parsed', '--output-parsed-mutations',
                        dest="output_parsed", action="store_true", default=False)
    parser.add_argument('--output-counted', '--output-counted-mutations',
                        dest="output_counted", action="store_true", default=False)
    parser.add_argument('--separate-ambig-counts', action="store_true", default=False)
    parser.add_argument('--output-classified', '--output-classified-mutations',
                        dest="output_classified", action="store_true", default=False)

    parser.add_argument('--render-mutations', action="store_true", default=False)

    parser.add_argument('--calc-correlations', action="store_true", default=False)
    parser.add_argument('--corr-args', type=str, default='')

    # TODO: optional syntax to auto-generate sample names from file/folder names

    fileparser = argparse.ArgumentParser()
    # TODO: support single-dash versions of --R1, --R2, --U?
    fileparser.add_argument("--R1", type=str, nargs='+')
    fileparser.add_argument("--R2", type=str, nargs='+')
    fileparser.add_argument("--U", type=str, nargs='+')
    fileparser.add_argument("--folder", type=str, nargs='+')
    fileparser.add_argument("--unpaired-folder", type=str, nargs='+')

    # give a specific message if no args provided
    if len(args) == 0:
        print("No arguments provided.")
        print("Usage: "+usage.lstrip())
        sys.exit(1)

    # first parse "global" options
    p, rest = parser.parse_known_args(args)

    # set some dependent flags for debug output option
    if p.render_mutations:
        p.output_aligned = True
        p.output_classified = True
        p.output_parsed = True
        p.preserve_order = True

    if p.calc_correlations:
        p.output_classified = True
        p.corr_args = p.corr_args.split(' ')

    if p.right_align_ambig:
        p.right_align_ambig_dels = True
        p.right_align_ambig_ins = True

    if p.target_raw == "" and (p.target is None or len(p.target) == 0):
        raise RuntimeError("Error: must provide at least one target sequence (--target)")

    if p.target is None:
        p.target = []

    # check sensible file extensions
    if p.target:
        for fa in p.target:
            check_fasta(fa)
    if p.out:
        check_folder(p.out)

    # then parse samples
    groups, rest = split_sample_args(rest)
    if len(rest) > 0:
        raise RuntimeError("Error: unrecognized argument(s): {}".format(rest))

    fastqs = {}
    def store_args(s, s_args):
        if sum([x is not None for x in [s_args.R1 and s_args.R2,
                                        s_args.U,
                                        s_args.folder,
                                        s_args.unpaired_folder]]) > 1:
            raise RuntimeError("Error: too many input arguments specified for sample {}.".format(s))
        if ( (bool(s_args.R1 is None) != bool(s_args.R2 is None)) or 
             ( ((s_args.R1 is not None) and (s_args.R2 is not None)) 
               and len(s_args.R1) != len(s_args.R2) ) ):
            raise RuntimeError("Error: specified R1 or R2 without matching paired file.")
        if sample_args.R1 and sample_args.R2:
            for f in sample_args.R1:
                check_fastq(f)
            for f in sample_args.R2:
                check_fastq(f)
            fastqs[s] = {"R1":sample_args.R1,
                         "R2":sample_args.R2}
        elif sample_args.U:
            for f in sample_args.U:
                check_fastq(f)
            fastqs[s] = {"U":sample_args.U}
        elif sample_args.folder:
            fastqs[s] = {"R1":[], "R2":[]}
            for f in sample_args.folder:
                R1, R2 = parse_paired_input_folder(f)
                fastqs[s]["R1"] += R1
                fastqs[s]["R2"] += R2
        elif sample_args.unpaired_folder:
            fastqs[s] = {"U":[]}
            for f in sample_args.unpaired_folder:
                U = parse_unpaired_input_folder(f)
                fastqs[s]["U"] += U
        else:
            msg = "Error: must provide either paired fastq files (--R1, --R2)"
            msg += ", folder of paired fastqs (--folder), "
            msg += "or folder of unpaired fastqs (--unpaired-folder) for each sample "
            msg += "provided."
            raise RuntimeError(msg)

    for sample in groups:
        sample_args, rest = fileparser.parse_known_args(groups[sample])
        if len(rest) > 0:
            raise RuntimeError("Error: unrecognized argument(s): {}".format(rest))
        store_args(sample, sample_args)

    if "modified" not in fastqs and "correct_seq" not in fastqs:
        msg = "Error: must provide FASTQ files or folders from "
        msg += "at least one modified sample (--modified --R1 <f1.fastq> --R2 <f2.fastq>)"
        msg += " or a sample to identify sequence variants (--correct-seq)."
        raise RuntimeError(msg)

    if "denatured" in fastqs and "untreated" not in fastqs:
        msg = "Error: denatured control specified without untreated control."
        raise RuntimeError(msg)

    possible_mutation_types = ["", "mismatch", "gap", "insert", 
                               "gap_multi", "insert_multi", 
                               "complex"]
    if p.mutation_type_to_count not in possible_mutation_types:
        msg = 'Unrecognized argument "{}" to "--mutation-type-to-count". Possible values: {}'
        msg = msg.format(p.mutation_type_to_count, 
                         ', '.join(['"{}"'.format(x) for x in possible_mutation_types]))
        raise RuntimeError(msg)

    return p, fastqs

def parse_paired_input_folder(input_folder):
    check_folder_exists(input_folder)
    R1 = []
    R2 = []
    file_list = [f for f in os.listdir(input_folder) if not os.path.isdir(f)]
    exts = [".fastq", ".fq", ".fastq.gz"]
    file_list = [f for f in file_list if any([f.endswith(ext) for ext in exts])]
    for f in file_list:
        # try to locate "R1" or "R2" in filename, separated from other fields
        # by underscores or periods
        fields = os.path.splitext(os.path.split(f)[1])[0].replace('.','_').split('_')
        if "R1" in fields or "r1" in fields:
            R1.append(f)
        elif "R2" in fields or "r2" in fields:
            R2.append(f)
        else:
            msg = "Error: FASTQ file(s) present in folder \"{}\" that do not contain".format(input_folder)
            msg += " underscore-separated R1 and R2 fields in the filename. Please rename these files to "
            msg += "identify paired reads, or provide unpaired reads with --unpaired-folder."
            raise RuntimeError(msg)
    R1.sort()
    R2.sort()
    if len(R1)==0 and len(R2)==0:
        msg = "Error: no fastq reads found in folder \"{}\"".format(input_folder)
        raise RuntimeError(msg)
    for f1, f2 in zip(R1, R2):
        if string_distance(f1, f2) > 1:
            msg = "Error: unable to identify paired read FASTQ files in folder \"" + input_folder + "\""
            msg += ". Ensure that paired files contain underscore-separated R1 and R2 fields in the "
            msg += "filename, and that their filenames are otherwise identical."
            raise RuntimeError(msg)
    R1 = [os.path.join(input_folder, f) for f in R1]
    R2 = [os.path.join(input_folder, f) for f in R2]

    return R1, R2

def parse_unpaired_input_folder(input_folder):
    check_folder_exists(input_folder)
    file_list = [f for f in os.listdir(input_folder) if not os.path.isdir(f)]
    exts = [".fastq", ".fq", ".fastq.gz"]
    file_list = [f for f in file_list if any([f.endswith(ext) for ext in exts])]
    file_list.sort()

    if len(file_list)==0:
        msg = "Error: no fastq reads found in folder {}".format(input_folder)
        raise RuntimeError(msg)
    U = [os.path.join(input_folder, f) for f in file_list]

    return U


def check_help_version(args):
    """
    Check if --help or --version is present in args. Display
    appropriate meesage and exit if so.
    """
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--version", action="version", 
                        version="ShapeMapper v{}".format(version()))
    p, rest = parser.parse_known_args(args)
    return


def check_no_whitespace(args):
    """
    Don't allow whitespace characters in any input arguments (currently causes problems for
    bowtie2 wrapper, which doesn't quote it's -x parameter).
    FIXME: remove once bowtie2 PR accepted and new release propagates to anaconda
    """
    for arg in args:
        for char in arg:
            if char in string.whitespace:
                raise RuntimeError("No whitespace characters are currently allowed in input arguments. Replace spaces in file and folder names with underscores ('_').")
    return

def get_log_path(args):
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--log", type=str)
    parser.add_argument("--name", type=str)

    p, rest = parser.parse_known_args(args)
    log = p.log
    if p.log is None:
        log = "shapemapper_log.txt"
        if p.name is not None:
            folder, filename = os.path.split(log)
            log = os.path.join(folder, p.name+'_'+filename)
            # put name back in arg list, since parsed later
            rest += ['--name', p.name]
    return log, rest


def get_out_path(args):
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--out", type=str, default="shapemapper_out")

    p, rest = parser.parse_known_args(args)
    return p.out, rest


def get_temp_path(args):
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("--temp", type=str, default="shapemapper_temp")

    p, rest = parser.parse_known_args(args)
    return p.temp, rest

def get_paths(args):
    """
    Get logfile, output folder, and temporary folder paths
    from commandline args. Also returns args with log argument
    removed.

    """
    log, rest = get_log_path(args)
    out, _ = get_out_path(args)
    temp, _ = get_temp_path(args)
    return log, out, temp, rest


def construct(args,
              **kwargs):
    """
    Construct a pipeline from commandline arguments.

    Args:
        args: Commandline arguments, typically sys.argv[1:]

    Returns:
        SerialPipeline or ParallelPipeline
    """
    params, fastqs = parse_args(args)

    kw = vars(params)
    kw.update(kwargs)

    return build_pipeline(fastqs=fastqs,
                          **kw)


