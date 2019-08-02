"""
Specific wrappers for pipeline executables, including
the Bowtie2 sequence aligner, BBmerge, and compiled shapemapper
modules for read trimming, mutation parsing, and mutation
counting.
"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import subprocess
from math import log2

from pyshapemap.component import *
from pyshapemap.util import require_explicit_kwargs, sanitize
from pyshapemap.nodes import connect as connect_nodes
from pyshapemap.nodes import disconnect as disconnect_nodes
from pyshapemap.nodes import connect_shared_input_nodes
from pyshapemap.connect import connect

this_dir = os.path.dirname(os.path.realpath(__file__))
bin_dir = os.path.join(this_dir, "../../bin")
pyexe = "python3"

# hacks for debugging
DISABLE_MERGING = False


# TODO: expose more module-specific params: could just forward custom args, kwargs to each process

class FastaFormatChecker(Component):
    def __init__(self,
                 fasta=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="fasta",
                           parallel=False))
        self.add(OutputNode(name="corrected",
                            extension="passthrough",
                            parallel=False))
        if fasta is not None:
            if isinstance(fasta, str):
                self.fasta.set_file(fasta)
            elif isinstance(fasta, Node):
                connect(fasta, self.fasta)
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "check_fasta_format.py"),
               "{fasta}",
               "{corrected}"]
        return cmd


class QualityTrimmer(Component):
    def __init__(self,
                 fastq=None,
                 min_qual=None,
                 window=None,
                 min_length=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.min_qual = min_qual
        self.min_length = min_length
        self.window = window
        self.add(InputNode(name="fastq"))
        if fastq is not None:
            self.fastq.set_file(fastq)
        self.add(OutputNode(name="trimmed",
                            extension="fastq"))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = "shapemapper_read_trimmer -i {fastq} -o {trimmed}"
        if self.min_qual is not None:
            cmd += " -p {min_qual}"
        if self.min_length is not None:
            cmd += " -l {min_length}"
        if self.window is not None:
            cmd += " -w {window}"
        return cmd

# NOTE: this is only used so bbmerge doesn't crash with pipe inputs
class Interleaver(Component):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="R1"))
        self.add(InputNode(name="R2"))
        self.add(OutputNode(name="interleaved",
                            extension="fastq"))
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "interleave_fastq.py"),
               "{R1}", "{R2}", "{interleaved}"]
        return cmd


class Tab6Interleaver(Component):
    def __init__(self,
                 separate_files=False,
                 **kwargs):
        self.separate_files = separate_files
        super().__init__(**kwargs)

        if self.separate_files:
            self.add(InputNode(name="R1", extension="fastq"))
            self.add(InputNode(name="R2", extension="fastq"))
        else:
            self.add(InputNode(name="fastq")) # mixed paired/unpaired
        self.add(OutputNode(name="tab6",
                            extension="tab6"))
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "tab6_interleave.py")]
        if self.separate_files:
            cmd += ["--R1", "{R1}",
                    "--R2", "{R2}"]
        else:
            cmd += ["--input", "{fastq}"]

        cmd += ["--output", "{tab6}"]

        return cmd


class Deinterleaver(Component):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="interleaved"))
        self.add(OutputNode(name="R1", extension="fastq"))
        self.add(OutputNode(name="R2", extension="fastq"))
        self.add(OutputNode(name="unpaired", extension="fastq"))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        #cmd = "paste - - - - - - - - < {interleaved} "
        #cmd += "| tee >(cut -f 1-4 | tr '\\t' '\\n' > {R1}) "
        #cmd += "| cut -f 5-8 | tr '\\t' '\\n' > {R2}"
        cmd = [pyexe,
               os.path.join(bin_dir, "deinterleave_fastq.py"),
               "--input", "{interleaved}",
               "--R1-out", "{R1}",
               "--R2-out", "{R2}",
               "--unpaired-out", "{unpaired}"]
        return cmd


class Appender(Component):
    """
    combine multiple (usually fastq) files into one by appending
    each file end-to-end

    """

    def __init__(self,
                 inputs=None,
                 add_extra_newline=False,
                 **kwargs):
        super().__init__(**kwargs)
        self.add(OutputNode(name="appended",
                            extension="passthrough"))
        self.add(StderrNode())
        self.add_extra_newline = add_extra_newline
        if add_extra_newline:
            # TODO: maybe hide temp output nodes or color differently in flowchart
            self.add(OutputNode(name="linebreak",
                                extension="txt",
                                parallel=False))

        # Option to set input files using args to constructor
        if inputs is not None:
            for i in range(len(inputs)):
                f = inputs[i]
                name = "f{}".format(i + 1)
                if isinstance(f, str):
                    self.add(InputNode(name=name,
                                       filename=f))
                elif isinstance(f, Node):
                    node = InputNode(name=name)
                    self.add(node)
                    connect(f, node)

    def cmd(self):
        if self.add_extra_newline:
            # for concatenating some files (like FASTA), add an extra linebreak
            # between files to ensure headers appear on their own lines
            cmd = "echo -e '\n' > {linebreak}; cat".split(' ')
            for i in range(len(self.input_nodes)):
                cmd += ["{{{}}}".format(self.input_nodes[i].get_name())]
                if i < len(self.input_nodes) - 1:
                    cmd += ["{linebreak}"]
        else:
            cmd = ["cat"]
            for node in self.input_nodes:
                cmd += ["{{{}}}".format(node.get_name())]
        cmd += [">", "{appended}"]
        return cmd


class Merger(Component):
    # TODO: expose more parameters
    # Note: "merged" output stream now includes both merged and
    #       unmerged reads
    def __init__(self,
                 preserve_order=None,
                 nproc=4,
                 **kwargs):
        self.nproc = nproc
        self.preserve_order = preserve_order
        super().__init__(**kwargs)
        self.add(InputNode(name="interleaved_fastq"))
        self.add(StdoutNode(name="output",
                            extension="fastq",
                            parallel=True))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = ["bbmerge.sh",
               "vstrict=t",
               "in=stdin",
               "out=stdout",
               "outu=stdout",
               #"out={merged}",
               #"outu={unmerged}",
               "interleaved=t",
               "usejni=t", # FIXME: autodetect whether JNI components are compiled
               "t={}".format(self.nproc), # number of threads
               #"-eoom",
               ]
        if self.preserve_order:
            cmd += ["ordered=t"]
        if DISABLE_MERGING:
            cmd += ["minoverlap=2000"] # quick hack to disable merging
        cmd += [">", "{output}"]
        cmd += ["<", "{interleaved_fastq}"]
        return cmd

    # Process doesn't seem to always exit on error, so make a special error wrapper that
    # just checks whether certain bad words appear in stderr
    def proc_status(self):
        status = Component.proc_status(self)
        if status == "failed":
            return status
        r = self.read_stderr()
        try:
            r = "\n".join(r.splitlines()[5:])  # first few lines are typically parameters and filename details
        except IndexError:
            r = ""
        # Ignore broken pipe errors, since those can be caused by a propagating downstream failure
        if ("Exception" in r or "Error" in r) and "Broken pipe" not in r:
            # print("\n\n\n######Error found, returning self#########\n\n\n\n")
            return "failed"
        else:
            return status

    def after_run_message(self):
        # stats are written to stderr
        return self.read_stderr()


class SamMixer(Component):
    '''
    Mix two SAM alignment streams, skipping header lines.
    '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="sam1", extension='sam'))
        self.add(InputNode(name="sam2", extension='sam'))
        self.add(OutputNode(name="mixed", extension='sam'))
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "mix_sam.py"),
               "{sam1}", "{sam2}", "{mixed}"]
        return cmd


class BowtieIndexBuilder(Component):
    out_extension = ""

    def __init__(self,
                 target=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="target",
                           filename=target,
                           parallel=False))
        self.add(OutputNode(name="index",
                            parallel=False,
                            extension=''))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = "bowtie2-build {target} {index}"
        return cmd


class BowtieAligner(Component):
    def __init__(self,
                 reorder=False,
                 disable_soft_clipping=False,
                 nproc=None,
                 maxins=None,
                 max_search_depth=None,
                 max_reseed=None,
                 **kwargs):
        self.reorder = reorder
        self.nproc = nproc
        self.disable_soft_clipping = disable_soft_clipping
        self.maxins = 800

        if maxins is not None:
            self.maxins = maxins
        self.max_search_depth = max_search_depth
        self.max_reseed = max_reseed
        if self.max_search_depth < 0:
            self.max_search_depth = None
        if self.max_reseed < 0:
            self.max_reseed = None
        if max_reseed is None and max_search_depth is not None:
            max_reseed = 2 # bowtie2 -R
        if max_search_depth is None and max_reseed is not None:
            max_search_depth = 15 # bowtie2 -D
        super().__init__(**kwargs)
        self.add(InputNode(name="index",
                           parallel=False))
        self.add(StdinNode(name="tab6",
                           extension="tab6",
                           parallel=True))
        self.add(OutputNode(name="aligned",
                            extension="sam",
                            parallel=True))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = ["bowtie2_wrapper.sh"]
        #cmd = ["bowtie2"]
        cmd += ["-p", str(self.nproc)]

        if not self.disable_soft_clipping:
            cmd += ["--local"]
            if self.max_search_depth is None:
                cmd += ["--sensitive-local"]
                # corresponds to -D 15 -R 2 -N 0 -L 20 -i S,1,0.75
            else:
                cmd += ["-D", str(self.max_search_depth), "-R", str(self.max_reseed), "-N", "0", "-L", "20", "-i", "S,1,0.75"]
        else:
            cmd += ["--end-to-end"]
            if self.max_search_depth is None:
                cmd += ["--sensitive"]
                # corresponds to -D 15 -R 2 -N 0 -L 22 -i S,1,1.15
            else:
                cmd += ["-D", str(self.max_search_depth), "-R", str(self.max_reseed), "-N", "0", "-L", "22", "-i", "S,1,0.75"]
        cmd += [
               "--mp", "3,1", # try to be slightly more consistent with STAR's alignment parameters
               "--rdg", "5,1", # multinuc deletions are a large part of the signal, so don't penalize
                               # as much as bowtie's default setting (open penalty 5, extension penalty 3)
               "--rfg", "5,1",
               "--dpad", "30", # this seems reasonable to allow fairly large deletions, but not too crazy (although these do exist)
               "--maxins", str(self.maxins),
               "--ignore-quals",
               "--no-unal",  # don't produce SAM records for reads that didn't map

               #"--quiet", # suppress warnings about e.g. R1 or R2 in a pair being zero-length
                           # - unfortunately this seems to also eliminate all alignment stats output
                           # - bin/bowtie2_wrapper.sh hacks around this to filter out spurious warnings
                           # - could pass reads missing mate pair to bowtie as unpaired input, but then
                           #   downstream stages can't distinguish merged reads from reads missing a mate pair
               #"--dovetail", # allow dovetailed paired reads to be considered concordantly aligned
               ]

        if self.reorder:
            cmd += ["--reorder"] # output in same order as input for debugging purposes

        cmd += ["--tab6", "-"]
        #cmd += ["-U", "{fastq}"]

        cmd += ["-x", "{index}"]
        cmd += ["<", "{tab6}"]
        cmd += ["-S", "{aligned}"]
        #cmd += [">", "{aligned}"]
        cmd = ' '.join(cmd)
        return cmd

    def after_run_message(self):
        # Note: Bowtie2 alignment stats are written to its stderr,
        # but the wrapper hack redirects stderr into stdout to enable suppressing some warnings
        return self.read_stdout()
        #return self.read_stderr()


class StarIndexBuilder(Component):
    out_extension = ""

    def __init__(self,
                 target=None,
                 num_targets=None,
                 total_target_length=None,
                 nproc=None,
                 genomeSAindexNbase=None,
                 **kwargs):
        self.num_targets = num_targets
        self.total_target_length = total_target_length
        self.nproc = nproc
        super().__init__(**kwargs)
        self.add(InputNode(name="target",
                           filename=target,
                           parallel=False))
        # - index folder must exist before running STAR,
        #   but temp folder must not exist
        self.add(OutputNode(name="index",
                            isfolder=True))
        self.add(OutputNode(name="temp",
                            isfolder=True,
                            make_parent=True,
                            error_on_existing=True))
        self.add(OutputNode(name="logs",
                            isfolder=True))
        self.add(StdoutNode())
        self.add(StderrNode())
        # The STAR manual recommends setting genomeSAindexNbases to
        # min(14, log2(GenomeLength)/2 - 1) to handle "short" reference
        # sequences (e.g. bacterial genomes or smaller), so the called script calculates
        # this value automatically.
        # STAR seems to segfault if genomeSAindexNbases is too high, even by a small amount
        # The manual also recommends setting genomeChrBinNbits to
        # min(18, log2(GenomeLength/NumberOfReferences))
        if genomeSAindexNbase == 0:
            self.genomeSAindexNbase = int(round(min(14, log2(self.total_target_length) / 2.0 - 1)))
        else:
            # allow setting directly (small target sequences may still segfault and require a
            # lower value (3 or 2))
            self.genomeSAindexNbase = genomeSAindexNbase
        self.genomeChrBinNbits = int(min(18, log2(self.total_target_length)/self.num_targets))

    def cmd(self):
        cmd = ["STAR",
               "--runMode", "genomeGenerate",
               "--genomeDir", "{index}",
               "--outTmpDir", "{temp}",
               "--outFileNamePrefix", "{logs}/",
               "--genomeFastaFiles", "{target}",
               "--runThreadN", str(self.nproc),
               "--genomeSAindexNbases", str(self.genomeSAindexNbase),
               "--genomeChrBinNbits", str(self.genomeChrBinNbits)]
        return cmd

# WARNING: unlike Bowtie2, STAR will align to all indices in a directory
# WARNING: STAR's performance degrades significantly if target sequences
#          are incomplete (for example, only providing the small subunit 
#          ribosomal sequence when both subunits are present in the data)
# FIXME: STAR is reporting a sizable fraction of multimappers even for a simple
#        amplicon dataset with a single target. Might be a result of ambiguous
#        indel alignment - need to figure out if that influences STAR's MAPQ values.
class StarAligner(Component):
    def __init__(self,
                 reorder=False,
                 disable_soft_clipping=False,
                 nproc=None,
                 paired=False,
                 shared_index=False,
                 fixed_index=None, # for debugging
                 **kwargs):
        #require_explicit_kwargs(locals())
        self.reorder = reorder
        self.nproc = nproc
        self.disable_soft_clipping = disable_soft_clipping
        self.shared_index = shared_index # use shared-memory index (seems buggy)
        self.fixed_index = fixed_index # use index located at given path
        # Note: STAR does not natively support mixed paired/unpaired input or
        #       interleaved fastq
        self.paired = paired
        super().__init__(**kwargs)
        if fixed_index is None:
            self.add(InputNode(name="index",
                               isfolder=True))
        if paired:
            self.add(InputNode(name="R1"))
            self.add(InputNode(name="R2"))
        else:
            self.add(InputNode(name="fastq", parallel=True))
        self.add(OutputNode(name="temp",
                            isfolder=True,
                            make_parent=True,
                            error_on_existing=True))
        self.add(OutputNode(name="logs",
                            isfolder=True))
        self.add(OutputNode(name="aligned",
                            extension="sam",
                            parallel=True))
        self.add(StderrNode())

    def cmd(self):
        cmd = ["STAR"]
        if self.paired:
            cmd += ["--readFilesIn", "{R1}", "{R2}"]
        else:
            cmd += ["--readFilesIn", "{fastq}"]
        cmd += ["--runThreadN", str(self.nproc)]
        if not self.disable_soft_clipping:
            cmd += ["--alignEndsType", "Local"]
        else:
            cmd += ["--alignEndsType", "EndToEnd"]
        cmd += [
               "--scoreGap", "-1000000",  # disable splice junction detection
               # reduce gap extension penalties (multinuc gaps are
               # a large part of the signal)
               "--scoreDelBase", "-1",
               "--scoreInsBase", "-1",
               # disable filtering by mismatch counts
               "--outFilterMismatchNmax", "999",
               "--outFilterMismatchNoverLmax", "999",
               ]
        if self.reorder:
            cmd += ["--outSAMorder", "PairedKeepInputOrder"]
        if self.shared_index:
            cmd += ["--genomeLoad", "LoadAndRemove"]
        if self.fixed_index is not None:
            cmd += ["--genomeDir", self.fixed_index]
        else:
            cmd += ["--genomeDir", "{index}"]
        cmd += ["--runMode", "alignReads",

                # only report one of the top alignments for a read that maps to multiple locations
                "--outMultimapperOrder", "Random",
                "--outSAMmultNmax", "1",
                "--outStd", "SAM",
                "--outSAMattributes", "MD",
                "--outTmpDir", "{temp}",
                "--outFileNamePrefix", "{logs}/",

                "> ", "{aligned}"]
        return cmd

    def after_run_message(self):
        return open(os.path.join(self.logs.output_nodes[0].foldername, "Log.final.out"), "rU").read()


class StarAlignerMixedInput(Component):
    '''
    Wrapper for STAR aligner to support mixed paired/unpaired input
    (spawns two STAR instances and mixes output streams)
    '''
    def __init__(self,
                 reorder=False,
                 disable_soft_clipping=False,
                 nproc=None,
                 star_shared_index=False,
                 **kwargs):
        super().__init__(**kwargs)
        deinterleaver = Deinterleaver()
        aligner_kwargs = dict(reorder=reorder,
                              disable_soft_clipping=disable_soft_clipping,
                              nproc=nproc)
        aligner_paired = StarAligner(name="StarAligner_paired",
                                     paired=True,
                                     shared_index=star_shared_index,
                                     **aligner_kwargs)
        aligner_unpaired = StarAligner(name="StarAligner_unpaired",
                                       paired=False,
                                       shared_index=star_shared_index,
                                       **aligner_kwargs)
        sam_mixer = SamMixer()
        connect(deinterleaver.R1, aligner_paired.R1)
        connect(deinterleaver.R2, aligner_paired.R2)

        connect(deinterleaver.unpaired, aligner_unpaired.fastq)
        connect(aligner_paired.aligned, sam_mixer.sam1)
        connect(aligner_unpaired.aligned, sam_mixer.sam2)

        self.add(SharedInputNode(name="index"))
        connect(self.index, [aligner_paired.index, aligner_unpaired.index])

        self.add([deinterleaver,
                  aligner_paired,
                  aligner_unpaired,
                  sam_mixer])
        # allow access to some nodes from this component's namespace
        self.add_node(aligner_paired.index)
        self.add(deinterleaver.interleaved, alias="interleaved_fastq")
        self.add(sam_mixer.mixed, alias="aligned")


class BowtieAlignerMixedInput(Component):
    '''
    Wrapper for Bowtie2 aligner to support interleaved fastq input with
    mixed paired/unpaired.
    '''
    def __init__(self,
                 reorder=False,
                 disable_soft_clipping=False,
                 nproc=None,
                 maxins=None,
                 max_search_depth=None,
                 max_reseed=None,
                 **kwargs):
        super().__init__(**kwargs)
        tab6interleaver = Tab6Interleaver()
        aligner_kwargs = dict(reorder=reorder,
                              disable_soft_clipping=disable_soft_clipping,
                              nproc=nproc,
                              maxins=maxins,
                              max_reseed=max_reseed,
                              max_search_depth=max_search_depth)
        aligner = BowtieAligner(name="BowtieAligner",
                                **aligner_kwargs)
        connect(tab6interleaver.tab6, aligner.tab6)
        self.add([tab6interleaver,
                  aligner])
        self.add(tab6interleaver.fastq, alias="interleaved_fastq")
        self.add(aligner.index, alias="index")
        self.add(aligner.aligned, alias="aligned")




class MutationParser(Component):
    def __init__(self,
                 min_mapq=None,
                 right_align_ambig_dels=None,
                 right_align_ambig_ins=None,
                 random_primer_len=None,
                 min_mutation_separation=None,
                 min_qual=None,
                 mutation_type=None,
                 variant_mode=None,
                 maxins=None,
                 amplicon=None,
                 input_is_unpaired=None,
                max_primer_offset=None,
                require_forward_primer_mapped=None,
                require_reverse_primer_mapped=None,
                trim_primers=None,
                debug_out=None,
                 **kwargs):
        self.min_mapq = 35 # FIXME: clarify or remove this, since this gets used for
                           # sequence variant correction instead of the lower default
                           # value of 30 for mutation counting
        super().__init__(**kwargs)
        self.min_mapq = min_mapq
        self.right_align_ambig_dels = right_align_ambig_dels
        self.right_align_ambig_ins = right_align_ambig_ins
        self.random_primer_len = random_primer_len
        self.min_mutation_separation = min_mutation_separation
        self.min_qual = min_qual
        self.mutation_type = mutation_type
        self.maxins = 800
        self.input_is_unpaired = input_is_unpaired
        if maxins is not None:
            self.maxins = maxins # analogous to bowtie2's --maxins param
        self.variant_mode = variant_mode

        self.amplicon = amplicon
        self.max_primer_offset = max_primer_offset
        self.require_forward_primer_mapped = require_forward_primer_mapped
        self.require_reverse_primer_mapped = require_reverse_primer_mapped
        self.trim_primers = trim_primers

        self.write_debug_out = debug_out

        self.add(InputNode(name="input"))
        if self.amplicon:
            self.add(InputNode(name="primers"))
        self.add(OutputNode(name="parsed_mutations",
                            extension="mut",
                            parallel=True))
        if self.write_debug_out:
            self.add(OutputNode(name="debug_out", parallel=True))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = ["shapemapper_mutation_parser",
               "-i", "{input}",
               "-o", "{parsed_mutations}",
               "-w"]
        if self.min_mapq is not None:
            cmd += ["-m", "{}".format(self.min_mapq)]
        if self.right_align_ambig_dels is not None and self.right_align_ambig_dels:
            cmd += ["--right_align_ambig_dels"]
        if self.right_align_ambig_ins is not None and self.right_align_ambig_ins:
            cmd += ["--right_align_ambig_ins"]
        if self.random_primer_len is not None:
            cmd += ["--exclude_3prime", str(self.random_primer_len+1)]
        if self.min_mutation_separation is not None:
            cmd += ["--max_internal_match", str(self.min_mutation_separation-1)]
        if self.min_qual is not None:
            cmd += ["--min_qual", str(self.min_qual)]
        if self.mutation_type != '':
            cmd += ["--use_only_mutation_type", "{mutation_type}"]
        if self.maxins is not None and self.maxins:
            cmd += ["--max_paired_fragment_length", str(self.maxins)]
        if self.input_is_unpaired is not None and self.input_is_unpaired:
            cmd += ["--input_is_unpaired"]
        if self.variant_mode is not None and self.variant_mode:
            cmd += ["--variant_mode"]
        if self.amplicon:
            cmd += ["--primers", "{primers}"]
        cmd += ["--max_primer_offset", str(self.max_primer_offset)]
        if self.require_forward_primer_mapped:
            cmd += ["--require_forward_primer_mapped"]
        if self.require_reverse_primer_mapped:
            cmd += ["--require_reverse_primer_mapped"]
        if self.trim_primers:
            cmd += ["--trim_primers"]
        if self.write_debug_out:
            cmd += ["--debug_out", "{debug_out}"]
        #cmd += ["--debug"] # FIXME: remove
        return cmd


class MutationCounter(Component):
    def __init__(self,
                 target_length=None,
                 primer_pairs=None,
                 variant_out=None,
                 mutations_out=None,
                 per_read_histograms=False,
                 separate_ambig_counts=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.per_read_histograms = per_read_histograms
        self.separate_ambig_counts = separate_ambig_counts
        self.add(InputNode(name="mut"))

        # target_length is either an int parameter, or an OutputNode outputting 
        # a file containing the target_length parameter (this is to allow
        # this component to get updated sequence length after sequence correction)
        if isinstance(target_length, int):
            self.target_length = target_length
        elif isinstance(target_length, OutputNode):
            self.add(ParameterNode(name="target_length"))
            connect_nodes(target_length, self.target_length)

        # number of amplicon primer pairs
        if primer_pairs is None or isinstance(primer_pairs, int):
            self.primer_pairs = primer_pairs
        if isinstance(primer_pairs, OutputNode):
            self.add(ParameterNode(name="primer_pairs"))
            connect_nodes(primer_pairs, self.primer_pairs)

        # a little convoluted, since this component can produce
        # variable numbers of outputs, and want to also either accept
        # explicit filepath or bool indicating filepath should be generated
        # automatically
        if variant_out is not None and variant_out:
            kw = {"name": "variants"}
            if isinstance(variant_out, str):
                kw["filename"] = variant_out
            self.add(OutputNode(parallel=False,
                                **kw))
        if mutations_out is not None and mutations_out:
            kw = {"name": "mutations"}
            if isinstance(mutations_out, str):
                kw["filename"] = mutations_out
            self.add(OutputNode(parallel=False,
                                **kw))

        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        # FIXME: support arbitrary number of input files
        cmd = ["shapemapper_mutation_counter", "-i"]
        cmd += ["{mut}"]
        cmd += ["-w"]
        node_names = [n.get_name() for n in self.output_nodes]
        if "variants" in node_names:
            cmd += ["-v", "{variants}"]
        if "mutations" in node_names:
            cmd += ["-c", "{mutations}"]

        if self.target_length is not None:
            cmd += ["--length", "{target_length}"]
        if self.primer_pairs is not None:
            cmd += ["--n_primer_pairs", "{primer_pairs}"]
        if self.per_read_histograms:
            cmd += ["--hist"]
        if self.separate_ambig_counts is not None and self.separate_ambig_counts:
            cmd += ["--separate_ambig_counts"]
        return cmd

    def after_run_message(self):
        lines = self.read_stdout().splitlines()
        # just display the histogram tables if present
        i = 0
        start_i = 0
        found_table = False
        for i in range(len(lines)):
            if lines[i] == "Read lengths":
                start_i = i
                found_table = True
                break
        end_i = 0
        if found_table:
            for i in range(len(lines)-1,-1,-1):
                if lines[i] == "--------------------":
                    end_i = i
                    break
            return "\n".join(lines[start_i:end_i+1])
        else:
            return ""


class ProgressMonitor(Component):
    def __init__(self,
                 input=None,
                 **kwargs):
        self.expected_bytes = None
        super().__init__(**kwargs)
        self.add(InputNode())
        if input is not None:
            self.input.set_file(input)
        self.add(OutputNode(extension="passthrough"))
        self.add(StderrNode())

    def cmd(self):
        cmd = "pv -f -p -e -b"
        term_width, term_height = shutil.get_terminal_size()
        cmd += " -w " + str(term_width - 6)  # leave some room for indents
        if self.expected_bytes is not None:
            cmd += " -s " + str(self.expected_bytes)
        cmd += " < {input} > {output}"
        return cmd

    # need direct access to stderr pipe for things to work as expected
    # FIXME: see if I can get things working using stderr output file instead
    def start_process(self,
                      verbose=False):
        kwargs = {"shell": True,
                  "executable": "/bin/bash",
                  "stderr": subprocess.PIPE,
                  "preexec_fn": os.setsid}
        self.proc = sp.Popen(self.format_command(self.cmd()),
                             **kwargs)


class CalcProfile(Component):
    def __init__(self,
                 mindepth=None,
                 maxbg=None,
                 random_primer_len = None,
                 num_samples=3,
                 target=None,
                 target_name=None,
                 **kwargs):
        self.mindepth = mindepth
        self.maxbg = maxbg
        self.random_primer_len = random_primer_len
        self.target_name = target_name
        super().__init__(**kwargs)
        self.add(InputNode(name="target",
                           parallel=False))
        if target is not None:
            if isinstance(target, Node):
                connect(target, self.target)
            elif isinstance(target, str):
                self.target.set_file(target)
        self.add(OutputNode(name="profile",
                            parallel=False))

        assert num_samples > 0
        samples = ["modified", "untreated", "denatured"]
        for i in range(num_samples):
            self.add(InputNode(name="counts_{}".format(samples[i]),
                               parallel=False))

        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "make_reactivity_profiles.py"),
               "--fa", "{target}"]
        if self.target_name is not None:
            cmd += ["--rna", '"{}"'.format(self.target_name)]
        cmd += ["--counts"]
        cmd += ["{{{}}}".format(n.get_name()) for n in self.input_nodes
                if n.get_name().startswith("counts")]
        cmd += ["--out", "{profile}"]
        if self.mindepth is not None:
            cmd += ["--mindepth", str(self.mindepth)]
        if self.maxbg is not None:
            cmd += ["--maxbg", str(self.maxbg)]
        if self.random_primer_len is not None and self.random_primer_len > 0:
            cmd += ["--random-primer-len", str(self.random_primer_len)]

        return cmd


class NormProfile(Component):
    def __init__(self,
                 profile=None,
                 profiles=None,
                 target_name=None,
                 target_names=None,
                 **kwargs):
        super().__init__(**kwargs)
        if profile is not None and profiles is not None:
            raise RuntimeError(
                "Error: for NormProfile component __init__(), can specify either profile or profiles, but not both.")
        if profile is not None:
            self.add(InputNode(name="profile",
                               parallel=False))
            if isinstance(profile, str):
                self.profile.set_file(profile)
            else:
                connect(profile, self.profile)
            self.add(OutputNode(name="normed",
                                parallel=False,
                                assoc_rna=target_name))
        elif profiles is not None:
            for i in range(len(profiles)):
                name = "profile_{}".format(i + 1)
                outname = "normed_{}".format(i + 1)
                self.add(InputNode(name=name,
                                   parallel=False))
                if isinstance(profiles[i], str):
                    self.__getattr__(name).set_file(profiles[i])
                else:
                    connect(profiles[i], self.__getattr__(name))
                assoc_rna = None
                try:
                    assoc_rna = target_names[i]
                except IndexError:
                    pass
                self.add(OutputNode(name=outname,
                                    parallel=False,
                                    assoc_rna=assoc_rna))
        else:
            # if no profiles provided, default to single profile input/output
            self.add(InputNode(name="profile",
                               parallel=False))
            self.add(OutputNode(name="normed",
                                parallel=False,
                                assoc_rna=target_name))

        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "normalize_profiles.py"),
               "--warn-on-error", # don't crash if not enough data to normalize
               "--tonorm"]
        for node in self.input_nodes:
            name = node.get_name()
            cmd += ["{{{}}}".format(name)]
        cmd += ["--normout"]
        for node in self.output_nodes:
            if isinstance(node, (StdoutNode, StderrNode)):
                continue
            name = node.get_name()
            cmd += ["{{{}}}".format(name)]
        return cmd

# FIXME: add --primers input for primers file if provided
class RenderFigures(Component):
    def __init__(self,
                 amplicon=False,
                 do_profiles=True,
                 do_histograms=True,
                 mindepth=5000,
                 maxbg=0.05,
                 **kwargs):
        self.amplicon = amplicon
        # TODO: expose the params min_depth_pass_frac, max_high_bg_frac, min_positive
        self.mindepth = mindepth
        self.maxbg = maxbg
        super().__init__(**kwargs)
        self.add(InputNode(name="profile",
                           parallel=False))
        if self.amplicon:
            self.add(InputNode(name="primers",
                               parallel=False))
        if do_profiles:
            self.add(OutputNode(name="profiles_fig",
                                extension="pdf",
                                parallel=False))
        if do_histograms:
            self.add(OutputNode(name="histograms_fig",
                                extension="pdf",
                                parallel=False))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "render_figures.py"),
               "--infile", "{profile}",
               "--mindepth", str(self.mindepth),
               "--maxbg", str(self.maxbg)]
        if self.amplicon:
            cmd += ["--primers", "{primers}"]
        node_names = [n.get_name() for n in self.output_nodes]
        if "profiles_fig" in node_names:
            cmd.extend(["--plot", "{profiles_fig}"])
        if "histograms_fig":
            cmd.extend(["--hist", "{histograms_fig}"])
        if self.assoc_rna is not None:
            cmd.extend(["--title", '"RNA: {}"'.format(self.assoc_rna)])
        return cmd

    def after_run_message(self):
        return self.read_stdout()


class RenderMappedDepths(Component):
    def __init__(self,
                 amplicon=False,
                 **kwargs):
        self.amplicon=amplicon
        super().__init__(**kwargs)
        self.add(InputNode(name="profile",
                           parallel=False))
        if self.amplicon:
            self.add(InputNode(name="primer_locations",
                               parallel=False))
            self.add(OutputNode(name="est_abundances",
                                extension="txt",
                                parallel=False))
        self.add(OutputNode(name="depth_fig",
                            extension="pdf",
                            parallel=False))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "render_mapped_depths.py"),
               "--rna-name", self.assoc_rna,
               "--tsv", "{profile}"]
        if self.amplicon:
            cmd += ["--primer-locations", "{primer_locations}"]
            cmd += ["--estimated-abundances", "{est_abundances}"]
        cmd += [ "--out", "{depth_fig}"]
        return cmd

    def after_run_message(self):
        return self.read_stdout()


class TabToShape(Component):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="profile",
                           parallel=False))
        self.add(OutputNode(name="shape",
                            parallel=False,
                            extension="shape"))
        self.add(OutputNode(name="map",
                            parallel=False,
                            extension="map"))
        self.add(OutputNode(name="varna",
                            parallel=False,
                            extension="txt"))
        self.add(OutputNode(name="ribosketch",
                            parallel=False,
                            extension="txt"))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "tab_to_shape.py"),
               "--infile", "{profile}",
               "--shape", "{shape}",
               "--map", "{map}",
               "--varna", "{varna}",
               "--ribosketch", "{ribosketch}"]
        return cmd


class ProfileHandler(Component):
    """
    Compute reactivity profile from 1-3 samples,
    normalize profile,
    convert to .map and .shape files,
    render pdf summary figures

    """

    def __init__(self,
                 target=None,
                 target_name=None,
                 mindepth=None,
                 maxbg=None,
                 random_primer_len=None,
                 counts=None,
                 norm=None,
                 amplicon=None,
                 **kwargs):
        require_explicit_kwargs(locals())
        super().__init__(**kwargs)
        profilemaker = CalcProfile(target=target,
                                   target_name=target_name,
                                   maxbg=maxbg,
                                   mindepth=mindepth,
                                   random_primer_len=random_primer_len,
                                   num_samples=len(counts))
        self.add(profilemaker)

        self.add(profilemaker.target)

        profilenode = profilemaker.profile
        if norm:
            normer = NormProfile(target_name=target_name)
            self.add(normer)
            connect(profilemaker.profile, normer.profile)
            profilenode = normer.normed

        try:
            connect(counts[0], profilemaker.counts_modified)
            connect(counts[1], profilemaker.counts_untreated)
            connect(counts[2], profilemaker.counts_denatured)
        except IndexError:
            pass

        tabtoshaper = TabToShape()
        self.add(tabtoshaper)
        connect(profilenode, tabtoshaper.profile)

        renderer = RenderFigures(assoc_rna=target_name,
                                 mindepth=mindepth,
                                 maxbg=maxbg,
                                 amplicon=amplicon)
        mapped_depth_renderer = RenderMappedDepths(assoc_rna=target_name,
                                                   amplicon=amplicon)
        self.add([renderer, mapped_depth_renderer])
        connect(profilenode, renderer.profile)
        connect(profilenode, mapped_depth_renderer.profile)


class SequenceCorrector(Component):
    # FIXME: clarify class names (right now very similar names for nested components)

    def __init__(self,
                 target=None,
                 target_name=None,
                 mindepth=None,
                 minfreq=None,
                 **kwargs):
        self.mindepth = mindepth
        self.minfreq = minfreq
        self.target_name = target_name
        super().__init__(**kwargs)
        self.add(InputNode(name="target",
                           parallel=False))
        if target is not None:
            self.target.set_file(target)
        self.add(InputNode(name="variants",
                           parallel=False))
        self.add(OutputNode(name="corrected",
                            extension="fa",
                            parallel=False))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "make_variant.py"),
               "--fa", "{target}"]
        # must provide name of sequence if multiple seqs present
        # in the input file
        if self.target_name is not None:
            cmd += ["--rna", '"{}"'.format(self.target_name)]
        cmd += ["--variants", "{variants}",
               "--out", "{corrected}"]

        if self.mindepth is not None:
            cmd += ["--mindepth", str(self.mindepth)]
        if self.minfreq is not None:
            cmd += ["--minfreq", str(self.minfreq)]
        return cmd

    def after_run_message(self):
        return self.read_stdout()

        # TODO: option to only require a certain depth to identify sequence variants
        # - (need far fewer reads to accurately sequence than accurately SHAPE)
        # - will require implementing some mechanism for detecting depth threshold
        #   met, then terminating upstream modules without triggering downstream
        #   failures


class SplitByTarget(Component):
    def __init__(self,
                 target_names=None,
                 **kwargs):
        require_explicit_kwargs(locals())
        super().__init__(**kwargs)
        self.add(InputNode(name="input"))
        self.add(StdoutNode())
        self.add(StderrNode())
        self.target_names = target_names
        for i in range(len(target_names)):
            self.add(OutputNode(name="rna_{}".format(i + 1),
                                extension="passthrough",
                                assoc_rna=target_names[i]))
            # TODO: store dict of output nodes indexed by target name? less fragile than int index

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "split_by_target.py")]
        cmd += ['-i', "{input}"]
        cmd += ['-n']
        for n in self.target_names:
            cmd += ['"'+n+'"']
        cmd += ['-o']
        for i in range(len(self.target_names)):
            node_name = "rna_{}".format(i + 1)
            cmd += ["{{{}}}".format(node_name)]
        return cmd


class GetSequenceLengths(Component):
    def __init__(self,
                 target_names=None,
                 **kwargs):
        require_explicit_kwargs(locals())
        super().__init__(**kwargs)
        self.add(InputNode(name="fasta",
                           parallel=False))
        self.add(StdoutNode())
        self.add(StderrNode())
        self.target_names = target_names
        for i in range(len(target_names)):
            self.add(OutputNode(name="L{}".format(i + 1),
                                extension="",
                                assoc_rna=target_names[i],
                                parallel=False))

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "get_sequence_lengths.py")]
        cmd += ['--fa', "{fasta}"]
        cmd += ['--out']
        for i in range(len(self.target_names)):
            node_name = "L{}".format(i + 1)
            cmd += ["{{{}}}".format(node_name)]
        return cmd


class PrimerLocator(Component):
    def __init__(self,
                 fastas=None,
                 target_names=None,
                 primer_files=None,
                 primers_in_sequence=None,
                 **kwargs):
        require_explicit_kwargs(locals())
        super().__init__(**kwargs)
        self.fastas = fastas
        self.target_names = target_names
        self.primer_files = primer_files
        self.primers_in_sequence = primers_in_sequence
        # FIXME: maybe move some of this boilerplate multinode input/output stuff to Component
        for i in range(len(fastas)):
            self.add(InputNode(name="fasta_{}".format(i+1),
                               filename=fastas[i]))
        for i in range(len(primer_files)):
            self.add(InputNode(name="primers_{}".format(i+1),
                               filename=primer_files[i]))
        for i in range(len(target_names)):
            self.add(OutputNode(name="locs_{}".format(i+1),
                                extension="txt",
                                assoc_rna=target_names[i],
                                parallel=False))
            # will be linked to ParameterNode of MutationCounter
            self.add(OutputNode(name="n_pairs_{}".format(i+1),
                                extension="txt",
                                assoc_rna=target_names[i],
                                parallel=False))
        #self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(this_dir, "locate_primers.py")]
        cmd += ["--fastas"]
        for i in range(len(self.fastas)):
            node_name = "fasta_{}".format(i + 1)
            cmd += ["{{{}}}".format(node_name)]
        if len(self.primer_files) > 0:
            cmd += ["--primer-files"]
            for i in range(len(self.primer_files)):
                node_name = "primers_{}".format(i+1)
                cmd += ["{{{}}}".format(node_name)]
        if self.primers_in_sequence:
            cmd += ["--primers-in-sequence"]
        cmd += ["--target-names"]
        for n in self.target_names:
            cmd += ['"' + n + '"']
        cmd += ["--locations-out"]
        for i in range(len(self.target_names)):
            node_name = "locs_{}".format(i + 1)
            cmd += ["{{{}}}".format(node_name)]
        cmd += ["--n-pairs-out"]
        for i in range(len(self.target_names)):
            node_name = "n_pairs_{}".format(i + 1)
            cmd += ["{{{}}}".format(node_name)]
        return cmd


class SplitToFile(Component):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add(StdinNode())
        self.add(StdoutNode(extension="passthrough"))
        self.add(OutputNode(extension="passthrough",
                            parallel=False,
                            name="to_file"))
        self.add(StderrNode())

    def cmd(self):
        out_node = None
        # parallel output node name gets changed before process started,
        # so locate it by elimination
        for node in self.output_nodes:
            if node.get_name() not in ["to_file", "stderr"]:
                out_node = node
        out_node_name = out_node.get_name()
        cmd = "tee {to_file} >"
        cmd += "{{{}}}".format(out_node_name)
        cmd += " <{stdin}"
        return cmd


# FIXME: clean this up somehow, maybe need more reflection to eliminate
#        ambiguities when trying to access specific components and/or nodes
#        - sometimes want to access nodes by name regardless of internal structure or
#          presence of wrapper
#        - sometimes need to know node location within wrapper
def split_to_file_wrapper(component, selected_out_names=None):
    """
    "Factory" for creating a wrapper component that
    will split all parallel outputs of a Component
    instance into two: one pipe and one file. The new
    component will have the same input and output nodes
    as the original component.

    """



    wrapper = Component()
    wrapper.name = component.get_name()
    wrapper.assoc_rna = component.assoc_rna
    wrapper.assoc_sample = component.assoc_sample
    wrapper.parent_component = component.parent_component
    wrapper.add([component])
    for node in component.input_nodes:
        wrapper.add(node)
    split_nodes = []
    for node in component.output_nodes:
        if selected_out_names is not None:
            if node.get_name() in selected_out_names:
                split_nodes.append(node)
        elif node.parallel and node.get_name() not in ["stdout", "stderr"]:
            split_nodes.append(node)
        else:
            wrapper.add(node)
    for i, node in enumerate(split_nodes):
        comp = SplitToFile(name="SplitToFile{}".format(i+1))
        wrapper.add(comp)
        try:
            connected_node = node.output_nodes[0]
        except IndexError:
            connected_node = None
        # TODO: update the disconnect() func to support single node arg?
        node.output_nodes = []
        connect(node, comp.stdin)
        # attempt to rewire/preserve existing connections
        if connected_node is not None:
            connect(comp.stdout, connected_node)
        wrapper.add(comp.stdout)
        comp.stdout.parallel = True
        # comp.stdout.name = "parallel_out"
        comp.stdout.name = node.get_name()

    # set assoc_rna property for all children
    for c in wrapper.collect_components():
        c.assoc_rna = wrapper.assoc_rna
        c.assoc_sample = wrapper.assoc_sample
    for n in wrapper.collect_component_nodes():
        n.assoc_rna = wrapper.assoc_rna
        n.assoc_sample = wrapper.assoc_sample
    return wrapper


class AlignPrep(Component):
    def __init__(self,
                 target=None,
                 num_targets=None,
                 total_target_length=None,
                 star_aligner=None,
                 genomeSAindexNbase=None,
                 nproc=None,
                 **kwargs):
        require_explicit_kwargs(locals())
        if isinstance(target, list) and len(target)==1:
            # allow init with list of 1 target by re-initializing with that target
            self.__init__(target=target[0],
                          num_targets=num_targets,
                          total_target_length=total_target_length,
                          star_aligner=star_aligner,
                          nproc=nproc,
                          genomeSAindexNbase=genomeSAindexNbase,
                          **kwargs)
            return
        super().__init__(**kwargs)

        if star_aligner is not None and star_aligner:
            indexbuilder = StarIndexBuilder(num_targets=num_targets,
                                            total_target_length=total_target_length,
                                            nproc=nproc,
                                            genomeSAindexNbase=genomeSAindexNbase)
        else:
            indexbuilder = BowtieIndexBuilder()
        if isinstance(target, str):
            fastachecker = FastaFormatChecker(fasta=target)
            self.add(fastachecker)
            connect(fastachecker.fasta.input_node, indexbuilder.target)
        elif isinstance(target, Node):
            fastachecker = FastaFormatChecker(fasta=target)
            connect(fastachecker.fasta.input_node, indexbuilder.target)
            self.add(fastachecker)
        elif isinstance(target, list):
            # if multiple target files, combine into single file before Bowtie index build
            fastacombine = Appender(inputs=target,
                                    add_extra_newline=True)
            for i in range(len(fastacombine.input_nodes)):
                fastachecker = FastaFormatChecker(name="FastaFormatChecker_{}".format(i+1))
                connect(fastacombine.input_nodes[i].input_node, fastachecker.fasta)
                self.add(fastachecker)
            self.add(fastacombine)
            connect(fastacombine.appended, indexbuilder.target)

        self.add(indexbuilder)
        self.add(indexbuilder.target)
        self.add(indexbuilder.index)


class Sample(Component):
    def __init__(self,
                 R1=None, R2=None,
                 U=None,
                 o1=None, o2=None,
                 assoc_sample=None,
                 assoc_rna=None,
                 min_qual_to_trim=None,
                 window_to_trim=None,
                 min_length_to_trim=None,
                 preserve_order=None,
                 total_target_length=None,
                 star_aligner=None,
                 star_shared_index=None,
                 min_mapq=None,
                 disable_soft_clipping=None,
                 nproc=None,
                 maxins=None,
                 max_search_depth=None,
                 max_reseed=None,
                 **kwargs):
        """
        Note: (does not perform aligner index building, mutation parsing/counting,
         or reactivity profile creation)

        Args:
            R1: path to FASTQ file for forward reads or list of
                paths to FASTQ files
            R2: path to FASTQ file for reverse reads or list of
                paths to FASTQ files
            U:  path to unpaired FASTQ file or list of paths to
                unpaired FASTQ files (must provide either U or R1+R2,
                but not both)

        """
        if not ((R1 is not None and R2 is not None) or
                U is not None):
            msg = "Error: must provide either R1+R2 or U, not both"
            raise RuntimeError(msg)

        super().__init__(**kwargs)

        self.total_target_length = total_target_length
        self.o1 = o1
        self.o2 = o2
        self.min_mapq = min_mapq

        # check if single fastq or multiple fastqs
        fastq_list = False
        if R1 is not None and isinstance(R1, list):
            fastq_list = True
            if len(R1)==1:
                fastq_list = False
                R1 = R1[0]
                R2 = R2[0]
        elif U is not None and isinstance(U, list):
            fastq_list = True
            if len(U)==1:
                fastq_list = False
                U = U[0]

        aligner = None
        if U is not None:
            # unpaired reads
            if fastq_list:
                # add component to concatenate input files into single streams
                append = Appender(inputs=U)
                progmonitor = ProgressMonitor()
                qtrimmer = QualityTrimmer(min_qual=min_qual_to_trim,
                                          window=window_to_trim,
                                          min_length=min_length_to_trim)
                connect(append, progmonitor)
                self.add([append,
                          progmonitor,
                          qtrimmer])
            else:
                progmonitor = ProgressMonitor(input=U)
                qtrimmer = QualityTrimmer(min_qual=min_qual_to_trim,
                                          window=window_to_trim,
                                          min_length=min_length_to_trim)
                self.add([progmonitor,
                          qtrimmer])
            connect(progmonitor, qtrimmer)

            aligner_params = {"reorder": preserve_order,
                              "assoc_rna": assoc_rna,
                              "disable_soft_clipping": disable_soft_clipping,
                              "nproc": nproc}

            if star_aligner is not None and star_aligner:
                aligner = StarAligner(name="StarAligner",
                                      paired=False,
                                      **aligner_params)
                connect(qtrimmer.trimmed, aligner.fastq)
            else:
                aligner = BowtieAlignerMixedInput(maxins=maxins,
                                                  max_search_depth=max_search_depth,
                                                  max_reseed=max_reseed,
                                                  **aligner_params)
                connect(qtrimmer.trimmed, aligner.interleaved_fastq)
            self.add(aligner)

        else:
            # paired reads
            # TODO: might be simpler to do this with nested components
            if fastq_list:
                # add components to concatenate input files into single streams
                append1 = Appender(name="Appender1", inputs=R1)
                append2 = Appender(name="Appender2", inputs=R2)
                progmonitor = ProgressMonitor()
                qtrimmer1 = QualityTrimmer(name="QualityTrimmer1",
                                           min_qual=min_qual_to_trim,
                                           window=window_to_trim,
                                           min_length=min_length_to_trim)
                qtrimmer2 = QualityTrimmer(name="QualityTrimmer2",
                                           min_qual=min_qual_to_trim,
                                           window=window_to_trim,
                                           min_length=min_length_to_trim)
                connect(append1, progmonitor)
                connect(progmonitor, qtrimmer1)
                connect(append2, qtrimmer2)
                self.add([append1,
                          append2,
                          progmonitor,
                          qtrimmer1,
                          qtrimmer2])
            else:
                progmonitor = ProgressMonitor(input=R1)
                qtrimmer1 = QualityTrimmer(name="QualityTrimmer1",
                                           min_qual=min_qual_to_trim,
                                           window=window_to_trim,
                                           min_length=min_length_to_trim)
                connect(progmonitor, qtrimmer1)
                qtrimmer2 = QualityTrimmer(name="QualityTrimmer2",
                                           min_qual=min_qual_to_trim,
                                           window=window_to_trim,
                                           min_length=min_length_to_trim,
                                           fastq=R2)
                self.add([progmonitor,
                          qtrimmer1,
                          qtrimmer2])
            

            interleaver = Interleaver()
            merger = Merger(preserve_order=preserve_order)
            connect(qtrimmer1.trimmed, interleaver.R1)
            connect(qtrimmer2.trimmed, interleaver.R2)
            connect(interleaver, merger)
            self.add([interleaver])
            self.add([merger])

            aligner_params = {"reorder":preserve_order,
                              "assoc_rna":assoc_rna,
                              "disable_soft_clipping":disable_soft_clipping,
                              "nproc":nproc}

            if star_aligner is not None and star_aligner:
                aligner = StarAlignerMixedInput(star_shared_index=star_shared_index,
                                                **aligner_params)
                connect(merger.output, aligner.interleaved_fastq)
                self.add(aligner)
            else:
                aligner = BowtieAlignerMixedInput(maxins=maxins,
                                                  max_search_depth=max_search_depth,
                                                  max_reseed=max_reseed,
                                                  **aligner_params)
                connect(merger.output, aligner.interleaved_fastq)
                self.add(aligner)

        if fastq_list:
            # sum sizes of input files so pipeviewer knows how much total data
            # will be coming through
            s = 0
            if U is not None:
                if isinstance(U, list):
                    for f in U:
                        s += os.path.getsize(f)
                else:
                    s += os.path.getsize(U)
            else:
                if isinstance(R1, list):
                    for f in R1:
                        s += os.path.getsize(f)
                else:
                    s += os.path.getsize(R1)
            progmonitor.expected_bytes = s

        # allow access to some nodes from top level object scope
        self.add(aligner.index, alias="index")
        self.add(aligner.aligned, alias="aligned")

        self.progmon = progmonitor

        # set assoc_sample for all children
        for c in self.collect_components():
            c.assoc_sample = assoc_sample
        for n in self.collect_component_nodes():
            n.assoc_sample = assoc_sample



class PostAlignment(Component):
    """
    Contains components to parse and count mutations and create reactivity
    profiles for a single RNA

    Args:
        target_name:
        target_length:
        num_samples:
        output_variant_counts: True or False (enable file output)
        output_mutation_counts:
        min_mapq:
        min_depth:
        max_bg:
        norm: normalize individual profile enabled
        right_align_ambig_dels:
        right_align_ambig_ins:
        min_mutation_separation:
        min_qual_to_count:
        random_primer_len:
    """

    def __init__(self,
                 target_name=None,
                 target_length=None,
                 primer_pairs=None,
                 target=None,
                 num_samples=None,
                 output_variant_counts=None,
                 output_mutation_counts=None,
                 min_mapq=None,
                 maxins=None,
                 input_is_unpaired=None,
                 min_depth=None,
                 max_bg=None,
                 norm=None,
                 separate_ambig_counts=None,
                 right_align_ambig_dels=None,
                 right_align_ambig_ins=None,
                 min_mutation_separation=None,
                 min_qual_to_count=None,
                 mutation_type_to_count=None,
                 random_primer_len=None,
                 amplicon=None,
                 max_primer_offset=None,
                 require_forward_primer_mapped=None,
                 require_reverse_primer_mapped=None,
                 trim_primers=None,
                 render_mutations=None,
                 render_must_span=None,
                 max_pages=None,
                 per_read_histograms=None,
                 **kwargs):
        require_explicit_kwargs(locals())
        assert isinstance(num_samples, int)
        super().__init__(**kwargs)
        self.assoc_rna = target_name
        samples = ["Modified", "Untreated", "Denatured"]
        # TODO: add kwargs to connect input nodes here?
        counts = []
        ambig_counts = []

        debug_out = False
        if render_mutations:
            debug_out = render_mutations

        for i in range(num_samples):
            sample = samples[i]

            parser = MutationParser(name="MutationParser_" + sample,
                                    min_mapq=min_mapq,
                                    maxins=maxins,
                                    input_is_unpaired=input_is_unpaired,
                                    right_align_ambig_dels=right_align_ambig_dels,
                                    right_align_ambig_ins=right_align_ambig_ins,
                                    min_mutation_separation=min_mutation_separation,
                                    min_qual=min_qual_to_count,
                                    random_primer_len=random_primer_len,
                                    mutation_type=mutation_type_to_count,
                                    variant_mode=output_variant_counts,
                                    amplicon=amplicon,
                                    max_primer_offset=max_primer_offset,
                                    require_forward_primer_mapped=require_forward_primer_mapped,
                                    require_reverse_primer_mapped=require_reverse_primer_mapped,
                                    trim_primers=trim_primers,
                                    assoc_sample=sample,
                                    debug_out=debug_out)

            self.add(parser)

            if debug_out:
                comp_name = "MutationRenderer_{}_{}".format(
                    sample,
                    self.assoc_rna)
                renderer = MutationRenderer(name=comp_name,
                                            assoc_rna=self.assoc_rna,
                                            assoc_sample=sample,
                                            maxins=maxins,
                                            amplicon=amplicon,
                                            max_pages=max_pages,
                                            span=render_must_span)
                connect(parser.debug_out, renderer.input)
                self.add(renderer)

            counter = MutationCounter(name="MutationCounter_" + sample,
                                      target_length=target_length,
                                      primer_pairs=primer_pairs,
                                      variant_out=output_variant_counts,
                                      mutations_out=output_mutation_counts,
                                      assoc_sample=sample,
                                      per_read_histograms=per_read_histograms,
                                      separate_ambig_counts=separate_ambig_counts)

            connect_nodes(parser.parsed_mutations, counter.mut)
            self.add(counter)

            counts.append(counter.mutations)

        profilehandler = ProfileHandler(target=target,
                                        target_name=target_name,
                                        mindepth=min_depth,
                                        maxbg=max_bg,
                                        random_primer_len=random_primer_len,
                                        counts=counts,
                                        norm=norm,
                                        amplicon=amplicon)
        self.add(profilehandler)

        # set assoc_rna property for all children
        for c in self.collect_components():
            c.assoc_rna = target_name
        for n in self.collect_component_nodes():
            n.assoc_rna = target_name


class CorrectSequence(Component):
    def __init__(self,
                 target=None,
                 target_names=None,
                 target_lengths=None,
                 U=None,
                 R1=None,
                 R2=None,
                 preserve_order=None,
                 total_target_length=None,
                 disable_soft_clipping=None,
                 min_mapq=None,
                 min_qual_to_trim=None,
                 window_to_trim=None,
                 min_length_to_trim=None,
                 min_qual_to_count=None,
                 min_seq_depth=None,
                 min_freq=None,
                 random_primer_len=None,
                 star_aligner=None,
                 genomeSAindexNbase=None,
                 star_shared_index=None,
                 nproc=None,
                 maxins=None,
                 max_search_depth=None,
                 max_reseed=None,
                 amplicon=None,
                 max_primer_offset=None,
                 require_forward_primer_mapped=None,
                 require_reverse_primer_mapped=None,
                 trim_primers=None,
                 **kwargs):
        locs = locals()
        for x in ["U", "R1", "R2"]:
            locs.pop(x)
        # FIXME: require either U or R1+R2 but not both
        require_explicit_kwargs(locs)
        super().__init__(**kwargs)

        prep = AlignPrep(target=target,
                         num_targets=len(target_names),
                         total_target_length=total_target_length,
                         star_aligner=star_aligner,
                         genomeSAindexNbase=genomeSAindexNbase,
                         nproc=nproc)

        sample = Sample(U=U,
                        R1=R1, R2=R2,
                        o1=True, o2=True,
                        disable_soft_clipping=disable_soft_clipping,
                        min_mapq=min_mapq,
                        assoc_sample="sequence correction",
                        preserve_order=preserve_order,
                        star_aligner=star_aligner,
                        star_shared_index=star_shared_index,
                        nproc=nproc,
                        maxins=maxins,
                        max_search_depth=max_search_depth,
                        max_reseed=max_reseed,
                        min_qual_to_trim=min_qual_to_trim,
                        window_to_trim=window_to_trim,
                        min_length_to_trim=min_length_to_trim,
                        )
        self.add([prep,
                  sample])
        connect(prep.index, sample.index)
        self.add(prep.target) # simplify access to main sequence input node

        # split aligned output by mapped RNA (if needed)
        if len(target_names)>1:
            splitter = SplitByTarget(target_names=target_names)
            connect(sample.aligned, splitter.input)
            self.add(splitter)
            corrected_nodes = []
            for i in range(len(target_names)):
                mapped_node = splitter["rna_{}".format(i+1)]
                parser = MutationParser(name="MutationParser_{}".format(i+1),
                                        assoc_rna=target_names[i],
                                        min_mapq=min_mapq,
                                        min_qual=min_qual_to_count,
                                        random_primer_len=random_primer_len,
                                        maxins=maxins,
                                        amplicon=amplicon,
                                        max_primer_offset=max_primer_offset,
                                        require_forward_primer_mapped=require_forward_primer_mapped,
                                        require_reverse_primer_mapped=require_reverse_primer_mapped,
                                        trim_primers=trim_primers,
                                        )

                counter = MutationCounter(name="MutationCounter_{}".format(i+1),
                                          assoc_rna=target_names[i],
                                          target_length=target_lengths[i],
                                          variant_out=True,
                                          mutations_out=False,
                                          )
                connect(splitter["rna_{}".format(i+1)], parser.input)
                connect(parser, counter)
                self.add([parser,
                          counter])

                sequencefixer = SequenceCorrector(name="SequenceCorrector_{}".format(i+1),
                                                  assoc_rna=target_names[i],
                                                  mindepth=min_seq_depth,
                                                  minfreq=min_freq,
                                                  target_name=target_names[i])  # FIXME: keep param names consistent across codebase
                self.add(sequencefixer)
                connect(prep.target.input_node, sequencefixer.target)
                connect(counter.variants, sequencefixer.variants)
                corrected_nodes.append(sequencefixer.corrected)
            # combine all corrected fasta sequences into a single file
            appender = Appender(inputs=corrected_nodes,
                                add_extra_newline=True)
            self.add(appender)
            self.add(appender.appended, alias="corrected")

        else:
            parser = MutationParser(min_mapq=min_mapq,
                                    min_qual=min_qual_to_count,
                                    random_primer_len=random_primer_len,
                                    maxins=maxins,
                                    amplicon=amplicon,
                                    max_primer_offset=max_primer_offset,
                                    require_forward_primer_mapped=require_forward_primer_mapped,
                                    require_reverse_primer_mapped=require_reverse_primer_mapped,
                                    trim_primers=trim_primers,
                                    )
            counter = MutationCounter(variant_out=True,
                                      mutations_out=False,
                                      target_length=target_lengths[0])
            connect(sample.aligned, parser.input)
            connect(parser, counter)
            self.add([parser,
                      counter])

            sequencefixer = SequenceCorrector(mindepth=min_seq_depth,
                                              minfreq=min_freq)  # FIXME: keep param names consistent across codebase
            self.add(sequencefixer)
            connect(prep.target.input_node, sequencefixer.target)
            connect(counter.variants, sequencefixer.variants)
            self.add(sequencefixer.corrected)

        lengthgetter = GetSequenceLengths(target_names=target_names)
        connect(self.corrected, lengthgetter.fasta)
        self.add(lengthgetter)
        self.add(lengthgetter.collect_component_nodes(name="L*"))


class Mangler(Component):
    def __init__(self,
                 **kwargs):
        gv_props = {"color": "red",
                    "penwidth": "8"}
        super().__init__(gv_props=gv_props,
                         **kwargs)
        self.add(InputNode())
        self.add(OutputNode(extension="passthrough"))

    def cmd(self):
        ext = self.input.get_extension()
        if any([ext.endswith(x) for x in ["bam", "gz"]]):
            mangler = "mangle_binary.py"
        else:
            mangler = "mangle_text_fixed.py"
        cmd = [pyexe,
               os.path.join(bin_dir, mangler),
               '<',
               '{input}',
               '>',
               '{output}']
        return cmd


class PsToPdf(Component):
    """
    Convert a PostScript file to PDF.
    """
    def __init__(self,
                 maxins=200,
                 **kwargs):
        super().__init__(**kwargs)
        self.maxins=maxins
        self.add(InputNode(name="ps"))
        self.add(OutputNode(name="pdf",
                            extension="pdf",
                            parallel=False))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [os.path.join(bin_dir, "ps2pdf_rescaled.sh"),
               str(self.maxins),
               "{ps}",
               "{pdf}"]
        return cmd


class MutationRendererPs(Component):
    """
    For debugging, make a postscript file with parsed and processed
    reads and mutations
    """
    def __init__(self,
                 maxins=200,
                 amplicon=False,
                 max_pages=100,
                 span=None,
                 **kwargs):
        super().__init__(**kwargs)
        self.maxins = maxins
        self.amplicon = amplicon
        self.max_pages = max_pages
        self.span = span
        self.add(InputNode(name="input",
                           parallel=True))
        if self.amplicon:
            self.add(InputNode(name="primers"))
        self.add(OutputNode(name="ps",
                            extension="ps",
                            parallel=True))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "render_mutations_ps.py"),
               "--max-length", str(self.maxins),
               "--max-pages", str(self.max_pages),
               "--input", "{input}",
               "--output", "{ps}"]
        if self.amplicon:
            cmd += ["--primers", "{primers}"]
        if self.span is not None and self.span != '':
            cmd += ["--must-span", self.span]
        return cmd


class MutationRenderer(Component):
    def __init__(self,
                 maxins=200,
                 amplicon=False,
                 max_pages=100,
                 span=None,
                 **kwargs):
        super().__init__(**kwargs)
        renderer = MutationRendererPs(maxins=maxins, amplicon=amplicon, max_pages=max_pages, span=span)
        converter = PsToPdf(maxins=maxins)
        connect(renderer.ps, converter.ps)

        self.add([renderer,
                  converter])

        # allow access to some nodes from this object scope
        self.add(renderer.input, alias="input")
        self.add(converter.pdf, alias="pdf")




