"""
Specific wrappers for pipeline executables, including
the Bowtie2 sequence aligner, BBmerge, and compiled shapemapper
modules for read trimming, mutation parsing, and mutation
counting.
"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import subprocess
from math import log2

from pyshapemap.component import *
from pyshapemap.util import require_explicit_kwargs, sanitize
from pyshapemap.nodes import connect as connect_nodes
from pyshapemap.nodes import disconnect as disconnect_nodes

this_dir = os.path.dirname(os.path.realpath(__file__))
bin_dir = os.path.join(this_dir, "../../bin")
pyexe = "python3"


# TODO: may be a better way to handle these func names - import order might confuse things in other files
# TODO: move these functions to component.py?
def connect(x, y):
    """
    Convenience function to link two Components with only
    a single input/output node (excluding stdout and stderr). Will also
    connect two nodes using connect() defined in nodes.py.
    """
    if isinstance(x, Node) and isinstance(y, Node):
        connect_nodes(x, y)
    elif isinstance(x, Component) and isinstance(y, Component):
        output_nodes = [node for node in x.output_nodes if node.get_name() not in ["stdout", "stderr"]]
        if len(output_nodes) != 1 or len(y.input_nodes) != 1:
            raise TypeError(
                "To connect() two Components, the first must have only one output node (excluding stdout and stderr), and the second must have only one input node.")
        connect_nodes(output_nodes[0], y.input_nodes[0])
    else:
        raise TypeError("Objects to connect() must both be Nodes or Components")


def disconnect(x, y):
    if isinstance(x, Node) and isinstance(y, Node):
        disconnect_nodes(x, y)
    elif isinstance(x, Component) and isinstance(y, Component):
        output_nodes = [node for node in x.output_nodes if node.get_name() not in ["stdout", "stderr"]]
        if len(output_nodes) != 1 or len(y.input_nodes) != 1:
            raise TypeError(
                "To disconnect() two Components, the first must have only one output node (excluding stdout and stderr), and the second must have only one input node.")
        disconnect_nodes(output_nodes[0], y.input_nodes[0])
    else:
        raise TypeError("Objects to connect() must both be Nodes or Components")


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


class Interleaver(Component):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="R1"))
        self.add(InputNode(name="R2"))
        self.add(OutputNode(name="interleaved",
                            extension="fastq"))
        self.add(StderrNode())

    def cmd(self):
        cmd = "paste <(paste - - - - < {R1}) "
        cmd += "<(paste - - - - < {R2}) "
        cmd += "| tr '\\t' '\\n' "
        cmd += "> {interleaved}"
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
    def __init__(self,
                 preserve_order=None,
                 **kwargs):
        self.preserve_order = preserve_order
        self.use_unmerged = True
        super().__init__(**kwargs)
        self.add(InputNode(name="interleaved_fastq"))
        self.add(OutputNode(name="merged",
                            extension="fastq"))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = ["bbmerge.sh",
               "vstrict=t",
               "in=stdin",
               "out={merged}",
               "interleaved=t",
               "usejni=t"]  # FIXME: autodetect whether JNI components are compiled
        if self.use_unmerged:
            cmd += ["mix=t"]
        if self.preserve_order:
            cmd += ["ordered=t"]
        cmd += ["<", "{interleaved_fastq}"]
        return cmd

    # ugly workaround: for paths containing spaces, bbmerge apparently requires them input as
    #                  quoted strings also containing escaped spaces
    # - This means we need to override the default command formatting just for the merged output node
    def format_command(self,
                       command):
        # rename {merged} temporarily so no error when calling super format
        if isinstance(command, list):
            for i in range(len(command)):
                command[i] = command[i].replace("{merged}","{{merged}}")
        else:
            command = command.replace("{merged}","{{merged}}")
        formatted = super().format_command(command)
        escaped_filename = '"'+self.merged.output_nodes[0].filename.replace(' ', '\ ')+'"'
        formatted = formatted.format(merged=escaped_filename)
        return formatted

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


class LengthFilter(Component):
    """
    Workaround to throw away 1-length "reads" that are
    making it through BBmerge for some reason
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="fastq"))
        self.add(OutputNode(name="filtered",
                            extension="fastq"))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = "shapemapper_read_length_filter -i {fastq} -o {filtered}"
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
                 **kwargs):
        self.reorder = reorder
        self.nproc = nproc
        self.disable_soft_clipping = disable_soft_clipping
        super().__init__(**kwargs)
        self.add(InputNode(name="index",
                           parallel=False))
        self.add(StdinNode(name="fastq"))
        self.add(OutputNode(name="aligned",
                            extension="sam"))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = ["bowtie2"]
        cmd += ["-p", str(self.nproc)]
        if not self.disable_soft_clipping:
            cmd += ["--local", "--sensitive-local"]
        else:
            cmd += ["--end-to-end", "--sensitive"]
        cmd += [
               "--mp", "3,1", # try to be slightly more consistent with STAR's alignment parameters
               "--rdg", "5,1", # multinuc deletions are a large part of the signal, so don't penalize
                               # as much as bowtie's default setting (open penalty 5, extension penalty 3)
               "--rfg", "5,1",
               "--dpad", "30", # this seems reasonable to allow fairly large deletions, but not too crazy (although these do exist)
               "--maxins", "800",
               "--ignore-quals",
               "--no-unal",  # don't produce SAM records for reads that didn't map
               ]

        if self.reorder:
            cmd += ["--reorder"] # output in same order as input for debugging purposes
        cmd += ["-x", "{index}",
               "-U", "{fastq}",
               "-S", "{aligned}"]
        cmd = ' '.join(cmd)
        return cmd

    def after_run_message(self):
        # Bowtie2 alignment stats are written to its stderr
        return self.read_stderr()

class StarIndexBuilder(Component):
    out_extension = ""

    def __init__(self,
                 target=None,
                 num_targets=None,
                 total_target_length=None,
                 nproc=None,
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

    def cmd(self):
        cmd = ["STAR",
               "--runMode", "genomeGenerate",
               "--genomeDir", "{index}",
               "--outTmpDir", "{temp}",
               "--outFileNamePrefix", "{logs}/",
               "--genomeFastaFiles", "{target}",
               "--runThreadN", str(self.nproc)]
        # The STAR manual recommends setting genomeSAindexNbases to
        # min(14, log2(GenomeLength)/2 - 1) to handle "short" reference
        # sequences (e.g. bacterial genomes or smaller), so the called script calculates
        # this value automatically.
        # STAR seems to segfault if genomeSAindexNbases is too high, even by a small amount
        # The manual also recommends setting genomeChrBinNbits to 
        # min(18, log2(GenomeLength/NumberOfReferences))
        index_n = int(round(min(14, log2(self.total_target_length)/2.0-1)))
        cmd += ["--genomeSAindexNbases", str(index_n)]
        bin_n = int(min(18, log2(self.total_target_length)/self.num_targets))
        cmd += ["--genomeChrBinNbits", str(bin_n)]
        return cmd

# WARNING: unlike Bowtie2, STAR will align to all indices in a directory
# WARNING: STAR's performance degrades significantly if target sequences
#          are incomplete (for example, only providing the small subunit 
#          ribosomal sequence when both subunits are present in the data)
class StarAligner(Component):
    def __init__(self,
                 reorder=False,
                 min_target_length=None,
                 disable_soft_clipping=False,
                 nproc=None,
                 **kwargs):
        require_explicit_kwargs(locals())
        self.reorder = reorder
        self.nproc = nproc
        self.min_target_length = min_target_length
        self.disable_soft_clipping = disable_soft_clipping
        super().__init__(**kwargs)
        self.add(InputNode(name="index",
                           isfolder=True))
        self.add(StdinNode(name="fastq"))
        self.add(OutputNode(name="temp",
                            isfolder=True,
                            make_parent=True,
                            error_on_existing=True))
        self.add(OutputNode(name="logs",
                            isfolder=True))
        self.add(OutputNode(name="aligned",
                            extension="sam"))
        self.add(StderrNode())

    def cmd(self):
        cmd = ["STAR"]
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
        cmd += ["--runMode", "alignReads",
                # only report one of the top alignments for a read that maps to multiple locations
                "--outMultimapperOrder", "Random",
                "--outSAMmultNmax", "1",
                "--outStd", "SAM",
                "--outSAMattributes", "MD",
                "--genomeDir", "{index}",
                "--readFilesIn", "{fastq}",
                 "--outTmpDir", "{temp}",
                "--outFileNamePrefix", "{logs}/",
                "> ", "{aligned}"]
        return cmd

    def after_run_message(self):
        return open(os.path.join(self.logs.output_nodes[0].foldername, "Log.final.out"), "rU").read()

class SamToBam(Component):
    # NOTE: currently unused
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add(InputNode(name="sam"))
        self.add(OutputNode(name="bam",
                            extension="bam"))

        self.add(StderrNode())

    def cmd(self):
        cmd = "samtools view -b -o {bam} {sam}"
        return cmd


class MutationParser(Component):
    def __init__(self,
                 min_mapq=None,
                 **kwargs):
        self.min_mapq = 35
        super().__init__(**kwargs)
        self.min_mapq = min_mapq
        self.add(InputNode(name="input"))
        self.add(OutputNode(name="parsed_mutations",
                            extension="mut"))
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = "shapemapper_mutation_parser -i {input} -o {parsed_mutations} -w"
        if self.min_mapq is not None:
            cmd += " -m {}".format(self.min_mapq)
        return cmd


class MutationCounter(Component):
    def __init__(self,
                 target_length=None,
                 variant_out=None,
                 mutations_out=None,
                 classified_out=None,
                 separate_ambig_counts=None,
                 right_align_ambig_dels=None,
                 right_align_ambig_ins=None,
                 random_primer_len=None,
                 min_mutation_separation=None,
                 min_qual=None,
                 mutation_type=None,
                 **kwargs):
        self.separate_ambig_counts = separate_ambig_counts
        self.right_align_ambig_dels = right_align_ambig_dels
        self.right_align_ambig_ins = right_align_ambig_ins
        self.random_primer_len = random_primer_len
        self.min_mutation_separation = min_mutation_separation
        self.min_qual = min_qual
        self.mutation_type = mutation_type
        super().__init__(**kwargs)
        self.add(InputNode(name="mut"))

        # target_length is either an int parameter, or an OutputNode outputting 
        # a file containing the target_length parameter (this is to allow
        # this component to get updated sequence length after sequence correction)
        if isinstance(target_length, int):
            self.target_length = target_length
        elif isinstance(target_length, OutputNode):
            self.add(ParameterNode(name="target_length"))
            connect_nodes(target_length, self.target_length)

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
        if classified_out is not None and classified_out:
            kw = {"name": "classified_mutations"}
            if isinstance(classified_out, str):
                kw["filename"] = classified_out
            self.add(OutputNode(**kw))

        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = ["shapemapper_mutation_counter",
               "-i", "{mut}", "-w"]
        node_names = [n.get_name() for n in self.output_nodes]
        if "variants" in node_names:
            cmd += ["-v", "{variants}"]
        if "mutations" in node_names:
            cmd += ["-c", "{mutations}"]
        if "classified_mutations" in node_names:
            cmd += ["-l", "{classified_mutations}"]
        if self.separate_ambig_counts is not None and self.separate_ambig_counts:
            cmd += ["--separate_ambig_counts"]
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
        if self.target_length is not None:
            cmd += ["--length", "{target_length}"]
        if self.mutation_type != '':
            cmd += ["--use_only_mutation_type", "{mutation_type}"]
        return cmd

    def after_run_message(self):
        lines = self.read_stdout().splitlines()
        # just display the histogram tables
        i = 0
        start_i = 0
        for i in range(len(lines)):
            if lines[i] == "Read lengths":
                start_i = i
                break
        end_i = 0
        for i in range(len(lines)-1,-1,-1):
            if lines[i] == "--------------------":
                end_i = i
                break
        return "\n".join(lines[start_i:end_i+1])


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


class RenderFigures(Component):
    def __init__(self,
                 do_profiles=True,
                 do_histograms=True,
                 mindepth=5000,
                 maxbg=0.05,
                 **kwargs):
        # TODO: expose the params min_depth_pass_frac, max_high_bg_frac, min_positive
        self.mindepth = mindepth
        self.maxbg = maxbg
        super().__init__(**kwargs)
        self.add(InputNode(name="profile",
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
        self.add(StdoutNode())
        self.add(StderrNode())

    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "tab_to_shape.py"),
               "--infile", "{profile}",
               "--shape", "{shape}",
               "--map", "{map}"]
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
                                 maxbg=maxbg)
        self.add(renderer)
        connect(profilenode, renderer.profile)


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


# TODO: also support splitting BAM files
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


def split_to_file_wrapper(component):
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
    for node in component.output_nodes:
        if node.parallel and node.get_name() not in ["stdout", "stderr"]:
            comp = SplitToFile()
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
        else:
            wrapper.add(node)
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
                          **kwargs)
            return
        super().__init__(**kwargs)

        if star_aligner is not None and star_aligner:
            indexbuilder = StarIndexBuilder(num_targets=num_targets,
                                            total_target_length=total_target_length,
                                            nproc=nproc)
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
                 min_target_length=None,
                 total_target_length=None,
                 star_aligner=None,
                 min_mapq=None,
                 disable_soft_clipping=None,
                 nproc=None,
                 **kwargs):
        """
        Note: (does not perform bowtie index building, mutation parsing/counting,
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

        if U is not None:
            # unpaired reads
            if fastq_list:
                # add component to concatenate input files into single streams
                append = Appender(inputs=U)
                progmonitor = ProgressMonitor()
                qtrimmer = QualityTrimmer(min_qual=min_qual_to_trim,
                		                  window=window_to_trim,
                	                      min_length=min_length_to_trim)
                filter = LengthFilter()
                connect(append, progmonitor)
                self.add([append,
                          progmonitor,
                          qtrimmer,
                          filter])
            else:
                progmonitor = ProgressMonitor(input=U)
                qtrimmer = QualityTrimmer(min_qual=min_qual_to_trim,
                                          window=window_to_trim,
                	                      min_length=min_length_to_trim)
                filter = LengthFilter()
                self.add([progmonitor,
                          qtrimmer,
                          filter])
            connect(progmonitor, qtrimmer)
            connect(qtrimmer, filter)
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
            filter = LengthFilter()
            connect(qtrimmer1.trimmed, interleaver.R1)
            connect(qtrimmer2.trimmed, interleaver.R2)
            connect(interleaver, merger)
            connect(merger, filter)
            self.add([interleaver,
                      merger,
                      filter])

        if star_aligner is not None and star_aligner:
            aligner = StarAligner(reorder=preserve_order,
                                  assoc_rna=assoc_rna,
                                  min_target_length=min_target_length,
                                  disable_soft_clipping=disable_soft_clipping,
                                  nproc=nproc)
        else:
           aligner = BowtieAligner(reorder=preserve_order,
                                   assoc_rna=assoc_rna,
                                   disable_soft_clipping=disable_soft_clipping,
                                   nproc=nproc)
        connect(filter.filtered, aligner.fastq)
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
        self.add(aligner.index)
        self.add(aligner.aligned)

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
        output_classified: output classified mutations (for debugging)
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
                 target=None,
                 num_samples=None,
                 output_variant_counts=None,
                 output_mutation_counts=None,
                 min_mapq=None,
                 min_depth=None,
                 max_bg=None,
                 norm=None,
                 output_classified=None,
                 separate_ambig_counts=None,
                 right_align_ambig_dels=None,
                 right_align_ambig_ins=None,
                 min_mutation_separation=None,
                 min_qual_to_count=None,
                 mutation_type_to_count=None,
                 random_primer_len=None,
                 **kwargs):
        require_explicit_kwargs(locals())
        assert isinstance(num_samples, int)
        super().__init__(**kwargs)
        self.assoc_rna = target_name
        samples = ["Modified", "Untreated", "Denatured"]
        # TODO: add kwargs to connect input nodes here?
        counts = []
        ambig_counts = []

        for i in range(num_samples):
            sample = samples[i]
            parser = MutationParser(name="MutationParser_" + sample,
                                    min_mapq=min_mapq,
                                    assoc_sample=sample)
            counter = MutationCounter(name="MutationCounter_" + sample,
                                      target_length=target_length,
                                      variant_out=output_variant_counts,
                                      mutations_out=output_mutation_counts,
                                      classified_out=output_classified,
                                      assoc_sample=sample,
                                      separate_ambig_counts=separate_ambig_counts,
                                      right_align_ambig_dels=right_align_ambig_dels,
                                      right_align_ambig_ins=right_align_ambig_ins,
                                      min_mutation_separation=min_mutation_separation,
                                      min_qual=min_qual_to_count,
                                      random_primer_len=random_primer_len,
                                      mutation_type=mutation_type_to_count)
            connect_nodes(parser.parsed_mutations, counter.mut)
            self.add([parser,
                      counter])
            counts.append(counter.mutations)

        profilehandler = ProfileHandler(target=target,
                                        target_name=target_name,
                                        mindepth=min_depth,
                                        maxbg=max_bg,
                                        random_primer_len=random_primer_len,
                                        counts=counts,
                                        norm=norm)
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
                 min_target_length=None,
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
                 nproc=None,
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
                         nproc=nproc)

        sample = Sample(U=U,
                        R1=R1, R2=R2,
                        o1=True, o2=True,
                        min_target_length=min_target_length,
                        disable_soft_clipping=disable_soft_clipping,
                        min_mapq=min_mapq,
                        assoc_sample="sequence correction",
                        preserve_order=preserve_order,
                        star_aligner=star_aligner,
                        nproc=nproc,
                        min_qual_to_trim=min_qual_to_trim,
                        window_to_trim=window_to_trim,
                        min_length_to_trim=min_length_to_trim
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
                                        min_mapq=min_mapq)
                counter = MutationCounter(name="MutationCounter_{}".format(i+1),
                                          assoc_rna=target_names[i],
                                          target_length=target_lengths[i],
                                          variant_out=True,
                                          mutations_out=False,
                                          min_qual=min_qual_to_count,
                                          random_primer_len=random_primer_len)
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
            parser = MutationParser(min_mapq=min_mapq)
            counter = MutationCounter(variant_out=True,
                                      mutations_out=False,
                                      target_length=target_lengths[0],
                                      min_qual=min_qual_to_count,
                                      random_primer_len=random_primer_len)
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


class MutationRenderer(Component):
    """
    For debugging, make a postscript file with parsed mutations
    and classified mutations rendered on each read
    """
    def __init__(self,
                 min_mapq=None,
                 target_name=None,
                 **kwargs):
        self.min_mapq=min_mapq
        self.target_name=target_name
        super().__init__(**kwargs)
        self.add(InputNode(name="sam",
                           parallel=False))
        self.add(InputNode(name="parsed_mutations",
                           parallel=False))
        self.add(InputNode(name="classified_mutations",
                           parallel=False))
        self.add(OutputNode(name="ps",
                            extension="ps",
                            parallel=False))
        self.add(StdoutNode())
        self.add(StderrNode())


    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "render_mutations_ps.py"),
               "{sam}",
               '"{}"'.format(sanitize(self.target_name)),
               str(self.min_mapq),
               "{parsed_mutations}",
               "{classified_mutations}",
               "{ps}"]
        return cmd

    def after_run_message(self):
        return self.read_stdout()


class MutationCorrelator(Component):
    """
    Compute correlated mutation pairs
    """
    def __init__(self,
                 mut=None,
                 mut_bg=None,
                 length=None,
                 arbitrary_args=None,
                 **kwargs):
        self.length = length
        self.arbitrary_args = arbitrary_args
        super().__init__(**kwargs)
        self.add(InputNode(name="mut"))
        if mut is not None:
            connect(mut, self.mut)
        if mut_bg is not None:
            self.add(InputNode(name="mut_bg"))
            connect(mut_bg, self.mut_bg)
        self.add(OutputNode(name="correlated",
                            extension=".txt"))
        self.add(OutputNode(name="matrix",
                            parallel=False,
                            extension=''))
        self.add(StdoutNode())
        self.add(StderrNode())


    def cmd(self):
        cmd = [pyexe,
               os.path.join(bin_dir, "single_molecule/RINGMaPcorrelations.py"),
               "{mut}",
                "{correlated}"]
        cmd += ["--molsize", str(self.length)]
        cmd += ["--writematrixfile", "{matrix}"]
        if self.mut_bg is not None:
            cmd += ["--bgsubtract", "{mut_bg}"]
        if self.arbitrary_args is not None and len(self.arbitrary_args) > 0:
            cmd += self.arbitrary_args
        return cmd

    def after_run_message(self):
        return self.read_stdout()
