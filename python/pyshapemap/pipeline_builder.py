"""
Set up a single end-to-end workflow for sequence correction and/or MaP analysis.

Note: The current workflow implementation does not allow for changes to the
workflow topology after partial execution.

"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os

from pyshapemap.component import *
from pyshapemap.components import *
from pyshapemap.pipeline import Pipeline

from pyshapemap.util import require_explicit_kwargs, read_fasta_names_lengths, sanitize
from pyshapemap.flowchart import draw_flowchart

def build_pipeline(fastqs,
                   out="shapemapper_out",
                   temp="shapemapper_temp",
                   structured_output=None,
                   overwrite=None,
                   name=None,
                   serial=None,
                   verbose=None,
                   target=None,
                   target_raw=None,
                   nproc=None,
                   min_depth=None,
                   max_bg=None,
                   preserve_order=None,
                   min_mapq=None,
                   min_seq_depth=None,
                   min_freq=None,
                   indiv_norm=None,
                   output_processed_reads=None,
                   output_aligned=None,
                   output_parsed=None,
                   output_counted=None,
                   output_classified=None,
                   render_mutations=None,
                   separate_ambig_counts=None,
                   right_align_ambig_dels=None,
                   right_align_ambig_ins=None,
                   mutation_type_to_count=None,
                   min_mutation_separation=None,
                   min_qual_to_count=None,
                   random_primer_len=None,
                   min_qual_to_trim=None,
                   window_to_trim=None,
                   min_length_to_trim=None,
                   star_aligner=None,
                   disable_soft_clipping=None,
                   calc_correlations=None,
                   corr_args=None,
                   render_flowchart=None,
                   **kwargs):
    require_explicit_kwargs(locals())

    if name is not None and len(name)>0:
        name = sanitize(name)
    else:
        name = "Pipeline"

    pipeline = Pipeline(name=name)
    pipeline.overwrite = overwrite
    pipeline.render_flowchart = render_flowchart

    # write raw fasta sequences to file. bit of a hack to allow
    # easily passing in short sequences from JS frontend without
    # shell access
    if target_raw != "":
        os.makedirs(temp, exist_ok=True)
        fa_name = temp+"/"+name+"_raw_target.fa"
        open(fa_name, "w").write(target_raw.replace('\\n','\n'))
        target.append(fa_name)

    target_names, target_lengths = read_fasta_names_lengths(target)
    pipeline.target_names = target_names
    pipeline.target_lengths = target_lengths
    pipeline.random_primer_len = random_primer_len
    pipeline.star_aligner = star_aligner

    # suggest STAR aligner for long sequences
    if max(target_lengths) > 2000 and not pipeline.star_aligner:
        msg = "Warning: Bowtie2 is slower than STAR for long sequences."
        msg += " Consider using STAR with the --star-aligner option."
        print(msg)

    # warn if no random primer length specified
    # (will also repeat this warning at end of run)
    if max(target_lengths) > 800 and random_primer_len==0:
        msg = "Warning: no random primer length was specified, "
        msg += "but at least one RNA is longer than a typical "
        msg += "directed-primer amplicon. Use --random-primer-len "
        msg += "to exclude mutations within primer binding regions."
        print(msg)

    # Reference sequence correction
    if "correct_seq" in fastqs:
        kw = {}
        if "U" in fastqs["correct_seq"]:
            kw["U"] = fastqs["correct_seq"]["U"]
        if "R1" in fastqs["correct_seq"]:
            kw["R1"] = fastqs["correct_seq"]["R1"]
            kw["R2"] = fastqs["correct_seq"]["R2"]
        seqcorrector = CorrectSequence(target=target,
                                       target_names=target_names,
                                       target_lengths=target_lengths,
                                       min_target_length=min(target_lengths),
                                       total_target_length=sum(target_lengths),
                                       nproc=nproc,
                                       min_qual_to_trim=min_qual_to_trim,
                                       window_to_trim=window_to_trim,
                                       min_length_to_trim=min_length_to_trim,
                                       disable_soft_clipping=False,
                                       min_mapq=10, # FIXME: put this default value somewhere more visible
                                       min_qual_to_count=min_qual_to_count,
                                       min_seq_depth=min_seq_depth,
                                       min_freq=min_freq,
                                       random_primer_len=random_primer_len,
                                       preserve_order=preserve_order,
                                       star_aligner=star_aligner,
                                       **kw)
        pipeline.add(seqcorrector)

    num_samples = sum([1 for sample in ["modified","untreated","denatured"]
                       if sample in fastqs])

    if num_samples > 0:
        # Preparation for main alignments
        if "correct_seq" in fastqs:
            alignprep = AlignPrep(target=seqcorrector.corrected,
                                  num_targets=len(target_names),
                                  total_target_length=sum(target_lengths),
                                  star_aligner=star_aligner,
                                  nproc=nproc)
        else:
            alignprep = AlignPrep(target=target,
                                  num_targets=len(target_names),
                                  total_target_length=sum(target_lengths),
                                  star_aligner=star_aligner,
                                  nproc=nproc)
        pipeline.add(alignprep)

        if name != '':
            n = name + '_'
        else:
            n = ''

        out_prefix = os.path.join(out, n)

        num_samples = 0
        # Main alignments
        # - also collect some output nodes to connect to next stages
        mapped_nodes = {}
        for sample in ["modified", "untreated", "denatured"]:
            if sample not in fastqs:
                continue
            num_samples += 1

            kw = {}
            if "U" in fastqs[sample]:
                kw["U"] = fastqs[sample]["U"]
            if "R1" in fastqs[sample]:
                kw["R1"] = fastqs[sample]["R1"]
                kw["R2"] = fastqs[sample]["R2"]
            if len(target_names) == 1:
                kw["assoc_rna"] = target_names[0]
            p = Sample(name=sample.capitalize(),
                       assoc_sample=sample.capitalize(),
                       nproc=nproc,
                       preserve_order=preserve_order,
                       min_target_length=min(target_lengths),
                       total_target_length=sum(target_lengths),
                       star_aligner=star_aligner,
                       disable_soft_clipping=disable_soft_clipping,
                       min_qual_to_trim=min_qual_to_trim,
                       window_to_trim=window_to_trim,
                       min_length_to_trim=min_length_to_trim,
                       **kw)

            connect(alignprep.index,
                    p.index)
            pipeline.add(p)

            # split aligned output by mapped RNA (if needed)
            if len(target_names)>1:
                splitter = SplitByTarget(name="SplitByTarget_"+sample.capitalize(),
                                         target_names=target_names)
                connect(p.aligned, splitter.input)
                pipeline.add(splitter)
                mapped_nodes[sample] = [node for node in splitter.output_nodes if not isinstance(node, (StdoutNode, StderrNode))]
            else:
                mapped_nodes[sample] = [p.aligned]

        # allow MutationCounter components to adjust to dynamic changes in 
        # sequence lengths
        updated_target_lengths = list(target_lengths)
        if "correct_seq" in fastqs:
            for i in range(len(target_names)):
                updated_target_lengths[i] = seqcorrector["L{}".format(i+1)]

        # Post-alignment steps, for each RNA target
        # - Mutation parsing, counting, and profile creation
        # - also collect profile output nodes for normalization (if normalizing as a group)
        if len(target_names) == 1:
            indiv_norm = True
        profile_nodes = []
        ambig_profile_nodes = []
        for i in range(len(target_names)):
            p = PostAlignment(name="PostAlignment_RNA_{}".format(i+1),
                              target_name=target_names[i],
                              target_length=updated_target_lengths[i],
                              target=alignprep.target.input_node,
                              num_samples=num_samples,
                              output_variant_counts=False,
                              output_mutation_counts=True,
                              min_mapq=min_mapq,
                              min_depth=min_depth,
                              max_bg=max_bg,
                              norm=indiv_norm,
                              output_classified=output_classified,
                              separate_ambig_counts=separate_ambig_counts,
                              right_align_ambig_dels=right_align_ambig_dels,
                              right_align_ambig_ins=right_align_ambig_ins,
                              min_mutation_separation=min_mutation_separation,
                              mutation_type_to_count=mutation_type_to_count,
                              min_qual_to_count=min_qual_to_count,
                              random_primer_len=random_primer_len)
            profile_nodes.append(p.ProfileHandler.CalcProfile.profile)
            # connect aligned reads nodes to post-alignment inputs
            # FIXME: gotta be a cleaner way to do this
            for sample in mapped_nodes:
                from_node = mapped_nodes[sample][i]
                to_node = p["MutationParser_"+sample.capitalize()].input
                connect(from_node, to_node)

            pipeline.add(p)

        # Normalize profiles as a group
        if not indiv_norm:
            normer = NormProfile(name="NormalizeAsGroup",
                                 profiles=profile_nodes,
                                 target_names=target_names)
            pipeline.add(normer)
            # rewire so normalized profiles get used for downstream steps
            for i in range(len(target_names)):
                node = profile_nodes[i]
                # need a list copy, since output_nodes will change during iteration
                connected_nodes = list(node.output_nodes)
                for cnode in connected_nodes:
                    # skip newly added connections
                    if cnode.parent_component is not normer:
                        disconnect(node, cnode)
                        connect(normer["normed_{}".format(i+1)],
                                cnode)


    pipeline.out = out
    if structured_output:
        pipeline.temp = out
    else:
        pipeline.temp = temp
    pipeline.structured_output = structured_output
    pipeline.serial = serial
    pipeline.verbose = verbose

    pipeline.flowchart_path = os.path.join(pipeline.out,
                                           name+"_flowchart.svg")


    # add wrapper components to split selected intermediate streams to files
    # TODO: could move some of this repeated code to function
    if ( (output_aligned or
          output_parsed or
          output_processed_reads)
         and not serial ):
        #if "skip_flowchart" not in kwargs:
        #    draw_flowchart(pipeline,
        #                   os.path.join(pipeline.temp, name+"_flowchart_pre-wrapper.svg"),
        #                   path=pipeline.temp,
        #                   highlight_dir=pipeline.out)
        if output_processed_reads:
            prealn_components = [c.fastq.input_node.parent_component
                                 for c in pipeline.collect_low_level_components(name="*Aligner*")
                                  if "CorrectSequence" not in c.get_parent_names()]
            for c in prealn_components:
                p = c.parent_component
                newc = split_to_file_wrapper(c)
                p.replace(c, newc)
        if output_aligned:
            aligner_components = [c for c in pipeline.collect_low_level_components(name="*Aligner*")
                                  if "CorrectSequence" not in c.get_parent_names()]
            for c in aligner_components:
                p = c.parent_component
                newc = split_to_file_wrapper(c)
                p.replace(c, newc)
        if output_parsed:
            parser_components = [c for c in pipeline.collect_low_level_components(name="MutationParser*")
                                 if "CorrectSequence" not in c.get_parent_names()]
            for c in parser_components:
                p = c.parent_component
                newc = split_to_file_wrapper(c)
                p.replace(c, newc)
        if output_classified:
           counter_components = [c for c in pipeline.collect_low_level_components(name="MutationCounter*")
                                 if "CorrectSequence" not in c.get_parent_names()]
           for c in counter_components:
               p = c.parent_component
               newc = split_to_file_wrapper(c)
               p.replace(c, newc)
        #if "skip_flowchart" not in kwargs:
        #    # FIXME: this render doesn't show certain connections correctly, even if no code is
        #    #        executed between the previous draw_flowchart() call and here
        #    draw_flowchart(pipeline,
        #                   os.path.join(pipeline.temp, name+"_flowchart_post-wrapper.svg"),
        #                   path=pipeline.temp,
        #                   highlight_dir=pipeline.out)

    # add components to render debug mutation info to postscript files
    # FIXME: skip sequence correction components
    if render_mutations:
        if not serial:
            comps = pipeline.collect_low_level_components(name="MutationParser*")
            parsed_nodes = []
            for c in comps:
                n = c.parent_component.SplitToFile.to_file
                parsed_nodes.append(n)
            comps = pipeline.collect_low_level_components(name="*Aligner*")
            sam_nodes = []
            for c in comps:
                n = c.parent_component.SplitToFile.to_file
                sam_nodes.append(n)
            classified_nodes = []
            comps = pipeline.collect_low_level_components(name="MutationCounter*")
            for c in comps:
                # use this if node is parallel
                n = c.parent_component.SplitToFile.to_file
                # use this if node is serial (need to rework parallel wrapper syntax)
                #n = c.classified_mutations
                classified_nodes.append(n)
        else:
            parsed_nodes = pipeline.collect_component_nodes(name="parsed_mutations")
            classified_nodes = pipeline.collect_component_nodes(name="classified_mutations")
            sam_nodes = pipeline.collect_component_nodes(name="aligned")
        # match up nodes with same RNA and sample
        for parsed_node in parsed_nodes:
            p = parsed_node.parent_component
            if p.assoc_sample == "sequence correction":
                continue

            sam_node = None
            for node in sam_nodes:
                p2 = node.parent_component
                if p.assoc_sample == p2.assoc_sample:
                    sam_node = node
                    break
            classified_node = None
            for node in classified_nodes:
                p2 = node.parent_component
                if ( p.assoc_rna == p2.assoc_rna and
                     p.assoc_sample == p2.assoc_sample ):
                    classified_node = node
                    break
            if sam_node is None or classified_node is None:
                raise RuntimeError("Error: unable to create mutation renderer component")
            comp_name = "MutationRenderer_{}_{}".format(
                p.assoc_sample,
                sanitize(p.assoc_rna)
            )
            renderer = MutationRenderer(name=comp_name,
                                        target_name=sanitize(p.assoc_rna),
                                        min_mapq=min_mapq,
                                        assoc_rna=p.assoc_rna,
                                        assoc_sample=p.assoc_sample)
            connect(sam_node, renderer.sam)
            connect(parsed_node, renderer.parsed_mutations)
            connect(classified_node, renderer.classified_mutations)
            pipeline.add(renderer)

    # component to perform single-molecule mutation correlation calcs
    if calc_correlations:
        for i in range(len(target_names)):
            rna = target_names[i]
            length = target_lengths[i]
            nodes = pipeline.collect_component_nodes(name="classified_mutations",
                                                     parent_name="MutationCounter*")
            # find output nodes with no existing connections
            nodes = [n for n in nodes if isinstance(n, OutputNode) and len(n.output_nodes) == 0]
            mut = None
            mut_bg = None
            for node in nodes:
                p = node.parent_component
                if p.assoc_sample == "Modified" and p.assoc_rna == rna:
                    mut = node
                elif p.assoc_sample == "Untreated" and p.assoc_rna == rna:
                    mut_bg = node
            if mut is None:
                raise RuntimeError("Error: unable to create mutation correlator component")
            comp_name = "MutationCorrelator_{}".format(sanitize(rna))
            correlator = MutationCorrelator(name=comp_name,
                                            mut=mut,
                                            mut_bg=mut_bg,
                                            length=length,
                                            arbitrary_args=corr_args,
                                            assoc_rna=rna)
            pipeline.add(correlator)

    #if "skip_flowchart" not in kwargs:
    #    draw_flowchart(pipeline,
    #                   os.path.join(pipeline.temp, name+"_flowchart_pre-setup.svg"),
    #                   path=pipeline.temp,
    #                   highlight_dir=pipeline.out)

    # option to skip setup is used by some testing modules
    if "skip_setup" not in kwargs:
        pipeline.setup()

    if not structured_output:
        move_output_files(pipeline,
                          output_processed_reads=output_processed_reads,
                          output_aligned=output_aligned,
                          output_parsed=output_parsed,
                          output_counted=output_counted,
                          output_classified=output_classified,
                          calc_correlations=calc_correlations,
                          render_mutations=render_mutations)

    return pipeline


def move_output_files(pipeline,
                      output_processed_reads=None,
                      output_aligned=None,
                      output_parsed=None,
                      output_counted=None,
                      output_classified=None,
                      calc_correlations=None,
                      render_mutations=None):
    """
    Change paths and simplify filenames for main output files
    """
    # A bit ugly doing this post-hoc, but the alternative (setting filenames
    # in component constructors) is really cluttered

    # corrected sequence FASTA
    try:
        pipeline.CorrectSequence.corrected.set_file(os.path.join(pipeline.out,
                                                    pipeline.name+"_corrected.fa"))
    except AttributeError:
        pass

    # tab-delimited profiles
    for node in pipeline.collect_component_nodes(name="normed*"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)
                                   +"_profile.txt"))

    # SHAPE files
    for node in pipeline.collect_component_nodes(name="shape"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)
                                   +".shape"))

    # MAP files
    for node in pipeline.collect_component_nodes(name="map"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)
                                   +".map"))

    # profile PDFs
    for node in pipeline.collect_component_nodes(name="profiles_fig"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)
                                   +"_profiles.pdf"))

    # histogram PDFs
    for node in pipeline.collect_component_nodes(name="histograms_fig"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)
                                   +"_histograms.pdf"))

    # whatever processed reads end up passed on to alignment
    if output_processed_reads:
        aligner_components = [c for c in pipeline.collect_low_level_components(name="*Aligner*")
                              if "CorrectSequence" not in c.get_parent_names()]
        for comp in aligner_components:
            sample = comp.assoc_sample
            extension = comp.fastq.get_extension()
            filename = os.path.join(pipeline.out,
                                    pipeline.name+'_'+sample+'_processed_reads.'+extension)
            if not pipeline.serial:
                # locate the split to file node, not the in memory pipe
                c = comp.fastq.input_node.input_node.parent_component.to_file
                c.set_file(filename)
            else:
                comp.fastq.input_node.set_file(filename)

    # aligned reads
    if output_aligned:
        for comp in pipeline.collect_low_level_components(name="*Aligner*"):
            if "CorrectSequence" in comp.get_parent_names():
                continue
            sample = comp.assoc_sample
            extension = comp.aligned.get_extension()
            filename = os.path.join(pipeline.out,
                                    pipeline.name+'_'+sample+'_aligned.'+extension)
            if not pipeline.serial:
                # locate the split to file node, not the in memory pipe
                c = comp.aligned.output_nodes[0].output_nodes[0].parent_component.to_file
                c.set_file(filename)
            else:
                # FIXME: this fails if called with --output-aligned and --serial
                comp.aligned.set_file(filename)

    # parsed mutations
    if output_parsed:
        for comp in pipeline.collect_low_level_components(name="MutationCounter*"):
            if "CorrectSequence" in comp.get_parent_names():
                continue
            sample = comp.assoc_sample
            rna = sanitize(comp.assoc_rna)
            filename = os.path.join(pipeline.out,
                                    pipeline.name+'_'+sample+'_'+rna+'_parsed.mut')
            if not pipeline.serial:
                # locate the split to file node, not the in memory pipe
                c = comp.mut.input_node.input_node.parent_component.to_file
                c.set_file(filename)
            else:
                comp.mut.set_file(filename)

    # mutation/variant/depth counts
    if output_counted or output_classified:
        for comp in pipeline.collect_low_level_components(name="MutationCounter*"):
            if "CorrectSequence" in comp.get_parent_names():
                continue
            sample = comp.assoc_sample
            rna = sanitize(comp.assoc_rna)
            if output_counted:
                comp.mutations.set_file(os.path.join(pipeline.out,
                                                     pipeline.name+'_'+sample+'_'+rna+'_mutation_counts.txt'))
                try:
                    comp.ambig_mutations.set_file(os.path.join(pipeline.out,
                                                               pipeline.name+'_'+sample+'_'+rna+
                                                               '_ambig_mutation_counts.txt'))
                except AttributeError:
                    pass
            if output_classified:
                filename = os.path.join(pipeline.out,
                                        pipeline.name+'_'+sample+'_'+rna+
                                        '_classified_mutations.txt')
                try:
                    # FIXME: now that the classified_mutations node is parallel by default, need
                    #        to rethink this.

                    #if not pipeline.serial:
                    #    # locate the split to file node, not the in memory pipe
                    #    # FIXME: this is pretty ugly
                    #    c = comp.classified_mutations.output_nodes[0].output_nodes[0].parent_component.to_file
                    #    c.set_file(filename)
                    #else:

                    if True:
                        comp.classified_mutations.set_file(filename)
                except AttributeError:
                    pass

    # correlated mutation pairs and associated matrices
    if calc_correlations:
        for node in pipeline.collect_component_nodes(name="correlated"):
            node.set_file(os.path.join(pipeline.out,
                                       pipeline.name+'_'+sanitize(node.parent_component.assoc_rna)+'_'+
                                       'correlated_mutations.txt'))
        for node in pipeline.collect_component_nodes(name="matrix"):
            node.set_file(os.path.join(pipeline.out,
                                       pipeline.name+'_'+sanitize(node.parent_component.assoc_rna)))

    # rendered mutation debug info
    if render_mutations:
        for node in pipeline.collect_component_nodes(name="ps"):
            node.set_file(os.path.join(pipeline.out,
                                       pipeline.name+'_'+node.parent_component.assoc_sample+
                                       '_'+sanitize(node.parent_component.assoc_rna)+'_'+
                                       'mutations.ps'))

