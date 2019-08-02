"""
Set up a single end-to-end workflow for sequence correction and/or MaP analysis.

Note: The current workflow implementation does not allow for changes to the
workflow topology after partial execution.

"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
from pprint import pprint

from pyshapemap.connect import *
from pyshapemap.components import *
from pyshapemap.pipeline import Pipeline


from pyshapemap.util import \
    require_explicit_kwargs, \
    read_fasta_names_lengths, \
    sanitize
from pyshapemap.flowchart import draw_flowchart

def build_pipeline(fastq=None,
                   out="shapemapper_out",
                   temp="shapemapper_temp",
                   structured_output=None,
                   overwrite=None,
                   name=None,
                   serial=None,
                   verbose=None,
                   target=None,
                   amplicon=None,
                   primers=None,
                   primers_in_sequence=None,
                   trim_primers = None,
                   require_forward_primer_mapped = None,
                   require_reverse_primer_mapped = None,
                   max_primer_offset = None,
                   target_raw=None,
                   nproc=None,
                   max_paired_fragment_length=None,
                   max_search_depth=None,
                   max_reseed=None,
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
                   render_mutations=None,
                   render_must_span=None,
                   max_pages=None,
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
                   genomeSAindexNbase=None,
                   star_shared_index=None,
                   rerun_on_star_segfault=None,
                   disable_soft_clipping=None,
                   render_flowchart=None,
                   per_read_histograms=None,
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
    primerlocator = None
    if amplicon:
        # FIXME: make clear on flowchart that fasta inputs are shared with AlignPrep
        primerlocator = PrimerLocator(fastas=target,
                                       target_names=target_names,
                                       primer_files=primers,
                                       primers_in_sequence=primers_in_sequence)
        pipeline.add(primerlocator)

    # boolean options - pass through to MutationParser
    pipeline.max_primer_offset = max_primer_offset
    pipeline.require_forward_primer_mapped = require_forward_primer_mapped
    pipeline.require_reverse_primer_mapped = require_reverse_primer_mapped
    pipeline.trim_primers = trim_primers

    pipeline.target_names = target_names
    pipeline.target_lengths = target_lengths
    pipeline.random_primer_len = random_primer_len
    pipeline.star_aligner = star_aligner
    pipeline.rerun_on_star_segfault = rerun_on_star_segfault

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
    if "correct_seq" in fastq:
        kw = {}
        if "U" in fastq["correct_seq"]:
            kw["U"] = fastq["correct_seq"]["U"]
        if "R1" in fastq["correct_seq"]:
            kw["R1"] = fastq["correct_seq"]["R1"]
            kw["R2"] = fastq["correct_seq"]["R2"]
        seqcorrector = CorrectSequence(target=target,
                                       target_names=target_names,
                                       target_lengths=target_lengths,
                                       total_target_length=sum(target_lengths),
                                       nproc=nproc,
                                       maxins=max_paired_fragment_length,
                                       max_search_depth=max_search_depth,
                                       max_reseed=max_reseed,
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
                                       genomeSAindexNbase=genomeSAindexNbase,
                                       star_shared_index=star_shared_index,
                                       amplicon=amplicon,
                                       max_primer_offset=max_primer_offset,
                                       require_forward_primer_mapped=require_forward_primer_mapped,
                                       require_reverse_primer_mapped=require_reverse_primer_mapped,
                                       trim_primers=trim_primers,
                                       **kw)
        pipeline.add(seqcorrector)


    num_samples = sum([1 for sample in ["modified","untreated","denatured"]
                       if sample in fastq])

    if num_samples > 0:
        # Preparation for main alignments
        if "correct_seq" in fastq:
            alignprep = AlignPrep(target=seqcorrector.corrected,
                                  num_targets=len(target_names),
                                  total_target_length=sum(target_lengths),
                                  star_aligner=star_aligner,
                                  genomeSAindexNbase=genomeSAindexNbase,
                                  nproc=nproc)
        else:
            alignprep = AlignPrep(target=target,
                                  num_targets=len(target_names),
                                  total_target_length=sum(target_lengths),
                                  star_aligner=star_aligner,
                                  genomeSAindexNbase=genomeSAindexNbase,
                                  nproc=nproc)
        pipeline.add(alignprep)

        if name != '':
            n = name + '_'
        else:
            n = ''

        out_prefix = os.path.join(out, n)

        input_is_unpaired = True
        for sample in ["modified", "untreated", "denatured"]:
            if (sample in fastq
                and "R1" in fastq[sample]
                and fastq[sample]["R1"] is not None):
                input_is_unpaired = False

        num_samples = 0
        # Main alignments
        # - also collect some output nodes to connect to next stages
        mapped_nodes = {}
        for sample in ["modified", "untreated", "denatured"]:
            if sample not in fastq:
                continue
            num_samples += 1

            kw = {}
            if "U" in fastq[sample]:
                kw["U"] = fastq[sample]["U"]
            if "R1" in fastq[sample]:
                kw["R1"] = fastq[sample]["R1"]
                kw["R2"] = fastq[sample]["R2"]
            if len(target_names) == 1:
                kw["assoc_rna"] = target_names[0]
            p = Sample(name=sample.capitalize(),
                       assoc_sample=sample.capitalize(),
                       nproc=nproc,
                       maxins=max_paired_fragment_length,
                       max_search_depth=max_search_depth,
                       max_reseed=max_reseed,
                       preserve_order=preserve_order,
                       total_target_length=sum(target_lengths),
                       star_aligner=star_aligner,
                       star_shared_index=star_shared_index,
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
        # sequence lengths, and also to read the number of primer pairs
        # associated with each target RNA
        updated_target_lengths = list(target_lengths)
        if "correct_seq" in fastq:
            for i in range(len(target_names)):
                updated_target_lengths[i] = seqcorrector["L{}".format(i+1)]

        num_primer_pairs = [-999]*len(target_names) # negative value indicates no amplicon filtering performed
        if require_forward_primer_mapped and require_reverse_primer_mapped:
            for i in range(len(target_names)):
                num_primer_pairs[i] = primerlocator["n_pairs_{}".format(i+1)]

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
                              primer_pairs=num_primer_pairs[i],
                              target=alignprep.target.input_node,
                              num_samples=num_samples,
                              output_variant_counts=False,
                              output_mutation_counts=True,
                              min_mapq=min_mapq,
                              maxins=max_paired_fragment_length,
                              input_is_unpaired=input_is_unpaired,
                              min_depth=min_depth,
                              max_bg=max_bg,
                              norm=indiv_norm,
                              separate_ambig_counts=separate_ambig_counts,
                              right_align_ambig_dels=right_align_ambig_dels,
                              right_align_ambig_ins=right_align_ambig_ins,
                              min_mutation_separation=min_mutation_separation,
                              mutation_type_to_count=mutation_type_to_count,
                              min_qual_to_count=min_qual_to_count,
                              random_primer_len=random_primer_len,
                              amplicon=amplicon,
                              max_primer_offset=max_primer_offset,
                              require_forward_primer_mapped=require_forward_primer_mapped,
                              require_reverse_primer_mapped=require_reverse_primer_mapped,
                              trim_primers=trim_primers,
                              render_mutations=render_mutations,
                              render_must_span=render_must_span,
                              max_pages=max_pages,
                              per_read_histograms=per_read_histograms,
                              )
            profile_nodes.append(p.ProfileHandler.CalcProfile.profile)
            # connect aligned reads nodes to post-alignment inputs
            # FIXME: gotta be a cleaner way to do this
            try:
                for sample in mapped_nodes:
                    from_node = mapped_nodes[sample][i]
                    to_node = p["MutationParser_"+sample.capitalize()].input
                    connect(from_node, to_node)
            except AttributeError:
                pass

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

    mp_comps = pipeline.collect_low_level_components(name="MutationParser*")
    mr_comps = pipeline.collect_low_level_components(name="MutationRendererPs")
    md_comps = pipeline.collect_low_level_components(name="RenderMappedDepths")
    rf_comps = pipeline.collect_low_level_components(name="RenderFigures")
    for i in range(len(target_names)):
        # connect amplicon primer pair location files to MutationParser,
        # RenderFigures, MutationRendererPs, and RenderMappedDepths components (if any)
        if primerlocator is not None:
            for i in range(len(target_names)):
                from_node = primerlocator["locs_{}".format(i + 1)]
                assoc_rna = from_node.assoc_rna
                for comp in mp_comps + mr_comps + rf_comps:
                    if comp.assoc_rna == assoc_rna:
                        connect(from_node, comp.primers)
                for comp in md_comps:
                    if comp.assoc_rna == assoc_rna:
                        connect(from_node, comp.primer_locations)



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
            prealn_components = []
            for c in pipeline.collect_low_level_components(name="*Aligner*"):
                if "CorrectSequence" not in c.get_parent_names():
                    for node in c.input_nodes:
                        if node.get_name() != "index":
                            prealn_components.append(node.input_node.parent_component)
            for c in list(set(prealn_components)):
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
                newc = split_to_file_wrapper(c, selected_out_names=["parsed_mutations"])
                p.replace(c, newc)
        #if "skip_flowchart" not in kwargs:
        #    # FIXME: this render doesn't show certain connections correctly, even if no code is
        #    #        executed between the previous draw_flowchart() call and here
        #    draw_flowchart(pipeline,
        #                   os.path.join(pipeline.temp, name+"_flowchart_post-wrapper.svg"),
        #                   path=pipeline.temp,
        #                   highlight_dir=pipeline.out)

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
                          render_mutations=render_mutations)

    return pipeline


def move_output_files(pipeline,
                      output_processed_reads=None,
                      output_aligned=None,
                      output_parsed=None,
                      output_counted=None,
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

    # simplified reactivity profiles for VARNA or Ribosketch
    for node in pipeline.collect_component_nodes(name="varna"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)
                                   +"_varna_colors.txt"))
    for node in pipeline.collect_component_nodes(name="ribosketch"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)
                                   +"_ribosketch_colors.txt"))

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

    # mapped depth PDFs
    for node in pipeline.collect_component_nodes(name="depth_fig"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)+
                                   "_mapped_depths.pdf"))

    # per primer-pair estimated abundances
    for node in pipeline.collect_component_nodes(name="est_abundances"):
        node.set_file(os.path.join(pipeline.out,
                                   pipeline.name+"_"+
                                   sanitize(node.assoc_rna)+
                                   "_per-amplicon_abundance.txt"))

    # whatever processed reads end up passed on to alignment
    if output_processed_reads:
        comps = [c for c in pipeline.collect_low_level_components(name="*Aligner*")
                 if "CorrectSequence" not in c.get_parent_names()]
        for i, comp in enumerate(comps):
            sample = comp.assoc_sample
            # for bowtie2, tab6
            # for STAR, R1+R2 on one component and fastq (unpaired) on another
            extension=".txt"
            node = None
            for n in ["tab6", "R1", "R2", "fastq"]:
                try:
                    node = getattr(comp, n)
                    extension = node.get_extension()
                except AttributeError:
                    continue

                if "Bowtie" in comp.get_name():
                    filename = os.path.join(pipeline.out,
                                            pipeline.name + '_' + sample + '_processed_reads.' + extension)
                else:
                    id = n
                    if n == "fastq":
                        id = "unpaired-or-merged"
                    else:
                        id = "merged-"+id
                    filename = os.path.join(pipeline.out,
                                            pipeline.name + '_' + sample + '_processed_reads_' + id + '.' + extension)

                if not pipeline.serial:
                    # locate the split to file node, not the in memory pipe
                    c = node.input_node.input_node.parent_component.to_file
                    c.set_file(filename)
                else:
                    node.input_node.input_node.set_file(filename)

    # aligned reads
    if output_aligned:
        comps = [c for c in pipeline.collect_low_level_components(name="*Aligner*")
                 if "CorrectSequence" not in c.get_parent_names()]
        for i, comp in enumerate(comps):
            sample = comp.assoc_sample
            extension = comp.aligned.get_extension()

            if "Bowtie" in comp.get_name():
                filename = os.path.join(pipeline.out,
                                        pipeline.name + '_' + sample + '_aligned.' + extension)
            else:
                paired = "paired" if "unpaired" not in comp.get_name() else "unpaired"
                filename = os.path.join(pipeline.out,
                                        pipeline.name+'_'+sample+'_aligned_'+paired+'.'+extension)

            if not pipeline.serial:
                # locate the split to file node, not the in memory pipe
                c = comp.aligned.output_nodes[0].output_nodes[0].parent_component.to_file
                c.set_file(filename)
            else:
                comp.aligned.set_file(filename)

    # parsed mutations
    if output_parsed:
        for comp in pipeline.collect_low_level_components(name="MutationParser*"):
            if "CorrectSequence" in comp.get_parent_names():
                continue
            sample = comp.assoc_sample
            rna = sanitize(comp.assoc_rna)
            filename = os.path.join(pipeline.out,
                                    pipeline.name+'_'+sample+'_'+rna+'_parsed.mut')
            if not pipeline.serial:
                # locate the split to file node, not the in memory pipe
                # - jeebus this is ugly. need a cleaner stream splitter interface
                c = comp.parsed_mutations.output_nodes[0].output_nodes[0].parent_component.to_file
                c.set_file(filename)
            else:
                comp.parsed_mutations.set_file(filename)

    # mutation/variant/depth counts
    if output_counted:
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
        for node in pipeline.collect_component_nodes(name="pdf"):
            comp = node.parent_component.parent_component
            if "MutationRenderer" not in comp.get_name():
                continue
            node.set_file(os.path.join(pipeline.out,
                                       pipeline.name+'_'+comp.assoc_sample+
                                       '_'+sanitize(comp.assoc_rna)+'_'+
                                       'debug_mutations.pdf'))

