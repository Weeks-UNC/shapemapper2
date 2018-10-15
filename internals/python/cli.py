"""
High-level shapemapper commandline interface

"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys
import os
import shutil

import pyshapemap.pipeline_arg_parser as ap
from pyshapemap.pipeline_builder import build_pipeline
from pyshapemap.util import Logger, timestamp, version
from pyshapemap.flowchart import draw_flowchart

def run(args):
    assert isinstance(args, list)
    try:
        ap.check_help_version(args[1:]) # print version or usage and exit if those args provided
        ap.check_no_whitespace(args[1:]) # check for disallowed whitespace characters in input arguments
        logpath, rest_args = ap.get_log_path(args[1:])

        log = Logger(filepath=logpath,
                     stream=sys.stdout)
        ts = timestamp()
        term_width, term_height = shutil.get_terminal_size()
        print("\n" + '#' * term_width)
        print("\nStarted ShapeMapper v{} at {}\nOutput will be logged to {}".format(version(),
                                                                                    ts,
                                                                                    logpath))
        log.file.write("\n" + '#' * term_width)
        log.file.write("\nStarted ShapeMapper v{} at {}\n".format(version(),
                                                                  ts))
        # override stdout and stderr globally to redirect through logger
        sys.stdout=log
        sys.stderr=log

        print("Running from directory: "+os.getcwd())

        msg = "args: "
        for arg in args[1:]:
            if len(arg) > 100:
                msg += ' '+arg[:100]+'...'
            else:
                msg += ' '+arg
        print(msg)

        pipeline, arg_dict = ap.construct(rest_args)
        print("Created pipeline at {}".format(timestamp()))

        if pipeline.render_flowchart:
            p = pipeline.flowchart_path
            s = draw_flowchart(pipeline,
                               p,
                               path=pipeline.temp,
                               highlight_dir=pipeline.out)
            if s:
                print("Rendered pipeline flowchart to {} at {}".format(p,
                                                                       timestamp()))
            else:
                msg = "Attempted to render pipeline flowchart to {} at {}, "
                msg += "but graphviz gave a warning or error."
                msg = msg.format(p,
                                 timestamp())
                print(msg)

        success = pipeline.run()

        if not success:
            if arg_dict["rerun_on_star_segfault"]:
                if pipeline.any_star_segfault():
                    # rerun with STAR --genomeSAindexNbase set to defined value (default 3)
                    g = arg_dict["rerun_genomeSAindexNbase"]
                    print("\nSTAR segfault detected. Rerunning pipeline with --genomeSAindexNbase {} and --overwrite.\n".format(g))
                    term_width, term_height = shutil.get_terminal_size()
                    print("\n" + '-' * term_width)
                    new_args = arg_dict
                    new_args["genomeSAindexNbase"] = g
                    new_args["overwrite"] = True
                    pipeline2 = build_pipeline(**new_args)
                    print("Created pipeline at {}".format(timestamp()))
                    success = pipeline2.run()
            if not success:
                print("ShapeMapper run failed at {}.".format(timestamp()))
                sys.exit(1)

        if success:
            failed_targets = pipeline.get_failed_targets()
            if len(failed_targets) == 0:
                print("ShapeMapper run successfully completed at {}.".format(timestamp()))
            else:
                print("ShapeMapper run completed at {}.".format(timestamp()))
                if len(failed_targets) == 1:
                    msg = "WARNING: This RNA has a possible poor quality reactivity profile: "+failed_targets[0]
                    print(msg)
                else:
                    msg = "WARNING: These RNAs have possible poor quality reactivity profiles: "
                    msg += ', '.join(failed_targets)
                    print(msg)
                print("See quality control checks above for details.")

            if min(pipeline.target_lengths) < 1000 and not pipeline.trim_primers:
                msg = "WARNING: amplicon primer trimming and read mapping "
                msg += "filtering is not enabled. To enable, set primer "
                msg += "sequence to lowercase and add the --amplicon option, "
                msg += "or provide a file with primer pairs with the --primers "
                msg += "option (see docs/primer_filtering.md)."
                print(msg)

            if max(pipeline.target_lengths) > 2000 and not pipeline.star_aligner:
                msg = "WARNING: Bowtie2 is slower than STAR for long sequences."
                msg += " Consider using STAR with the --star-aligner option."
                print(msg)

            # FIXME: don't show this warning if tiled amplicon primer pairs used
            if max(pipeline.target_lengths) > 800 and pipeline.random_primer_len==0:
                msg = "WARNING: no random primer length was specified, "
                msg += "but at least one RNA is longer than a typical "
                msg += "directed-primer amplicon. If analyzing a randomly "
                msg += "primed experiment that was not subjected to a "
                msg += "Nextera prep, use the --random-primer-len option "
                msg += "to exclude mutations within primer binding regions."
                print(msg)

            sys.exit(0)

    finally:
        try:
            pipeline.clean_up()
        except (NameError, AttributeError):
            pass
        try:
            pipeline2.clean_up()
        except (NameError, AttributeError):
            pass

if __name__ == "__main__":
    run(sys.argv)

