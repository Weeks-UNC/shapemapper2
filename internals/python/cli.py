"""
High-level shapemapper commandline interface

"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2018 Steven Busan.                     #
# --------------------------------------------------------------------- #

import sys
import re
import os
import shutil

import pyshapemap.pipeline_arg_parser as ap
from pyshapemap.pipeline_builder import build_pipeline
from pyshapemap.util import Logger, timestamp, version
from pyshapemap.flowchart import draw_flowchart
from pyshapemap.nodes import FileNode


def file_check(args):
   """
      Checks if all .fastq/.fastq.gz's and .fa files
      actually exist.
      
      If they do not exist, this function notifies 
      user and halts execution.
   """

   files = []
   for arg in args:
      if arg.find(".gz") != -1 or arg.find(".fastq") != -1 or arg.find(".fa") != -1:
         files.append(arg)

   for f in files:
      if not os.path.exists(f):
         print("Error, {} not found. Please verify name of file and resubmit.".format(f))
         sys.exit()

    

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


        #If N7 and parallel makes additional FIFO pipes so script may be run in parallel. 
        #If dms, makes a fifo in accordance with new dms processing naming scheme.
        if "--N7" in args and "--serial" not in args:
           success = pipeline.run(N7=True, dms=True)
        elif "--dms" in args and "--serial" not in args:
           success = pipeline.run(dms=True)
        else:
           success = pipeline.run()

        # Fix to aid downstream file moving stuff
        if "--star-aligner" in args and "--serial" not in args:
            args.append("--serial")

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

            #Transfers .mut* and .sam files if specified during a parallel run. Necessary due to changes
            #required by N7 parallel processing
            if ('--output-aligned-reads' in args or '--output-aligned' in args) and '--serial' not in args:
               samFiles = []
               for c in pipeline.collect_low_level_components():
                  for node in c.get_component_nodes():

                     if "SplitToFile" in c.get_name() and "--dms" not in args:
                        for out_node in node.output_nodes:
                           if isinstance(out_node, FileNode):
                              fname = out_node.filename
                              if "to_file.sam" in fname:
                                 samFiles.append(fname)


                     if isinstance(node.input_node, FileNode) and "--dms" in args:
                        fname = node.input_node.filename

                        if '--star-aligner' in args:
                           if ".sam" in fname and "SamMixer" in fname:
                              if fname not in samFiles:
                                 samFiles.append(fname)

                        else:
                           if ".sam" in fname:
                              if fname not in samFiles:
                                 samFiles.append(fname)


               for sam in samFiles:
                  splFile = sam.split("/")
                  truncFile = splFile[-1]
                  out_name = pipeline.out + "/" + truncFile
                  if "--dms" in args:
                     shutil.move(sam + "move", out_name)
                  else:
                     shutil.move(sam , out_name)

            # Transfer processed reads at end of run
            if "--output-processed-reads" in args and "--star-aligner" in args and "--serial" not in args:
                if pipeline.temp == "shapemapper_temp":
                    t = pipeline.temp + f"/{arg_dict['name']}"
                else:
                    t = pipeline.temp
                extracted_files = []
                for root, _, files in os.walk(pipeline.temp):
                    for f in files:
                        if f.endswith(".processed_read_move_extension"):
                            extracted_files.append(root + "/" + f)
                for f in extracted_files:
                    os.rename(f, pipeline.out + "/" + f.removesuffix(".processed_read_move_extension").split("/")[-1])
                    


            #If N7 QC is triggered, removes all .mutga and .txtga files so user must 
            #at least read warning (in order to obtain flag) before gaining access
            #to these datastreams for subsequent processing. (Like RINGMAP or PAIRMAP)
            #print("Removing files from {}".format(pipeline.out))
            sub_name=arg_dict['name']
            try:
               files = os.listdir(pipeline.out)
               n7_message_file_regex =  "\." + sub_name + ".*_n7message.txt" 
               for f in files:
                  if re.findall(n7_message_file_regex, f) != []:
                     n7_message_file = pipeline.out + "/" + f
               with open(n7_message_file, "r") as n7_file:
                  warnings = [line for line in n7_file][0].split(",")
                  if "low N7" in warnings:
                     to_remove = []

                     mod_ga_regex = sub_name + "_Modified.*\.mutga"
                     unt_ga_regex = sub_name + "_Untreated.*\.mutga"
                     txt_ga_regex =   "^[^.]*" + sub_name + ".*" + "profile\.txtga"
                     for f in files:
                        #if f.find(".mutga") != -1:
                        if re.findall(mod_ga_regex, f) != []:
                           to_remove.append(f)
                        elif re.findall(unt_ga_regex, f) != []:
                           to_remove.append(f)
                        elif re.findall(txt_ga_regex, f) != []:
                           to_remove.append(f)
                     for r in to_remove:
                        os.remove(pipeline.out + "/" + r)
               os.remove(n7_message_file)

            except UnboundLocalError as e:
                if str(e).find("n7_message_file") != -1: #If n7_message_file not bound, it simply wasn't generated and we can skip.
                   pass
                else:
                   raise # Else re raise error as it is unexpected.

            if '--output-temp' not in args:
               print("Deleting temp files. Use --output-temp to retain temp files.")
               shutil.rmtree(pipeline.temp + "/" + sub_name)

               if os.path.exists(pipeline.temp): # Workaround to avoid error message if run is un named
                   subdirs = os.listdir(pipeline.temp)
                   if len(subdirs) == 0:
                       shutil.rmtree(pipeline.temp)

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
    file_check(sys.argv)
    run(sys.argv)

