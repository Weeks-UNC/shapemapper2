"""
Construct a complete pipeline and map endpoints (PipelineModules)
Then construct new pipelines with error injectors and test for
module failure detection at each endpoint.

"""
# --------------------------------------------------------------------- #
#  This file is a part of ShapeMapper, and is licensed under the terms  #
#  of the MIT license. Copyright 2017 Steven Busan.                     #
# --------------------------------------------------------------------- #

import os
import sys
import shutil
import traceback

# so we can import from python files in pyshapemap
# folder, even though this script is a bit isolated
this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(this_dir, '../python'))

import pyshapemap.pipeline_arg_parser as ap
from pyshapemap.util import Logger, timestamp

from pyshapemap.pipeline import run_fail_test
#from pyshapemap.pipeline_templates import run_fail_test
#from pyshapemap.pipeline_base import map_pipeline_tree, get_module


# FIXME: suppress misleading BrokenPipeError warnings

if __name__ == "__main__":

    try:
        args = sys.argv[1:]
        out_folder, _ = ap.get_out_path(args)
        temp_folder, _ = ap.get_temp_path(args)

        reference_pipeline = ap.construct(args,
                                          skip_flowchart=True,
                                          skip_setup=True)

        names, _, component_locations =reference_pipeline.map_pipeline_tree()

        tests = []
        for n in range(len(component_locations)):
            # skip some components that don't do any format checking
            skip = False
            name = names[n]
            for s in ["Interleaver",
                      "Append",
                      "ProgressMonitor"]:
                if s in name:
                    skip = True
                    break
            if not skip:
                tests.append(n)

        passed = [False for n in tests]
        fail_msgs = []

        #print("[==========]")
        print("[----------] {} tests for component failure detection".format(len(tests)))

        k = -1
        for n in tests:
            k += 1
            p = component_locations[n]

            #term_width, term_height = shutil.get_terminal_size()
            #pad = ''.join(['-' for x in range(term_width)])+'\n'

            name = names[n]

            pass_flag = False
            fail_msg = "ERROR: unable to construct pipeline."
            pipeline = None
            try:
                pipeline = ap.construct(args,
                                        skip_flowchart=True,
                                        skip_setup=True)
                pipeline.name = pipeline.get_name()+"_"+p.replace('.','-')
                pipeline.flowchart_path = os.path.join(pipeline.out,
                                                       pipeline.get_name()+"_flowchart.svg")

                if pipeline is not None:
                    fail_msg = "ERROR: unable to run pipeline fail test."
                    pass_flag, success, err = run_fail_test(pipeline=pipeline,
                                                            failing_module_loc=p,
                                                            failing_module_name=name)

                    if pass_flag:
                        print("[       OK ] Fail test for component at {} ({})".format(p,
                                                                                       name))
                    else:
                        print("Expected failure at component \"{}\"".format(name))
                        if success:
                            print("Instead got no failure.")
                        else:
                            print("Instead got failure:\n"+err)
                        fail_msg = "[  FAILED  ] Fail test for component at {} ({})".format(p,
                                                                                            name)
                        print(fail_msg)
                    sys.stdout.flush()

            except Exception as e:
                if isinstance(e, KeyboardInterrupt):
                    raise KeyboardInterrupt(e)
                fail_msg += " {}".format(traceback.format_exc())
                fail_msg += " {}".format(e)
                print(fail_msg)

            passed[k] = pass_flag
            if not pass_flag:
                fail_msgs.append(fail_msg)


        print("[----------] {} tests for component failure detection".format(len(tests)))
        print("[==========]")

        total_passed = sum(passed)
        if total_passed == len(tests):
            print("[  PASSED  ] {} tests for component failure detection".format(len(tests)))
        else:
            print("[  FAILED  ] {} tests, listed below:".format(len(tests)-total_passed))
            for msg in fail_msgs:
                print(msg)


    finally:
        try:
            reference_pipeline.clean_up()
        except (NameError, AttributeError):
            pass
