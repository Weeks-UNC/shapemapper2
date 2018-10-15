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

# hack so we can import from python files in pyshapemap
# folder, even though this script is a bit isolated
this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(this_dir, '../python'))

import pyshapemap.pipeline_arg_parser as ap
from pyshapemap.util import Logger, timestamp

from pyshapemap.pipeline import run_fail_test



# FIXME: suppress misleading BrokenPipeError warnings
#Exception ignored in: <_io.TextIOWrapper name='<stdout>' mode='w' encoding='UTF-8'>
#BrokenPipeError: [Errno 32] Broken pipe


if __name__ == "__main__":

    latency = 0.1 # default is 0.1. decreasing latency results in more test failures
    success_term_pause = 0.05 # default 0.05
    fail_term_pause = 0.5 # default 0.5
    timeout = 30 # default is 30
    exit_early = True # whether to run all failure tests or quit after
                       # first test giving incorrect results
    quiet = True # if False, display log messages while running, and
                 # append --verbose to args, otherwise
                 # only display error messages on test failure)

    limit_to_component = None # only run fail test for specific component, e.g. '1.7.1'
    # intermittent test failures for various StarAligner components
    # - sometimes times out, sometimes broken pipe
    # - eliminating shared memory index seems to get rid of timeouts,
    #   probably deadlocking


    try:
        args = sys.argv[1:]
        if not quiet:
            args.append('--verbose')
        out_folder, _ = ap.get_out_path(args)
        temp_folder, _ = ap.get_temp_path(args)

        reference_pipeline, _ = ap.construct(args,
                                             skip_flowchart=True,
                                             skip_setup=True)

        names, _, component_locations =reference_pipeline.map_pipeline_tree()

        tests = []
        for n in range(len(component_locations)):
            # skip some components that don't do much input validation
            skip = False
            name = names[n]
            for s in ["Interleaver",
                      "Deinterleaver",
                      "Append",
                      "ProgressMonitor"]:
                if s in name:
                    skip = True
                    break
            if limit_to_component is not None and component_locations[n] != limit_to_component:
                skip = True
            #if "StarAligner" not in name:
            #    skip = True

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

            name = names[n]

            pass_flag = False
            fail_msg = "ERROR: unable to construct pipeline."
            pipeline = None
            try:
                pipeline, _ = ap.construct(args,
                                           skip_flowchart=True,
                                           skip_setup=True)
                pipeline.name = pipeline.get_name()+"_"+p.replace('.','-')
                pipeline.flowchart_path = os.path.join(pipeline.out,
                                                       pipeline.get_name()+"_flowchart.svg")

                if pipeline is not None:
                    fail_msg = "ERROR: unable to run pipeline fail test."
                    pass_flag, success, err = run_fail_test(pipeline=pipeline,
                                                            failing_module_loc=p,
                                                            failing_module_name=name,
                                                            quiet=quiet,
                                                            latency=latency,
                                                            timeout=timeout,
                                                            success_term_pause=success_term_pause,
                                                            fail_term_pause=fail_term_pause)

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
                if exit_early:
                    break


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
