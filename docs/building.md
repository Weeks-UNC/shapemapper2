<!---
NOTE:
If you're reading this, instead try opening README.html in a web browser 
or view this file from within the github repository website.

This is a github-flavored markdown file not meant to be easily readable.
-->

Building
========

Note: All required executables are included with the release tarball,
available on the [release page](https://github.com/Weeks-UNC/shapemapper2/releases).
Download using the `shapemapper-2.1.5.tar.gz` link, _not_ the link for source code only.
Rebuilding should only be needed in very specialized cases.

Runtime dependencies
--------------------
All third-party executables are included in the main release, within
the subdirectory `internals/thirdparty`. If attempting to run from a
source-only release, these will need to be installed at the system level.
 - bash
 - python >= 3.5.3 (`python3` must be visible in the PATH)
 - threaded perl >= 5.22.0 (for bowtie2 wrapper)
 - bowtie2 >= 2.3.0
 - bbmap >= 37.78 (includes bbmerge), requires building JNI components
      - jdk >= 8.0.45 (required for compiling bbmerge JNI components)
 - java >= 8.0.45 (for bbmerge)
 - matplotlib >= 1.5.1 (very recent versions might break some things)
 - pipeviewer >= 1.6.0
 - STAR aligner >= 2.5.2a (Optional. Only useful to speed up alignment against
                           large target sequences.)
 - scikit-learn >= 0.18.1 (Optional. Only needed if running ROC_tests.sh)
 - graphviz >= 2.38.0 (Optional. only needed if running shapemapper with
                       `--render-flowchart` option)
 - ghostscript >= 9.25 (Optional. For debugging use with --render-mutations)


Build requirements
------------------
- cmake >= 3.4.3
- gcc/g++ >= 5.3.0 (must support C++11)
- zlib
- boost >= 1.60.0 libs: filesystem, program_options, iostreams, system
- doxygen >= 1.8.10 (Optional. For rendering c++ documentation)


Building ShapeMapper binaries
-----------------------------

To rebuild ShapeMapper executables
- Delete the folder `build` if present
- Run `internals/bin/build_binaries.sh` or run the following commands,
  starting from inside the top-level shapemapper-2.1.5 folder:
    - `mkdir build`
    - `cd build`
    - `cmake ..`
    - `make`


Source code documentation
-------------------------
- C++ documentation at `internals/cpp-src/doc/html/index.html`
- Python documentation in source code comments (folder `internals/python`)


&nbsp;&nbsp;&nbsp;&nbsp;

[← back to README](../README.md)
