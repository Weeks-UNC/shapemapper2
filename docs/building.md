<!---
NOTE:
If you're reading this, instead try opening README.html in a web browser 
or view this file from within the github repository website.

This is a github-flavored markdown file not meant to be easily readable.
-->

Building
========

Note: All required executables are included with the release tarball.
Rebuilding should only be needed in very specialized cases.

Third-party dependencies
------------------------
All third-party executables are included in the release, within
the subdirectory `internals/thirdparty`. Deleting this folder and running
`internals/install/build_thirdparty.sh` will regenerate these files (requires 
an active Internet connection). Most packages will be downloaded
in executable form using conda or from the package author's preferred
release URL, and a few components will be built from source.

- Bowtie2>=2.3.0 (2.3.4.3 included)
- STAR>=2.5.2a
- Python>=3.5 (3.5.3 included)
- matplotlib>=1.5.1
- numpy
- BBmerge>=37.78 (requires java)
- Pipe Viewer>=1.6.0
- Graphviz (optional: for debugging use with --render-flowchart)
- Ghostscript (optional: for debugging use with --render-mutations)
- scikit-learn>=0.18.1 (optional: for area under ROC curve tests)

In the rare case that the provided third-party executables are not 
compatible with your platform (for example, if you get strange "GLIBC-2.XX" 
or "ABI incompatible" errors), you're on your own. Delete the `internals/thirdparty` 
folder entirely, and then install the listed packages at the system level 
through your package manager or by building from source.


Source code documentation
-------------------------
- C++ documentation at `internals/cpp-src/doc/html/index.html`
- Python documentation in source code comments (folder `internals/python`)


C++ build requirements
----------------------
- recent `cmake` (>=3.3)
- `gcc` compiler with C++11 support
- `git` with https support
- `build-essential` (should include `dpkg-dev`, `g++`, `libc-dev`, `libstdc++`)
- `libboost-dev` (>=1.60.0)
- `libboost-filesystem-dev`
- `libboost-iostreams-dev`
- `libboost-program-options-dev`
- `doxygen` (optional)


Building C++ modules
--------------------

To rebuild ShapeMapper executables
- Delete the folders `build` and `build_deps` if present
- Run `internals/install/build_binaries.sh --use-system-libs` 
- This assumes the tools and libraries listed above are already 
  present at the system level.

Alternatively, run `internals/install/build_all.sh`
- This will take substantial disk space and time
- Compiler and library dependencies listed above (with the exception
  of Doxygen) will be downloaded and built locally within the 
  `internals/build_deps` folder


&nbsp;&nbsp;&nbsp;&nbsp;

[← back to README](../README.md)
