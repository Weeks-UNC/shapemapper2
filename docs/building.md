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
Download using the `shapemapper-2.1.5.tar.gz` link.


Building from source
--------------------
Navigate to `shapemapper2/internals/thirdparty_helper/`

Run:
`get_miniconda.sh`

Upon completion run:
`build_thirdparty.sh`

Navigate to `shapemapper2/internals/bin`

Run:
`build_binaries.sh`


Building ShapeMapper binaries
-----------------------------

To build ShapeMapper executables
- Delete the folder `shapemapper2/build` if present

Navigate to `shapemapper2/internals/bin`
Run:
`build_binaries.sh`


Source code documentation
-------------------------
- C++ documentation at `internals/cpp-src/doc/html/index.html`
- Python documentation in source code comments (folder `internals/python`)


&nbsp;&nbsp;&nbsp;&nbsp;

[← back to README](../README.md)
