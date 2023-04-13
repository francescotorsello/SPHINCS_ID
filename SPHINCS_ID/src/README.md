# README.md file for the src directory of SPHINCS_ID

The src directory contains the source code for all of the MODULES of SPHINCS_ID.

The code is organized in the following hierarchical structure. Each MODULE has its own directory, named as the MODULE itself, that contains the source files for the MODULE and its SUBMODULES. If a MODULE contains the implementation of a TYPE foo that EXTENDS TYPE bar, then the directory foo is placed inside the directory bar.

In this way, the user gets an idea of the structure of the code just by looking at the directory structure.

The SConscript file in this directory is called by the SConstruct file in the root directory of SPHINCS_ID, during compilation.
