<!-- README.md file for tools directory of SPHINCS_ID -->

## tools

Tools meant to be used during the compilation of SPHINCS_ID, or during the generation of the documentation, are supposed to be placed in the tools directory.

The Python files \*.py are used by the src/SConscript file (they are not Python packages; they are literally simply included in src/SConscript). They allow for a cleaner and more intuitive structure for src/SConscript.

The bash script clean-build.sh can be sourced (not executed) in the root directory of SPHINCS_ID, to make a clean compilation of the code. See the file itself for details.

The file documentation_sphincs_id.md is the Project File needed by FORD to produce the documentation of SPHINCS_ID.
