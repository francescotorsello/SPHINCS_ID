#!/bin/bash
################################################################################
# File:       clean-build.sh
# Author:     Francesco Torsello
################################################################################
# Copyright (C) 2020-2023 Francesco Torsello
#
# This file is part of SPHINCS_ID
#
# SPHINCS_ID is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SPHINCS_ID is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SPHINCS_ID. If not, see <https://www.gnu.org/licenses/>.
# The copy of the GNU General Public License should be in the file
# 'COPYING'.
################################################################################
################################################################################
# BASH SCRIPT TO MAKE A CLEAN COMPILATION OF SPHINCS_ID
#
# IT TAKES A MANDATORY ARGUMENT WHICH SPECIFIES WHICH COMPILER TO USE.
# THE ARGUMENT MUST TAKE ONE OF THESE VALUES:
#
# - intel : USE ifort
# - gnu   : USE gfortran
#
# THE SCRIPT MUST BE SOURCED, NOT EXECUTED, IN THE ROOT DIRECTORY OF SPHINCS_ID,
# FOR EXAMPLE WITH THE LINE: 'source tools/clean-build.sh'
#
# NOTE THAT THE SPHINCS_BSSN LIBRARY MUST BE COMPILED WITH THE SAME COMPILER
# AS SPHINCS_ID
################################################################################

# Check that the argument is appropriate
if [ -z "$1" ]; then
    echo ""
    echo "Please specify the argument: the variable 'compilers', which can take 2 values: 'intel' or 'gnu'."
    echo ""
    return
fi

# Run SCons with the option --clean
if [[ $1 == 'intel' ]]; then

    (use ifort; scons compilers=intel --clean)

else

    scons compilers=gnu --clean
    
fi

# SCons does not remove the *.smod files when using the --clean option.
# Hence, remove them and make sure that also the object files *.o,
# the module files *.mod, and the optimization reports *.optrpt are removed.
find . -type f -name '*.o' -delete
find . -type f -name '*.mod' -delete
find . -type f -name '*.smod' -delete
find . -type f -name '*.optrpt' -delete

# Compile SPHINCS_ID
if [[ $1 == 'intel' ]]; then

    (use ifort; scons compilers=intel)

else

    scons compilers=gnu

fi
