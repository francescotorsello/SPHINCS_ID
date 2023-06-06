################################################################################
# SConstruct file to build SPHINCS_ID
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
# ACRONYMS
################################################################################
# ID          : Initial Data
# SPHINCS_ID  : Smoothed Particle Hydrodynamics IN Curved Spacetime - ID builder
################################################################################

import sys
import os
import glob
import socket

################################################################################
# DEFAULT OPTIONS
################################################################################

# A 'flavour' of SPHINCS_ID is defined as a set of its MODULES
# The following flavours are implemented

# Include all the MODULES
full_flavour        = '1'

# Include the MODULES that link to the LORENE library to read the ID, and
# interpolate the ID from a Cartesian, uniform grid
lorene_flavour      = '2'

# Include the MODULES that link to the FUKA library to read the ID, and
# interpolate the ID from a Cartesian, uniform grid
fuka_flavour        = '3'

# Include the MODULES that interpolate from a Cartesian, uniform grid
interpolate_flavour = '4'

# Set the default flavour
default_flavour = full_flavour

# 'TRUE'  : build in DEBUG MODE
# 'FALSE' : build in PRODUCTION MODE
default_debug = 'FALSE'

# Which compilers to use
#   Fortran : 'gfortran', 'ifort'
#   C++     : 'g++', 'gcc' or 'icpc'
default_fortran_compiler = 'ifort'
#default_cpp_compiler     = 'icpc'

# 'TRUE'  : enables verbose output
# 'FALSE' : disables verbose output
default_verbose = 'FALSE'

# Specifies the host machine where SPHINCS_ID is used
# 'r3x'     : The r30, r31, r32 machines @ Department of Astronomy,
#             Stockholm University
# 'Sunrise' : The Sunrise HPC cluster @ Fysikum, Stockholm University
#             The following line needs to be added to the user's .bash_profile
#             file, so that the necessary modules are loaded to the user's local
#             environment on Sunrise:
#   module load -openmpi3 +openmpi4 +gsl +fftw +scalapack +libgfortran +pgplot
#             Regarding compilation of FUKA on Sunrise, follow the
#             instructions at (on Sunrise):
#               /cfs/home/pg/CHAP/compile-kadath-scripts/
current_host= socket.gethostname()

if current_host == 'r30' or current_host == 'r31' or current_host == 'r32':
  default_host = 'r3x'
  MPI_ranks    = '40'

elif current_host == 'sol-login.fysik.su.se':
  default_host = 'Sunrise'
  MPI_ranks    = '128'

else:
  default_host = 'r3x'
  MPI_ranks    = '40'

# Version of SPHINCS_ID
version = 'v2.0'

################################################################################

################################################################################
# USEFUL VARIABLES
################################################################################

# Home directory for SPHINCS_ID
working_dir = os.getcwd()

################################################################################

################################################################################
# CUSTOM OPTIONS
################################################################################

flavour = ARGUMENTS.get('flavour')
if flavour is None: flavour = default_flavour

debug = ARGUMENTS.get('debug')
if debug is None: debug = default_debug

fortran_compiler = ARGUMENTS.get('fortran_compiler')
if fortran_compiler is None: fortran_compiler = default_fortran_compiler

#cpp_compiler = ARGUMENTS.get('cpp_compiler')
#if cpp_compiler is None: cpp_compiler = default_cpp_compiler

compilers = ARGUMENTS.get('compilers')
if compilers == 'intel':
    fortran_compiler = 'ifort'
    #cpp_compiler     = 'icpc'
if compilers == 'gnu':
    fortran_compiler = 'gfortran'
    #cpp_compiler     = 'g++'

verbose = ARGUMENTS.get('verbose')
if verbose is None: verbose = default_verbose

host = ARGUMENTS.get('host')
if host is None: host = default_host

################################################################################

################################################################################

################################################################################
# ENVIRONMENT
################################################################################

if fortran_compiler == 'ifort': #

  if host == 'r3x': #

    env = DefaultEnvironment(

      ENV       = { 'PATH'               : os.environ['PATH'],
                    'INTEL_LICENSE_FILE' : os.environ['INTEL_LICENSE_FILE'],
                    'LD_LIBRARY_PATH'    : os.environ['LD_LIBRARY_PATH']  },
      LINK      = 'ifort',
      LINKFLAGS = '-g -qopenmp',
      TOOLS     = ['default', 'ifort']

    )

  if host == 'Sunrise': #

    env = DefaultEnvironment(

      ENV       = { 'PATH'               : os.environ['PATH'],
                    'LD_LIBRARY_PATH'    : os.environ['LD_LIBRARY_PATH']  },
      LINK      = 'ifort',
      LINKFLAGS = '-g -traceback -qopenmp -static-intel -qopenmp-link=static \
                   -static-libgcc -static-libstdc++',
      TOOLS     = ['default', 'ifort']

    )

if fortran_compiler == 'gfortran': #

  env = DefaultEnvironment(

    ENV       = { 'PATH' : os.environ['PATH'] },
    LINK      = 'gfortran',
    LINKFLAGS = '-g -fopenmp',
    TOOLS     = ['default', 'gfortran']

  )

if verbose == 'FALSE':

  if sys.stdout.isatty():

    env['CXXCOMSTR']  = "\033[92mCompiling\033[0m $TARGET"
    env['F90COMSTR']  = "\033[92mCompiling\033[0m $TARGET"
    env['LINKCOMSTR'] = "\033[94mLinking\033[0m $TARGET"
    env['HDF5COMSTR'] = "\033[95mGenerating\033[0m $TARGET"

  else:

    env['CXXCOMSTR']  = "Compiling $TARGET"
    env['F90COMSTR']  = "Compiling $TARGET"
    env['LINKCOMSTR'] = "Linking $TARGET"
    env['HDF5COMSTR'] = "Generating $TARGET"

################################################################################

Export('env flavour full_flavour lorene_flavour fuka_flavour ' \
       + 'interpolate_flavour version MPI_ranks working_dir host ' \
       + 'fortran_compiler debug') # cpp_compiler
################################################################################

SConscript('src/SConscript', variant_dir='build', duplicate=False)
