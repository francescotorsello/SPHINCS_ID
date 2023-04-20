# File:       options_r3x_dbg.py
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
# IN THIS FILE, HOST-SPECIFIC (r3x) DEBUGGING OPTIONS FOR THE COMPILATION OF
# SPHINCS_ID ARE SET. THIS FILE IS SUPPOSED TO BE INCLUDED IN A SConscript FILE.
################################################################################
############################################################################
# FORTRAN DEBUG BUILD OPTIONS
# (NOT USED ANYMORE AS THE C++ SOURCES ARE NOW COMPILED IN SPHINCS_BSSN ONLY)
############################################################################
if fortran_compiler == 'ifort':

  env['F90FLAGS'] = [ '-O0', '-g', '-qopenmp', '-qoverride-limits', \
                      '-xHOST', '-heap-arrays', '-qopt-report', \
                      '-qopt-report-phase=vec,openmp', '-no-wrap-margin', \
                      '-CB', '-CS', '-traceback', '-fpe0', '-warn', \
                      '-debug extended', '-fpp', \
                      incmod_sphincs_bssn, build_flavour, \
                      build_mpi, build_dir, build_host, build_version ] \
                      + mod_dir

  # More F90 flags for debugging
  #env['F90FLAGS'] = [ '-r8','-g','-qopenmp','-O0', '-ftz', \
  #                    '-heap-arrays 5000', '-fno-alias', '-traceback', \
  #                    '-debug', '-debug all', '-nolib-inline', \
  #                    '-align all',\
  #                    '-check bounds', '-fno-inline-functions', \
  #                    '-prec-div',\
  #                    '-prec-sqrt', '-fp-stack-check', \
  #                    '-init=snan,arrays', \
  #                    '-traceback', '-fPIC', '-fpe3', '-heap-arrays 32', \
  #                    '-assume realloc_lhs', \
  #                    '-assume protect_parens,minus0',\
  #                    '-assume no old_maxminloc', '-warn unused', \
  #                    '-align dcommons',  '-xHOST', '-qopt-report', \
  #                    '-qopt-report-phase=vec,openmp', '-fp-model strict', \
  #                    '-no-wrap-margin', '-ftrapuv', '-mp' ]

if fortran_compiler == 'gfortran':

  env['F90FLAGS'] = [ '-O0', '-g3', '-fopenmp', '-ftree-vectorize', \
                      '-fopt-info-vec', '-fopt-info-loop', '-fbacktrace', \
                      '-ftrapping-math', '-fbounds-check', \
                      '-ffpe-trap=zero,overflow,underflow', \
                      '-Wuninitialized','-W','-Wall', '-Wextra', '-cpp', \
                      '-ffree-line-length-none', \
                      '-ffixed-line-length-none', '-fdollar-ok', \
                      mod_dir, incmod_sphincs_bssn, \
                      build_flavour, build_mpi, build_dir, build_host, \
                      build_version ] + mod_dir

############################################################################
# C++ DEBUG BUILD OPTIONS
############################################################################
#if cpp_compiler == 'icpc':
#
#  # C++ flags for debugging
#  # See also https://www.nas.nasa.gov/hecc/support/kb/recommended-intel-compiler-debugging-options_92.html
#  env['CXXFLAGS'] = [ '-O0', '-g', '-std=c++11', '-qopenmp', '-xHOST', \
#                      '-qopt-report', '-qopt-report-phase=vec,openmp', \
#                      '-Wall', '-m64', '-pedantic', '-Wall', '-Wundef',\
#                      '-Wshadow', '-Wcast-qual', '-Wcast-align', \
#                      '-Wconversion', '-Winline', '-Wabi=11', \
#                      '-Wold-style-cast', '-Woverloaded-virtual', \
#                      '-traceback', '-check-uninit', '-ftrapuv', '-debug', \
#                      '-debug extended', '-fpe3', '-mp', \
#                      '-fp-model strict', \
#                      '-align all', '-check bounds', '-assume realloc_lhs',\
#                      '-assume protect_parens,minus0', \
#                      '-assume no old_maxminloc', '-warn unused', \
#                      '-align dcommons' ]
#
#if cpp_compiler == 'gcc' or cpp_compiler == 'g++':
#
#  env['CXXFLAGS'] = [ '-O0', '-g3', '-std=c++11', '-fopenmp', '-Wall', \
#                      '-ftree-vectorize','-fopt-info-vec', \
#                      '-fopt-info-loop',\
#                      '-m64', '-pedantic', '-Wall', \
#                      '-Wundef', '-Wshadow', '-Wcast-qual', '-Wcast-align',\
#                      '-Wconversion', '-Winline', '-Wabi=11', \
#                      '-Wold-style-cast', '-Woverloaded-virtual', \
#                      '-Wfatal-errors' ]
