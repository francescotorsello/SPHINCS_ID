################################################################################
# File:       host_sunrise.py
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
# IN THIS FILE, HOST-SPECIFIC (Sunrise) PRODUCTION OPTIONS FOR THE COMPILATION
# OF SPHINCS_ID ARE SET.
# THIS FILE IS SUPPOSED TO BE INCLUDED IN A SConscript FILE.
################################################################################
############################################################################
# FORTRAN PRODUCTION BUILD OPTIONS
############################################################################
if fortran_compiler == 'ifort':

  env['F90FLAGS'] = [ '-O3', '-qopenmp', '-qoverride-limits', \
                      '-heap-arrays', '-qopt-report', \
                      '-qopt-report-phase=vec,openmp', '-no-wrap-margin', \
                      '-warn', '-CB', '-CS', '-fpp', \
                      '-diag-disable=10346', '-static-intel', \
                      '-qopenmp-link=static', '-static', '-xHOST', \
                      '-align array64byte', \
                      incmod_sphincs_bssn, \
                      build_flavour, build_mpi, build_dir, build_host, \
                      build_version ] + mod_dir
#, '-g', '-CB', '-CS', '-traceback'
if fortran_compiler == 'gfortran':

  env['F90FLAGS'] = [ '-O2', '-fopenmp', '-ftree-vectorize', \
                      '-fopt-info-vec', '-fdollar-ok', '-fbounds-check', \
                      '-fopt-info-loop', '-g', '-fbacktrace', '-cpp', \
                      '-ffree-line-length-none', \
                      '-ffixed-line-length-none', \
                      '-Wall', '-Wextra', '-Wimplicit-interface', \
                      '-Wimplicit-procedure', \
                      incmod_sphincs_bssn, build_flavour, \
                      build_mpi, build_dir, build_host, build_version ] \
                      + mod_dir

############################################################################
# C++ PRODUCTION BUILD OPTIONS
############################################################################
if cpp_compiler == 'icpc':

  env['CXXFLAGS'] = [ '-O3', '-g', '-std=c++11', '-qopenmp', '-xHOST', \
                      '-qopt-report', '-qopt-report-phase=vec,openmp', \
                      '-Wall', '-m64', '-DNDEBUG', '-pedantic',
                      '-traceback', '-diag-disable=10397', \
                      '-qoverride-limits', '-static-intel', \
                      '-qopenmp-link=static', '-static', '-fma']
  #-ip, -ipo, -ipo=n
  # icx compiler: best of icc and best of clang
  # -o prog

if cpp_compiler == 'gcc' or cpp_compiler == 'g++':

  env['CXXFLAGS'] = [ '-O3', '-g', '-std=c++11', '-fopenmp', \
                      '-ftree-vectorize','-fopt-info-vec', \
                      '-fopt-info-loop',\
                      '-m64', '-DNDEBUG', '-pedantic', '-Wall', \
                      '-Wundef', '-Wshadow', '-Wcast-qual', '-Wcast-align',\
                      '-Wconversion', '-Winline', '-Wabi=11', \
                      '-Wold-style-cast', '-Woverloaded-virtual', \
                      '-Wfatal-errors' ]
