################################################################################
# File:       flavours.py
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
# IN THIS FILE, FLAVOUR-SPECIFIC OPTIONS FOR THE COMPILATION OF SPHINCS_ID
# ARE SET. THIS FILE IS SUPPOSED TO BE INCLUDED IN A SConscript FILE.
################################################################################
# FLAVOUR= 1. Link to LORENE and Kadath libraries, and interpolate from
#             Cartesian grid

if flavour == full_flavour: #

  build_flavour = '-Dflavour=1'

  if host == 'r3x': #

    env['LIBPATH'] = liblorene_dir + libkadath_dir + libsphincs_bssn_dir

    if debug == 'FALSE': #

      env['LIBS'] = ['lorene_export', 'lorene', 'lorenef77',  \
                     'kadath', 'gsl', 'lapack', 'fftw3', 'blas', \
                     'gslcblas', 'gfortran', 'sphincs_bssn', 'stdc++', 'm']

    if debug == 'TRUE': #

      # It may happen that the *_g libraries by LORENE stop the execution due to
      # arithmetic exceptions in LORENE. This may not allow to debug SPHINCS_ID.
      # If that happens, link to the production LORENE libraries, listed above.
      env['LIBS'] = ['lorene_export_g', 'lorene_g', 'lorenef77_g', \
                     'kadath-debug', 'gsl', 'lapack', 'fftw3', 'blas', \
                     'gslcblas', 'gfortran', 'sphincs_bssn', 'stdc++', 'm']

  if host == 'Sunrise': #

    env['LIBPATH'] = liblorene_dir + libkadath_dir + libsphincs_bssn_dir + \
                     ['/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/scalapack/2.1.0/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/fftw/3.3.8/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openblas/0.3.12/lib']

    if debug == 'FALSE': #

      env['LIBS'] = ['lorene_export', 'lorene', 'lorenef77', 'kadath', \
                     'gsl', 'scalapack', 'openblas', 'gslcblas',  \
                     'gfortran', 'fftw3', 'sphincs_bssn', 'stdc++', 'm']

    if debug == 'TRUE': #

      # It may happen that the *_g libraries by LORENE stop the execution due to
      # arithmetic exceptions in LORENE. This may not allow to debug SPHINCS_ID.
      # If that happens, link to the production LORENE libraries, listed above.
      env['LIBS'] = ['lorene_export_g', 'lorene_g', \
                     'lorenef77_g', 'kadath-debug', 'gsl', 'scalapack', \
                     'openblas', 'gslcblas', 'gfortran', 'fftw3', \
                     'sphincs_bssn', 'stdc++', 'm']

  sources_flavour = module_bns_lorene + module_diffstar_lorene \
                  + module_bns_fuka \
                  + module_ejecta_generic + module_sphincs_id_full

################################################################################
# FLAVOUR= 2. Link to LORENE library, and interpolate from Cartesian grid

if flavour == lorene_flavour: #

  build_flavour = '-Dflavour=2'

  if host == 'r3x': #

    env['LIBPATH'] = liblorene_dir + libsphincs_bssn_dir

    if debug == 'FALSE': #

      env['LIBS'] = ['lorene_export', 'lorene', 'lorenef77', \
                     'gsl', 'lapack', 'fftw3', 'blas', \
                     'gslcblas', 'sphincs_bssn', \
                     'stdc++', 'm']

    if debug == 'TRUE': #

      # It may happen that the *_g libraries by LORENE stop the execution due to
      # arithmetic exceptions in LORENE. This may not allow to debug SPHINCS_ID.
      # If that happens, link to the production LORENE libraries, listed above.
      env['LIBS'] = ['lorene_export_g', 'lorene_g', \
                     'lorenef77_g', 'gsl', 'lapack', 'fftw3', 'blas', \
                     'gslcblas', 'gfortran', 'sphincs_bssn', 'stdc++', 'm']

  if host == 'Sunrise': #

    env['LIBPATH'] = liblorene_dir + libsphincs_bssn_dir + \
                     ['/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/scalapack/2.1.0/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/fftw/3.3.8/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openblas/0.3.12/lib']

    if debug == 'FALSE': #

      env['LIBS'] = ['lorene_export', 'lorene', 'lorenef77', \
                     'gsl', 'scalapack', 'openblas', 'gslcblas', \
                     'gfortran', 'fftw3', 'sphincs_bssn', 'stdc++', 'm']

    if debug == 'TRUE': #

      # It may happen that the *_g libraries by LORENE stop the execution due to
      # arithmetic exceptions in LORENE. This may not allow to debug SPHINCS_ID.
      # If that happens, link to the production LORENE libraries, listed above.
      env['LIBS'] = ['lorene_export_g', 'lorene_g', 'lorenef77_g', \
                     'gsl', 'scalapack', 'openblas', 'gslcblas', \
                     'gfortran', 'fftw3', 'sphincs_bssn', 'stdc++', 'm']

  sources_flavour = module_bns_lorene + module_diffstar_lorene \
                  + module_ejecta_generic + module_sphincs_id_lorene

################################################################################
# FLAVOUR= 3. Link to Kadath library, and interpolate from Cartesian grid

if flavour == fuka_flavour: #

  build_flavour = '-Dflavour=3'

  if host == 'r3x': #

    env['LIBPATH'] = libkadath_dir + libsphincs_bssn_dir

    if debug == 'FALSE': #

      env['LIBS'] = ['kadath', 'gsl', 'fftw3', 'lapack', \
                     'sphincs_bssn', 'stdc++', 'm']

    if debug == 'TRUE': #

      env['LIBS'] = ['kadath-debug', 'gsl', 'fftw3', 'lapack', \
                     'sphincs_bssn', 'stdc++', 'm']

  if host == 'Sunrise': #

    env['LIBPATH'] = libkadath_dir + libsphincs_bssn_dir + \
                     ['/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/scalapack/2.1.0/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/fftw/3.3.8/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openblas/0.3.12/lib']

    if debug == 'FALSE': #

      env['LIBS'] = ['kadath', 'gsl', 'scalapack', 'openblas', \
                     'fftw3', 'sphincs_bssn', 'stdc++', 'm']

    if debug == 'TRUE': #

      env['LIBS'] = ['kadath', 'gsl', 'scalapack', 'openblas', \
                     'fftw3', 'sphincs_bssn', 'stdc++', 'm']

  sources_flavour = module_bns_fuka \
                    + module_ejecta_generic + module_sphincs_id_fuka

################################################################################
# FLAVOUR= 4. Interpolate data from Cartesian grid

if flavour == interpolate_flavour: #

  build_flavour= '-Dflavour=4'

  env['LIBPATH'] = libsphincs_bssn_dir

  env['LIBS'] = ['sphincs_bssn', 'stdc++', 'm']

  sources_flavour = module_ejecta_generic + module_sphincs_id_interpolate

################################################################################
