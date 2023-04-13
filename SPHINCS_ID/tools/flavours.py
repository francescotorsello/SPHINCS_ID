################################################################################
# FLAVOUR= 1. Link to LORENE and Kadath libraries, and interpolate from
#             Cartesian grid

if flavour == full_flavour: #

  build_flavour = '-Dflavour=1'

  if host == 'r3x': #

    env['LIBPATH'] = liblorene_dir + libkadath_dir + libsphincs_bssn_dir

    if debug == 'FALSE': #

      env['LIBS'] = ['stdc++', 'm', \
                     'lorene_export', 'lorene', 'lorenef77',  \
                     'kadath', 'gsl', 'lapack', 'fftw3', 'blas', \
                     'gslcblas', 'gfortran', 'sphincs_bssn']

    if debug == 'TRUE': #

      env['LIBS'] = ['stdc++', 'm', 'lorene_export_g', 'lorene_g', \
                     'lorenef77_g', 'kadath-debug', 'gsl', 'lapack', 'fftw3', \
                     'blas', 'gslcblas', 'gfortran', 'sphincs_bssn']

  if host == 'Sunrise': #

    env['LIBPATH'] = liblorene_dir + libkadath_dir + \
                     ['/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/scalapack/2.1.0/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/fftw/3.3.8/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openblas/0.3.12/lib']

    if debug == 'FALSE': #

      env['LIBS'] = ['stdc++', 'm', 'lorene_export', 'lorene', \
                     'lorenef77', 'kadath', 'gsl', 'scalapack', 'openblas', \
                     'gslcblas', 'gfortran', 'fftw3', 'sphincs_bssn']

    if debug == 'TRUE': #

      env['LIBS'] = ['stdc++', 'm', 'lorene_export_g', 'lorene_g', \
                     'lorenef77_g', 'kadath-debug', 'gsl', 'scalapack', \
                     'openblas', 'gslcblas', 'gfortran', 'fftw3', \
                     'sphincs_bssn']

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

      env['LIBS'] = ['stdc++', 'm', \
                     'lorene_export', 'lorene', 'lorenef77', \
                     'gsl', 'lapack', 'fftw3', 'blas', \
                     'gslcblas', 'gfortran', 'sphincs_bssn']

    if debug == 'TRUE': #

      env['LIBS'] = ['stdc++', 'm', 'lorene_export_g', 'lorene_g', \
                     'lorenef77_g', 'gsl', 'lapack', 'fftw3', 'blas', \
                     'gslcblas', 'gfortran', 'sphincs_bssn']

  if host == 'Sunrise': #

    env['LIBPATH'] = liblorene_dir + \
                     ['/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/scalapack/2.1.0/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/fftw/3.3.8/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openblas/0.3.12/lib']

    if debug == 'FALSE': #

      env['LIBS'] = ['stdc++', 'm', 'lorene_export', 'lorene', \
                     'lorenef77', 'gsl', 'scalapack', 'openblas', \
                     'gslcblas', 'gfortran', 'fftw3', 'sphincs_bssn']

    if debug == 'TRUE': #

      env['LIBS'] = ['stdc++', 'm', 'lorene_export_g', 'lorene_g', \
                     'lorenef77_g', 'gsl', 'scalapack', 'openblas', \
                     'gslcblas', 'gfortran', 'fftw3', 'sphincs_bssn']

  sources_flavour = module_bns_lorene + module_diffstar_lorene \
                  + module_ejecta_generic + module_sphincs_id_lorene

################################################################################
# FLAVOUR= 3. Link to Kadath library, and interpolate from Cartesian grid

if flavour == fuka_flavour: #

  build_flavour = '-Dflavour=3'

  if host == 'r3x': #

    env['LIBPATH'] = libkadath_dir + libsphincs_bssn_dir

    if debug == 'FALSE': #

      env['LIBS'] = ['stdc++', 'm', 'kadath', 'gsl', 'fftw3', 'lapack', \
                     'sphincs_bssn']

    if debug == 'TRUE': #

      env['LIBS'] = ['stdc++', 'm', 'kadath-debug', 'gsl', 'fftw3', 'lapack', \
                     'sphincs_bssn']

  if host == 'Sunrise': #

    env['LIBPATH'] = libkadath_dir + \
                     ['/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/scalapack/2.1.0/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openmpi4/fftw/3.3.8/lib', \
                      '/opt/ohpc/pub/libs/gnu8/openblas/0.3.12/lib']

    if debug == 'FALSE': #

      env['LIBS'] = ['stdc++', 'm', 'kadath', 'gsl', 'scalapack', 'openblas', \
                     'fftw3', 'sphincs_bssn']

    if debug == 'TRUE': #

      env['LIBS'] = ['stdc++', 'm', 'kadath', 'gsl', 'scalapack', 'openblas', \
                     'fftw3', 'sphincs_bssn']

  sources_flavour = module_bns_fuka \
                    + module_ejecta_generic + module_sphincs_id_fuka

################################################################################
# FLAVOUR= 4. Interpolate data from Cartesian grid

if flavour == interpolate_flavour: #

  build_flavour= '-Dflavour=4'

  env['LIBPATH'] = libsphincs_bssn_dir

  env['LIBS'] = ['stdc++', 'm', 'sphincs_bssn']

  sources_flavour = module_ejecta_generic + module_sphincs_id_interpolate

################################################################################
