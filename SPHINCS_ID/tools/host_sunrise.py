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
