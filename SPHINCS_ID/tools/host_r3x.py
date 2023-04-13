############################################################################
# FORTRAN PRODUCTION BUILD OPTIONS
############################################################################
if fortran_compiler == 'ifort':

  env['F90FLAGS'] = [ '-O0', '-qopenmp', '-qoverride-limits', '-xHOST', \
                      '-heap-arrays', '-qopt-report', \
                      '-qopt-report-phase=vec,openmp', '-no-wrap-margin', \
                      '-warn', '-CB', '-CS', '-g', '-traceback', '-fpp', \
                      incmod_sphincs_bssn, \
                      build_flavour, build_mpi, build_dir, build_host, \
                      build_version ] + mod_dir
#, '-CB', '-CS', '-g', '-traceback'
if fortran_compiler == 'gfortran':

  env['F90FLAGS'] = [ '-O0', '-fopenmp', '-ftree-vectorize', \
                      '-fopt-info-vec', '-fdollar-ok', '-fbounds-check', \
                      '-fopt-info-loop', '-g', '-fbacktrace', '-cpp', \
                      '-ffree-line-length-none', \
                      '-ffixed-line-length-none', \
                      '-Wall', '-Wextra', '-Wimplicit-interface', \
                      '-Wimplicit-procedure', \
                      incmod_sphincs_bssn, \
                      build_flavour, build_mpi, \
                      build_dir, build_host, build_version ] + mod_dir

############################################################################
# C++ PRODUCTION BUILD OPTIONS
############################################################################
if cpp_compiler == 'icpc':

  env['CXXFLAGS'] = [ '-O0', '-g', '-std=c++11', '-qopenmp', '-xHOST', \
                      '-qopt-report', '-qopt-report-phase=vec,openmp', \
                      '-Wall', '-m64', '-DNDEBUG', '-pedantic', '-Wall', \
                      '-W', '-Wundef', '-Wshadow', '-Wcast-qual', \
                      '-Wconversion', '-Winline', '-Woverloaded-virtual' ]
  #-ip, -ipo, -ipo=n
  # icx compiler: best of icc and best of clang
  # -o prog

if cpp_compiler == 'gcc' or cpp_compiler == 'g++':

  env['CXXFLAGS'] = [ '-O0', '-g', '-std=c++11', '-fopenmp', '-Wall', \
                      '-ftree-vectorize','-fopt-info-vec', \
                      '-fopt-info-loop',\
                      '-m64', '-DNDEBUG', '-pedantic', '-Wall', \
                      '-Wundef', '-Wshadow', '-Wcast-qual', '-Wcast-align',\
                      '-Wconversion', '-Winline', '-Wabi=11', \
                      '-Wold-style-cast', '-Woverloaded-virtual', \
                      '-Wfatal-errors' ]
