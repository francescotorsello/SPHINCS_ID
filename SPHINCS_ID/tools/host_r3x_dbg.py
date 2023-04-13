############################################################################
# FORTRAN DEBUG BUILD OPTIONS
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
if cpp_compiler == 'icpc':

  # C++ flags for debugging
  # See also https://www.nas.nasa.gov/hecc/support/kb/recommended-intel-compiler-debugging-options_92.html
  env['CXXFLAGS'] = [ '-O0', '-g', '-std=c++11', '-qopenmp', '-xHOST', \
                      '-qopt-report', '-qopt-report-phase=vec,openmp', \
                      '-Wall', '-m64', '-pedantic', '-Wall', '-Wundef',\
                      '-Wshadow', '-Wcast-qual', '-Wcast-align', \
                      '-Wconversion', '-Winline', '-Wabi=11', \
                      '-Wold-style-cast', '-Woverloaded-virtual', \
                      '-traceback', '-check-uninit', '-ftrapuv', '-debug', \
                      '-debug extended', '-fpe3', '-mp', \
                      '-fp-model strict', \
                      '-align all', '-check bounds', '-assume realloc_lhs',\
                      '-assume protect_parens,minus0', \
                      '-assume no old_maxminloc', '-warn unused', \
                      '-align dcommons' ]

if cpp_compiler == 'gcc' or cpp_compiler == 'g++':

  env['CXXFLAGS'] = [ '-O0', '-g3', '-std=c++11', '-fopenmp', '-Wall', \
                      '-ftree-vectorize','-fopt-info-vec', \
                      '-fopt-info-loop',\
                      '-m64', '-pedantic', '-Wall', \
                      '-Wundef', '-Wshadow', '-Wcast-qual', '-Wcast-align',\
                      '-Wconversion', '-Winline', '-Wabi=11', \
                      '-Wold-style-cast', '-Woverloaded-virtual', \
                      '-Wfatal-errors' ]
