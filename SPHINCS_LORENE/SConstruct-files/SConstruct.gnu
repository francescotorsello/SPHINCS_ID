import sys
import os
import glob

env = DefaultEnvironment(ENV = {'PATH' : os.environ['PATH']},
                         LINK = 'gfortran',
                         LINKFLAGS = '-g -fopenmp',
                         TOOLS= ['default', 'gfortran'])

# pretty output
if ARGUMENTS.get('VERBOSE') != '1':
  if sys.stdout.isatty():
    env['CXXCOMSTR'] = "\033[92mCompiling\033[0m $TARGET"
    env['F90COMSTR'] = "\033[92mCompiling\033[0m $TARGET"
    env['LINKCOMSTR'] = "\033[94mLinking\033[0m $TARGET"
    env['HDF5COMSTR'] = "\033[95mGenerating\033[0m $TARGET"
  else:
    env['CXXCOMSTR'] = "Compiling $TARGET"
    env['F90COMSTR'] = "Compiling $TARGET"
    env['LINKCOMSTR'] = "Linking $TARGET"
    env['HDF5COMSTR'] = "Generating $TARGET"

# Build options
env['LIBS'] = ['stdc++', 'm', 'gsl', 'lapack', 'fftw3', 'blas', 'gslcblas', 'gfortran', 'lorene_export', 'lorene', 'lorenef77']
env['LIBPATH'] = ['../../../Lorene/Lib']
env['F90FILESUFFIXES']=['.f90','.f']
env['F90']     = 'gfortran'
#env['F90FLAGS'] = ['-O0', '-g3', '-fopenmp', '-ftree-vectorize','-fopt-info-vec', '-fopt-info-loop', '-fbacktrace', '-ftrapping-math', '-fbounds-check', '-ffpe-trap=zero,overflow,underflow','-Wuninitialized','-W','-Wall']
env['F90FLAGS'] = ['-O2', '-fopenmp', '-ftree-vectorize','-fopt-info-vec', '-fopt-info-loop']
env['CXXFLAGS'] = ['-O2', '-g', '-std=c++11', '-fopenmp', '-Wall', '-ftree-vectorize','-fopt-info-vec', '-fopt-info-loop']
env['CXX'] = 'g++'

env['F90PATH'] = [ '.', '../../BSSN', '../SPHINCS_fix_metric', "../SPHINCS_BSSN" ]
env['CPPPATH'] = [ '../../../Lorene/Export/C++/Include', '../../../Lorene/C++/Include/' ]
Progress('Evaluating $TARGET\n')

#lorene_sources = '../../../Lorene/C++/Source/**/*.C'
lorene_sources = [
Glob('../../../Lorene/C++/Source/Binaire/*.C'),
#Glob('../../../Lorene/C++/Source/Binary_xcts/*.C'),
#Glob('../../../Lorene/C++/Source/Star/*.C'),
Glob('../../../Lorene/C++/Source/Tenseur/*.C'),
Glob('../../../Lorene/C++/Source/Valeur/*.C'),
Glob('../../../Lorene/C++/Source/Tensor/*.C'),
Glob('../../../Lorene/C++/Source/Tensor/Scalar/*.C'),
Glob('../../../Lorene/C++/Source/Metric/*.C'),
Glob('../../../Lorene/C++/Source/Metrique/*.C'),
Glob('../../../Lorene/C++/Source/Connection/*.C'),
Glob('../../../Lorene/C++/Source/Diff/*.C'),
Glob('../../../Lorene/C++/Source/Itbl/*.C'),
Glob('../../../Lorene/C++/Source/Tbl/*.C'),
Glob('../../../Lorene/C++/Source/Mtbl/*.C'),
Glob('../../../Lorene/C++/Source/Grille3d/*.C'),
Glob('../../../Lorene/C++/Source/Map/*.C'),
Glob('../../../Lorene/C++/Source/Mtbl_cf/*.C'),
Glob('../../../Lorene/C++/Source/Etoile/*.C'),
Glob('../../../Lorene/C++/Source/Eos/*.C'),
Glob('../../../Lorene/C++/Source/Map_et/*.C'),
Glob('../../../Lorene/C++/Source/Mg3d/*.C'),
Glob('../../../Lorene/C++/Source/Cmp/*.C'),
Glob('../../../Lorene/C++/Source/Matrice/*.C'),
Glob('../../../Lorene/C++/Source/Param/*.C'),
Glob('../../../Lorene/C++/Source/Param_elliptic/*.C'),
Glob('../../../Lorene/C++/Source/Coord/*.C'),
Glob('../../../Lorene/C++/Source/Base_val/*.C'),
Glob('../../../Lorene/C++/Source/Base_vect/*.C'),
Glob('../../../Lorene/C++/Source/Non_class_members/Utilities/*.C'),
Glob('../../../Lorene/C++/Source/Non_class_members/PDE/*.C'),
Glob('../../../Lorene/C++/Source/Non_class_members/Operators/*.C'),
Glob('../../../Lorene/C++/Source/Non_class_members/Coef/*.C'),
Glob('../../../Lorene/C++/Source/Non_class_members/Coef/FFTW3/*.C'),
Glob('../../../Lorene/C++/Source/Non_class_members/Graphics/save_profile.C'),
#Glob('../../../Lorene/C++/Source/Non_class_members/Coef/FFT991/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_helmholtz_minus_2d/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_helmholtz_minus_pseudo_1d/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_poisson_2d/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_poisson_pseudo_1d/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_pois_vect_r/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_sec_order/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_sec_order_r2/*.C'),
Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_vorton/*.C')
]

lorene_sources_bin_ns = [#'../../../Lorene/Export/BinNS/read_bin_ns.C',
                         '../../../Lorene/Export/C++/Source/bin_ns.C',
                         '../../../Lorene/Export/C++/Source/bin_ns_aux.C',
                         '../../../Lorene/Export/C++/Source/write_lines.C']

sources = ['src/module_bindings.f90',
           'src/module_utility.f90',
           'src/module_SPHINCS_LORENE.f90',
           'src/module_bns_id.f90',
           'src/module_particles_id.f90',
           'src/module_formul_3p1_id.f90',
           'src/module_bssn_id.f90',
           'src/submodule_bns_constructor.f90',
           'src/submodule_bns_methods.f90',
           'src/submodule_particles_constructor.f90',
           'src/submodule_particles_methods.f90',
           'src/submodule_formul_3p1_methods.f90',
           'src/submodule_BSSN_id_constructor.f90',
           'src/submodule_BSSN_id_methods.f90',
           '../../BSSN/module_ADM.f90',
           '../../BSSN/module_BSSN.f90',
           '../../BSSN/module_BSSN_parameters.f90',
           '../../BSSN/module_GravityAcceleration.f90',
           '../../BSSN/module_McLachlan.f90',
           '../../BSSN/module_TOV.f90',
           '../../BSSN/module_Tmunu.f90',
           '../../BSSN/module_gravity_grid.f90',
           '../../BSSN/module_map_metric_2_particles.f90',
           '../../BSSN/module_tensor.f90',
           '../../BSSN/module_Numerics.f90',
           '../../BSSN/Hermite_3D_Weights.f',
           '../../BSSN/Hermite5_3D_Weights.f',
           '../../BSSN/module_WENO.f90',
           '../../BSSN/module_NaNChecker.f90',
           '../SPHINCS_BSSN/module_CP_distribute.f90',
           '../SPHINCS_BSSN/module_MLS.f90',
           '../SPHINCS_BSSN/module_SelfRegularization.f90',
           '../SPHINCS_BSSN/module_particle_mesh.f90',
           '../SPHINCS_BSSN/module_map_particles_2_grid.f90',
           '../SPHINCS_BSSN/module_integ_hydro_BSSN.f90',
           '../SPHINCS_BSSN/module_grid_output.f90',
           '../SPHINCS_BSSN/module_gradient.f90',
           '../SPHINCS_BSSN/module_timing.f90',
           '../SPHINCS_BSSN/module_debug_output.f90',
           '../SPHINCS_BSSN/module_grav_acc_check.f90',
           '../SPHINCS_BSSN/module_dyn_metric_on_particles.f90',
           '../../BSSN/ML_BSSN_NV_EvolutionInteriorSplitBy1.cc',
           '../../BSSN/ML_BSSN_NV_EvolutionInteriorSplitBy2.cc',
           '../../BSSN/ML_BSSN_NV_EvolutionInteriorSplitBy3.cc',
           '../../BSSN/ML_BSSN_NV_InitialADMBase1Everywhere.cc',
           '../../BSSN/ML_BSSN_NV_InitialADMBase2Interior.cc',
           '../../BSSN/ML_BSSN_NV_MatterAccelerationInterior.cc',
           '../../BSSN/ML_BSSN_NV_EnforceEverywhere.cc',
           '../../BSSN/ML_BSSN_NV_ADMDerivativesInterior.cc',
           '../../BSSN/ML_BSSN_NV_ConstraintsInterior.cc',
           '../../BSSN/module_Extract_Mass.f90',
           '../../BSSN/bssn_parameters.cc',
           '../SPHINCS_fix_metric/module_files.f90',
           '../SPHINCS_fix_metric/module_constants.f90',
           '../SPHINCS_fix_metric/module_units.f90',
           '../SPHINCS_fix_metric/module_SPHINCS_options.f90',
           '../SPHINCS_fix_metric/module_boundaries.f90',
           '../SPHINCS_fix_metric/module_Numerical_Recipes.f90',
           '../SPHINCS_fix_metric/module_matrix.f90',
           '../SPHINCS_fix_metric/module_3plus1.f90',
           '../SPHINCS_fix_metric/module_kernel_table.f90',
           '../SPHINCS_fix_metric/module_RCB_tree_neighbours.f90',
           '../SPHINCS_fix_metric/module_SPH_variables.f90',
           '../SPHINCS_fix_metric/module_analyze.f90',
           '../SPHINCS_fix_metric/module_EOS.f90',
           '../SPHINCS_fix_metric/module_set_h.f90',
           '../SPHINCS_fix_metric/module_BH_params.f90',
           '../SPHINCS_fix_metric/module_input_output.f90',
           '../SPHINCS_fix_metric/module_recovery_via_pressure.f90',
           '../SPHINCS_fix_metric/module_SPHINCS_SPH_v2.f90',
           '../SPHINCS_fix_metric/module_SPHINCS_version.f90',
           '../SPHINCS_fix_metric/module_SPHINCS_logbook.f90',
           '../SPHINCS_fix_metric/module_relativistic_Eigenvalues.f90',
           '../SPHINCS_fix_metric/module_reconstruction.f90',
           '../SPHINCS_fix_metric/module_deactivate_particles.f90']

#tov_executable = 'setup_TOV.x'
#executable = 'sphincs_bssn.x'
lorene_id_executable = 'setup_lorene_bns_id.x'

#Program(tov_executable, sources+['../SPHINCS_BSSN/Setup_TOV_star.f90'])
#Program(executable, sources+['../SPHINCS_BSSN/SPHINCS_BSSN.f90'])
Program(lorene_id_executable, lorene_sources_bin_ns + lorene_sources + sources + ['src/setup_lorene_id.f90'])

Decider('MD5-timestamp')
