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

env['F90']= 'gfortran'

env['CXX']= 'g++'

# F90 flags for debugging
#env['F90FLAGS'] = ['-O0', '-g3', '-fopenmp', '-ftree-vectorize','-fopt-info-vec', '-fopt-info-loop', '-fbacktrace', '-ftrapping-math', '-fbounds-check', '-ffpe-trap=zero,overflow,underflow','-Wuninitialized','-W','-Wall', '-cpp', '-ffree-line-length-none', '-ffixed-line-length-none']

# F90 flags for production
env['F90FLAGS'] = ['-O3', '-fopenmp', '-ftree-vectorize','-fopt-info-vec', '-fopt-info-loop', '-g', '-fbacktrace', '-cpp', '-ffree-line-length-none', '-ffixed-line-length-none']
#, '-cpp', '-ffree-line-length-none', '-ffixed-line-length-none', '-E', '-dM'

# C++ flags for production (union of the flags in the SPHINCS SConstruct files
# and the flags in the file local_settings in the LORENE home directory)
env['CXXFLAGS'] = ['-O3', '-g', '-std=c++11', '-fopenmp', '-Wall', '-ftree-vectorize','-fopt-info-vec', '-fopt-info-loop', '-m64', '-DNDEBUG', '-pedantic', '-Wall', '-W', '-Wundef', '-Wshadow', '-Wcast-qual', '-Wcast-align', '-Wconversion', '-Winline', '-Wabi=11', '-Wold-style-cast', '-Woverloaded-virtual']

env['F90PATH'] = [ '.', '../../BSSN', '../SPHINCS_fix_metric', "../SPHINCS_BSSN" ]
env['CPPPATH'] = [ '../../../Lorene/Export/C++/Include', '../../../Lorene/C++/Include/' ]
Progress('Evaluating $TARGET\n')

#lorene_sources = '../../../Lorene/C++/Source/**/*.C'
#lorene_sources = [
#Glob('../../../Lorene/C++/Source/Binaire/*.C'),
##Glob('../../../Lorene/C++/Source/Binary_xcts/*.C'),
##Glob('../../../Lorene/C++/Source/Star/*.C'),
#Glob('../../../Lorene/C++/Source/Tenseur/*.C'),
#Glob('../../../Lorene/C++/Source/Valeur/*.C'),
#Glob('../../../Lorene/C++/Source/Tensor/*.C'),
#Glob('../../../Lorene/C++/Source/Tensor/Scalar/*.C'),
#Glob('../../../Lorene/C++/Source/Metric/*.C'),
#Glob('../../../Lorene/C++/Source/Metrique/*.C'),
#Glob('../../../Lorene/C++/Source/Connection/*.C'),
#Glob('../../../Lorene/C++/Source/Diff/*.C'),
#Glob('../../../Lorene/C++/Source/Itbl/*.C'),
#Glob('../../../Lorene/C++/Source/Tbl/*.C'),
#Glob('../../../Lorene/C++/Source/Mtbl/*.C'),
#Glob('../../../Lorene/C++/Source/Grille3d/*.C'),
#Glob('../../../Lorene/C++/Source/Map/*.C'),
#Glob('../../../Lorene/C++/Source/Mtbl_cf/*.C'),
#Glob('../../../Lorene/C++/Source/Etoile/*.C'),
#Glob('../../../Lorene/C++/Source/Eos/*.C'),
#Glob('../../../Lorene/C++/Source/Map_et/*.C'),
#Glob('../../../Lorene/C++/Source/Mg3d/*.C'),
#Glob('../../../Lorene/C++/Source/Cmp/*.C'),
#Glob('../../../Lorene/C++/Source/Matrice/*.C'),
#Glob('../../../Lorene/C++/Source/Param/*.C'),
#Glob('../../../Lorene/C++/Source/Param_elliptic/*.C'),
#Glob('../../../Lorene/C++/Source/Coord/*.C'),
#Glob('../../../Lorene/C++/Source/Base_val/*.C'),
#Glob('../../../Lorene/C++/Source/Base_vect/*.C'),
#Glob('../../../Lorene/C++/Source/Non_class_members/Utilities/*.C'),
#Glob('../../../Lorene/C++/Source/Non_class_members/PDE/*.C'),
#Glob('../../../Lorene/C++/Source/Non_class_members/Operators/*.C'),
#Glob('../../../Lorene/C++/Source/Non_class_members/Coef/*.C'),
#Glob('../../../Lorene/C++/Source/Non_class_members/Coef/FFTW3/*.C'),
#Glob('../../../Lorene/C++/Source/Non_class_members/Graphics/save_profile.C'),
##Glob('../../../Lorene/C++/Source/Non_class_members/Coef/FFT991/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_helmholtz_minus_2d/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_helmholtz_minus_pseudo_1d/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_poisson_2d/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_poisson_pseudo_1d/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_pois_vect_r/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_sec_order/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_sec_order_r2/*.C'),
#Glob('../../../Lorene/C++/Source/Ope_elementary/Ope_vorton/*.C')
#]

lorene_sources_bin_ns = ['../../../Lorene/Export/C++/Source/bin_ns.C',
                         '../../../Lorene/Export/C++/Source/bin_ns_aux.C',
                         '../../../Lorene/Export/C++/Source/write_lines.C']

sources = ['src/module_utility.f90',
           'src/module_id_base.f90',
           'src/module_bns_base.f90',
           'src/module_bns_lorene.f90',
           'src/module_particles_id.f90',
           'src/module_formul_3p1_id.f90',
           'src/module_bssn_id.f90',
           #'src/submodule_idbase_mass_profile.f90',
           'src/submodule_bns_base_access.f90',
           'src/submodule_bns_base_mass_profile.f90',
           'src/submodule_bns_lorene_constructor.f90',
           'src/submodule_bns_lorene_import.f90',
           'src/submodule_bns_lorene_memory.f90',
           'src/submodule_bns_lorene_access.f90',
           'src/submodule_bns_lorene_params.f90',
           'src/submodule_particles_constructor.f90',
           'src/submodule_particles_memory.f90',
           'src/submodule_particles_lattices.f90',
           'src/submodule_particles_spherical_surfaces.f90',
           'src/submodule_particles_sph_variables.f90',
           'src/submodule_particles_access.f90',
           'src/submodule_particles_compose.f90',
           'src/submodule_particles_redistribute_nu.f90',
           'src/submodule_particles_apm.f90',
           'src/submodule_formul_3p1_constructor.f90',
           'src/submodule_formul_3p1_access.f90',
           'src/submodule_formul_3p1_analysis.f90',
           'src/submodule_bssn_id_constructor.f90',
           'src/submodule_bssn_id_variables.f90',
           'src/submodule_bssn_id_constraints.f90',
           'src/submodule_bssn_id_memory.f90',
           'src/module_SPHINCS_LORENE.f90',
           '../../BSSN/module_mesh_refinement.f90',
           #'../../BSSN/module_ADM.f90',
           '../../BSSN/module_ADM_refine.f90',
           #'../../BSSN/module_BSSN.f90',
           '../../BSSN/module_BSSN_refine.f90',
           '../../BSSN/module_BSSN_parameters.f90',
           #'../../BSSN/module_GravityAcceleration.f90',
           '../../BSSN/module_GravityAcceleration_refine.f90',
           #'../../BSSN/module_McLachlan.f90',
           '../../BSSN/module_McLachlan_refine.f90',
           #'../../BSSN/module_TOV.f90',
           '../../BSSN/module_TOV_refine.f90',
           #'../../BSSN/module_Tmunu.f90',
           '../../BSSN/module_Tmunu_refine.f90',
           #'../../BSSN/module_gravity_grid.f90',
           #'../../BSSN/module_map_metric_2_particles.f90',
           '../../BSSN/module_map_metric_2_particles_refine.f90',
           '../../BSSN/module_evolution_parameters.f90',
           '../../BSSN/submodule_prolongation.f90',
           '../../BSSN/submodule_restriction.f90',
           '../../BSSN/submodule_utilities.f90',
           '../../BSSN/submodule_parameter.f90',
           '../../BSSN/submodule_allocation.f90',
           '../../BSSN/submodule_output.f90',
           '../../BSSN/module_tensor.f90',
           '../../BSSN/module_Numerics.f90',
           '../../BSSN/module_Hermite_refine.f90',
           '../../BSSN/Hermite_3D_Weights.f',
           '../../BSSN/Hermite5_3D_Weights.f',
           '../../BSSN/module_WENO.f90',
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
           '../../BSSN/submodule__Extract_Mass_impl.f90',
           '../../BSSN/fermi_step_f.f90',
           '../../BSSN/bssn_parameters.cc',
           #'../../BSSN/module_NaNChecker.f90',
           '../SPHINCS_BSSN/module_CP_distribute.f90',
           '../SPHINCS_BSSN/module_MLS.f90',
           '../SPHINCS_BSSN/module_SelfRegularization.f90',
           '../SPHINCS_BSSN/module_particle_mesh.f90',
           '../SPHINCS_BSSN/module_particle_mesh_hash_grid.f90',
           '../SPHINCS_BSSN/module_map_particles_2_grid.f90',
           '../SPHINCS_BSSN/module_integ_hydro_BSSN.f90',
           #'../SPHINCS_BSSN/module_grid_output.f90',
           '../SPHINCS_BSSN/module_gradient.f90',
           '../SPHINCS_BSSN/module_timing.f90',
           '../SPHINCS_BSSN/module_debug_output.f90',
           #'../SPHINCS_BSSN/module_grav_acc_check.f90',
           '../SPHINCS_BSSN/module_dyn_metric_on_particles.f90',
           '../SPHINCS_BSSN/module_alive_flag.f90',
           '../SPHINCS_BSSN/module_quadrup_GW_extraction.f90',
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
           '../SPHINCS_fix_metric/module_deactivate_particles.f90',
           '../SPHINCS_fix_metric/module_sorting.f90',
           '../SPHINCS_fix_metric/module_Artificial_Pressure_Method.f90',
           '../SPHINCS_fix_metric/module_piecewise_polytrope.f90']

lorene_id_executable = 'sphincs_lorene_bns.x'
lorene_id_convtest_executable = 'convergence_test.x'
write_par_eos = 'write_par_eos.x'

#+ lorene_sources + lorene_sources_bin_ns

Program(lorene_id_executable, sources + ['src/sphincs_lorene_bns.f90'])
Program(lorene_id_convtest_executable, sources + ['src/convergence_test.f90'])
Program(write_par_eos, sources + ['src/write_par_eos.f90'])
#
Decider('MD5-timestamp')
