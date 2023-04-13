home_sphincs='../../../SPHINCS_20230411'
#env['F90PATH'] = [
#  '.',
#  home_sphincs + '/BSSN',
#  home_sphincs + '/sphincs_repository/SPHINCS_fix_metric',
#  home_sphincs + '/sphincs_repository/SPHINCS_BSSN'
#]
env['F90PATH'] = ['.']

mesh_dir         = [home_sphincs + '/BSSN/']
sphincs_bssn_dir = [home_sphincs + '/sphincs_repository/SPHINCS_BSSN/']
sphincs_fm_dir   = [home_sphincs + '/sphincs_repository/SPHINCS_fix_metric/']

sources_mesh = [
  'module_mesh_refinement.f90',
  'submodule_prolongation.f90',
  'submodule_restriction.f90',
  'submodule_utilities.f90',
  'submodule_parameter.f90',
  'submodule_allocation.f90',
  'submodule_output.f90',
  'submodule_refine.f90',
  'submodule_RK4.f90',
  #
  'module_ADM_refine.f90',
  'module_BSSN_refine.f90',
  'module_BSSN_parameters.f90',
  'module_GravityAcceleration_refine.f90',
  #
  'module_McLachlan_refine.f90',
  'submodule_ML_evolution.f90',
  'submodule_ML_init.f90',
  'submodule_ML_utilities.f90',
  'submodule_ML_constraints.f90',
  'submodule_ML_boundary.f90',
  'submodule_ML_allocation.f90',
  'submodule_ML_output.f90',
  #
  'module_TOV_refine.f90',
  'module_Tmunu_refine.f90',
  'module_map_metric_2_particles_refine.f90',
  'module_evolution_parameters.f90',
  'module_Numerics.f90',
  'module_Hermite_refine.f90',
  'Hermite_3D_Weights.f',
  'Hermite5_3D_Weights.f',
  'module_WENO_refine.f90',
  'ML_BSSN_NV_EvolutionInteriorSplitBy1.cc',
  'ML_BSSN_NV_EvolutionInteriorSplitBy2.cc',
  'ML_BSSN_NV_EvolutionInteriorSplitBy3.cc',
  'ML_BSSN_NV_InitialADMBase1Everywhere.cc',
  'ML_BSSN_NV_InitialADMBase2Interior.cc',
  'ML_BSSN_NV_MatterAccelerationInterior.cc',
  'ML_BSSN_NV_EnforceEverywhere.cc',
  'ML_BSSN_NV_ADMDerivativesInterior.cc',
  'ML_BSSN_NV_ConstraintsInterior.cc',
  'ML_BSSN_NV_ConstraintsInterior_mod.cc',
  'ML_BSSN_NV_RicciInterior.cc',
  #
  'module_Extract_Mass.f90',
  'submodule__Extract_Mass_impl.f90',
  #
  #
  'module_Collapse.f90',
  'submodule_Collapse.f90',
  #
  'fermi_step_f.f90',
  'bssn_parameters.cc'
]

sources_sphincs_bssn = [
  'module_CP_distribute.f90',
  'module_MLS.f90',
  'module_SelfRegularization.f90',
  'module_particle_mesh_hash_grid.f90',
  'module_map_particles_2_grid.f90',
  'module_integ_hydro_BSSN.f90',
  'module_gradient.f90',
  'module_timing.f90',
  'module_debug_output.f90',
  'module_dyn_metric_on_particles.f90',
  'module_alive_flag.f90',
  'module_quadrup_GW_extraction.f90',
  'module_LRE.f90'
]

sources_sphincs_fm = [
  'module_files.f90',
  'module_constants.f90',
  'module_units.f90',
  'module_SPHINCS_options.f90',
  'module_boundaries.f90',
  'module_Numerical_Recipes.f90',
  'module_matrix.f90',
  'module_3plus1.f90',
  'module_kernel_table.f90',
  'module_RCB_tree_neighbours.f90',
  'module_SPH_variables.f90',
  'module_analyze.f90',
  'submodule_analyze_adm_variables.f90',
  'module_EOS.f90',
  'module_set_h.f90',
  'module_tensor.f90',
  'module_BH_params.f90',
  'module_input_output.f90',
  'module_recovery_via_pressure.f90',
  'module_SPHINCS_SPH_v2.f90',
  'module_SPHINCS_version.f90',
  'module_SPHINCS_logbook.f90',
  'module_relativistic_Eigenvalues.f90',
  'module_reconstruction.f90',
  'module_deactivate_particles.f90',
  'module_sorting.f90',
  'module_Finite_Diff_1D.f90',
  'module_Artificial_Pressure_Method.f90',
  'module_piecewise_polytrope.f90',
  'module_particle_identity.f90',
  'submodule_particle_identity.f90',
]

for i in range(len(sources_mesh)):
 sources_mesh[i]           = mesh_dir[0]          + sources_mesh[i]

for i in range(len(sources_sphincs_bssn)):
 sources_sphincs_bssn[i]   = sphincs_bssn_dir[0]  + sources_sphincs_bssn[i]

for i in range(len(sources_sphincs_fm)):
 sources_sphincs_fm[i]     = sphincs_fm_dir[0]    + sources_sphincs_fm[i]
