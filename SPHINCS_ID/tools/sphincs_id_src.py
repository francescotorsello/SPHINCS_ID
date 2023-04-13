################################################################################

env['F90PATH'] = ['.']

module_id_base_dir = ['id_base/']
module_id_base = [
  'module_id_base.f90',
  'submodule_id_base_initialization.f90',
  'submodule_id_base_access.f90',
  'submodule_id_base_mass_profile.f90',
  'submodule_id_base_length_scale.f90'
]

module_sph_particles_dir = ['sph_particles/']
module_sph_particles = [
  'module_sph_particles.f90',
  'submodule_sph_particles_constructor_std.f90',
  'submodule_sph_particles_constructor_bin.f90',
  'submodule_sph_particles_memory.f90',
  'submodule_sph_particles_lattices.f90',
  'submodule_sph_particles_ellipsoidal_surfaces.f90',
  'submodule_sph_particles_sph_variables.f90',
  'submodule_sph_particles_access.f90',
  'submodule_sph_particles_io.f90',
  'submodule_sph_particles_compose.f90',
  #'submodule_sph_particles_redistribute_nu.f90',
  'submodule_sph_particles_apm.f90',
  'submodule_sph_particles_recovery.f90',
  'submodule_sph_particles_handle_positions.f90',
  #'submodule_sph_particles_adm_variables.f90',
  'submodule_sph_particles_quality_indicators.f90'
]

module_standard_tpo_formulation_dir = ['standard_tpo_formulation/']
module_standard_tpo_formulation = [
  'module_standard_tpo_formulation.f90',
  'submodule_standard_tpo_formulation_standard_tpo_variables.f90',
  'submodule_standard_tpo_formulation_access.f90',
  'submodule_standard_tpo_formulation_io.f90',
  'submodule_standard_tpo_formulation_recovery.f90',
  'submodule_standard_tpo_formulation_analysis.f90',
  'submodule_standard_tpo_formulation_sph_adm_variables.f90'
]

module_bssn_formulation_dir = \
  ['standard_tpo_formulation/bssn_formulation/']
module_bssn_formulation = [
  'module_bssn_formulation.f90',
  'submodule_bssn_formulation_constructor.f90',
  'submodule_bssn_formulation_memory.f90',
  'submodule_bssn_formulation_io.f90',
  'submodule_bssn_formulation_bssn_variables.f90',
  'submodule_bssn_formulation_constraints.f90',
  'submodule_bssn_formulation_ricci.f90',
  'submodule_bssn_formulation_landau_lifshitz.f90'
]

module_bns_base_dir = ['id_base/bns_base/']
module_bns_base = [
  'module_bns_base.f90',
  'submodule_bns_base_access.f90',
  'submodule_bns_base_io.f90',
  'submodule_bns_base_geometry.f90'
]

module_bns_lorene_dir = ['id_base/bns_base/bns_lorene/']
module_bns_lorene = [
  'module_bns_lorene.f90',
  'submodule_bns_lorene_constructor.f90',
  'submodule_bns_lorene_read.f90',
  'submodule_bns_lorene_memory.f90',
  'submodule_bns_lorene_access.f90',
  'submodule_bns_lorene_properties.f90',
  'submodule_bns_lorene_io.f90',
  'submodule_bns_lorene_finalize_id.f90'
]

module_diffstar_base_dir = ['id_base/diffstar_base/']
module_diffstar_base = [
  'module_diffstar_base.f90',
  'submodule_diffstar_base_access.f90',
  'submodule_diffstar_base_io.f90'
]

module_diffstar_lorene_dir = ['id_base/diffstar_base/diffstar_lorene/']
module_diffstar_lorene = [
  'module_diffstar_lorene.f90',
  'submodule_diffstar_lorene_constructor.f90',
  'submodule_diffstar_lorene_read.f90',
  'submodule_diffstar_lorene_memory.f90',
  'submodule_diffstar_lorene_access.f90',
  'submodule_diffstar_lorene_properties.f90',
  'submodule_diffstar_lorene_io.f90'
]

module_bns_fuka_dir = ['id_base/bns_base/bns_fuka/']
module_bns_fuka = [
  'module_bns_fuka.f90',
  'submodule_bns_fuka_constructor.f90',
  'submodule_bns_fuka_read.f90',
  'submodule_bns_fuka_interpolate.f90',
  'submodule_bns_fuka_memory.f90',
  'submodule_bns_fuka_access.f90',
  'submodule_bns_fuka_properties.f90',
  'submodule_bns_fuka_io.f90'
]

module_ejecta_generic_dir = ['id_base/ejecta_generic/']
module_ejecta_generic = [
  'module_ejecta_generic.f90',
  'submodule_ejecta_generic_constructor.f90',
  'submodule_ejecta_generic_access.f90',
  'submodule_ejecta_generic_memory.f90',
  'submodule_ejecta_generic_io.f90',
  'submodule_ejecta_generic_interpolate.f90'
]

module_sphincs_id_full_dir = ['sphincs_id_full/']
module_sphincs_id_full = [
  'module_sphincs_id_full.f90'
]

module_sphincs_id_lorene_dir = ['sphincs_id_lorene/']
module_sphincs_id_lorene = [
  'module_sphincs_id_lorene.f90'
]

module_sphincs_id_fuka_dir = ['sphincs_id_fuka/']
module_sphincs_id_fuka = [
  'module_sphincs_id_fuka.f90'
]

module_sphincs_id_interpolate_dir = ['sphincs_id_interpolate/']
module_sphincs_id_interpolate = [
  'module_sphincs_id_interpolate.f90'
]

module_cauchy_convergence_test_dir = ['cauchy_convergence_test/']
module_cauchy_convergence_test = [
  'module_cauchy_convergence_test.f90',
  'submodule_cauchy_convergence_test_shared_grid.f90',
  'submodule_cauchy_convergence_test_perform_test.f90'
]

module_lorentz_group_dir = ['lorentz_group/']
module_lorentz_group = [
  'module_lorentz_group.f90',
  'submodule_lorentz_group_constructors.f90',
  'submodule_lorentz_group_actions.f90'
]

module_wd_eos_dir = ['wd_eos/']
module_wd_eos = [
  'module_wd_eos.f90'
]

module_tabulated_eos_dir = ['tabulated_eos/']
module_tabulated_eos = [
  'module_tabulated_eos.f90'
]

module_utility_dir = ['utility/']
module_utility = [
  'module_utility.f90'
]

for i in range(len(module_id_base)):
  module_id_base[i] = module_id_base_dir[0] + module_id_base[i]

for i in range(len(module_sph_particles)):
  module_sph_particles[i] = module_sph_particles_dir[0] \
                          + module_sph_particles[i]

for i in range(len(module_standard_tpo_formulation)):
  module_standard_tpo_formulation[i] = module_standard_tpo_formulation_dir[0] \
                                     + module_standard_tpo_formulation[i]

for i in range(len(module_bssn_formulation)):
  module_bssn_formulation[i] = module_bssn_formulation_dir[0] \
                               + module_bssn_formulation[i]

for i in range(len(module_bns_base)):
  module_bns_base[i] = module_bns_base_dir[0] + module_bns_base[i]

for i in range(len(module_bns_lorene)):
  module_bns_lorene[i] = module_bns_lorene_dir[0] + module_bns_lorene[i]

for i in range(len(module_diffstar_base)):
  module_diffstar_base[i] = module_diffstar_base_dir[0] \
                          + module_diffstar_base[i]

for i in range(len(module_diffstar_lorene)):
  module_diffstar_lorene[i] = module_diffstar_lorene_dir[0] \
                            + module_diffstar_lorene[i]

for i in range(len(module_bns_fuka)):
  module_bns_fuka[i] = module_bns_fuka_dir[0] + module_bns_fuka[i]

for i in range(len(module_ejecta_generic)):
  module_ejecta_generic[i] = module_ejecta_generic_dir[0] \
                           + module_ejecta_generic[i]

for i in range(len(module_sphincs_id_full)):
  module_sphincs_id_full[i] = module_sphincs_id_full_dir[0] \
                            + module_sphincs_id_full[i]

for i in range(len(module_sphincs_id_lorene)):
  module_sphincs_id_lorene[i] = module_sphincs_id_lorene_dir[0] \
                              + module_sphincs_id_lorene[i]

for i in range(len(module_sphincs_id_fuka)):
  module_sphincs_id_fuka[i] = module_sphincs_id_fuka_dir[0] \
                            + module_sphincs_id_fuka[i]

for i in range(len(module_sphincs_id_interpolate)):
  module_sphincs_id_interpolate[i] = module_sphincs_id_interpolate_dir[0] \
                                   + module_sphincs_id_interpolate[i]

for i in range(len(module_cauchy_convergence_test)):
  module_cauchy_convergence_test[i] = module_cauchy_convergence_test_dir[0] \
                                    + module_cauchy_convergence_test[i]

for i in range(len(module_lorentz_group)):
  module_lorentz_group[i] = module_lorentz_group_dir[0] \
                          + module_lorentz_group[i]

for i in range(len(module_wd_eos)):
  module_wd_eos[i] = module_wd_eos_dir[0] + module_wd_eos[i]

for i in range(len(module_tabulated_eos)):
  module_tabulated_eos[i] = module_tabulated_eos_dir[0] \
                          + module_tabulated_eos[i]

for i in range(len(module_utility)):
  module_utility[i] = module_utility_dir[0] + module_utility[i]

sources_base = module_id_base + module_bns_base + module_diffstar_base \
             + module_sph_particles + module_standard_tpo_formulation \
             + module_bssn_formulation + module_cauchy_convergence_test \
             + module_lorentz_group + module_wd_eos + module_tabulated_eos \
             + module_utility

#all_modules = [ module_id_base, module_utility, module_sphincs_id_interpolate ]

################################################################################
