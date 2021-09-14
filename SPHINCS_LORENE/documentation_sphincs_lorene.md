Project: SPHINCS_LORENE
Version: 2.0
Project_Bitbucket: https://bitbucket.org/ftorsello/sphincs_repository_ft/src
Summary: Documentation of the FORTRAN 2018 code SPHINCS_LORENE
Author: Francesco Torsello
display: private
         protected
         public
docmark_alt: #
predocmark: >
predocmark_alt: &
source: true
graph: false
sort: alpha
src_dir: ./src
exclude_dir: ./src/prototypes
exclude: submodule_bns_bindings.f90
output_dir: ./doc
creation_date: %Y-%m-%d %H:%M:%S
print_creation_date: true
proc_internals: true
license: gfdl
github: https://github.com/francescotorsello
linkedin: https://www.linkedin.com/in/francescotorsello

SPHINCS_LORENE is a modular, object-oriented, OMP parallelized FORTRAN 2018 code to produce binary neutron stars initial data to be evolved in time with the FORTRAN 2018 code SPHINCS_BSSN ([1][1]{:target="_blank"}), using the C++ code LORENE ([2][2]{:target="_blank"},[3][3]{:target="_blank"}).

SPHINCS_LORENE acts as an interface between LORENE and SPHINCS_BSSN; it reads the spectral ID produced by LORENE and produces the SPH and BSSN ID to be read and evolved in time with SPHINCS_BSSN.

[1]: <https://iopscience.iop.org/article/10.1088/1361-6382/abee65>
[2]: <https://lorene.obspm.fr/>
[3]: <https://arxiv.org/abs/gr-qc/0007028>
