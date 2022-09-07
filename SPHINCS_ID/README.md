### **S**moothed **P**article **H**ydrodynamics **IN** **C**urved **S**pacetime &mdash; **I**nitial **D**ata builder
___

SPHINCS_ID is a modular, object-oriented, OMP parallelized Fortran 2018 code to produce initial data to be evolved in time with the General Relativistic, Lagrangian Hydrodynamics, Fortran 2018 code SPHINCS_BSSN ([1][1]{:target="_blank"}), and the Newtonian, Lagrangian Hydrodynamics, Fortran code MAGMA2 ([2][2]{:target="_blank"}).

Presently, SPHINCS_ID does not solve any equations for the initial data, but acts as an interface between an initial data solver and SPHINCS_BSSN or MAGMA2. It reads the data computed by the solver and produces the SPH and BSSN ID to be read and evolved in time with SPHINCS_BSSN or MAGMA2.

Currently, it produces initial data for:

  - binary neutron star mergers and differentially rotating stars, using the data computed by the solvers within the C++ library LORENE ([3][3]{:target="_blank"},[4][4]{:target="_blank"})
  - data on a Cartesian, uniform grid, representing a generic physical system

The modular and hierarchical structure of the code makes it easy to extend it to be able to set up initial data for other types of physical systems and other formulations of the Einstein equations. The code is currently under heavy development.

SPHINCS_ID needs SPHINCS_BSSN to be compiled.

The User Manual for SPHINCS_ID can be found in doc-pages/user_manual.pdf

---

[1]: <https://iopscience.iop.org/article/10.1088/1361-6382/abee65>
[2]: <https://academic.oup.com/mnras/article/498/3/4230/5897370>
[3]: <https://lorene.obspm.fr/>
[4]: <https://arxiv.org/abs/gr-qc/0007028>
[5]: <https://www.gnu.org/licenses/gpl-3.0.en.html>