Project: SPHINCS_ID
Version: 1.0
Project_Bitbucket: https://bitbucket.org/ftorsello/sphincs_repository_ft/src
Summary: Documentation of \(\mathrm{SPHINCS\_ID}\) <br /><br /> ![SPHINCS_ID](|media|/binary.PNG){: style="text-align: center", width=100%} <br /> <font size="2"> **Figure caption:** Projection of the SPH particles on the \(xy\) plane, for a binary neutron star system with gravitational masses of \(1.2M_\odot\) (left) and \(1.8M_\odot\) (right), with equation of state 'BBKF(DD2-SF) quark-hadron model RDF 1.8 (with electrons)' from the \(\texttt{CompOSE}\) database. The solution was computed with \(\texttt{LORENE}\) and the \(\texttt{CompOSE}\) software, and the particles were placed with \(\texttt{SPHINCS_ID}\) using the Artificial Pressure Method. The color bar shows the baryon mass density. The plot was made with \(\texttt{SPLASH}\) and \(\texttt{GIMP}\). </font>
Author: Francesco Torsello
display: private
         protected
         public
docmark_alt: #
predocmark: >
predocmark_alt: &
source: true
graph: true
sort: alpha
src_dir: ./src
media_dir: ./media
exclude_dir: ./src/prototypes
exclude: submodule_bns_lorene_bindings.f90
         submodule_diffstar_lorene_bindings.f90
         submodule_sph_particles_redistribute_nu.f90
output_dir: ./doc
page_dir: ./doc-pages
creation_date: %Y-%m-%d %H:%M:%S
print_creation_date: true
proc_internals: true
license: gfdl
github: https://github.com/francescotorsello
linkedin: https://www.linkedin.com/in/francescotorsello
alias: sphincsid = \(\texttt{SPHINCS_ID}\)
       sphincslorene = \(\texttt{SPHINCS_LORENE}\)
       lorene = \(\texttt{LORENE}\)
       compose = \(\mathrm{CompOSE}\)
       sphincsbssn = \(\texttt{SPHINCS_BSSN}\)
       sphincs = \(\texttt{SPHINCS}\)
       sph = \(\mathrm{SPH}\)
       id = \(\mathrm{ID}\)
       bssnok = \(\mathrm{BSSNOK}\)
       bssn = \(\mathrm{BSSNOK}\)
       fuka = \(\texttt{FUKA}\)
       binns = \(\texttt{Bin_NS}\)
       etrotdiff = \(\texttt{Et_rot_diff}\)
       etdiffrot = \(\texttt{Et_diffrot}\)
       eos = \(\mathrm{EOS}\)
       ee = Einstein equations

#### **S**moothed **P**article **H**ydrodynamics **IN** **C**urved **S**pacetime &mdash; **I**nitial **D**ata builder
___

SPHINCS_ID is a modular, object-oriented, OMP parallelized Fortran 2018 code to produce initial data to be evolved in time with the General Relativistic, Lagrangian Hydrodynamics, Fortran 2018 code SPHINCS_BSSN ([1][1]{:target="_blank"}), and the Newtonian, Lagrangian Hydrodynamics, Fortran code MAGMA2 ([4][4]{:target="_blank"}).

Currently, it produces initial data for binary neutron star mergers and differentially rotating stars, using the data computed by the solvers within the C++ library LORENE ([2][2]{:target="_blank"},[3][3]{:target="_blank"}).

Presently, SPHINCS_ID does not solve any equations for the initial data, but acts as an interface between an initial data solver and SPHINCS_BSSN or MAGMA2. It reads the data computed by the solver and produces the SPH and BSSN ID to be read and evolved in time with SPHINCS_BSSN or MAGMA2.

The modular and hierarchical structure of the code makes it easy to extend it to be able to set up initial data for other types of physical systems and other formulations of the Einstein equations. The code is currently under heavy development.

SPHINCS_ID needs SPHINCS_BSSN to be compiled.

---

Copyright (C) 2022 Francesco Torsello

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

[1]: <https://iopscience.iop.org/article/10.1088/1361-6382/abee65>
[2]: <https://lorene.obspm.fr/>
[3]: <https://arxiv.org/abs/gr-qc/0007028>
[4]: <https://academic.oup.com/mnras/article/498/3/4230/5897370>
[5]: <https://www.gnu.org/licenses/gpl-3.0.en.html>
___
