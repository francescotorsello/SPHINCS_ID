Project: SPHINCS_ID
Version: 1.8
linkedin: https://www.linkedin.com/in/francescotorsello
Summary: Documentation of \(\mathrm{SPHINCS\_ID}\) <br /><br /> ![SPHINCS_ID](|media|/binary.PNG){: style="text-align: center", width=100%} <br /> <font size="2"> **Figure caption:** Projection of the SPH particles on the \(xy\) plane, for a binary neutron star system with gravitational masses of \(1.2M_\odot\) (left) and \(1.8M_\odot\) (right), with equation of state 'BBKF(DD2-SF) quark-hadron model RDF 1.8 (with electrons)' from the \(\texttt{CompOSE}\) database. The solution was computed with \(\texttt{LORENE}\) and the \(\texttt{CompOSE}\) software, and the particles were placed with \(\texttt{SPHINCS_ID}\) using the Artificial Pressure Method. The color bar shows the baryon mass density. The plot was made with \(\texttt{SPLASH}\) [1][8]{:target="_blank"} and \(\texttt{GIMP}\). </font>
Author: Francesco Torsello
display: private
         protected
         public
docmark_alt: #
predocmark: >
predocmark_alt: &
source: true
incl_src: false
graph: true
sort: alpha
src_dir: ./src
media_dir: ./doc-media
exclude: ./src/sph_particles/submodule_sph_particles_redistribute_nu.f90
output_dir: ./doc
page_dir: ./doc-pages
creation_date: %Y-%m-%d %H:%M:%S
print_creation_date: true
proc_internals: true
favicon: ./doc-media/favicon.png
alias: sphincsid = \(\texttt{SPHINCS_ID}\)
       lorene = \(\texttt{LORENE}\)
       compose = \(\mathrm{CompOSE}\)
       sphincsbssn = \(\texttt{SPHINCS_BSSN}\)
       sphincs = \(\texttt{SPHINCS}\)
       sph = \(\mathrm{SPH}\)
       id = \(\mathrm{ID}\)
       bssnok = \(\mathrm{BSSNOK}\)
       bssn = \(\mathrm{BSSNOK}\)
       fuka = \(\texttt{FUKA}\)
       kadath = \(\texttt{Kadath}\)
       binns = \(\texttt{Bin_NS}\)
       etrotdiff = \(\texttt{Et_rot_diff}\)
       etdiffrot = \(\texttt{Et_diffrot}\)
       bnsexp = \(\texttt{bns_export}\)
       eos = \(\mathrm{EOS}\)
       tov = \(\mathrm{TOV}\)
       ee = Einstein equations
       bns = \(\mathrm{BNS}\)
       drs = \(\mathrm{DRS}\)
       nb = n_\mathrm{b}
parallel: 80
<!---graph_dir: doc-graphs--->

### **S**moothed **P**article **H**ydrodynamics **IN** **C**urved **S**pacetime &mdash; **I**nitial **D**ata builder
___

SPHINCS_ID is a modular, object-oriented, OMP parallelized Fortran 2018 code to produce initial data to be evolved in time with the General Relativistic, Lagrangian Hydrodynamics, Fortran 2018 code SPHINCS_BSSN ([2][1]{:target="_blank"}), and the Newtonian, Lagrangian Hydrodynamics, Fortran code MAGMA2 ([3][2]{:target="_blank"}).

Presently, SPHINCS_ID does not solve any equations for the initial data, but acts as an interface between an initial data solver and SPHINCS_BSSN or MAGMA2. It reads the data computed by the solver and produces the SPH and BSSN ID to be read and evolved in time with SPHINCS_BSSN or MAGMA2.

Currently, it produces initial data for:

  - binary systems of neutron stars and differentially rotating stars, using the data computed by the solvers within the C++ library LORENE ([4][3]{:target="_blank"},[5][4]{:target="_blank"})
  - binary systems of neutron stars, using the data computed by the FUKA solvers within the C++ library Kadath ([6][5]{:target="_blank"},[7][6]{:target="_blank"})
  - data on a Cartesian, uniform grid, representing a generic physical system
  - Newtonian binary systems of neutron stars and white dwarfs, using the data computed by the TOV solver within SPHINCS_BSSN; in other words, two TOV stars are placed on an orbit given by the Newtonian 2-body problem

The modular and hierarchical structure of the code makes it easy to extend it to be able to set up initial data for other types of physical systems and other formulations of the Einstein equations.

SPHINCS_ID needs SPHINCS_BSSN to be compiled.

The User Manual for SPHINCS_ID is <a href="page/SPHINCS_ID-User_Manual.pdf" target="_blank">SPHINCS_ID-User_Manual.pdf</a>

---

#### Acknowledgements

It is a pleasure to thank Peter Diener and Stephan Rosswog for all the help and support they gave me during the development of SPHINCS_ID.

---

Copyright (C) 2020-2023 Francesco Torsello.

SPHINCS_ID is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SPHINCS_ID. If not, see <https://www.gnu.org/licenses/>.
___

Copyright (C) 2020-2023 Francesco Torsello.

Permission is granted to copy, distribute and/or modify this documentation
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
A copy of the license is included in the section entitled "License", reachable by clicking "About" at the top of the webpage, or at <https://www.gnu.org/licenses/fdl-1.3.html/>.
---

[1]: <https://iopscience.iop.org/article/10.1088/1361-6382/abee65>
[2]: <https://academic.oup.com/mnras/article/498/3/4230/5897370>
[3]: <https://lorene.obspm.fr/>
[4]: <https://arxiv.org/abs/gr-qc/0007028>
[5]: <https://kadath.obspm.fr/fuka/>
[6]: <https://arxiv.org/abs/2103.09911>
[7]: <https://www.gnu.org/licenses/gpl-3.0.en.html>
[8]: <https://doi.org/10.1071/AS07022>
