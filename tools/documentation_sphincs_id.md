project: SPHINCS_ID
version: 2.0
linkedin: https://www.linkedin.com/in/francescotorsello
project_bitbucket: https://bitbucket.org/sphincsid/sphincs_id
summary: Documentation of \(\mathrm{SPHINCS\_ID}\) <br /><br /> ![SPHINCS_ID](|media|/binary.png){: style="text-align: center", width=100%} <br /> <font size="2"> **Figure caption:** Projection of \(5\times 10^6\) total SPH particles on the \(xy\) plane, for a binary neutron star system with gravitational masses of \(1.1M_\odot\) (left) and \(1.5M_\odot\) (right), with equation of state 'BBKF(DD2-SF) quark-hadron model RDF 1.8 (with electrons)' from the <a href="https://compose.obspm.fr/" target="_blank">\(\texttt{CompOSE}\) service</a>. The solution was computed with <a href="https://lorene.obspm.fr/" target="_blank">\(\texttt{LORENE}\)</a> and the <a href="https://compose.obspm.fr/software" target="_blank">\(\texttt{CompOSE}\) software</a>, and the particles were placed with \(\texttt{SPHINCS_ID}\) using the Artificial Pressure Method. The color bar shows the baryon mass density. The plot was made with <a href="https://doi.org/10.1071/AS07022" target="_blank">\(\texttt{SPLASH}\)</a> and <a href="https://www.gimp.org/" target="_blank">\(\texttt{GIMP}\)</a>. </font>
author: Francesco Torsello
display: private
         protected
         public
docmark_alt: #
predocmark: >
predocmark_alt: &
source: true
incl_src: true
graph: true
sort: alpha
src_dir: ../src
media_dir: ../res/doc-media
page_dir: ../res/doc-pages
exclude: ./src/sph_particles/submodule_sph_particles_redistribute_nu.f90
output_dir: ../doc
creation_date: %Y-%m-%d %H:%M:%S
print_creation_date: true
proc_internals: true
favicon: ../res/doc-media/favicon.png
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
<!---graph_dir: ../res/doc-graphs--->

### **S**moothed **P**article **H**ydrodynamics **IN** **C**urved **S**pacetime &mdash; **I**nitial **D**ata builder
___

SPHINCS_ID is a modular, object-oriented, OMP parallelized Fortran 2018 code to produce initial data to be evolved in time with the General Relativistic, Lagrangian Hydrodynamics, Fortran 2018 code SPHINCS_BSSN ([1][1]{:target="_blank"}), and the Newtonian, Lagrangian Hydrodynamics, Fortran code MAGMA2 ([2][2]{:target="_blank"}).

Presently, SPHINCS_ID does not solve any equations for the initial data, but acts as an interface between an initial data solver and SPHINCS_BSSN or MAGMA2. It reads the data computed by the solver and produces the SPH and BSSN ID to be read and evolved in time with SPHINCS_BSSN or MAGMA2.

Currently, it produces initial data for:

  - binary systems of neutron stars and differentially rotating stars, using the data computed by the solvers within the C++ library LORENE ([3][3]{:target="_blank"},[4][4]{:target="_blank"})
  - binary systems of neutron stars, using the data computed by the FUKA solvers within the C++ library Kadath ([5][5]{:target="_blank"},[6][6]{:target="_blank"})
  - data on a Cartesian, uniform grid, representing a generic physical system
  - Newtonian binary systems of neutron stars and white dwarfs, using the data computed by the TOV solver within SPHINCS_BSSN; in other words, two TOV stars are placed on an orbit given by the Newtonian 2-body problem

The modular and hierarchical structure of the code makes it easy to extend it to be able to set up initial data for other types of physical systems and other formulations of the Einstein equations.

SPHINCS_ID links to SPHINCS_BSSN.

The User Manual for SPHINCS_ID is <a href="page/SPHINCS_ID-User_Manual.pdf" target="_blank">SPHINCS_ID-User_Manual.pdf</a>

---

#### Required citations

The algorithms implemented in SPHINCS_ID were presented in the following references, so please cite them if you use SPHINCS_ID.

- Diener, P., Rosswog, S. & Torsello, F. Simulating neutron star mergers with the Lagrangian Numerical Relativity code SPHINCS_BSSN. Eur. Phys. J. A 58, 74 (2022). [https://doi.org/10.1140/epja/s10050-022-00725-7](https://doi.org/10.1140/epja/s10050-022-00725-7)

- Rosswog, S., Torsello, F. & Diener, P., The Lagrangian numerical relativity code SPHINCS_BSSN_v1.0. Front. Appl. Math. Stat. 9 (2023). [https://doi.org/10.3389/fams.2023.1236586](https://doi.org/10.3389/fams.2023.1236586)

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
along with SPHINCS_ID. If not, see <https://www.gnu.org/licenses/gpl-3.0.html>.

---

Copyright (C) 2020-2023 Francesco Torsello.

Permission is granted to copy, distribute and/or modify this documentation
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
A copy of the license is included in the section entitled "License", reachable by clicking "About" at the top of the webpage, or at <https://www.gnu.org/licenses/fdl-1.3.html>.

---

[1]: <https://iopscience.iop.org/article/10.1088/1361-6382/abee65>
[2]: <https://academic.oup.com/mnras/article/498/3/4230/5897370>
[3]: <https://lorene.obspm.fr/>
[4]: <https://arxiv.org/abs/gr-qc/0007028>
[5]: <https://kadath.obspm.fr/fuka/>
[6]: <https://arxiv.org/abs/2103.09911>
[7]: <https://www.gnu.org/licenses/gpl-3.0.en.html>
[8]: <https://doi.org/10.1071/AS07022>
