Project: SPHINCS_ID
Version: 1.0
Project_Bitbucket: https://bitbucket.org/ftorsello/sphincs_repository_ft/src
Summary: Documentation of \(\mathrm{SPHINCS\_ID}\) <br /><br /> ![SPHINCS_ID](|media|/binary-tr.PNG){: style="text-align: center", width=100%} <font size="2"> **Image caption:** Projection of the SPH particles on the \(xy\) plane, for a binary neutron star system with gravitational masses of \(1.2M_\odot\) (left) and \(1.8M_\odot\) (right), with equation of state DD2-SF RDF 1.8 from the CompOSE database. The solution was computed with \(\texttt{LORENE}\), and the particles placed with \(\texttt{SPHINCS_ID}\) using the Artificial Pressure Method. The color bar shows the baryon mass density. The plot was made with \(\texttt{SPLASH}\) and \(\texttt{GIMP}\). </font>
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
         submodule_particles_redistribute_nu.f90
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
       fuka = \(\texttt{FUKA}\)
       binns = \(\texttt{Bin_NS}\)
       etrotdiff = \(\texttt{Et_rot_diff}\)
       etdiffrot = \(\texttt{Et_diffrot}\)
       eos = \(\mathrm{EOS}\)

#### **S**moothed **P**article **H**ydrodynamics **IN** **C**urved **S**pacetime &mdash; **I**nitial **D**ata builder
___

SPHINCS_ID is a modular, object-oriented, OMP parallelized FORTRAN 2018 code to produce initial data to be evolved in time with the FORTRAN 2018 code SPHINCS_BSSN ([1][1]{:target="_blank"}). Currently, it produces initial data for binary neutron star mergers and differentially rotating stars, using the data computed by the solvers within the C++ library LORENE ([2][2]{:target="_blank"},[3][3]{:target="_blank"}).

Currently, SPHINCS_ID does not solve any equations for the initial data, but acts as an interface between an initial data solver and SPHINCS_BSSN. It reads the data computed by the solver and produces the SPH and BSSN ID to be read and evolved in time with SPHINCS_BSSN.

[1]: <https://iopscience.iop.org/article/10.1088/1361-6382/abee65>
[2]: <https://lorene.obspm.fr/>
[3]: <https://arxiv.org/abs/gr-qc/0007028>
___
