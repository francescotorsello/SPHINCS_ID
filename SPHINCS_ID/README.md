# **S**moothed **P**article **H**ydrodynamics **IN** **C**urved **S**pacetime &mdash; **I**nitial **D**ata builder
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

The User Manual for SPHINCS_ID is res/SPHINCS_ID-User_Manual.pdf

Please read the README.md files in each directory for more details.

# Compilation

SPHINCS_ID is compiled using SCons.

Follow the instructions at [the SCons documentation](https://scons.org/doc/production/HTML/scons-user/index.html){:target="_blank"} (or any later version of it) to install SCons. Once SCons is installed, go to the root directory of SPHINCS_ID where the SConstruct file is placed; if you are compiling SPHINCS_ID on a new host, open SConstruct and set up your local environment. After that, run 'scons' and the compilation will start. See res/SPHINCS_ID-User_Manual.pdf for more details.

The compilation of SPHINCS_ID will create a directory named build and place the object files *.o inside it, following the same directory structure as in src. The *.mod files will instead be placed inside the mod directory. The executable files *.x will be placed inside programs/bin. The configuration files (or parameter files) needed to run the executables are placed inside config. See res/SPHINCS_ID-User_Manual.pdf for more details.

# Producing the documentation

The documentation of SPHINCS_ID is produced with FORD.

To install FORD, run 'pip install ford' or follow the instructions at [its GitHub repository](https://github.com/Fortran-FOSS-Programmers/ford){:target="_blank"}. One FORD is installed, go to the root directory of SPHINCS_ID and run 'ford tools/documentation_sphincs_id.md'. The documentation will be generated into the doc directory as an HTML document. After it is produced, open the file doc/index.html with any browser, to read it. See res/SPHINCS_ID-User_Manual.pdf for more details.

---

[1]: <https://iopscience.iop.org/article/10.1088/1361-6382/abee65>
[2]: <https://academic.oup.com/mnras/article/498/3/4230/5897370>
[3]: <https://lorene.obspm.fr/>
[4]: <https://arxiv.org/abs/gr-qc/0007028>
[5]: <https://kadath.obspm.fr/fuka/>
[6]: <https://arxiv.org/abs/2103.09911>
[7]: <https://www.gnu.org/licenses/gpl-3.0.en.html>
