---
title: "The General Molecular Simulation Object (GMSO): A Data Structure for the Molecular Simulation Design Framework (MoSDeF)"

tags:
- python
- molecular-simulations
- data-structure
- MoSDeF

authors:
- name: Co D. Quach
  orcid: 0000-0002-1255-4161
  equal-contrib: true
  affiliation: "1, 2"
- name: Nicholas C. Craven
  orcid: 0000-0002-4607-4377
  affiliation: "2, 3"
- name: Umesh Timalsina
  orcid: 0000-0002-5430-3993
  affiliation: "4"
- name: Justin B. Gilmer
  orcid: 0000-0002-6915-5591
  affiliation: "2, 3"
- name: Matthew W. Thompson
  orcid: 0000-0002-1460-3983
  affiliation: "1, 2"
- name: Alexander Yang
  affiliation: "1, 2"
- name: Ray A. Matsumoto
  orcid: 0000-0002-9124-3512
  affiliation: "1, 2"
- name: Parashara Shamaprasad
  affiliation: "1, 2"
- name: Chris Jones
  affiliation: "5"
- name: Ryan S. DeFever
  affiliation: "6"
- name: Brad Crawford
  orcid: 0000-0003-0638-7333
  affiliation: "7, 8"
- name: Christopher R. Iacovella
  orcid: 0000-0003-0557-0427
  affiliation: "1, 2"
- name: Jeffrey Potoff
  affiliation: "7"
  orcid: "0000-0002-4421-8787"
-name: Eric Jankowski
  affiliation: "5"
  orcid: "0000-0002-3267-1410"
- name: Edward J. Maginn
  affiliation: "6"
- name: Clare McCabe
  orcid: 0000-0002-8552-9135
  affiliation: "1, 2, 9"
- name: Peter T. Cummings
  orcid: 0000-0002-9766-2216
  affiliation: "1, 2, 9"



affiliations:
- name: Department of Chemical and Biological Engineering, Vanderbilt University, Nashville, TN, USA
  index: 1
- name: Multiscale Modeling and Simulation (MuMS) Center, Vanderbilt University, Nashville, TN, USA
  index: 2
- name: Interdisciplinary Material Science Program, Vanderbilt University, Nashville, TN, USA
  index: 3
- name: Institute for Software Integrated Systems (ISIS), Vanderbilt University, Nashville, TN, USA
  index: 4
- name: Micron School of Materials Science and Engineering, Boise State University, Boise, ID, USA
  index: 5
- name: Department of Chemical and Biomolecular Engineering, University of Notre Dame, Notre Dame, IN, USA
  index: 6
- name: Atomfold LLC, PA, USA
  index: 7
- name: Department of Chemical Engineering, Wayne State University, Detroit, MI, USA
  index: 8
- name: School of Engineering and Physical Sciences, Heriot-Watt University, Edinburgh, Scotland, U.K
  index: 9


dates: 2 January, 2024

bibliography: paper.bib


# Summary
The General Molecular Simulation Object, or GMSO, stands as an open-source Python data structure, offering a versatile and expandable framework for handling chemical and biomolecular topologies. This library is an integral component of the Molecular Simulation Design Framework (MoSDeF), dedicated to streamlining the creation, parameterization, and representation of systems for molecule simulations. The GMSO library serves as a dynamic repository for storing chemical/biomolecular structures, encompassing metadata, coordinates, and interaction potentials. Moreover, the library includes routines for exporting stored structures into various file formats, facilitating compatibility with other software for visualization (e.g., VMD[@humphrey1996vmd] and OVITO[@]) or conducting molecular simulations (e.g., GROMACS [@abraham2015gromacs], LAMMPS[@thompson2022lammps], GOMC[@nejahi2021update]).


# Statement of need

The Molecular Simulation Design Framework (MoSDeF) is a software suite tailored to facilitate the initialization of chemical and biomolecular systems for computational simulations [@cummings2021opena]. These tools were developed to specifically address a critical aspect of the (ir)reproducibility issue within the molecular simulation community — namely, the insufficient documentation of the structure preparation process and force field parameter implementation [@thompson2020towards]. The initialization step, often performed through Graphical User Interfaces (GUI) or via the use of ad-hoc, unpublished, and unreviewed code, poses the risk of introducing irreproducible and untraceable errors[@baker2016reproducibility]. By providing general-purposed and standardized tools that build and parameterized molecular systems for molecular simulations, directly support  various molecular dynamics (MD) and Monte Carlo (MC) engines, MoSDeF aims to trivialize the describing and disseminating such processes without creating extra burdens for computational simulation researchers.[@cummings2021opena].


The initialization of chemical/biomolecular systems comprises of three key steps:
1.  Constructing structures: encompassing loading or creating molecules/structures that mirrors the phenomena under investigation.
2. Parameterizing: interactional parameters are assigned to all particles and connections within the structures.
3. Storing Structures and Output Generation: parameterized structures are stored, and output is generated in file formats compatible with various simulation software, including GROMACS, LAMMPS, HOOMD-Blue, Cassandra, and GOMC[@abraham2015gromacs; @thompson2022lammps; @anderson2010hoomd, @shah2017cassandra, @nejahi2021update].
Each of these steps necessitates distinct routines, and as such, is addressed by a series of specialized libraries — specifically, mBuild [@klein2016hierarchical], Foyer [@klein2019formalizing], and GMSO, which is introduced in this work. mBuild functions as a molecular builder, equipped with extensive utilities for creating, loading, and manipulating positions of atoms and molecules, along with managing their connectivity through bonds[@klein2016hierarchical]. Foyer assumes the role of parameterizing for the created structures, involving the identification and assignment of interaction parameters to each atom or group of atoms[@klein2019formalizing] and their associated connections (e.g., bonds, angles, and dihedrals). This process entails matching the connectivity (bond graph) of the provided structure with the SMARTS grammar of the corresponding atom type, defining the interactional parameters[@klein2019formalizing]. The use of a graph matching method, departing from the traditional approach of matching via atom indices, allows for a more flexible parameterization. This feature proves particularly advantageous in the study of functionalized polymers, whose structures consistently deviate slightly from the standard polymer[@summers2020mosdef, @quach2022high]. These utilities have been applied in various projects to explore a wide range of structures and applications [@thompson2019scalable; @summers2020mosdef; @quach2022high; @ma2022dynamics], as well as integrated into other scientific libraries [@albooyeh2023flowermd; @defever2021mosdef; @crawford2023mosdefgomc].


The parameterization step introduces additional information, requiring a more sophisticated data structure for representation. Beyond the initial details concerning positions and connectivity established during system construction, the new structure incorporates supplementary metadata and interaction parameters. These encompass not only atoms and bonds but also extend to angles, dihedrals, and improper dihedrals, along with interaction parameters associated with each of these objects. This requires the development of data structures that are capable of:
- Supporting a variety of models
- Providing flexibility for exotic potentials and unit systems
- Being compatible with existing community tools
- Being extensible (to support new simulation models/engines/workflows)
Currently, existing data structures, such as ParmEd and OpenMM[@shirts2016lessons; @eastmann2017openmm], fulfill many functionalities and are widely adopted [add citations]. However, their underlying structures are tailored to specific subsets of simulation workflows and ecosystems, as well as force field equation forms, sacrificing generality and broad applicability. This limitation includes hard-coding and assumptions about potential expressions and units. They lack the generality that MoSDeF and its users seek, such as the ability to define and store arbitrary potential expressions or unit systems. Integrating these new features into existing software, unfortunately, would require a major overhaul, potentially impacting existing simulation workflows and is not appealing to current project stakeholders.

Hence, we develop the General Molecular Simulation Object (GMSO) library, which is a lightweight, extensible data structure encapsulating chemical/biomolecular systems and their associated interaction parameters, i.e., force fields, to cater to MoSDeF ecosystem. The library is designed to accommodate a wide range of chemical/biomolecular models, offering the capability to support arbitrary potential expressions and unit systems. Generalizing these potential (force field) expressions also allows users to enter the force field in its native form and units, minimizing user error when setting up the force field file while providing the ability to easily auto-convert the potential form and units to the molecular engine's required form. GMSO satisfies the broader community's need for a general, extensible, and reproducible method of setting up molecular simulations. In addition to core data classes, the library includes routines for interacting/converting to and from other ecosystems, including ParmEd and OpenMM, enhancing interoperability without reinventing functionalities. GMSO supports output to multiple molecular simulation engine-specific file formats, currently including GROMACS, LAMMPS, HOOMD-Blue, GOMC, Cassandra, with plans for future expansion. When integrated with other MoSDeF software, GMSO facilitates large-scale automated molecular screening for diverse molecules/structures or state points, applying correct force field parameters using SMARTS strings[@klein2016hierarchical, @klein2019formalizing]. Combining MoSDeF software in a unified workflow using a manager like Signac [@adorf2018simple] enables large-scale molecular screenings critical for developing new materials, chemicals, and drugs[@quach2022high, @thompson2019scalable].

# Acknowledgements
This research was partially supported by the National Science Foundation OAC-1835713 and OAC-1835874. Atomfold also donated research and development time and computational resources for this research and software.
