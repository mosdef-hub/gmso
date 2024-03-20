---
title: "The General Molecular Simulation Object (GMSO): A Data Structure for the Molecular Simulation Design Framework (MoSDeF)"

tags:
- python
- molecular-simulations
- data-structure
- MoSDeF
- interoperability
- force fields

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
  orcid: 0000-0002-6196-5274
  affiliation: "5"
- name: Ryan S. DeFever
  orcid: 0000-0001-5311-6718
  affiliation: "6"
- name: Brad Crawford
  orcid: 0000-0003-0638-7333
  affiliation: "7, 8"
- name: Christopher R. Iacovella
  orcid: 0000-0003-0557-0427
  affiliation: "1, 2"
- name: Jeffrey Potoff
  affiliation: "7"
  orcid: 0000-0002-4421-8787
-name: Eric Jankowski
  affiliation: "5"
  orcid: 0000-0002-3267-1410
- name: Edward J. Maginn
  orcid: 0000-0002-6309-1347
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
- name: Atomfold, PA, USA
  index: 7
- name: Department of Chemical Engineering, Wayne State University, Detroit, MI, USA
  index: 8
- name: School of Engineering and Physical Sciences, Heriot-Watt University, Edinburgh, Scotland, U.K
  index: 9


dates: 17 March, 2024

bibliography: paper.bib


# Summary

The General Molecular Simulation Object, or GMSO, is an open-source Python package designed to supplement molecular simulation workflow. This library offers a versatile and expandable data structures crucial for storage of chemical and biomolecular topologies, along with utilities necessary for editing and outputting these systems. GMSO is a core component of the Molecular Simulation Design Framework (MoSDeF), dedicated to streamlining the creation, parameterization, and representation of systems for molecular simulations. The GMSO library serves as a dynamic repository for storing chemical/biomolecular structures, encompassing metadata, coordinates, and interaction potentials. Moreover, the library includes routines for editing and exporting stored structures into various file formats, which can be used with other software for visualization (e.g., VMD[@humphrey1996vmd] and OVITO[@]) or conducting molecular simulations (e.g., GROMACS[@abraham2015gromacs], LAMMPS[@thompson2022lammps], GOMC[@nejahi2021update], and HOOMD-blue[@anderson2020hoomd]).


# Statement of need

The General Molecular Simulation Object (GMSO) is a component of the Molecular Simulation Design Framework (MoSDeF), provides a framework and utilities for storing, manipulating, and outputting of molecular systems. MoSDeF is a suite of software tailored to facilitate the initialization of chemical and biomolecular systems for computational simulations[@cummings2021opena]. These tools were developed to specifically address a critical aspect of the (ir)reproducibility issue within the molecular simulation community — namely, the insufficient documentation of the structure preparation process and force field parameter implementation[@thompson2020towards]. The initialization step, often performed through Graphical User Interfaces (GUI) or via the use of ad-hoc, unpublished, and unreviewed code, poses the risk of introducing irreproducible and untraceable errors[@baker2016reproducibility]. By providing general-purposed and standardized tools that build and parameterize molecular systems for molecular simulations, directly support  various molecular dynamics and Monte Carlo engines, MoSDeF aims to trivialize the describing and disseminating such processes without creating extra burdens for computational simulation researchers[@cummings2021opena].


The initialization of chemical/biomolecular systems comprises of three key steps:
1.  Constructing structures: Encompassing loading and/or creating molecules/structures that mirrors the phenomena under investigation.
2. Parameterizing: Assigning interactional parameters to all particles and connections within the structures.
3. Storing Structures and Output Generation: Storing parameterized structures, and outputting to file formats compatible with various simulation software.

Each of these steps necessitates distinct routines, and as such, is addressed by a series of specialized libraries — specifically, mBuild[@klein2016hierarchical], Foyer[@klein2019formalizing], and GMSO, which is introduced in this work. mBuild functions as a molecular builder, equipped with extensive utilities for creating, loading, and manipulating positions of atoms and molecules, along with managing their connectivity through bonds[@klein2016hierarchical]. Foyer assumes the role of parameterizing for the created structures, involving the identification and assignment of interaction parameters to each atom or group of atoms and their associated connections (e.g., bonds, angles, and dihedrals)[@klein2019formalizing]. This process entails matching the connectivity (bond graph) of the provided structure with the SMARTS grammar of the corresponding atom type, defining the interactional parameters[@klein2019formalizing]. The use of a graph matching method, departing from the traditional approach of matching via atom indices, allows for a more flexible parameterization. This feature proves particularly advantageous in the study of functionalized polymers, whose structures consistently deviate slightly from the standard polymer[@summers2020mosdef, @quach2022high]. These utilities have been utilized in various projects to explore a wide range of structures and applications[@thompson2019scalable; @summers2020mosdef; @quach2022high; @ma2022dynamics], and integrated into other scientific libraries[@albooyeh2023flowermd; @defever2021mosdef; @crawford2023mosdefgomc].


The parameterization step introduces additional information, requiring a more sophisticated data structure for representation. Beyond the initial details concerning positions and connectivity established during system construction, the new structure incorporates supplementary metadata and interaction parameters. These encompass not only atoms and bonds but also extend to angles, dihedrals, and improper dihedrals, along with interaction parameters associated with each of these objects. This requires the development of data structures that are capable of:
- Supporting a variety of models
- Providing flexibility for exotic potentials and unit systems
- Being compatible with existing community tools
- Being extensible (to support new simulation models/engines/workflows)

Currently, existing data structures, such as ParmEd and OpenMM[@shirts2016lessons; @eastmann2017openmm], fulfill many functionalities and are widely adopted[@elenareal2023real; @kehrein2023unravel; @tesei2021accurate; @marrink2019computational]. However, their underlying structures are tailored to specific subsets of simulation workflows and ecosystems, as well as force field equation forms, sacrificing generality and broad applicability. This limitation includes hard-coding and assumptions about potential expressions and units. They lack the generality that MoSDeF and its users seek, such as the ability to define and store arbitrary potential expressions or unit systems. Integrating these new features into existing software, unfortunately, would require a major overhaul, potentially impacting existing simulation workflows and is not appealing to current project stakeholders.


Hence, we developed the General Molecular Simulation Object (GMSO) library, which is a lightweight and extensible data structure encapsulating chemical/biomolecular systems and their associated interaction parameters, i.e., force fields, to cater to the general force fields. The library is designed to accommodate a wide range of chemical/biomolecular models, offering the capability to support arbitrary potential expressions and unit systems. Generalizing these potential (force field) expressions allows users to enter the force field in its native form and units, minimizing user error when setting up the force field file while providing the ability to easily auto-convert the potential form and units to the molecular engine's required form. GMSO satisfies the broader community's need for a general, extensible, and reproducible method of setting up molecular simulations. In addition to core data classes, the library includes routines for interacting/converting to and from other ecosystems, including ParmEd and OpenMM, enhancing interoperability without reinventing functionalities. GMSO supports output to multiple engine-specific molecular file formats, currently including: GROMACS[@abraham2015gromacs], LAMMPS[@thompson2022lammps], HOOMD-Blue[@anderson2020hoomd], NAMD[@phillips2020scalable], CP2K[@kuhne2020cp2k], Cassandra[@shah2017cassandra], and GOMC[@nejahi2021update], with plans for future expansion. When integrated with other MoSDeF software and workflow manager like Signac[@adorf2018simple], GMSO facilitates large-scale automated molecular screening for diverse molecules/structures and state points, which is critical for developing new materials, chemicals and drugs[@craven2021examining, @quach2022high, @thompson2019scalable].


# Figures

<p align="center">
  <img src="summary_fig.png" alt="Summary of the GMSO Workflow" style="width:50%">
  <figcaption style="text-align: center";>
    <b>Figure 1. Summary of the GMSO workflow.</b>
  </figcaption>
</p>

# Acknowledgements
This research was partially supported by the National Science Foundation OAC-1835713 and OAC-1835874. Atomfold also donated research and development time and computational resources for this research and software.
