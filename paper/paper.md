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
- name: School of Engineering and Physical 551 Sciences, Heriot-Watt University, Edinburgh, Scotland, U.K
  index: 9


dates: 2 January, 2024

bibliography: paper.bib


# Summary
The General Molecular Simulation Object, or GMSO, stands as an open-source Python data structure, offering a versatile and expandable framework for handling chemical and biomolecular topologies. This library is an integral component of the Molecular Simulation Design Framework (MoSDeF), dedicated to streamlining the creation, parameterization, and representation of chemical systems for molecular simulations. The GMSO library serves as a dynamic repository for storing chemical/biomolecular structures, encompassing metadata, coordinates, and interaction potentials. Moreover, the library includes routines for exporting stored structures into various file formats, facilitating compatibility with other software for visualization (e.g., VMD and OVITO) or conducting molecular simulations (e.g., GROMACS, LAMMPS, GOMC).


# Statement of need

The Molecular Simulation Design Framework (MoSDeF) is a suite of software tailored to facilitate the initialization of chemical and biomolecular systems for computational simulations [@cummings2021opena]. These tools were developed to specifically address a critical aspect of the (ir)reproducibility issue within the molecular simulation community — namely, the insufficient documentation of the structure preparation process[@thompson2020towards]. The initialization step, often performed through Graphical User Interfaces (GUI) or via the use of ad-hoc, unpublished, and unreviewed code, poses the risk of introducing irreproducible and untraceable errors[@baker2016reproducibility]. By providing general-purposed and standardized tools, MoSDeF aims to trivialize the process of describing and dissiminating such process, without creating extra burden for computional simulation researcher.[@cummings2021opena]

The process of initializing chemical/biomolecular systems involves constructing structures, assigning interaction parameters, and generating output structures in file formats compatible with a multitude of simulation software, such as GROMACS, LAMMPS, or GOMC[@abraham2015gromacs; @thompson2022lammps; @nejahi2021update]. Each of these steps necessitates distinct routines, and as such, is addressed by a series of specialized libraries—specifically, mBuild [@klein2016hierarchical], Foyer [@klein2019formalizing], and GMSO, which will be elaborated upon in this work.

mBuild functions as a molecular builder, meticulously crafted with extensive utilities for creating, loading, and manipulating the positions of atoms and molecules, along with managing their connectivity through bonds[@klein2016hierarchical]. These utilities have been applied in various projects to explore a wide range of structures and applications, as well as integrated into other scientific libraries[@quach2022high; @albooyeh2023flowermd; @ma2022dynamics].

Foyer assumes the role of atom typing for the created structures, involving the identification and assignment of interaction parameters to each atom or group of atoms[@klein2019formalizing]. This process entails matching the connectivity (bond graph) of the provided structure with the SMARTS grammar of the corresponding atom type, defining the interactional parameters[@klein2019formalizing]. The use of a graph matching method, departing from the traditional approach of matching via atom indices, allows for a more flexible parameterization. This feature proves particularly advantageous in the study of functionalized polymers, whose structures consistently deviate slightly from the standard polymer[@quach2022high].

The parameterization step introduces additional information, requiring a more intricate data structure for representation. Beyond the initial details concerning positions and connectivity established during system construction, the new structure incorporates supplementary metadata and interaction parameters. These encompass not only atoms and bonds but also extend to angles, dihedrals, and improper dihedrals, along with specific interaction parameters associated with each of these objects. This necessitates the development of data structures that are capable of:
- Supporting a variety of models
- Providing flexibility for exotic potentials
- Being compatible with existing community tools
- Being extensible (to support new simulation models/engines/workflows)

Currently, there are data structures designed to represent these systems, such as ParmEd and OpenMM[@shirts2016lessons; @eastmann2017openmm]. However, these data structures are tailored to a specific subset of simulation workflows/ecosystems, sacrificing some generality. This limitation includes, but is not limited to, the hard-coding and assumption of potential (interactional) expressions and units. They lack the generality that MoSDeF seeks, such as the ability to define and store arbitrary potential expressions or unit systems. Integrating these new features would require a major overhaul of these data structures, potentially impacting existing simulation workflows and may not be appealing to current stakeholders in those projects. Hence, we have developed a new data structure called the General Molecular Simulation Object, or GMSO, specifically catering to the MoSDeF ecosystem. GMSO satisfies our needs for generality as well as extensibility.

The General Molecular Simulation Object (GMSO) library is a lightweight, extensible data structure encapsulating chemical/biomolecular systems and their associated interaction parameters, i.e., force fields. The library is designed to accommodate a wide range of chemical/biomolecular models, offering the capability to support arbitrary potential expressions and unit systems. In addition to the core data structure, the library includes routines essential for converting to and from other data structures, allowing the utilization of written parsers (enabling the creation of objects from disk). This improves interoperability with other ecosystems, while avoiding the reinvention of the wheel for well-established code. Furthermore, GMSO allows outputting to several molecular simulation engine-specific file formats, currently supporting GROMACS, LAMMPS, HOOMD-Blue, GOMC, Cassandra, with plans for future expansion.


# Acknowledgements
This research was partially supported by the National Science Foundation OAC-1835713 and OAC-1835874. Atomfold also donated research and development time and computational resources for this research and software.
