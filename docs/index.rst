.. GMSO documentation master file

GMSO: Flexible storage of chemical topology for molecular simulation
====================================================================

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

`GMSO` is a flexible storage of chemical topology for molecular simulation.
With a few lines of `GMSO` code, together with `mBuild <https://mbuild.mosdef.org>`_ and
`foyer <https://foyer.mosdef.org>`_, users can rapidly prototype arbitrary parameterized
chemical systems and generate data files for a wide variety of simulation engines.


GMSO is a part of the MoSDeF ecosystem
---------------------------------------

`GMSO` is designed to be a general and flexible representation of chemical topologies for
molecular simulation.  With an emphasis on assuming as little as possible about the chemical
system, model, or engine, `GMSO` can enable support for a variety of systems.

`GMSO` is a part of the `MoSDeF (Molecular Simulation and Design Framework)
<https://mosdef.org>`_ ecosystem, and is intended to be a generalized alternative for the
`foyer package <https://foyer.mosdef.org>`_.  Libraries in the MoSDeF ecosystem are
designed to provide utilities necessary to streamline a researcher's simulation workflow.
When setting up simulation studies we also recommend users follow the `TRUE
<https://www.tandfonline.com/doi/full/10.1080/00268976.2020.1742938>`_
(Transparent, Reproducible, Usable-by-others, and Extensible) standard.


Goals and Features
------------------

`GMSO`'s goal is to provide a flexible backend framework to store topological information
of a chemical system in a reproducible fashion.  **Topology** in this case is defined as
the information needed to initialise a molecular simulation. Depending on the type of
simulation performed, this ranges from:

* particle positions
* particle connectivity
* box information
* forcefield data

    * functional forms defined as `sympy <https://www.sympy.org>`_ expressions
    * parameters with defined units
    * partial charges
    * tabulated data
    * etc.

* Other optional data

    * particle mass
    * elemental data

With these driving goals for `GMSO`, the following features are enabled:

#.  **Supporting a variety of models**: no assumptions are made about an interaction site
    representing an atom or bead — these can be atomistic, united-atom/coarse-grained,
    polarisable, and other models.

#.  **Greater flexibility for exotic potentials**: :class:`~gmso.AtomType` (and analogous
    classes for intramolecular interactions) uses `sympy <https://www.sympy.org>`_ to store
    any potential that can be represented by a mathematical expression.

#.  **Adaptable for new engines**: by not being designed for compatibility with any
    particular simulation engine, it is more tractable for community developers to add
    support for engines not currently supported.

#.  **Compatibility with existing community tools**: GMSO includes functions to convert
    between various file formats and libraries.  Currently supported conversions include:

    * `ParmEd`
    * `OpenMM`
    * `mBuild`
    * `HOOMD-blue`
    * `NetworkX`

#.  **Native support for reading and writing many common file formats**:

    * XYZ
    * GRO / TOP / ITP
    * LAMMPS DATA
    * GSD
    * MOL2
    * MCF
    * JSON


.. toctree::
   :hidden:
   :caption: User Guide

   installation
   design_principles
   contributing
   docker

.. toctree::
   :hidden:
   :caption: API Reference

   api/index
