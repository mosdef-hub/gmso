.. GMSO documentation master file, created by
   sphinx-quickstart on Mon Feb 24 11:46:58 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GMSO: Flexible storage of chemical topology for molecular simulation
====================================================================

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

`GMSO` is a flexible storage of chemical topology for molecular simulation.
With a few lines of `GMSO` code, together with `mBuild <https://mbuild.mosdef.org>`_ and `foyer <https://foyer.mosdef.org>`_, users can rapidly prototype arbitrary parameterized chemical systems and generate data files for a wide variety of simulation engines.


GMSO is a part of the MoSDeF ecosystem
-----------------------------------------
`GMSO` is designed to be a general and flexible representation of chemical topolgies for molecular simulation.
With an emphasis on assuming as little as possible about the chemical system, model, or engine, `GMSO` can enable support for a variety of systems.
`GMSO` is a part of the `MoSDeF (Molecular Simulation and Design Framework) <https://mosdef.org>`_ ecosystem, and is intended to be a generalized alternative for the `foyer package <https://foyer.mosdef.org>`_.
Libraries in the MoSDeF ecosystem are designed to provide utilities neccessary to streamline
a researcher's simulation workflow. When setting up simulation studies,
we also recommend users to follow the `TRUE <https://www.tandfonline.com/doi/full/10.1080/00268976.2020.1742938>`_
(Transparent, Reproducible, Usable-by-others, and Extensible) standard, which is a set of common
practices meant to improve the reproducibility of computational simulation research.


Goals and Features
------------------

`GMSO`'s goal is to provide a flexible backend framework to store topological information of a chemical system in a reproducible fashion.
**Topology** in this case is defined as the information needed to initialize a molecular simulation.
Depending on the type of simulation performed, this ranges from:

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
    * etc.

With these driving goals for `GMSO`, the following features are enabled:

    #.  **Supporting a variety of models**: in the molecular simulation/computational chemistry community: No assumptions are made about an interaction site representing an atom or bead, instead these can be atomistic, united-atom/coarse-grained, polarizable, and other models!

    #.  **Greater flexibility for exotic potentials**: The `AtomType` (and analogue classes for intramolecular interactions) uses `sympy <https://www.sympy.org>`_ to store any potential that can be represented by a mathematical expression.

    #.  **Adaptable for new engines**: by not being designed for compatibility with any particular molecular simulation engine or ecosystem, it becomes more tractable for developers in the community to add glue for engines that are not currently supported.

    #.  **Compatibility with existing community tools**: No single molecular simulation tool will ever be a silver bullet, so ``GMSO`` includes functions to convert between various file formats and libraries. These can be used in their own right to convert between objects in-memory and also to support conversion to file formats not natively supported at any given time. Currently supported conversions include:

        * `ParmEd`
        * `OpenMM`
        * `mBuild`
        * more in the future!

    #.  **Native support for reading and writing many common file formats**:

        * `XYZ`
        * `GRO`
        * `TOP`
        * `LAMMPSDATA`
        * `GSD`
        * indirect support, through other libraries, for many more!


.. toctree::
   :hidden:

   design_principles
   data_structures
   formats
   external
   installation
   docker
   contributing
