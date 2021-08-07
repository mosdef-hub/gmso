=================
Design Principles
=================

Scope and Features of ``GMSO``
----------------------------------

``GMSO`` is designed to enable the flexible, general representation of
chemical topologies for molecular simulation. Efforts are made to enable
lossless, bias-free storage of data, without assuming particular chemistries,
models, or using any particular engine's ecosystem as a starting point. The
scope is generally restrained to the preparation, manipulation, and conversion
of and of input files for molecular simulation, i.e. before engines are called
to execute the simulations themselves. ``GMSO`` currently does not support
conversions between trajectory file formats for analysis codes. In the scope of
molecular simulation, we loosely define a chemical topology as everything
needed to reproducibly prepare a chemical system for simulation. This includes
particle coordinates and connectivity, box information, force field data
(functional forms, parameters tagged with units, partial charges, etc.) and
some optional information that may not apply to all systems (i.e. specification
of elements with each particle).

``GMSO`` enables the following features:

* Supporting a variety of models in the molecular simulation/computational
  chemistry community: No assumptions are made about an interaction site
  representing an atom or bead, instead supported atomistic,
  united-atom/coarse-grained, polarizable, and other models!

* Greater flexibility for exotic potentials: The ``AtomType`` (and analogue
  classes for intramolecular interactions) uses ``sympy`` to store any
  potential that can be represented by a mathematical expression. If you can
  write it down, it can be stored!

* Easier development for glue to new engines: by not being designed for
  compatibility with any particular molecular simulation engine or ecosystem,
  it becomes more tractable for developers in the community to add glue for
  engines that are not currently supported (and even ones that do not exist at
  present)!


* Compatibility with existing community tools: No single molecular simulation
  tool will be a silver bullet, so ``GMSO`` includes functions to convert
  objects. These can be used in their own right to convert between objects in-memory
  and also to support conversion to file formats not natively supported at
  any given time. Currently supported conversions include ``ParmEd``,
  ``OpenMM``, ``mBuild``, ``MDTraj``, with others coming in the future!


* Native support for reading and writing many common file formats (``XYZ``,
  ``GRO``, ``TOP``, ``LAMMPSDATA``) and indirect support, through other
  libraries, for many more!


Structure of ``GMSO``
-------------------------
There are three main modules within the Python package:

* ``gmso.core`` stores the classes that constitute the core data structures.
* ``gmso.formats`` stores readers and writers for (on-disk) file formats.
* ``gmso.external`` includes functions that convert core data structures between external libraries and their internal representation.
