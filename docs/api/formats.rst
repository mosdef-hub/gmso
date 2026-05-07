File Format Readers and Writers
================================

The ``gmso.formats`` module provides readers and writers for common molecular simulation
file formats.  All public functions are importable from ``gmso.formats``.

.. contents:: Contents
   :local:
   :depth: 2

GROMACS
-------

GROMACS coordinate files (`.gro`) and topology/parameter files (`.top`, `.itp`) are
supported for both reading and writing.

.. autosummary::
   :nosignatures:

   gmso.formats.read_gro
   gmso.formats.write_gro
   gmso.formats.write_top

read_gro
~~~~~~~~

.. autofunction:: gmso.formats.read_gro

write_gro
~~~~~~~~~

.. autofunction:: gmso.formats.write_gro

write_top
~~~~~~~~~

.. autofunction:: gmso.formats.write_top

GSD
---

GSD (General Simulation Data) is the native trajectory format for HOOMD-blue.

.. autosummary::
   :nosignatures:

   gmso.formats.write_gsd

write_gsd
~~~~~~~~~

.. autofunction:: gmso.formats.write_gsd

XYZ
---

The plain-text XYZ format stores atomic symbols and Cartesian coordinates.

.. autosummary::
   :nosignatures:

   gmso.formats.read_xyz
   gmso.formats.write_xyz

read_xyz
~~~~~~~~

.. autofunction:: gmso.formats.read_xyz

write_xyz
~~~~~~~~~

.. autofunction:: gmso.formats.write_xyz

LAMMPS DATA
-----------

LAMMPS data files encode particle positions, topology, and forcefield coefficients.

.. autosummary::
   :nosignatures:

   gmso.formats.write_lammpsdata
   gmso.formats.read_lammpsdata

write_lammpsdata
~~~~~~~~~~~~~~~~

.. autofunction:: gmso.formats.write_lammpsdata

read_lammpsdata
~~~~~~~~~~~~~~~

.. autofunction:: gmso.formats.read_lammpsdata

MOL2
----

Tripos MOL2 is a common format for small-molecule structures.

.. autosummary::
   :nosignatures:

   gmso.formats.read_mol2
   gmso.formats.write_mol2

read_mol2
~~~~~~~~~

.. autofunction:: gmso.formats.read_mol2

write_mol2
~~~~~~~~~~

.. autofunction:: gmso.formats.write_mol2

MCF
---

Monte Carlo forcefield (MCF) files are used by Cassandra Monte Carlo.

.. autosummary::
   :nosignatures:

   gmso.formats.write_mcf

write_mcf
~~~~~~~~~

.. autofunction:: gmso.formats.write_mcf

JSON
----

GMSO topologies can be serialized to and from JSON for storage and interoperability.

.. autosummary::
   :nosignatures:

   gmso.formats.write_json
   gmso.formats.load_json

write_json
~~~~~~~~~~

.. autofunction:: gmso.formats.write_json

load_json
~~~~~~~~~

.. autofunction:: gmso.formats.load_json

Format Registry
---------------

The format registry provides the decorator-based mechanism for registering new file
format readers and writers.

.. autoclass:: gmso.formats.formats_registry.LoadersRegistry
   :members:
   :show-inheritance:

.. autoclass:: gmso.formats.formats_registry.SaversRegistry
   :members:
   :show-inheritance:
