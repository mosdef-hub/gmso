External Library Conversions
============================

The ``gmso.external`` module provides functions to convert :class:`~gmso.Topology` objects
to and from in-memory representations used by other molecular simulation libraries.
All public functions are importable from ``gmso.external``.

.. contents:: Contents
   :local:
   :depth: 2

mBuild
------

`mBuild <https://mbuild.mosdef.org>`_ is a hierarchical molecule builder.  The conversions
below allow a built ``mBuild.Compound`` to be typed and written to simulation files via
GMSO, and vice-versa.

.. autosummary::
   :nosignatures:

   gmso.external.from_mbuild
   gmso.external.from_mbuild_box
   gmso.external.to_mbuild

from_mbuild
~~~~~~~~~~~

.. autofunction:: gmso.external.from_mbuild

from_mbuild_box
~~~~~~~~~~~~~~~

.. autofunction:: gmso.external.from_mbuild_box

to_mbuild
~~~~~~~~~

.. autofunction:: gmso.external.to_mbuild

ParmEd
------

`ParmEd <https://parmed.github.io/ParmEd>`_ is a library for editing and converting
molecular mechanics parameter files.

.. autosummary::
   :nosignatures:

   gmso.external.from_parmed
   gmso.external.to_parmed

from_parmed
~~~~~~~~~~~

.. autofunction:: gmso.external.from_parmed

to_parmed
~~~~~~~~~

.. autofunction:: gmso.external.to_parmed

OpenMM
------

`OpenMM <https://openmm.org>`_ is a GPU-accelerated molecular dynamics library.  GMSO can
export a typed topology as an OpenMM ``Topology`` object.

.. autosummary::
   :nosignatures:

   gmso.external.to_openmm

to_openmm
~~~~~~~~~

.. autofunction:: gmso.external.to_openmm

HOOMD-blue
----------

`HOOMD-blue <https://hoomd-blue.readthedocs.io>`_ is a particle simulation toolkit.
GMSO can export directly to HOOMD-blue's GSD snapshot and forcefield objects.

.. autosummary::
   :nosignatures:

   gmso.external.to_hoomd_snapshot
   gmso.external.to_hoomd_forcefield
   gmso.external.to_gsd_snapshot

to_hoomd_snapshot
~~~~~~~~~~~~~~~~~

.. autofunction:: gmso.external.to_hoomd_snapshot

to_hoomd_forcefield
~~~~~~~~~~~~~~~~~~~

.. autofunction:: gmso.external.to_hoomd_forcefield

to_gsd_snapshot
~~~~~~~~~~~~~~~

.. autofunction:: gmso.external.to_gsd_snapshot

NetworkX
--------

GMSO topologies can be converted to and from `NetworkX <https://networkx.org>`_ graphs for
graph-based analysis.

.. autosummary::
   :nosignatures:

   gmso.external.from_networkx
   gmso.external.to_networkx

from_networkx
~~~~~~~~~~~~~

.. autofunction:: gmso.external.from_networkx

to_networkx
~~~~~~~~~~~

.. autofunction:: gmso.external.to_networkx
