Core Classes
============

The ``gmso.core`` module contains the concrete data structures used to represent a full
chemical topology.  All classes are importable directly from the ``gmso`` namespace.

.. contents:: Contents
   :local:
   :depth: 2

Topology
--------

:class:`~gmso.Topology` is the central data structure in GMSO.  It acts as the container
for all sites, connections, and forcefield parameters that together describe a molecular
system ready for simulation.

Summary
~~~~~~~

.. autosummary::
   :nosignatures:

   gmso.Topology

Reference
~~~~~~~~~

.. autoclass:: gmso.Topology
   :members:
   :show-inheritance:

Sites
-----

Sites are the fundamental interaction points in a topology.  GMSO ships two concrete site
types: :class:`~gmso.Atom` for element-bearing particles, and
:class:`~gmso.VirtualSite` for massless interaction centres (e.g., TIP4P water models).

.. autosummary::
   :nosignatures:

   gmso.Atom
   gmso.VirtualSite

Atom
~~~~

.. autoclass:: gmso.Atom
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

VirtualSite
~~~~~~~~~~~

.. autoclass:: gmso.VirtualSite
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Connections
-----------

Connections encode the covalent topology of a system.  They are built from ordered sequences
of :class:`~gmso.Atom` members and an optional connection type that carries the forcefield
parameters.

Bonded connections must be defined, either in the file used to create the `gmso.Topology` or
in the conversion from another object into a `gmso.Topology`. If the bond connections are
correctly identified, you can use :func:`gmso.utils.connectivity.identify_connections()` to infer
and populate the rest of the connections (i.e., angles, dihedrals and improper dihedrals.)

.. autosummary::
   :nosignatures:

   gmso.Bond
   gmso.Angle
   gmso.Dihedral
   gmso.Improper

Bond
~~~~

.. autoclass:: gmso.Bond
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Angle
~~~~~

.. autoclass:: gmso.Angle
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Dihedral
~~~~~~~~

.. autoclass:: gmso.Dihedral
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Improper
~~~~~~~~

.. autoclass:: gmso.Improper
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Potential Types
---------------

Potential types carry the functional form and parameters for a given interaction class.
They attach to the corresponding connection (or site) via the ``*_type`` attributes.

.. autosummary::
   :nosignatures:

   gmso.AtomType
   gmso.BondType
   gmso.AngleType
   gmso.DihedralType
   gmso.ImproperType
   gmso.PairPotentialType
   gmso.VirtualType

AtomType
~~~~~~~~

.. autoclass:: gmso.AtomType
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

BondType
~~~~~~~~

.. autoclass:: gmso.BondType
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

AngleType
~~~~~~~~~

.. autoclass:: gmso.AngleType
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

DihedralType
~~~~~~~~~~~~

.. autoclass:: gmso.DihedralType
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

ImproperType
~~~~~~~~~~~~

.. autoclass:: gmso.ImproperType
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

PairPotentialType
~~~~~~~~~~~~~~~~~

.. autoclass:: gmso.PairPotentialType
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

VirtualType
~~~~~~~~~~~

.. autoclass:: gmso.VirtualType
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

ForceField
----------

:class:`~gmso.ForceField` parses and stores the potential types and parameters defined in a
forcefield XML file (Foyer/OpenMM format).  It is the primary input to
:func:`gmso.parameterization.apply`.

.. autosummary::
   :nosignatures:

   gmso.ForceField

.. autoclass:: gmso.ForceField
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Box
---

:class:`~gmso.Box` stores periodic boundary information for a simulation cell.

.. autoclass:: gmso.Box
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Element
-------

:class:`~gmso.Element` provides periodic-table data (mass, symbol, atomic number, etc.).

.. autoclass:: gmso.Element
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields
