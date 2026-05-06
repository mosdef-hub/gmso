Abstract Base Classes
=====================

The ``gmso.abc`` module defines the abstract base classes that underpin all data structures in
GMSO. These classes are not meant to be instantiated directly; instead, they provide a shared
interface and common behaviour that concrete classes (e.g., :class:`gmso.Atom`,
:class:`gmso.Bond`) build upon.

.. contents:: Contents
   :local:
   :depth: 2

GMSOBase
--------

:class:`~gmso.abc.gmso_base.GMSOBase` is the root Pydantic ``BaseModel`` that every GMSO
object derives from. It enforces strict validation, forbids extra fields, and provides
JSON-serialisation helpers.

.. autoclass:: gmso.abc.gmso_base.GMSOBase
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Site
----

:class:`~gmso.abc.abstract_site.Site` is the abstract base for all interaction sites
(atoms, virtual sites, coarse-grained beads, etc.).  It stores a 3-D Cartesian position,
an optional name, and optional molecule/residue labels.

.. autoclass:: gmso.abc.abstract_site.Site
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Molecule
~~~~~~~~

A lightweight label that groups sites into named, numbered molecules.

.. autoclass:: gmso.abc.abstract_site.Molecule
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Residue
~~~~~~~

A lightweight label that groups sites into named, numbered residues.

.. autoclass:: gmso.abc.abstract_site.Residue
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

AbstractPotential
-----------------

:class:`~gmso.abc.abstract_potential.AbstractPotential` is the abstract base for all
potential-energy functions.  It stores the functional form as a :mod:`sympy` expression and
the parameters—with physical units—as a dictionary.

.. autoclass:: gmso.abc.abstract_potential.AbstractPotential
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields

Connection
----------

:class:`~gmso.abc.abstract_connection.Connection` is the abstract base for all topological
connections between sites (bonds, angles, dihedrals, impropers).

.. autoclass:: gmso.abc.abstract_connection.Connection
   :members:
   :show-inheritance:
   :exclude-members: model_config, model_fields
