Parameterization
================

The ``gmso.parameterization`` module provides functions and classes for applying a
:class:`~gmso.ForceField` to an un-typed :class:`~gmso.Topology`.  The primary entry point
is :func:`~gmso.parameterization.apply`.

.. contents:: Contents
   :local:
   :depth: 2

apply
-----

:func:`~gmso.parameterization.apply` is the main public API for parameterizing a topology.
It accepts a :class:`~gmso.ForceField` (or a dictionary of forcefields keyed by molecule
name) and writes atom types, bond types, angle types, dihedral types, and improper types
onto the topology in-place.

.. autofunction:: gmso.parameterization.parameterize.apply
