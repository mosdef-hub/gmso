External
=========

This submodule includes functions that convert core data structures between external libraries and their internal representation.

mBuild
-------
The following methods are available for converting `mBuild <https://mbuild.mosdef.org>`_ objects to and from ``GMSO``.

    .. autofunction:: gmso.external.from_mbuild
    .. autofunction:: gmso.external.to_mbuild


Parmed
------
Conversion methods for `Parmed <https://parmed.github.io/ParmEd/html/index.html>`_ objects to and from ``GMSO``.

    .. autofunction:: gmso.external.from_parmed

OpenMM
-------
Conversion methods for `OpenMM <http://openmm.org/>`_ objects to and from ``GMSO``.

    .. autofunction:: gmso.external.to_openmm
