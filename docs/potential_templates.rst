.. _potential_templates:

Potential Templates
===================

GMSO ships a library of built-in potential templates stored as JSON files in
``gmso/lib/jsons/``.  Each template defines a named functional form, its
mathematical expression (as a SymPy-compatible string), and the physical
dimensions expected for each parameter.

Templates are loaded via :class:`~gmso.utils.io.PotentialTemplateLibrary` and
are used internally by forcefield parsers and writers to validate that
parameters carry the correct units.

The more common use case will be to use the expressions of these templates in
a GMSO `.xml` file, where the independent parameters for the expression are defined
for each unique interaction type.

.. contents:: Contents
   :local:
   :depth: 2

----

Non-bonded Potentials
---------------------

LennardJonesPotential
~~~~~~~~~~~~~~~~~~~~~

The classic 12-6 Lennard-Jones pair potential.  Widely used for van der Waals
interactions in atomistic and united-atom forcefields.

.. math::

   U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}
          - \left(\frac{\sigma}{r}\right)^{6}\right]

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``epsilon``
     - energy
     - Well depth; controls the strength of the interaction.
   * - ``sigma``
     - length
     - Finite distance at which the potential is zero; related to particle size.

----

BuckinghamPotential
~~~~~~~~~~~~~~~~~~~

Exponential-6 (Buckingham) pair potential.  Replaces the steep repulsive
:math:`r^{-12}` wall of Lennard-Jones with a physically motivated exponential
repulsion term.

.. math::

   U(r) = a\exp(-b\,r) - c\,r^{-6}

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``a``
     - energy
     - Prefactor for the repulsive exponential term.
   * - ``b``
     - 1/length
     - Exponent controlling the range of repulsion.
   * - ``c``
     - energy·length\ :sup:`6`
     - Prefactor for the attractive dispersion term.

----

MiePotential
~~~~~~~~~~~~

Generalised Mie pair potential.  Lennard-Jones is the special case
:math:`n = 12`, :math:`m = 6`.

.. math::

   U(r) = \frac{n}{n-m}\left(\frac{n}{m}\right)^{m/(n-m)}
          \epsilon\left[\left(\frac{\sigma}{r}\right)^n
          - \left(\frac{\sigma}{r}\right)^m\right]

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``epsilon``
     - energy
     - Well depth.
   * - ``sigma``
     - length
     - Length-scale parameter.
   * - ``n``
     - dimensionless
     - Repulsive exponent.
   * - ``m``
     - dimensionless
     - Attractive exponent.

----

Bond Potentials
---------------

HarmonicBondPotential
~~~~~~~~~~~~~~~~~~~~~

Standard harmonic (quadratic) bond-stretching potential.

.. math::

   U(r) = \frac{1}{2}\,k\,(r - r_{eq})^2

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy/length\ :sup:`2`
     - Force constant.
   * - ``r_eq``
     - length
     - Equilibrium bond length.

----

LAMMPSHarmonicBondPotential
~~~~~~~~~~~~~~~~~~~~~~~~~~~

LAMMPS convention for the harmonic bond.  The factor of ½ is absorbed into
``k``, so the force constant value is twice the physical spring constant
compared to :ref:`HarmonicBondPotential`.

.. math::

   U(r) = k\,(r - r_{eq})^2

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy/length\ :sup:`2`
     - Force constant (note: equals 2× the physical spring constant).
   * - ``r_eq``
     - length
     - Equilibrium bond length.

----

FixedBondPotential
~~~~~~~~~~~~~~~~~~

Rigid constraint that fixes the bond length to exactly ``r_eq``.  Expressed
as a Dirac delta to be consistent with the potential template framework.

.. math::

   U(r) = \delta(r - r_{eq})

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``r_eq``
     - length
     - Fixed (constrained) bond length.

----

FENEBondPotential
~~~~~~~~~~~~~~~~~

Finitely Extensible Nonlinear Elastic (FENE) bond.  Used in coarse-grained
polymer models to prevent bond extension beyond a maximum length.

.. math::

   U(r) = -\frac{1}{2}\,k\,r_{eq}^2\,\ln\!\left[1
          - \left(\frac{r}{r_{eq}}\right)^2\right]

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy/length\ :sup:`2`
     - Spring constant.
   * - ``r_eq``
     - length
     - Maximum extensible bond length.

----

LAMMPSFENEBondPotential
~~~~~~~~~~~~~~~~~~~~~~~

FENE bond with an embedded Weeks-Chandler-Andersen (WCA) repulsion, as
implemented in LAMMPS.

.. math::

   U(r) = -\frac{1}{2}\,K\,R_0^2\,\ln\!\left[1 - \left(\frac{r}{R_0}\right)^2\right]
          + 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}
          - \left(\frac{\sigma}{r}\right)^{6}\right] + \epsilon

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``K``
     - energy/length\ :sup:`2`
     - FENE spring constant.
   * - ``R0``
     - length
     - Maximum extensible bond length.
   * - ``epsilon``
     - energy
     - WCA well depth.
   * - ``sigma``
     - length
     - WCA length-scale parameter.

----

HOOMDFENEWCABondPotential
~~~~~~~~~~~~~~~~~~~~~~~~~

FENE + WCA bond as implemented in HOOMD-blue.  Adds a displacement offset
``delta`` relative to the LAMMPS form.

.. math::

   U(r) = -\frac{1}{2}\,k\,r_0^2\,\ln\!\left[1
          - \left(\frac{r-\delta}{r_0}\right)^2\right]
          + 4\epsilon\left[\left(\frac{\sigma}{r-\delta}\right)^{12}
          - \left(\frac{\sigma}{r-\delta}\right)^{6}\right] + \epsilon

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy/length\ :sup:`2`
     - FENE spring constant.
   * - ``r0``
     - length
     - Maximum extensible bond length.
   * - ``epsilon``
     - energy
     - WCA well depth.
   * - ``sigma``
     - length
     - WCA length-scale parameter.
   * - ``delta``
     - length
     - Position offset applied before evaluating the potential.

----

Angle Potentials
----------------

HarmonicAnglePotential
~~~~~~~~~~~~~~~~~~~~~~

Standard harmonic angle-bending potential.

.. math::

   U(\theta) = \frac{1}{2}\,k\,(\theta - \theta_{eq})^2

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy/angle\ :sup:`2`
     - Bending force constant.
   * - ``theta_eq``
     - angle
     - Equilibrium bond angle.

----

LAMMPSHarmonicAnglePotential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

LAMMPS convention for the harmonic angle.  The factor of ½ is absorbed into
``k``.

.. math::

   U(\theta) = k\,(\theta - \theta_{eq})^2

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy/angle\ :sup:`2`
     - Force constant (equals 2× the physical bending constant).
   * - ``theta_eq``
     - angle
     - Equilibrium bond angle.

----

FixedAnglePotential
~~~~~~~~~~~~~~~~~~~

Rigid angle constraint.

.. math::

   U(\theta) = \delta(\theta - \theta_{eq})

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``theta_eq``
     - angle
     - Fixed (constrained) bond angle.

----

Torsion / Dihedral Potentials
------------------------------

PeriodicTorsionPotential
~~~~~~~~~~~~~~~~~~~~~~~~~

Single-term periodic (cosine) dihedral.  The default torsion form in many
AMBER- and CHARMM-derived forcefields.

.. math::

   U(\phi) = k\left[1 + \cos(n\phi - \phi_{eq})\right]

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy
     - Torsion barrier height.
   * - ``n``
     - dimensionless
     - Periodicity (multiplicity) of the dihedral.
   * - ``phi_eq``
     - angle
     - Phase offset.

----

HarmonicTorsionPotential
~~~~~~~~~~~~~~~~~~~~~~~~~

Harmonic dihedral restraint.  Penalises deviations from a reference torsion
angle.

.. math::

   U(\phi) = \frac{1}{2}\,k\,(\phi - \phi_{eq})^2

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy/angle\ :sup:`2`
     - Torsion force constant.
   * - ``phi_eq``
     - angle
     - Reference dihedral angle.

----

OPLSTorsionPotential
~~~~~~~~~~~~~~~~~~~~~

Four-term cosine expansion used in the OPLS forcefield family.

.. math::

   U(\phi) = \frac{1}{2}k_1(1+\cos\phi)
           + \frac{1}{2}k_2(1-\cos 2\phi)
           + \frac{1}{2}k_3(1+\cos 3\phi)
           + \frac{1}{2}k_4(1-\cos 4\phi)

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k1``
     - energy
     - Coefficient for the first cosine term.
   * - ``k2``
     - energy
     - Coefficient for the second cosine term.
   * - ``k3``
     - energy
     - Coefficient for the third cosine term.
   * - ``k4``
     - energy
     - Coefficient for the fourth cosine term.

----

FourierTorsionPotential
~~~~~~~~~~~~~~~~~~~~~~~~

Five-term Fourier dihedral expansion. Extends the OPLS form with a constant
offset term ``k0``.

.. math::

   U(\phi) = \frac{1}{2}k_0
           + \frac{1}{2}k_1(1+\cos\phi)
           + \frac{1}{2}k_2(1-\cos 2\phi)
           + \frac{1}{2}k_3(1+\cos 3\phi)
           + \frac{1}{2}k_4(1-\cos 4\phi)

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k0``
     - energy
     - Constant offset.
   * - ``k1``
     - energy
     - Coefficient for the first cosine term.
   * - ``k2``
     - energy
     - Coefficient for the second cosine term.
   * - ``k3``
     - energy
     - Coefficient for the third cosine term.
   * - ``k4``
     - energy
     - Coefficient for the fourth cosine term.

----

RyckaertBellemansTorsionPotential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Polynomial expansion in :math:`\cos\phi`, used in the GROMOS and early AMBER
forcefields.  Related to the Fourier form by a trigonometric identity.

.. math::

   U(\phi) = \sum_{n=0}^{5} c_n \cos^n\!\phi

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``c0``
     - energy
     - 0th-order coefficient (constant).
   * - ``c1``
     - energy
     - 1st-order coefficient.
   * - ``c2``
     - energy
     - 2nd-order coefficient.
   * - ``c3``
     - energy
     - 3rd-order coefficient.
   * - ``c4``
     - energy
     - 4th-order coefficient.
   * - ``c5``
     - energy
     - 5th-order coefficient.

----

LAMMPSHarmonicDihedralPotential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Single-term periodic dihedral in LAMMPS convention.

.. math::

   U(\phi) = K\left[1 + d\cos(n\phi)\right]

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``K``
     - energy
     - Barrier height.
   * - ``d``
     - dimensionless
     - Phase sign (typically +1 or −1).
   * - ``n``
     - dimensionless
     - Periodicity (multiplicity).

----

HOOMDPeriodicDihedralPotential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Periodic dihedral as implemented in HOOMD-blue.  Equivalent to the LAMMPS
form but adds a phase offset ``phi0``.

.. math::

   U(\phi) = \frac{1}{2}\,k\left[1 + d\cos(n\phi - \phi_0)\right]

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy
     - Barrier height.
   * - ``d``
     - dimensionless
     - Phase sign (typically +1 or −1).
   * - ``n``
     - dimensionless
     - Periodicity (multiplicity).
   * - ``phi0``
     - angle
     - Phase offset.

----

Improper Potentials
-------------------

HarmonicImproperPotential
~~~~~~~~~~~~~~~~~~~~~~~~~~

Harmonic improper dihedral restraint.  Keeps a central atom and its three
bonded neighbours coplanar (or at a fixed out-of-plane angle).

.. math::

   U(\phi) = \frac{1}{2}\,k\,(\phi - \phi_{eq})^2

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy/angle\ :sup:`2`
     - Improper force constant.
   * - ``phi_eq``
     - angle
     - Reference out-of-plane angle (0 for planar groups).

----

PeriodicImproperPotential
~~~~~~~~~~~~~~~~~~~~~~~~~~

Periodic (cosine) improper dihedral.

.. math::

   U(\phi) = k\left[1 + \cos(n\phi - \phi_{eq})\right]

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``k``
     - energy
     - Barrier height.
   * - ``n``
     - dimensionless
     - Periodicity.
   * - ``phi_eq``
     - angle
     - Phase offset.

----

Virtual Sites
-------------

TIP4PPotential (MSite)
~~~~~~~~~~~~~~~~~~~~~~

Lennard-Jones interaction centred on the virtual M-site used in 4-point water
models such as TIP4P.  The M-site itself carries no charge in this template;
electrostatics are handled separately.

.. math::

   U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}
          - \left(\frac{\sigma}{r}\right)^{6}\right]

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``epsilon``
     - energy
     - Well depth.
   * - ``sigma``
     - length
     - Length-scale parameter.

----

Type3fdVirtualPosition
~~~~~~~~~~~~~~~~~~~~~~~

Position rule for a type 3fd virtual site: placed along the bisector of two
bond vectors, at a fixed distance from the central atom.  Used to construct
out-of-plane virtual sites in 4-point and 5-point water models.

.. math::

   \mathbf{r}_v = \mathbf{r}_i
     + b\,\frac{\mathbf{r}_j - \mathbf{r}_i
       + a\,(\mathbf{r}_k - \mathbf{r}_j)}
      {\|\mathbf{r}_j - \mathbf{r}_i + a\,(\mathbf{r}_k - \mathbf{r}_j)\|}

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``a``
     - dimensionless
     - Fractional weight along the :math:`\mathbf{r}_j - \mathbf{r}_i` vector.
   * - ``b``
     - dimensionless
     - Distance from :math:`\mathbf{r}_i` to the virtual site (in units of the normalised vector length).

----

Engine-Specific Potentials
---------------------------

HOOMDDPDForce
~~~~~~~~~~~~~

Dissipative Particle Dynamics (DPD) conservative force as implemented in
HOOMD-blue.  The full DPD force also includes dissipative and random terms,
which are set at the integrator level and are not encoded here.

.. math::

   U(r) = A\left(1 - \frac{r}{r_{cut}}\right) - \gamma

.. list-table::
   :header-rows: 1
   :widths: 15 25 60

   * - Parameter
     - Dimensions
     - Description
   * - ``A``
     - force
     - Conservative force amplitude.
   * - ``r_cut``
     - length
     - Cutoff radius; the force is zero for :math:`r \geq r_{cut}`.
   * - ``γ``
     - mass/(length·time)
     - Friction coefficient for the dissipative term.
