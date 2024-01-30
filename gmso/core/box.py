"""Represent box regions in gmso."""

import warnings

import numpy as np
import unyt as u
from unyt.array import allclose_units


def _validate_lengths(lengths):
    """Ensure the lengths of the box are positive and check dimension."""
    if not isinstance(lengths, u.unyt_array):
        if all(isinstance(length, u.unyt_quantity) for length in lengths):
            print("Converting list of unyt quantities to a unyt array")
            lengths = u.unyt_array([l for l in lengths], str(lengths[0].units))
        else:
            warnings.warn("Lengths are assumed to be in nm")
            lengths *= u.nm
    input_unit = lengths.units

    lengths = np.asarray(lengths, dtype=float, order="C")
    np.reshape(lengths, newshape=(3,), order="C")

    lengths *= input_unit
    if input_unit != u.Unit("dimensionless"):
        lengths.convert_to_units(u.nm)

    if np.any(
        np.less(
            lengths,
            [0, 0, 0],
        )
    ):
        raise ValueError(
            "Negative length(s) passed. Lengths must be a value "
            "greater than 0.0. You passed {}".format(lengths)
        )

    if np.any(
        np.equal(
            lengths,
            [0, 0, 0],
        )
    ):
        if lengths[0] > 0 and lengths[1] > 0:
            warnings.warn(
                "A c value of 0 was passed. This will be "
                "interpreted as a 2-D box."
            )
        else:
            raise ValueError(
                "Length(s) of value 0 were passed. Lengths must "
                "be a value greater than 0.0. You passed "
                "{}".format(lengths)
            )
    return lengths


def _validate_angles(angles):
    """Convert angles to degree units and reshape to expected input."""
    if angles is None:
        angles = np.asarray([90, 90, 90], dtype=float, order="C")
        angles *= u.degree
    else:
        if not isinstance(angles, u.unyt_array):
            warnings.warn("Angles are assumed to be in degrees")
            angles *= u.degree

        input_unit = angles.units

        angles = np.asarray(angles, dtype=float, order="C")
        np.reshape(angles, newshape=(3, 1), order="C")

        angles *= input_unit
        angles.convert_to_units(u.degree)

    return angles


class Box(object):
    """A box that bounds a `Topology`.

    The `Box` data structure contains the relevant information to fully
    describe a simulation box in three dimensions, including lengths, angles,
    and vectors

    It is based on the Bravais lattice concept. A 3-dimensional
    prism can be fully described with 6 parameters, the lengths of the 3
    edges of the prism (`a`,`b`,`c`); and the 3 interplanar angles that
    describe the tilt of the prism edges (alpha, beta, gamma). For example,
    the Bravais parameters where a=b=c, and alpha=beta=gamma=90 define a
    cubic prism (cube).

    Parameters
    ----------
    lengths : array-like, shape(3,), dtype=float
        Lengths of the box [a, b, c]. Units are assumed to be in nm; if
        passed in as a `unyt_array` it will be converted to nm; if passed in
        as floats, nm is assumed.
    angles : array-like, optional, shape(3,), dtype=float
        Interplanar angles, [alpha, beta, gamma], that describe the box shape.
        Units are assumed to be in degrees; if passed in as a `unyt_array` it
        will be converted to degrees; if passed in as floats, degrees are
        assumed.

    Attributes
    ----------
    a, b, c : float
        Lengths of the box.
    alpha, beta, gamma : float
        Interplanar angles fully describing the `Box`.

    Methods
    -------
    get_vectors()
        Output the vectors describing the shape of the `Box`.
    get_unit_vectors()
        Output the unit vectors describing the shape of the `Box`.

    """

    def __init__(self, lengths, angles=None):
        """Construct a `Box` based on lengths and angles."""
        self._lengths = _validate_lengths(lengths)
        self._angles = _validate_angles(angles)

    @property
    def lengths(self):
        """Return edge lengths of the box."""
        return self._lengths

    @property
    def angles(self):
        """Return angles of the box."""
        return self._angles

    @lengths.setter
    def lengths(self, lengths):
        """Set the lengths of the box."""
        self._lengths = _validate_lengths(lengths)

    @angles.setter
    def angles(self, angles):
        """Set the angles of the box."""
        self._angles = _validate_angles(angles)

    def _unit_vectors_from_angles(self):
        """Return unit vectors describing prism from angles."""
        (alpha, beta, gamma) = self.angles

        cosa = np.cos(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)
        cosg = np.cos(gamma)
        sing = np.sin(gamma)
        mat_coef_y = (cosa - cosb * cosg) / sing
        mat_coef_z = np.power(sinb, 2, dtype=float) - np.power(
            mat_coef_y, 2, dtype=float
        )

        if mat_coef_z > 0.0:
            mat_coef_z = np.sqrt(mat_coef_z)
        else:
            raise Warning(
                "Non-positive z-vector. Angles {} "
                "do not generate a box with the z-vector in the"
                "positive z direction".format(self._angles)
            )

        # Note that our box vectors are always aligned along the x-axis
        #    and then the xy plane, with the z-axis
        box_vec = [[1, 0, 0], [cosg, sing, 0], [cosb, mat_coef_y, mat_coef_z]]

        return u.unyt_array(box_vec, u.dimensionless, dtype=float)

    def get_vectors(self):
        """Return the vectors of the box."""
        return (self._lengths * self.get_unit_vectors().T).T

    def get_unit_vectors(self):
        """Return the normalized vectors of the box."""
        return self._unit_vectors_from_angles()

    def json_dict(self):
        from gmso.abc.serialization_utils import unyt_to_dict

        return {
            "lengths": unyt_to_dict(self._lengths),
            "angles": unyt_to_dict(self._angles),
        }

    def __repr__(self):
        """Return formatted representation of the box."""
        return "Box(a={}, b={}, c={}, alpha={}, beta={}, gamma={})".format(
            *self._lengths, *self._angles
        )

    def __eq__(self, other):
        """Compare two boxes for equivalence."""
        if self is other:
            return True

        if not isinstance(other, Box):
            return False

        if not allclose_units(
            self.lengths, other.lengths, rtol=1e-5, atol=1e-8
        ):
            return False

        if not allclose_units(self.angles, other.angles, rtol=1e-5, atol=1e-8):
            return False

        return True
