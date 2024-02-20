"""Base data structure for GMSO chemical systems."""

import itertools
import warnings
from pathlib import Path

import numpy as np
import unyt as u
from boltons.setutils import IndexedSet

import gmso
from gmso.abc.abstract_site import Site
from gmso.abc.serialization_utils import unyt_to_dict
from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.pairpotential_type import PairPotentialType
from gmso.core.views import TopologyPotentialView
from gmso.exceptions import GMSOError
from gmso.utils.connectivity import (
    identify_connections as _identify_connections,
)
from gmso.utils.conversions import (
    convert_params_units,
    convert_topology_expressions,
)
from gmso.utils.units import GMSO_UnitRegistry as UnitReg

scaling_interaction_idxes = {"12": 0, "13": 1, "14": 2}


class Topology(object):
    """A topology.

    A topology represents a chemical structure wherein lie the collection
    of sites which together form a chemical structure containing connections
    (gmso.Bond, gmso.Angle and gmso.Dihedral (along with their associated types).
    A topology is the fundamental data structure in GMSO, from which we can gather
    various information about the chemical structure and apply a forcefield
    before converting the structure into a format familiar to various simulation
    engines.

    Parameters
    ----------
    name : str, optional, default='Topology'
        A name for the Topology.
    box : gmso.Box, optional, default=None
        A gmso.Box object bounding the topology

    Attributes
    ----------
    typed : bool
        True if the topology is typed

    combining_rule : str, ['lorentz', 'geometric']
        The combining rule for the topology, can be either 'lorentz' or 'geometric'

    scaling_factors : dict
        A collection of scaling factors used in the forcefield

    n_sites : int
        Number of sites in the topology

    n_connections : int
        Number of connections in the topology (Bonds, Angles, Dihedrals, Impropers)

    n_bonds : int
        Number of bonds in the topology

    n_angles: int
        Number of angles in the topology

    n_dihedrals : int
        Number of dihedrals in the topology

    n_impropers : int
        Number of impropers in the topology

    connections : tuple of gmso.Connection objects
        A collection of bonds, angles, dihedrals, and impropers in the topology

    bonds : tuple of gmso.Bond objects
        A collection of bonds in the topology

    angles : tuple of gmso.Angle objects
        A collection of angles in the topology

    dihedrals : tuple of gmso.Dihedral objects
        A collection of dihedrals in the topology

    impropers : tuple of gmso.Improper objects
        A collection of impropers in the topology

    connection_types : tuple of gmso.Potential objects
        A collection of BondTypes, AngleTypes, DihedralTypes, and ImproperTypes in the topology

    atom_types : tuple of gmso.AtomType objects
        A collection of AtomTypes in the topology

    bond_types : tuple of gmso.BondType objects
        A collection of BondTypes in the topology

    angle_types : tuple of gmso.AngleType objects
        A collection go AngleTypes in the topology

    dihedral_types : tuple of gmso.DihedralType objects
        A collection of DihedralTypes in the topology

    improper_types : tuple of gmso.ImproperType objects
        A collection of ImproperTypes in the topology

    pairpotential_types : tuple of gmso.PairPotentialType objects
        A collection of PairPotentialTypes in the topology

    atom_type_expressions : list of gmso.AtomType.expression objects
        A collection of all the expressions for the AtomTypes in topology

    connection_type_expressions : list of gmso.Potential.expression objects
        A collection of all the expressions for the Potential objects in the topology that represent a connection type

    bond_type_expressions : list of gmso.BondType.expression objects
        A collection of all the expressions for the BondTypes in topology

    angle_type_expressions : list of gmso.AngleType.expression objects
        A collection of all the expressions for the AngleTypes in topology

    dihedral_type_expressions : list of gmso.DihedralType.expression objects
        A collection of all the expression for the DihedralTypes in the topology

    improper_type_expressions : list of gmso.ImproperType.expression objects
        A collection of all the expression for the ImproperTypes in the topology

    pairpotential_type_expressions : list of gmso.PairPotentialType.expression objects
        A collection of all the expression for the PairPotentialTypes in the topology
    """

    def __init__(self, name="Topology", box=None):
        self.name = name
        self._box = box
        self._sites = IndexedSet()
        self._typed = False
        self._bonds = IndexedSet()
        self._angles = IndexedSet()
        self._dihedrals = IndexedSet()
        self._impropers = IndexedSet()
        self._combining_rule = "lorentz"
        self._pairpotential_types = IndexedSet()
        self._scaling_factors = np.array(
            [
                [0.0, 0.0, 0.5],  # lj scales
                [0.0, 0.0, 0.5],  # electrostatics scale
            ],
            dtype=float,
        )
        self._molecule_scaling_factors = {}
        self.is_updated = True
        self._potentials_count = {
            "atom_types": 0,
            "bond_types": 0,
            "angle_types": 0,
            "dihedral_types": 0,
            "improper_types": 0,
            "pairpotential_types": 0,
        }

        self._unique_connections = {}
        self._unit_system = None

    @property
    def unit_system(self):
        """Return the unyt system of the topology."""
        return self._unit_system

    @unit_system.setter
    def unit_system(self, unit_system):
        """Set the unyt system of the topology."""
        self._name = unit_system

    @property
    def name(self):
        """Return the name of the topology."""
        return self._name

    @name.setter
    def name(self, name):
        """Set the name of the topology."""
        self._name = str(name) if name else "Topology"

    @property
    def box(self):
        """Return the Box of the topology."""
        return self._box

    @box.setter
    def box(self, box):
        """Set the box of the topology."""
        self._box = box

    @property
    def typed(self):
        """Check if the system is parametrized."""
        return self._typed

    @typed.setter
    def typed(self, typed):
        """Set the status of the topology if parametrized."""
        self._typed = typed

    @property
    def combining_rule(self):
        """Return the combining rule for the topology."""
        return self._combining_rule

    @combining_rule.setter
    def combining_rule(self, rule):
        """Set the combining rule for the topology."""
        if rule not in ["lorentz", "geometric"]:
            raise GMSOError("Combining rule must be `lorentz` or `geometric`")
        self._combining_rule = rule

    @property
    def scaling_factors(self):
        return self._scaling_factors.copy()

    @property
    def molecule_scaling_factors(self):
        return {k: v.copy() for k, v in self._molecule_scaling_factors.items()}

    @property
    def positions(self):
        """Return the positions of the sites in the topology."""
        xyz = np.empty(shape=(self.n_sites, 3)) * u.nm
        for i, site in enumerate(self._sites):
            xyz[i, :] = site.position
        return xyz

    @property
    def n_sites(self):
        """Return the number of sites in the topology."""
        return len(self._sites)

    @property
    def n_connections(self):
        """Return the number of connections in the topology."""
        return len(self.connections)

    @property
    def n_bonds(self):
        """Return the number of bonds in the topology."""
        return len(self._bonds)

    @property
    def n_angles(self):
        """Return the amount of angles in the topology."""
        return len(self._angles)

    @property
    def n_dihedrals(self):
        """Return the amount of dihedrals in the topology."""
        return len(self._dihedrals)

    @property
    def n_impropers(self):
        """Return the number of impropers in the topology."""
        return len(self._impropers)

    @property
    def sites(self):
        """Return all sites in the topology."""
        return self._sites

    @property
    def connections(self):
        """Return all connections in topology."""
        return IndexedSet(
            [*self._bonds, *self._angles, *self._dihedrals, *self._impropers]
        )

    @property
    def bonds(self):
        """Return all bonds in the topology."""
        return self._bonds

    @property
    def angles(self):
        """Return all angles in the topology."""
        return self._angles

    @property
    def dihedrals(self):
        """Return all dihedrals in the topology."""
        return self._dihedrals

    @property
    def impropers(self):
        """Return all impropers in the topology."""
        return self._impropers

    def unique_site_labels(self, label_type="molecule", name_only=False):
        """Return a list of all molecule/residue labels in the Topology."""
        # Not super happy with this method name, open for suggestion.
        unique_tags = IndexedSet()
        if name_only and label_type in ("molecule", "residue"):
            for site in self.sites:
                label = getattr(site, label_type)
                unique_tags.add(label.name if label else None)
        else:
            for site in self.sites:
                unique_tags.add(getattr(site, label_type))
        return unique_tags

    @property
    def atom_types(self):
        """Return all atom_types in the topology.

        Notes
        -----
        This returns a TopologyPotentialView object which can be used as
        an iterator. By default, this will return a view with all the atom_types
        in the topology (if multiple sites point to the same atom_type, only a
        single reference is returned/iterated upon). Use, different filters(builtin or custom) to suit your needs.
        See examples below.

        Examples
        --------
        >>> from gmso.core.atom import Atom
        >>> from gmso.core.atom_type import AtomType
        >>> from gmso.core.topology import Topology
        >>> from gmso.core.views import PotentialFilters
        >>> top = Topology(name="my_top")
        >>> atom_type = AtomType(name="my_atom_type")
        >>> for j in range(100):
        ...     atom = Atom(name=f"atom_{j+1}")
        ...     atom.atom_type = atom_type
        ...     top.add_site(atom)
        >>> len(top.atom_types)
        1
        >>> len(top.atom_types(filter_by=PotentialFilters.REPEAT_DUPLICATES))
        100
        >>> len(top.atom_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS))
        1

        See Also
        --------
        gmso.core.views.TopologyPotentialView
            An iterator/filter based view of Potentials in a gmso Topology.

        gmso.core.views.PotentialFilters
            Builtin filters for viewing potentials in a Topology.

        Returns
        -------
        gmso.core.views.TopologyPotentialView
            An iterator of the atom_types in the system filtered according to the
            filter function supplied.
        """
        return TopologyPotentialView(self._sites)

    @property
    def connection_types(self):
        """Return all connection_types in the topology.

        Notes
        -----
        This returns a TopologyPotentialView object which can be used as
        an iterator.

        See Also
        --------
        gmso.core.views.TopologyPotentialView
            An iterator/filter based view of Potentials in a gmso Topology.
        """

        return TopologyPotentialView(
            itertools.chain(
                self.bonds, self.angles, self.dihedrals, self.impropers
            )
        )

    @property
    def bond_types(self):
        """Return all bond_types in the topology.

        Notes
        -----
        This returns a TopologyPotentialView object which can be used as
        an iterator.By default, this will return a view with all the bond_types
        in the topology (if multiple bonds point to the same bond_type, only a
        single reference is returned/iterated upon). Use, different filters(builtin or custom) to suit your needs.
        See examples below.

        Examples
        --------
        >>> from gmso.core.atom import Atom
        >>> from gmso.core.bond import Bond
        >>> from gmso.core.bond_type import BondType
        >>> from gmso.core.topology import Topology
        >>> from gmso.core.views import PotentialFilters
        >>> top = Topology(name="my_top")
        >>> for j in range(100):
        ...     atom1 = Atom(name=f"atom_A_{j+1}")
        ...     atom2 = Atom(name=f"atom_B_{j+1}")
        ...     bond = Bond(connection_members=[atom1, atom2])
        ...     bond.bond_type = BondType(name=f"bond_type", member_types=('atom_A', 'atom_B'))
        ...     conn = top.add_connection(bond)
        >>> len(top.bond_types)
        100
        >>> len(top.bond_types(filter_by=PotentialFilters.UNIQUE_ID))
        100
        >>> len(top.bond_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS))
        1

        See Also
        --------
        gmso.core.views.TopologyPotentialView
            An iterator/filter based view of Potentials in a gmso Topology.

        gmso.core.views.PotentialFilters
            Builtin filters for viewing potentials in a Topology.

        Returns
        -------
        gmso.core.views.TopologyPotentialView
            An iterator of the bond_types in the system filtered according to the
            filter function supplied.
        """
        return TopologyPotentialView(self._bonds)

    @property
    def angle_types(self):
        """Return all angle_types in the topology.

        Notes
        -----
        This returns a TopologyPotentialView object which can be used as
        an iterator. By default, this will return a view with all the angle_types
        in the topology (if multiple angles point to the same angle_type, only a
        single reference is returned/iterated upon). Use, different filters(builtin or custom) to suit
        your needs. See examples below.

        Examples
        --------
        >>> from gmso.core.atom import Atom
        >>> from gmso.core.angle import Angle
        >>> from gmso.core.angle_type import AngleType
        >>> from gmso.core.topology import Topology
        >>> from gmso.core.views import PotentialFilters
        >>> for j in range(100):
        ...     atom1 = Atom(name=f"atom_A_{j+1}")
        ...     atom2 = Atom(name=f"atom_B_{j+1}")
        ...     atom3 = Atom(name=f"atom_C_{j+1}")
        ...     angle = Angle(connection_members=[atom1, atom2, atom3])
        ...     angle.angle_type = AngleType(name=f"angle_type", member_types=('atom_A', 'atom_B', 'atom_C'))
        ...     conn = top.add_connection(angle)
        >>> len(top.angle_types)
        100
        >>> len(top.angle_types(filter_by=PotentialFilters.UNIQUE_ID))
        100
        >>> len(top.angle_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS))
        1


        See Also
        --------
        gmso.core.views.TopologyPotentialView
            An iterator/filter based view of Potentials in a gmso Topology.

        gmso.core.views.PotentialFilters
            Builtin filters for viewing potentials in a Topology.

        Returns
        -------
        gmso.core.views.TopologyPotentialView
            An iterator of the angle_types in the system filtered according to the
            filter function supplied.
        """
        return TopologyPotentialView(self._angles)

    @property
    def dihedral_types(self):
        """Return all dihedral_types in the topology.

        Notes
        -----
        This returns a TopologyPotentialView object which can be used as
        an iterator. By default, this will return a view with all the dihedral_types
        in the topology (if multiple dihedrals point to the same dihedral types, only a
        single reference is returned/iterated upon). Use, different filters(builtin or custom)
        to suit your needs. See examples below.

        Examples
        --------
        >>> from gmso.core.atom import Atom
        >>> from gmso.core.dihedral import Dihedral
        >>> from gmso.core.dihedral_type import DihedralType
        >>> from gmso.core.topology import Topology
        >>> from gmso.core.views import PotentialFilters
        >>> for j in range(100):
        ...     atom1 = Atom(name=f"atom_A_{j+1}")
        ...     atom2 = Atom(name=f"atom_B_{j+1}")
        ...     atom3 = Atom(name=f"atom_C_{j+1}")
        ...     atom4 = Atom(name=f"atom_D_{j+1}")
        ...     dihedral = Dihedral(connection_members=[atom1, atom2, atom3, atom4])
        ...     dihedral.dihedral_type = DihedralType(
        ...         name=f"dihedral_type",
        ...         member_types=('atom_A', 'atom_B', 'atom_C', 'atom_D')
        ...     )
        ...     conn = top.add_connection(dihedral)
        >>> len(top.dihedral_types)
        100
        >>> len(top.dihedral_types(filter_by=PotentialFilters.UNIQUE_ID))
        100
        >>> len(top.dihedral_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS))
        1

        See Also
        --------
        gmso.core.views.TopologyPotentialView
            An iterator/filter based view of Potentials in a gmso Topology.

        gmso.core.views.PotentialFilters
            Builtin filters for viewing potentials in a Topology.

        Returns
        -------
        gmso.core.views.TopologyPotentialView
            An iterator of the dihedral_types in the system filtered according to the
            filter function supplied.
        """
        return TopologyPotentialView(self._dihedrals)

    @property
    def improper_types(self):
        """Return all improper_types in the topology.

        Notes
        -----
        This returns a TopologyPotentialView object which can be used as
        an iterator. By default, this will return a view with all the improper_types
        in the topology (if multiple impropers point to the same improper_type, only a
        single reference is returned/iterated upon). Use, different filters(builtin or custom) to
        suit your needs. See examples below.

        Examples
        --------
        >>> from gmso.core.atom import Atom
        >>> from gmso.core.improper import Improper
        >>> from gmso.core.improper_type import ImproperType
        >>> from gmso.core.topology import Topology
        >>> from gmso.core.views import PotentialFilters
        >>> for j in range(100):
        ...     atom1 = Atom(name=f"atom_A_{j+1}")
        ...     atom2 = Atom(name=f"atom_B_{j+1}")
        ...     atom3 = Atom(name=f"atom_C_{j+1}")
        ...     atom4 = Atom(name=f"atom_D_{j+1}")
        ...     improper = Improper(connection_members=[atom1, atom2, atom3, atom4])
        ...     improper.improper_type = ImproperType(
        ...         name=f"dihedral_type",
        ...         member_types=('atom_A', 'atom_B', 'atom_C', 'atom_D')
        ...     )
        ...     conn = top.add_connection(improper)
        >>> len(top.improper_types)
        100
        >>> len(top.improper_types(filter_by=PotentialFilters.UNIQUE_ID))
        100
        >>> len(top.improper_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS))
        1

        See Also
        --------
        gmso.core.views.TopologyPotentialView
            An iterator/filter based view of Potentials in a gmso Topology.

        gmso.core.views.PotentialFilters
            Builtin filters for viewing potentials in a Topology.

        Returns
        -------
        gmso.core.views.TopologyPotentialView
            An iterator of the dihedral_types in the system filtered according to the
            filter function supplied.
        """
        return TopologyPotentialView(self._impropers)

    @property
    def pairpotential_types(self):
        return self._pairpotential_types

    @property
    def atom_type_expressions(self):
        """Return all atom_type expressions in the topology."""
        return list(set([atype.expression for atype in self.atom_types]))

    @property
    def connection_type_expressions(self):
        """Return all connection_type expressions in the topology."""
        return list(
            set([contype.expression for contype in self.connection_types])
        )

    @property
    def bond_type_expressions(self):
        """Return all bond_type expressions in the topology."""
        return list(set([btype.expression for btype in self.bond_types]))

    @property
    def angle_type_expressions(self):
        """Return all angle_type expressions in the topology."""
        return list(set([atype.expression for atype in self.angle_types]))

    @property
    def dihedral_type_expressions(self):
        """Return all dihedral_type expressions in the topology."""
        return list(set([dtype.expression for dtype in self.dihedral_types]))

    @property
    def improper_type_expressions(self):
        """Return all improper_type expressions in the topology."""
        return list(set([itype.expression for itype in self.improper_types]))

    @property
    def pairpotential_type_expressions(self):
        return list(
            set([ptype.expression for ptype in self._pairpotential_types])
        )

    def get_lj_scale(self, *, molecule_id=None, interaction=None):
        """Return the selected lj_scales defined for this topology."""
        return self._get_scaling_factor(molecule_id, interaction, "lj_scale", 0)

    def set_lj_scale(self, value, *, molecule_id=None, interaction=None):
        """Set the correct lj_scaling factors for this topology."""
        self._set_scaling_factor(value, molecule_id, interaction, "lj_scale", 0)

    def get_scaling_factors(self, *, molecule_id=None):
        """Get the scaling factor of this topology or for a particular molecule"""
        return np.vstack(
            [
                self.get_lj_scale(molecule_id=molecule_id),
                self.get_electrostatics_scale(molecule_id=molecule_id),
            ]
        )

    def remove_site(self, site):
        """Remove a site from the topology.

        Parameters
        ----------
        site : gmso.core.Site
            The site to be removed.

        Notes
        -----
        When a site is removed, any connections that site belonged
        to are also removed.

        See Also
        --------
        gmso.core.topology.Topology.iter_connections_by_site
            The method that shows all connections belonging to a specific site
        """
        if site not in self._sites:
            raise ValueError(
                f"Site {site} is not currently part of this topology."
            )
        site_connections = [
            conn for conn in self.iter_connections_by_site(site)
        ]
        for conn in site_connections:
            self.remove_connection(conn)
        self._sites.remove(site)

    def remove_connection(self, connection):
        """Remove a connection from the topology.

        Parameters
        ----------
        connection : gmso.abc.abstract_conneciton.Connection
            The connection to be removed from the topology

        Notes
        -----
        The sites that belong to this connection are
        not removed from the topology.
        """
        if connection not in self.connections:
            raise ValueError(
                f"Connection {connection} is not currently part of this topology."
            )
        if isinstance(connection, gmso.core.bond.Bond):
            self._bonds.remove(connection)
        elif isinstance(connection, gmso.core.angle.Angle):
            self._angles.remove(connection)
        elif isinstance(connection, gmso.core.dihedral.Dihedral):
            self._dihedrals.remove(connection)
        elif isinstance(connection, gmso.core.improper.Improper):
            self._impropers.remove(connection)

    def set_scaling_factors(self, lj, electrostatics, *, molecule_id=None):
        """Set both lj and electrostatics scaling factors."""
        self.set_lj_scale(
            lj,
            molecule_id=molecule_id,
            interaction=None,
        )

        self.set_electrostatics_scale(
            electrostatics,
            molecule_id=molecule_id,
        )

    def get_electrostatics_scale(self, *, molecule_id=None, interaction=None):
        """Return the selected electrostatics_scale defined for this topology.

        Parameters
        ----------
        molecule_id: str, default=None
            The molecule id that this scaling factor applies to, if None
            this will return the Topology's global scaling factors

        interaction: str, one of {'12', '13', '14'}, default=None
            The interaction for which to return the scaling factor for, if None
            a 3 tuple

        Raises
        ------
        GMSOError
            If the specified parameters can't return a scaling factor
        """
        return self._get_scaling_factor(
            molecule_id, interaction, "electrostatics_scale", 1
        )

    def set_electrostatics_scale(
        self, value, *, molecule_id=None, interaction=None
    ):
        """Set the correct lj_scaling factors for this topology.

        Parameters
        ----------
        value: float, numpy.ndarray, list, or tuple of floats
            The value to set for this scale

        molecule_id: str, default=None
            The molecule id that this scaling factor applies to, if None
            this will return the Topology's global scaling factors

        interaction: str, one of {'12', '13', '14'}, default=None
            The interaction for which to return the scaling factor for, if None
            a 3 tuple

        Raises
        ------
        GMSOError
            If the specified parameters can't return a scaling factor
        """
        self._set_scaling_factor(
            value, molecule_id, interaction, "electrostatics_scale", 1
        )

    def _get_scaling_factor(self, molecule_id, interaction, name, index):
        """Get the scaling factor according to molecule_id, interaction, and name."""
        if molecule_id is None:
            all_scales = self._scaling_factors
        else:
            if molecule_id not in self._molecule_scaling_factors:
                warnings.warn(
                    f"Scaling factors for molecule `{molecule_id}` is not defined "
                    f"in the topology. Returning None."
                )
                return None
            all_scales = self._molecule_scaling_factors[molecule_id]

        if interaction is None:
            return all_scales[index].copy()
        else:
            if interaction not in scaling_interaction_idxes:
                raise GMSOError(f"Unknown `{name}` interaction `{interaction}`")
            return all_scales[index][scaling_interaction_idxes[interaction]]

    def _set_scaling_factor(self, value, molecule_id, interaction, name, index):
        """Set the scaling factor according to molecule_id, interaction, and name."""
        org_value = value
        value = np.array(value, dtype=float).reshape(-1)

        if any(np.isnan(value)):
            raise ValueError(
                f"Cannot assign a nan/NoneType to `{name}`. "
                f"Provided value: {org_value}"
            )

        if value.shape != (1,) and value.shape != (3,):
            raise ValueError(
                f"Cannot determine the appropriate shape for {org_value} to "
                f"assign it to `{name}`"
            )

        if molecule_id is None:
            all_scales = self._scaling_factors
        else:
            if molecule_id not in self._molecule_scaling_factors:
                self._molecule_scaling_factors[molecule_id] = (
                    self._scaling_factors.copy()
                )
            all_scales = self._molecule_scaling_factors[molecule_id]

        if interaction is None:
            all_scales[index] = value
        else:
            if interaction not in scaling_interaction_idxes:
                raise GMSOError(f"Unknown `{name}` interaction `{interaction}`")
            all_scales[index][scaling_interaction_idxes[interaction]] = value[0]

    def add_site(self, site, update_types=False):
        """Add a site to the topology.

        This method will add a site to the existing topology, since
        sites are stored in an indexed set, adding redundant site
        will have no effect. If the update_types parameter is set to
        true (default behavior), this method will also check if there
        is an gmso.AtomType associated with the site and it to the
        topology's AtomTypes collection.

        Parameters
        ----------
        site : gmso.core.Site
            Site to be added to this topology
        update_types : (bool), default=True
            If true, add this site's atom type to the topology's set of AtomTypes
        """
        self._sites.add(site)
        self.is_updated = False
        if update_types:
            self.update_topology()

    def add_connection(self, connection, update_types=False):
        """Add a gmso.Connection object to the topology.

        This method will add a gmso.Connection object to the
        topology, it can be used to generically include any
        Connection object i.e. Bond or Angle or Dihedral to
        the topology. According to the type of object added,
        the equivalent collection in the topology is updated.
        For example- If you add a Bond, this method will update
        topology.connections and topology.bonds object. Additionally,
        if update_types is True (default behavior), it will also
        update any Potential objects associated with the connection.

        Parameters
        ----------
        connection : one of gmso.Connection, gmso.Bond, gmso.Angle, gmso.Dihedral, or gmso.Improper object
        update_types : bool, default=True
            If True also add any Potential object associated with connection to the
            topology.

        Returns
        -------
        gmso.Connection
            The Connection object or equivalent Connection object that
            is in the topology
        """
        # Check if an equivalent connection is in the topology
        equivalent_members = connection.equivalent_members()
        if equivalent_members in self._unique_connections:
            warnings.warn(
                "An equivalent connection already exists. "
                "Providing the existing equivalent Connection."
            )
            connection = self._unique_connections[equivalent_members]

        for conn_member in connection.connection_members:
            self.add_site(conn_member)

        self._unique_connections.update({equivalent_members: connection})

        connections_sets = {
            Bond: self._bonds,
            Angle: self._angles,
            Dihedral: self._dihedrals,
            Improper: self._impropers,
        }
        connections_sets[type(connection)].add(connection)
        if update_types:
            self.update_topology()

        return connection

    def identify_connections(self):
        """Identify all connections in the topology."""
        _identify_connections(self)

    def update_atom_types(self):
        """Keep an up-to-date length of all the connection types."""
        self.update_topology()

    def update_connection_types(self):
        """Keep an up-to-date length of all the connection types."""
        self.update_topology()

    def update_topology(self):
        """Update the entire topology."""
        self._bookkeep_potentials()
        self.is_updated = True
        self.is_typed(updated=True)

    def _bookkeep_potentials(self):
        self._potentials_count = {
            "atom_types": len(self.atom_types),
            "bond_types": len(self.bond_types),
            "angle_types": len(self.angle_types),
            "dihedral_types": len(self.dihedral_types),
            "improper_types": len(self.improper_types),
            "pairpotential_types": len(self._pairpotential_types),
        }

    def add_pairpotentialtype(self, pairpotentialtype, update=True):
        """add a PairPotentialType to the topology

        This method checks whether the member_type of a PairPotentialType
        object is already stored in pairpotential_types. If so, update the
        pair potential between the member_type, and if not, add the
        PairPotentialType to the topology.

        Parameters
        ----------
        pairpotentialtype: gmso.core.PairPotentialType
            The PairPotentialType object to be added
        update: Boolean, default=True

        See Also:
        --------
        gmso.core.pairpotential_type: Pairwise potential that does not follow
        combination rules
        """
        if not isinstance(pairpotentialtype, PairPotentialType):
            raise GMSOError(
                "Non-PairPotentialType {} provided".format(pairpotentialtype)
            )
        for atype in pairpotentialtype.member_types:
            if atype not in {t.name for t in self.atom_types}:
                if atype not in {t.atomclass for t in self.atom_types}:
                    raise GMSOError(
                        "There is no name/atomclass of AtomType {} in current topology".format(
                            atype
                        )
                    )
        self._pairpotential_types.add(pairpotentialtype)

    def remove_pairpotentialtype(self, pair_of_types):
        """Remove the custom pairwise potential between two AtomTypes/Atomclasses

        Parameters
        ----------
        pair_of_types: list-like of strs
            The pair (or set) of names or atomclasses of gmso.AtomTypes of which
            the custom pairwise potential should be removed
        """
        to_delete = []
        for t in self._pairpotential_types:
            if t.member_types == tuple(pair_of_types):
                to_delete.append(t)
        if len(to_delete) > 0:
            for t in to_delete:
                self._pairpotential_types.remove(t)
        else:
            warnings.warn(
                "No pair potential specified for such pair of AtomTypes/atomclasses"
            )

    def is_typed(self, updated=False):
        """Verify if the topology is parametrized."""
        if not updated:
            self.update_topology()
        self._typed = any(self._potentials_count.values())
        return self._typed

    def is_fully_typed(self, group="topology", updated=False):
        """Check if the topology or a specifc group of objects that make up the topology are fully typed

        Parameters
        ----------
        updated : bool, optional, default=False
            If False, will update the atom type and connection type list of the
            Topology before the checking step.
        group : str, optional, default='topology'
            Specific objects to be checked. Options include:
            'topology'  : check for status of all the topology constituents
            'sites'     : check for status of all topology.sites
            'bonds'     : check for status of all topology.bonds
            'angles'    : check for status of all topology.angles
            'dihedrals' : check for status of all topology.dihedrals
            'impropers' : check for status of all topology.impropers

        Returns
        -------
        bool
            Status of the check
        Notes
        -----
        `self._type` is set to True as long as the Topology is at least
        partially typed.
        """

        if not updated:
            self.update_topology()

        typed_status = {
            "sites": lambda top: all(site.atom_type for site in top._sites),
            "bonds": lambda top: all(bond.bond_type for bond in top._bonds),
            "angles": lambda top: all(
                angle.angle_type for angle in top._angles
            ),
            "dihedrals": lambda top: all(
                dihedral.dihedral_type for dihedral in top._dihedrals
            ),
            "impropers": lambda top: all(
                improper.improper_type for improper in top._impropers
            ),
        }

        if group == "topology":
            result = list()
            for subgroup in typed_status:
                result.append(typed_status[subgroup](self))
            return all(result)
        elif group in typed_status:
            return typed_status[group](self)
        else:
            raise ValueError(
                f"Could not check typing status of {group}. "
                "Available options: 'topology', 'sites', 'bonds', "
                "'angles', 'dihedrals', 'impropers'."
            )

    def get_untyped(self, group):
        """Get the untyped (non-parametrized) objects of the Topology.

        Parameters
        ----------
        group : {'sites', 'bonds', 'angles', 'dihedrals', 'impropers', 'topology'}
            The group of objects to be checked. The 'topology' option will return
            all untyped object of the topology.

        Returns
        -------
        untyped : dict
            Dictionary of all untyped object, key of the dictionary corresponds to
            object group names define above.
        """
        untyped = dict()
        untyped_extractors = {
            "sites": self._get_untyped_sites,
            "bonds": self._get_untyped_bonds,
            "angles": self._get_untyped_angles,
            "dihedrals": self._get_untyped_dihedrals,
            "impropers": self._get_untyped_impropers,
        }
        if group == "topology":
            for subgroup in untyped_extractors:
                untyped.update(untyped_extractors[subgroup]())
        elif isinstance(group, (list, tuple, set)):
            for subgroup in group:
                untyped.update(untyped_extractors[subgroup]())
        elif isinstance(group, str) and group in untyped_extractors:
            untyped = untyped_extractors[group]()
        else:
            raise ValueError(
                f"Cannot get untyped {group}. "
                f"Available options: {[untyped_extractors.keys()]}."
            )
        return untyped

    def _get_untyped_sites(self):
        "Return a list of untyped sites"
        untyped = {"sites": list()}
        for site in self._sites:
            if not site.atom_type:
                untyped["sites"].append(site)
        return untyped

    def _get_untyped_bonds(self):
        "Return a list of untyped bonds"
        untyped = {"bonds": list()}
        for bond in self._bonds:
            if not bond.bond_type:
                untyped["bonds"].append(bond)
        return untyped

    def _get_untyped_angles(self):
        "Return a list of untyped angles"
        untyped = {"angles": list()}
        for angle in self._angles:
            if not angle.angle_type:
                untyped["angles"].append(angle)
        return untyped

    def _get_untyped_dihedrals(self):
        "Return a list of untyped dihedrals"
        untyped = {"dihedrals": list()}
        for dihedral in self._dihedrals:
            if not dihedral.dihedral_type:
                untyped["dihedrals"].append(dihedral)
        return untyped

    def _get_untyped_impropers(self):
        "Return a list of untyped impropers"
        untyped = {"impropers": list()}
        for improper in self._impropers:
            if not improper.improper_type:
                untyped["impropers"].append(improper)
        return untyped

    def update_angle_types(self):
        """Use gmso.Topology.update_connection_types to update AngleTypes in the topology.

        This method is an alias for gmso.Topology.update_connection_types.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology.
        """
        self.update_connection_types()

    def update_bond_types(self):
        """Use gmso.Topology.update_connection_types to update BondTypes in the topology.

        This method is an alias for gmso.Topology.update_connection_types.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology.
        """
        self.update_connection_types()

    def update_dihedral_types(self):
        """Use gmso.Topology.update_connection_types to update DihedralTypes in the topology.

        This method is an alias for gmso.Topology.update_connection_types.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology.
        """
        self.update_connection_types()

    def update_improper_types(self):
        """Use gmso.Topology.update_connection_types to update ImproperTypes in the topology.

        This method is an alias for gmso.Topology.update_connection_types.

        See Also
        --------
        gmso.Topology.update_connection_types :
            Update the connection types based on the connection collection in the topology.
        """
        self.update_connection_types()

    def _get_bonds_for(self, site):
        """Return a list of bonds in this Topology that the site is a part of."""
        bonds = []
        for bond in self.bonds:
            if site in bond.connection_members:
                bonds.append(bond)
        return bonds

    def _get_angles_for(self, site):
        """Return a list of angles in this Topology that the site is a part of."""
        angles = []
        for angle in self.angles:
            if site in angle.connection_members:
                angles.append(angle)
        return angles

    def _get_dihedrals_for(self, site):
        """Return a list of dihedrals in this Topology that the site is a part of."""
        dihedrals = []
        for dihedral in self.dihedrals:
            if site in dihedral.connection_members:
                dihedrals.append(dihedral)
        return dihedrals

    def get_index(self, member):
        """Get index of a member in the topology.

        Parameters
        ----------
        member : gmso Topology objects
            The member to for which to return index for.
            `member` can be of type gmso.Site, gmso.Bond, gmso.Angle, gmso.Dihedral, gmso.Improper,
            gmso.AtomType, gmso.BondType, gmso.AngleType, gmso.DihedralType or gmso.ImproperType.

        Returns
        -------
        int
            The index of the member in the topology's collection objects
        """
        refs = {
            Atom: self._sites,
            Bond: self._bonds,
            Angle: self._angles,
            Dihedral: self._dihedrals,
            Improper: self._impropers,
            AtomType: self.atom_types,
            BondType: self.bond_types,
            AngleType: self.angle_types,
            DihedralType: self.dihedral_types,
            ImproperType: self.improper_types,
            PairPotentialType: self.pairpotential_types,
        }

        member_type = type(member)

        if member_type not in refs.keys():
            raise TypeError(
                f"Cannot index member of type {member_type.__name__}"
            )

        index = refs[member_type].index(member)

        return index

    def write_forcefield(self, filename, overwrite=False):
        """Save an xml file for all parameters found in the topology.

        Parameters
        ----------
        filename: Union[str, pathlib.Path], default=None
            The filename to write the XML file to
        overwrite: bool, default=False
            If True, overwrite an existing file if it exists

        Notes
        -----
        This method can be used to save a small, trimmed down forcefield
        from a larger forcefield (e.g. oplsaa). This is useful for
        editing, saving, and sharing forcefield parameters.

        Raises
        ------
        GMSOError
            If the topology is untyped
        """
        ff = self.get_forcefield()
        ff.to_xml(filename=filename, overwrite=overwrite)

    def to_dataframe(self, parameter="sites", site_attrs=None, unyts_bool=True):
        """Return a pandas dataframe object for the sites in a topology

        Parameters
        ----------
        parameter : str, default='sites'
            A string determining what aspects of the gmso topology will be reported.
            Options are: 'sites', 'bonds', 'angles', 'dihedrals', and 'impropers'. Defaults to 'sites'.
        site_attrs : list of str, default=None
             List of strings that are attributes of the topology site and can be included as entries in the pandas dataframe.
            Examples of these can be found by printing `topology.sites[0].__dict__`.
            See https://gmso.mosdef.org/en/stable/data_structures.html#gmso.Atom for additional information on labeling.
        unyts_bool: bool, default=True
            Determine if numerical values are saved as unyt quantities or floats. See
            https://unyt.readthedocs.io/en/stable/usage.html
            for more information about manipulating unyt quantities.
            Default is True.

        Returns
        -------
        Pandas Dataframe
            A pandas.Dataframe object, see https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
            for further information.

        Examples
        ________
        >>> topology.to_dataframe(parameter = 'sites', site_attrs = ['charge'])
            This will return a dataframe with a listing of the sites and include the charges that correspond to each site.
        >>> topology.to_dataframe(parameter = 'dihedrals', site_attrs = ['positions'])
            This will return a dataframe with a listing of the sites that make up each dihedral, the positions of each of
            those sites, and the parameters that are associated with the dihedrals.

        Notes
        ____
        A dataframe is easily manipulated. In order to change the rounding to two decimals places for a column named `label`:
            >>> df['label'] = df['label'].round(2)
        The column labels can also be easily modified. This line can take a dataframe `df` and rename a column labeled
        `Atom0` to `newname` using a dictionary.
            >>> df.rename(columns = {'Atom0':'newname'})
        See https://pandas.pydata.org/pandas-docs/stable/getting_started/intro_tutorials/index.html for further information.
        """
        from gmso.utils.io import import_

        pd = import_("pandas")
        if not site_attrs:
            site_attrs = []
        df = pd.DataFrame()
        if not self.is_typed():
            raise GMSOError(
                "This topology is not typed, please type this object before converting to a pandas dataframe"
            )
        if parameter == "sites":
            df["atom_types"] = list(site.atom_type.name for site in self.sites)
            df["names"] = list(site.name for site in self.sites)
            for attr in site_attrs:
                df = self._parse_dataframe_attrs(
                    df, attr, parameter, unyts_bool
                )
        elif parameter in ["bonds", "angles", "dihedrals", "impropers"]:
            if len(getattr(self, parameter)) == 0:
                raise GMSOError(
                    f"There arent any {parameter} in the topology. The dataframe would be empty."
                )
            df = self._pandas_from_parameters(
                df,
                parameter=parameter,
                site_attrs=site_attrs,
                unyts_bool=unyts_bool,
            )
            df = self._parse_parameter_expression(df, parameter, unyts_bool)
        else:
            raise AttributeError(
                "{} is not yet supported for outputting parameters to a dataframe. \
            Please use  one of 'sites', 'bonds', 'angles', 'dihedrals', or \
            'impropers'".format(
                    str(parameter)
                )
            )

        return df

    def get_forcefield(self):
        """Get an instance of gmso.ForceField out of this topology

        Raises
        ------
        GMSOError
            If the topology is untyped
        """
        if not self.is_typed():
            raise GMSOError(
                "Cannot create a ForceField from an untyped topology."
            )
        else:
            from gmso import ForceField
            from gmso.utils._constants import FF_TOKENS_SEPARATOR

            ff = ForceField()
            ff.name = self.name + "_ForceField"
            ff.scaling_factors = {
                "electrostatics14Scale": self.scaling_factors[1, 2],
                "nonBonded14Scale": self.scaling_factors[0, 2],
            }
            for atom_type in self.atom_types:
                ff.atom_types[atom_type.name] = atom_type.copy(
                    deep=True, exclude={"topology", "set_ref"}
                )

            ff_conn_types = {
                BondType: ff.bond_types,
                AngleType: ff.angle_types,
                DihedralType: ff.dihedral_types,
                ImproperType: ff.improper_types,
            }

            for connection_type in self.connection_types:
                ff_conn_types[type(connection_type)][
                    FF_TOKENS_SEPARATOR.join(connection_type.member_types)
                ] = connection_type.copy(
                    deep=True, exclude={"topology", "set_ref"}
                )

        return ff

    def iter_sites(self, key, value):
        """Iterate through this topology's sites based on certain attribute and their values.

        Parameters
        ----------
        key: str
            The attribute of the site to look for
        value:
            The value that the given attribute should be equal to

        Yields
        ------
        gmso.abc.abstract_site.Site
            The site where getattr(site, key) == value
        """
        if key not in Site.__iterable_attributes__:
            raise ValueError(
                f"`{key}` is not an iterable attribute for Site. "
                f"To check what the iterable attributes are see gmso.abc.abstract_site module."
            )

        if value is None:
            raise ValueError(
                "Expected `value` to be something other than None. Provided None."
            )
        if key in ("molecule", "residue") and isinstance(value, str):
            for site in self._sites:
                if getattr(site, key) and getattr(site, key).name == value:
                    yield site
        for site in self._sites:
            if getattr(site, key) == value:
                yield site

    def iter_sites_by_residue(self, residue_tag):
        """Iterate through this topology's sites which contain this specific residue name.

        See Also
        --------
        gmso.core.topology.Topology.iter_sites
            The method to iterate over Topology's sites
        """
        if isinstance(residue_tag, str):
            for site in self._sites:
                if (
                    site.residue
                    and getattr(site, "residue").name == residue_tag
                ):
                    yield site
        else:
            return self.iter_sites("residue", residue_tag)

    def iter_sites_by_molecule(self, molecule_tag):
        """Iterate through this topology's sites which contain this specific molecule name.

        See Also
        --------
        gmso.core.topology.Topology.iter_sites
            The method to iterate over Topology's sites
        """
        if isinstance(molecule_tag, str):
            for site in self._sites:
                if (
                    site.molecule
                    and getattr(site, "molecule").name == molecule_tag
                ):
                    yield site
        else:
            return self.iter_sites("molecule", molecule_tag)

    def iter_connections_by_site(self, site, connections=None):
        """Iterate through this topology's connections which contain
        this specific site.

        Parameters
        ----------
        site : gmso.core.Site
            Site to limit connections search to.
        connections : set or list or tuple, optional, default=None
            The connection types to include in the search.
            If None, iterates through all of a site's connections.
            Options include "bonds", "angles", "dihedrals", "impropers"

        Yields
        ------
        gmso.abc.abstract_conneciton.Connection
            Connection where site is in Connection.connection_members

        """
        if site not in self._sites:
            raise ValueError(
                f"Site {site} is not currently part of this topology."
            )
        if connections is None:
            connections = ["bonds", "angles", "dihedrals", "impropers"]
        else:
            connections = set([option.lower() for option in connections])
            for option in connections:
                if option not in ["bonds", "angles", "dihedrals", "impropers"]:
                    raise ValueError(
                        "Valid connection types are limited to: "
                        '"bonds", "angles", "dihedrals", "impropers"'
                    )
        for conn_str in connections:
            for conn in getattr(self, conn_str):
                if site in conn.connection_members:
                    yield conn

    def create_subtop(self, label_type, label):
        """Create a new Topology object from a molecule or graup of the current Topology.

        Parameters
        ----------
        label_type: str
            Category of the label ("group" or "molecule")
        label: str (group) or tuple (molecule)
            The label of molecule or group that need to be cloned.

        Returns
        -------
        gmso.Topology
        """
        from gmso.parameterization.molecule_utils import (
            molecule_angles,
            molecule_bonds,
            molecule_dihedrals,
            molecule_impropers,
        )

        of_group = True if label_type == "group" else False
        sites_dict = {
            site: (idx, site.clone())
            for idx, site in enumerate(self.iter_sites(label_type, label))
        }
        bonds_dict = {
            bond: tuple(
                sites_dict[bond.connection_members[i]][0] for i in range(2)
            )
            for bond in molecule_bonds(self, label, of_group)
        }

        angles_dict = {
            angle: tuple(
                sites_dict[angle.connection_members[i]][0] for i in range(3)
            )
            for angle in molecule_angles(self, label, of_group)
        }

        dihedrals_dict = {
            dihedral: tuple(
                sites_dict[dihedral.connection_members[i]][0] for i in range(4)
            )
            for dihedral in molecule_dihedrals(self, label, of_group)
        }

        impropers_dict = {
            improper: tuple(
                sites_dict[improper.connection_members[i]][0] for i in range(4)
            )
            for improper in molecule_impropers(self, label, of_group)
        }

        new_top = gmso.Topology(
            name=label if isinstance(label, str) else label[0]
        )

        for ref_site, new_site in sites_dict.items():
            new_top.add_site(new_site[1])
        for ref_conn, conn_idx in bonds_dict.items():
            bond = gmso.Bond(
                connection_members=[
                    new_top.sites[conn_idx[i]] for i in range(2)
                ],
                bond_type=(
                    None
                    if not ref_conn.connection_type
                    else ref_conn.connection_type.clone()
                ),
            )
            new_top.add_connection(bond)
        for ref_conn, conn_idx in angles_dict.items():
            angle = gmso.Angle(
                connection_members=[
                    new_top.sites[conn_idx[i]] for i in range(3)
                ],
                angle_type=(
                    None
                    if not ref_conn.connection_type
                    else ref_conn.connection_type.clone()
                ),
            )
            new_top.add_connection(angle)
        for ref_conn, conn_idx in dihedrals_dict.items():
            dihedral = gmso.Dihedral(
                connection_members=[
                    new_top.sites[conn_idx[i]] for i in range(4)
                ],
                dihedral_type=(
                    None
                    if not ref_conn.connection_type
                    else ref_conn.connection_type.clone()
                ),
            )
            new_top.add_connection(dihedral)
        for ref_conn, conn_idx in impropers_dict.items():
            improper = gmso.Improper(
                connection_members=[
                    new_top.sites[conn_idx[i]] for i in range(4)
                ],
                improper_type=(
                    None
                    if not ref_conn.connection_type
                    else ref_conn.connection_type.clone()
                ),
            )
            new_top.add_connection(improper)

        new_top.update_topology()
        return new_top

    def save(self, filename, overwrite=False, **kwargs):
        """Save the topology to a file.

        Parameters
        ----------
        filename: str, pathlib.Path
            The file to save the topology as
        overwrite: bool, default=True
            If True, overwrite the existing file if it exists
        **kwargs:
            The arguments to specific file savers listed below(as extensions):
            * json: types, update, indent
            * gro: precision
            * lammps/lammpsdata: atom_style
        """
        if not isinstance(filename, Path):
            filename = Path(filename).resolve()

        if filename.exists() and not overwrite:
            raise FileExistsError(
                f"The file {filename} exists. Please set "
                f"overwrite=True if you wish to overwrite the existing file"
            )

        from gmso.formats import SaversRegistry

        saver = SaversRegistry.get_callable(filename.suffix)
        saver(self, filename, **kwargs)

    def __repr__(self):
        """Return custom format to represent topology."""
        if not self.is_updated:
            self.update_topology()
        return (
            f"<Topology {self.name}, {self.n_sites} sites,\n "
            f"{self.n_connections} connections,\n "
            f"{sum(self._potentials_count.values())} potentials,\n "
            f"id: {id(self)}>"
        )

    def __str__(self):
        """Return custom format to represent topology as a string."""
        return f"<Topology {self.name}, {self.n_sites} sites, id: {id(self)}>"

    def _pandas_from_parameters(
        self, df, parameter, site_attrs=None, unyts_bool=True
    ):
        """Add to a pandas dataframe the site indices for each connection member in a
        multimember topology attribute such as a bond. Also include information about
        those sites in the site_attrs list"""
        if site_attrs is None:
            site_attrs = []
        sites_per_connection = len(
            getattr(self, parameter)[0].connection_members
        )
        for site_index in np.arange(sites_per_connection):
            df["Atom" + str(site_index)] = list(
                str(connection.connection_members[site_index].name)
                + f"({self.get_index(connection.connection_members[site_index])})"
                for connection in getattr(self, parameter)
            )
        for attr in site_attrs:
            df = self._parse_dataframe_attrs(
                df, attr, parameter, sites_per_connection, unyts_bool
            )
        return df

    def _parse_dataframe_attrs(
        self, df, attr, parameter, sites_per_connection=1, unyts_bool=True
    ):
        """Parses an attribute string to correctly format and return the topology attribute
        into a pandas dataframe"""
        if parameter == "sites":
            if "." in attr:
                try:
                    attr1, attr2 = attr.split(".")
                    df[attr] = list(
                        _return_float_for_unyt(
                            getattr(getattr(site, attr1), attr2),
                            unyts_bool,
                        )
                        for site in self.sites
                    )
                except AttributeError:
                    raise AttributeError(
                        f"The attribute {attr} is not in this gmso object."
                    )
            elif attr == "positions" or attr == "position":
                for i, dimension in enumerate(["x", "y", "z"]):
                    df[dimension] = list(
                        _return_float_for_unyt(
                            getattr(site, "position")[i], unyts_bool
                        )
                        for site in self.sites
                    )
            elif attr == "charge" or attr == "charges":
                df["charge (e)"] = list(
                    site.charge.in_units(
                        u.Unit(
                            "elementary_charge", registry=UnitReg.default_reg()
                        )
                    ).to_value()
                    for site in self.sites
                )
            else:
                try:
                    df[attr] = list(
                        _return_float_for_unyt(getattr(site, attr), unyts_bool)
                        for site in self.sites
                    )
                except AttributeError:
                    raise AttributeError(
                        f"The attribute {attr} is not in this gmso object."
                    )

        elif parameter in ["bonds", "angles", "dihedrals", "impropers"]:
            for site_index in np.arange(sites_per_connection):
                if "." in attr:
                    try:
                        attr1, attr2 = attr.split(".")
                        df[attr + " Atom" + str(site_index)] = list(
                            _return_float_for_unyt(
                                getattr(
                                    getattr(
                                        connection.connection_members[
                                            site_index
                                        ],
                                        attr1,
                                    ),
                                    attr2,
                                ),
                                unyts_bool,
                            )
                            for connection in getattr(self, parameter)
                        )
                    except AttributeError:
                        raise AttributeError(
                            f"The attribute {attr} is not in this gmso object."
                        )
                elif attr == "positions" or attr == "position":
                    df["x Atom" + str(site_index) + " (nm)"] = list(
                        _return_float_for_unyt(
                            getattr(
                                connection.connection_members[site_index],
                                "position",
                            )[0],
                            unyts_bool,
                        )
                        for connection in getattr(self, parameter)
                    )
                    df["y Atom" + str(site_index) + " (nm)"] = list(
                        _return_float_for_unyt(
                            getattr(
                                connection.connection_members[site_index],
                                "position",
                            )[1],
                            unyts_bool,
                        )
                        for connection in getattr(self, parameter)
                    )
                    df["z Atom" + str(site_index) + " (nm)"] = list(
                        _return_float_for_unyt(
                            getattr(
                                connection.connection_members[site_index],
                                "position",
                            )[2],
                            unyts_bool,
                        )
                        for connection in getattr(self, parameter)
                    )
                elif attr == "charge" or attr == "charges":
                    df["charge Atom" + str(site_index) + " (e)"] = list(
                        getattr(
                            connection.connection_members[site_index],
                            "charge",
                        )
                        .in_units(
                            u.Unit(
                                "elementary_charge",
                                registry=UnitReg.default_reg(),
                            )
                        )
                        .value
                        for connection in getattr(self, parameter)
                    )
                else:
                    try:
                        df[f"{attr} Atom {site_index}"] = list(
                            _return_float_for_unyt(
                                getattr(
                                    connection.connection_members[site_index],
                                    attr,
                                ),
                                unyts_bool,
                            )
                            for connection in getattr(self, parameter)
                        )
                    except AttributeError:
                        raise AttributeError(
                            f"The attribute {attr} is not in this gmso object."
                        )
        else:
            raise AttributeError(
                f"{parameter} is not yet supported for adding labels to a dataframe. \
                 Please use  one of 'sites', 'bonds', 'angles', 'dihedrals', or 'impropers'"
            )
        return df

    def _parse_parameter_expression(self, df, parameter, unyts_bool):
        """Take a given topology attribute and return the parameters associated with it"""
        for i, param in enumerate(
            getattr(
                getattr(self, parameter)[0], parameter[:-1] + "_type"
            ).parameters
        ):
            df[
                f"Parameter {i} ({param}) {getattr(getattr(self, parameter)[0], parameter[:-1]+'_type').parameters[param].units}"
            ] = list(
                _return_float_for_unyt(
                    getattr(connection, parameter[:-1] + "_type").parameters[
                        param
                    ],
                    unyts_bool,
                )
                for connection in getattr(self, parameter)
            )
        return df

    @classmethod
    def load(cls, filename, **kwargs):
        """Load a file to a topology"""
        filename = Path(filename).resolve()
        from gmso.formats import LoadersRegistry

        loader = LoadersRegistry.get_callable(filename.suffix)
        return loader(filename, **kwargs)

    def convert_potential_styles(self, expressionMap={}):
        """Convert from one parameter form to another.

        Parameters
        ----------
        expressionMap : dict, default={}
            Map where the keys represent the current potential
            type and the corresponding values represent the desired
            potential type. The desired potential style can be
            either a string with the corresponding name, or
            a gmso.utils.expression.PotentialExpression type.

        Examples
        ________
        # Convert from RB torsions to OPLS torsions
        top.convert_potential_styles({"dihedrals": "OPLSTorsionPotential"})
        # TODO: convert_potential_styles with PotentialExpression
        """
        # TODO: raise warnings for improper values or keys in expressionMap

        return convert_topology_expressions(self, expressionMap)

    def convert_unit_styles(self, unitsystem, exp_unitsDict):
        """Convert from one set of base units to another.

        Parameters
        ----------
        unitsystem : unyt.UnitSystem
            set of base units to use for all expressions of the topology
            in `unyt package <https://unyt.readthedocs.io/en/stable/>_`
        exp_unitsDict : dict
            keys with topology attributes that should be converted and
            values with dictionary of parameter: expected_dimension

        Examples
        ________
        top.convert_unit_styles(
            u.UnitSystem(
                "lammps_real", "Å", "amu", "fs", "K", "rad",
            ),
            {"bond":{"k":"energy/length**2", "r_eq":"length"}},
        )
        """

        ref_values = {"energy": "kJ/mol", "length": "nm", "angle": "radians"}

        # all potContainer ["atom", "bond", "angle", "dihedral", "improper"]
        for potStr in exp_unitsDict:
            potContainer = getattr(self, potStr + "_types")
            convert_params_units(
                potContainer,
                expected_units_dim=exp_unitsDict[potStr],
                base_units=unitsystem,
                ref_values=ref_values,
            )


def _return_float_for_unyt(unyt_quant, unyts_bool):
    try:
        return unyt_quant if unyts_bool else unyt_to_dict(unyt_quant)["array"]
    except TypeError:
        return unyt_quant
