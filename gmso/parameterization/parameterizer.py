"""The parameterizer module for a gmso Topology."""

import warnings
from typing import Dict, Union

from pydantic import Field

from gmso.abc.gmso_base import GMSOBase
from gmso.core.forcefield import ForceField
from gmso.core.topology import Topology
from gmso.exceptions import GMSOError
from gmso.parameterization.foyer_utils import (
    get_atomtyping_rules_provider,
    get_topology_graph,
    get_topology_graph_from_subtop,
    typemap_dict,
)
from gmso.parameterization.isomorph import partition_isomorphic_topology_graphs
from gmso.parameterization.subtopology_utils import (
    assert_no_boundary_bonds,
    subtop_angles,
    subtop_bonds,
    subtop_dihedrals,
    subtop_impropers,
)
from gmso.parameterization.utils import POTENTIAL_GROUPS


class GMSOParameterizationError(GMSOError):
    """Raise when parameterization fails."""


class Parameterizer(GMSOBase):
    """Utility class to parameterize a topology with gmso Forcefield."""

    topology: Topology = Field(..., description="The gmso topology.")

    forcefields: Union[ForceField, Dict[str, ForceField]] = Field(
        ...,
        description="The gmso forcefield/ a dictionary of gmso "
        "forcefields per sub-topology, where the keys "
        "should match the subtopology names",
    )

    def get_ff(self, key=None):
        """Return the forcefield of choice by looking up the forcefield dictionary."""
        if isinstance(self.forcefields, Dict):
            return self.forcefields.get(key)
        else:
            return self.forcefields

    def _parameterize_sites(self, top_or_subtop, typemap, ff):
        """Parameterize sites with appropriate atom-types from the forcefield."""
        for j, site in enumerate(top_or_subtop.sites):
            site.atom_type = ff.get_potential(
                "atom_type", typemap[j]["atom_type"]
            ).clone()  # Always properly indexed or not?

    def _parameterize_connections(
        self,
        top_or_subtop,
        ff,
        is_subtop=False,
        assert_bond_params=True,
        assert_angle_params=True,
        assert_dihedral_params=True,
        assert_improper_params=False,
    ):
        """Parameterize connections with appropriate potentials from the forcefield."""
        if is_subtop:
            bonds = subtop_bonds(top_or_subtop)
            angles = subtop_angles(top_or_subtop)
            dihedrals = subtop_dihedrals(top_or_subtop)
            impropers = subtop_impropers(top_or_subtop)
        else:
            bonds = top_or_subtop.bonds
            angles = top_or_subtop.angles
            dihedrals = top_or_subtop.dihedrals
            impropers = top_or_subtop.impropers

        self._apply_connection_parameters(bonds, ff, assert_bond_params)
        self._apply_connection_parameters(angles, ff, assert_angle_params)
        self._apply_connection_parameters(dihedrals, ff, assert_dihedral_params)
        self._apply_connection_parameters(impropers, ff, assert_improper_params)

    def _apply_connection_parameters(
        self, connections, ff, error_on_missing=True
    ):
        """Find and assign potentials from the forcefield for the provided connections."""
        visited = dict()

        for connection in connections:
            group, connection_identifiers = self.connection_identifier(
                connection
            )
            match = None
            for identifier_key in connection_identifiers:
                if tuple(identifier_key) in visited:
                    match = visited[tuple(identifier_key)]
                    break

                match = ff.get_potential(
                    group=group, key=identifier_key, warn=True
                )
                if match:
                    visited[tuple[identifier_key]] = match
                    break

            if not match and error_on_missing:
                raise GMSOParameterizationError(
                    f"No parameters found for connection {connection} in the Forcefield."
                )
            setattr(connection, group, match.clone())

    def _parameterize(
        self,
        subtop_or_top,
        typemap,
        is_subtop=False,
        assert_bond_params=True,
        assert_angle_params=True,
        assert_dihedral_params=True,
        assert_improper_params=False,
    ):
        """Parameterize a topology/subtopology based on an atomtype map."""
        forcefield = self.get_ff(subtop_or_top.name)
        self._parameterize_sites(subtop_or_top.sites, typemap, forcefield)
        self._parameterize_connections(
            subtop_or_top,
            forcefield,
            is_subtop=is_subtop,
            assert_bond_params=assert_bond_params,
            assert_angle_params=assert_angle_params,
            assert_dihedral_params=assert_dihedral_params,
            assert_improper_params=assert_improper_params,
        )

    def apply(
        self,
        identify_connected_components=True,
        use_residue_info=False,
        assert_bond_params=True,
        assert_angle_params=True,
        assert_dihedral_params=True,
        assert_improper_params=False,
    ):
        """Apply the current forcefield(s) to the topology/ various subtopologies."""
        parameterize_kwargs = dict(
            assert_bond_params=assert_bond_params,
            assert_angle_params=assert_angle_params,
            assert_dihedral_params=assert_dihedral_params,
            assert_improper_params=assert_improper_params,
        )
        if self.topology.is_typed():
            raise GMSOParameterizationError(
                "Cannot parameterize a typed topology. Please provide a topology without any types"
            )

        if isinstance(self.forcefields, Dict):
            for subtop in self.topology.subtops:
                if subtop.name not in self.forcefields:
                    warnings.warn(
                        f"Subtopology {subtop.name} will not be parameterized, as the forcefield to parameterize it "
                        f"is missing."
                    )  # FixMe: Will warning be enough?
                else:
                    assert_no_boundary_bonds(subtop)
                    typemap = self._get_atomtypes(
                        self.get_ff(subtop.name),
                        subtop,
                        identify_connected_components,
                        is_subtop=True,
                    )
                    self._parameterize(
                        typemap,
                        subtop,
                        is_subtop=True,  # This will be removed from the future iterations
                        **parameterize_kwargs,
                    )
        else:
            typemap = self._get_atomtypes(
                self.get_ff(),
                self.topology,
                identify_connected_components,
                is_subtop=False,
            )
            self._parameterize(
                self.topology,
                typemap,
                is_subtop=False,  # This will be removed from the future iterations
                **parameterize_kwargs,
            )

    @staticmethod
    def connection_identifier(
        connection,
    ):  # This can extended to incorporate a pluggable object from the forcefield.
        """Return the group and list of identifiers for a connection to query the forcefield for its potential."""
        group = POTENTIAL_GROUPS[type(connection)]
        return group, [
            list(
                member.atom_type.atomclass
                for member in connection.connection_members
            ),
            list(
                member.atom_type.name
                for member in connection.connection_members
            ),
        ]

    @staticmethod
    def _get_atomtypes(
        forcefield, topology, use_isomprohic_checks=False, is_subtop=False
    ):
        """Run atom-typing in foyer and return a typemap."""
        atom_typing_rules_provider = get_atomtyping_rules_provider(forcefield)

        if is_subtop:
            foyer_topology_graph = get_topology_graph_from_subtop(topology)
        else:
            foyer_topology_graph = get_topology_graph(topology)

        if use_isomprohic_checks:
            isomorphic_substructures = partition_isomorphic_topology_graphs(
                foyer_topology_graph
            )
            for graph, mirrors in isomorphic_substructures.items():
                typemap = typemap_dict(
                    atomtyping_rules_provider=atom_typing_rules_provider,
                    topology_graph=graph,
                )
                for mirror, mapping in mirrors:
                    for node in mirror:
                        typemap[node] = typemap[mapping[node]]
            return typemap

        else:
            return typemap_dict(
                topology_graph=foyer_topology_graph,
                atomtyping_rules_provider=atom_typing_rules_provider,
            )
