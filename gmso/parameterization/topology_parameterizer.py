"""The parameterizer module for a gmso Topology."""

import warnings
from typing import Dict, Union

import networkx as nx
from boltons.setutils import IndexedSet
from pydantic import Field

from gmso.abc.gmso_base import GMSOBase
from gmso.core.forcefield import ForceField
from gmso.core.topology import Topology
from gmso.exceptions import GMSOError
from gmso.parameterization.foyer_utils import (
    get_atomtyping_rules_provider,
    get_topology_graph,
    typemap_dict,
)
from gmso.parameterization.isomorph import (
    partition_isomorphic_topology_graphs,
    top_node_match,
)
from gmso.parameterization.molecule_utils import (
    assert_no_boundary_bonds,
    molecule_angles,
    molecule_bonds,
    molecule_dihedrals,
    molecule_impropers,
)
from gmso.parameterization.utils import POTENTIAL_GROUPS


class ParameterizationError(GMSOError):
    """Raise when parameterization fails."""


class TopologyParameterizationConfig(GMSOBase):
    """Configuration options for parameterizing a topology."""

    clone_topology: bool = Field(
        default=False,
        description="If true, clone the topology and apply parameters to the cloned one.",
    )  # Unused

    identify_connections: bool = Field(
        default=False,
        description="If true, add connections identified using networkx graph matching to match"
        "the topology's bonding graph to smaller sub-graphs that correspond to an "
        "angle, dihedral, improper etc",
    )

    identify_connected_components: bool = Field(
        default=False,
        description="A flag to determine whether or not to search the topology"
        " for repeated disconnected structures, otherwise known as "
        "molecules and type each molecule only once.",
    )

    use_molecule_info: bool = Field(
        default=False,
        description="A flag to determine whether or not to look at site.molecule "
        "to look parameterize each molecule only once. Will only be used if "
        "identify_connected_components=False",
    )  # Unused

    assert_bond_params: bool = Field(
        default=True,
        description="If True, an error is raised if parameters are not found for "
        "all system bonds.",
    )

    assert_angle_params: bool = Field(
        default=True,
        description="If True, an error is raised if parameters are not found for "
        "all system angles",
    )

    assert_dihedral_params: bool = (
        Field(
            default=True,
            description="If True, an error is raised if parameters are not found for "
            "all system dihedrals.",
        ),
    )

    assert_improper_params: bool = Field(
        default=False,
        description="If True, an error is raised if parameters are not found for "
        "all system impropers.",
    )


class TopologyParameterizer(GMSOBase):
    """Utility class to parameterize a topology with gmso Forcefield."""

    topology: Topology = Field(..., description="The gmso topology.")

    forcefields: Union[ForceField, Dict[str, ForceField]] = Field(
        ...,
        description="The gmso forcefield/ a dictionary of gmso "
        "forcefields per molecule/group, where the keys "
        "should match the molecule/group names",
    )

    config: TopologyParameterizationConfig = Field(
        ..., description="The configuration options for the parameterizer."
    )

    def get_ff(self, key=None):
        """Return the forcefield of choice by looking up the forcefield dictionary."""
        if isinstance(self.forcefields, Dict):
            return self.forcefields.get(key)
        else:
            return self.forcefields

    def _parameterize_sites(self, sites, typemap, ff, use_molecule_info=None):
        """Parameterize sites with appropriate atom-types from the forcefield."""
        for j, site in enumerate(sites):
            site.atom_type = ff.get_potential(
                "atom_type", typemap[j]["atomtype"]
            ).clone()  # Always properly indexed or not?

    def _parameterize_connections(
        self,
        top,
        ff,
        of_group=None,
    ):
        """Parameterize connections with appropriate potentials from the forcefield."""
        if of_group:
            bonds = molecule_bonds(top, of_group)
            angles = molecule_angles(top, of_group)
            dihedrals = molecule_dihedrals(top, of_group)
            impropers = molecule_impropers(top, of_group)
        else:
            bonds = top.bonds
            angles = top.angles
            dihedrals = top.dihedrals
            impropers = top.impropers

        self._apply_connection_parameters(
            bonds, ff, self.config.assert_bond_params
        )
        self._apply_connection_parameters(
            angles, ff, self.config.assert_angle_params
        )
        self._apply_connection_parameters(
            dihedrals, ff, self.config.assert_dihedral_params
        )
        self._apply_connection_parameters(
            impropers, ff, self.config.assert_improper_params
        )

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
                    visited[tuple(identifier_key)] = match
                    break

            if not match and error_on_missing:
                raise ParameterizationError(
                    f"No parameters found for connection {connection}, group: {group}, "
                    f"identifiers: {connection_identifiers} in the Forcefield."
                )
            elif match:
                setattr(connection, group, match.clone())

    def _parameterize(
        self, top, typemap, of_group=None, use_molecule_info=False
    ):
        """Parameterize a topology/subtopology based on an atomtype map."""
        if of_group:
            forcefield = self.get_ff(of_group.name)
            sites = top.iter_sites("molecule", of_group)
        else:
            forcefield = self.get_ff(top.name)
            sites = top.sites

        self._parameterize_sites(
            sites, typemap, forcefield, use_molecule_info=use_molecule_info
        )
        self._parameterize_connections(
            top,
            forcefield,
            of_group=of_group,
        )

    def _verify_forcefields_metadata(self):
        """Verify all the provided forcefields have the same scaling factors and combining rule."""
        if isinstance(self.forcefields, dict):
            ffs = list(self.forcefields.values())
            init_scaling_factors = ffs[0].scaling_factors
            init_combining_rule = ffs[0].combining_rule
            for ff in ffs[1:]:
                if ff.scaling_factors != init_scaling_factors:
                    raise ParameterizationError(
                        "Scaling factors of the provided forcefields do not"
                        "match, please provide forcefields with same scaling"
                        "factors that apply to a Topology"
                    )

                if ff.combining_rule != init_combining_rule:
                    raise ParameterizationError(
                        "Combining rules of the provided forcefields do not"
                        "match, please provide forcefields with same scaling"
                        "factors that apply to a Topology"
                    )
            return init_scaling_factors, init_combining_rule
        else:
            return (
                self.forcefields.scaling_factors,
                self.forcefields.combining_rule,
            )

    def run_parameterization(self):
        """Run parameterization of the topology with give forcefield(s) and configuration."""
        scaling_factors, combining_rule = self._verify_forcefields_metadata()
        if self.topology.is_typed():
            raise ParameterizationError(
                "Cannot parameterize a typed topology. Please provide a topology without any types"
            )

        if self.config.identify_connections:
            """ToDo: This mutates the topology and is agnostic to downstream
            errors. So, here we should use index only option"""
            self.topology.identify_connections()

        if isinstance(self.forcefields, Dict):
            group_labels = self.topology.unique_site_labels("molecule")
            if not group_labels or group_labels == IndexedSet([None]):
                raise ParameterizationError(
                    f"The provided gmso topology doesn't have any group."
                    f"Either use a single forcefield to apply to to whole topology "
                    f"or provide an appropriate topology whose molecule names are "
                    f"the keys of the `forcefields` dictionary. Provided Forcefields: "
                    f"{self.forcefields}, Topology: {self.topology}"
                )
            assert_no_boundary_bonds(self.topology)
            for group in group_labels:
                if group.name not in self.forcefields:
                    warnings.warn(
                        f"Group {group.name} will not be parameterized, as the forcefield to parameterize it "
                        f"is missing."
                    )  # FixMe: Will warning be enough?
                else:
                    typemap = self._get_atomtypes(
                        self.get_ff(group.name),
                        self.topology,
                        self.config.identify_connected_components,
                        of_group=group,
                        use_molecule_info=self.config.use_molecule_info,
                    )
                    self._parameterize(
                        self.topology,
                        typemap,
                        of_group=group,
                        use_molecule_info=self.config.use_molecule_info,  # This will be removed from the future iterations
                    )
        else:
            typemap = self._get_atomtypes(
                self.get_ff(),
                self.topology,
                self.config.identify_connected_components,
                use_molecule_info=self.config.use_molecule_info,
            )
            self._parameterize(
                self.topology,
                typemap,
                use_molecule_info=self.config.use_molecule_info,
            )

        self.topology.scaling_factors.update(scaling_factors)
        self.topology.combining_rule = combining_rule
        self.topology.update_topology()

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
        forcefield,
        topology,
        use_isomorphic_checks=False,
        of_group=None,
        use_molecule_info=False,
    ):
        """Run atom-typing in foyer and return the typemap."""
        #if use_isomorphic_checks and use_molecule_info:
        #    raise GMSOError(
        #        "Cannot set `use_isomorphic_checks=True` "
        #        "and `use_molecule_info=True` at the same time."
        #    )
        atom_typing_rules_provider = get_atomtyping_rules_provider(forcefield)
        foyer_topology_graph = get_topology_graph(topology, of_group)

        if use_isomorphic_checks:
            isomorphic_substructures = partition_isomorphic_topology_graphs(
                foyer_topology_graph
            )
            typemap = {}
            for graph, mirrors in isomorphic_substructures.items():
                typemap.update(
                    typemap_dict(
                        atomtyping_rules_provider=atom_typing_rules_provider,
                        topology_graph=graph,
                    )
                )
                for mirror, mapping in mirrors:
                    for node in mirror:
                        typemap[node] = typemap[mapping[node]]
            return typemap

        elif use_molecule_info:
            typemap = {}
            reference = dict()
            for connected_component in nx.connected_components(
                foyer_topology_graph
            ):
                subgraph = foyer_topology_graph.subgraph(connected_component)
                nodes_idx = tuple(subgraph.nodes)
                molecule = subgraph.nodes[nodes_idx[0]]["atom_data"].molecule
                if molecule.name not in reference:
                    reference[molecule.name] = {
                        "typemap": typemap_dict(
                            atomtyping_rules_provider=atom_typing_rules_provider,
                            topology_graph=subgraph,
                        ),
                        "graph": subgraph,
                    }
                    typemap.update(reference[molecule.name]["typemap"])
                else:
                    if use_isomorphic_checks:
                        # Check for isomorphism submatching to typemap
                        matcher = nx.algorithms.isomorphism.GraphMatcher(
                            subgraph,
                            reference[molecule.name]["graph"],
                            node_match=top_node_match,
                        )
                        assert matcher.is_isomorphic()
                        for node in subgraph.nodes:
                            typemap[node] = reference[molecule.name]["typemap"][
                                matcher.mapping[node]
                            ]
                    else:
                        # Assume molecule sites are in the same order
                        for node in subgraph.nodes:
                            typemap[node] = reference[molecule.name]["typemap"][node]

            return typemap
        else:
            foyer_topology_graph = get_topology_graph(topology, of_group)
            return typemap_dict(
                topology_graph=foyer_topology_graph,
                atomtyping_rules_provider=atom_typing_rules_provider,
            )
