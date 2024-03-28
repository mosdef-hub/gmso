"""The parameterizer module for a gmso Topology."""

import warnings
from typing import Dict, Union

import networkx as nx
from boltons.setutils import IndexedSet

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

try:
    from pydantic.v1 import Field
except ImportError:
    from pydantic import Field


class ParameterizationError(GMSOError):
    """Raise when parameterization fails."""


class TopologyParameterizationConfig(GMSOBase):
    """Configuration options for parameterizing a topology."""

    clone_topology: bool = Field(
        default=False,
        description="If true, clone the topology and apply parameters to the cloned one.",
    )  # Unused

    match_ff_by: str = Field(
        default=None,
        description="The site's' label used to matched with the provided dictionary.",
    )

    identify_connections: bool = Field(
        default=False,
        description="If true, add connections identified using networkx graph matching to match"
        "the topology's bonding graph to smaller sub-graphs that correspond to an "
        "angle, dihedral, improper etc",
    )

    speedup_by_molgraph: bool = Field(
        default=False,
        description="A flag to determine whether or not to search the topology"
        " for repeated disconnected structures, otherwise known as "
        "molecules and type each molecule only once.",
    )

    speedup_by_moltag: bool = Field(
        default=False,
        description="A flag to determine whether or not to look at site.molecule "
        "to look parameterize each molecule only once. Will only be used if "
        "speedup_by_molgraph=True",
    )

    ignore_params: list = Field(
        default=[],
        description="Skipping the checks that make sure all connections (in the list) "
        "have a connection types.",
    )

    remove_untyped: bool = Field(
        default=False,
        description="If True, after the atomtyping and parameterization step, "
        "remove all connection that has no connection_type",
    )

    fast_copy: bool = Field(
        default=True,
        description="If True, don't deepcopy sympy expression and sympy independent, "
        "variables to save time on parameterization step.",
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

    def _parameterize_sites(self, sites, typemap, ff, speedup_by_moltag=None):
        """Parameterize sites with appropriate atom-types from the forcefield."""
        for j, site in enumerate(sites):
            site.atom_type = ff.get_potential(
                "atom_type", typemap[j]["atomtype"]
            ).clone(self.config.fast_copy)
            assert site.atom_type, site

    def _parameterize_connections(
        self,
        top,
        ff,
        label_type=None,
        label=None,
    ):
        """Parameterize connections with appropriate potentials from the forcefield."""
        if label_type and label:
            bonds = molecule_bonds(
                top, label, True if label_type == "group" else False
            )
            angles = molecule_angles(
                top, label, True if label_type == "group" else False
            )
            dihedrals = molecule_dihedrals(
                top, label, True if label_type == "group" else False
            )
            impropers = molecule_impropers(
                top, label, True if label_type == "group" else False
            )
        else:
            bonds = top.bonds
            angles = top.angles
            dihedrals = top.dihedrals
            impropers = top.impropers

        self._apply_connection_parameters(
            bonds, ff, False if "bond" in self.config.ignore_params else True
        )
        self._apply_connection_parameters(
            angles, ff, False if "angle" in self.config.ignore_params else True
        )
        self._apply_connection_parameters(
            dihedrals,
            ff,
            False if "dihedral" in self.config.ignore_params else True,
        )
        self._apply_connection_parameters(
            impropers,
            ff,
            False if "improper" in self.config.ignore_params else True,
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
                    group=group,
                    key=identifier_key,
                    return_match_order=True,
                    warn=True,
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
                setattr(
                    connection, group, match[0].clone(self.config.fast_copy)
                )
                matched_order = [
                    connection.connection_members[i] for i in match[1]
                ]
                connection.connection_members = matched_order
                if not match[0].member_types:
                    connection.connection_type.member_types = tuple(
                        member.atom_type.name for member in matched_order
                    )
                if not match[0].member_classes:
                    connection.connection_type.member_classes = tuple(
                        member.atom_type.atomclass for member in matched_order
                    )

    def _parameterize(
        self, top, typemap, label_type=None, label=None, speedup_by_moltag=False
    ):
        """Parameterize a topology/subtopology based on an atomtype map."""
        if label and label_type:
            forcefield = self.get_ff(label)
            sites = top.iter_sites(label_type, label)
        else:
            forcefield = self.get_ff(top.name)
            sites = top.sites

        self._parameterize_sites(
            sites, typemap, forcefield, speedup_by_moltag=speedup_by_moltag
        )
        self._parameterize_connections(
            top,
            forcefield,
            label_type,
            label,
        )

    def _set_combining_rule(self):
        """Verify all the provided forcefields have the same combining rule and set it for the Topology."""
        if isinstance(self.forcefields, dict):
            all_comb_rules = set(
                ff.combining_rule for ff in self.forcefields.values()
            )
        else:
            all_comb_rules = {self.forcefields.combining_rule}

        if not len(all_comb_rules) == 1:
            raise ParameterizationError(
                "Combining rules of the provided forcefields do not"
                "match, please provide forcefields with same scaling"
                "factors that apply to a Topology"
            )
        self.topology.combining_rule = all_comb_rules.pop()

    def _set_scaling_factors(self):
        """Set either per-molecule or global scaling factors for the topology based on the forcefields provided."""
        # ToDo: Set other scaling factors by extending the forcefield schema
        # ToDo: What to do when all the scaling factors matchup? Should we promote them to be global?
        # ToDo: Do we want to also parse other interaction if provided?
        lj_scales = {
            f"nonBonded{interaction}Scale": interaction
            for interaction in ["12", "13", "14"]
        }
        electrostatics_scales = {
            f"electrostatics{interaction}Scale": interaction
            for interaction in ["12", "13", "14"]
        }

        if isinstance(self.forcefields, Dict):
            for group_or_molecule, ff in self.forcefields.items():
                for name, interaction in lj_scales.items():
                    if ff.scaling_factors.get(name) is not None:
                        self.topology.set_lj_scale(
                            ff.scaling_factors[name],
                            interaction=interaction,
                            molecule_id=group_or_molecule,
                        )
                for name, interaction in electrostatics_scales.items():
                    if ff.scaling_factors.get(name) is not None:
                        self.topology.set_electrostatics_scale(
                            ff.scaling_factors[name],
                            interaction=interaction,
                            molecule_id=group_or_molecule,
                        )
        else:
            for name, interaction in lj_scales.items():
                if self.forcefields.scaling_factors.get(name) is not None:
                    self.topology.set_lj_scale(
                        self.forcefields.scaling_factors[name],
                        interaction=interaction,
                    )
            for name, interaction in electrostatics_scales.items():
                if self.forcefields.scaling_factors.get(name) is not None:
                    self.topology.set_electrostatics_scale(
                        self.forcefields.scaling_factors[name],
                        interaction=interaction,
                    )

    def run_parameterization(self):
        """Run parameterization of the topology with give forcefield(s) and configuration."""
        if self.topology.is_typed():
            raise ParameterizationError(
                "Cannot parameterize a typed topology. Please provide a topology without any types"
            )

        self._set_combining_rule()  # Fail Early if no match

        if self.config.identify_connections:
            """ToDo: This mutates the topology and is agnostic to downstream
            errors. So, here we should use index only option"""
            self.topology.identify_connections()

        if isinstance(self.forcefields, Dict):
            labels = self.topology.unique_site_labels(
                self.config.match_ff_by, name_only=True
            )
            if not labels or labels == IndexedSet([None]):
                # raise ParameterizationError(
                warnings.warn(
                    f"The provided gmso topology doesn't have any group/molecule."
                    f"Either use a single forcefield to apply to to whole topology "
                    f"or provide an appropriate topology whose molecule names are "
                    f"the keys of the `forcefields` dictionary. Provided Forcefields: "
                    f"{self.forcefields}, Topology: {self.topology}"
                )

            assert_no_boundary_bonds(self.topology)
            for label in labels:
                if label not in self.forcefields:
                    warnings.warn(
                        f"Group/molecule {label} will not be parameterized, as the forcefield to parameterize it "
                        f"is missing."
                    )  # FixMe: Will warning be enough?
                else:
                    typemap = self._get_atomtypes(
                        self.get_ff(label),
                        self.topology,
                        self.config.match_ff_by,
                        label,
                        self.config.speedup_by_moltag,
                        self.config.speedup_by_molgraph,
                    )
                    self._parameterize(
                        self.topology,
                        typemap,
                        label_type=self.config.match_ff_by,
                        label=label,
                        speedup_by_moltag=self.config.speedup_by_moltag,  # This will be removed from the future iterations
                    )
        else:
            typemap = self._get_atomtypes(
                self.get_ff(),
                self.topology,
                speedup_by_moltag=self.config.speedup_by_moltag,
                use_isomorphic_checks=self.config.speedup_by_molgraph,
            )
            self._parameterize(
                self.topology,
                typemap,
                speedup_by_moltag=self.config.speedup_by_moltag,
            )

        self._set_scaling_factors()  # Set global or per molecule scaling factors
        self.topology.update_topology()

        if self.config.remove_untyped:
            # TEMP CODE: copied from foyer/general_forcefield.py, will update later
            for i in range(self.topology.n_bonds - 1, -1, -1):
                if not self.topology.bonds[i].bond_type:
                    self.topology._bonds.pop(i)
            for i in range(self.topology.n_angles - 1, -1, -1):
                if not self.topology.angles[i].angle_type:
                    self.topology._angles.pop(i)
            for i in range(self.topology.n_dihedrals - 1, -1, -1):
                if not self.topology.dihedrals[i].dihedral_type:
                    self.topology._dihedrals.pop(i)
            for i in range(self.topology.n_impropers - 1, -1, -1):
                if not self.topology.impropers[i].improper_type:
                    self.topology._impropers.pop(i)

    @staticmethod
    def connection_identifier(
        connection,
    ):  # This can extended to incorporate a pluggable object from the forcefield.
        """Return the group and list of identifiers for a connection to query the forcefield for its potential."""
        group = POTENTIAL_GROUPS[type(connection)]
        return group, [
            list(
                member.atom_type.name
                for member in connection.connection_members
            ),
            list(
                member.atom_type.atomclass
                for member in connection.connection_members
            ),
        ]

    @staticmethod
    def _get_atomtypes(
        forcefield,
        topology,
        label_type=None,
        label=None,
        speedup_by_moltag=False,
        use_isomorphic_checks=False,
    ):
        """Run atom-typing in foyer and return the typemap."""
        atom_typing_rules_provider = get_atomtyping_rules_provider(forcefield)
        foyer_topology_graph = get_topology_graph(
            topology,
            label_type,
            label,
        )

        if speedup_by_moltag:
            # Iterate through foyer_topology_graph, which is a subgraph of label_type
            typemap, reference = dict(), dict()
            for connected_component in nx.connected_components(
                foyer_topology_graph
            ):
                subgraph = foyer_topology_graph.subgraph(connected_component)
                nodes_idx = tuple(subgraph.nodes)
                molecule = subgraph.nodes[nodes_idx[0]]["atom_data"].molecule
                if molecule not in reference:
                    reference[molecule] = {
                        "typemap": typemap_dict(
                            atomtyping_rules_provider=atom_typing_rules_provider,
                            topology_graph=subgraph,
                        ),
                        "graph": subgraph,
                    }
                    typemap.update(reference[molecule]["typemap"])
                else:
                    if use_isomorphic_checks:
                        # Check for isomorphism submatching to typemap
                        matcher = nx.algorithms.isomorphism.GraphMatcher(
                            subgraph,
                            reference[molecule]["graph"],
                            node_match=top_node_match,
                        )
                        assert matcher.is_isomorphic()
                        for node in subgraph.nodes:
                            typemap[node] = reference[molecule]["typemap"][
                                matcher.mapping[node]
                            ]
                    else:
                        # Assume nodes in repeated structures are in the same order
                        for node, ref_node in zip(
                            sorted(subgraph.nodes),
                            sorted(reference[molecule]["typemap"]),
                        ):
                            typemap[node] = reference[molecule]["typemap"][
                                ref_node
                            ]
            return typemap
        elif use_isomorphic_checks:
            # Iterate through each isomorphic connected component
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

        else:
            return typemap_dict(
                topology_graph=foyer_topology_graph,
                atomtyping_rules_provider=atom_typing_rules_provider,
            )
