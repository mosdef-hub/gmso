"""Functions used to atomtype a gmso.Topology."""
from gmso.parameterization.topology_parameterizer import (
    TopologyParameterizationConfig,
    TopologyParameterizer,
)

__all__ = ["apply"]


def apply(
    top,
    forcefields,
    identify_connections=False,
    identify_connected_components=True,
    use_residue_info=False,
    assert_bond_params=True,
    assert_angle_params=True,
    assert_dihedral_params=True,
    assert_improper_params=False,
):
    """Set Topology parameter types from GMSO Forcefields.

    Parameters
    ----------
    top: gmso.core.topology.Topology, required
        The GMSO topology on which to apply forcefields
    forcefields: dict, required
        The keys are labels that match the subtopology label or site residue_name, and the
        values are gmso Forcefield objects that gets applied to the specified molecule

    identify_connections: bool, optional, default=False
        If true, add connections identified using networkx graph matching to match
        the topology's bonding graph to smaller sub-graphs that correspond to an angle,
        dihedral, improper etc

    identify_connected_components: bool, optional, default=True
        A flag to determine whether or not to search the topology for repeated disconnected
        structures, otherwise known as molecules and type each molecule only once.
    use_residue_info: bool, optional, default=False
        A flag to determine whether or not to look at site.residue_name to look parameterize
        each molecule only once. Will only be used if identify_connected_components=False
    assert_bond_params : bool, optional, default=True
        If True, an error is raised if parameters are not found for all system
        bonds.
    assert_angle_params : bool, optional, default=True
        If True, an error is raised if parameters are not found for all system
        angles.
    assert_dihedral_params : bool, optional, default=True
        If True, an error is raised if parameters are not found for all system
        proper dihedrals.
    assert_improper_params : bool, optional, default=False
        If True, an error is raised if parameters are not found for all system
        improper dihedrals.
    """
    config = TopologyParameterizationConfig.parse_obj(
        dict(
            identify_connections=identify_connections,
            identify_connected_components=identify_connected_components,
            use_residue_info=use_residue_info,
            assert_bond_params=assert_bond_params,
            assert_angle_params=assert_angle_params,
            assert_dihedral_params=assert_dihedral_params,
            assert_improper_params=assert_improper_params,
        )
    )
    parameterizer = TopologyParameterizer(
        topology=top, forcefields=forcefields, config=config
    )

    parameterizer.run_parameterization()

    return parameterizer.topology
