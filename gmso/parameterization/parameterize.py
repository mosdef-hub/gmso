"""Functions used to atomtype a gmso.Topology."""
from gmso.parameterization.topology_parameterizer import (
    TopologyParameterizationConfig,
    TopologyParameterizer,
)

__all__ = ["apply"]


def apply(
    top,
    forcefields,
    match_ff_by="molecule",
    identify_connections=False,
    identify_connected_components=True,
    use_molecule_info=False,
    assert_bond_params=True,
    assert_angle_params=True,
    assert_dihedral_params=True,
    assert_improper_params=False,
    remove_untyped=False,
    fast_copy=True,
):
    """Set Topology parameter types from GMSO ForceFields.

    Parameters
    ----------
    top: gmso.core.topology.Topology, required
        The GMSO topology on which to apply forcefields

    forcefields: ForceField or dict, required
        The forcefield to apply. If a dictionary is used the keys are labels that match
        the molecule name (specified as a label of site), and the values are gmso ForceField objects that gets applied
        to the specified molecule.
        Note: if a Topology with no molecule is provided, this option will only take
        a ForceField object. If a dictionary of ForceFields is provided, this method will
        fail.

    match_ff_by: str, optional, default="molecule"
        They site's tag used to match the forcefields provided above to the Topology.
        Options include "molecule" and "group". This option is only valid if forcefields are provided
        as a dict.

    identify_connections: bool, optional, default=False
        If true, add connections identified using networkx graph matching to match
        the topology's bonding graph to smaller sub-graphs that correspond to an angle,
        dihedral, improper etc

    identify_connected_components: bool, optional, default=True
        A flag to determine whether or not to search the topology for repeated disconnected
        structures, otherwise known as molecules and type each molecule only once.

    use_molecule_info: bool, optional, default=False
        A flag to determine whether or not to look at site.residue_name to look parameterize
        each molecule only once. Currently unused.

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

    remove_untyped : bool, optional, default=False
        If True, after the atomtyping and parameterization step, remove all connection
        that has no connection_type.

    fast_copy : bool, optional, default=True
        If True, sympy expressions and parameters will not be deep copied during replicated
        parameterization. This can lead to the potentials for multiple sites/connections
        to be changed if a single parameter_type independent variable or expression is
        modified after the topology is parameterized. However, this leads to much faster
        application of forcefield parameters, and so is defaulted to True. Note that
        this should be changed to False if further modification of expressions are
        necessary post parameterization.
    """
    config = TopologyParameterizationConfig.parse_obj(
        dict(
            match_ff_by=match_ff_by,
            identify_connections=identify_connections,
            identify_connected_components=identify_connected_components,
            use_molecule_info=use_molecule_info,
            assert_bond_params=assert_bond_params,
            assert_angle_params=assert_angle_params,
            assert_dihedral_params=assert_dihedral_params,
            assert_improper_params=assert_improper_params,
            remove_untyped=remove_untyped,
            fast_copy=True,
        )
    )
    parameterizer = TopologyParameterizer(
        topology=top, forcefields=forcefields, config=config
    )

    parameterizer.run_parameterization()

    return parameterizer.topology
