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
    speedup_by_molgraph=False,
    speedup_by_moltag=False,
    ignore_params=["improper"],
    remove_untyped=True,
    fast_copy=True,
):
    """Set Topology parameter types from GMSO ForceFields.

    Parameters
    ----------
    top : gmso.core.topology.Topology, required
        The GMSO topology on which to apply forcefields

    forcefields : ForceField or dict, required
        The forcefield to apply. If a dictionary is used the keys are labels that match
        the molecule name (specified as a label of site), and the values are gmso ForceField objects that gets applied
        to the specified molecule.
        Note: if a Topology with no molecule is provided, this option will only take
        a ForceField object. If a dictionary of ForceFields is provided, this method will
        fail.

    match_ff_by : str, optional, default="molecule"
        They site's tag used to match the forcefields provided above to the Topology.
        Options include "molecule" and "group". This option is only valid if forcefields are provided
        as a dict.

    identify_connections : bool, optional, default=False
        If true, add connections identified using networkx graph matching to match
        the topology's bonding graph to smaller sub-graphs that correspond to an angle,
        dihedral, improper etc

    speedup_by_molgraph: bool, optional, default=False
        A flag to determine whether or not to search the topology for repeated disconnected
        structures, otherwise known as molecules and type each molecule only once.
        This option will be usefult to handle systems with many repeated small molecules,
        but may slow down system with large molecule, e.g., monolayer.

    speedup_by_moltag : bool, optional, default=False
        A flag to determine whether or not to look at site.molecule_name to try to parameterize
        each molecule only once. This option provides speedup for topologies with properly
        assigned molecule and residue labels.

    ignore_params : set or list or tuple, optional, default=["impropers"]
        Skipping the checks that make sure all connections (in the list) have a connection types.
        Available options includes "bonds", "angles", "dihedrals", and "impropers".
        If you wish to have all connection types checks, provides an empty set/list/tuple.

    remove_untyped : bool, optional, default=True
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
    ignore_params = set([option.lower() for option in ignore_params])
    config = TopologyParameterizationConfig.model_validate(
        dict(
            match_ff_by=match_ff_by,
            identify_connections=identify_connections,
            speedup_by_molgraph=speedup_by_molgraph,
            speedup_by_moltag=speedup_by_moltag,
            ignore_params=ignore_params,
            remove_untyped=remove_untyped,
            fast_copy=fast_copy,
        )
    )
    parameterizer = TopologyParameterizer(
        topology=top, forcefields=forcefields, config=config
    )

    parameterizer.run_parameterization()

    return parameterizer.topology
