"""Functions used to atomtype a gmso.Topology."""

from typing import Dict, List, Set, Tuple, Union

from gmso import Topology, ForceField
from gmso.parameterization.topology_parameterizer import (
    TopologyParameterizationConfig,
    TopologyParameterizer,
)

__all__ = ["apply"]


def apply(
    top: Topology,
    forcefields: Union["gmso.ForceField", Dict[str, "gmso.ForceField"]],
    match_ff_by: str = "molecule",
    identify_connections: bool = False,
    speedup_by_molgraph: bool = False,
    speedup_by_moltag: bool = False,
    ignore_params: Union[List[str], Set[str], Tuple[str, ...]] = ["improper"],
    remove_untyped: bool = True,
    fast_copy: bool = True,
) -> Topology:
    """Apply forcefield parameters to a :class:`~gmso.Topology`.

    Atom-types all sites in *top* using the supplied forcefield(s) and
    writes the resulting potential types onto the topology.

    Parameters
    ----------
    top : gmso.Topology
        The un-typed topology to parameterize.
    forcefields : gmso.ForceField or dict
        The forcefield(s) to apply.  When a single :class:`~gmso.ForceField`
        is supplied it is applied to every site in *top*.  When a
        ``dict`` is supplied, each key must match a molecule name (or group
        name, see *match_ff_by*) present in the topology, and the
        corresponding :class:`~gmso.ForceField` is applied only to sites
        belonging to that molecule/group.  A topology with no molecule labels
        can only accept a single :class:`~gmso.ForceField`; passing a dict
        in that case will raise an error.

    match_ff_by : str, optional, default="molecule"
        Site attribute used to map forcefields when *forcefields* is a dict.
        Accepted values: ``"molecule"`` or ``"group"``.
    identify_connections : bool, optional, default=False
        When ``True``, use NetworkX graph matching to detect angles,
        dihedrals, and impropers from the bond graph before parameterizing.
    speedup_by_molgraph : bool, optional, default=False
        When ``True``, detect repeated disconnected sub-graphs (molecules) and
        parameterize each unique molecule only once.  Useful for systems with
        many identical small molecules; may slow down systems with large unique
        molecules (e.g., monolayers).
    speedup_by_moltag : bool, optional, default=False
        When ``True``, use ``site.molecule`` labels to identify repeated
        molecules and parameterize each unique molecule only once.  Requires
        that molecule and residue labels have been properly assigned.
    ignore_params : list, set, or tuple of str, optional, default=["improper"]
        Connection types whose completeness check is skipped after
        parameterization.  Valid entries: ``"bond"``, ``"angle"``,
        ``"dihedral"``, ``"improper"``.  Pass an empty collection to enforce
        checks for all connection types.
    remove_untyped : bool, optional, default=True
        When ``True``, remove all connections that have no associated
        connection type after parameterization.
    fast_copy : bool, optional, default=True
        When ``True``, sympy expressions and parameter dictionaries are
        *not* deep-copied when the same potential type is assigned to multiple
        sites or connections.  This is significantly faster but means that
        modifying a parameter on one site/connection will affect all others
        that share the same :class:`~gmso.AtomType` (or connection-type)
        object.  Set to ``False`` if you need to modify expressions or
        parameters after parameterization.

    Returns
    -------
    gmso.Topology
        The parameterized topology (modified in-place and returned).

    Examples
    --------
    Apply a single forcefield to an entire topology:

    >>> from gmso import ForceField
    >>> from gmso.parameterization import apply
    >>> ff = ForceField("oplsaa.xml")
    >>> typed_top = apply(top, ff)

    Apply different forcefields to different molecule types:

    >>> ff_water = ForceField("spce.xml")
    >>> ff_ethanol = ForceField("oplsaa.xml")
    >>> typed_top = apply(top, {"water": ff_water, "ethanol": ff_ethanol})
    """
    ignore_params = set([option.lower().rstrip("s") for option in ignore_params])
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
