from gmso.core.topology import Topology


def slice_topology_by_molecule(topology, molecule_tag, molecule_number=None):
    """Create a Topology that contains a subset of molecules from another Topology.

    Parameters
    ----------
    topology : gmso.core.topology.Topology
        The gmso Topology to perform the slice on
    molecule_tag : str
        The name of the gmso.abstract_site.Molecule object to include in the slice
    molecule_number : int, default None
        If given, only include a single molecule's sites
        If None, then all sites in every molecule matching `molecule_tag` are included
        in the sliced topology.

    Returns
    -------
    gmso.core.topology.Topology
        A new Topology instance containing only sites and connections from matching molecules.
    """
    sites = [
        s
        for s in topology.iter_sites_by_molecule(
            molecule_tag=molecule_tag, molecule_number=molecule_number
        )
    ]
    return slice_by_sites(topology=topology, sites=sites)


def slice_topology_by_residue(topology, residue_tag, residue_number=None):
    """Create a Topology that contains a subset of residues from another Topology.

    Parameters
    ----------
    topology : gmso.core.topology.Topology
        The gmso Topology to perform the slice on.
    residue_tag : str
        The name of the gmso.abstract_site.Residue object to include in the slice.
    residue_number : int, default None
        If given, only include a single residue's sites
        If None, then all sites in every residue matching `residue_tag` are included
        in the sliced topology.

    Returns
    -------
    gmso.core.topology.Topology
        A new Topology instance containing only sites and connections from matching molecules.
    """
    sites = [
        s
        for s in topology.iter_sites_by_residue(
            residue_tag=residue_tag, residue_number=residue_number
        )
    ]
    return slice_by_sites(topology=topology, sites=sites)


def slice_by_sites(topology, sites):
    """Used by slice_topology_by_molecule() and slice_topology_by_residue()

    Parameters
    ----------
    sites : list of gmso.core.atom.Atom
        List of sites to include in the sub-topology.
    topology : gmso.core.topology.Topology
        The topology being sliced.
    """
    new_topology = Topology()

    connections = set()
    for site in sites:
        new_topology.add_site(site)
        for connection in topology.iter_connections_by_site(site):
            connections.add(connection)

    for connection in connections:
        new_topology.add_connection(connection)

    return new_topology
