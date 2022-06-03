"""Convert to and from an mbuild.Compound."""
from warnings import warn

import numpy as np
import unyt as u

from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.box import Box
from gmso.core.element import (
    element_by_atomic_number,
    element_by_mass,
    element_by_name,
    element_by_symbol,
)
from gmso.core.subtopology import SubTopology
from gmso.core.topology import Topology
from gmso.utils.io import has_mbuild

if has_mbuild:
    import mbuild as mb


def from_mbuild(
    compound, box=None, search_method=element_by_symbol, parse_label=False
):
    """Convert an mbuild.Compound to a gmso.Topology.

    This conversion makes the following assumptions about the inputted `Compound`:

        * All positional and box dimension values in compound are in nanometers.

        * The hierarchical structure of the Compound will be flattened a translated to a
        system of label in GMSO Sites. The directly supported labels include `Site.group`,
        `Site.molecule_name`, and `Site.residue_name`. `group` is determined as the
        second-highest level Compound; `molecule_name` is determined by traversing through
        hierarchy of the mb.Compound, starting from the particle level, until the lowest
        independent mb.Compound is reached (determined as an mb.Compound that does not have
        any bond outside its boundary); `residue_name` is the `mb.Compound` level right above
        particle level. `molecule_number` and `residue_number` can be used to distinguish between
        different molecules and residues of the same name. There may be overlap between these labels.

        * Only `Bonds` are added for each bond in the `Compound`. If `Angles`\
          and `Dihedrals` are desired in the resulting `Topology`, they must be\
          added separately from this function.

    Parameters
    ----------
    compound : mbuild.Compound
        mbuild.Compound instance that need to be converted
    box : mbuild.Box, optional, default=None
        Box information to be loaded to a gmso.Topologyl
    search_method : function, optional, default=element_by_symbol
        Searching method used to assign element from periodic table to
        particle site.
        The information specified in the `search_method` argument is extracted
        from each `Particle`'s `name` attribute.
        Valid functions are element_by_symbol, element_by_name,
        element_by_atomic_number, and element_by_mass, which can be imported
        from `gmso.core.element'
    parse_label : bool, optional, default=False
        Option to parse hierarchy info of the compound into system of top label,
        including, group, molecule and residue labels.

    Returns
    -------
    top : gmso.Topology
    """
    msg = "Argument compound is not an mbuild.Compound"
    assert isinstance(compound, mb.Compound), msg

    top = Topology()
    top.typed = False

    site_map = dict()
    if not parse_label:

        for particle in compound.particles():
            pos = particle.xyz[0] * u.nanometer
            if particle.element:
                ele = search_method(particle.element.symbol)
            else:
                ele = search_method(particle.name)
            site = Atom(
                name=particle.name,
                position=pos,
                element=ele,
            )
            site_map[particle] = site
            top.add_site(site)
    else:
        connected_subgraph = compound.bond_graph.connected_components()
        molecule_tracker = dict()
        for molecule in connected_subgraph:
            ancestors = IndexedSet(molecule[0].ancestors())
            for particle in molecule:
                ancestors.intersection(IndexedSet(particle.ancestors()))

            # This works because particle.ancestors traversed, and hence
            # the lower level will be in the front.
            # The intersection will only remove, not add to the IndexedSet
            # of the first particle we used as reference.
            # Hence, this called will return the lowest-level Compound
            # that is a molecule
            """Parse molecule information"""
            molecule_tag = ancestors[0]
            if molecule_tag in molecule_tracker:
                molecule_tracker[molecule_tag] += 1
            else:
                molecule_tracker[molecule_tag] = 0
            molecule_number = molecule_tracker[molecule_tag]
            """End of molecule parsing"""

            residue_tracker = dict()
            for particle in molecule:
                pos = particle.xyz[0] * u.nanometer
                if particle.element:
                    ele = search_method(particle.element.symbol)
                else:
                    ele = search_method(particle.name)

                """Parse residue information"""
                residue_tag = particle.parent
                if residue_tag in residue_tracker:
                    residue_tracker[residue_tag] += 1
                else:
                    residue_tracker[residue_tag] = 0
                residue_number = residue_tracker[residue_tag]

                site = Atom(
                    name=particle.name,
                    position=pos,
                    element=ele,
                    molecule=(molecule_tag.name, molecule_number),
                    residue=(residue_tag.name, residue_number),
                )
                site_map[particle] = site
                top.add_site(site)

    for b1, b2 in compound.bonds():
        assert site_map[b1].molecule == site_map[b2].molecule
        new_bond = Bond(
            connection_members=[site_map[b1], site_map[b2]],
            bond_type=None,
        )
        top.add_connection(new_bond, update_types=False)

    if box:
        top.box = from_mbuild_box(box)
    # Assumes 2-D systems are not supported in mBuild
    # if compound.periodicity is None and not box:
    else:
        if compound.box:
            top.box = from_mbuild_box(compound.box)
        else:
            top.box = from_mbuild_box(compound.get_boundingbox())
    top.periodicity = compound.periodicity

    return top


def to_mbuild(topology):
    """Convert a gmso.Topology to mbuild.Compound.

    Parameters
    ----------
    topology : gmso.Topology
        topology instance that need to be converted

    Returns
    -------
    compound : mbuild.Compound

    """
    msg = "Argument topology is not a Topology"
    assert isinstance(topology, Topology), msg

    compound = mb.Compound()
    if topology.name is Topology().name:
        compound.name = "Compound"
    else:
        compound.name = topology.name

    particle_map = dict()
    for site in topology.sites:
        if site.element:
            element = site.element.symbol
        else:
            element = None
        particle = mb.Compound(
            name=site.name, pos=site.position, element=element
        )
        particle_map[site] = particle
        compound.add(particle)

    for connect in topology.connections:
        if isinstance(connect, Bond):
            compound.add_bond(
                (
                    particle_map[connect.connection_members[0]],
                    particle_map[connect.connection_members[1]],
                )
            )

    return compound


def from_mbuild_box(mb_box):
    """Convert an mBuild box to a GMSO box.

    Assumes that the mBuild box dimensions are in nanometers

    Parameters
    ----------
    mb_box : mbuild.Box
        mBuild box object to be converted to a gmso.core.Box object

    Returns
    -------
    box : gmso.core.Box

    """
    # TODO: Unit tests

    if not isinstance(mb_box, mb.Box):
        raise ValueError("Argument mb_box is not an mBuild Box")

    if np.allclose(mb_box.lengths, [0, 0, 0]):
        warn("No box or boundingbox information detected, setting box to None")
        return None

    box = Box(
        lengths=np.asarray(mb_box.lengths) * u.nm,
        angles=np.asarray(mb_box.angles) * u.degree,
    )

    return box
