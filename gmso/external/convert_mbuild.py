"""Convert to and from an mbuild.Compound."""

from warnings import warn

import mbuild as mb
import numpy as np
import unyt as u
from boltons.setutils import IndexedSet
from unyt import Unit

from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.box import Box
from gmso.core.element import (
    element_by_atomic_number,
    element_by_mass,
    element_by_name,
    element_by_symbol,
)
from gmso.core.topology import Topology
from gmso.exceptions import GMSOError
from gmso.utils.io import has_mbuild

if has_mbuild:
    import mbuild as mb


def from_mbuild(
    compound,
    box=None,
    search_method=element_by_symbol,
    parse_label=True,
    custom_groups=None,
    infer_elements=False,
):
    """Convert an mbuild.Compound to a gmso.Topology.

    This conversion makes the following assumptions about the inputted `Compound`:

        * All positional and box dimension values in compound are in nanometers.

        * The hierarchical structure of the Compound will be flattened and translated to labels\
         in GMSO Sites. The directly supported labels include `Site.group`,\
        `Site.molecule_name`, and `Site.residue_name`.

            * `group` is determined as te second-highest level Compound and is automatically generated;\
            * `molecule` is determined by traversing through\
            hierarchy of the mb.Compound, starting from the particle level, until the lowest\
            independent mb.Compound is reached (determined as an mb.Compound that does not have\
            any bond outside its boundary);
            * `residue` is the `mb.Compound` level right above particle level.
            * `molecule` and `residue` take the format of (name, index), where the latter can be used\
            to distinguish between molecule/residue of the same name. These two labels are only generated\
            if parse_label=True.

        * Only `Bonds` are added for each bond in the `Compound`. If `Angles`\
          and `Dihedrals` are desired in the resulting `Topology`, they must be\
          added separately from this function.

    Parameters
    ----------
    compound : mbuild.Compound
        mbuild.Compound instance that need to be converted
    box : mbuild.Box, optional, default=None
        Box information to be loaded to a gmso.Topology
    search_method : function, optional, default=element_by_symbol
        Searching method used to assign element from periodic table to particle site.
        The information specified in the `search_method` argument is extracted
        from each `Particle`'s `name` attribute.
        Valid functions are element_by_symbol, element_by_name,
        element_by_atomic_number, and element_by_mass, which can be imported
        from `gmso.core.element`
    parse_label : bool, optional, default=True
        Option to parse hierarchy info of the compound into system of top label,
        including, group, molecule and residue labels.
    custom_groups : list or str, optional, default=None
        Allows user to identify the groups assigned to each site in the topology
        based on the compound.name attributes found traversing down the hierarchy. Be
        sure to supply names such that every particle will be pass through one
        matching name on the way down from compound.children. Only the first match
        while moving downwards will be assigned to the site. If parse_label=False,
        this argument does nothing.
    infer_elements : bool, default=False
        Allows the reader to try to load element info from the mbuild Particle.name
        instead of only from the populated Particle.element

    Returns
    -------
    top : gmso.Topology
    """
    msg = "Argument compound is not an mbuild.Compound"
    assert isinstance(compound, mb.Compound), msg
    msg = "Compound is not a top level compound. Make a copy to pass to the `compound` \
    argument that has no parents"
    assert not compound.parent, msg

    top = Topology()
    top.typed = False

    site_map = {
        particle: {
            "site": None,
            "residue": None,
            "molecule": None,
            "group": None,
        }
        for particle in compound.particles()
    }
    if parse_label:
        _parse_molecule_residue(site_map, compound)
        _parse_group(site_map, compound, custom_groups)

    # Use site map to apply Compound info to Topology.
    for part in compound.particles():
        site = _parse_site(
            site_map, part, search_method, infer_element=infer_elements
        )
        top.add_site(site)

    for b1, b2 in compound.bonds():
        assert site_map[b1]["site"].molecule == site_map[b2]["site"].molecule
        new_bond = Bond(
            connection_members=[site_map[b1]["site"], site_map[b2]["site"]],
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


def to_mbuild(topology, infer_hierarchy=True):
    """Convert a gmso.Topology to mbuild.Compound.

    Parameters
    ----------
    topology : gmso.Topology
        topology instance that need to be converted
    infer_hierarchy : bool, optional, default=True
        Option to infer the hierarchy from Topology's labels
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
    if not infer_hierarchy:
        particle_list = []
        for site in topology.sites:
            particle = _parse_particle(particle_map=particle_map, site=site)
            particle_list.append(particle)
        compound.add(particle_list)

    else:
        molecule_list = []
        for molecule_tag in topology.unique_site_labels(label_type="molecule"):
            mb_molecule = mb.Compound()
            mb_molecule.name = (
                molecule_tag.name if molecule_tag else "DefaultMolecule"
            )
            residue_dict = dict()
            residue_dict_particles = dict()

            if molecule_tag:
                sites_iter = topology.iter_sites("molecule", molecule_tag)
            else:
                sites_iter = (
                    site for site in topology.sites if not site.molecule
                )

            for site in sites_iter:
                particle = _parse_particle(particle_map, site)
                # Try to add the particle to a residue level
                residue_tag = (
                    site.residue if site.residue else ("DefaultResidue", 0)
                )  # the 0 idx is placeholder and does nothing
                if residue_tag in residue_dict:
                    residue_dict_particles[residue_tag] += [particle]
                else:
                    residue_dict[residue_tag] = mb.Compound(name=residue_tag[0])
                    residue_dict_particles[residue_tag] = [particle]

            for key, item in residue_dict.items():
                item.add(residue_dict_particles[key])
                molecule_list.append(item)

        compound.add(molecule_list)
    for connect in topology.bonds:
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


def _parse_particle(particle_map, site):
    """Parse information for a mb.Particle from a gmso.Site and add it to particle map."""
    element = site.element.symbol if site.element else None
    charge = (
        site.charge.in_units(u.elementary_charge).value if site.charge else None
    )
    mass = site.mass.in_units(u.amu).value if site.mass else None

    particle = mb.Compound(
        name=site.name,
        pos=site.position.to_value(u.nm),
        element=element,
        charge=charge,
        mass=mass,
    )
    particle_map[site] = particle
    return particle


def _parse_site(site_map, particle, search_method, infer_element=False):
    """Parse information for a gmso.Site from a mBuild.Compound and add it to the site map."""
    pos = particle.xyz[0] * u.nm
    if particle.element:
        ele = search_method(particle.element.symbol)
    else:
        ele = search_method(particle.name) if infer_element else None

    charge = particle.charge * u.elementary_charge if particle.charge else None
    mass = particle.mass * Unit("amu") if particle.mass else None

    site = Atom(
        name=particle.name,
        position=pos,
        element=ele,
        charge=charge,
        mass=mass,
        molecule=site_map[particle]["molecule"],
        residue=site_map[particle]["residue"],
        group=site_map[particle]["group"],
    )
    site_map[particle]["site"] = site
    return site


def _parse_molecule_residue(site_map, compound):
    """Parse information necessary for residue and molecule labels when converting from mbuild."""
    connected_subgraph = compound.bond_graph.connected_components()
    molecule_tracker = dict()
    residue_tracker = dict()
    for molecule in connected_subgraph:
        if len(molecule) == 1:
            ancestors = [molecule[0]]
        else:
            ancestors = IndexedSet(molecule[0].ancestors())
            for particle in molecule[1:]:
                # This works because particle.ancestors traversed, and hence
                # the lower level will be in the front.
                # The intersection will be left at the end,
                # ancestor of the first particle is used as reference.
                # Hence, this called will return the lowest-level Compound]
                # that is a molecule
                ancestors = ancestors.intersection(
                    IndexedSet(particle.ancestors())
                )

        """Parse molecule information"""
        molecule_tag = ancestors[0]
        # The duplication is determined solely by name
        if molecule_tag.name in molecule_tracker:
            molecule_tracker[molecule_tag.name] += 1
        else:
            molecule_tracker[molecule_tag.name] = 0
        molecule_number = molecule_tracker[molecule_tag.name]
        """End of molecule parsing"""

        for particle in molecule:
            """Parse residue information"""
            residue_tag = (
                particle if not particle.n_direct_bonds else particle.parent
            )
            if residue_tag.name in residue_tracker:
                if residue_tag not in residue_tracker[residue_tag.name]:
                    residue_tracker[residue_tag.name][residue_tag] = len(
                        residue_tracker[residue_tag.name]
                    )
            else:
                residue_tracker[residue_tag.name] = {residue_tag: 0}

            residue_number = residue_tracker[residue_tag.name][residue_tag]
            site_map[particle]["residue"] = (residue_tag.name, residue_number)
            site_map[particle]["molecule"] = (
                molecule_tag.name,
                molecule_number,
            )


def _parse_group(site_map, compound, custom_groups):
    """Parse group information."""
    if custom_groups:
        if isinstance(custom_groups, str):
            custom_groups = [custom_groups]
        elif not hasattr(custom_groups, "__iter__"):
            raise TypeError(
                f"Please pass groups {custom_groups} as a list of strings."
            )
        elif not np.all([isinstance(g, str) for g in custom_groups]):
            raise TypeError(
                f"Please pass groups {custom_groups} as a list of strings."
            )
        for part in _traverse_down_hierarchy(compound, custom_groups):
            for particle in part.particles():
                site_map[particle]["group"] = part.name
        try:
            applied_groups = set(map(lambda x: x["group"], site_map.values()))
            assert applied_groups == set(custom_groups)
        except AssertionError:
            warn(
                f"""Not all custom groups ({custom_groups}, is are being used when
            traversing compound hierachy. Only {applied_groups} are used.)"""
            )
    elif not compound.children:
        for particle in compound.particles():
            site_map[particle]["group"] = compound.name
    elif not np.any(
        list(map(lambda c: len(c.children), compound.children))
    ):  # compound is a 2 level hierarchy
        for particle in compound.particles():
            site_map[particle]["group"] = compound.name
    else:  # set compund name to se
        for child in compound.children:
            for particle in child.particles():
                site_map[particle]["group"] = child.name


def _traverse_down_hierarchy(compound, group_names):
    if compound.name in group_names:
        yield compound
    elif compound.children:
        for child in compound.children:
            yield from _traverse_down_hierarchy(child, group_names)
    else:
        raise GMSOError(
            f"""A particle named {compound.name} cannot be associated with the
        custom_groups {group_names}. Be sure to specify a list of group names that will cover
        all particles in the compound. This particle is one level below {compound.parent.name}."""
        )
