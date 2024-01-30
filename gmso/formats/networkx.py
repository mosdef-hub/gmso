"""NetworkX methods for operating with GMSO topologies."""

from gmso.utils import nx_utils
from gmso.utils.io import has_ipywidgets, import_, run_from_ipython

widgets = import_("ipywidgets")
plt = import_("matplotlib.pyplot")
if has_ipywidgets:
    from ipywidgets import interact, fixed

from gmso.external.convert_networkx import to_networkx


def interactive_networkx_atomtypes(topology, list_of_labels=None):
    """Get an interactive networkx plot showing the atom types of a topology object.

    Parameters
    ----------
    topology : gmso.Topology
        This should be a gmso topology object that you want to visualize the atom types
        that have been parameterized.
    list_of_labels : List[str]
         Default labels are ['atom_type.name','charge','mass','element','label','position'].
         Any additonal labels can be appended from the list of:
         `dir(list(topology.sites)[0])` or `dir(list(topology.sites)[0].atom_type)`

    Notes
    -----
    This function will output interactive ipywidgets. An ipywidgets dropdown object will
        select which atomtypes that should be shown on the networkx graph object below. This
        graph object is a visual representation for the connections within a gmso.Topology
        object.
    Select from a list of labels labels to be included on the graph.

    See Also
    --------
    gmso.formats.networkx.interactive_networkx_bonds, gmso.formats.networkx.interactive_networkx_angles,
    gmso.formats.networkx.interactive_networkx_dihedrals
    for other interactive visualization methods.

    ipywidgets.Dropdown()
        Two dropdown options, corresponding to:
            Label = What labels are shown on the networkx graph.
            Atom_Name = Which sites will show these labels.
    """
    if not run_from_ipython():
        raise RuntimeError(
            "Unsupported visualization outside of jupyter notebooks."
        )
    networkx_graph = to_networkx(topology)
    # get a unique list of site names
    site_names = []
    for node in networkx_graph.nodes:
        site_names.append(node.name)
    site_names = set(site_names)

    # create a tuple of selectable names for each site
    names_tuple = []
    names_tuple.append(("All Sites", None))
    for name in site_names:
        names_tuple.append((name, name))

    # Create list of labels to put on plot
    if not list_of_labels:
        list_of_labels = []
    base_list = [
        "atom_type.name",
        "charge",
        "mass",
        "element",
        "label",
        "position",
    ]
    list_of_labels += base_list

    @interact
    def get_interactive_sites(Label=list_of_labels, Atom_Name=names_tuple):
        """Return the interactive sites."""
        # Plot atom types for the widget inputs
        nx_utils.plot_networkx_atomtypes(topology, Atom_Name, [Label])

        return

    return


def interactive_networkx_bonds(topology, additional_labels=None):
    """Get an interactive networkx plot showing the bond types of a topology object.

    Parameters
    ----------
    topology : A gmso.core.topology.Topology object
        This should be a gmso topology object that you want to visualize the atom types
        that have been parameterized.
    additonal_labels : Labels at each site to be included on the plot.
        Default labels are ['atom_type.name'].
        Any additonal labels can be appended from the list of:
        `dir(list(topology.sites)[0])` or `dir(list(topology.sites)[0].atom_type)`

    Notes
    -----
    This function will output interactive ipywidgets. An ipywidgets dropdown object will
        select which atom names that will make up the bonds shown on the networkx graph object below.
        This graph object is a visual representation for the connections within a gmso.Topology
        object, and the bonds that makeup that bondtype are highlighted in red.
    Select from a list of available sites, the sites that make up bonds to be highlighted on the graph.

    See Also
    --------
    gmso.formats.networkx.interactive_networkx_atomtypes, gmso.formats.networkx.interactive_networkx_angles,
    gmso.formats.networkx.interactive_networkx_dihedrals
    for other interactive visualization methods.

    ipywidgets.Dropdown()
        Two dropdown options, corresponding to:
            Atom1 (req) = Filters through all bonds to find bonds that contain this site. If not specified,
                only bonds with missing bondtypes will be highlighted (in red).
            Atom2 (opt) = A second site to specify which bonds can be selected on the graph.
        A third dropdown option:
            list_of_bonds = a list of the bonds labeled by the types that makeup the two sites.

    matplotlib.pyplot.figure()
        A figure showing the networkx.graph object of the molecule with highlighted edges corresponding to the
            selected bonds

    ipywidgets.checkbox()
        A checkbox that will allow you to print the parameters associated with the selected bond in the list_of_bonds
            dropdown option.
    """
    if not run_from_ipython():
        raise RuntimeError(
            "Unsupported visualization outside of jupyter notebooks."
        )
    networkx_graph = to_networkx(topology)

    # Create a list of labels to go on the nodes
    if not additional_labels:
        additional_labels = []
    base_list = ["atom_type.name"]
    list_of_labels = base_list + additional_labels

    # Create list of nodes to plot
    site_names = []
    for node in networkx_graph.nodes:
        site_names.append(node.name)
    site_names = set(site_names)

    # Create a tuple of keys for each selected site
    names_tuple = []
    names_tuple.append(("All Sites", None))
    for name in site_names:
        names_tuple.append((name, name))
    atom_selection = []
    descriptions = ["Atom1 (req)", "Atom2 (opt)"]
    for i in [0, 1]:
        atom_selection.append(
            widgets.Dropdown(
                options=names_tuple,
                layout=widgets.Layout(width="30%"),
                style=dict(description_width="initial"),
                description=descriptions[i],
            )
        )
    interact(
        nx_utils.call_interactive_sites,
        Atom1=atom_selection[0],
        Atom2=atom_selection[1],
        list_of_labels=fixed(list_of_labels),
        networkx_graph=fixed(networkx_graph),
        topology=fixed(topology),
    )

    return


def interactive_networkx_angles(topology):
    """Get an interactive networkx plot showing the angle types of a topology object.

    Parameters
    ----------
    topology : A gmso.core.topology.Topology object
        This should be a gmso topology object that you want to visualize the angle types
        that have been parameterized.

    Notes
    -----
    This function will output interactive ipywidgets. An ipywidgets dropdown object will
        select which atom names that will make up the angles shown on the networkx graph object below.
        This graph object is a visual representation for the connections within a gmso.Topology
        object.
    Select from a list of available sites, the sites that make up angles to be highlighted (in red) on the graph.

    See Also
    --------
    gmso.formats.networkx.interactive_networkx_atomtypes, gmso.formats.networkx.interactive_networkx_bonds,
       gmso.formats.networkx.interactive_networkx_dihedrals
       for other interactive visualization methods.

    ipywidgets.Dropdown()
        Three dropdown options, corresponding to:
            Central Atom1 (req) = Filters through all angles to find angles that contain this site as the center of the angle.
                If not specified only angles with missing angletypes will be highlighted.
            Atom2 (opt) = A second site to filter which angles can be selected on the graph.
            Atom3 (opt) = A third site to filter which angles can be selected on the graph.
        A fourth dropdown option:
            Selected Edge = a list of the angles labeled by the types that makeup the three sites.

    matplotlib.pyplot.figure()
        A figure showing the networkx.graph object of the molecule with highlighted (in red) edges corresponding to the
            selected angles.

    ipywidgets.checkbox()
        A checkbox that will allow you to print the parameters associated with the selected angle in the Selected Edge
            dropdown option.
    """
    if not run_from_ipython():
        raise RuntimeError(
            "Unsupported visualization outside of jupyter notebooks."
        )

    networkx_graph = to_networkx(topology)

    # Create list of nodes to plot
    site_names = []
    for node in networkx_graph.nodes:
        site_names.append(node.name)
    site_names = set(site_names)

    # Create a tuple of keys for each selected site
    names_tuple = []
    names_tuple.append(("All Sites", None))
    for name in site_names:
        names_tuple.append((name, name))

    # Call recursive interacts. The top level determines what bonds can be selected
    atom_selection = []
    descriptions = ["Central Atom1 (req)", "Atom2 (opt)", "Atom3 (opt)"]
    for i in [0, 1, 2]:
        atom_selection.append(
            widgets.Dropdown(
                options=names_tuple,
                layout=widgets.Layout(width="30%"),
                style=dict(description_width="initial"),
                description=descriptions[i],
            )
        )
    interact(
        nx_utils.select_angles_from_sites,
        networkx_graph=fixed(networkx_graph),
        top=fixed(topology),
        Atom1=atom_selection[0],
        Atom2=atom_selection[1],
        Atom3=atom_selection[2],
    )

    return


def interactive_networkx_dihedrals(topology):
    """Get an interactive networkx plot showing the dihedral types of a topology object.

    Parameters
    ----------
    topology : A gmso.core.topology.Topology object
        This should be a gmso topology object that you want to visualize the dihedral types
        that have been parameterized.

    Notes
    -----
    This function will output interactive ipywidgets. An ipywidgets dropdown object will
        select which atom names that will make up the dihedrals shown on the networkx graph object below.
        This graph object is a visual representation for the connections within a gmso.Topology
        object.
    Select from a list of available sites, the sites that make up dihedrals to be highlighted (in red) on the graph.

    See Also
    --------
        gmso.formats.networkx.interactive_networkx_atomtypes, gmso.formats.networkx.interactive_networkx_bonds,
        gmso.formats.networkx.interactive_networkx_angles
        for other interactive visualization methods.

    ipywidgets.Dropdown()
        Four dropdown options, corresponding to:
            Central Atom1 (req) = Filters through all bonds to find dihedrals that contain this site in one of its central
                two positions. If not specified, only dihedrals with missing dihedraltypes will be highlighted.
            Central Atom2 (req) = Filters through all bonds to find dihedrals that contain this site in one of its central
                two positions. If not specified, only dihedrals with missing dihedraltypes will be highlighted.
            Atom3 (opt) = A third site to filter which dihedrals can be selected on the graph.
            Atom4 (opt) = A fourth site to filter which dihedrals can be selected on the graph.
        A fourth dropdown option:
            Selected Edge = a list of the angles labeled by the types that makeup the three sites.

    matplotlib.pyplot.figure()
        A figure showing the networkx.graph object of the molecule with highlighted (in red) edges corresponding to the
            selected dihedrals.

    ipywidgets.checkbox()
        A checkbox that will allow you to print the parameters associated with the selected dihedral in the Selected Edge
            dropdown option.
    """
    if not run_from_ipython():
        raise RuntimeError(
            "Unsupported visualization outside of jupyter notebooks."
        )

    networkx_graph = to_networkx(topology)
    # Create list of nodes to plot
    site_names = []
    for node in networkx_graph.nodes:
        site_names.append(node.name)
    site_names = set(site_names)

    # Create a tuple of keys for each selected site
    names_tuple = []
    names_tuple.append(("All Sites", None))
    for name in site_names:
        names_tuple.append((name, name))

    # Call recursive interacts. The top level determines what bonds can be selected
    atom_selection = []
    descriptions = [
        "Central Atom1 (req)",
        "Central Atom2 (req)",
        "Atom3 (opt)",
        "Atom4 (opt)",
    ]
    for i in [0, 1, 2, 3]:
        atom_selection.append(
            widgets.Dropdown(
                options=(names_tuple),
                layout=widgets.Layout(width="30%"),
                style=dict(description_width="initial"),
                description=descriptions[i],
            )
        )
    interact(
        nx_utils.select_dihedrals_from_sites,
        networkx_graph=fixed(networkx_graph),
        top=fixed(topology),
        Atom1=atom_selection[0],
        Atom2=atom_selection[1],
        Atom3=atom_selection[2],
        Atom4=atom_selection[3],
    )

    return
