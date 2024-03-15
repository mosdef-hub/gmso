"""NetowrkX support for GMSO."""

import networkx as nx
import unyt

from gmso.utils.io import has_ipywidgets, import_

widgets = import_("ipywidgets")
plt = import_("matplotlib.pyplot")
if has_ipywidgets:
    from ipywidgets import interact, fixed

from gmso.external.convert_networkx import to_networkx


def plot_networkx_params(networkx_graph, list_of_edges):
    """Return a networkx plot showing the specified edges in a networkx graph object.

    Parameters
    ----------
    networkx_graph : A networkx.Graph object
        This is a networkx graph object that can be created from a topology using the to_networkx
        function found in gmso.external.convert_networkx
    list_of_edges : a list of edges that should be shown on the plot
        Will be of the shape [(node1,node2), (node2,node3), etc...]

    Returns
    -------
    matplotlib.pyplot.figure
        The drawn networkx plot of the topology on a matplotlib editable figures

    matplotlib.pyplot.axis
        The axis information that corresponds to a networkx drawn figure. This output can be
        shown using
        matplotlib.pyplot.show()
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    edge_weights, edge_colors = highlight_networkx_edges(
        networkx_graph, list_of_edges
    )
    ax = plot_networkx_nodes(
        networkx_graph, ax, edge_weights=edge_weights, edge_colors=edge_colors
    )

    return (fig, ax)


def call_interactive_sites(
    Atom1, Atom2, networkx_graph, topology, list_of_labels
):
    """Return a static plot or interactive display for visualization."""
    list_of_edges = get_edges(networkx_graph, Atom1, Atom2)
    if list_of_edges:
        # Call the second level, which takes the selected bonds and shows them on the figure
        interact(
            select_edges,
            list_of_bonds=list_of_edges,
            list_of_labels=fixed(list_of_labels),
            networkx_graph=fixed(networkx_graph),
            topology=fixed(topology),
        )
    else:
        plot_networkx_bonds(networkx_graph, list_of_labels=list_of_labels)

    return


def select_edges(networkx_graph, topology, list_of_bonds, list_of_labels):
    """Select and view edges in the topology."""
    plot_networkx_bonds(
        networkx_graph,
        list_of_labels=list_of_labels,
        list_of_bonds=list_of_bonds,
    )
    # The final level prints the parameters of the selected bond using a checkbox.
    checkbox = widgets.Checkbox(value=False, description="Show parameters")
    interact(
        show_bond_info,
        w=checkbox,
        topology=fixed(topology),
        list_of_bonds=fixed(list_of_bonds),
    )

    return


def show_bond_info(w, topology, list_of_bonds):
    """Return detailed information about the bonds."""
    if w:
        report_bond_parameters(topology, list_of_bonds)
    else:
        # TODO: Should be able to remove this blank print statement so that deselecting the
        # checkbox removes the listed parameters.
        print("")

    return


def select_angles_from_sites(networkx_graph, top, Atom1, Atom2, Atom3):
    """Return angles based on interactive selection."""
    params_list = select_params_on_networkx(
        networkx_graph, [Atom1, Atom2, Atom3]
    )
    if params_list:
        edges_widget = widgets.Dropdown(
            options=params_list,
            layout=widgets.Layout(width="60%"),
            style=dict(description_width="initial"),
            description="Selected Edge",
        )
        interact(
            select_edges_on_networkx,
            networkx_graph=fixed(networkx_graph),
            top=fixed(top),
            list_of_params=edges_widget,
        )
    else:
        plot_networkx_params(networkx_graph, list_of_edges=[])

    return


def select_dihedrals_from_sites(
    networkx_graph, top, Atom1=None, Atom2=None, Atom3=None, Atom4=None
):
    """Select dihedrals based on interactive selection."""
    params_list = select_params_on_networkx(
        networkx_graph, [Atom1, Atom2, Atom3, Atom4]
    )
    if params_list:
        edges_widget = widgets.Dropdown(
            options=params_list,
            layout=widgets.Layout(width="60%"),
            style=dict(description_width="initial"),
            description="Selected Edge",
        )
        interact(
            select_edges_on_networkx,
            networkx_graph=fixed(networkx_graph),
            top=fixed(top),
            list_of_params=edges_widget,
        )
    else:
        plot_networkx_params(networkx_graph, list_of_edges=[])

    return


def select_params_on_networkx(networkx_graph, atoms):
    """Create a list of the edges that have been selected."""
    # Create a list of the edges that have been selected corresponding to the selected
    # list of sites names. If all sites are None, then select edges that are missing parameter types.
    # Works for selecting angles or dihedals based on the number of atoms sent in the `atoms` variable.
    list_of_params = []
    list_of_names = []
    mia_angle_flag = 0
    if len(atoms) == 3:
        if all(atoms):
            for node, angles in networkx_graph.nodes(data="angles"):
                for angle in angles:
                    names = [member.name for member in angle.connection_members]
                    if atoms == [names[i] for i in [1, 0, 2]] or atoms == [
                        names[i] for i in [1, 2, 0]
                    ]:
                        list_of_params.append(angle.connection_members)
                        list_of_names.append(
                            _get_formatted_atom_types_names_for(angle)
                        )
        elif all(atoms[0:2]):
            for node, angles in networkx_graph.nodes(data="angles"):
                for angle in angles:
                    names = [member.name for member in angle.connection_members]
                    if atoms[0:2] == names[1:3] or atoms[0:2] == list(
                        reversed(names[0:2])
                    ):
                        list_of_params.append(angle.connection_members)
                        list_of_names.append(
                            _get_formatted_atom_types_names_for(angle)
                        )
        elif atoms[0]:
            for node, angles in networkx_graph.nodes(data="angles"):
                for angle in angles:
                    names = [member.name for member in angle.connection_members]
                    if atoms[0] == names[1]:
                        list_of_params.append(angle.connection_members)
                        list_of_names.append(
                            _get_formatted_atom_types_names_for(angle)
                        )
        else:
            for node, angles in networkx_graph.nodes(data="angles"):
                if not angles:
                    print(
                        "No angle parameters have been applied to this topology"
                    )
                    return
                for angle in angles:
                    if angle.angle_type is None:
                        list_of_params.append(angle.connection_members)
                        list_of_names.append(
                            _get_formatted_atom_types_names_for(angle)
                        )
                        mia_angle_flag = 1
            if not mia_angle_flag:
                print(
                    "All angles are typed. Select a central atom to look at the different angle_types."
                )
            else:
                print(
                    "Since no sites are input, angles with missing types are shown."
                )

    elif len(atoms) == 4:
        # select a list of dihedrals
        if all(atoms):
            for node, dihedrals in networkx_graph.nodes(data="dihedrals"):
                for dihedral in dihedrals:
                    names = [
                        member.name for member in dihedral.connection_members
                    ]
                    if (
                        atoms == [names[i] for i in [2, 1, 0, 3]]
                        or atoms == [names[i] for i in [1, 2, 3, 0]]
                        or atoms == [names[i] for i in [1, 2, 0, 3]]
                        or atoms == [names[i] for i in [2, 1, 3, 0]]
                    ):
                        list_of_params.append(dihedral.connection_members)
                        list_of_names.append(
                            _get_formatted_atom_types_names_for(dihedral)
                        )
        elif all(atoms[0:3]):
            for node, dihedrals in networkx_graph.nodes(data="dihedrals"):
                for dihedral in dihedrals:
                    names = [
                        member.name for member in dihedral.connection_members
                    ]
                    if (
                        atoms[0:3] == [names[i] for i in [1, 2, 3]]
                        or atoms[0:3] == [names[i] for i in [2, 1, 3]]
                        or atoms[0:3] == [names[i] for i in [1, 2, 0]]
                        or atoms[0:3] == [names[i] for i in [2, 1, 0]]
                    ):
                        list_of_params.append(dihedral.connection_members)
                        list_of_names.append(
                            _get_formatted_atom_types_names_for(dihedral)
                        )
        elif all(atoms[0:2]):
            for node, dihedrals in networkx_graph.nodes(data="dihedrals"):
                for dihedral in dihedrals:
                    names = [
                        member.name for member in dihedral.connection_members
                    ]
                    if atoms[0:2] == names[1:3] or atoms[0:2] == [
                        names[2],
                        names[1],
                    ]:
                        list_of_params.append(dihedral.connection_members)
                        list_of_names.append(
                            _get_formatted_atom_types_names_for(dihedral)
                        )
        else:
            for node, dihedrals in networkx_graph.nodes(data="dihedrals"):
                if not dihedrals:
                    print(
                        "No dihedral parameters have been applied to this topology"
                    )
                    return
                for dihedral in dihedrals:
                    if dihedral.dihedral_type is None:
                        list_of_params.append(dihedral.connection_members)
                        list_of_names.append(
                            _get_formatted_atom_types_names_for(dihedral)
                        )
                        mia_angle_flag = 1

            if not mia_angle_flag:
                print(
                    "All dihedrals are typed. Select two central atoms to see associated dihedrals."
                )
            else:
                print(
                    "Since no sites are input, dihedrals with missing types are shown."
                )

    else:
        print("invalid atom selections")

    labeled_params = zip(list_of_names, list_of_params)

    # create a dict so each selected bond selects the proper edges on the networkx.graph object.
    selectable_list = {}
    for label, param in labeled_params:
        if label in selectable_list.keys():
            selectable_list[label].append(param)
        else:
            selectable_list[label] = []
            selectable_list[label].append(param)

    # turn the dict selectable list into a list of tuples.
    list_of_edges = []
    for key in selectable_list:
        list_of_edges.append((key, selectable_list[key]))

    return list_of_edges


def _get_formatted_atom_types_names_for(connection):
    """Return formatted atom_type names for a connection."""
    names = []
    for member in connection.connection_members:
        if not member.atom_type:
            label = ""
        else:
            label = member.atom_type.name
        names.append(label)

    return " --- ".join(names)


def select_edges_on_networkx(networkx_graph, top, list_of_params):
    """Plot or display the selected connection parameters."""
    list_of_edges = get_networkx_edges(list_of_params)
    plot_networkx_params(networkx_graph, list_of_edges=list_of_edges)

    # The interact prints the parameters of the selected bond using a checkbox.
    checkbox = widgets.Checkbox(value=False, description="Show parameters")
    interact(
        show_parameter_values,
        topology=fixed(top),
        list_of_params=fixed(list_of_params),
        checkbox=checkbox,
    )

    return


def get_networkx_edges(list_of_params):
    """Return a list of edges within a given dihedral or angle."""
    # Return a list of edges within a given dihedral or angle
    # Both orientations of every edge are saved, to guarentee the edge can be
    # found on a networkx_graph.edge object.
    list_of_edges = []
    if len(list_of_params[0]) == 4:
        for param in list_of_params:
            edge1 = (param[0], param[1])
            edge2 = (param[1], param[2])
            edge3 = (param[1], param[0])
            edge4 = (param[2], param[1])
            edge5 = (param[2], param[3])
            edge6 = (param[3], param[2])
            list_of_edges.append(edge1)
            list_of_edges.append(edge2)
            list_of_edges.append(edge3)
            list_of_edges.append(edge4)
            list_of_edges.append(edge5)
            list_of_edges.append(edge6)
    elif len(list_of_params[0]) == 3:
        for param in list_of_params:
            edge1 = (param[0], param[1])
            edge2 = (param[1], param[2])
            edge3 = (param[1], param[0])
            edge4 = (param[2], param[1])
            list_of_edges.append(edge1)
            list_of_edges.append(edge2)
            list_of_edges.append(edge3)
            list_of_edges.append(edge4)
    else:
        raise ValueError(
            "The parameters are not proper angles or dihedrals. Connection members are missing."
        )

    return list_of_edges


def plot_networkx_nodes(
    networkx_graph,
    ax,
    atom_name=None,
    edge_weights=None,
    edge_colors=None,
    node_sizes=None,
    list_of_labels=["atom_type.name"],
):
    """Plot the nodes of the networkX graph."""
    # Place nodes at 2D positions related to position in the topology
    layout = nx.drawing.layout.kamada_kawai_layout(networkx_graph)

    # Use this dictionary to color specific atoms
    node_color_dict = {
        "C": "grey",
        "H": "silver",
        "O": "red",
        "N": "blue",
        "Cl": "green",
    }
    node_colors = []
    for node in networkx_graph.nodes:
        if node.name in list(node_color_dict.keys()):
            node_colors.append(node_color_dict[node.name])
        else:
            node_colors.append("black")

    # Node sizes determines if looking at just sites. Bonds, angles, dihedrals have edge_weights
    # and edge_colors as identifiers
    if node_sizes:
        nx.draw(
            networkx_graph,
            layout,
            ax,
            node_color=node_colors,
            node_size=node_sizes,
        )
    else:
        nx.draw(
            networkx_graph,
            layout,
            ax,
            node_color=node_colors,
            width=list(edge_weights.values()),
            edge_color=list(edge_colors.values()),
        )

    # Offset positions to place labels
    for atom, pos in layout.items():
        layout[atom] = pos + [0.09, 0]

    # Get a list of labels to plot
    labels = identify_labels(networkx_graph, list_of_labels, atom_name)
    # Plot labels on current figure
    nx.draw_networkx_labels(
        networkx_graph, layout, labels, horizontalalignment="left"
    )
    ax.margins(0.3, 0.3)

    return ax


def identify_labels(networkx_graph, list_of_labels, atom_name=None):
    """Show labels for all or a selected site."""
    # If atom_name specified, only show labels on that site.
    # Otherwise, show labels for every atom from the label list.
    if atom_name:
        list_of_nodes = []
        for node in networkx_graph.nodes:
            if node.name == atom_name:
                list_of_nodes.append(node)
        labels = return_labels_for_nodes(list_of_nodes, list_of_labels)
    else:
        labels = return_labels_for_nodes(
            list(networkx_graph.nodes), list_of_labels
        )

    return labels


def return_labels_for_nodes(list_of_nodes, list_of_labels):
    """Return the label values for specified sites."""
    # Get the label values for the sites specified.
    # labels is a dict of each node and the labels to put there
    labels = {}
    for i, node in enumerate(list_of_nodes):
        node.label = str(i) + ":" + str(node.name)
        for label in list_of_labels:
            if "." in label:
                label1, label2 = label.split(".")
                try:
                    node.label = (
                        node.label
                        + "\n"
                        + str(getattr(getattr(node, label1), label2))
                    )
                except AttributeError:
                    node.label = node.label + "\nNoneType"
            elif label == "charge":
                if isinstance(getattr(node, label), unyt.array.unyt_quantity):
                    node.label = (
                        node.label
                        + "\n"
                        + str(
                            (getattr(node, label) / unyt.electron_charge).round(
                                4
                            )
                        )
                        + " e"
                    )
                else:
                    node.label = node.label + "\nNone"
            elif label == "position":
                if isinstance(
                    getattr(node, label)[0], unyt.array.unyt_quantity
                ):
                    node.label = (
                        node.label
                        + "\n"
                        + str(
                            getattr(node, label).to("angstrom").round(2)
                            * unyt.angstrom
                        )
                    )
                else:
                    node.label = node.label + "\nNone"
            else:
                try:
                    node.label = node.label + "\n" + str(getattr(node, label))
                except AttributeError:
                    node.label = node.label + "\nNoneType"
                if len(node.label) > 12:
                    node.label = "".join(
                        [line + "\n" for line in node.label.split()]
                    )
        labels[node] = node.label

    return labels


def show_parameter_values(topology, list_of_params, checkbox):
    """Show value of a parameter."""
    if checkbox:
        try:
            report_parameter_expression(topology, list_of_params[0])
        except AttributeError:
            print("There are no values for the parameter expression")
    else:
        # TODO: Should be able to remove this blank print statement so that deselecting the
        # checkbox removes the listed parameters.
        print("")

    return


def report_parameter_expression(topology, param):
    """Return formatted parameters for a given edge."""
    # return nicely printed parameters for a given edge.
    if len(param) == 4:
        try:
            for dihedral in list(topology.dihedrals):
                if dihedral.connection_members == (
                    param[0],
                    param[1],
                    param[2],
                    param[3],
                ):
                    print(dihedral.dihedral_type.expression, "\n")
                    print("{:<12} {:<15}".format("Parameter", "Value"))
                    for k, v in dihedral.dihedral_type.parameters.items():
                        print("{:<12} {:<15}".format(k, v))
        except AttributeError:
            print("Dihedral not typed")
    elif len(param) == 3:
        try:
            for angle in list(topology.angles):
                if angle.connection_members == (param[0], param[1], param[2]):
                    print(angle.angle_type.expression, "\n")
                    print("{:<12} {:<15}".format("Parameter", "Value"))
                    for k, v in angle.angle_type.parameters.items():
                        print("{:<12} {:<15}".format(k, v))
        except AttributeError:
            print("Angle not typed")
    else:
        raise ValueError(
            "Parameters are not proper angles or dihedrals. Connection members are missing"
        )

    return


def get_edges(networkx_graph, atom_name1, atom_name2):
    """Create a list of edges that have been selected."""
    # Create a list of the edges that have been selected corresponding to the selected atom_name1
    # and atom_name2. If both are None, then select edges that are missing bond types.
    selectable_dict = {}
    mia_bond_flag = 0
    if atom_name1 and atom_name2:
        for edge in list(networkx_graph.edges):
            if (edge[0].name == atom_name1 and edge[1].name == atom_name2) or (
                edge[0].name == atom_name2 and edge[1].name == atom_name1
            ):
                selectable_dict = create_dict_of_labels_for_edges(
                    selectable_dict, edge
                )
    elif atom_name1:
        for edge in list(networkx_graph.edges):
            if edge[0].name == atom_name1 or edge[1].name == atom_name1:
                selectable_dict = create_dict_of_labels_for_edges(
                    selectable_dict, edge
                )
    else:
        for nodes in list(networkx_graph.edges.items()):
            if nodes[1]["connection"].bond_type is None:
                selectable_dict = create_dict_of_labels_for_edges(
                    selectable_dict, nodes[0]
                )
                mia_bond_flag = 1
        if not mia_bond_flag:
            return print("All bonds are typed")
    # turn the dic selectable dict into a list of tuples.
    return list(selectable_dict.items())


def create_dict_of_labels_for_edges(selectable_dict, edge):
    """Create list of labels for edges."""
    label0 = ""
    label1 = ""
    try:
        label0 = edge[0].atom_type.name
    except AttributeError:
        print("An atomtype for {} is missing".format(edge[0].label))
    try:
        label1 = edge[1].atom_type.name
    except AttributeError:
        print("An atomtype for {} is missing".format(edge[1].label))
    label = label0 + " --- " + label1
    if label in selectable_dict.keys():
        selectable_dict[label].append(edge)
    else:
        selectable_dict[label] = []
        selectable_dict[label].append(edge)
    return selectable_dict


def plot_networkx_bonds(
    networkx_graph,
    atom_name1=None,
    atom_name2=None,
    list_of_labels=["atom_type.name"],
    list_of_bonds=[],
):
    """Plot the bonds of the networkX graph."""
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    # Create dictionaries of edges that correspond red thick lines for the selected bonds
    edge_weights = {}
    edge_colors = {}
    for bond in networkx_graph.edges:
        if bond in list_of_bonds:
            edge_weights[bond] = 5
            edge_colors[bond] = "red"
        else:
            edge_weights[bond] = 1
            edge_colors[bond] = "k"

    ax = plot_networkx_nodes(
        networkx_graph,
        ax,
        edge_weights=edge_weights,
        edge_colors=edge_colors,
        list_of_labels=list_of_labels,
        atom_name=atom_name1,
    )
    return fig, ax


def report_bond_parameters(topology, edge):
    """Return formatted bond parameters."""
    # return nicely printed bond parameters for a given edge.
    for bond in topology.bonds:
        try:
            if bond.connection_members == edge[
                0
            ] or bond.connection_members == (
                edge[0][1],
                edge[0][0],
            ):
                print(bond.bond_type.expression, "\n")
                print("{:<12} {:<15}".format("Parameter", "Value"))
                for k, v in bond.bond_type.parameters.items():
                    print("{:<12} {:<15}".format(k, v))
        except AttributeError:
            print("The bond between {} is missing parameters".format(edge))

    return


def plot_networkx_atomtypes(
    topology, atom_name=None, list_of_labels=["atom_type.name"]
):
    """Get a networkx plot showing the atom types in a topology object.

    Parameters
    ----------
    topology : A gmso.core.topology.Topology object
        This should be a gmso topology object that you want to visualize the atom types
        that have been parameterized.
    atom_name : The atom name which will have larger node sizes.
        When drawing the networkx graph, all atoms with this name will be 3X as large.
        This input will be of type string. To see what atom names are available, use
        for site in topology.sites:
            print(site.name)

    Returns
    -------
    matplotlib.pyplot.figure
        The drawn networkx plot of the topology on a matplotlib editable figures

    matplotlib.pyplot.axis
        The axis information that corresponds to that figure. This output can be
        shown using
        matplotlib.pyplot.show()
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    networkx_graph = to_networkx(topology)

    # larger nodes for the atom selected
    node_sizes = []
    for node in networkx_graph.nodes:
        if node.name == atom_name:
            node_sizes.append(900)
        else:
            node_sizes.append(300)

    ax = plot_networkx_nodes(
        networkx_graph,
        ax,
        edge_weights=None,
        edge_colors=None,
        node_sizes=node_sizes,
        list_of_labels=list_of_labels,
        atom_name=atom_name,
    )

    return (fig, ax)


def highlight_networkx_edges(networkx_graph, list_of_edges):
    """Return edge parameters for the selected networkX edge list."""
    # return edge parameters for the specified networkx edge list.
    edge_weights = {}
    edge_colors = {}
    for edge in networkx_graph.edges:
        edge_weights[edge] = 1
        edge_colors[edge] = "k"
    for edge in networkx_graph.edges:
        if edge in list_of_edges:
            edge_weights[edge] = 5
            edge_colors[edge] = "r"

    return edge_weights, edge_colors
