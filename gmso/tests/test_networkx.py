import networkx as nx
import pytest

from gmso.external.convert_networkx import to_networkx
from gmso.tests.base_test import BaseTest
from gmso.utils.io import has_ipywidgets, has_matplotlib, import_
from gmso.utils.nx_utils import (
    _get_formatted_atom_types_names_for,
    call_interactive_sites,
    get_edges,
    get_networkx_edges,
    highlight_networkx_edges,
    identify_labels,
    plot_networkx_atomtypes,
    plot_networkx_bonds,
    plot_networkx_nodes,
    plot_networkx_params,
    report_bond_parameters,
    report_parameter_expression,
    return_labels_for_nodes,
    select_angles_from_sites,
    select_dihedrals_from_sites,
    select_edges,
    select_edges_on_networkx,
    select_params_on_networkx,
    show_bond_info,
    show_parameter_values,
)

if has_matplotlib:
    plt = import_("matplotlib.pyplot")
if has_ipywidgets:
    widgets = import_("ipywidgets")


@pytest.mark.skipif(not has_ipywidgets, reason="Ipywidgets is not installed")
@pytest.mark.skipif(not has_matplotlib, reason="Matplotlib is not installed")
class TestNetworkx(BaseTest):
    def test_highlight_networkx_edges(self, typed_ethane):
        list(typed_ethane.angles)[0].angle_type = None
        list(typed_ethane.dihedrals)[0].dihedral_type = None
        graph = to_networkx(typed_ethane)
        list_edges = list(graph.edges)[0:3]
        test_edge_weights, test_edge_colors = highlight_networkx_edges(
            graph, list_edges[0:2]
        )
        assert test_edge_weights[list_edges[0]] == 5
        assert test_edge_weights[list_edges[1]] == 5
        assert test_edge_weights[list_edges[2]] == 1
        assert test_edge_colors[list_edges[0]] == "r"
        assert test_edge_colors[list_edges[1]] == "r"
        assert test_edge_colors[list_edges[2]] == "k"

    def test_plot_networkx_atomtypes(self, typed_ethane):
        fig, ax = plot_networkx_atomtypes(typed_ethane, atom_name=None)
        test_fig, test_ax = plt.subplots(1)
        assert isinstance(fig, test_fig.__class__)
        assert isinstance(ax, test_ax.__class__)

    def test_plot_networkx_bonds(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        fig, ax = plot_networkx_bonds(graph)
        test_fig, test_ax = plt.subplots(1)
        assert isinstance(fig, test_fig.__class__)
        assert isinstance(ax, test_ax.__class__)

    def test_select_params_on_networkx(self, typed_ethane, capsys):
        graph = to_networkx(typed_ethane)
        assert (
            len(select_params_on_networkx(graph, [None, None, None, None])) == 0
        )
        assert len(select_params_on_networkx(graph, ["C", "H", "H"])) == 1
        assert len(select_params_on_networkx(graph, ["C", "C", "H", "H"])) == 1
        assert len(select_params_on_networkx(graph, [None, None, None])) == 0
        assert len(select_params_on_networkx(graph, ["C", "H", None])) == 3
        assert len(select_params_on_networkx(graph, ["C", None, None])) == 3
        assert len(select_params_on_networkx(graph, ["C", "C", "H", None])) == 1
        assert (
            len(select_params_on_networkx(graph, ["C", "C", None, None])) == 1
        )

    def test_select_params_on_networkx_output(self, typed_ethane, capsys):
        graph = to_networkx(typed_ethane)
        select_params_on_networkx(graph, [None, None, None])
        captured, err = capsys.readouterr()
        assert captured.startswith("All angles")
        select_params_on_networkx(graph, [None, None, None, None])
        captured, err = capsys.readouterr()
        assert captured.startswith("All dihedrals")
        for node, angles in graph.nodes(data="angles"):
            if angles[0]:
                angles[0].angle_type = None
        result = select_params_on_networkx(graph, [None, None, None])
        captured, err = capsys.readouterr()
        assert captured.startswith("Since no sites are input, angles")
        for node, dihedrals in graph.nodes(data="dihedrals"):
            if dihedrals[0]:
                dihedrals[0].dihedral_type = None
        result = select_params_on_networkx(graph, [None, None, None, None])
        captured, err = capsys.readouterr()
        assert captured.startswith("Since no sites are input, dihedrals")
        nx.set_node_attributes(graph, None, name="angles")
        result = select_params_on_networkx(graph, [None, None, None])
        captured, err = capsys.readouterr()
        assert captured.startswith("No angle")
        nx.set_node_attributes(graph, None, name="dihedrals")
        result = select_params_on_networkx(graph, [None, None, None, None])
        captured, err = capsys.readouterr()
        assert captured.startswith("No dihedral")
        result = select_params_on_networkx(graph, [None, None])
        captured, err = capsys.readouterr()
        assert captured.startswith("invalid")

    def test__get_formatted_atom_types_names_for(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        for node, dihedrals in graph.nodes(data="angles"):
            assert isinstance(
                _get_formatted_atom_types_names_for(dihedrals[0]), str
            )

    def test_get_networkx_edges(self, typed_ethane, capsys):
        assert (
            len(
                get_networkx_edges(
                    [list(typed_ethane.dihedrals[0].connection_members)]
                )
            )
            == 6
        )
        assert (
            len(
                get_networkx_edges(
                    [list(typed_ethane.angles[0].connection_members)]
                )
            )
            == 4
        )
        with pytest.raises(ValueError):
            assert get_networkx_edges(["C", "C"])

    def test_identify_labels(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        assert len(identify_labels(graph, ["name"], atom_name="C")) == 2

    def test_show_parameter_values(self, typed_ethane):
        parameters = list(typed_ethane.angles[0].connection_members[0:2])
        with pytest.raises(ValueError):
            show_parameter_values(typed_ethane, [parameters], True)

    def test_select_dihedrals_from_sites(self, typed_ethane, capsys):
        graph = to_networkx(typed_ethane)
        select_dihedrals_from_sites(graph, typed_ethane)
        select_dihedrals_from_sites(graph, top, "C", "C", "H", "H")
        captured, err = capsys.readouterr()
        assert isinstance(err, str)

    def test_select_dihedrals_from_sites(self, typed_ethane, capsys):
        graph = to_networkx(typed_ethane)
        select_dihedrals_from_sites(graph, typed_ethane)
        captured, err = capsys.readouterr()
        assert isinstance(err, str)

    def test_plot_networkx_nodes(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        fig, ax = plt.subplots(1, 1)
        plot_networkx_nodes(
            graph, ax, edge_weights={1: 5}, edge_colors={1: "r"}
        )

    def test_plot_networkx_params(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        list_of_edges = get_networkx_edges(
            list_of_params=[
                typed_ethane.dihedrals[0].connection_members,
                typed_ethane.dihedrals[1].connection_members,
            ]
        )
        fig, ax = plot_networkx_params(graph, list_of_edges)
        test_fig, test_ax = plt.subplots(1)

        assert isinstance(fig, test_fig.__class__)
        assert isinstance(ax, test_ax.__class__)

    def test_select_edges_on_networkx(self, typed_ethane, capsys):
        graph = to_networkx(typed_ethane)
        edges = select_params_on_networkx(graph, ["C", "C", "H"])
        select_edges_on_networkx(graph, typed_ethane, edges[0][1])
        captured, err = capsys.readouterr()
        assert isinstance(err, str)

    def test_report_parameter_expression(self, typed_ethane, capsys):
        report_parameter_expression(
            typed_ethane, list(typed_ethane.dihedrals[0].connection_members)
        )
        report_parameter_expression(
            typed_ethane, list(typed_ethane.angles[0].connection_members)
        )
        captured, err = capsys.readouterr()
        assert isinstance(err, str)

    def test_get_edges(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        assert len(get_edges(graph, "C", "C")[0][1]) == 1
        assert len(get_edges(graph, "H", None)[0][1]) == 6
        assert not get_edges(graph, None, None)

    def test_report_bond_parameters(self, typed_ethane, capsys):
        report_bond_parameters(
            typed_ethane, [typed_ethane.bonds[0].connection_members]
        )
        captured, err = capsys.readouterr()
        assert isinstance(err, str)

    def test_return_labels_for_nodes(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        assert (
            len(
                return_labels_for_nodes(
                    graph.nodes, ["atom_type.name", "charge", "positions"]
                )
            )
            == 8
        )
        assert (
            list(
                return_labels_for_nodes(
                    graph.nodes, ["atom_type.error"]
                ).values()
            )[0][-8:]
            == "NoneType"
        )
        assert (
            list(return_labels_for_nodes(graph.nodes, ["error"]).values())[0][
                -8:
            ]
            == "NoneType"
        )

    def test_select_angles_from_sites(self, typed_ethane, capsys):
        graph = to_networkx(typed_ethane)
        select_angles_from_sites(
            graph, typed_ethane, Atom1="C", Atom2="H", Atom3="C"
        )
        select_angles_from_sites(
            graph, typed_ethane, Atom1="O", Atom2="H", Atom3="C"
        )
        captured, err = capsys.readouterr()
        assert isinstance(err, str)

    def test_call_interactive_sites(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        assert not call_interactive_sites(
            "C", "C", graph, typed_ethane, list_of_labels=["name"]
        )
        assert not call_interactive_sites(
            "C", "O", graph, typed_ethane, list_of_labels=["name"]
        )

    def test_select_edges(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        list_of_edges = get_edges(graph, "C", "C")
        assert not select_edges(graph, typed_ethane, list_of_edges, ["name"])

    def test_show_bond_info(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        list_of_edges = get_edges(graph, "C", "C")
        assert not show_bond_info(True, typed_ethane, list_of_edges)
        assert not show_bond_info(False, typed_ethane, list_of_edges)
