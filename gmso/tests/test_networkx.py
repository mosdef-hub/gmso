import unyt as u
import numpy as np
#import matplotlib.pyplot as plt
#import networkx as nx
import pytest

from gmso.formats.networkx import *
from gmso.tests.base_test import BaseTest
from unyt.testing import assert_allclose_units
from gmso.utils.io import import_, has_networkx, has_pyplot
from gmso.external.convert_networkx import to_networkx

if has_networkx:
    networkx = import_('networkx')

if has_pyplot:
    plt = import_('matplotlib.pyplot')

@pytest.mark.skipif(not has_networkx, reason="Networkx is not installed")
@pytest.mark.skipif(not has_pyplot, reason="Matplotlib.pyplot is not installed")
class TestNetworkx(BaseTest):
    def test_highlight_networkx_edges(self, typed_ethane):
        list(typed_ethane.angles)[0].angle_type = None
        list(typed_ethane.dihedrals)[0].dihedral_type = None
        
        graph = to_networkx(typed_ethane)
        test_edge_weights, test_edge_colors = highlight_networkx_edges(graph, 'angles')
        nodes = list(graph.nodes)
        assert test_edge_weights[nodes[0],nodes[4]] == 5
        assert test_edge_weights[nodes[4],nodes[5]] == 5
        assert test_edge_weights[nodes[0],nodes[1]] == 1
  
        test_edge_weights, test_edge_colors = highlight_networkx_edges(graph, 'dihedrals')
        assert test_edge_weights[nodes[0],nodes[4]] == 5
        assert test_edge_weights[nodes[4],nodes[5]] == 5
        assert test_edge_weights[nodes[0],nodes[1]] == 5
        assert test_edge_weights[nodes[0],nodes[3]] == 1

    def test_plot_networkx_atomtypes(self,typed_ethane):
        graph = to_networkx(typed_ethane)
        fig, ax = plot_networkx_atomtypes(typed_ethane,atom_name=None)
        test_fig, test_ax = plt.subplots(1)
        
        assert isinstance(fig, test_fig.__class__)
        assert isinstance(ax, test_ax.__class__)

    def test_plot_networkx_bonds(self,typed_ethane):
        graph = to_networkx(typed_ethane)
        fig, ax = plot_networkx_bonds(graph)
        test_fig, test_ax = plt.subplots(1)

        assert isinstance(fig, test_fig.__class__)
        assert isinstance(ax, test_ax.__class__)

    def test_select_params_on_networkx(self,typed_ethane):
        graph = to_networkx(typed_ethane)
        assert len(select_params_on_networkx(graph,[None,None,None,None])) == 0
        assert len(select_params_on_networkx(graph,['C','H','H'])) == 1
        assert len(select_params_on_networkx(graph,['C','C','H','H'])) == 1
        assert len(select_params_on_networkx(graph,[None,None,None])) == 0

    def test__get_formatted_atom_types_names_for(self,typed_ethane):
        graph = to_networkx(typed_ethane)
        for node, dihedrals in graph.nodes(data='angles'):
            assert isinstance(_get_formatted_atom_types_names_for(dihedrals[0]),str)

    def test_get_networkx_edges(self):
        with pytest.raises(ValueError):
            get_networkx_edges(list_of_params = ['C','C'])

    def test_identify_labels(self,typed_ethane):
        graph = to_networkx(typed_ethane)
        assert len(identify_labels(graph,['name'],atom_name = 'C')) == 2

    def test_show_parameter_values(self,typed_ethane):
        parameters = list(typed_ethane.angles[0].connection_members[0:2])
        with pytest.raises(ValueError):
            show_parameter_values(typed_ethane, [parameters], True) 
