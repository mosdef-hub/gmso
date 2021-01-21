import unyt as u
import numpy as np
#import matplotlib.pyplot as plt
#import networkx as nx
import pytest

from gmso.formats.networkx import (plot_networkx_atomtypes,
plot_networkx_bonds, plot_networkx_angles, plot_networkx_dihedrals,
highlight_bonds, sel_input_params, initialize_edge_params)
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
    def test_highlight_bonds(self, typed_ethane):
        list(typed_ethane.angles)[0].angle_type = None
        list(typed_ethane.dihedrals)[0].dihedral_type = None
        
        graph = to_networkx(typed_ethane)
        test_edge_weights, test_edge_colors = highlight_bonds(graph, 'angles')
        nodes = list(graph.nodes)
        assert test_edge_weights[nodes[0],nodes[4]] == 5
        assert test_edge_weights[nodes[4],nodes[5]] == 5
        assert test_edge_weights[nodes[0],nodes[1]] == 1
  
        test_edge_weights, test_edge_colors = highlight_bonds(graph, 'dihedrals')
        assert test_edge_weights[nodes[0],nodes[4]] == 5
        assert test_edge_weights[nodes[4],nodes[5]] == 5
        assert test_edge_weights[nodes[0],nodes[1]] == 5
        assert test_edge_weights[nodes[0],nodes[3]] == 1

             
   
    def test_sel_input_params(self,typed_ethane):
        graph = to_networkx(typed_ethane)
        edge_weights, edge_colors = initialize_edge_params(graph)
        test_edge_weights, test_edge_colors = sel_input_params(graph, 0, 'angles', edge_weights, edge_colors)

        assert all([a == b for a, b in zip(list(test_edge_weights.values()), [5, 5, 5, 5, 1, 1, 1])])
        assert all([a == b for a, b in zip(list(test_edge_colors.values()), ['red', 'red', 'red', 'red', 'k', 'k', 'k'])])

    def test_initialize_edge_params(self, typed_ethane):
        graph = to_networkx(typed_ethane)
        test_edge_weights, test_edge_colors = initialize_edge_params(graph)

        assert all([a == b for a, b in zip(list(test_edge_weights.values()), [1, 1, 1, 1, 1, 1, 1])])
        assert all([a == b for a, b in zip(list(test_edge_colors.values()), ['k', 'k', 'k', 'k', 'k', 'k', 'k'])])

    def test_plot_networkx_atomtypes(self,typed_ethane):
        graph = to_networkx(typed_ethane)
        fig, ax = plot_networkx_atomtypes(typed_ethane,atom_name=None)
        test_fig, test_ax = plt.subplots(1)
        
        assert isinstance(fig, test_fig.__class__)
        assert isinstance(ax, test_ax.__class__)

    def test_plot_networkx_bonds(self,typed_ethane):
        graph = to_networkx(typed_ethane)
        fig, ax = plot_networkx_bonds(typed_ethane)
        test_fig, test_ax = plt.subplots(1)

        assert isinstance(fig, test_fig.__class__)
        assert isinstance(ax, test_ax.__class__)
