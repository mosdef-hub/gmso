import numpy as np
import unyt as u
import matplotlib.pyplot as plt
import networkx as nx

from gmso.core.topology import Topology
from gmso.external.convert_networkx import to_networkx

def plot_networkx_atomtypes(topology,atom_name=None):
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

    fig,ax = plt.subplots(1,1,figsize=(8,8))
    networkx_graph = to_networkx(topology)
    node_sizes = []
    for node in networkx_graph.nodes:
        if node.name == atom_name:
            node_sizes.append(900)
        else:
            node_sizes.append(300)
    ax = plot_nodes(networkx_graph,ax,edge_weights=None,edge_colors=None,node_sizes = node_sizes)

    return(fig,ax)

def plot_networkx_bonds(topology,atom_name1=None,atom_name2=None):
    """Get a networkx plot showing the bonds in a topology object.

    Parameters
    ----------
    topology : A gmso.core.topology.Topology object
        This should be a gmso topology object that you want to visualize the atom types
        that have been parameterized.
    atom_name1 : The atom name to of which bonds you want to have selected
        When drawing the networkx graph, all bonds connected to atoms with
            this name will be indicated.
        This input will be of type string. To see what atom names are available, use
            for site in topology.sites:
                print(site.name)
        If no atom_name is given, then only bonds missing bond information will be 
            indicated
    atom_name2 : A second atom name to restrict what bonds will be selected. only bonds
                     between the two given atoms will be indicated.
                 This input will be of type string.

    Returns
    -------
    matplotlib.pyplot.figure
        The drawn networkx plot of the topology on a matplotlib editable figures

    matplotlib.pyplot.axis
        The axis information that corresponds to that figure. This output can be
        shown using
        matplotlib.pyplot.show()
    """

    #Color and weight edges between particular atoms. If both atom_names are none, plot missing bond types.
    from gmso.external import convert_networkx
    
    fig,ax = plt.subplots(1,1,figsize=(8,8))
    networkx_graph = to_networkx(topology)

    edge_weights = {}
    edge_colors = {}
    mia_bond_ind = 0
    if atom_name1 and atom_name2:
        for edge in networkx_graph.edges:
            if edge[0].name == atom_name1 and edge[1].name == atom_name2:
                edge_weights[edge] = 5
                edge_colors[edge] = 'red'
            elif edge[0].name == atom_name2 and edge[1].name == atom_name1:
                edge_weights[edge] = 5
                edge_colors[edge] = 'red'
            else:
                edge_weights[edge] = 1
                edge_colors[edge] = 'k'
    elif atom_name1:
        for edge in networkx_graph.edges:
            if edge[0].name == atom_name1 or edge[1].name == atom_name1:
                edge_weights[edge] = 5
                edge_colors[edge] = 'red'
            else:
                edge_weights[edge] = 1
                edge_colors[edge] = 'k'
    else:
        for bond in list(networkx_graph.edges.items()):
            if bond[1]['connection'].bond_type == None:
                edge_weights[bond[0]] = 5
                edge_colors[bond[0]] = 'red'
                mia_bond_ind = 1
            else:
                edge_weights[bond[0]] = 1
                edge_colors[bond[0]] = 'k'
        if not mia_bond_ind:
            print('All bonds are typed')

    ax = plot_nodes(networkx_graph,ax,edge_weights,edge_colors)
    return fig, ax


def plot_networkx_angles(topology,center_atom_index = None):
    """Get a networkx plot showing the angles in a topology object.
    
    Parameters 
    ----------
    topology : A gmso.core.topology.Topology object
        This should be a gmso topology object that you want to visualize the angle types
        that have been parameterized.
    center_atom_index : The central atom to of which angles you want to have selected
        When drawing the networkx graph, all angles connected to atoms with
            this index will be indicated.
        This input will be of type int. To see what atoms correspond to these indices, see 
            gmso.formats.networkx.plot_networkx_atomtypes
        If no atom_name is given, then only atoms missing angle information will be 
            indicated
    
    Returns
    ------- 
    matplotlib.pyplot.figure
        The drawn networkx plot of the topology on a matplotlib editable figures
        
    matplotlib.pyplot.axis
        The axis information that corresponds to a networkx drawn figure. This output can be
        shown using
        matplotlib.pyplot.show()
    """
    networkx_graph = to_networkx(topology)
    
    fig,ax = plt.subplots(1,1,figsize=(8,8))
    
    edge_weights, edge_colors = highlight_bonds(networkx_graph,'angles',atom_index=center_atom_index)

    ax = plot_nodes(networkx_graph,ax,edge_weights,edge_colors)

    return(fig,ax)


def plot_networkx_dihedrals(topology,center_atom_index = None):
    """Get a networkx plot showing the dihedrals in a topology object.
    
    Parameters
    ----------
    topology : A gmso.core.topology.Topology object
        This should be a gmso topology object that you want to visualize the dihedral types
        that have been parameterized.
    center_atom_index : The second or third atom of the dihedrals you want to have selected
        When drawing the networkx graph, all dihedrals connected to atoms with
            this index will be indicated. 
        This input will be of type int. To see what atoms correspond to these indices, see 
            gmso.formats.networkx.plot_networkx_atomtypes
        If no atom_name is given, then only atoms missing dihedral information will be
            indicated
    
    Returns
    -------
    matplotlib.pyplot.figure
        The drawn networkx plot of the topology on a matplotlib editable figures

    matplotlib.pyplot.axis
        The axis information that corresponds to that figure. This output can be
        shown using
        matplotlib.pyplot.show()    
    """
    networkx_graph = to_networkx(topology)
    
    fig,ax = plt.subplots(1,1,figsize=(8,8))
    
    edge_weights, edge_colors = highlight_bonds(networkx_graph,'dihedrals',atom_index=center_atom_index)
    ax = plot_nodes(networkx_graph,ax,edge_weights,edge_colors)

    return(fig,ax)


def highlight_bonds(networkx_graph,attribute,atom_index=None):
    edge_weights={};edge_colors={}
    for edge in networkx_graph.edges:
        edge_weights[edge] = 1; edge_colors[edge] = 'k'
   
    def_param = 1
    if atom_index == None:
        for node in networkx_graph.nodes:
            for parameter in networkx_graph.nodes[node][attribute]:
                var = attribute[:-1] + '_type'
                if getattr(parameter,var) == None:
                    def_param = 0
                    members = list(parameter.connection_members)
                    for i in np.arange(len(members)-1):
                        edge_weights[(members[i],members[i+1])] = 5; edge_weights[(members[i+1],members[i])] = 5
                        edge_colors[(members[i],members[i+1])] = 'red'; edge_colors[(members[i+1],members[i])] = 'red'

        if def_param:
            print('No {} selected, and all {} typed'.format(attribute,attribute))

    elif isinstance(atom_index,int) and len(list(networkx_graph.nodes)) > atom_index:
        node = list(networkx_graph.nodes)[atom_index]
        for parameter in networkx_graph.nodes[node][attribute]:
            if (node == parameter.connection_members[1] or 
               node == parameter.connection_members[
                            len(networkx_graph.nodes[node][attribute][0].connection_members)-2]):
                for i in np.arange(len(networkx_graph.nodes[node][attribute][0].connection_members)-1):
                    node1 = list(parameter.connection_members)[i+1]
                    node0 = list(parameter.connection_members)[i]
                    edge = (node0,node1)
                    edge_weights[edge] = 5
                    edge_colors[edge] = 'red'
                    edge = (node1,node0)
                    edge_weights[edge] = 5
                    edge_colors[edge] = 'red'
    else:
        print('Invalid input for atom or node index')
        
    return edge_weights, edge_colors

def plot_nodes(networkx_graph,ax,edge_weights=None,edge_colors=None,node_sizes = None):
    pos={}
    for node in networkx_graph.nodes:
        pos[node] = node.position.value[0:2]

    layout = nx.drawing.layout.spring_layout(networkx_graph,k=.5,pos=pos)
    
    node_color_dict = {'C':'grey','H':'silver','O':'red','N':'blue','Cl':'green'}
    node_colors = []
    for node in networkx_graph.nodes:
        if node.name in list(node_color_dict.keys()):
            node_colors.append(node_color_dict[node.name])
        else:
            node_colors.append('black')

    if node_sizes:
        nx.draw(networkx_graph,layout,ax,node_color=node_colors,node_size=node_sizes)
    else:
        nx.draw(networkx_graph,layout,ax,node_color=node_colors,
            width=list(edge_weights.values()),edge_color=list(edge_colors.values()))
    labels = {}
    i=0
    for node in list(networkx_graph.nodes()):
        node.label = str(i) + ': ' + node.name + '\n' + node.atom_type.name
        labels[node] = node.label
        i+=1

    for atom,pos in layout.items():
        layout[atom] = pos + [0.09,0]
    nx.draw_networkx_labels(networkx_graph,layout,labels,horizontalalignment='left')
    ax.margins(.3,.3)
    
    return ax
