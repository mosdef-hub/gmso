"""Convert to and from a TRIPOS mol2 file."""
#TODO add sources of mol2 files
import os

import unyt as u
import warnings

from gmso import Topology, Atom, Bond
from gmso.core.element import element_by_name

def from_mol2(filename): #TODO add flags for information to return
    #TODO: descriptions and examples
    #TODO: Be sure to be clear about how to read in to mbuild compound using gmso.external.to_mbuild function 
    
    msg = "Provided path to file that does not exist"
    if not os.path.isfile(filename):
        raise OSError(msg)
    #TODO: write a function that verifies a file path is a mol2 file

    # Initialize topology
    topology = Topology(name=os.path.splitext(os.path.basename(filename))[0])
    #save the name from the filename
    f = open(filename, 'r')
    line = f.readline()
    while f:
        # check for header character in line
        if line.startswith('@<TRIPOS>'):
            # if header character in line, send to a function that will direct it properly
            line, topology = parse_record_type_indicator(f, line, topology)
        elif line == "":
            break
        else: 
            # else, skip to next line
            line = f.readline()
    f.close()

    # return warnings if any sections are not covered
    # save sections to a list for each
    # Iterate through list of <ATOM> to save to Topology.sites
    # Iterate through list of <BOND> to save to Topology.bonds
    # Save box dimensions to Topology.box
    # Make sure other attributes of the topology are updated accordingly
    #TODO: read in parameters to correct attribute as well
    return topology

def load_top_sites(f, topology):
    """Take a mol2 file section with the heading @<TRIPOS>ATOM and save to the topology.sites attribute"""
    while True:
        line = f.readline()
        if '@' not in line:
            line = line.split()
            position = [float(x) for x in line[2:5]] * u.Å
            # TODO: make sure charges are also saved as a unyt value
            # TODO: add validation for element names
            atom = Atom(name=line[1], position=position.to('nm'), charge=float(line[8]), element=element_by_name(line[5]))
            topology.add_site(atom)
        else:
            break
    return line, topology

def load_top_bonds(f, topology):
    """Take a mol2 file section with the heading @<TRIPOS>BOND and save to the topology.bonds attribute"""
    while True:
        line = f.readline()
        if '@' not in line:
            line = line.split()
            bond = Bond(connection_members=(topology.sites[int(line[1])-1], topology.sites[int(line[2])-1]))
            topology.add_connection(bond)
        else:
            break
    return line, topology

def load_top_box(f, topology):
    """Take a mol2 file section with the heading @<TRIPOS>FF_PBC and save to a topology"""
    while True:
        line = f.readline()
        if '@' not in line:
            line = line.split()
            #TODO: write to box information
        else:
            break
    return line, topology

def parse_record_type_indicator(f, line, topology):
    """Take a specific record type indicator from a mol2 file format and save to the proper attribute of a gmso topology.
    Supported record type indicators include Atom, Bond, FF_PBC."""
    supported_rti = {'@<TRIPOS>ATOM\n':load_top_sites,
                     '@<TRIPOS>BOND\n':load_top_bonds,
                     '@<TRIPOS>FF_PBC\n':load_top_box}
    #read in to atom attribute
    try:
        return supported_rti[line](f, topology)
    except KeyError:
        warnings.warn('The record type indicator {} is not supported'.format(line))
        line = f.readline()
        return line, topology
