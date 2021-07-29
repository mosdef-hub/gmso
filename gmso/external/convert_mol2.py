"""Convert to and from a TRIPOS mol2 file."""
#TODO add sources of mol2 files
import unyt as u
import warnings

from gmso import Topology, Atom, Bond
from gmso.core.element import element_by_name

def from_mol2(filename): #TODO add flags for information to return
    #TODO descriptions and examples
    msg = "Provided path to file is not a support mol2 file"
    #TODO write a function that verifies a file path is a mol2 file
    #TODO implement reading in and saving to an empty topology
        # Initialize topology
    topology = Topology(name='Mol2_top')
        # Read in filename, return error if filename not found
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
    #TODO read in parameters to correct attribute as well
    return topology
    #TODO Be sure to be clear about how to read in to mbuild compound using gmso.external.to_mbuild function 

def load_top_sites(f, topology):
    """Take a mol2 file section with the heading @<TRIPOS>ATOM and save to the topology.sites attribute"""
    while True:
        line = f.readline()
        if '@' not in line:
            line = line.split()
            atom = Atom(name=line[1], position=line[2:5], charge=line[8], element=element_by_name(line[5]))
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
        
        

def to_mol2(topology, openmm_object="topology"):
    """Convert an untyped topology object to an untyped OpenMM modeller or topology.

    This is useful if it's preferred to atom-type a system within OpenMM.
    See http://openmm.org for more information.

    Parameters
    ----------
    topology : `Topology` object
        An untyped topology object
    open_mm_object : 'topology' or 'modeller' OpenMM object, default='topology'
        Untyped OpenMM object to convert to

    Returns
    -------
    open_mm_object : Untyped `topology` or `modeller` object

    """
    openmm_top = app.Topology()

    # Get topology.positions into OpenMM form
    openmm_unit = 1 * simtk_unit.nanometer
    topology.positions.convert_to_units(openmm_unit.unit.get_symbol())
    value = [i.value for i in topology.positions]
    openmm_pos = simtk_unit.Quantity(value=value, unit=openmm_unit.unit)

    # Adding a default chain and residue temporarily
    chain = openmm_top.addChain()
    residue = openmm_top.addResidue(name="RES", chain=chain)

    for site in topology.sites:
        openmm_top.addAtom(
            name=site.name, element=site.element.name, residue=residue
        )

    # Set box
    box = topology.box
    box.lengths.convert_to_units(u.nanometer)
    lengths = box.lengths.value
    openmm_top.setUnitCellDimensions(lengths)

    # TODO: Figure out how to add residues
    # TODO: Convert connections to OpenMM Bonds

    if openmm_object == "topology":

        return openmm_top

    else:
        modeller = app.Modeller(openmm_top, openmm_pos)

        return modeller


