"""Convert to and from a TRIPOS mol2 file."""
#TODO add sources of mol2 files
import unyt as u

from gmso.core import Topology

def from_mol2(filename): #TODO add flags for information to return
    #TODO descriptions and examples
    #TODO implement reading in and saving to an empty topology
        # Initialize topology
        # Read in filename, return error if filename not found
        with open(filename, 'r') as f :
            line = f.readline()
            # check for header character in line
            # if header character in line, send to a function that will direct it properly
            parse_record_type_indicator(f, line, topology)
            # else, skip the line
	    # return warnings if any sections are not covered
	    # save sections to a list for each
        # Iterate through list of <ATOM> to save to Topology.sites
        # Iterate through list of <BOND> to save to Topology.bonds
        # Save box dimensions to Topology.box
        # Make sure other attributes of the topology are updated accordingly
        #TODO read in parameters to correct attribute as well
        # Return gmso.Topology
        #TODO Be sure to be clear about how to read in to mbuild compound using gmso.external.to_mbuild function 

def parse_record_type_indicator(
        
        

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


