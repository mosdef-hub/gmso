import unyt as u

import gmso
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType


def read_itp(itp_file):
    """Create a topology from a provided gro file.
    The Gromos87 (gro) format is a common plain text structure file used
    commonly with the GROMACS simulation engine.  This file contains the
    simulation box parameters, number of atoms, the residue and atom number for
    each atom, as well as their positions and velocities (velocity is
    optional).  This method will receive a filepath representation either as a
    string, or a file object and return a `topology`.
    Parameters
    ----------
    filename : str or file object
        The path to the gro file either as a string, or a file object that
        points to the gro file.
    Returns
    -------
    gmso.core.topology
        A `topology` object containing site information
    Notes
    -----
    Gro files do not specify connections between atoms, the returned topology
    will not have connections between sites either.
    Currently this implementation does not support parsing velocities from a gro file o
    with more than 1 frame.
    All residues and resid information from the gro file are currently lost
    when converting to `topology"""

    # for different functional forms for bonds, angles, dihedrals include an if statement for angle
    # pdb=pmd.load_file('liquid.pdb')
    # pdb=from_parmed(pdb,refer_type=True)

    # for different functional forms for bonds, angles, dihedrals include an if statement for angle
    # pdb=pmd.load_file('liquid.pdb')
    # pdb=from_parmed(pdb,refer_type=True)

    mass = []
    charge = []
    atype = []
    epsilon = []
    sigma = []

    dicts_atype = []

    # Bond parameters
    b_1 = []
    b_2 = []
    b_K = []
    b_d = []

    # Get atom types, epsilon and sigma in a dictionary
    with open(itp_file, "r") as file:
        for line in file:
            if "atomtypes" in line:
                for line in file:
                    if "molecule" in line:
                        break
                    elif line.split() and "name" not in line.split():
                        dicts_atype.append(
                            {
                                "type": line.split()[0],
                                "epsilon": float(line.split()[5]),
                                "sigma": float(line.split()[6]),
                            }
                        )
                        pass
                        # print(line)

    # Get charges, mass and atomtype in a list
    with open("LIQ.itp", "r") as file:
        for line in file:
            # Atoms

            if "atoms" in line:
                k = 1
                for line in file:
                    if "[" in line:
                        break
                    elif line.split() and line.split()[0] == str(k):
                        # print(line.split())
                        # print(line)
                        # print(line.split())
                        atype.append(line.split()[1])
                        mass.append(line.split()[7])
                        charge.append(line.split()[6])
                        k = k + 1
            # bonds

            if "bonds" in line:
                kbonds = 1
                for line in file:
                    if "[" in line:
                        break
                    elif line.split() and "funct" not in line.split():
                        # print(line)

                        kbonds = kbonds + 1

            # Angles

            if "angles" in line:
                kang = 1
                for line in file:
                    if "[" in line:
                        break
                    elif line.split() and "funct" not in line.split():
                        # print(line)

                        kang = kang + 1

            # Dihedrals

            if "dihedrals" in line:
                kdih = 1
                for line in file:
                    if "[" in line:
                        break
                    elif line.split() and "funct" not in line.split():
                        # print(line)

                        kdih = kdih + 1
        pass
    print(k, kbonds, kang, kdih)
    pdb = gmso.Topology()
    for i in range(len(mass)):
        site = Atom()
        site.mass = float(mass[i])
        site.charge = float(charge[i])

        for index in range(len(dicts_atype)):
            # Loop over dictionary to assign epsilon
            for key in dicts_atype[index]:
                # print(index,key,dicts_atype[index][key])
                if dicts_atype[index]["type"] == atype[i]:
                    sigma = dicts_atype[index]["sigma"]
                    epsilon = dicts_atype[index]["epsilon"]
                    break
                break
        pass
        site.atom_type = AtomType(
            name=atype[i],
            charge=float(charge[i]),
            mass=float(mass[i]),
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            parameters={
                "sigma": sigma * u.angstrom,
                "epsilon": epsilon * u.Unit("kcal / mol"),
            },
            independent_variables={"r"},
        )
        pdb.add_site(site)
        print(i, site.atom_type.parameters)

    for i in range(kbonds):
        bonds = Bond(connection_members=(pdb.sites[0], pdb.sites[1]))
        bonds.bond_type = (
            BondType(
                name="2",
                expression="0.5 * k * (r-r_eq)**2",
                independent_variables={"r"},
                parameters={"k": 1000 * u.Unit("kJ / (nm**2)"), "r_eq": 0.14 * u.nm},
            ),
        )

        pdb.add_connection(bonds)

    for i in range(kbonds):
        bonds = Bond(connection_members=(pdb.sites[0], pdb.sites[1]))
        bonds.bond_type = (
            BondType(
                name="2",
                expression="0.5 * k * (r-r_eq)**2",
                independent_variables={"r"},
                parameters={"k": 1000 * u.Unit("kJ / (nm**2)"), "r_eq": 0.14 * u.nm},
            ),
        )

        pdb.add_connection(bonds)

    # bond=gmso.core.bond.Bond(connection_members=(pdb.sites[0],pdb.sites[1]))
    return pdb
