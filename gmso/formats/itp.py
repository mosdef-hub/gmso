import unyt as u

import gmso
from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.forcefield import *


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

    # for different functional forms for bonds, angles, dihedrals include an if statement for angle
    # pdb=pmd.load_file('liquid.pdb')
    # pdb=from_parmed(pdb,refer_type=True)

    mass = []
    charge = []
    atype = []
    epsilon = []
    sigma = []

    dicts_atype = []

    # collect bond parameters
    bind_1 = []
    bind_2 = []
    b_type = []
    b_K = []
    b_req = []

    # collect angle parameters
    aind_1 = []
    aind_2 = []
    aind_3 = []
    a_type = []
    a_K = []
    a_thetaeq = []

    # collect dihedral parameters
    dind_1 = []
    dind_2 = []
    dind_3 = []
    dind_4 = []
    d_type = []
    d_K = []
    d_phi = []
    d_n = []

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
    with open(itp_file, "r") as file:
        for line in file:
            # Atoms

            if "atoms" in line:
                natoms = 1
                for line in file:
                    if "[" in line:
                        break
                    elif line.split() and line.split()[0] == str(natoms):
                        # print(line.split())
                        # print(line)
                        # print(line.split())
                        atype.append(line.split()[1])
                        mass.append(line.split()[7])
                        charge.append(line.split()[6])
                        natoms = natoms + 1
            # bonds

            if "bonds" in line:
                nbonds = 0
                for line in file:
                    if "[" in line:
                        break
                    elif line.split() and "funct" not in line.split():
                        print(line)
                        bind_1.append(int(line.split()[0]) - 1)
                        bind_2.append(int(line.split()[1]) - 1)
                        b_type.append(int(line.split()[2]))
                        b_K.append(float(line.split()[3]))
                        b_req.append(float(line.split()[4]))
                        nbonds = nbonds + 1

            # Angles

            if "angles" in line:
                nang = 0
                for line in file:
                    if "[" in line:
                        break
                    elif line.split() and "funct" not in line.split():
                        print(line)
                        aind_1.append(int(line.split()[0]) - 1)
                        aind_2.append(int(line.split()[1]) - 1)
                        aind_3.append(int(line.split()[2]) - 1)
                        a_type.append(int(line.split()[3]))
                        a_K.append(float(line.split()[4]))
                        a_thetaeq.append(float(line.split()[5]))
                        nang = nang + 1

            # Dihedrals
            if "dihedrals" in line:
                ndih = 0
                for line in file:
                    if "[" in line:
                        break
                    elif line.split() and "funct" not in line.split():
                        # print(line)
                        dind_1.append(int(line.split()[0]) - 1)
                        dind_2.append(int(line.split()[1]) - 1)
                        dind_3.append(int(line.split()[2]) - 1)
                        dind_4.append(int(line.split()[3]) - 1)
                        d_type.append(int(line.split()[4]))
                        d_K.append(float(line.split()[5]))
                        d_phi.append(float(line.split()[6]))
                        d_n.append(float(line.split()[7]))
                        ndih = ndih + 1

    print(natoms, nbonds, nang, ndih)
    pdb = gmso.Topology()
    for i in range(natoms - 1):
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
        # print(i,site.atom_type.parameters)

    print(len(pdb.sites))
    print(bind_1)
    print(bind_2)

    for i in range(nbonds):
        site_1 = pdb.sites[bind_1[i]]
        site_2 = pdb.sites[bind_2[i]]
        bonds = Bond(connection_members=(site_1, site_2))

        K = b_K[i]
        req = b_req[i]

        bonds.bond_type = BondType(
            name="HarmonicBondPotential",
            expression="0.5 * k * (r-r_eq)**2",
            independent_variables={"r"},
            member_types=(site_1.atom_type.name, site_2.atom_type.name),
            member_classes=(site_1.atom_type.atomclass, site_2.atom_type.atomclass),
            parameters={"k": K * u.Unit("kJ / (nm**2)"), "r_eq": req * u.nm},
        )

        pdb.add_connection(bonds)

    for i in range(nang):
        site_1 = pdb.sites[aind_1[i]]
        site_2 = pdb.sites[aind_2[i]]
        site_3 = pdb.sites[aind_3[i]]
        angles = Angle(connection_members=(site_1, site_2, site_3))

        K = a_K[i]
        thetaeq = a_thetaeq[i]

        angles.angle_type = AngleType(
            name="HarmonicAnglePotential",
            expression="0.5 * k * (theta-theta_eq)**2",
            independent_variables={"theta"},
            member_types=(
                site_1.atom_type.name,
                site_2.atom_type.name,
                site_3.atom_type.name,
            ),
            member_classes=(
                site_1.atom_type.atomclass,
                site_2.atom_type.atomclass,
                site_3.atom_type.atomclass,
            ),
            parameters={"k": K * u.Unit("kJ / (deg**2)"), "theta_eq": thetaeq * u.deg},
        )

        pdb.add_connection(angles)

    for i in range(ndih):
        site_1 = pdb.sites[dind_1[i]]
        site_2 = pdb.sites[dind_2[i]]
        site_3 = pdb.sites[dind_3[i]]
        site_4 = pdb.sites[dind_4[i]]
        diehdrals = Dihedral(connection_members=(site_1, site_2, site_3, site_4))
        K = d_K[i]
        req = d_phi[i]
        n = d_n[i]
        diehdrals.dihedral_type = DihedralType(
            name="HarmonicTorsionPotential",
            expression="k * (1 + cos(n * phi - phi_eq))**2",
            independent_variables={"phi"},
            member_types=(
                site_1.atom_type.name,
                site_2.atom_type.name,
                site_3.atom_type.name,
                site_4.atom_type.name,
            ),
            member_classes=(
                site_1.atom_type.atomclass,
                site_2.atom_type.atomclass,
                site_3.atom_type.atomclass,
                site_4.atom_type.atomclass,
            ),
            parameters={
                "k": K * u.Unit("kJ / (deg**2)"),
                "phi_eq": 180 * u.deg,
                "n": n * u.dimensionless,
            },
        )

        pdb.add_connection(diehdrals)
    #
    return pdb 