import unyt as u

import gmso
from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType


def read_itp(itp_file):
    """Create a topology from a provided gro file.
    The Gromos87 (gro) format is a common plain text structure file used
    commonly with the GROMACS simulation engine.  This file contains the
    simulation box parameters, number of atoms, the residue and atom number for
    each atom, as well as their positions and velocities (velocity is
    optional).  This only takes a string, to be passed into with open(filename,"r") and returns gmso `topology`.

    Parameters
    ----------
    filename : str
        The path to the gro file either as a string, or a file object that
        points to the gro file.
    Returns
    -------
    gmso.core.topology
        A `topology` object containing site information
    Notes
    -----
    """
    topology = gmso.Topology()
    with open(itp_file, "r") as file:
        # build topology sequentially
        _parse_atoms(file, topology)
        _parse_bonds(file, topology)
        _parse_angles(file, topology)
        _parse_torsions(file, topology)
    return topology



def _parse_atoms(file, pdb):
    """Get a dictionary of parsable atomtypes."""
    mass = []
    charge = []
    atype = []
    dicts_atype = []

    for line in file:
        if "atomtypes" in line:
            for line_1 in file:
                if "molecule" in line_1:
                    break
                elif line_1.split() and ";" not in line_1.split():
                    dicts_atype.append(
                        {
                            "type": line_1.split()[0],
                            "epsilon": float(line_1.split()[5]),
                            "sigma": float(line_1.split()[6]),
                        }
                    )

        if "atoms" in line:
            natoms = 1
            for line_1 in file:
                if "[" in line_1:
                    break
                elif line_1.split() and line_1.split()[0] == str(
                    natoms
                ):  # ignore comment and read any non-empty line
                    # print(line.split())
                    # print(line)
                    # print(line.split())
                    atype.append(line_1.split()[1])
                    mass.append(line_1.split()[7])
                    charge.append(line_1.split()[6])
                    natoms = natoms + 1
            break

    for i in range(natoms - 1):
        site = Atom()
        site.mass = float(mass[i])
        site.charge = float(charge[i])
        for index in range(len(dicts_atype)):  # change the for loop
            # Loop over dictionary to assign epsilon
            for key in dicts_atype[index]:
                # print(index,key,dicts_atype[index][key])
                if dicts_atype[index]["type"] == atype[i]:
                    sigma = dicts_atype[index]["sigma"]
                    epsilon = dicts_atype[index]["epsilon"]
                    break
                break

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

    file.seek(0)
    return pdb


def _parse_bonds(file, pdb):
    bind_1 = []
    bind_2 = []
    b_type = []
    b_K = []
    b_req = []
    for line in file:
        if "bonds" in line:
            nbonds = 0
            for line_1 in file:
                if "[" in line_1:
                    break
                elif line_1.split() and ";" not in line_1.split():
                    print(line)
                    bind_1.append(int(line_1.split()[0]) - 1)
                    bind_2.append(int(line_1.split()[1]) - 1)
                    b_type.append(int(line_1.split()[2]))
                    b_K.append(float(line_1.split()[3]))
                    b_req.append(float(line_1.split()[4]))
                    nbonds = nbonds + 1
            break

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
    file.seek(0)
    return pdb


def _parse_angles(file, pdb):
    # collect angle parameters in a list
    aind_1 = []
    aind_2 = []
    aind_3 = []
    a_type = []
    a_K = []
    a_thetaeq = []

    for line in file:
        if "angles" in line:
            nang = 0
            for line_1 in file:
                if "[" in line_1:
                    break
                elif line_1.split() and ";" not in line_1.split():
                    print(line)
                    aind_1.append(int(line_1.split()[0]) - 1)
                    aind_2.append(int(line_1.split()[1]) - 1)
                    aind_3.append(int(line_1.split()[2]) - 1)
                    a_type.append(int(line_1.split()[3]))
                    a_K.append(float(line_1.split()[4]))
                    a_thetaeq.append(float(line_1.split()[5]))
                    nang = nang + 1
            break

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
    file.seek(0)
    return pdb


def _parse_RB(pdb, line):
    dind_1, dind_2, dind_3, dind_4 = (
        int(line.split()[0]) - 1,
        int(line.split()[1]) - 1,
        int(line.split()[2]) - 1,
        int(line.split()[3]) - 1,
    )
    c0, c1, c2, c3, c4, c5 = (
        float(line.split()[5]),
        float(line.split()[6]),
        float(line.split()[7]),
        float(line.split()[8]),
        float(line.split()[9]),
        float(line.split()[9]),
    )

    site_1 = pdb.sites[dind_1]
    site_2 = pdb.sites[dind_2]
    site_3 = pdb.sites[dind_3]
    site_4 = pdb.sites[dind_4]
    diehdrals = Dihedral(connection_members=(site_1, site_2, site_3, site_4))
    diehdrals.dihedral_type = DihedralType(
        name="RyckaertBellemansTorsionPotential",
        expression="c0 * cos(phi)**0 + c1 * cos(phi)**1 + c2 * cos(phi)**2 + c3 * cos(phi)**3 + c4 * cos(phi)**4 + c5 * cos(phi)**5",
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
            "c0": c0 * u.Unit("kJ / (deg**2)"),
            "c1": c1 * u.Unit("kJ / (deg**2)"),
            "c2": c2 * u.Unit("kJ / (deg**2)"),
            "c3": c3 * u.Unit("kJ / (deg**2)"),
            "c4": c4 * u.Unit("kJ / (deg**2)"),
            "c5": c5 * u.Unit("kJ / (deg**2)"),
        },
    )
    pdb.add_connection(diehdrals)


def _parse_improper_harmonic(pdb,line):
    dind_1, dind_2, dind_3, dind_4 = (
    int(line.split()[0]) - 1,
    int(line.split()[1]) - 1,
    int(line.split()[2]) - 1,
    int(line.split()[3]) - 1,
    )
    d_phi,d_K,d_n = (
    float(line.split()[5]),
    float(line.split()[6]),
    float(line.split()[7]),
    )
    site_1 = pdb.sites[dind_1]
    site_2 = pdb.sites[dind_2]
    site_3 = pdb.sites[dind_3]
    site_4 = pdb.sites[dind_4]
    impropers = Improper(connection_members=(site_1, site_2, site_3, site_4))
    K = d_K
    phi_eq = d_phi
    n = d_n
    impropers.improper_type = ImproperType(
    name="PeriodicImproperPotential",
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
        "phi_eq": phi_eq * u.deg,
        "n": n * u.dimensionless,
    },
    )
    pdb.add_connection(impropers)


def _parse_harmonic(pdb, line):
    dind_1, dind_2, dind_3, dind_4 = (
        int(line.split()[0]) - 1,
        int(line.split()[1]) - 1,
        int(line.split()[2]) - 1,
        int(line.split()[3]) - 1,
    )
    d_K, d_phi, d_n = (
        float(line.split()[5]),
        float(line.split()[6]),
        float(line.split()[7]),
    )

    site_1 = pdb.sites[dind_1]
    site_2 = pdb.sites[dind_2]
    site_3 = pdb.sites[dind_3]
    site_4 = pdb.sites[dind_4]
    diehdrals = Dihedral(connection_members=(site_1, site_2, site_3, site_4))
    K = d_K
    phi_eq = d_phi
    n = d_n
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
            "phi_eq": phi_eq * u.deg,
            "n": n * u.dimensionless,
        },
    )
    pdb.add_connection(diehdrals)


def _parse_torsions(file, pdb):
    
    #This function parses both proper (dihedrals) and improper torsions
    # dry run to get the dihedral type
    for i, line in enumerate(file):
        if "dihedrals" in line:
            for i, line_1 in enumerate(file):
                if line_1.split() and ";" not in line_1.split():
                    d_type = int(line_1.split()[4])
                    break
            break
    
    file.seek(0)

    dicts_dtype = {
        1: _parse_harmonic,
        3: _parse_RB,
        4: _parse_improper_harmonic
    }  # This one should come from the GROMACS page
    # https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html
    parser = dicts_dtype[d_type]

    for line in file:
        if "dihedrals" in line:
            for line_1 in file:
                if "[" in line_1:
                    break
                elif line_1.split() and ";" not in line_1.split():
                    print(line_1)
                    parser(pdb, line_1)
            break
        
    
    # dry run to get the dihedral type
    # for i, line in enumerate(file):
        # if "dihedrals" in line:
            # for i, line_1 in enumerate(file):
                # if line_1.split() and ";" not in line_1.split():
                    # d_type = int(line_1.split()[4])
                    # break
        # break
    
    #file.seek(0)

    d_type=3
    parser = dicts_dtype[d_type]
    # #return to previous line command not working (will fix it later) so using file.seek(0)#
    for line in file:
        if  line.split() and ";" not in line.split():
            parser(pdb, line)

   

    return pdb
    

    # https://github.com/mosdef-hub/gmso/tree/main/gmso/lib/jsons
