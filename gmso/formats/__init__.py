"""Readers and writers for various file formats."""

from gmso.utils.io import has_ipywidgets

from .formats_registry import LoadersRegistry, SaversRegistry
from .gro import read_gro, write_gro
from .gsd import write_gsd
from .json import write_json
from .lammpsdata import write_lammpsdata
from .mcf import write_mcf
from .mol2 import read_mol2
from .top import write_top
from .xyz import read_xyz, write_xyz

if has_ipywidgets:
    from .networkx import (
        interactive_networkx_angles,
        interactive_networkx_atomtypes,
        interactive_networkx_bonds,
        interactive_networkx_dihedrals,
    )
