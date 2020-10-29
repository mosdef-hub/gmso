![](https://anaconda.org/conda-forge/gmso/badges/license.svg)
[![](https://anaconda.org/conda-forge/gmso/badges/version.svg)](https://anaconda.org/conda-forge/gmso)
[![Build Status](https://dev.azure.com/mosdef/mosdef/_apis/build/status/mosdef-hub.gmso?branchName=master)](https://dev.azure.com/mosdef/mosdef/_build/latest?definitionId=9&branchName=master)
[![codecov](https://codecov.io/gh/mosdef-hub/gmso/branch/master/graph/badge.svg?token=rqPGwmXDzu)](undefined)

## GMSO: General Molecular Simulation Object

Flexible storage of chemical topology for molecular simulation.
With few lines of `GMSO` code, together with [`mBuild`](https://mbuild.mosdef.org) and [`foyer`](https://foyer.mosdef.org), users can rapidly prototype arbitrary parameterized chemical systems and generate data files for a wide variety of simulation engines.

```python
import foyer
from mbuild.lib.molecules import Ethane
from gmso.external.convert_parmed import from_parmed
from gmso.formats.lammpsdata import write_lammpsdata
# Start with a mBuild compound
mb_ethane = Ethane()
oplsaa = foyer.Forcefield(name='oplsaa')
# atomtype the system with foyer, and convert the resulting structure to a topology
typed_ethane = from_parmed(oplsaa.apply(mb_ethane))
typed_ethane.name = 'ethane'
# Write out lammps datafile
write_lammpsdata(typed_ethane, filename='ethane.lammps', atom_style='full')
```


Introduction
------------

`GMSO` is designed to be a general and flexible representation of chemical topolgies for molecular simulation.
With an emphasis on assuming as little as possible about the chemical system, model, or engine, `GMSO` can enable support for a variety of systems.
`GMSO` is a part of the [MoSDeF (Molecular Simulation and Design Framework)](https://mosdef.org) ecosystem, and is intended to be the backend replacement for the [`foyer` package](https://foyer.mosdef.org).


Goals and Features
------------------

`GMSO`'s goal is to provide a flexible backend framework to store topological information of a chemical system in a reproducible fashion.
**Topology** in this case is defined as the information needed to initialize a molecular simulation.
Depending on the type of simulation performed, this ranges from:
* particle positions
* particle connectivity
* box information
* forcefield data
    - functional forms defined as [`sympy` expressions](https://www.sympy.org)
    - parameters with defined units
    - partial charges
    - tabulated data
    - etc.
* Other optional data
    - particle mass
    - elemental data
    - etc.

With these driving goals for `GMSO`, the following features are enabled:
1. __Supporting a variety of models__ in the molecular simulation/computational
  chemistry community_:
  No assumptions are made about an interaction site
  representing an atom or bead, instead these can be atomistic,
  united-atom/coarse-grained, polarizable, and other models!

1. __Greater flexibility for exotic potentials__: The [`AtomType`](./gmso/core/atom_type.py) (and [analogue
  classes for intramolecular interactions](./gmso/core)) uses [`sympy`](https://www.sympy.org) to store any
  potential that can be represented by a mathematical expression.

1. __Adaptable for new engines__: by not being designed for
  compatibility with any particular molecular simulation engine or ecosystem,
  it becomes more tractable for developers in the community to add glue for
  engines that are not currently supported.

1. __Compatibility with existing community tools__: No single molecular simulation
  tool will ever be a silver bullet, so ``GMSO`` includes functions to convert
  between various file formats and libraries. These can be used in their own right to convert between objects in-memory
  and also to support conversion to file formats not natively supported at
  any given time. Currently supported conversions include:
    * [`ParmEd`](./gmso/external/convert_parmed.py)
    * [`OpenMM`](./gmso/external/convert_openmm.py)
    * [`mBuild`](./gmso/external/convert_mbuild.py)
    * more in the future!

1. __Native support for reading and writing many common file formats__: We natively have support for:
    * [`XYZ`](./gmso/formats/xyz.py)
    * [`GRO`](./gmso/formats/gro.py)
    * [`TOP`](gmso/formats/top.py)
    * [`LAMMPSDATA`](gmso/formats/lammpsdata.py)
    * indirect support, through other libraries, for many more!


Installation
------------
For full, detailed instructions, refer to the [documentation for installation](https://gmso.mosdef.org/en/latest/installation.html)

### `conda` installation quickstart
`GMSO` is available on `conda` and can be installed as:
```bash
conda install -c conda-forge gmso
```

### `pip` installation quickstart
_Note: `GMSO` is not on `pypi` currently, but its dependencies are._

```bash
git clone  https://github.com/mosdef-hub/gmso.git
cd gmso
pip install -r requirements.txt
pip install -e .
```

`pip` quickstart will install `GMSO` in [`editable` mode](https://pip.pypa.io/en/stable/reference/pip_install/#editable-installs), which means that as you edit the source code of `GMSO` those edits will be reflected in your installation.


Documentation
-------------

The full documentation can be found at [gmso.mosdef.org](https://gmso.mosdef.org).
