Data Structures in GMSO
---------------------------
Following data structures are available within GMSO.

Core Classes
============
.. autosummary::
    :nosignatures:

    gmso.Topology
    gmso.Atom
    gmso.Bond
    gmso.Angle
    gmso.Dihedral
    gmso.Improper


Topology
********
    .. autoclass:: gmso.Topology
        :members: add_site, add_connection, update_topology

Atom
****
    .. autoclass:: gmso.Atom
        :members:
        :exclude-members: model_config, model_fields

Bond
****
    .. autoclass:: gmso.Bond
        :members:
        :exclude-members: model_config, model_fields

Angle
*****
    .. autoclass:: gmso.Angle
        :members:
        :exclude-members: model_config, model_fields

Dihedral
********
    .. autoclass:: gmso.Dihedral
        :members:
        :exclude-members: model_config, model_fields

Improper
********
    .. autoclass:: gmso.Improper
        :members:
        :exclude-members: model_config, model_fields


Potential Classes
=================
.. autosummary::
    :nosignatures:

    gmso.AtomType
    gmso.BondType
    gmso.AngleType
    gmso.DihedralType
    gmso.ImproperType
    gmso.PairPotentialType

AtomType
********
    .. autoclass:: gmso.AtomType
        :members:
        :exclude-members: model_config, model_fields

BondType
********
    .. autoclass:: gmso.BondType
        :members:
        :exclude-members: model_config, model_fields

AngleType
**********
    .. autoclass:: gmso.AngleType
        :members:
        :exclude-members: model_config, model_fields

DihedralType
************
    .. autoclass:: gmso.DihedralType
        :members:
        :exclude-members: model_config, model_fields

ImproperType
************
    .. autoclass:: gmso.ImproperType
        :members:
        :exclude-members: model_config, model_fields

PairPotentialType
*****************
    .. autoclass:: gmso.PairPotentialType
        :members:
        :exclude-members: model_config, model_fields

ForceField
==========
    .. autoclass:: gmso.ForceField
        :members:
        :exclude-members: model_config, model_fields
