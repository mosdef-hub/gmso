from topology.core.atom_type import AtomType


def test_new_atom_type():
    new_type = AtomType(name='mytype')
    assert new_type.name == 'mytype'
