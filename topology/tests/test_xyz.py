from topology.formats.xyz import read_xyz


def test_read_xyz():
    top = read_xyz('ethane.xyz')

    assert top.n_sites == 8
