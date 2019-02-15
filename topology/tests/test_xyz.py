import unyt as u

from topology.formats.xyz import read_xyz
from topology.utils.io import get_fn


def test_read_xyz():
    top = read_xyz(get_fn('ethane.xyz'))
    assert top.n_sites == 8
    assert top.n_connections == 0
    assert set([type(site.position) for site in top.site_list]) == {u.unyt_array}
    assert set([site.position.units for site in top.site_list]) == {u.nm}

    top = read_xyz(get_fn('cu_block.xyz'))
    assert top.n_sites == 108
    assert top.n_connections == 0
    assert set([type(site.position) for site in top.site_list]) == {u.unyt_array}
    assert set([site.position.units for site in top.site_list]) == {u.nm}
