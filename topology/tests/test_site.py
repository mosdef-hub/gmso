import numpy as np
import unyt as u

from topology.core.site import Site


def test_new_site():
    site = Site(name='site')
    assert site.name == 'site'

def test_dtype():
    site = Site(name='site', position=np.zeros(3))
    assert site.position.dtype == float
    assert isinstance(site.position, u.unyt_array)
    assert isinstance(site.position, np.ndarray)
