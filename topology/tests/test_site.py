from topology.core.site import Site


def test_new_site():
    site = Site(name='site')
    assert site.name == 'site'