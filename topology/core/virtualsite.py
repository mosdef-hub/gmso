import warnings
from topology.core.site import Site
from topology.core.potential import Potential


class VirtualSite(object):
    """A virtual site construction

    Partners
    --------
    virtual_site: topology.site
        An instance of topology.site that is the virtual site
    virtual_site_members: list of topology.Site
        Should be length 2 or 3
    virtual_site_type: topology.Potential
        An instance of topology.Potential that describes
        the virtual site function and parameters of this interaction
    name: name of the virtual site
    """

    def __init__(self, virtual_site, virtual_site_members=[], virtual_site_type=None, name="VirtualSite"):
        if virtual_site_members is None:
            virtual_site_members = tuple()
        
        # TODO: validate virtual site members
        #self._virtual_site_members = _validate_virtual_site_members(virtual_site_members)
        self._virtual_site = virtual_site
        self._virtual_site_members = virtual_site_members
        self._virtual_site_type = _validate_vs_type(virtual_site_type)
        self._name = _validate_name(name)
        # TODO: update members
        #self._update_members()

    @property
    def virtual_site(self):
        return self._virtual_site

    @property
    def virtual_site_members(self):
        return self._virtual_site_members

    @property
    def virtual_site_type(self):
        return self._virtual_site_type

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, vsname):
        self._name = _validate_name(vsname)

    def __repr__(self):
        descr = '<{} Virtual Site, id {}'.format(
                self.virtual_site_type, id(self))
        if self.name:
            descr += ', name {}'.format(self.name)
        descr += '>'

        return descr

def _validate_vs_type(c_type):
    if c_type is None:
        warnings.warn("Non-parametrized virtual site detected")
    elif not isinstance(c_type, Potential):
        raise TopologyError("Supplied non-Potential {}".format(c_type))
    return c_type

def _validate_name(conname):
    if not isinstance(conname, str):
        raise TopologyError("Supplied name {} is not a string".format(conname))
    return conname
