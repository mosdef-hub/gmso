from topology.core.potential import Potential


class VirtualSiteTemplate(Potential):
    def __init__(self,
                 name='VirtualSiteTemplate',
                 expression='1 - a',
                 independent_variables=None,
                 template=True):

        super(VirtualSiteTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=template,
        )

class TwoSiteConstruction(VirtualSiteTemplate):
    def __init__(self,
                 name='TwoSiteConstruction',
                 expression='1 - a)',
                 independent_variables=None):
        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )

class ThreeSiteConstruction(VirtualSiteTemplate):
    def __init__(self,
                 name='ThreeSiteConstruction',
                 expression='1 - a - b)',
                 independent_variables=None):
        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )
