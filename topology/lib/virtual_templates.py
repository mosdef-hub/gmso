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

class ThreeSiteFDConstruction(VirtualSiteTemplate):
    def __init__(self,
                 name='ThreeSiteFDConstruction',
                 # This is not right
                 expression='1 - a - b)',
                 independent_variables=None):
        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )

class ThreeSiteFADConstruction(VirtualSiteTemplate):
    def __init__(self,
                 name='ThreeSiteFDConstruction',
                 # This is not right
                 expression='1 - a - b)',
                 independent_variables=None):
        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )

class ThreeSiteOutConstruction(VirtualSiteTemplate):
    def __init__(self,
                 name='ThreeSiteOutConstruction',
                 # This is not right
                 expression='1 - a - b)',
                 independent_variables=None):
        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )

class FourSiteConstruction(VirtualSiteTemplate):
    def __init__(self,
                 name='FourSiteConstruction',
                 # This is not right
                 expression='1 - a - b)',
                 independent_variables=None):
        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )
