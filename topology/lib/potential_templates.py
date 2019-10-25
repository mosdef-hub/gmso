from topology.core.potential import Potential


class PotentialTemplate(Potential):
    def __init__(self,
                 name='PotentialTemplate',
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 independent_variables={'r'},
                 template=True):

        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=template,
        )

class LennardJonesPotential(PotentialTemplate):
    def __init__(self,
                 name='LennardJonesPotential',
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 independent_variables=None):
        if independent_variables is None:
            independent_variables = {'r'}
        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )

class BuckinghamPotential(PotentialTemplate):
    def __init__(self,
                 name='BuckinghamPotential',
                 expression='a*exp(-b*r) - c*r**-6',
                 independent_variables=None):
        if independent_variables is None:
            independent_variables = {'r'}

        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )
