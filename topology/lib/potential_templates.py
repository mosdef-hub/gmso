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

    def __hash__(self):
        return hash(
            tuple(
                (
                    self.name,
                    self.expression,
                    tuple(self.independent_variables),
                )
            )
        )


class LennardJonesPotential(PotentialTemplate):
    def __init__(self,
                 name='LennardJonesPotential',
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 independent_variables={'r'}):

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
                 independent_variables={'r'}):

        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )

class HarmonicBondPotential(PotentialTemplate):
    def __init__(self,
                 name='HarmonicBondPotential',
                 expression='0.5 * k * (r-r_eq)**2',
                 independent_variables={'r'}):

        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )

class HarmonicAnglePotential(PotentialTemplate):
    def __init__(self,
                 name='HarmonicAnglePotential',
                 expression='0.5 * k * (theta-theta_eq)**2',
                 independent_variables={'theta'}):

        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=True,
        )

