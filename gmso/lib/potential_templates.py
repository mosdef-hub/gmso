from gmso.core.potential import Potential


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

class MiePotential(PotentialTemplate):
    def __init__(self,
                 name='MiePotential',
                 expression=('(n/(n-m)) * (n/m)**(m/(n-m)) * '
                            'epsilon * ((sigma/r)**n - (sigma/r)**m)'),
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

class HarmonicTorsionPotential(PotentialTemplate):
    def __init__(self,
                 name='HarmonicTorsionPotential',
                 expression='0.5 * k * (phi - phi_eq)**2',
                 independent_variables={'phi'}):

        super(PotentialTemplate, self).__init__(
                name=name,
                expression=expression,
                independent_variables=independent_variables,
                template=True
        )

class PeriodicTorsionPotential(PotentialTemplate):
    def __init__(self,
                 name='PeriodicTorsionPotential',
                 expression='k * (1 + cos(n * phi - phi_eq))**2',
                 independent_variables={'phi'}):

        super(PotentialTemplate, self).__init__(
                name=name,
                expression=expression,
                independent_variables=independent_variables,
                template=True
        )

class OPLSTorsionPotential(PotentialTemplate):
    def __init__(self,
                 name='OPLSTorsionPotential',
                 expression=('0.5 * k0 + '
                             '0.5 * k1 * (1 + cos(phi)) + '
                             '0.5 * k2 * (1 - cos(2*phi)) + '
                             '0.5 * k3 * (1 + cos(3*phi)) + '
                             '0.5 * k4 * (1 - cos(4*phi))'),
                 independent_variables={'phi'}):

        super(PotentialTemplate, self).__init__(
                name=name,
                expression=expression,
                independent_variables=independent_variables,
                template=True
        )

class RyckaertBellemansTorsionPotential(PotentialTemplate):
    def __init__(self,
                 name='RyckaertBellemansTorsionPotential',
                 expression=('c0 * cos(phi)**0 + c1 * cos(phi)**1 + '
                             'c2 * cos(phi)**2 + c3 * cos(phi)**3 + '
                             'c4 * cos(phi)**4 + c5 * cos(phi)**5'),
                 independent_variables={'phi'}):

        super(PotentialTemplate, self).__init__(
                name=name,
                expression=expression,
                independent_variables=independent_variables,
                template=True
        )

class HarmonicImproperPotential(PotentialTemplate):
    def __init__(self,
                 name='HarmonicImproperPotetial',
                 expression='0.5 * k * (phi - phi_eq)**2',
                 independent_variables={'phi'}):

        super(PotentialTemplate, self).__init__(
                name=name,
                expression=expression,
                independent_variables=independent_variables,
                template=True
        )

