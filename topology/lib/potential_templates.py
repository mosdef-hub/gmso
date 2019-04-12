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
