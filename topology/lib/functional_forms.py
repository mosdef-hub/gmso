from abc import ABC, abstractmethod, abstractclassmethod

import sympy


class NonbondedPotential(ABC):
    def __init__(self, symbols=None, potential=None):
        super().__init__()
        self._symbols = symbols
        self._potential = potential

    @property
    def symbols(self):
        return self._symbols

    @symbols.setter
    def symbols(self, symbols):
        self._symbols = symbols

    @property
    def potential(self):
        return self._potential

    @potential.setter
    def potential(self, potential):
        self._potential = potential


class LennardJonesPotential(NonbondedPotential):

    def __init__(self,
                 symbols=sympy.symbols('epsilon sigma r'),
                 potential='4*epsilon*((sigma/r)**12 - (sigma/r)**6)'
                 ):
        super().__init__()

        self._symbols = symbols
        self._potential = potential
