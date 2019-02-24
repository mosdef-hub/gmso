from abc import ABC, abstractmethod, abstractclassmethod

import sympy


class NonbondedPotential(ABC):

    def __init__(self, symbols, potential):
        self._symbols = symbols
        self._potential = potential
        super().__init__()

    @property
    @abstractmethod
    def symbols(self):
        return self._symbols

    @property
    @abstractmethod
    def potential(self):
        return self._potential


class LennardJonesPotential(NonbondedPotential):

    def __init__(self):
                 symbols=sympy.symbols('epsilon sigma r'),
                 potential='4*epsilon*((sigma/r)**12 - (sigma/r)**6)'

    def set_potential(self):
        epsilon, sigma, r = self.symbols
        self._potential = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)
