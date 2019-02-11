import numpy as np
import sympy

class AtomType(object):
    """An atom type."""
    def __init__(self, name="AtomType", charge=0.0, 
            nb_function='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
            parameters={'sigma':1, 'epsilon':100}):

        self._name = name
        self._charge = charge

        if isinstance(parameters, dict):
            self._parameters = parameters
        else:
            raise ValueError("Please enter dictionary for parameters")


        if nb_function is None:
            self._nb_function = None
        elif isinstance(nb_function, str):
            self._nb_function = sympy.sympify(nb_function)
        elif isinstance(nb_function, sympy.Expr):
            self._nb_function = nb_function
        else:
            raise ValueError("Please enter a string, sympy expression, "
                            "or None for nb_function")

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, val):
        self._name = val

    @property
    def charge(self):
        return self._charge

    @charge.setter
    def charge(self, val):
        self._charge = val

    @property
    def parameters(self):
        return self._parameters 

    @property 
    def nb_function(self):
        return self._nb_function

    def set_nb_function(self, function=None, parameters=None):
        # Check valid function type (string or sympy expression)
        # If func is undefined, just keep the old one
        if function is None:
            pass
        elif isinstance(function, str):
            self._nb_function = sympy.sympify(function)
        elif isinstance(function, sympy.Expr):
            self._nb_function = function
        else:
            raise ValueError("Please enter a string or sympy expression")

        # If params is undefined, keep the old one
        if parameters is None:
            symbols = sympy.symbols(set(self.parameters.keys()))
            parameters = self.parameters
        else:
            symbols = sympy.symbols(set(parameters.keys()))

        # Now verify that the parameters and nb_function have consistent symbols
        if symbols.issubset(self.nb_function.free_symbols):
            # Rebuild the parameters, eliminate unnecessary parameters
            self._parameters.update(parameters)
            self._parameters = {key: val for key, val in self._parameters.items() 
                    if key in set(str(sym) for sym in self.nb_function.free_symbols)}
        else:
            raise ValueError("NB function and parameter symbols do not agree")

    def __eq__(self, other):
        return ((self.name == other.name) & 
                (np.isclose(self.charge, other.charge)) & 
                (self.parameters == other.parameters) & 
                (self.nb_function == other.nb_function))

    def __repr__(self):
        desc = "<AtomType {}, id {}>".format(self._name, id(self)) 
        return desc
