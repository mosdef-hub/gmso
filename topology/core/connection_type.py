import numpy as np
import sympy
import unyt as u

class ConnectionType(object):
    """A connection type."""
    def __init__(self, 
            potential_function='0.5 * k * (r-r_eq)**2',
            parameters={'k': 1000 * u.joule / (u.mol * u.nm**2), 
                'r_eq':1 * u.nm}):

        if isinstance(parameters, dict):
            self._parameters = parameters
        else:
            raise ValueError("Please enter dictionary for parameters")

        if potential_function is None:
            self._potential_function = None
        elif isinstance(potential_function, str):
            self._potential_function = sympy.sympify(potential_function)
        elif isinstance(potential_function, sympy.Expr):
            self._potential_function = potential_function
        else:
            raise ValueError("Please enter a string, sympy expression, "
                            "or None for potential_function")
    
    @property
    def parameters(self):
        return self._parameters

    @property
    def potential_function(self):
        return self._potential_function

    def set_potential_function(self, function=None, parameters=None):
        # Check valid function type (string or sympy expression)
        # If func is undefined, just keep the old one
        if function is None:
            pass
        elif isinstance(function, str):
            self._potential_function = sympy.sympify(function)
        elif isinstance(function, sympy.Expr):
            self._potential_function = function
        else:
            raise ValueError("Please enter a string or sympy expression")

         # If params is undefined, keep the old one
        if parameters is None:
            parameters = self.parameters
        symbols = sympy.symbols(set(parameters.keys()))

         # Now verify that the parameters and potential_function have consistent symbols
        if symbols.issubset(self.potential_function.free_symbols):
            # Rebuild the parameters, eliminate unnecessary parameters
            self._parameters.update(parameters)
            self._parameters = {key: val for key, val in self._parameters.items()
                    if key in set(str(sym) 
                        for sym in self.potential_function.free_symbols)}
        else:
            extra_syms = symbols - self.potential_function.free_symbols
            raise ValueError("Potential function and parameter"
                            " symbols do not agree,"
                            " you supplied extraneous symbols:"
                            " {}".format(extra_syms))

    def __eq__(self, other):
        return ((self.parameters == other.parameters) &
                (self.potential_function == other.potential_function))
