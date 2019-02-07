import sympy

class AtomType(object):
    """An atom type."""
    def __init__(self, name="AtomType", charge=0.0, 
            nb_function='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
            parameters={'sigma':1, 'epsilon':100}):

        self._name = name
        self._charge = charge
        self._parameters = parameters
        if nb_function is None:
            self._nb_function = None
        elif isinstance(nb_function, str):
            self._nb_function = sympy.sympify(nb_function)
        elif isinstance(nb_function, sympy.Expr):
            self._nb_function = nb_function
        else:
            raise ValueError("Please enter a string, sympy expression, or None")

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

    @parameters.setter
    def parameters(self, val):
        self._parameters = val

    @property 
    def nb_function(self):
        return self._nb_function

    @nb_function.setter
    def nb_function(self, val):
        if val is None:
            self._nb_function = None
        elif isinstance(val, str):
            self._nb_function = sympy.sympify(val)
        elif isinstance(val, sympy.Expr):
            self._nb_function = val
        else:
            raise ValueError("Please enter a string or sympy expression")

    def __eq__(self, other):
        return (self.name == other.name)

    def __repr__(self):
        desc = "<AtomType {}, id {}>".format(self._name, id(self)) 
        return desc
