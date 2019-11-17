from functools import wraps


def confirm_set_existence(setter_function):
    """This decorator confirms that any core type
     member is in the topology's set (if it is used to
    wrap setters of the core type member class)
    """
    @wraps(setter_function)
    def setter_function_with_set_removal(self, *args, **kwargs):
        if self._topology and (self in self._topology[self._set_ref]):
            self._topology['set_ref'].discard(self)
            setter_function(self, *args, **kwargs)
            self._topology['set_ref'].add(self)
        else:
            setter_function(self, *args, **kwargs)
    return setter_function_with_set_removal
