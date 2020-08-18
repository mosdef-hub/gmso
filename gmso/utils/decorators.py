from functools import wraps


def confirm_dict_existence(setter_function):
    """This decorator confirms that any core type
     member is in the topology's set (if it is used to
    wrap setters of the core type member class)
    """
    @wraps(setter_function)
    def setter_with_dict_removal(self, *args, **kwargs):
        if self.topology:
            self.topology._set_refs[self._set_ref].pop(self, None)
            setter_function(self, *args, **kwargs)
            self.topology._set_refs[self._set_ref][self] = (self)
            self.topology._reindex_connection_types(self._set_ref)
        else:
            setter_function(self, *args, **kwargs)
    return setter_with_dict_removal
