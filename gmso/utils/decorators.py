"""Various decorators for GMSO."""
from functools import wraps


def confirm_dict_existence(setter_function):
    """Confirm that any core type member is in the topology's set.

    If it is used to wrap setters of the core type member class
    """

    @wraps(setter_function)
    def setter_with_dict_removal(self, *args, **kwargs):
        if self.topology:
            self.topology._set_refs[self.set_ref].pop(self, None)
            setter_function(self, *args, **kwargs)
            self.topology._set_refs[self.set_ref][self] = self
            self.topology._reindex_connection_types(self.set_ref)
        else:
            setter_function(self, *args, **kwargs)

    return setter_with_dict_removal
