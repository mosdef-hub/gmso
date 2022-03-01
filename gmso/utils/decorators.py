"""Various decorators for GMSO."""
from functools import update_wrapper, wraps

from gmso.abc import GMSOJSONHandler
from gmso.utils.io import import_


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


class register_pydantic_json(object):
    """Provides a way to register json encoders for a non-JSON serializable class."""

    def __init__(self, method="json"):
        self.method = method

    def __call__(self, cls):
        """Register this class's json encoder method to GMSOJSONHandler."""
        json_method = getattr(cls, self.method)
        GMSOJSONHandler.register(cls, json_method)
        return cls


class requires_package:
    """Decorator to use when a function requires a package to function."""

    def __init__(self, pkgname):
        self.pkgname = pkgname

    def __call__(self, func):
        """Return a wrapped function by checking if the package is importable."""

        def wrapped(*args, **kwargs):
            import_(self.pkgname)
            return func(*args, **kwargs)

        update_wrapper(wrapped, func)
        return wrapped
