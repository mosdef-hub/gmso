from functools import wraps

from gmso.abc import GMSOJSONHandler

def confirm_dict_existence(setter_function):
    """This decorator confirms that any core type
     member is in the topology's set (if it is used to
    wrap setters of the core type member class)
    """
    @wraps(setter_function)
    def setter_with_dict_removal(self, *args, **kwargs):
        if self.topology:
            self.topology._set_refs[self.set_ref].pop(self, None)
            setter_function(self, *args, **kwargs)
            self.topology._set_refs[self.set_ref][self] = (self)
            self.topology._reindex_connection_types(self.set_ref)
        else:
            setter_function(self, *args, **kwargs)
    return setter_with_dict_removal


class register_pydantic_json(object):
    """Provides a way to register json encoders for a non-JSON serializable class"""
    def __init__(self, method='json'):
        self.method = method

    def __call__(self, cls):
        """Register this class's json encoder method to GMSOJSONHandler"""
        json_method = getattr(cls, self.method)
        GMSOJSONHandler.register(cls, json_method)
        return cls
