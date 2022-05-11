"""Various decorators for GMSO."""
from gmso.abc import GMSOJSONHandler


class register_pydantic_json(object):
    """Provides a way to register json encoders for a non-JSON serializable class."""

    def __init__(self, method="json"):
        self.method = method

    def __call__(self, cls):
        """Register this class's json encoder method to GMSOJSONHandler."""
        json_method = getattr(cls, self.method)
        GMSOJSONHandler.register(cls, json_method)
        return cls
