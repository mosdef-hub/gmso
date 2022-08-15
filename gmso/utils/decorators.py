"""Various decorators for GMSO."""
import functools

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

def deprecate_kwargs(deprecated_kwargs=None):
    if deprecated_kwargs is None:
        deprecated_kwargs = set()

    def decorate_deprecate_kwargs(func):
        @functools.wraps(func)
        def wrapper(self_or_cls, *args, **kwargs):
            _deprecate_kwargs(kwargs, deprecated_kwargs)
            return func(self_or_cls, *args, **kwargs)

        return wrapper

    return decorate_deprecate_kwargs


def _deprecate_kwargs(kwargs, deprecated_kwargs):
    added_args = []
    for kwarg in kwargs:
        if kwarg in deprecated_kwargs:
            added_args.append(kwarg)
    if len(added_args) > 1:
        message = (
            "Keyword arguments `{dep_args}` are deprecated and will be removed in the "
            "next minor release of the package. Please update your code accordingly"
        )
    else:
        message = (
            "Keyword argument `{dep_args}` is deprecated and will be removed in the "
            "next minor release of the package. Please update your code accordingly"
        )
    if added_args:
        warnings.warn(
            message.format(dep_args=", ".join(added_args)),
            DeprecationWarning,
            3,
        )
