"""Various decorators for GMSO."""

import functools
from warnings import warn
import logging

logger = logging.getLogger(__name__)


def deprecate_kwargs(deprecated_kwargs=None):
    """Decorate functions with deprecated/deprecating kwargs."""
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
    added_params = []
    deprecated_args = [kwarg[0] for kwarg in deprecated_kwargs]
    deprecated_params = [kwarg[1] for kwarg in deprecated_kwargs]
    for kwarg in kwargs:
        if kwarg in deprecated_args and kwargs[kwarg] in deprecated_params:
            added_args.append(kwarg[0])
            added_params.append(kwarg[1])
    if len(added_args) > 1:
        message = (
            "Keyword arguments `{dep_args}={dep_params}` are deprecated and will be removed in the "
            "next minor release of the package. Please update your code accordingly"
        )
    else:
        message = (
            "Keyword argument `{dep_args}={dep_params}` is deprecated and will be removed in the "
            "next minor release of the package. Please update your code accordingly"
        )
    if added_args:
        warn(
            message.format(
                dep_args=", ".join(added_args),
                dep_params=", ".join(added_params),
            ),
            DeprecationWarning,
            3,
        )


def mark_WIP(message=""):
    """Decorate functions with WIP marking."""

    def _function_wrapper(function):
        @functools.wraps(function)
        def _inner(*args, **kwargs):
            logger.info(
                "Call to function {} is WIP.".format(function.__name__),
                category=UserWarning,
                stacklevel=2,
            )
            return function(*args, **kwargs)

        return _inner

    return _function_wrapper


def deprecate_function(msg, klass=PendingDeprecationWarning):
    """Raise a warning that a given function will be deprecated soon."""

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            logger.info(msg, klass, stacklevel=2)
            return func(*args, **kwargs)

        return wrapper

    return decorator
