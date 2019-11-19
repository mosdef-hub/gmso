

TYPE_ERROR_STRING = "Error: Expected {0} to be of type {1}, found {2}."

class NotYetImplementedWarning(Warning):
    """Warning for behavior that is planned but not currently supported."""


class TopologyError(Exception):
    """Base class for all non-trivial errors raised by `topology`."""
