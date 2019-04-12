class NotYetImplementedWarning(Warning):
    """Warning for behavior that is planned but not currently supported."""

class TopologyError(Exception):
    """Base class for all non-trivial errors raised by `topology`."""

class RedundancyError(Exception):
    """Error indicating there would be undesired redundancy."""
