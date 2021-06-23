"""Custom exception module for GMSO."""
TYPE_ERROR_STRING = "Error: Expected {0} to be of type {1}, found {2}."


class NotYetImplementedWarning(Warning):
    """Warning for behavior that is planned but not currently supported."""


class GMSOError(Exception):
    """Base class for all non-trivial errors raised by `gmso`."""


class ForceFieldError(Exception):
    """Base class for forcefield related errors."""


class ForceFieldParseError(Exception):
    """Base class for forcefield parsing errors."""


class EngineIncompatibilityError(GMSOError):
    """Error for engine incompatibility when writing or converting."""


class MissingAtomTypesError(ForceFieldParseError):
    """Error for missing AtomTypes when creating a ForceField from an XML file."""


class MissingPotentialError(ForceFieldError):
    """Error for missing Potential when searching for Potentials in a ForceField."""
