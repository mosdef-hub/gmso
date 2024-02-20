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


class MixedClassAndTypesError(ForceFieldParseError):
    """Error for missing AtomTypes when creating a ForceField from an XML file."""


class MissingPotentialError(ForceFieldError):
    """Error for missing Potential when searching for Potentials in a ForceField."""


class ParameterError(GMSOError):
    """Errors related to parameters."""

    def __init__(self, param, expected):
        self.param = param
        self.params = expected


class UnknownParameterError(ParameterError):
    """Errors to be raised when a parameter is unknown."""

    def __str__(self):
        """Error message."""
        err = f"Parameter {self.param} is not one of the expected parameters {self.params}"
        return err


class MissingParameterError(ParameterError):
    """Error to be raised when a parameter is missing."""

    def __str__(self):
        """Error message."""
        err = f"Parameter '{self.param}' missing from the provided parameters {self.params}"
        return err
