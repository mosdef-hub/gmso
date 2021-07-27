"""Registry utilities to handle formats for gmso Topology."""


class UnsupportedFileFormatError(Exception):
    """Exception to be raised whenever the file loading or saving is not supported."""


class Registry:
    """A registry to incorporate a callable with a file extension."""

    def __init__(self):
        self.handlers = {}

    def _assert_can_process(self, extension):
        if extension not in self.handlers:
            raise UnsupportedFileFormatError(
                f"Extension {extension} cannot be processed as no utility "
                f"is defined in the current API to handle {extension} files."
            )

    def get_callable(self, extension):
        """Get the callable associated with extension."""
        self._assert_can_process(extension)
        return self.handlers[extension]


SaversRegistry = Registry()
LoadersRegistry = Registry()


class saves_as:
    """Decorator to aid saving."""

    def __init__(self, extension):
        self.extension = extension

    def __call__(self, method):
        """Register the method as saver for an extension."""
        SaversRegistry.handlers[self.extension] = method


class loads_as:
    """Decorator to aid loading."""

    def __init__(self, extension):
        self.extension = extension

    def __call__(self, method):
        """Register the method as loader for an extension."""
        LoadersRegistry.handlers[self.extension] = method
