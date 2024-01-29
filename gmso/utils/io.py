"""General module for IO operations in GMSO."""

import importlib
import inspect
import os
import sys
import textwrap
from tempfile import TemporaryFile
from unittest import SkipTest

from pkg_resources import resource_filename

MESSAGES = dict()
MESSAGES[
    "matplotlib.pyplot"
] = """
The code at {filename}:{line_number} requires the "matplotlib" package
matplotlib can be installed using:
# conda install -c conda-forge matplotlib
or
# pip install matplotlib
"""

MESSAGES[
    "matplotlib"
] = """
The code at {filename}:{line_number} requires the "matplotlib" package
matplotlib can be installed using:f
# conda install -c conda-forge matplotlib
or
# pip install matplotlib
"""


class DelayImportError(ImportError, SkipTest):
    """Delay an import error for testing utilities."""

    pass


def get_fn(filename):
    """Get the full path to one of the reference files shipped for utils.

    In the source distribution, these files are in ``gmso/utils/files``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    filename : str
        Name of the file to load (with respect to the files/ folder).

    Returns
    -------
    fn : str
        Full path to filename
    """
    fn = resource_filename("gmso", os.path.join("utils", "files", filename))
    if not os.path.exists(fn):
        raise IOError("Sorry! {} does not exists.".format(fn))
    return fn


def import_(module):
    """Import a module.

    Attempt to import a module and if it isn't installed, print a message to to STDERR
    instead of erroring out. This code is copied from foyer in foyer/utils/io.py.

    Parameters
    ----------
    module : str
        The module you'd like to import, as a string

    Returns
    -------
    module : {module, object}
        The module object

    Examples
    --------
    >>> # the following two lines are equivalent. the difference is that the
    >>> # second will check for an ImportError and print you a very nice
    >>> # user-facing message about what's wrong (where you can install the
    >>> # module from, etc) if the import fails
    >>> import tables
    >>> tables = import_('tables')
    """
    try:
        return importlib.import_module(module)
    except ImportError as e:
        try:
            message = MESSAGES[module]
        except KeyError:
            message = (
                "The code at {filename}:{line_number} requires the "
                + module
                + " package"
            )
            e = ImportError("No module named %s" % module)

        (
            frame,
            filename,
            line_number,
            function_name,
            lines,
            index,
        ) = inspect.getouterframes(inspect.currentframe())[1]

        m = message.format(
            filename=os.path.basename(filename), line_number=line_number
        )
        m = textwrap.dedent(m)

        bar = (
            "\033[91m"
            + "#" * max(len(line) for line in m.split(os.linesep))
            + "\033[0m"
        )

        print("", file=sys.stderr)
        print(bar, file=sys.stderr)
        print(m, file=sys.stderr)
        print(bar, file=sys.stderr)
        raise DelayImportError(m)


try:
    import mbuild

    has_mbuild = True
    del mbuild
except ImportError:
    has_mbuild = False

try:
    import gsd

    has_gsd = True
    del gsd
except ImportError:
    has_gsd = False

try:
    import hoomd

    has_hoomd = True
    del hoomd
except ImportError:
    has_hoomd = False

try:
    import parmed

    has_parmed = True
    del parmed
except ImportError:
    has_parmed = False

try:
    import mdtraj

    has_mdtraj = True
    del mdtraj
except ImportError:
    has_mdtraj = False

try:
    import openmm

    has_openmm = True
    del openmm
except ImportError:
    has_openmm = False

try:
    from openmm import unit

    has_openmm_unit = True
    del unit
except ImportError:
    has_openmm_unit = False

try:
    import ipywidgets

    has_ipywidgets = True
    del ipywidgets
except ImportError:
    has_ipywidgets = False

try:
    import matplotlib

    has_matplotlib = True
    del matplotlib
except ImportError:
    has_matplotlib = False

try:
    import pandas

    has_pandas = True
    del pandas
except ImportError:
    has_pandas = False

try:
    import mosdef_cassandra

    has_cassandra = True
    del mosdef_cassandra
except ImportError:
    has_cassandra = False


def run_from_ipython():
    """Verify that the code is running in an ipython kernel."""
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
