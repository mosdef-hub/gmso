# ruff: noqa: F401
"""GMSO: General Molecular Simulation Object."""

from .core.angle import Angle
from .core.angle_type import AngleType
from .core.atom import Atom
from .core.atom_type import AtomType
from .core.bond import Bond
from .core.bond_type import BondType
from .core.box import Box
from .core.dihedral import Dihedral
from .core.dihedral_type import DihedralType
from .core.element import Element
from .core.forcefield import ForceField
from .core.improper import Improper
from .core.improper_type import ImproperType
from .core.pairpotential_type import PairPotentialType
from .core.topology import Topology
from .core.virtual_site import VirtualSite
from .core.virtual_type import VirtualType

__version__ = "0.14.0"

import logging
import sys
from logging.handlers import RotatingFileHandler


class DeduplicationFilter(logging.Filter):
    """A logging filter that suppresses duplicate messages."""

    def __init__(self):
        super().__init__()
        self.logged_messages = set()

    def filter(self, record):
        log_entry = (record.name, record.levelno, record.msg)
        if log_entry not in self.logged_messages:
            self.logged_messages.add(log_entry)
            return True
        return False


class GMSOLogger:
    def __init__(self):
        self.library_logger = logging.getLogger("gmso")
        self.library_logger.setLevel(logging.INFO)

        # Create handlers
        self.console_handler = logging.StreamHandler(sys.stdout)
        self.console_handler.setLevel(logging.INFO)

        self.file_handler = RotatingFileHandler(
            "gmso.log", maxBytes=10**6, backupCount=3
        )
        self.file_handler.setLevel(logging.DEBUG)

        # Create a formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )

        # Add formatter to handlers
        self.console_handler.setFormatter(formatter)
        self.file_handler.setFormatter(formatter)

        # Initialize and add the deduplication filter
        self.dedup_filter = DeduplicationFilter()
        self.console_handler.addFilter(self.dedup_filter)
        self.file_handler.addFilter(self.dedup_filter)

        # Clear any previous handlers to avoid duplicates in Jupyter
        self._clear_handlers()

        # Add handlers to the library logger
        self.library_logger.addHandler(self.console_handler)
        # self.library_logger.addHandler(self.file_handler)

    def _clear_handlers(self):
        handlers = self.library_logger.handlers[:]
        for handler in handlers:
            self.library_logger.removeHandler(handler)


# Example usage in __init__.py
# gmso_logger = GMSOLogger()
