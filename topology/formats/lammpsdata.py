from __future__ import division

from collections import OrderedDict
from warnings import warn
import itertools as it

import numpy as np


from topology.core.topology import Topology
from topology.core.site import Site


def write_lammpsdata(topology, filename, atom_style='full',
        nbfix_in_data_file=True):
