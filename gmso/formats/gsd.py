"""Write GSD files from GMSO topologies."""

from __future__ import division

from gmso.external.convert_hoomd import to_gsd_snapshot
from gmso.formats.formats_registry import saves_as
from gmso.utils.io import has_gsd

__all__ = ["write_gsd"]

if has_gsd:
    import gsd.hoomd


@saves_as(".gsd")
def write_gsd(
    top,
    filename,
    base_units=None,
    shift_coords=True,
    write_special_pairs=True,
):
    """Output a GSD file (HOOMD v3 default data format).

    The `GSD` binary file format is the native format of HOOMD-Blue. This file
    can be used as a starting point for a HOOMD-Blue simulation, for analysis,
    and for visualization in various tools.

    Parameters
    ----------
    top : gmso.Topology
        gmso.Topology object
    filename : str
        Path of the output file.
    shift_coords : bool, optional, default=True
        Shift coordinates from (0, L) to (-L/2, L/2) if necessary.
    write_special_pairs : bool, optional, default=True
        Writes out special pair information necessary to correctly use the OPLS
        fudged 1,4 interactions in HOOMD.

    Notes
    -----
    Force field parameters are not written to the GSD file and must be included
    manually in a HOOMD input script. Work on a HOOMD plugin is underway to
    read force field parameters from a Foyer XML file.

    """
    gsd_snapshot = to_gsd_snapshot(
        top=top,
        base_units=base_units,
        shift_coords=shift_coords,
        parse_special_pairs=write_special_pairs,
    )[0]
    with gsd.hoomd.open(filename, mode="w") as gsd_file:
        gsd_file.append(gsd_snapshot)
