"""Write GSD files from GMSO topologies."""

from __future__ import division

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from gmso.external.convert_hoomd import to_gsd_snapshot
from gmso.formats.formats_registry import saves_as
from gmso.utils.io import has_gsd

if TYPE_CHECKING:
    from gmso.core.topology import Topology

__all__ = ["write_gsd"]

if has_gsd:
    import gsd.hoomd


@saves_as(".gsd")
def write_gsd(
    top: "Topology",
    filename: Union[str, Path],
    base_units: Optional[dict] = None,
    shift_coords: bool = True,
    write_special_pairs: bool = True,
) -> None:
    """Output a GSD file (HOOMD default data format).

    The GSD binary file format is the native format of HOOMD-blue.  It can be
    used as a starting point for a HOOMD-blue simulation, for trajectory
    analysis, and for visualisation in various tools.

    To write a HOOMD snapshot, see
    :func:`gmso.external.convert_hoomd.to_hoomd_snapshot` and
    :func:`gmso.external.convert_hoomd.to_gsd_snapshot`

    Parameters
    ----------
    top : gmso.Topology
        Typed topology to write.
    filename : str or pathlib.Path
        Path of the output file.
    base_units : dict, optional, default=None
        Dictionary of base units to use when writing the GSD snapshot.
        If ``None``, HOOMD-blue default units are used.
    shift_coords : bool, optional, default=True
        Shift coordinates from ``(0, L)`` to ``(-L/2, L/2)`` if necessary.
    write_special_pairs : bool, optional, default=True
        Write special pair information needed to correctly apply the OPLS
        fudged 1-4 interactions in HOOMD-blue.

    Returns
    -------
    None
        Writes the GSD snapshot to *filename* in place.

    Notes
    -----
    Force field parameters are not stored in the GSD file.
    You can use GMSO to create the HOOMD force objects. See :func:`gmso.external.convert_hoomd.to_hoomd_forcefield`

    """
    gsd_snapshot = to_gsd_snapshot(
        top=top,
        base_units=base_units,
        shift_coords=shift_coords,
        parse_special_pairs=write_special_pairs,
    )[0]
    with gsd.hoomd.open(filename, mode="w") as gsd_file:
        gsd_file.append(gsd_snapshot)
