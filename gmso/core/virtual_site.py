
class VirtualSite(Site):
    """A generalized virtual site class in GMSO.

    Virtual sites are massless particles that represet represent off-atom charge sites, lone pairs, or other non-physical sites.
    
    Attributes
    ----------
    charge : float
        The charge of the virtual site in elementary charge units.
    parent_atoms : List[Site]
        The real constituent atoms that define the virtual site's position.
    """
