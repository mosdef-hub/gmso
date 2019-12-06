import unyt as u

UNITS_MAP = {
    'amu': u.gram / u.mol
}

def validate(xml, schema):
    """Validate a given xml file with a refrence schema"""
    pass

def mass_to_unyt(mass):
    try:
        unit = str(mass).split(' ')[1]
        return float(str(mass)).split(' ')[0] * UNITS_MAP[unit]
    except IndexError:
        return float(mass) * UNITS_MAP['amu']
