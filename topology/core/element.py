import unyt as u
from collections import OrderedDict
class Element(object):
    """An element."""
    def __init__(self, atomic_number=None,name=None, symbol=None, mass=None):
        self.atomic_number = atomic_number
        self.name = name
        self.symbol = symbol
        self.mass = mass



def element_by_symbol(symbols):
    # In: symbol or list of symbol (str or list of str)
    # Out: Element or list of Element
    # Search by symbol of the element
    # First trial: create dict to map symbol (key) to element (value)
    symbol_dict = {element.symbol:element for element in Periodict}
    element_list = []
    if isinstance(symbols, str):
        symbols = list(symbols)
    for symbol in symbols:
        element_list.append(symbol_dict[symbol])
    return element_list

def element_by_name(names):
    # In: name or list of name (str or list of str)
    # Out: Element or list of Element
    # Search by name of the element
    name_dict = {element.name:element for element in Periodict}
    element_list = []
    if isinstance(names, str):
        names = list(names)
    for name in names:
        name = name.lower()
        element_list.append(name_dict[name])
    return element_list

def element_by_atomic_number(atomic_numbers):
    # In: atomic number or list of atomic number (int or list of int)
    # Out: Element or list of Element
    # Search by atomic number of the element
    atomic_dict = {element.atomic_number:element for element in Periodict}
    element_list = []
    if isinstance(atomic_numbers, int):
        atomic_numbers = list(atomic_numbers)
    for number in atomic_numbers:
        element_list.append(atomic_dict[number])
    return element_list

def element_by_mass(masses):
    # In: mass or list of mass (int/float or list of int/float)
    # Out: Element or list of Element
    # Search by mass
    return None



Hydrogen = 	Element(atomic_number=1, name='hydrogen', symbol='H', mass=1.0079 * u.amu)
Helium = 	Element(atomic_number=2, name='helium', symbol='He', mass=4.0026 * u.amu)
Lithium = 	Element(atomic_number=3, name='lithium', symbol='Li', mass=6.941 * u.amu)
Beryllium = 	Element(atomic_number=4, name='beryllium', symbol='Be', mass=9.0122 * u.amu)
Boron = 	Element(atomic_number=5, name='boron', symbol='B', mass=10.811 * u.amu)
Carbon = 	Element(atomic_number=6, name='carbon', symbol='C', mass=12.0107 * u.amu)
Nitrogen = 	Element(atomic_number=7, name='nitrogen', symbol='N', mass=14.0067 * u.amu)
Oxygen = 	Element(atomic_number=8, name='oxygen', symbol='O', mass=15.9994 * u.amu)
Fluorine = 	Element(atomic_number=9, name='fluorine', symbol='F', mass=18.9984 * u.amu)
Neon = 		Element(atomic_number=10, name='neon', symbol='Ne', mass=20.1797 * u.amu)
Sodium = 	Element(atomic_number=11, name='sodium', symbol='Na', mass=22.9898 * u.amu)
Magnesium = 	Element(atomic_number=12, name='magnesium', symbol='Mg', mass=24.305 * u.amu)
Aluminum = 	Element(atomic_number=13, name='aluminum', symbol='Al', mass=26.9815 * u.amu)
Silicon = 	Element(atomic_number=14, name='silicon', symbol='Si', mass=28.0855 * u.amu)
Phosphorus = 	Element(atomic_number=15, name='phosphorus', symbol='P', mass=30.9738 * u.amu)
Sulfur = 	Element(atomic_number=16, name='sulfur', symbol='S', mass=32.065 * u.amu)
Chlorine = 	Element(atomic_number=17, name='chlorine', symbol='Cl', mass=35.453 * u.amu)
Argon = 	Element(atomic_number=18, name='argon', symbol='Ar', mass=39.948 * u.amu)
Potassium = 	Element(atomic_number=19, name='potassium', symbol='K', mass=39.0983 * u.amu)
Calcium = 	Element(atomic_number=20, name='calcium', symbol='Ca', mass=40.078 * u.amu)
Scandium = 	Element(atomic_number=21, name='scandium', symbol='Sc', mass=44.9559 * u.amu)
Titanium = 	Element(atomic_number=22, name='titanium', symbol='Ti', mass=47.867 * u.amu)
Vanadium = 	Element(atomic_number=23, name='vanadium', symbol='V', mass=50.9415 * u.amu)
Chromium = 	Element(atomic_number=24, name='chromium', symbol='Cr', mass=51.9961 * u.amu)
Manganese = 	Element(atomic_number=25, name='manganese', symbol='Mn', mass=54.938 * u.amu)
Iron = 		Element(atomic_number=26, name='iron', symbol='Fe', mass=55.845 * u.amu)
Cobalt = 	Element(atomic_number=27, name='cobalt', symbol='Co', mass=58.9331 * u.amu)
Nickel = 	Element(atomic_number=28, name='nickel', symbol='Ni', mass=58.6934 * u.amu)
Copper = 	Element(atomic_number=29, name='copper', symbol='Cu', mass=63.546 * u.amu)
Zinc = 		Element(atomic_number=30, name='zinc', symbol='Zn', mass=65.409 * u.amu)
Gallium = 	Element(atomic_number=31, name='gallium', symbol='Ga', mass=69.723 * u.amu)
Germanium = 	Element(atomic_number=32, name='germanium', symbol='Ge', mass=72.64 * u.amu)
Arsenic = 	Element(atomic_number=33, name='arsenic', symbol='As', mass=74.9216 * u.amu)
Selenium = 	Element(atomic_number=34, name='selenium', symbol='Se', mass=78.96 * u.amu)
Bromine = 	Element(atomic_number=35, name='bromine', symbol='Br', mass=79.904 * u.amu)
Krypton = 	Element(atomic_number=36, name='krypton', symbol='Kr', mass=83.798 * u.amu)
Rubidium = 	Element(atomic_number=37, name='rubidium', symbol='Rb', mass=85.4678 * u.amu)
Strontium = 	Element(atomic_number=38, name='strontium', symbol='Sr', mass=87.62 * u.amu)
Yttrium = 	Element(atomic_number=39, name='yttrium', symbol='Y', mass=88.9059 * u.amu)
Zirconium = 	Element(atomic_number=40, name='zirconium', symbol='Zr', mass=91.224 * u.amu)
Niobium = 	Element(atomic_number=41, name='niobium', symbol='Nb', mass=92.9064 * u.amu)
Molybdenum = 	Element(atomic_number=42, name='molybdenum', symbol='Mo', mass=95.94 * u.amu)
Technetium = 	Element(atomic_number=43, name='technetium', symbol='Tc', mass=98.0 * u.amu)
Ruthenium = 	Element(atomic_number=44, name='ruthenium', symbol='Ru', mass=101.07 * u.amu)
Rhodium = 	Element(atomic_number=45, name='rhodium', symbol='Rh', mass=102.9055 * u.amu)
Palladium = 	Element(atomic_number=46, name='palladium', symbol='Pd', mass=106.42 * u.amu)
Silver = 	Element(atomic_number=47, name='silver', symbol='Ag', mass=107.8682 * u.amu)
Cadmium = 	Element(atomic_number=48, name='cadmium', symbol='Cd', mass=112.411 * u.amu)
Indium = 	Element(atomic_number=49, name='indium', symbol='In', mass=114.818 * u.amu)
Tin = 		Element(atomic_number=50, name='tin', symbol='Sn', mass=118.71 * u.amu)
Antimony = 	Element(atomic_number=51, name='antimony', symbol='Sb', mass=121.76 * u.amu)
Tellurium = 	Element(atomic_number=52, name='tellurium', symbol='Te', mass=127.6 * u.amu)
Iodine = 	Element(atomic_number=53, name='iodine', symbol='I', mass=126.9045 * u.amu)
Xenon = 	Element(atomic_number=54, name='xenon', symbol='Xe', mass=131.293 * u.amu)
Cesium = 	Element(atomic_number=55, name='cesium', symbol='Cs', mass=132.9055 * u.amu)
Barium = 	Element(atomic_number=56, name='barium', symbol='Ba', mass=137.327 * u.amu)
Lanthanum = 	Element(atomic_number=57, name='lanthanum', symbol='La', mass=138.9055 * u.amu)
Cerium = 	Element(atomic_number=58, name='cerium', symbol='Ce', mass=140.116 * u.amu)
Praseodymium = 	Element(atomic_number=59, name='praseodymium', symbol='Pr', mass=140.9077 * u.amu)
Neodymium = 	Element(atomic_number=60, name='neodymium', symbol='Nd', mass=144.242 * u.amu)
Promethium = 	Element(atomic_number=61, name='promethium', symbol='Pm', mass=145.0 * u.amu)
Samarium = 	Element(atomic_number=62, name='samarium', symbol='Sm', mass=150.36 * u.amu)
Europium = 	Element(atomic_number=63, name='europium', symbol='Eu', mass=151.964 * u.amu)
Gadolinium = 	Element(atomic_number=64, name='gadolinium', symbol='Gd', mass=157.25 * u.amu)
Terbium = 	Element(atomic_number=65, name='terbium', symbol='Tb', mass=158.9254 * u.amu)
Dysprosium = 	Element(atomic_number=66, name='dysprosium', symbol='Dy', mass=162.5 * u.amu)
Holmium = 	Element(atomic_number=67, name='holmium', symbol='Ho', mass=164.9303 * u.amu)
Erbium = 	Element(atomic_number=68, name='erbium', symbol='Er', mass=167.259 * u.amu)
Thulium = 	Element(atomic_number=69, name='thulium', symbol='Tm', mass=168.9342 * u.amu)
Ytterbium = 	Element(atomic_number=70, name='ytterbium', symbol='Yb', mass=173.04 * u.amu)
Lutetium = 	Element(atomic_number=71, name='lutetium', symbol='Lu', mass=174.967 * u.amu)
Hafnium = 	Element(atomic_number=72, name='hafnium', symbol='Hf', mass=178.49 * u.amu)
Tantalum = 	Element(atomic_number=73, name='tantalum', symbol='Ta', mass=180.9479 * u.amu)
Tungsten = 	Element(atomic_number=74, name='tungsten', symbol='W', mass=183.84 * u.amu)
Rhenium = 	Element(atomic_number=75, name='rhenium', symbol='Re', mass=186.207 * u.amu)
Osmium = 	Element(atomic_number=76, name='osmium', symbol='Os', mass=190.23 * u.amu)
Iridium = 	Element(atomic_number=77, name='iridium', symbol='Ir', mass=192.217 * u.amu)
Platinum = 	Element(atomic_number=78, name='platinum', symbol='Pt', mass=195.084 * u.amu)
Gold = 		Element(atomic_number=79, name='gold', symbol='Au', mass=196.9666 * u.amu)
Mercury = 	Element(atomic_number=80, name='mercury', symbol='Hg', mass=200.59 * u.amu)
Thallium = 	Element(atomic_number=81, name='thallium', symbol='Tl', mass=204.3833 * u.amu)
Lead = 		Element(atomic_number=82, name='lead', symbol='Pb', mass=207.2 * u.amu)
Bismuth = 	Element(atomic_number=83, name='bismuth', symbol='Bi', mass=208.9804 * u.amu)
Polonium = 	Element(atomic_number=84, name='polonium', symbol='Po', mass=209.0 * u.amu)
Astatine = 	Element(atomic_number=85, name='astatine', symbol='At', mass=210.0 * u.amu)
Radon = 	Element(atomic_number=86, name='radon', symbol='Rn', mass=222.0 * u.amu)
Francium = 	Element(atomic_number=87, name='francium', symbol='Fr', mass=223.0 * u.amu)
Radium = 	Element(atomic_number=88, name='radium', symbol='Ra', mass=226.0 * u.amu)
Actinium = 	Element(atomic_number=89, name='actinium', symbol='Ac', mass=227.0 * u.amu)
Thorium = 	Element(atomic_number=90, name='thorium', symbol='Th', mass=232.0381 * u.amu)
Proactinium = 	Element(atomic_number=91, name='proactinium', symbol='Pa', mass=231.0359 * u.amu)
Uranium = 	Element(atomic_number=92, name='uranium', symbol='U', mass=238.0289 * u.amu)
Neptunium = 	Element(atomic_number=93, name='neptunium', symbol='Np', mass=237.0 * u.amu)
Plutonium = 	Element(atomic_number=94, name='plutonium', symbol='Pu', mass=244.0 * u.amu)
Americium = 	Element(atomic_number=95, name='americium', symbol='Am', mass=243.0 * u.amu)
Curium = 	Element(atomic_number=96, name='curium', symbol='Cm', mass=247.0 * u.amu)
Berkelium = 	Element(atomic_number=97, name='berkelium', symbol='Bk', mass=247.0 * u.amu)
Californium = 	Element(atomic_number=98, name='californium', symbol='Cf', mass=251.0 * u.amu)
Einsteinium = 	Element(atomic_number=99, name='einsteinium', symbol='Es', mass=252.0 * u.amu)
Fermium = 	Element(atomic_number=100, name='fermium', symbol='Fm', mass=257.0 * u.amu)
Mendelevium = 	Element(atomic_number=101, name='mendelevium', symbol='Md', mass=258.0 * u.amu)
Nobelium = 	Element(atomic_number=102, name='nobelium', symbol='No', mass=259.0 * u.amu)
Lawrencium = 	Element(atomic_number=103, name='lawrencium', symbol='Lr', mass=262.0 * u.amu)
Rutherfordium = Element(atomic_number=104, name='rutherfordium', symbol='Rf', mass=261.0 * u.amu)
Dubnium = 	Element(atomic_number=105, name='dubnium', symbol='Db', mass=262.0 * u.amu)
Seaborgium = 	Element(atomic_number=106, name='seaborgium', symbol='Sg', mass=266.0 * u.amu)
Bohrium = 	Element(atomic_number=107, name='bohrium', symbol='Bh', mass=264.0 * u.amu)
Hassium = 	Element(atomic_number=108, name='hassium', symbol='Hs', mass=277.0 * u.amu)
Meitnerium = 	Element(atomic_number=109, name='meitnerium', symbol='Mt', mass=268.0 * u.amu)
Darmstadtium = 	Element(atomic_number=110, name='darmstadtium', symbol='Ds', mass=281.0 * u.amu)
Roentgenium = 	Element(atomic_number=111, name='roentgenium', symbol='Rg', mass=272.0 * u.amu)
Copernicium = 	Element(atomic_number=112, name='copernicium', symbol='Cn', mass=285.0 * u.amu)
Ununtrium = 	Element(atomic_number=113, name='ununtrium', symbol='Uut', mass=284.0 * u.amu)
Ununquadium = 	Element(atomic_number=114, name='ununquadium', symbol='Uuq', mass=289.0 * u.amu)
Ununpentium = 	Element(atomic_number=115, name='ununpentium', symbol='Uup', mass=288.0 * u.amu)
Ununhexium = 	Element(atomic_number=116, name='ununhexium', symbol='Uuh', mass=292.0 * u.amu)
Ununseptium = 	Element(atomic_number=117, name='ununseptium', symbol='Uus', mass=291.0 * u.amu)
Ununoctium = 	Element(atomic_number=118, name='ununoctium', symbol='Uuo', mass=294.0 * u.amu)

Periodict = {  Hydrogen:None, Helium:None, Lithium:None, Beryllium:None, Boron:None,
	       Carbon:None, Nitrogen:None, Oxygen:None, Fluorine:None, Neon:None,
               Sodium:None, Magnesium:None, Aluminum:None, Silicon:None, Phosphorus:None,
               Sulfur:None, Chlorine:None, Argon:None, Potassium:None, Calcium:None,
               Scandium:None, Titanium:None, Vanadium:None, Chromium:None, Manganese:None,
               Iron:None, Cobalt:None, Nickel:None, Copper:None, Zinc:None, Gallium:None,
               Germanium:None, Arsenic:None, Selenium:None, Bromine:None, Krypton:None,
               Rubidium:None, Strontium:None, Yttrium:None, Zirconium:None, Niobium:None,
               Molybdenum:None, Technetium:None, Ruthenium:None, Rhodium:None, Palladium:None,
               Silver:None, Cadmium:None, Indium:None, Tin:None, Antimony:None, Tellurium:None,
               Iodine:None, Xenon:None, Cesium:None, Barium:None, Lanthanum:None, Cerium:None,
               Praseodymium:None, Neodymium:None, Promethium:None, Samarium:None, Europium:None,
               Gadolinium:None, Terbium:None, Dysprosium:None, Holmium:None, Erbium:None, Thulium:None,
               Ytterbium:None, Lutetium:None, Hafnium:None, Tantalum:None, Tungsten:None, Rhenium:None,
               Osmium:None, Iridium:None, Platinum:None, Gold:None, Mercury:None, Thallium:None, Lead:None,
               Bismuth:None, Polonium:None, Astatine:None, Radon:None, Francium:None, Radium:None, Actinium:None,
               Thorium:None, Proactinium:None, Uranium:None, Neptunium:None, Plutonium:None, Americium:None,
               Curium:None, Berkelium:None, Californium:None, Einsteinium:None, Fermium:None, Mendelevium:None,
               Nobelium:None, Lawrencium:None, Rutherfordium:None, Dubnium:None, Seaborgium:None, Bohrium:None,
               Hassium:None, Meitnerium:None, Darmstadtium:None, Roentgenium:None, Copernicium:None, Ununtrium:None,
               Ununquadium:None, Ununpentium:None, Ununhexium:None, Ununseptium:None, Ununoctium:None}

