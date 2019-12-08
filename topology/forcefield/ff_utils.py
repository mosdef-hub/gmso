import os


from lxml import etree
from lxml.etree import DocumentInvalid

import unyt as u


UNITS_MAP = {
    'amu': u.gram / u.mol
}


def validate(xml_path, schema=None):
    """Validate a given xml file with a refrence schema"""
    if schema is None:
        schema_path = os.path.join(os.path.abspath(__file__), 'schema', 'article-schema.xsd')

    xml_doc = etree.parse(schema_path)
    xmlschema = etree.XMLSchema(xml_doc)

    xmlschema.assertValid
