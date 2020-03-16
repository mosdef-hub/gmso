import glob
import os
import json

from gmso.core.potential import Potential
from gmso.utils.singleton import Singleton
from gmso.exceptions import GMSOError

POTENTIAL_JSONS = glob.glob(os.path.join(os.path.dirname(__file__), 'jsons', '*.json'))
JSON_DIR = os.path.join(os.path.dirname(__file__), 'jsons')


def _create_hidden_dir(new_dir=None):
    gmso_dir = os.path.join(os.path.expanduser('~'), '.gmso')
    if os.path.exists(gmso_dir) and os.path.isdir(gmso_dir):
        return
    else:
        os.makedirs(gmso_dir)
    if new_dir is not None:
        dir_loc = os.path.join(gmso_dir, new_dir)
        if os.path.exists(dir_loc) and os.path.isdir(dir_loc):
            os.makedirs(dir_loc)
        return dir_loc
    return gmso_dir


def _verify_potential_template_keys(_dict, name):
    assert 'name' in _dict, f"Key name not found in the potential template {name}.json"
    assert 'expression' in _dict, f"Key expression not found in the potential template {name}.json"
    assert 'independent_variables' in _dict, f"Key independent_variables not found in the potential template {name}.json"
    if str(name) != _dict['name']:
        raise GMSOError(f'Mismatch between Potential name {name} and {_dict["name"]}')


class PotentialTemplate(Potential):
    def __init__(self,
                 name='PotentialTemplate',
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 independent_variables='r',
                 template=True):
        if not isinstance(independent_variables, set):
            independent_variables = set(independent_variables.split(','))

        super(PotentialTemplate, self).__init__(
            name=name,
            expression=expression,
            independent_variables=independent_variables,
            template=template,
        )


class PotentialTemplates(Singleton):
    """A singleton collection of all the potential templates"""

    def __init__(self):
        try:
            self.json_refs
        except AttributeError:
            self.json_refs = POTENTIAL_JSONS
            potential_names = [(os.path.split(filename)[-1]).replace('.json', '') for filename in POTENTIAL_JSONS]
            self._ref_dict = {potential.replace('.json', ''): None for potential in potential_names}

    def __getitem__(self, item):
        if item not in self._ref_dict:
            raise KeyError(f'No Potential named {item} found in the library')
        if self._ref_dict[item] is None:
            potential_dict = self._load_json(item)
            self._ref_dict[item] = PotentialTemplate(**potential_dict)

        return self._ref_dict[item]

    def _load_json(self, item):
        with open(os.path.join(JSON_DIR, f'{item}.json')) as json_file:
            potential_dict = json.load(json_file)
            _verify_potential_template_keys(potential_dict, item)
            return potential_dict

    @staticmethod
    def save_potential_template(name, template_dict, update=True, user_template=True):
        """Add a new potential template to the PotentialTemplates Library"""
        parent_dir = JSON_DIR
        if user_template:
            parent_dir = _create_hidden_dir('potential_templates')

        _verify_potential_template_keys(template_dict, name)
        with open(os.path.join(parent_dir, f'{name}.json')) as json_file:
            json.dump(template_dict, json_file)

        if update:
            delattr(PotentialTemplates(), 'json_refs')
            PotentialTemplates()
