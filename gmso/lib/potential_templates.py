import glob
import os
import json

from gmso.core.potential import Potential
from gmso.utils.singleton import Singleton
from gmso.exceptions import GMSOError

POTENTIAL_JSONS = glob.glob(os.path.join(os.path.dirname(__file__), 'jsons', '*.json'))
JSON_DIR = os.path.join(os.path.dirname(__file__), 'jsons')
USER_TEMPLATE_PREFIX = 'potential_templates'


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


def _load_template_json(item, json_dir=JSON_DIR):
    with open(os.path.join(json_dir, f'{item}.json')) as json_file:
        potential_dict = json.load(json_file)
        _verify_potential_template_keys(potential_dict, item)
        return potential_dict


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
            self.user_jsons = []
            potential_names = [(os.path.split(filename)[-1]).replace('.json', '') for filename in POTENTIAL_JSONS]
            self._ref_dict = {potential.replace('.json', ''): None for potential in potential_names}
            self._user_ref_dict = {}

    def get_available_template_names(self):
        return tuple(self._ref_dict.keys())

    def load_user_templates(self):
        gmso_dir = _create_hidden_dir(USER_TEMPLATE_PREFIX)
        user_jsons = glob.glob(os.path.join(gmso_dir, '*.json'))
        self.user_jsons = list(set(self.user_jsons.extend(user_jsons)))
        user_potential_names = [(os.path.split(filename)[-1]).replace('.json', '') for filename in user_jsons]
        self._user_ref_dict.update({potential.replace('.json', ''): self._user_ref_dict.get(potential, None)
                                    for potential in user_potential_names})

    def __getitem__(self, item):
        exists_user = True
        exists_package = True
        if item not in self._ref_dict:
            exists_package = False
        if item not in self._user_ref_dict:
            exists_user = False

        if exists_user:
            return self._load_from_json(item, ref='user')
        if exists_package:
            return self._load_from_json(item, ref='package')

        raise KeyError(f'Potential Template {item} not found in user as well as package templates'
                       f'Please make sure the user templates are loaded by using load_user_templates')

    def _load_from_json(self, item, ref='package'):
        if ref == 'package':
            ref_dict = self._ref_dict
            _dir = JSON_DIR
        elif ref == 'user':
            ref_dict = self._user_ref_dict
            _dir = _create_hidden_dir(USER_TEMPLATE_PREFIX)
        else:
            raise GMSOError(f'Cannot load template from {ref}. Please use either package or user.')

        if ref_dict[item] is None:
            potential_dict = _load_template_json(item, _dir)
            ref_dict[item] = PotentialTemplate(**potential_dict)
        return ref_dict[item]

    @staticmethod
    def save_potential_template(name, template_dict, update=True, user_template=True):
        """Add a new potential template to the PotentialTemplates Library"""
        parent_dir = JSON_DIR
        if user_template:
            parent_dir = _create_hidden_dir(USER_TEMPLATE_PREFIX)

        _verify_potential_template_keys(template_dict, name)
        with open(os.path.join(parent_dir, f'{name}.json'), 'w') as json_file:
            json.dump(template_dict, json_file)

        if update:
            delattr(PotentialTemplates(), 'json_refs')
            PotentialTemplates()
