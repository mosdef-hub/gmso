import re
from typing import Type, Any, Dict, Tuple, Union, Optional, List
from pydantic import BaseModel


BASE_DOC_ATTR = '__base_doc__'
FIELDS_IN_DOCSTRING = '__alias_to_fields__'
DOCS_GENERATED = '__docs_generated__'
FIELDS_KEY = '__fields__'

__all__ = ['AutoDocGenerator', 'apply_docs']


def _infer_type(type_: Type) -> str:
    """Infer types based on typing or builtin"""
    # FIXME: Inferred types clunky pydantic.main.ModelMetaclass
    if type(type_) is type:
        return type_.__name__
    else:
        return repr(type_).replace('typing.', '')


def _find_extended_summary(doc_lines: List) -> str:
    """Find extended summary from the docstring if it exists"""
    if len(doc_lines) == 2:
        return doc_lines[1]

    first_section_index = None
    for i, line in enumerate(doc_lines):
        if re.match(r'-+', line):
            first_section_index = i - 1
            break

    if first_section_index:
        return '\n'.join(doc_lines[1: first_section_index])
    else:
        return ''


def _find_section_at(doc_lines: List, name: str) -> str:
    try:
        sec_index_start = doc_lines.index(name)
        sec_index_end = -1
        for j in range(sec_index_start+2, len(doc_lines)):
            if re.match(r'-+', doc_lines[j]):
                sec_index_end = j-1
                break

        if sec_index_end != -1:
            return '\n'.join(doc_lines[sec_index_start: sec_index_end])
        else:
            return '\n'.join(doc_lines[sec_index_start:])

    except ValueError:
        return ''


def _base_doc_to_sections(base_doc: str) -> dict:
    lines = list(filter(lambda x: x.strip() != '', base_doc.split('\n')))
    doc_lines = [line.strip() for line in lines]
    leading_spaces = [(len(line) - len(line.lstrip()) - 4) for line in lines]
    # Hack but works
    doc_lines = [f'\t{doc_lines[j]}' if leading_spaces[j] > 0 else doc_lines[j]
                      for j in range(len(doc_lines))]

    sections = {
            'Parameters': '',
            'Attributes': '',
            'Notes': '',
            'Examples': '',
            'See Also': '',
            'References': '',
            'Warnings': ''
        }

    summary = ''
    if len(doc_lines) > 0:
        summary = doc_lines[0]

    extended_summary = _find_extended_summary(doc_lines)

    for section in sections:
        sections[section] = _find_section_at(doc_lines, section)

    sections.update(
        {
            'Summary': summary,
            'Extended Summary': extended_summary
        }
    )
    return sections


def _inject_parameters_from_fields(fields: Dict[str, Any],
                                   by_alias: bool = True,
                                   name_map: Optional[Tuple[str, str]] = None,
                                   inject_init: bool = True) -> Union[Tuple[str, str], str]:
    """Inject the parameters name from the fields"""
    parameters = ['Parameters', '----------']
    init_signature = None
    if inject_init:
        init_signature = ['__init__(', ]

    for key, value in fields.items():
        default = value.default
        if default == '':
            default = "''"

        if by_alias:
            name = value.alias
        else:
            name = key

        parameters.append(f'{name} : '
                          f'{_infer_type(value.type_)}'
                          f', default={default}')
        if init_signature:
            init_signature.append(f'{name} : {_infer_type(value.type_)} = {default}, ')

        if value.field_info.description:
            parameters.append(f'\t{value.field_info.description}')

    if len(parameters) > 2:
        parameters = '\n'.join(parameters)
    else:
        parameters = ''

    if init_signature:
        init_signature[-1] = init_signature[-1].replace(',', '')
        init_signature.append(')')
        init_signature = ''.join(init_signature)

    # TODO: Will this always be lower? Could use `re`
    if name_map:
        parameters = parameters.replace(
            name_map[1].lower(),
            name_map[0].lower()
        )
    if init_signature:
        return parameters, init_signature

    return parameters


class AutoDocGenerator:
    """Generates __doc__ attribute for pydantic base model classses and descendants

    Parameters
    ----------
    target: pydantic.BaseModel or its children
        The basemodel class to generate __doc__ attribute for
    map_names: bool, default=True
        If true, change field descriptions by replacing superclass name by subclass
    """
    def __init__(self, target: Type[BaseModel], map_names=True) -> None:
        self.should_apply = True
        if target is BaseModel or issubclass(target, BaseModel):
            if hasattr(target, BASE_DOC_ATTR):
                if getattr(target, DOCS_GENERATED):
                    self.should_apply = False
            else:
                raise AttributeError(f'No mating attribute {BASE_DOC_ATTR} found in {target.__name__}')

        else:
            raise TypeError('Cannot generate documentation for non-basemodel descendants')
        self.base_doc = getattr(target, BASE_DOC_ATTR)
        self.target = target
        setattr(self.target, DOCS_GENERATED, True)
        self.name_map = self._name_map(map_names=map_names)
        self.should_inject_init = self._should_inject_init()

    def _should_inject_init(self) -> bool:
        immediate_super_class = self.target.mro()[1]
        self_init = self.target.__init__
        super_init = immediate_super_class.__init__
        return id(self_init) == id(super_init)

    def _name_map(self, map_names: bool = True) -> Union[Tuple[str, str], None]:
        super_classes = self.target.mro()
        if map_names and len(super_classes) > 2:
            return (
                self.target.__name__,
                super_classes[1].__name__
            )
        else:
            return None

    def __get__(self, *args: Any, **kwargs: Any) -> str:
        """Return the __doc__ attribute"""
        return self.get_docstring(
            getattr(self.target, FIELDS_KEY),
            self.base_doc,
            self.name_map,
            self.should_inject_init
        )

    def __del__(self) -> None:
        self.target.__doc__ = self.base_doc
        setattr(self.target, DOCS_GENERATED, False)
        self.should_apply = True

    @staticmethod
    def get_docstring(fields: Dict[str, Any],
                      base_doc: str,
                      name_map: Tuple[str, str],
                      inject_init_signature: bool = True) -> str:
        """Get the docstring names based on fields"""
        sections = _base_doc_to_sections(base_doc)
        params_or_params_init = _inject_parameters_from_fields(
            fields,
            by_alias=True,
            name_map=name_map,
            inject_init=inject_init_signature
        )
        docstring = []

        if inject_init_signature:
            sections['Parameters'] = params_or_params_init[0]
            docstring.append(params_or_params_init[1])
        else:
            sections['Parameters'] = params_or_params_init

        docstring.append(sections.pop('Summary'))
        extended_summary = sections.pop('Extended Summary')

        if extended_summary:
            docstring.append(extended_summary)

        for section in sections:
            sec_str = sections.get(section)
            if sec_str:
                docstring.append(sec_str)

        return '\n\n'.join(docstring)


def apply_docs(target_class: Type[BaseModel],
               map_names: bool = True,
               silent: bool = True) -> None:
    """Apply __doc__ attribute to the BaseModel and its descendants"""
    if hasattr(target_class, DOCS_GENERATED) and getattr(target_class, DOCS_GENERATED):
        return
    try:
        target_class.__doc__ = AutoDocGenerator(target_class, map_names=map_names)
    except (TypeError, ValueError, AttributeError) as e:
        if not silent:
            raise e
