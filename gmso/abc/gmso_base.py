import warnings
from abc import ABC
from typing import Any, ClassVar, Type


from pydantic import BaseModel
from pydantic.validators import dict_validator

from gmso.abc import GMSOJSONHandler
from gmso.abc.auto_doc import apply_docs
from gmso.abc.serialization_utils import dict_to_unyt


class GMSOBase(BaseModel, ABC):
    """A BaseClass to all abstract classes in GMSO"""
    __base_doc__: ClassVar[str] = """A base class to all abstract base classes in gmso."""

    __docs_generated__: ClassVar[bool] = False

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def __setattr__(self, name: Any, value: Any) -> None:
        if name in self.__config__.alias_to_fields:
            name = self.__config__.alias_to_fields[name]
        else:
            warnings.warn(
                'Use of internal fields is discouraged. '
                'Please use external fields to set attributes.'
            )

        super().__setattr__(name, value)

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__()
        setattr(cls, '__docs_generated__', False)
        for super_class in cls.mro()[1:]:
            if hasattr(super_class, 'Config') and \
               hasattr(super_class.Config, 'alias_to_fields') and \
               hasattr(cls.Config, 'alias_to_fields'):
                cls.Config.alias_to_fields.update(super_class.Config.alias_to_fields)
        apply_docs(cls, map_names=True, silent=False)

    @classmethod
    def parse_obj(cls: Type['Model'], obj: Any) -> 'Model':
        dict_to_unyt(obj)
        return super(GMSOBase, cls).parse_obj(obj)

    @classmethod
    def validate(cls, value):
        if isinstance(value, cls):
            return value
        else:
            return cls(**dict_validator(value))

    @classmethod
    def __get_validators__(cls) -> 'CallableGenerator':
        yield cls.validate

    class Config:
        arbitrary_types_allowed = True
        alias_to_fields = dict()
        extra = 'forbid'
        json_encoders = GMSOJSONHandler.json_encoders
