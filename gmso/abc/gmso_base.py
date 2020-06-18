import warnings
from abc import ABC
from typing import Any, ClassVar

from pydantic import BaseModel

from gmso.abc.auto_doc import apply_docs


class GMSOBase(BaseModel, ABC):
    """A BaseClass to all abstract classes in GMSO"""
    __base_doc__: ClassVar[str] = """A base class to all abstract base classes in gmso."""

    __docs_generated__: ClassVar[bool] = False

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
        for super_class in cls.mro():
            if hasattr(super_class, 'Config') and hasattr(super_class.Config, 'alias_to_fields'):
                cls.Config.alias_to_fields.update(super_class.Config.alias_to_fields)
        apply_docs(cls, map_names=True, silent=False)

    class Config:
        arbitrary_types_allowed = True
        alias_to_fields = dict()
        extra = 'forbid'





