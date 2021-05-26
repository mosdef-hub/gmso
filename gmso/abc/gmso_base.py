"""Base model all classes extend."""
import warnings
from abc import ABC
from typing import Any, ClassVar

from pydantic import BaseModel
from pydantic.validators import dict_validator

from gmso.abc.auto_doc import apply_docs


class GMSOBase(BaseModel, ABC):
    """A BaseClass to all abstract classes in GMSO."""

    __base_doc__: ClassVar[
        str
    ] = """A base class to all abstract base classes in gmso."""

    __docs_generated__: ClassVar[bool] = False

    def __hash__(self):
        """Return the unique hash of the object."""
        return id(self)

    def __eq__(self, other):
        """Test if two objects are equivalent."""
        return self is other

    def __setattr__(self, name: Any, value: Any) -> None:
        """Set the attributes of the object."""
        if name in self.__config__.alias_to_fields:
            name = self.__config__.alias_to_fields[name]
        elif name in self.__config__.alias_to_fields.values():
            warnings.warn(
                "Use of internal fields is discouraged. "
                "Please use external fields to set attributes."
            )

        super().__setattr__(name, value)

    @classmethod
    def __init_subclass__(cls, **kwargs):
        """Initialize the subclass of the object."""
        super().__init_subclass__()
        setattr(cls, "__docs_generated__", False)
        for super_class in cls.mro()[1:]:
            if (
                hasattr(super_class, "Config")
                and hasattr(super_class.Config, "alias_to_fields")
                and hasattr(cls.Config, "alias_to_fields")
            ):
                cls.Config.alias_to_fields.update(
                    super_class.Config.alias_to_fields
                )
        apply_docs(cls, map_names=True, silent=False)

    @classmethod
    def validate(cls, value):
        """Ensure that the object is validated before use."""
        if isinstance(value, cls):
            return value
        else:
            return cls(**dict_validator(value))

    @classmethod
    def __get_validators__(cls) -> "CallableGenerator":
        """Get the validators of the object."""
        yield cls.validate

    class Config:
        """Pydantic configuration for base object."""

        arbitrary_types_allowed = True
        alias_to_fields = dict()
        extra = "forbid"
