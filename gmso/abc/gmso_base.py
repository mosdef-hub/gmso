"""Base model all classes extend."""
import warnings
from abc import ABC
from typing import Any, ClassVar, Type

from pydantic import BaseModel, ConfigDict, validators

from gmso.abc.auto_doc import apply_docs
from gmso.abc.serialization_utils import dict_to_unyt

dict_validator = validators.getattr_migration("dict_validator")


class GMSOBase(BaseModel, ABC):
    """A BaseClass to all abstract classes in GMSO."""

    __base_doc__: ClassVar[
        str
    ] = """A base class to all abstract base classes in gmso."""

    __docs_generated__: ClassVar[bool] = False

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        validate_assignment=True,
        extra="forbid",
        populate_by_name=True,
    )

    def __hash__(self):
        """Return the unique hash of the object."""
        return id(self)

    def __eq__(self, other):
        """Test if two objects are equivalent."""
        return self is other

    def __setattr__(self, name: Any, value: Any) -> None:
        """Set the attributes of the object."""
        if name in self.model_config.get("alias_to_fields"):
            name = self.model_config.get("alias_to_fields")[name]
        elif name in self.model_config.get("alias_to_fields").values():
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
    def model_validate(cls: Type["Model"], obj: Any) -> "Model":
        dict_to_unyt(obj)
        return super(GMSOBase, cls).model_validate(obj)

    def model_dump(self, **kwargs) -> "DictStrAny":
        kwargs["by_alias"] = True

        additional_excludes = set()
        if "exclude" in kwargs:
            for term in kwargs["exclude"]:
                if term in self.model_config["alias_to_fields"]:
                    additional_excludes.add(
                        self.model_config["alias_to_fields"][term]
                    )
            kwargs["exclude"] = kwargs["exclude"].union(additional_excludes)
        super_dict = super(GMSOBase, self).model_dump(**kwargs)
        return super_dict

    def model_dump_json(self, **kwargs):
        kwargs["by_alias"] = True

        additional_excludes = set()
        if "exclude" in kwargs:
            for term in kwargs["exclude"]:
                if term in self.model_config["alias_to_fields"]:
                    additional_excludes.add(
                        self.model_config["alias_to_fields"][term]
                    )
            kwargs["exclude"] = kwargs["exclude"].union(additional_excludes)
        super_dict = super(GMSOBase, self).model_dump_json(**kwargs)

        return super_dict

    def _iter(self, **kwargs) -> "TupleGenerator":
        exclude = kwargs.get("exclude")
        include = kwargs.get("include")
        include_alias = set()
        exclude_alias = set()

        if include:
            for included in include:
                if included in self.Config.alias_to_fields:
                    include_alias.add(self.Config.alias_to_fields[included])
                else:
                    include_alias.add(included)
            kwargs["include"] = include_alias

        if exclude:
            for excluded in exclude:
                if excluded in self.Config.alias_to_fields:
                    exclude_alias.add(self.Config.alias_to_fields[excluded])
                else:
                    exclude_alias.add(excluded)
            kwargs["exclude"] = exclude_alias

        yield from super()._iter(**kwargs)

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
