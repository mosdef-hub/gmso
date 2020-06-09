from typing import ClassVar
from pydantic import BaseModel

import pytest

from gmso.abc.auto_doc import apply_docs
from gmso.tests.base_test import BaseTest


class InheritedBaseModelNoBaseDocAttribute(BaseModel):
    name: str = 'BaseModelNoBaseDocAttribute'


class InheritedBaseModelWithBaseDocButNoFlag(BaseModel):
    __base_doc__: ClassVar[str] = """Base Documentation"""
    name: str = 'BaseModelBaseDocButNoFlag'


class InheritedBaseModelEmptyBaseDoc(BaseModel):
    __base_doc__: ClassVar[str] = ''
    __docs_generated__: ClassVar[bool] = False
    name: str = 'InheritedBaseModelEmptyBaseDoc'


class InheritedBaseModelWithBaseDoc(BaseModel):
    __base_doc__: ClassVar[str] = """Inherited BaseModel with base documentation
    Notes
    -----
    This is a sample note
    
    Warnings
    --------
    This is a sample warning.
    
    Examples
    --------
        >>> model = InheritedBaseModelWithBaseDoc()
    """

    __docs_generated__: bool = False


class TestAutoDoc(BaseTest):

    def test_no_base_doc_attribute(self):
        with pytest.raises(AttributeError):
            apply_docs(
                InheritedBaseModelNoBaseDocAttribute,
                map_names=False,
                silent=False
            )

    def test_base_doc_but_no_flag(self):
        with pytest.raises(AttributeError):
            apply_docs(
                InheritedBaseModelNoBaseDocAttribute,
                map_names=False,
                silent=False
            )

    def test_apply_docs_silent(self):
        apply_docs(
            InheritedBaseModelNoBaseDocAttribute,
            map_names=True,
            silent=True
        )

    def test_empty_base_doc(self):
        apply_docs(
            InheritedBaseModelEmptyBaseDoc,
            map_names=False,
            silent=False
        )
        print(InheritedBaseModelEmptyBaseDoc.__doc__)

    def test_inherited_base_model_with_base_doc(self):
        apply_docs(
            InheritedBaseModelWithBaseDoc,
            map_names=False,
            silent=True
        )
        assert '__init__' in InheritedBaseModelWithBaseDoc.__doc__
        assert getattr(InheritedBaseModelWithBaseDoc, '__docs_generated__')
