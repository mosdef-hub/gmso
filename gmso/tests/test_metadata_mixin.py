import pytest

from gmso.abc.metadata_mixin import MetadataMixin
from gmso.tests.base_test import BaseTest


class TestMetadataMixin(BaseTest):
    @pytest.fixture(scope="session")
    def metadata_mixin(self):
        return MetadataMixin()

    def test_metadata_empty_tags(self, metadata_mixin):
        assert metadata_mixin.tag_names == []
        assert list(metadata_mixin.tag_names_iter) == []

    def test_metadata_add_tags(self, metadata_mixin):
        metadata_mixin.add_tag("tag1", dict([("tag_name_1", "value_1")]))
        metadata_mixin.add_tag("tag2", dict([("tag_name_2", "value_2")]))
        metadata_mixin.add_tag("int_tag", 1)
        assert len(metadata_mixin.tag_names) == 3

    def test_metadata_mixin_add_tags_overwrite(self, metadata_mixin):
        with pytest.raises(ValueError):
            metadata_mixin.add_tag("tag2", "new_value", overwrite=False)
        metadata_mixin.add_tag("tag2", "new_value", overwrite=True)
        assert metadata_mixin.get_tag("tag2") == "new_value"
        assert len(metadata_mixin.tag_names) == 3

    def test_metadata_get_tags(self, metadata_mixin):
        assert metadata_mixin.get_tag("tag1").get("tag_name_1") == "value_1"
        assert metadata_mixin.get_tag("int_tag") == 1
        assert metadata_mixin.get_tag("non_existent_tag") is None
        with pytest.raises(KeyError):
            metadata_mixin.get_tag("non_existent_tag", throw=True)

    def test_metadata_mixin_all_tags(self, metadata_mixin):
        assert "int_tag" in metadata_mixin.tags

    def test_metadata_mixin_delete_tags(self, metadata_mixin):
        with pytest.raises(KeyError):
            metadata_mixin.delete_tag("non_existent_tag")
        assert metadata_mixin.pop_tag("non_existent_tag") is None
        metadata_mixin.delete_tag("int_tag")
        assert len(metadata_mixin.tag_names) == 2
