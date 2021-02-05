from typing import Dict, Any, List, Iterator

from pydantic import Field, BaseModel, validator


class MetadataMixin(BaseModel):
    tags: Dict[str, Any] = Field(
        {},
        description='Tags associated with the metadata'
    )

    @property
    def tag_names(self) -> List[str]:
        return list(self.__dict__.get('tags'))

    @property
    def tag_names_iter(self) -> Iterator[str]:
        return iter(self.__dict__.get('tags'))

    def add_tag(self, tag: str, value: Any, overwrite=True) -> None:
        """Add metadata for a particular tag"""
        if self.tags.get(tag) is None and not overwrite:
            raise ValueError(
                f'Tag {tag} already exists. '
                f'Please use overwrite=True to overwrite'
            )
        self.tags[tag] = value

    def delete_tag(self, tag: str) -> None:
        del self.tags[tag]

    def pop_tag(self, tag: str) -> Any:
        return self.tags.pop(tag)

    @validator('tags', pre=True)
    def validate_tags(cls, value):
        if value is None:
            value = dict()
        return value

    class Config:
        validate_assignment = True
