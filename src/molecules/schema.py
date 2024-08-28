import re
from pydantic import BaseModel, Field, field_validator


class MoleculeResponse(BaseModel):
    id: int
    smiles: str = Field(
        ..., min_length=1, max_length=100, description="structure of chemical molecules"
    )

    @field_validator("smiles")
    def validate_smiles(cls, value):
        if not re.match(
            r"^(?:[A-Z][a-z]?|[a-z])(?:(?:[1-9]\d*)?(?:\[(?:(?:[A-Z][a-z]?(?:@[@]?)?)"
            r"|[#+-]|\d+)?\])?|(?:[-=#$:/\\])?(?:[A-Z][a-z]?|[a-z])|[().\[\]])*((?:[1-9]\d*)?)$",
            value,
        ):
            raise ValueError("This SMILES has an invalid structure")
        return value

class MoleculeAdd(BaseModel):
    smiles: str = Field(
        ..., min_length=1, max_length=100, description="structure of chemical molecules"
    )

    @field_validator("smiles")
    def validate_smiles(cls, value):
        if not re.match(
            r"^(?:[A-Z][a-z]?|[a-z])(?:(?:[1-9]\d*)?(?:\[(?:(?:[A-Z][a-z]?(?:@[@]?)?)"
            r"|[#+-]|\d+)?\])?|(?:[-=#$:/\\])?(?:[A-Z][a-z]?|[a-z])|[().\[\]])*((?:[1-9]\d*)?)$",
            value,
        ):
            raise ValueError("This SMILES has an invalid structure")
        return value
    
class MoleculeUpdate(BaseModel):
    smiles: str = Field(
        ..., min_length=1, max_length=100, description="structure of chemical molecules"
    )

    @field_validator("smiles")
    def validate_smiles(cls, value):
        if not re.match(
            r"^(?:[A-Z][a-z]?|[a-z])(?:(?:[1-9]\d*)?(?:\[(?:(?:[A-Z][a-z]?(?:@[@]?)?)"
            r"|[#+-]|\d+)?\])?|(?:[-=#$:/\\])?(?:[A-Z][a-z]?|[a-z])|[().\[\]])*((?:[1-9]\d*)?)$",
            value,
        ):
            raise ValueError("This SMILES has an invalid structure")
        return value