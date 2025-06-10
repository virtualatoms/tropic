"""Beanie documents for the API."""

from beanie import Document, Link
from pydantic import BaseModel, Field

from tropic.core.models import Monomer, Polymerisation


class MonomerDocument(Document, Monomer):
    """Document representing a monomer."""

    class Settings:
        """Settings for the Monomer document."""

        name = "monomers"
        indexes = ("monomer_id", "smiles")


class PolymerisationDocument(Document, Polymerisation):
    """Document representing a polymerisation record."""

    monomer: Link[MonomerDocument] = Field(
        description="Monomer of the polymerisation",
    )

    class Settings:
        """Settings for the Polymerisation document."""

        name = "polymerizations"
        indexes = ("polymerisation_id", "monomer.smiles")


class MonomerSummary(BaseModel):
    monomer_id: str = Field(alias="_id")
    smiles: str
    has_exp: bool
    has_comp: bool
