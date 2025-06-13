"""Beanie documents for the API."""

from beanie import Document, Link
from pydantic import Field

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
