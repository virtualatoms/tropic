"""Beanie documents for the API."""

from beanie import Document, Link
from pydantic import Field

from tropic.core.models import Monomer, Reaction


class MonomerDocument(Document, Monomer):
    """Document representing a monomer."""

    class Settings:
        """Settings for the Monomer document."""

        name = "monomers"
        indexes = ("monomer_id", "smiles")


class ReactionDocument(Document, Reaction):
    """Document representing a reaction record."""

    monomer: Link[MonomerDocument] = Field(
        description="Monomer of the reaction",
    )

    class Settings:
        """Settings for the Reaction document."""

        name = "reactions"
        indexes = ("reaction_id", "monomer.smiles")
