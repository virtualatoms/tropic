"""Beanie documents for the API."""

from beanie import Document

from tropic.core.models import MonomerSummary, Polymerisation


class MonomerSummaryDocument(Document, MonomerSummary):
    """Document representing a summary of monomers."""

    class Settings:
        """Settings for the MonomerSummary document."""

        name = "monomers_summary"
        indexes = ("monomer_id", "monomer.smiles")


class PolymerisationDocument(Document, Polymerisation):
    """Document representing a polymerisation record."""

    class Settings:
        """Settings for the Polymerisation document."""

        name = "polymerizations"
        indexes = ("polymerisation_id", "monomer.smiles")
