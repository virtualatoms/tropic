from beanie import Document

from tropic.core.models import MonomerSummary, Polymerisation


class MonomerSummaryDocument(Document, MonomerSummary):
    class Settings:
        name = "monomers_summary"
        indexes = ["monomer_id", "monomer.smiles"]


class PolymerisationDocument(Document, Polymerisation):
    class Settings:
        name = "polymerizations"
        indexes = ["polymerisation_id", "monomer.smiles"]
