from pydantic import Field
from beanie import Document, Link

from roppy.core.models import MonomerSummary, PolymerisationSummary, Polymerisation


class PolymerisationSummaryDocument(Document, PolymerisationSummary):
    class Settings:
        name = "polymerisations_summary"
        indexes = [
            "polymerisation_summary_id",
            "reaction_smiles",
        ]


class MonomerSummaryDocument(Document, MonomerSummary):
    polymerisations: list[Link[PolymerisationSummary]] = Field(
        default_factory=list,
        description="Links to polymerization data involving the monomer",
    )

    class Settings:
        name = "monomers_summary"
        indexes = [
            "monomer_id",
            "smiles",
        ]


class PolymerizationDocument(Document, Polymerisation):
    class Settings:
        name = "polymerizations"
        indexes = [
            "polymerisation_id",
            "monomer.smiles",
        ]
