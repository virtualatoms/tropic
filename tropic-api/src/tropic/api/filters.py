"""Filters for the API."""

from typing import ClassVar

from fastapi_filter import FilterDepends, with_prefix
from fastapi_filter.contrib.beanie import Filter

from tropic.api.documents import MonomerSummaryDocument
from tropic.core.models import Monomer


class MonomerFilter(Filter):
    """Filter for monomers based on various criteria."""

    smiles: str | None = None
    smiles__in: list[str] | None = None
    molecular_weight__lte: float | None = None
    molecular_weight__gte: float | None = None
    functional_groups__in: list[str] | None = None
    ring_size__lte: int | None = None
    ring_size__gte: int | None = None
    has_exp: bool | None = None
    has_calc: bool | None = None
    search: str | None = None
    order_by: ClassVar[list[str]] = ["smiles"]

    class Constants(Filter.Constants):
        """Settings for the Monomer filter."""

        model = Monomer


class MonomerSummaryFilter(Filter):
    """Filter for monomer summaries based on various criteria."""

    monomer: MonomerFilter | None = FilterDepends(with_prefix("monomer", MonomerFilter))
    monomer_id: str | None = None
    monomer_id__in: list[str] | None = None
    has_exp: bool | None = None
    has_calc: bool | None = None
    search: str | None = None
    order_by: ClassVar[list[str]] = ["monomer_id"]

    class Constants(Filter.Constants):
        """Settings for the MonomerSummary filter."""

        model = MonomerSummaryDocument
        search_model_fields = ("monomer.smiles",)
