from fastapi_filter.contrib.beanie import Filter
from roppy.api.documents import MonomerSummaryDocument
from typing import Optional


class MonomerSummaryFilter(Filter):
    monomer_id: Optional[str] = None 
    monomer_id__in: Optional[list[str]] = None
    smiles: Optional[str] = None
    smiles__in: Optional[list[str]] = None
    ring_size: Optional[int] = None
    ring_size__lt: Optional[int] = None
    ring_size__gte: Optional[int] = None
    has_exp: Optional[bool] = None
    has_calc: Optional[bool] = None
    search: Optional[str] = None
    order_by: list[str] = ["monomer_id"]

    class Constants(Filter.Constants):
        model = MonomerSummaryDocument
        search_model_fields = ["smiles"]
