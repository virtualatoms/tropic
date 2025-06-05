from fastapi_filter.contrib.beanie import Filter
from fastapi_filter import FilterDepends, with_prefix
from roppy.api.documents import MonomerSummaryDocument
from roppy.core.models import Monomer
from typing import Optional



class MonomerFilter(Filter):
    smiles: Optional[str] = None
    smiles__in: Optional[list[str]] = None
    molecular_weight__lt: Optional[float] = None
    molecular_weight__gt: Optional[float] = None
    functional_groups__in: Optional[list[str]] = None
    ring_size__lt: Optional[int] = None
    ring_size__gte: Optional[int] = None
    has_exp: Optional[bool] = None
    has_calc: Optional[bool] = None
    search: Optional[str] = None
    order_by: list[str] = ["smiles"]

    class Constants(Filter.Constants):
        model = Monomer

class MonomerSummaryFilter(Filter):
    monomer: Optional[MonomerFilter] = FilterDepends(with_prefix("monomer", MonomerFilter))
    monomer_id: Optional[str] = None
    monomer_id__in: Optional[list[str]] = None
    has_exp: Optional[bool] = None
    has_calc: Optional[bool] = None
    search: Optional[str] = None
    order_by: list[str] = ["monomer_id"]

    class Constants(Filter.Constants):
        model = MonomerSummaryDocument
        search_model_fields = ["monomer.smiles"]