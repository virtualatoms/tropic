"""Filters for the API."""

from typing import ClassVar

from fastapi_filter import FilterDepends, with_prefix
from fastapi_filter.contrib.beanie import Filter

from tropic.core.models import Metadata, Monomer, Parameters, Product, Reaction, Thermo


class MonomerFilter(Filter):
    """Filter for monomers based on various criteria."""

    smiles: str | None = None
    smiles__in: list[str] | None = None
    monomer_id: str | None = None
    monomer_id__in: list[str] | None = None
    inchi: str | None = None
    inchi__in: list[str] | None = None
    molecular_weight__lte: float | None = None
    molecular_weight__gte: float | None = None
    functional_groups__in: list[str] | None = None
    iupac_name: str | None = None
    iupac_name__in: list[str] | None = None
    pubchem_cid: str | None = None
    pubchem_cid__in: list[str] | None = None
    ring_size: int | None = None
    ring_size__lte: int | None = None
    ring_size__gte: int | None = None

    class Constants(Filter.Constants):
        """Settings for the Monomer filter."""

        model = Monomer


class ProductFilter(Filter):
    """Filter for monomers based on various criteria."""

    smiles: str | None = None
    smiles__in: list[str] | None = None
    repeating_units__lte: int | None = None
    repeating_units__gte: int | None = None
    repeating_units: int | None = None
    deg_of_poly__lte: float | None = None
    deg_of_poly__gte: float | None = None
    dispersity__lte: float | None = None
    dispersity__gte: float | None = None
    n_avg_molar_mass__lte: float | None = None
    n_avg_molar_mass__gte: float | None = None
    m_avg_molar_mass__lte: float | None = None
    m_avg_molar_mass__gte: float | None = None

    class Constants(Filter.Constants):
        """Settings for the Monomer filter."""

        model = Product


class ParametersFilter(Filter):
    """Filter for reaction parameters based on various criteria."""

    is_experimental: bool | None = None
    is_experimental__in: list[bool] | None = None
    temperature__lte: float | None = None
    temperature__gte: float | None = None
    pressure__lte: float | None = None
    pressure__gte: float | None = None
    monomer_state: str | None = None
    monomer_state__in: list[str] | None = None
    polymer_state: str | None = None
    polymer_state__in: list[str] | None = None
    initiator_smiles: str | None = None
    initiator_smiles__in: list[str] | None = None
    initial_monomer_conc__lte: float | None = None
    initial_monomer_conc__gte: float | None = None
    bulk_monomer_conc__lte: float | None = None
    bulk_monomer_conc__gte: float | None = None
    medium: str | None = None
    medium__in: list[str] | None = None
    solvent: str | None = None
    solvent__in: list[str] | None = None
    cosolvent: str | None = None
    cosolvent__in: list[str] | None = None
    method: str | None = None
    method__in: list[str] | None = None
    functional: str | None = None
    functional__in: list[str] | None = None
    basis_set: str | None = None
    basis_set__in: list[str] | None = None
    dispersion: str | None = None
    dispersion__in: list[str] | None = None
    forcefield: str | None = None
    forcefield__in: list[str] | None = None
    solvent_model: str | None = None
    solvent_model__in: list[str] | None = None
    topology: str | None = None
    topology__in: list[str] | None = None

    class Constants(Filter.Constants):
        """Settings for the Parameters filter."""

        model = Parameters


class ThermoFilter(Filter):
    """Filter for reaction thermodynamic data based on various criteria."""

    delta_h__lte: float | None = None
    delta_h__gte: float | None = None
    delta_s__lte: float | None = None
    delta_s__gte: float | None = None
    delta_g__lte: float | None = None
    delta_g__gte: float | None = None
    ceiling_temperature__lte: float | None = None
    ceiling_temperature__gte: float | None = None

    class Constants(Filter.Constants):
        """Settings for the Thermo filter."""

        model = Thermo


class MetadataFilter(Filter):
    """Filter for reaction metadata based on various criteria."""

    year__lte: int | None = None
    year__gte: int | None = None
    doi: str | None = None
    doi__in: list[str] | None = None

    class Constants(Filter.Constants):
        """Settings for the Metadata filter."""

        model = Metadata


class ReactionFilter(Filter):
    """Filter for reactions based on various criteria."""

    reaction_id: str | None = None
    reaction_id__in: list[str] | None = None
    type: str | None = None
    type__in: list[str] | None = None
    monomer: MonomerFilter | None = FilterDepends(with_prefix("monomer", MonomerFilter))
    product: ProductFilter | None = FilterDepends(with_prefix("product", ProductFilter))
    parameters: ParametersFilter | None = FilterDepends(
        with_prefix("parameters", ParametersFilter),
    )
    thermo: ThermoFilter | None = FilterDepends(with_prefix("thermo", ThermoFilter))
    metadata: MetadataFilter | None = FilterDepends(
        with_prefix("metadata", MetadataFilter),
    )
    order_by: ClassVar[list[str]] = ["reaction_id"]

    class Constants(Filter.Constants):
        """Settings for the Reaction filter."""

        model = Reaction


class MonomerSummariesFilter(Filter):
    """Filter for monomer summaries based on various criteria."""

    monomer: MonomerFilter | None = FilterDepends(with_prefix("monomer", MonomerFilter))

    class Constants(Filter.Constants):
        """Settings for the MonomerSummaries filter."""

        model = Reaction
