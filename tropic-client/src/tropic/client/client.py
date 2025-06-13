"""Module for Tropic API client."""

from dataclasses import dataclass, field
from typing import Any, ClassVar

import requests

from tropic.core.models import Polymerisation


@dataclass
class TropicClient:
    """Client for interacting with the Tropic API.

    Example
    -------

    The client can be used with pythons "with" statement to ensure proper
    resource management::

        from tropic.client import TropicClient

        with TropicClient() as client:
            polymerisations = client.get_polymerisations(type="addition")
    """

    ENDPOINT: str = "http://127.0.0.1:8000/"
    ALLOWED_FIELDS: ClassVar[set[str]] = {
        "polymerisation_id",
        "polymerisation_id__in",
        "type",
        "type__in",
        "monomer__smiles",
        "monomer__smiles__in",
        "monomer__monomer_id",
        "monomer__monomer_id__in",
        "monomer__inchi",
        "monomer__inchi__in",
        "monomer__molecular_weight__lte",
        "monomer__molecular_weight__gte",
        "monomer__functional_groups__in",
        "monomer__iupac_name",
        "monomer__iupac_name__in",
        "monomer__pubchem_cid",
        "monomer__pubchem_cid__in",
        "monomer__ring_size__lte",
        "monomer__ring_size__gte",
        "product__smiles",
        "product__smiles__in",
        "product__repeating_units__lte",
        "product__repeating_units__gte",
        "product__repeating_units",
        "product__deg_of_poly__lte",
        "product__deg_of_poly__gte",
        "product__dispersity__lte",
        "product__dispersity__gte",
        "product__n_avg_molar_mass__lte",
        "product__n_avg_molar_mass__gte",
        "product__m_avg_molar_mass__lte",
        "product__m_avg_molar_mass__gte",
        "thermo__delta_h__lte",
        "thermo__delta_h__gte",
        "thermo__delta_s__lte",
        "thermo__delta_s__gte",
        "thermo__delta_g__lte",
        "thermo__delta_g__gte",
        "thermo__ceiling_temperature__lte",
        "thermo__ceiling_temperature__gte",
        "metadata__year__lte",
        "metadata__year__gte",
        "metadata__doi",
        "metadata__doi__in",
        "order_by",
    }
    session: requests.Session = field(
        default_factory=requests.Session,
        init=False,
        repr=False,
    )

    def __enter__(self) -> "TropicClient":
        """Support for "with" context."""
        return self

    def __exit__(self, *_: object) -> None:
        """Support for "with" context."""
        self.session.close()

    def request(self, sub_url: str, method: str = "GET") -> list | dict[str, Any]:
        """
        Make a request to the Tropic API.

        Parameters
        ----------
        sub_url : str
            The sub-path of the API endpoint to request.
        method : str
            The HTTP method to use for the request (default is "GET").

        Returns
        -------
        Any
            The JSON response from the API.
        """
        url = f"{self.ENDPOINT}{sub_url}"
        response = requests.request(method, url, timeout=5, verify=True)
        response.raise_for_status()
        return response.json()

    def get_polymerisations(
        self,
        **kwargs: dict[str, str | list[str] | None | float | int],
    ) -> list[Polymerisation]:
        """Retrieve polymerisations from the Tropic API.

        Parameters
        ----------
        polymerisation_id : str | None
            Polymerisation ID to filter by.
        polymerisation_id__in : list[str] | None
            List of polymerisation IDs to filter by.
        type : str | None
            Type of polymerisation to filter by.
        type__in : list[str] | None
            List of types to filter by.
        monomer__smiles : str | None
            SMILES string of the monomer to filter by.
        monomer__smiles__in : list[str] | None
            List of SMILES strings to filter by.
        monomer__monomer_id : str | None
            Monomer ID to filter by.
        monomer__monomer_id__in : list[str] | None
            List of monomer IDs to filter by.
        monomer__inchi : str | None
            InChI string of the monomer to filter by.
        monomer__inchi__in : list[str] | None
            List of InChI strings to filter by.
        monomer__molecular_weight__lte : float | None
            Maximum molecular weight of the monomer to filter by.
        monomer__molecular_weight__gte : float | None
            Minimum molecular weight of the monomer to filter by.
        monomer__functional_groups__in : list[str] | None
            List of functional groups to filter by.
        monomer__iupac_name : str | None
            IUPAC name of the monomer to filter by.
        monomer__iupac_name__in : list[str] | None
            List of IUPAC names to filter by.
        monomer__pubchem_cid : str | None
            PubChem CID of the monomer to filter by.
        monomer__pubchem_cid__in : list[str] | None
            List of PubChem CIDs to filter by.
        monomer__ring_size__lte : int | None
            Maximum ring size of the monomer to filter by.
        monomer__ring_size__gte : int | None
            Minimum ring size of the monomer to filter by.
        monomer__ring_size : int | None
            Ring size of the monomer to filter by.
        product__smiles : str | None
            SMILES string of the product to filter by.
        product__smiles__in : list[str] | None
            List of SMILES strings to filter by.
        product__repeating_units__lte : int | None
            Maximum number of repeating units in the product to filter by.
        product__repeating_units__gte : int | None
            Minimum number of repeating units in the product to filter by.
        product__repeating_units : int | None
            Number of repeating units in the product to filter by.
        product__deg_of_poly__lte : float | None
            Maximum degree of polymerisation in the product to filter by.
        product__deg_of_poly__gte : float | None
            Minimum degree of polymerisation in the product to filter by.
        product__dispersity__lte : float | None
            Maximum dispersity of the product to filter by.
        product__dispersity__gte : float | None
            Minimum dispersity of the product to filter by.
        product__n_avg_molar_mass__lte : float | None
            Maximum number average molar mass of the product to filter by.
        product__n_avg_molar_mass__gte : float | None
            Minimum number average molar mass of the product to filter by.
        product__m_avg_molar_mass__lte : float | None
            Maximum mass average molar mass of the product to filter by.
        product__m_avg_molar_mass__gte : float | None
            Minimum mass average molar mass of the product to filter by.
        thermo__delta_h__lte : float | None
            Maximum enthalpy change of the polymerisation to filter by.
        thermo__delta_h__gte : float | None
            Minimum enthalpy change of the polymerisation to filter by.
        thermo__delta_s__lte : float | None
            Maximum entropy change of the polymerisation to filter by.
        thermo__delta_s__gte : float | None
            Minimum entropy change of the polymerisation to filter by.
        thermo__delta_g__lte : float | None
            Maximum Gibbs free energy change of the polymerisation to filter by.
        thermo__delta_g__gte : float | None
            Minimum Gibbs free energy change of the polymerisation to filter by.
        thermo__ceiling_temperature__lte : float | None
            Maximum ceiling temperature of the polymerisation to filter by.
        thermo__ceiling_temperature__gte : float | None
            Minimum ceiling temperature of the polymerisation to filter by.
        metadata__year__lte : int | None
            Maximum year of the polymerisation to filter by.
        metadata__year__gte : int | None
            Minimum year of the polymerisation to filter by.
        metadata__doi : str | None
            DOI of the polymerisation to filter by.
        metadata__doi__in : list[str] | None
            List of DOIs to filter by.
        order_by : str | None
            Field to order the results by.

        Returns
        -------
        list[Polymerisation]
            List of polymerisation records matching the filters.
        """
        if invalid_keys := {key for key in kwargs if key not in self.ALLOWED_FIELDS}:
            raise ValueError(
                f"Invalid query fields: {invalid_keys}."
                f" Must be in: {self.ALLOWED_FIELDS}",
            )

        sub_url = "polymerisations"
        if kwargs:
            query_params = "&".join(
                f"{key}={value}" for key, value in kwargs.items() if value is not None
            )
            sub_url += f"?{query_params}"

        documents = self.request(sub_url=sub_url, method="GET")
        return [Polymerisation(**doc) for doc in documents]

    def get_polymerisation(
        self,
        polymerisation_id: str,
    ) -> Polymerisation:
        """Retrieve a specific polymerisation by its ID.

        Parameters
        ----------
        polymerisation_id : str
            The ID of the polymerisation to retrieve.

        Returns
        -------
        Polymerisation
            The polymerisation record with the specified ID.
        """
        sub_url = f"polymerisations/{polymerisation_id}"
        document = self.request(sub_url=sub_url, method="GET")
        return Polymerisation(**document)
