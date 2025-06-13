# TROPIC: Thermodynamic Ring-Opening Polymerisation Information Collection

**TROPIC** is a curated database and web platform for collecting, organising, and analysing thermodynamic data related to the ring-opening polymerisation (ROP) of polar cyclic monomers. It supports both experimental and computational data, enabling structured comparison of reaction parameters such as enthalpy and entropy of polymerisation across monomers and conditions.

## Features

- 📚 **Database**: Manually curated thermodynamic data from the literature for ROP of cyclic monomers (e.g. lactones, cyclic carbonates).
- 🔬 **Metadata-rich schema**: Includes reaction conditions, methods of determination, and provenance.
- 🔎 **REST API**: Searchable interface for programmatic access to data (e.g. by SMILES, ring size, method).
- 📊 **Web Interface**: Interactive tools for monomer search, visualisation, filtering, and data export.
- ⚙️ **Built with**: Pydantic, FastAPI, MongoDB, Beanie ODM, Plotly Dash, Dash Mantine Components, Dash AG Grid.

## Getting Started

This repository contains the packages for running the tropic website, API server, and API client. The repository is structured as:
- `tropic-core`: Pydantic document models representing reactions, used by the API server and client.
- `tropic-client`: The TROPIC Python API client.
- `tropic-api`: The REST API implementation.
- `tropic-web`: The interactive website implementation.

### Installation

The TROPIC software is available as namespace packages on PyPI. This enables you to install only the specific portion of the software that you need.
For example, the TROPIC API client can be installed using.

```bash
pip install tropic-client
```

### Using the API

The easiet way to use the TROPIC API is through the Python client:

```python
from tropic.client import TropicClient

with TropicClient() as client:
    reactions = client.get_reactions(type="ROR", monomer__ring_size__gte=10)
```

## Documentation

- 📄 The full API documentation available on the [polytropic website](https://polytropic.org/api).
- 📁 See [`roppy-core/models.py`](https://github.com/virtualatoms/tropic/blob/main/tropic-core/src/tropic/core/models.py) for full schema details.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request to contribute new features, bug fixes, or data entries.

## License

This project is released under the MIT License. See [LICENSE](https://github.com/virtualatoms/tropic/blob/main/LICENSE) for details.
