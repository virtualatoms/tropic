# TROPIC 🌴 Thermodynamic Ring-Opening Polymerisation Information Collection

[![Website](https://img.shields.io/badge/website-polytropic.org-3b82f6?style=flat-square)](https://polytropic.org)
[![Docs](https://img.shields.io/badge/docs-API%20Reference-6366f1?style=flat-square)](https://polytropic.org/docs)
[![License](https://img.shields.io/badge/license-MIT-10b981?style=flat-square)](LICENSE)

> Curated thermodynamic intelligence for polar cyclic monomers — bridging experimental and computational insights for ring-opening polymerisation.

- 🌐 Visit the live platform at [polytropic.org](https://polytropic.org)
- 🧭 Explore the [interactive documentation](https://polytropic.org/docs)
- 💬 Join the discussion by opening an [issue](https://github.com/virtualatoms/tropic/issues)

## Table of Contents

- [Overview](#overview)
- [Highlights](#highlights)
- [Platform Architecture](#platform-architecture)
- [Quick Start](#quick-start)
- [Python API Example](#python-api-example)
- [Documentation & Resources](#documentation--resources)
- [Contributing](#contributing)
- [License](#license)

## Overview

**TROPIC** is a curated database and web platform for collecting, organising, and analysing thermodynamic data related to the ring-opening polymerisation (ROP) of polar cyclic monomers. It supports both experimental and computational submissions so you can rapidly compare enthalpy, entropy, and reaction conditions across monomers and methodologies.

## Highlights

- 📚 **Rich database** of literature-sourced thermodynamic measurements for lactones, cyclic carbonates, and more.
- 🔬 **Metadata-first schema** capturing reaction provenance, methodology, and experimental notes.
- 🔎 **Searchable REST API** for programmatic access via SMILES, ring size, method, and custom filters.
- 📊 **Interactive web interface** with visual exploration, filtering, and data export utilities.
- 🧠 **Insight-ready analytics** to juxtapose experimental and computational datasets side-by-side.
- ⚙️ **Modern stack** powered by Pydantic, FastAPI, MongoDB, Beanie, React, Mantine, and AG Grid.

## Platform Architecture

This repository is a mono-repo containing all TROPIC packages:

- 📦 `tropic-core` — Pydantic models for reactions and thermodynamic entities shared across services.
- 🤝 `tropic-api` — FastAPI service that exposes the REST API and connects to the database.
- 🐍 `tropic-client` — Python client for streamlined programmatic consumption of the API.
- 🕸️ `tropic-web` — React front end delivering the interactive exploration experience.

## Quick Start

The TROPIC ecosystem is published to PyPI as namespace packages so you can install only what you need. For instance, install the Python client with:

```bash
pip install tropic-client
```

Looking to self-host the full stack? Use the provided `docker-compose.yml` to stand up the API and web interface locally.

## Python API Example

Interact with the platform programmatically in just a few lines:

```python
from tropic.client import TropicClient

with TropicClient() as client:
    reactions = client.get_reactions(
        type="ROR",
        monomer__ring_size__gte=10,
        dataset="computational"
    )
    for reaction in reactions:
        print(reaction.monomer.common_name, reaction.delta_h, reaction.delta_s)
```

## Documentation & Resources

- 📄 Dive into the [API reference](https://polytropic.org/docs) for endpoint details and usage notes.
- 🧬 Browse [`tropic-core/src/tropic/core/models.py`](https://github.com/virtualatoms/tropic/blob/main/tropic-core/src/tropic/core/models.py) to inspect the full data schema.
- 🧪 Check the [`data/`](data) directory for curated datasets and examples.

## Contributing

We welcome contributions of all shapes and sizes — from new datasets to tooling improvements.

1. 🧠 Start a conversation by filing an [issue](https://github.com/virtualatoms/tropic/issues) outlining the change or addition.
2. 🔄 Fork the repository, create a feature branch, and make your updates.
3. ✅ Add or update tests and documentation where relevant.
4. 📬 Open a pull request describing the motivation, approach, and any follow-up ideas.

Not sure where to begin? Look for issues labelled `good first issue` or propose a dataset you would like to contribute.

## License

This project is released under the MIT License. See [LICENSE](https://github.com/virtualatoms/tropic/blob/main/LICENSE) for details.

---

Made with ❤️ by the TROPIC team — contributions and feedback keep the rainforest of data thriving! 🌿
