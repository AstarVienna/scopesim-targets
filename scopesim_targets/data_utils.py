# -*- coding: utf-8 -*-
"""Auxilliary functions for downloading and caching static data."""

from pathlib import Path

import pooch

from . import PKG_DIR, DATA_DIR

RETRIEVER = pooch.create(
    path=Path.home() / ".astar/scopesim-targets",
    base_url="https://github.com/AstarVienna/scopesim-targets/data/",
    # env=None,
    registry=None,  # load afterwards
    retry_if_failed=5,
)
RETRIEVER.load_registry(PKG_DIR / "_static/registry.txt")


def fetch_data_file(filename: str):
    """Load local data or fetch via pooch.

    Data will be available locally if this is running from a cloned repo or
    editable install from such. Pooch will try to cache the downloaded file.
    """
    if DATA_DIR.is_dir() and (path := DATA_DIR / filename).exists():
        return path
    RETRIEVER.fetch(filename, progressbar=True)
