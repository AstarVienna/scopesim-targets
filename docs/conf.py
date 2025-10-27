# -*- coding: utf-8 -*-
"""Sphinx configuration."""

project = "ScopeSim-Targets"

extensions = [
    # "sphinx.ext.todo",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    # "sphinx.ext.intersphinx",
    # "sphinx.ext.inheritance_diagram",
    # "sphinx.ext.mathjax",
    # "sphinx.ext.extlinks",
    # "sphinx.ext.doctest",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "myst_nb",
]

templates_path = ["_templates"]
autosummary_generate = True
autoclass_content = "class"
autodoc_default_flags = ["members", "inherited-members"]
autodoc_docstring_signature = False
napoleon_numpy_docstring = True
napoleon_use_admonition_for_references = True

source_encoding = "utf-8"
source_suffix = {
    ".rst": "restructuredtext",
    ".myst": "myst-nb",
    ".md": "myst-nb",
}

nb_execution_timeout = 3600  # [s]
nb_execution_mode = "auto"

html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": "https://github.com/AstarVienna/scopesim-targets",
    "use_repository_button": True,
    "use_download_button": True,
    "home_page_in_toc": True,
}
