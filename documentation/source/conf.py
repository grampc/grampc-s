# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'grampc-s'
copyright = '2025, Daniel Landgraf, Andreas Völz, Knut Graichen'
author = 'Daniel Landgraf, Andreas Völz, Knut Graichen'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinxcontrib.bibtex']
bibtex_bibfiles = ['literature.bib']

templates_path = ['_templates']
exclude_patterns = []
numfig = True
math_numfig = True
numfig_secnum_depth = 2
math_eqref_format = "({number})"

# define macros for html math output
mathjax3_config = {
    "tex": {
        "macros": {
            "vm": ["\\boldsymbol{#1}", 1],
            "trans": ["^{\mathrm{T}}"],
            "inv": ["^{-1}"],
            "diag": ["\DeclareMathOperator*{\diag}{diag}"],
            "gp": ["{\mathcal{GP}"],
            "ex": ["\mathds{E}\left[#1\right]", 1],
        },
    }
}

latex_elements = {
    "preamble": r"\input{include.tex.txt}"
}
latex_documents = [("index", "manual.tex", "GRAMPC-S documentation", "Daniel Landgraf, Andreas Völz, Knut Graichen", "manual", False)]
latex_additional_files = ["include.tex.txt"]
latex_domain_indices = [False]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_theme_options = {
    "show_toc_level": 2
}