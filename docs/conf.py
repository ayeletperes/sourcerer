#!/usr/bin/env python3
#
# sourcerer documentation build configuration file
#
# This file is execfile()d with the current directory set to its
# containing dir.

import datetime
import os

# Python 3.14's argparse colorizes usage/help text by default when it thinks
# stdout is a terminal (e.g. when `make html` is run interactively). That
# color is just ANSI escape codes, and sphinxcontrib-autoprogram embeds
# argparse's format_usage()/format_help() output verbatim into the docs, so
# without this the escape codes show up as literal text in the rendered
# usage blocks. NO_COLOR is checked before argparse's tty check, so setting
# it here forces plain text regardless of where the build is invoked from.
os.environ.setdefault('NO_COLOR', '1')

# Sourcerer imports
import sourcerer.Version

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.6'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.intersphinx',
              'sphinx.ext.napoleon',
              'sphinx.ext.todo',
              'sphinxcontrib.autoprogram',
              'sphinx_rtd_theme']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = 'sourcerer'
copyright = 'Kleinstein Lab, Yale University, ' + str(datetime.datetime.now().year)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = sourcerer.Version.__version__
# The full version, including alpha/beta/rc tags.
release = '%s-%s' % (sourcerer.Version.__version__, sourcerer.Version.__date__)

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# The name of the Pygments (syntax highlighting) style to use.
highlight_language = 'bash'
pygments_style = 'vs'

# If True, show todo entries
todo_include_todos = True

# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.
html_theme = 'sphinx_rtd_theme'

html_theme_options = {}

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_static/immcantation.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = ['overrides.css']

# Output file base name for HTML help builder.
htmlhelp_basename = 'sourcererdoc'


# -- Options for LaTeX output ---------------------------------------------

latex_documents = [
  ('index', 'sourcerer.tex', 'sourcerer Documentation',
   'Susanna Marquez', 'manual'),
]


# -- Options for manual page output ---------------------------------------

man_pages = [
    ('index', 'sourcerer', 'sourcerer Documentation',
     ['Susanna Marquez'], 1)
]


# -- Options for Texinfo output -------------------------------------------

texinfo_documents = [
  ('index', 'sourcerer', 'sourcerer Documentation',
   'Susanna Marquez', 'sourcerer',
   'Download data from online immune repertoire databases and format it for '
   'Immcantation.',
   'Miscellaneous'),
]

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'python': ('https://docs.python.org/3', None),
                       'changeo': ('https://changeo.readthedocs.io/en/stable', None),
                       'immcantation': ('https://immcantation.readthedocs.io/en/stable', None)}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
