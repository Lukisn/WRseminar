How to build this documentation
===============================

This documentation is created and built using sphinx.
The initial structure was initialized using ``sphinx-quickstart``.
Further tweaking can be done in the ``conf.py`` file.

The automatic api documentation files in the ``apidoc/`` folder were created
using the ``sphinx-apidoc`` command:

.. code-block:: bash

    .../doc/$ sphinx-apidoc -o apidoc/ ..

To build the actual documentation sphinx offers an easy entry point using by
offering a Makefile.
So building the documentation is as simple as using the ``make`` command, for
example the html docs are created by:

.. code-block:: bash

    .../doc/$ make html

This is all there is about it.

The plain text of the documentation is organized within the restructuredText
files in the ``docs/`` folder.