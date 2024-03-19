============
Installation
============

Installing with `conda <https://repo.anaconda.com/miniconda>`__
---------------------------------------------------------------

Starting from ``GMSO`` version ``0.3.0``, you can use `conda <https//repo.anaconda.com/miniconda>`_ to install ``GMSO`` in your preferred environment. This will also install the dependencies of ``GMSO``.

.. code-block:: bash

    (your-env) $ conda install -c conda-forge gmso


Installing from source `conda <https://repo.anaconda.com/miniconda>`__
----------------------------------------------------------------------

Dependencies of GMSO are listed in the files ``environment.yml`` (lightweight environment specification containing minimal dependencies) and ``environment-dev.yml`` (comprehensive environment specification including optional and testing packages for developers).
The ``gmso`` or ``gmso-dev`` conda environments can be created with

.. code-block:: bash

    $ git clone https://github.com/mosdef-hub/gmso.git
    $ cd gmso
    # for gmso conda environment
    $ conda env create -f environment.yml
    $ conda activate gmso

    # for gmso-dev
    $ conda env create -f environment-dev.yml
    $ conda activate gmso

    # install a non-editable version of gmso
    $ pip install .



Install an editable version from source
---------------------------------------

Once all dependencies have been installed and the ``conda`` environment has been created, the ``GMSO`` itself can be installed.

.. code-block:: bash

    $ cd gmso
    $ conda activate gmso-dev # or gmso depending on your installation
    $ pip install -e .


Supported Python Versions
-------------------------

Python 3.9-3.11 is the recommend version for users. It is the only version on which
development and testing consistently takes place.  Older (3.6-3.9) and newer (3.12+)
versions of Python 3 are likely to work but no guarantee is made and, in
addition, some dependencies may not be available for other versions.  No effort
is made to support Python 2 because it is considered obsolete as of early 2020.

Testing your installation
-------------------------

``GMSO`` uses ``py.test`` to execute its unit tests. To run them, first install the ``gmso-dev`` environment from above as well as ``gmso`` itself

.. code-block:: bash

    $ conda activate gmso-dev
    $ pip install -e .

And then run the tests with the ``py.test`` executable:

.. code-block:: bash

    $ py.test -v

Install pre-commit
------------------

We use [pre-commit](https://pre-commit.com/) to automatically handle our code formatting and this package is included in the dev environment.
With the ``gmso-dev`` conda environment active, pre-commit can be installed locally as a git hook by running

.. code-block:: bash

    $ pre-commit install

And (optional) all files can be checked by running

.. code-block:: bash

    $ pre-commit run --all-files



Building the documentation
--------------------------

``GMSO`` uses `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ to build its documentation. To build the docs locally, run the following while in the ``docs`` directory::

    $ conda env create -f docs-env.yml
    $ conda activate gmso-docs
    $ make html
