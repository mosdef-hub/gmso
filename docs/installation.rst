============
Installation
============

Installing with `conda <https://repo.anaconda.com/miniconda>`_
--------------------------------------------------------

Starting from ``GMSO`` version ``0.3.0``, you can use `conda <https//repo.anaconda.com/miniconda>`_ to install ``GMSO`` in your preferred environment. This will also install the dependencies of ``GMSO``.
::
    (your-env) $ conda install -c conda-forge gmso


Installing dependencies with `conda <https://repo.anaconda.com/miniconda>`_
---------------------------------------------------------------------

Dependencies of ``GMSO`` are listed in the file ``environment.yml``. They
can be installed in a ``conda`` environment named ``gmso`` with::

    $ conda env create -f environment.yml
    $ conda activate gmso


Install an editable version from source
---------------------------------------

Once all dependencies are installed, the ``GMSO`` itself can be installed.
::

    $ git clone https://github.com/mosdef-hub/gmso.git
    $ cd gmso
    $ conda env create -f environment-dev.yml
    $ conda activate gmso-dev
    $ pip install -e .


Supported Python Versions
-------------------------

Python 3.7 is the recommend version for users. It is the only version on which
development and testing consistently takes place.  Older (3.6) and newer (3.8+)
versions of Python 3 are likely to work but no guarantee is made and, in
addition, some dependencies may not be available for other versions.  No effort
is made to support Python 2 because it is considered obsolete as of early 2020. 

Testing your installation
-------------------------

``GMSO`` uses ``py.test`` to execute its unit tests. To run them, first install the dev environment
::

    $ git clone https://github.com/mosdef-hub/gmso.git
    $ cd gmso
    $ conda env create -f environment-dev.yml
    $ conda activate gmso-dev
    $ pip install -e .

And then run the tests with the ``py.test`` executable:
::
    $ py.test -v


Building the documentation
--------------------------

``GMSO`` uses `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ to build its documentation. To build the docs locally, run the following while in the ``docs`` directory::
    
    $ pip install -r requirements.txt
    $ make html

