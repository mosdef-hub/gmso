============
Installation
============

Installing dependencies with `conda <http://continuum.io/downloads>`_
---------------------------------------------------------------------

Dependencies of ``Topology`` are listed in the file ``requirements.txt``. They
can be installed in one line:
::

    $ conda install -c omnia -c mosdef -c conda-forge --file requirements.txt

Alternatively you can add all the required channels to your ``.condarc`` file
and then install dependenices.

    $ conda config --add channels omnia
    $ conda config --add channels mosdef
    $ conda config --add channels conda-forge
    $ conda install --file requirements.txt

.. note::
    These commands will likely change a configuration file on your computer and
    may affect installation of other packages in other projects you are working
    on. However, the channel priority recommended is fairly common
    (in particle, having ``conda-forge`` having the highest priority) and should
    work well for most installations.

Installing dependencies with `pip <https://pypi.org/project/pip/>`_
-------------------------------------------------------------------
::

    $ pip install --file requirements.txt

.. note::
    Compared to ``conda`` installation, this is less tested. Some upstream
    dependencies may not be available on ``PyPI`` but can be installed via
    source or ``conda``.

Install an editable version from source
---------------------------------------

Once all dependencies are installed, the ``Topology`` itself can be installed.
It is currently only available through its source code. It will be available
through ``pip`` and ``conda`` in the future.
::

    $ git clone https://github.com/mattwthompson/topology
    $ cd topology
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

``Topology`` uses ``py.test`` to execute its unit tests. To run them, first install some extra depdencies:
::
    $ conda install --file requirements-test.txt


And then run the tests with the ``py.test`` executable:
::
    $ py.test -v
