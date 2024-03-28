Using GMSO with Docker
========================

As much of scientific software development happens in unix platforms, to avoid the quirks of development dependent on system you use, a recommended way is to use docker or other containerization technologies. This section is a how to guide on using ``GMSO`` with docker.

Prerequisites
-------------
A docker installation in your machine. Follow this `link <https://docs.docker.com/get-docker/>`_ to get a docker installation working on your machine. If you are not familiar with docker and want to get started with docker, the Internet is full of good tutorials like the ones `here <https://docker-curriculum.com/>`__ and `here <https://www.youtube.com/watch?v=zJ6WbK9zFpI&feature=youtu.be>`__.

Quick Start
-----------
After you have a working docker installation, please use the following command to use run a jupyter-notebook with all the dependencies for `GMSO` installed:

.. code-block:: bash

    $ docker pull mosdef/gmso:latest
    $ docker run -it --name gmso -p 8888:8888 mosdef/gmso:latest

If no command is provided to the container (as above), the container starts a ``jupyter-notebook`` at the (container) location
``/home/anaconda/data``.
To access the notebook, paste the notebook URL into a web browser on your computer.  When you are finished, you can control-C to
exit the notebook as usual.  The docker container will exit upon notebook shutdown.


Alternatively, you can also start a Bourne shell to use python from the container's terminal:

.. code-block:: bash

    $ docker run -it --mount type=bind,source=$(pwd),target="/home/anaconda/data" mosdef/gmso:latest sh
    ~ $ source .profile
    (gmso-dev) ~ $

.. Warning::
    Containers by nature are ephemeral, so filesystem changes (e.g., adding a new notebook) only persist until the end of the
    container's lifecycle.  If the container is removed, any changes or code additions will not persist.  See the section below for
    persistent data.

.. note::

    The -it flags connect your keyboard to the terminal running in the container.
    You may run the prior command without those flags, but be aware that the container will not respond to any keyboard input.
    In that case, you would need to use the docker ``ps`` and ``docker kill`` commands to shut down the container.


Persisting User Volumes
-----------------------
If you will be using `GMSO` from a docker container, a recommended way is to mount what are called user volumes in the container. User volumes will provide a way to persist all filesystem/code additions made to a container regardless of the container lifecycle. For example, you might want to create a directory called `gmso-notebooks` in your local system, which will store all your `GMSO` notebooks/code. In order to make that accessible to the container(where the notebooks will be created/edited), use the following steps:

.. code-block:: bash

    $ mkdir -p /path/to/gmso-notebooks
    $ cd /path/to/gmso-notebooks
    $ docker run -it --name gmso --mount type=bind,source=$(pwd),target=/home/anaconda/data -p 8888:8888 mosdef/gmso:latest

You can easily mount a different directory from your local machine by changing ``source=$(pwd)`` to
``source=/path/to/my/favorite/directory``.

.. note::

    The ``--mount`` flag mounts a volume into the docker container. Here we
    use a ``bind`` mount to bind the current directory on our local filesystem
    to the ``/home/anaconda/data`` location in the container. The files you see
    in the ``jupyter-notebook`` browser window are those that exist on your
    local machine.

.. warning::

    If you are using the container with jupyter notebooks you should use
    the ``/home/anaconda/data`` location as the mount point inside the container;
    this is the default notebook directory.

Running Python scripts in the container
---------------------------------------
Jupyter notebooks are a great way to explore new software and prototype
code. However, when it comes time for production science, it is often
better to work with python scripts. In order to execute a python script
(``example.py``) that exists in the current working directory of your
local machine, run:

.. code-block:: bash

    $ docker run --mount type=bind,source=$(pwd),target=/home/anaconda/data mosdef/gmso:latest "python data/test.py"

Note that once again we are ``bind`` mounting the current working directory
to ``/home/anaconda/data``. The command we pass to the container is
``python data/test.py``. Note the prefix ``data/`` to the script; this is because
we enter the container in the home folder (``/home/anaconda``), but our script
is located under ``/home/anaconda/data``.

.. warning::

    Do not bind mount to ``target=/home/anaconda``. This will cause errors.


If you don't require a Jupyter notebook, but just want a Python interpreter,
you can run:

.. code-block:: bash

    $ docker run --mount type=bind,source=$(pwd),target=/home/anaconda/data mosdef/gmso:latest python

If you don't need access to any local data, you can of course drop the
``--mount`` command:

.. code-block:: bash

    $ docker run mosdef/gmso:latest python

Cleaning Up
-----------
You can remove the created container by using the following command:

.. code-block:: bash

    $ docker container rm gmso

.. note::

    Instead of using `latest`, you can use the image :code:`mosdef/gmso:stable` for most recent stable release of ``GMSO`` and run the tutorials.
