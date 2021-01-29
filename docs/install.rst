Installation
============


Requirements
^^^^^^^^^^^^

GECCO requires additional libraries that can be installed directly from PyPI,
the Python Package Index. Contrary to other tools in the field
(such as DeepBGC or AntiSMASH), it does not require any external binary.


.. PyPi
.. ^^^^
..
.. GECCO is hosted on the EMBL Git server, but the easiest way to install it is
.. to download the latest release from its `PyPi repository <https://pypi.python.org/pypi/gecco>`_.
.. It will install all dependencies then install the ``gecco`` module:
..
.. .. code:: console
..
.. 	$ pip install gecco

.. Conda
.. ^^^^^
..
.. GECCO is also available as a `recipe <https://anaconda.org/bioconda/GECCO>`_
.. in the `bioconda <https://bioconda.github.io/>`_ channel. To install, simply
.. use the `conda` installer:
..
.. .. code:: console
..
.. 	 $ conda install -c bioconda GECCO
..

Git + ``pip``
^^^^^^^^^^^^^

Until GECCO is released on PyPI, you can install it from the GitHub repository
directly with ``pip``:

.. If, for any reason, you prefer to download the library from the git repository,
.. you can clone the repository and install the repository by running:

.. code:: console

	$ pip install https://github.com/zellerlab/GECCO


Keep in mind this will install the development version of the library, so not
everything may work as expected compared to a stable versioned release.


GitHub + ``setuptools``
^^^^^^^^^^^^^^^^^^^^^^^

If you do not have ``pip`` installed, you can do the following (after
having properly installed all the dependencies):

.. code:: console

	$ git clone https://git.embl.de/grp-zeller/GECCO/
	$ cd GECCO
	# python setup.py install
