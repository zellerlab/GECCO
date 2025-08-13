Installation
============


Requirements
------------

GECCO requires additional libraries that can be installed directly from PyPI,
the Python Package Index. Contrary to other tools in the field
(such as DeepBGC or AntiSMASH), it does not require any external binary.


Installing GECCO locally
------------------------

PyPi
^^^^

The GECCO source is hosted on the EMBL Git server and mirrored on GitHub, but
the easiest way to install it is to download the latest release from
`PyPi <https://pypi.python.org/pypi/gecco>`_ with ``pip``.

.. code:: console

    $ pip install gecco


Conda
^^^^^

GECCO is also available as a `recipe <https://anaconda.org/bioconda/GECCO>`_
in the `Bioconda <https://bioconda.github.io/>`_ channel. To install, simply
use the ``conda`` installer:

.. code:: console

	 $ conda install -c bioconda GECCO


.. Git + ``pip``
.. ^^^^^^^^^^^^^
..
.. If, for any reason, you prefer to download the library from the git repository,
.. you can clone the repository and install the repository by running:
..
.. .. code:: console
..
.. 	$ pip install https://github.com/zellerlab/GECCO/archive/master.zip
..
..
.. Keep in mind this will install the development version of the library, so not
.. everything may work as expected compared to a stable versioned release.
..
..
.. GitHub + ``setuptools``
.. ^^^^^^^^^^^^^^^^^^^^^^^
..
.. If you do not have ``pip`` installed, you can do the following (after
.. having properly installed all the dependencies):
..
.. .. code:: console
..
.. 	$ git clone https://github.com/zellerlab/GECCO/
.. 	$ cd GECCO
.. 	# python setup.py install


Using GECCO in Galaxy
---------------------

GECCO is available as a Galaxy tool in the `Toolshed <https://toolshed.g2.bx.psu.edu/>`_.
It can be used directly on the `Galaxy Europe server <https://usegalaxy.eu/>`_. To
add it to your local Galaxy server, get in touch with the admin and ask them
to add the `Toolshed repository for GECCO <https://toolshed.g2.bx.psu.edu/view/althonos/gecco/88dc16b4f583>`_
to the server tools.
