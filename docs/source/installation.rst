Installation
============

Requirements
------------

* Python 3.8+
* numpy
* astropy
* scipy
* matplotlib

Using pip
---------

The package is available on PyPI and can be installed with pip:

.. code-block:: bash

    pip install a3cosmos-gas-evolution

Using conda
-----------

You can also install the package using conda with the following command:

.. code-block:: bash

    conda install -c conda-forge a3cosmos-gas-evolution

From source
-----------

To install from source, clone the repository and run:

.. code-block:: bash

    git clone https://github.com/1054/a3cosmos-gas-evolution.git
    cd a3cosmos-gas-evolution
    pip install -e .

Development installation
------------------------

To install in development mode with all dev dependencies:

.. code-block:: bash

    git clone https://github.com/1054/a3cosmos-gas-evolution.git
    cd a3cosmos-gas-evolution
    pip install -e ".[dev]"
    pip install pytest
