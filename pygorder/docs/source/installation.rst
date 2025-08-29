Installation
============

Using ``uv``
++++++++++++

The recommended way to install ``gorder`` as a Python package is to use the
`uv <https://github.com/astral-sh/uv>`_ package manager. To add ``gorder`` to your
``uv`` project, run:

.. code-block:: bash

    uv add git+https://github.com/Ladme/gorder.git#subdirectory=pygorder

Using ``conda``
+++++++++++++++

``gorder`` is also available on
`conda-forge <https://anaconda.org/conda-forge/pygorder>`_
(credit to `@RubenChM <https://github.com/RubenChM>`_):

.. code-block:: bash

    conda create -n gorder
    conda activate gorder
    conda install conda-forge::pygorder

Using ``pip``
+++++++++++++

If you use neither ``uv`` nor ``conda``, you can also install ``gorder`` using ``pip``:

.. code-block:: bash

    pip install git+https://github.com/Ladme/gorder.git#subdirectory=pygorder