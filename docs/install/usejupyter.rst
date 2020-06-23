
Using Jupyter notebook
======================

`Jupyter notebooks <http://jupyter-notebook.readthedocs.io/en/latest/>`__ provide an interactive way to work with Netgen/NGSolve.

If you are using Anaconda, Jupyter is already installed. Make sure to have Anaconda with Python 3.7 installed.

Using the Python package manager pip you can install Jupyter via

.. code:: bash

   # Mac/Linux
   pip3 install jupyter

On Windows you need to use pip instead of pip3. Due to a bug in ipykernel on Windows you also need to downgrade it to version 4.10 (see https://github.com/ipython/ipykernel/issues/358 ):

.. code:: bash

   # Windows
   pip install jupyter
   pip install ipykernel==4.10.0

To install the NGSolve Jupyter notebook extension, run

.. code:: bash

   jupyter nbextension install --user --py ngsolve

Now, download the first NGSolve Jupyter notebook :download:`poisson.ipynb</../py_tutorials/poisson.ipynb>`, and start

.. code:: bash

   jupyter notebook poisson.ipynb

Step though the notebook by pressing "Shift-Enter" for every cell. A separate visualization window will pop up showing the generated mesh and computed solution.


Using Jupyter lab
======================

The command needed to install the NGSolve Jupyter lab extension depends on the system and is generated with the following code:

.. code:: bash

   python3 -c 'import ngsolve.webgui; ngsolve.webgui.howtoInstallJupyterLabextension()'
