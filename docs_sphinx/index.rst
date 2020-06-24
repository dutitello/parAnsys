.. parAnsys documentation master file, created by
   sphinx-quickstart on Tue May 14 22:08:02 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PARANSYS's documentation!
====================================

This is part of my master degree thesis, which is currently under development.
------------------------------------------------------------------------------

This module calls ANSYS by Probabilistic Design System (PDS) as an FEM tool for
running a couple of parametric analysis and get some variables back to Python.

The class **paransys.ANSYS** is the heart of the process, this class receive the
ANSYS APDL model, the parameters names and values, runs everything and return
the control parameters values for each simulation back to Python.

This class makes easy to run a lot of analysis with the same model just changing
some parameters values, in this scope we have Reliability Analysis,
Optimizations and Mesh Convergence Studies, for example.


For Reliability Analysis we have two others classes implemented using **paransys.ANSYS**:

* **paransys.MonteCarlo**

  This class performs Monte Carlo reliability analyses, with and without using
  importance sampling. When using importance sampling it's possible to fix the
  sampling point or search it based on the center of failures weights.
  It's also possible to run Monte Carlo simulations for a limit state equation,
  without using ANSYS.

  **For now this class just accept one limit state function/failure mode.**

  |

* **paransys.FORM**

  This class applies the First Order Reliability Method (FORM) inside Python
  using ParAnsys as a connection with ANSYS for evaluate FEM models.
  It is possible to run simulations without using ANSYS, just
  defining the limit state equation and all variables.

  |

This module makes no claim to own any rights to ANSYS, it's just an interface
to call the program owned by ANSYS.


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   ANSYS
   MonteCarlo
   FORM
   Examples
   Tutorial


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
