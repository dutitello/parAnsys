Welcome to PARANSYS's documentation!
====================================
**PARANSYS**: *Python pArametric Reliability Analysis on ANSYS*, but it's not limited just for that!

PARANSYS is a module that can connect the ANSYS software to Python scripts using APDL 
scripts. It has one connection class, that can be used with many types of APDL analysis 
(more than reliability analysis), and two other classes for reliability analysis, 
that can be used with implicit and explicit Python functions, not just with ANSYS.

Introduction
------------

This module calls ANSYS by Probabilistic Design System (PDS) as an FEM tool for
running a couple of parametric analysis and get some variables back to Python.

The class **paransys.ANSYS** is the heart of the process, this class receive the
ANSYS APDL model, the parameters names and values, runs everything and return
the control parameters values for each simulation back to Python. 
This class makes easy to run a lot of analysis with the same model just changing
some parameters values, in this scope we have Reliability Analysis,
Optimizations and Mesh Convergence Studies, for example.

For Reliability Analysis there are two classes implemented:

* **paransys.MonteCarlo**

  This class performs Monte Carlo reliability analyses, with and without using
  importance sampling. When using importance sampling it's possible to fix the
  sampling point or search it based on the center of failures weights. *It's 
  currently limited to one limit state.*


* **paransys.FORM**

  This class applies the First Order Reliability Method (FORM) inside Python
  using ParAnsys as a connection with ANSYS for evaluate FEM models.
  It is also possible to run simulations without using ANSYS, just
  defining the limit state equation and all variables.


How to install
--------------

PARANSYS is available at 
`github.com/dutitello/parAnsys <https://github.com/dutitello/parAnsys>`_.
The installer isn't done yet, but you can just copy the `paransys` dir inside
the repository to your Python modules folder, or have it in your working directory. 


How to cite
-----------

This module is published at `Zenodo <https://github.com/dutitello/parAnsys>`_. 

You can cite it as follows, or you can get another format in it's Zenodo page. 
    
    Titello, E.P., Campos Filho, A., & Real, M.V. (2020). 
    PARANSYS: Python pArametric Reliability Analysis on ANSYS. Version 1.0. 
    Zenodo. http://doi.org/10.5281/zenodo.xxxxxxxx.



Final notes
-----------
This is part of my master degree thesis, available `here <https://lume.ufrgs.br/handle/10183/213006>`_.

**This module makes no claim to own any rights to ANSYS, it's just an interface
to call the program owned by ANSYS company.**


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   ANSYS
   MonteCarlo
   FORM
   Examples


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
