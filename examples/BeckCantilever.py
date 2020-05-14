"""
This is and PARANSYS-ANSYS problem that is really quick to run, just for validation.

Example 1 from Beck, A.T. and da Rosa, E., 2006, "Structural Reliability Analysis Using Deterministic Finite Element Programs", Latin American Journal of Solids and Structures, 3, pp. 197-222.
	- https://www.lajss.org/index.php/LAJSS/article/download/101/95/
	- https://www.researchgate.net/publication/228570899_Structural_reliability_analysis_using_deterministic_finite_element_programs

Here some parameters need to be set:
 - form.ANSYS:
   - exec_loc = ansys executable location
   - run_location = where ansys will work (an empty dir)
 - form.SetANSYSModel:
   - directory = where is the file BeckCantileverAPDL.inp

"""

import paransys
import numpy as np

# Call ParAnsys
form = paransys.FORM()

form.ANSYS(exec_loc='C:\\Program Files\\ANSYS Inc\\v194\\ansys\\bin\\winx64\\ansys194.exe',
           run_location='C:\\Temp\\wk',
		   jobname='file', nproc=4, override=True, cleardir=False, add_flags='')

# Console log On
form.Info(True)

# Set ANSYS Model
form.SetANSYSModel(inputname='BeckCantileverAPDL.inp', extrafiles=[], directory='C:\\ANSYS\\Conf_Dissert\\quali_parAnsys\\_tests\\')

# Create random variables
form.CreateVar('q', 'norm', 1.15, cv=0.029)
form.CreateVar('l', 'norm', 60, cv=0.01)
form.CreateVar('b', 'norm', 4, cv=0.03)
form.CreateVar('h', 'normal', 1, cv=0.03)
form.CreateVar('sy', 'norm', 3600, cv=0.083)

# ANSYS Random variables
form.SetANSYSVar('q')
form.SetANSYSVar('l')
form.SetANSYSVar('b')
form.SetANSYSVar('h')

# Set ANSYS output
form.SetANSYSOutVar('stress')

# Create limit state
form.SetLimState('sy-stress')

# Run
values = form.Run(maxIter=12, tolRel=0.001, dh=0.005, diff='forward', meth='iHLRF')

# Exports to form69.csv
form.ExportDataCSV('BeckCantilever')
