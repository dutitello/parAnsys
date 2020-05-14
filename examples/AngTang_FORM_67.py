"""
Running example 6.7 from Ang & Tang 1984
"""

import paransys

# Call ParAnsys
form = paransys.FORM()

# Console log On
form.Info(True)

# Create random variables
form.CreateVar('y', 'gauss', 40, cv=0.125)
form.CreateVar('z', 'gauss', 50, cv=0.050)
form.CreateVar('m', 'gauss', 1000, cv=0.200)

# Create limit state
form.SetLimState('y*z-m')

# Run
values = form.Run(dh=0.01, meth='iHLRF')
form.ExportDataCSV('AngTang-FORM-67', 'Comentarios?')
