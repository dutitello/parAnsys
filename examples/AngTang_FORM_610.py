"""
Running example 6.10 from Ang & Tang 1984
"""

import paransys

# Call ParAnsys
form = paransys.FORM()

# Console log On
form.Info(True)

# Create random variables
form.CreateVar('y', 'logn', 40, cv=0.125)
form.CreateVar('z', 'logn', 50, cv=0.050)
form.CreateVar('m', 'gumbel', 1000, cv=0.200)

# Set correlation
form.SetCorrel('y', 'z', 0.40)

# Create limit state
form.SetLimState('y*z-m')

# Run
values = form.Run(dh=0.001, meth='iHLRF', tolRel=1E-6)
form.ExportDataCSV('AngTang-FORM-610', 'Comentarios?')
