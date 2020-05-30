"""
Running example 6.9 from Ang & Tang 1984
"""

import paransys
import numpy as np

# Call ParAnsys
mc = paransys.MonteCarlo()

# Console log On
mc.Info(True)

# Create random variables
mc.CreateVar('y', 'gauss', 40, cv=0.125)
mc.CreateVar('z', 'gauss', 50, cv=0.050)
mc.CreateVar('m', 'gauss', 1000, std=0.200*1000)


# An Python function that will be the Implicit function
def myFunc(y, z, m):
    # Just to show that we can do a lot of things here...

    # Some conditional
    if m > 1e27:
        m = 1e27

    # Determine the result
    result = y * z - m

    # Print the result externally to PARANSYS
    print('~~~The LS result is: {}'.format(result))

    # Send return
    return result


# There are two ways to use implicit functions:
# Setting a Python function in the place of the LS string, as next line:
mc.CreateLimState(myFunc)
# Or setting the Python function as userf and using it inside LS string:
#mc.CreateLimState("userf(y=y, z=z, m=m)", userf=myFunc)


# Define correlation
mc.SetCorrel('y', 'z', 0.40)

# Sampling for first limit state (0)
# It create failures before =)
k = 2
# For GAUSS and GUMBEL is better to use STD and for LOGN CV ;)
mc.SetRandomVarSampl('y', 0, 'gauss', 40*(1-k*0.125), std=0.125*40)
mc.SetRandomVarSampl('z', 0, 'gauss', 50*(1-k*0.050), std=0.050*50)
mc.SetRandomVarSampl('m', 0, 'gauss', 1000*(1+k*0.200), std=0.200*1000)

# Running
values = mc.Run(100, 1000, 0.05, 0.005)

# Export
#mc.ExportDataCSV('AngTang-MCAI-69', 'Comentarios?')

# Figures
#mc.Graph(['N_Pf', 'N_Beta', 'N_CVPf'], show=True)
