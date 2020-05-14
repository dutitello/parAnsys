"""
Running example 6.8 from Ang & Tang 1984
"""

import paransys
import numpy as np

# Call ParAnsys
mc = paransys.MonteCarlo()

# Console log On
mc.Info(True)

# Create random variables
mc.CreateVar('y', 'logn', 40, cv=0.125)
mc.CreateVar('z', 'logn', 50, cv=0.050)
mc.CreateVar('m', 'gumbel', 1000, cv=0.200)

# First limit state (0)
mc.CreateLimState('y*z-m')

# Sampling for first limit state (0)
# It create failures before =)
k = 2
# For GAUSS and GUMBEL is better to use ORIGINAL STD and for LOGN CV ;)
mc.SetRandomVarSampl('y', 0, 'logn', 40*(1-k*0.125), cv=0.125)
mc.SetRandomVarSampl('z', 0, 'logn', 50*(1-k*0.050), cv=0.050)
mc.SetRandomVarSampl('m', 0, 'gumbel', 1000*(1+k*0.200), std=0.200*1000)

# Running
values = mc.Run(100, 1000, 0.05, 0.005)

# Export
mc.ExportDataCSV('AngTang-MCAI-68', 'Comentarios?')

# Figures
mc.Graph(['N_Pf', 'N_Beta', 'N_CVPf'], show=True)
