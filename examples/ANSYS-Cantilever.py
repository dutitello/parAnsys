"""
In this example we just connect ANSYS to Python by PARANSYS.

Here some parameters need to be set:
 - form.ANSYS:
   - exec_loc = ansys executable location
   - run_location = where ansys will work (an empty dir)
 - form.SetANSYSModel:
   - directory = where is the file BeckCantileverAPDL.inp

"""

import paransys
import numpy as np
import pandas as pd

#===========================================================================================
# Analysis setup
#
# Call ANSYS
ans = paransys.ANSYS(exec_loc='C:\\Program Files\\ANSYS Inc\\v194\\ansys\\bin\\winx64\\ansys194.exe',
                    run_location='C:\\Temp\\wk',
                    jobname='file', nproc=4, override=True, cleardir=False, add_flags='')

# Activate the output
ans.Info(True)

# Set APDL model
ans.SetModel(inputname='BeckCantileverAPDL.inp', extrafiles=[], directory='.')

# Set input parameters
ans.CreateVarIn('q')
ans.CreateVarIn('l')
ans.CreateVarIn('b')
ans.CreateVarIn('h')

# Set output parameters
ans.CreateVarOut('stress')

# Set analysis lenght (number of tests)
ans.SetLength(3)

#===========================================================================================
# Set values and Run! If the model need to run many times just this part need to change!
#
# We are doing 3 analysis using 4 variables, so we need to create 4 NumPy arrays with the 
# values that we will be testing at each simulation, like in this table:
#  ____________________________________________
# |  Sim\Var |   q   |    l   |   b   |   h   |
# |  Sim #1  |  1.15 |  60.0  |  4.0  |  1.0  |
# |  Sim #2  |  1.15 |  62.0  |  4.0  |  1.0  |
# |  Sim #3  |  1.15 |  60.0  |  4.0  |  1.2  |
# |__________|_______|________|_______|_______|
#
# Values at (Sim#1, Sim#2, Sim#3)
q = np.array([1.15,  1.15,  1.15])
l = np.array([60.0,  62.0,  60.0])
b = np.array([4.00,  4.00,  4.00])
h = np.array([1.00,  1.00,  1.20])

# Passing it to ANSYS
ans.SetVarInValues('q', q)
ans.SetVarInValues('l', l)
ans.SetVarInValues('b', b)
ans.SetVarInValues('h', h)

# Run baby, RUN!
ans.Run()

# Now we get the results in a dictionary
results = ans.GetVarOutValues()

# Print all results
print(results)
# The output parameters always came in UPPER CASE because of ANSYS is made in FORTRAN77 (I think it's because of that!)

#
# NOW IT'S JUST FUN!
# Print just stresses
print(results['STRESS'])

# Put all in a Pandas DataFrame:
df = pd.DataFrame({'q': q, 'l': l, 'b': b, 'h': h, 'stress': results['STRESS']})
print(df)

# Derivatives of dStress/dl and dStress/dh:
dSdl = (df['stress'][1]-df['stress'][0])/(df['l'][1]-df['l'][0])
print(f'dStress/dl = {dSdl}')
dSdh = (df['stress'][2]-df['stress'][0])/(df['h'][2]-df['h'][0])
print(f'dStress/dh = {dSdh}')


#===========================================================================================
# Bonus:
#   If you have to change variables names, APDL model or anything beside the values
#   you may have to use:
# 
# To clear everything after paransys.ANSYS():
#   It will erase the model, variables names and values.
#ans.ClearAll()
#
# To clear just the variables values:
#ans.ClearValues()
#
#===========================================================================================
