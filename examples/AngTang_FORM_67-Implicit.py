"""
Running example 6.7 from Ang & Tang 1984 with an Implicit limit state function
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
form.SetLimState(myFunc)
# Or setting the Python function as userf and using it inside LS string:
#form.SetLimState("userf(y=y, z=z, m=m)", userf=myFunc)

# Run
values = form.Run(dh=0.01, meth='iHLRF')
