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

# A clean Implicit limit state function
def myFunc(y, z, m):
    return y * z - m

# But it's a function, so the parameters doesn't need to have the same name of random variables and we can have some fun, like this. 
#   Of course, here I'm trying keep the result, but you can see the idea... 
#   It's possible to connect complex functions here or create a function that do the connection (calling a lot of functions inside it)
def CrazyFunction(x1, x2, x3):
    if(x1 < -1000): x1 = 1000
    if(x2 >= 1e28): x2 = 0
    res = x1 * x2 - x3
    if(res <= -100000): res = -100000

    return res


# Create limit state
#   You can hange the function name at userf=myFunc/CrazyFunction, but it's always called as userf() inside the LS
#       Please, see that myFunc() uses x,y,m and CrazyFunction() uses x1,x2,x3, but inside userf() is used the random variables names
form.SetLimState('userf(x3=m,x2=z,x1=y)', userf=CrazyFunction)

# Run
values = form.Run(dh=0.01, meth='iHLRF')
