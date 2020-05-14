# parAnsys
A Python module to do some parametric analysis on ANSYS software, it includes some reliability methods.

**Docs are available in https://dutitello.github.io/paransys/**

### ANSYS is a paid software!

# ToDo:
1) Create aliases for variables, like fy4=Normal(50, 0.05) and fy2=fy3=fy4. (same random value)
   
2) FORM Method: the final gradient could not be the correct gradient! It's because of tolRel? 
   
3) Add some R and S to FORM-ANSYS variables, this way it could know if grad is right or not (is it real? R is always positive?)
   
4) Correct this: ``Running ANSYS.
Solution is done. It took 5.424501 minutes.
Importing results from PDS results file ("C:\ANSYS\Analises\paransys-wkdir\wk1\file_current.pdrs").
iHLRF step size is 1.000000.C:\Programacao\WPy64-3720\python-3.7.2.amd64\lib\site-packages\paransys\form.py:1383: RuntimeWarning: divide by zero encountered in true_divide
  relErrorPoint = max(abs((curVecRedPts-newVecRedPts)/newVecRedPts))

Maximum relative error on design point = inf.
Absolute error betwen current and next beta = 0.5169.``

