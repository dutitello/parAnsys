## -*- coding: UTF-8 -*-
"""
This module performns Reliability Analysis with First Order Reliability Method using Python.
Please read the class docstring for more. 

Docs are available at https://dutitello.github.io/parAnsys/ 
"""

import os
import time
from paransys.ansys import ANSYS
import numpy as np
import scipy.stats
from math import * # It's necessry to evaluate limit states
import math
from inspect import isfunction


class FORM(object):
	"""
	This class applies the First Order Reliability Method (FORM) inside Python
	using ParAnsys as a connection with ANSYS for evaluate FEM models.

	It is possible to run simulations without using ANSYS, just
	defining the limit state equation and all variables.

	This code was made following the ideia of ANSYS being a tool for getting the
	ultimate load of the structure. This works applying a displacement in the
	loaded node, and then getting the biggest reaction force on that node,
	following this way the limit state defined here is 'R-S', where R are the
	values get from ANSYS and S the values generated in Python. It's also
	possible to work applying the true load on ANSYS, it's just necessary to
	formulate a valid limit state equation.

		
	|
	|

	**Class methods:**

	
	"""

	#---------------------------------------------------------------------------



	def __init__(self):
		"""

		"""

		# Before
		self._ANSYS = False
		self.PrintR = True
		self.limstate = None
		self._userf = None
		self.variableDistrib = {}
		self.variableConst = {}
		self.variableStartPt = {}
		self.controls = {}
		self.corlist = {}
		self._last_ck = 0

		# Options
		self._options = {}
		self._options['iHLRF_forced_lambdk'] = 'auto'
		self._options['iHLRF_prod_ck'] = 2.00
		self._options['iHLRF_add_ck'] = 0.00
		self._options['iHLRF_par_a'] = 0.10
		self._options['iHLRF_par_b'] = 0.50
		self._options['iHLRF_step_lambdk_test'] = 4
		self._options['rHLRF_relax'] = 0.50
		self._options['APDLdebug'] = False

		# After
		self.variableDesPt = {}
		self.results = {}
		self._stnumb = 99


	def ANSYS(self, exec_loc=None, run_location=os.getcwd()+'\\ansys_anl\\', jobname='file',
			     nproc=2, override=False, cleardir=False, add_flags=''):
		"""
		If ANSYS will be used it defines ANSYS properties, for initialize the
		paransys.ANSYS class.

		Parameters
		----------
		exec_loc : str, obligatory
			Location of ANSYS executable file.

		run_location : str, optional
			ANSYS working directory. Recomended to be a separated directory.
			Defaults to ansys_anl on current directory.

		jobname : str, optional
			ANSYS jobname. Defaults to 'file'.

		nproc : int, optional
			Number of processors. Defaults to 2.

		override : bool, optional
			Attempts to delete the .lock file at working directory.
			It's useful when ANSYS was interrupted.
			Defaults to False

		cleardir : bool, optional
			Delete all the files from ANSYS working directory when call the Run command.
			Defaults to False

		add_flags : str, optional
			Additional flags to be called with ANSYS.
			If it's an academic version use add_flags='-aa_r'
			Do not use '-b -i -o'
			Flags can be found at https://www.sharcnet.ca/Software/Ansys/16.2.3/en-us/help/ans_ope/Hlp_G_OPE3_1.html

		"""
		self.ansys = ANSYS(exec_loc=exec_loc, run_location=run_location,
								jobname=jobname, nproc=nproc, override=override,
								cleardir=cleardir, add_flags=add_flags)

		self._ANSYS = True

	def Info(self, act=False):
		"""
		Turn on/off the return of the commands to Python.

		Parameters
		----------
		act : bool, obligatory
			True turn On and False turn Off the return of the commands to Python.
		"""
		self.PrintR = act
		if self._ANSYS:
			self.ansys.Info(act)
		#return self._PrintR('Now the commands will send a return to Python (like this).')


	def _PrintR(self, value):
		"""
		Internal function to print or not the return of commands based on command Info()
		"""
		if self.PrintR:
			return print(value, flush=True)
		else:
			pass

	def SetANSYSModel(self, inputname, extrafiles=[], directory=os.getcwd()):
		"""
		Set the input script file to be used on ANSYS and extra files that should
		be copied together.
		All this files must be in the same directory set in parameter directory.

		Parameters
		----------
		inputname : str, obligatory
			Name with extension of the script that will be executed in the analysis.
			The script must be done in function of the INPUT variables defined here,
			(as parameters of ANSYS), and must define/calculate ANSYS parameters with
			the results using the names defined here.

		extrafiles : list of strings, optional
			A list of strings containing extra files (with extension) that are necessary to
			run the script analys, could be an MODEL with the MESH already generated,
			for example.
			An example of extrafiles list is:
			``extrafiles = ['temps.txt', 'model1.ans', 'file.db']``

		directory : str, optional
			If the script is not in the current running Python directory you should
			place the entire location, if it's in a subdirectory of current directory
			you can use '/dirname/filename.ext'.
			Defaults to current running Python directory.
		"""
		if self._ANSYS:
			self.ansys.SetModel(inputname, extrafiles, directory)
		else:
			exception = Exception('ANSYS not declared yet. Before set ANSYS '+
				'model you must define ANSYS properties with ANSYS(...).')
			raise exception


	def SetANSYSOutVar(self, name):
		"""
		Defines a parameter/variable from ANSYS APDL script as an variable to
		return values for Python.

		Parameters
		----------
		name : str, obligatory
			Variable/Parameter name, as defined in APDL script.
		"""
		if self._ANSYS:
			self.ansys.CreateVarOut(name)
		else:
			exception = Exception('ANSYS not declared yet. Before set ANSYS '+
				'variables you must define ANSYS properties with ANSYS(...).')
			raise exception


	# Setting distribution of variables
	def CreateVar(self, name, distrib, mean, std=0, cv=None, par1=None, par2=None):
		"""
		Create a Variable, random or not.

		If it's used on ANSYS it need to be told, so after this use:

		>>> form.SetANSYSVar(name)

		Parameters
		----------
		name : str, obligatory
			Name of variable.

		distrib : str, obligatory
			Probabilistic variable distribution type.

			For all distributions Mean and Std are related to Normal distribution
			(the code determines the parameters for the desired distribution).

			Available types are:
			* gaussian (or gauss, normal);
			* lognormal (or log, logn, ln, lognorm);
			* gumbel (or gumb, type1);
			* constant (or const) - Constant value (doesn't need std).

		mean : float, obligatory
			Standard mean of variable values.

		std : float, optional
			Standard deviation of variable. You must define it or cv for variables
			that aren't constant, if both (cv and std) declared std will be used.

			For LogNormal variables it's recommend to use CV!


		cv : float, optional
			Coeficient of Variation of variable. You must define it or std for variables
			that aren't constant, if both (cv and std) declared std will be used.

			For LogNormal variables it's recommend to use CV!


		par1 and par2 : float, optional
			Parameters for future implementations.

		"""

		# Verify the variable name
		if name in ['', None, ' ']:
			exception = Exception('You must define the name of the variable that will receive the values.')
			raise exception

		# Set distribution name as lower case
		name = name.lower()
		distrib = distrib.lower()

		# CV or STD?
		if distrib not in ['constant', 'const', 'cons', 'c']:
			if cv == None:
				cv = std/mean
			else:
				std = cv*mean
		

		# Verify the distribution and then store the parameters
		# Gaussian variable
		if distrib in ['gauss', 'gaus', 'gaussian', 'normal', 'norm']:
			# Store values
			self.variableDistrib[name] = ['gauss', mean, std, cv]
			# Set the start point as the mean
			self.variableStartPt[name] = mean
			return self._PrintR('Variable \"%s\" defined as Gaussian with mean=%f, std. dev.=%f, CV=%f.'
									% (name, mean, std, cv))

		# Lognormal distribution
		elif distrib in ['lognormal', 'logn', 'ln', 'log', 'lognorm']:
			# Store values
			self.variableDistrib[name] = ['logn', mean, std, cv]
			# Set the start point as the mean
			self.variableStartPt[name] = mean
			return self._PrintR('Variable \"%s\" defined as LogNormal with mean=%f, std. dev.=%f, CV=%f.'
									% (name, mean, std, cv))

		# Gumbel distribution
		elif distrib in ['gumbel', 'gumb', 'type1']:
			# Store values
			self.variableDistrib[name] = ['gumbel', mean, std, cv]
			# Set the start point as the mean
			self.variableStartPt[name] = mean
			return self._PrintR('Variable \"%s\" defined as Gumbel with mean=%f, std. dev.=%f, CV=%f.'
									% (name, mean, std, cv))

		# Constant value
		elif distrib in ['constant', 'const', 'cons', 'c']:
			# Store value
			self.variableConst[name] = mean
			return self._PrintR('Variable \"%s\" defined as Constant with value=%f'
									% (name, mean))

		# or what?
		else:
			exception = Exception('Distribution \"%s\" set on variable \"%s\" is not recognized.'
									% (distrib.upper(), name))
			raise exception


	def _SetControls(self, maxIter, tolRel, tolLS, dh, diff, meth):
		"""
		Internal from Run()
		Set the controls of process.

		Parameters
		----------
		maxIter : integer
			Maximum of iterations that can be performed. After this the process
			will stop with error.

		tolRel : float


		tolLS : float or string


		dh : float
			delta_h step when applying derivatives, value applied in reduced
			space.

		diff : str


		meth : str
			FORM method used. Available methods are:

				* HLRF: Hasofer Lind Rackwitz and Fiessler method.

				* iHLRF: improved Hasofer Lind Rackwitz and Fiessler method.
		"""

		if diff not in ['center', 'forward', 'backward']:
			exception = Exception('Invalid derivative method.')
			raise exception

		if min(maxIter, tolRel, dh) >= 0 and meth in ['HLRF', 'iHLRF', 'rHLRF'] and \
			diff in ['center', 'forward', 'backward']:

			# Save controls variable
			self.controls['maxIter'] = maxIter
			self.controls['tolRel'] = tolRel
			self.controls['tolLS'] = tolLS
			self.controls['dh'] = dh
			self.controls['diff'] = diff
			self.controls['meth'] = meth


			self._PrintR('Process controls set as:')
			self._PrintR('   Limit of iterations: %d.' % maxIter)
			self._PrintR('   Relative error tolerance: %2.3E.' % tolRel)
			if isinstance(tolLS, str):
				self._PrintR('   Absolute LS error tolerance: auto.')
			else:
				self._PrintR('   Absolute LS error tolerance: %2.3E.' % tolLS)
			self._PrintR('   deltah (for derivatives): %2.3E.' % dh)
			self._PrintR('   Finite difference method: %s.' % diff)
			self._PrintR('   FORM Method: %s.' % meth)

		else:
			exception = Exception('Error while setting simulation controls. Please verify the set values.')
			raise exception


	def SetLimState(self, equat, userf=None):
		"""
		Set the limit state equation.

		Parameters
		----------
		equat : obligatory
			1) It could be a string with the equation of the limit state. It must be write as a
			function of defined variables (In and Out).
			2) It could be a Python function created by the user, just passing the function name in 
			the place of the string.

		userf : function, optional
			An user defined function that could be used inside the limit state
			equation (string), called inside equat as ``userf()``. 

			It's similar to use a function instead of a string in the equat parameter.

			For example, you can create a complex Python function with loops, ifs
			and whatever for evaluate the R part of your limit state function
			for a concrete beam. An example is showed after.

		First example: if ANSYS returns the maximum load on a truss as variable
		FxMAX, and applied loads to be tested are ``(g+q)*sin(theta)``, where
		``g``, ``q``, theta are defined random variables created with ``CreateVar()``.

		.. code-block:: python

			form.SetLimState(equat='FxMAX-(g+q)*sin(theta)')

		Note that you can use math expressions as ``sin()``, ``cos()``, ``tan()``, ``sqrt()``
		from Python math module inside the equation.

		|

		Second example: you have a steel bar in tension that hasn't hardening.
		It's stress is a function of ``(eps, fy, E)``, where ``eps`` is current
		deformation, ``fy`` is yield stress and ``E`` the elastic moduli,
		you can create inside your code	an function like:

		.. code-block:: python

			def stress(eps, fy, E):
				if eps > fy/E:
					return fy
				else:
					return eps*E

		And now defining ``userf=stress`` we can:

		.. code-block:: python

			form.SetLimState(equat='userf(eps,fy,E)-q', userf=stress)

		where ``eps``, ``fy``, ``E`` and ``q`` are random variables.
		Note that the function inside the limit state equation should be
		called as ``userf()`` with the parameters from ``stress``.

		Or we can do the same using the functions instead of the string:

		.. code-block:: python

			form.SetLimState(equat=stress)

		"""
		# Store it
		self.limstate = equat
		
		if(type(self.limstate) == str):
			# String equation
			# Change equation to lowcase
			self.limstate = self.limstate.lower()
			self._userf = userf

			return self._PrintR('Limit state defined as "{}".'.format(equat))
		else:
			# LS is a Python Function 
			return self._PrintR('Limit state defined as "{}" function.'.format(self.limstate.__name__))

	def SetANSYSVar(self, name):
		"""
		Set a variable as ANSYS variable.

		Parameters
		----------
		name : str, obligatory
			Name of variable.

		"""

		# Verify ANSYS object
		if self._ANSYS == False:
			exception = Exception('ANSYS not declared yet. Before set ANSYS '+
				'variables you must define ANSYS properties with ANSYS(...).')
			raise exception


		# Verify the variable name
		if name in ['', None, ' ']:
			exception = Exception('You must define the name of the variable.')
			raise exception
		else:
			name = name.lower()
			if name not in self.variableDistrib and name not in self.variableConst:
				exception = Exception('This variable name is not declared yet. '+
						'Only Random variables can be set as ANSYS variables.\n'+
						'Please use CreateVar() to declare it.')
				raise exception

		# Declare it on ansys object
		self.ansys.CreateVarIn(name)


	def SetStartPoint(self, name, value):
		"""
		Set the point, for each variable, that process will start.
		If it's not declared it will start with the mean value.

		Parameters
		----------
		name : str, obligatory
			Variable name.

		value : float, obligatory
			Starting point for this variable.

		"""
		# Set lower
		name = name.lower()
		# Verify if it's already created
		if name not in self.variableDistrib:
			exception = Exception('This variable name is not declared yet. '+
					'Before set the start value you must create the random '+
					'variable with CreateVar().')
			raise exception

		# If no error: store the value
		self.variableStartPt[name] = value

		#
		self._PrintR('Variable \"%s\" will start the process at point %f.' % (name, value))


	def Options(self, option, value=None):
		"""
		Set extra options values.


		Parameters
		----------
		option : str, obligatory
			Name of option, listed next.

		value : optional
			Value to be set, type varies with option.
			If not defined it will return current value.


		** Valid options:**
		For iHLRF method:
			* iHLRF_forced_lambdk : float
			  Forced value when line search doesnt found a valid ``lambdak``.
			  Being ``lambdak`` the step size.

			  It could be set as ``'auto'``, when it
			  is the complement of ``|y*.gradG|/(|y*|.|gradG|)``.
			  Defaults to 'auto'.

			* iHLRF_prod_ck : float
			  Scalar value that will be multiplied by calculated ``ck`` value. For a
			  fix ``ck`` value turn it to 0 and then use 'iHLRF_add_ck'.

			* iHLRF_add_ck : float
			  Scalar value that will be added to ``ck`` value.

			* iHLRF_par_a : float
			  Value presented as ``a`` in line search equation for iHLRF.

			* iHLRF_par_b : float
			  Value presented as ``b`` in line search equation for iHLRF,
			  ``lambdak`` value is ``b**nk``.

			* iHLRF_step_lambdk_test : float
			  Size of ``lambdak`` test block, after each block convergence is checked.

		For analyses using ANSYS:
			* APDLdebug: bool
			  If it's true it will be print the dict with results imported from ANSYS
			  at each call. Great use for APDL debug.

		** If an invalid option or value is set the process could stop (or not).**

		"""

		if value == None:
			# Return current value
			return self._options[option]

		else:
			self._options[option] = value
			self._PrintR('Option \"%s\" set as \"%s\".' % (option, str(value)))



	def SetCorrel(self, var1, var2, correl):
		"""
		Set the correlation betwen two variables. The values will be transformed
		by the Nataf process before running.



		Parameters
		----------
		var1 : str, obligatory
			First variable name.

		var2 : str, obligatory
			Second variable name.

		correl : float, obligatory
			Correlation betwen var1 and var2.
		"""

		# Set lower
		var1 = var1.lower()
		var2 = var2.lower()

		# Verify if it's already created
		if var1 not in self.variableDistrib:
			exception = Exception('Variable "%s" is not declared yet. ' % var1 +
					'Before set the correlation you must create the random '+
					'variable with CreateVar().')
			raise exception

		if var2 not in self.variableDistrib:
			exception = Exception('Variable "%s" is not declared yet. ' % var2 +
					'Before set the correlation you must create the random '+
					'variable with CreateVar().')
			raise exception

		if var1 == var2:
			exception = Exception('You cannot change the correlation from a variable with itself.')
			raise exception

		if correl < -1 or correl > 1:
			exception = Exception('Correlation must be a value betwen -1 and +1.')
			raise exception

		# Store correlations on correlation list self.corlist
		self.corlist[var1, var2] = correl
		self.corlist[var2, var1] = correl

		self._PrintR('Correlation betwen \"%s\" and \"%s\" set as %f.' %(var1, var2, correl))


	def _EquivNormal(self, name, pt):
		"""
		Internal function that determines the equivalent normal distribution parameters

		Parameters
		----------
		name : str, obligatory
			Name of variable that is being transformed.

		pt : float, obligatory
			Point that is being transformed.

		Returns
		-------
		[mean, std] where:
			* mean : float
			  Equivalent normal mean

			* std : float
		      Equivalent normal standard deviation.
		"""

		# Copy the data to var
		name = name.lower()
		var = self.variableDistrib[name]

		# Gaussian
		if var[0] == 'gauss':
			# Doesn't need to transform
			mean = var[1]
			std = var[2]

		# Lognormal
		elif var[0] == 'logn':
			# LogNormal parameters
			qsi = math.sqrt(math.log(1 + (var[3])**2))
			lmbd = math.log(var[1]) - 0.5*qsi**2

			# Equiv parametes: analitycal since lognormal is derivated from normal
			std = pt*qsi
			mean = pt*(1-math.log(pt)+lmbd)

		# Gumbel
		elif var[0] == 'gumbel':
			# Gumbel parameters
			scl = math.sqrt(6)*var[2]/math.pi
			loc = var[1] - 0.57721*scl

			# Equiv with PDF/CDF functions
			# Gumbel pdf: scipy.stats.gumbel_r.pdf(x, loc, scl)
			# Gumbel cdf: scipy.stats.gumbel_r.cdf(x, loc, scl)
			# Normal pdf: scipy.stats.norm.pdf(x, mu, std)
			# Normal cdf: scipy.stats.norm.cdf(x, mu, std)
			# Normal icdf: scipy.stats.norm.ppf(q, mu, std)

			# invCum = icdf(gumbel_cdf())
			invCum = scipy.stats.gumbel_r.cdf(pt, loc, scl)
			invCum = scipy.stats.norm.ppf(invCum, 0, 1)

			std = scipy.stats.norm.pdf(invCum, 0, 1) / scipy.stats.gumbel_r.pdf(pt, loc, scl)
			mean = pt - invCum*std

		# What??
		else:
			exception = Exception('This couldnt happen!')
			raise exception

		# Return
		return [mean, std]



	def Run(self, maxIter=50, tolRel=0.01, tolLS='auto', dh=0.05, diff='forward', meth='iHLRF'):
		"""
		Run the FORM process.

		Parameters
		----------
		maxIter : integer, optional
			Maximum of iterations that can be performed. After this the process
			will stop with error.
			Defaults to 50.

		tolRel : float, optional
			Maximum **relative** error tolerance, for example on search for X point
			``|X_k - X_(k-1)|/|X_(k-1)|<=tolRel``. Defaults to 0.005.

		tolLS : float, optional
			Maximum **absolute** error tolerance for limit state function,
			``|G(X)|~=tolLS``. It should be calibrated based on the magnitude of
			limit state function. 
			
			It's possible to automatically determine it using tolLS='auto', it will be set
			as (tolRel)*(first cycle limit state value).

			Defaults to 'auto'.

		dh : float, optional
			delta_h step when applying derivatives, value applied over means, as
			h=mean(x)*dh, so, ``f'(x) = (f(x+mean(x)*dh)-f(x)) / (mean(x)*dh)``.
			Defaults to 0.05.

		diff : str, optional
			Numeric derivative calcultation method. The possible mehtods are:

				* center: for finite difference method with central difference,
				``f'(x) = (f(x+h)-f(x-h)) / (2h)``, it needs ``1 + 2*Nvars``
				evaluations of the limit state function.

				* forward: for finite difference method with forward difference,
				``f'(x) = (f(x+h)-f(x)) / h``, it needs ``1 + Nvars``
				evaluations of the limit state function.

				* backward: for finite difference method with backward difference,
				``f'(x) = (f(x)-f(x-h)) / h``, it needs ``1 + Nvars``
				evaluations of the limit state function.

				Defaults to forward.

		meth : str, optional
			FORM method used. Available methods are:

				* HLRF: Hasofer Lind Rackwitz and Fiessler method.

				* iHLRF: improved Hasofer Lind Rackwitz and Fiessler method.

			Defaults to iHLRF.

		**Returns a dictionary with:**

			* status : integer
			  Status of solution, values can be found after this list.

			* Pf : float
			  Probability of failure.

			* Beta : float
			  Reliability index.

			* {DesignPoint} : dictionary of values
			  Dictionary with the design points for each variable.

			* {gradG} : dictionary of values
			  Dictionary with the final gradient for each variable.

			* {alpha} : dictionary of values
			  Dictionary with the final director cossines for each variable.

			* cycles : int
			  Number of iterations performed to obtain the solution.


		**Status values:**

			* 0: no problem;
			* 1: warning, maximum of cycles reached with no convergence of CVPf;
			* 99: undefined error!

		"""

		#-----------------------------------------------------------------------
		# Set controls
		self._SetControls(maxIter=maxIter, tolRel=tolRel, tolLS=tolLS, dh=dh, diff=diff, meth=meth)
		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Verify the limit state equation
		if self.limstate == None:
			exception = Exception('Before Run you must define the limit state'+
								  'equation with SetLimState().')
			raise exception
		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Verify if ANSYS is being used
		#
		if self._ANSYS:
			# Verify if the model is set
			if self.ansys.Model == {}:
				exception = Exception('Before running ANSYS you must define the model that will be analysed with SetANSYSModel().')
				raise exception
		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Copy the start point to self.variableDesPt()
		#
		self.variableDesPt = self.variableStartPt.copy()
		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Create the correlation Matrix, VarId list and Jyz/Jzy matrices
		#
		NInRandVars = len(self.variableDistrib)
		NInConstVars = len(self.variableConst)
		# NOutVars is ANSYS out vars
		NOutVars = 0

		# correlMat start as an eye matrix, just with RandomVars!
		self.correlMat = np.eye(NInRandVars)

		# Create var id list
		varId = {}
		idx = 0

		# For random variables
		for eachVar in self.variableDistrib:
			varId[eachVar] = idx
			idx += 1

		# For constant variables
		for eachVar in self.variableConst:
			varId[eachVar] = idx
			idx += 1

		# and Now ANSYS output variables, if using ANSYS
		if self._ANSYS:
			NOutVars = len(self.ansys.varOutNames)
			for eachVar in self.ansys.varOutNames:
				eachVar = eachVar.lower()
				varId[eachVar] = idx
				idx += 1

		# Save varId in the object
		self.varId = varId

		# Run all the declareted correls
		for each in self.corlist:
			i = varId[each[0]]
			j = varId[each[1]]
			cor = self.corlist[each]

			# Apply Nataf to transform the correlation
			var1props = self.variableDistrib[each[0]]
			var2props = self.variableDistrib[each[1]]

			# Both are gauss
			if (var1props[0] == 'gauss' and var2props[0] == 'gauss'):
				cor = cor

			# Both are LN
			elif var1props[0] == 'logn' and var2props[0] == 'logn':
				cv1 = var1props[3]
				cv2 = var2props[3]
				cor = cor*(math.log(1+cor*cv1*cv2)/(cor*math.sqrt(math.log(1+cv1**2)*math.log(1+cv2**2))))

			# Both are Gumbel
			elif var1props[0] == 'gumbel' and var2props[0] == 'gumbel':
				cor = cor*(1.064 - 0.069*cor + 0.005*cor**2)

			# One is gauss and other is logn
			elif (var1props[0] == 'gauss' and var2props[0] == 'logn') \
				or (var2props[0] == 'gauss' and var1props[0] == 'logn'):

				# who is logn?
				if var1props[0] == 'logn':
					cv = var1props[3]
				else:
					cv = var2props[3]

				# cor is
				cor = cor*cv/math.sqrt(math.log(1+cv**2))

			# One is gauss and other is gumbel
			elif (var1props[0] == 'gauss' and var2props[0] == 'gumbel') \
				or (var2props[0] == 'gauss' and var1props[0] == 'gumbel'):
				cor = 1.031*cor

			# One is logn and other is gumbel
			elif (var1props[0] == 'logn' and var2props[0] == 'gumbel') \
				or (var2props[0] == 'logn' and var1props[0] == 'gumbel'):
				# who is logn?
				if var1props[0] == 'logn':
					cv = var1props[3]
				else:
					cv = var2props[3]

				# cor is
				cor = cor*(1.029 + 0.001*cor + 0.014*cv + 0.004*cor**2 + 0.233*cv**2 - 0.197*cor*cv)

			# Forbiden zone
			else:
				exception = Exception('When applying NATAF on variables \"%s\" and \"%s\" the variables '+
						'conditions wasn\'t possible.')
				raise exception

			# Save it!
			self.correlMat[i, j] = cor
			self.correlMat[j, i] = cor

		# Obtain the transformation matrices Jyz/Jzy
		# Jzy is Lower matrix obtained by Cholesky at correlMat
		Jzy = np.linalg.cholesky(self.correlMat)
		# and Jyz is the inverse matrix
		Jyz = np.linalg.inv(Jzy)

		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		#  CYCLEs
		#
		self._PrintR('\nStarting FORM process.')
		timei = time.time()

		# Start values
		lastcycle = False
		curVecRedPts = np.zeros(NInRandVars)
		valG = 0
		curBeta = 0
		self.results['Beta'] = []
		self.results['Beta'].append(0)

		# Cycle (or iteration) start at 1, so +1
		for cycle in range(1, 1+self.controls['maxIter']):
			#self._PrintR('----------------------------------------------------------------------------\n')
			self._PrintR('____________________________________________________________________________')
			self._PrintR('Iteration cycle %d.' % cycle)

			#-------------------------------------------------------------------
			# Mount mean vector, stddev matrix and transformations (Jxy/Jyx) matrices
			#
			matStd = np.zeros([NInRandVars, NInRandVars])
			vecMean = np.zeros(NInRandVars)
			# Vector with current design point (and constant values)
			vecPts = np.zeros(NInRandVars+NInConstVars)
			for eachVar in self.variableDistrib:
				# Variable id
				id = varId[eachVar]
				# Current X point
				pt = self.variableDesPt[eachVar]
				# Transform
				[vecMean[id], matStd[id, id]] = self._EquivNormal(eachVar, pt)
				# Put point in vecPts
				vecPts[id] = pt

			for eachVar in self.variableConst:
				# Get variable id
				id = varId[eachVar]
				# Put point in vecPts
				vecPts[id] = self.variableConst[eachVar]

			# In case of correlations we need to apply that over Jxy/Jyx
			Jxy = matStd.dot(Jzy)
			Jyx = Jyz.dot(np.linalg.inv(matStd))
			#-------------------------------------------------------------------


			#-------------------------------------------------------------------
			# Evaluate G and grad(G)
			#

			# Get dh
			dh = self.controls['dh']
			if diff == 'center':
				# size is 1+2*NInRandVars
				matEvalPts = np.zeros([(1+2*NInRandVars+NInConstVars), (NInRandVars+NInConstVars+NOutVars)])
				# All lines starts with design point/constant values, after random variables replace it
				matEvalPts[:, 0:(NInRandVars+NInConstVars)] = vecPts

				# Run all lines from [1,end] with step=2
					# curId is current varId
				curId = 0
				for eachLine in range(1, 1+2*NInRandVars, 2):
					# vecdh is not zero on current var item
					vecdh = np.zeros(NInRandVars)
					vecdh[curId] = dh
					# Odd line is +h
					matEvalPts[eachLine, 0:NInRandVars] = vecPts[0:NInRandVars] + vecMean*vecdh
					# and now -h
					matEvalPts[eachLine+1, 0:NInRandVars] = vecPts[0:NInRandVars] - vecMean*vecdh
					curId += 1

			elif diff == 'forward':
				# size is 1+NInRandVars
				matEvalPts = np.zeros([(1+NInRandVars+NInConstVars), (NInRandVars+NInConstVars+NOutVars)])
				# All lines starts with design point/constant values, after random variables replace it
				matEvalPts[:, 0:(NInRandVars+NInConstVars)] = vecPts
				curId = 0
				for eachLine in range(1, 1+NInRandVars):
					# vecdh is not zero on current var item
					vecdh = np.zeros(NInRandVars)
					vecdh[curId] = dh
					matEvalPts[eachLine, 0:NInRandVars] = vecPts[0:NInRandVars] + vecMean*vecdh
					curId += 1

			elif diff == 'backward':
				# size is 1+NInRandVars
				matEvalPts = np.zeros([(1+NInRandVars+NInConstVars), (NInRandVars+NInConstVars+NOutVars)])
				# All lines starts with design point/constant values, after random variables replace it
				matEvalPts[:, 0:(NInRandVars+NInConstVars)] = vecPts
				curId = 0
				for eachLine in range(1, 1+NInRandVars):
					# vecdh is not zero on current var item
					vecdh = np.zeros(NInRandVars)
					vecdh[curId] = dh
					matEvalPts[eachLine, 0:NInRandVars] = vecPts[0:NInRandVars] - vecMean*vecdh
					curId += 1


			#-------------------------------------------------------------------
			# If using ANSYS send it's variables dependents to simulation.
			#
			if self._ANSYS:
				# Mount a list of matEvalPts lines that should be simunlated
				# Line 0 with X value is used for all not calculated lines + G(x)
				ansysSendingList = [0]

				if diff == 'center':
					for eachVar in self.ansys.varInNames:
						# Current variable == random or constant?
						if eachVar.lower() in self.variableDistrib:
							eachVar = eachVar.lower()
							ansysSendingList.append(2*varId[eachVar]+1)
							ansysSendingList.append(2*varId[eachVar]+2)

				elif diff == 'forward' or diff == 'backward':
					for eachVar in self.ansys.varInNames:
						# Current variable is random or constant?
						if eachVar.lower() in self.variableDistrib:
							eachVar = eachVar.lower()
							ansysSendingList.append(varId[eachVar]+1)

				# Set length as Ns
				Ns = len(ansysSendingList)
				self.ansys.ClearValues()
				self.ansys.SetLength(Ns)

				# Send from the list to ANSYS varInValues
				for eachVar in self.ansys.varInNames:
					eachVar = eachVar.lower()
					self.ansys.SetVarInValues(eachVar, matEvalPts[ansysSendingList, varId[eachVar]])

				# Run ANSYS
				self.ansys.Run()

				# Get results from ANSYS
				resANSYS = self.ansys.GetVarOutValues()

				# If control APDL debug is true!
				if self._options['APDLdebug'] == True:
					self._PrintR('---\nPrinting \'resANSYS\' dict:')
					self._PrintR(resANSYS)
					self._PrintR('---\n')

				# Add ANSYS results to matEvalPts
				for eachVar in self.ansys.varOutNames:
					# First all lines has value from X
					matEvalPts[:, varId[eachVar.lower()]] = resANSYS[eachVar][0]
					# Now values are stored as in ansysSendingList
					matEvalPts[ansysSendingList, varId[eachVar.lower()]] = resANSYS[eachVar]

			#-------------------------------------------------------------------


			#-------------------------------------------------------------------
			# Evaluate Limit State (valG) and Gradient (in x) (gradGx)
			#
			self._PrintR('Evaluating limit state.')
			# dic of values in each evaluation
			varVal = {}

			# Eval Limit State:
			for eachVar in varId:
				varVal[eachVar] = matEvalPts[0, varId[eachVar]]

			# Test if limstate is a string or a function
			if(type(self.limstate) == str):
				varVal['userf'] = self._userf
				valG = eval(self.limstate, globals(), varVal)
			else:
				valG = self.limstate(**varVal)

			# tolLS 'auto' is tolRel*(initial valG)
			if cycle == 1 and tolLS == 'auto':
				tolLS = self.controls['tolRel']*valG
				self.controls['tolLS'] = tolLS
				self._PrintR('Limit state tolerance set to %f.' % tolLS)

			# Print limit state value:
			self._PrintR('Limit state value = %f (tolerance = %f).' % (valG, self.controls['tolLS']))

			# Eval Gradient
			self._PrintR('Evaluating gradient.')

			curId = 0
			gradGx = 0
			gradGx = np.zeros(NInRandVars)

			#Derivatives only for random variables/input var

			if diff == 'center':
				for eachLine in range(1, (1+2*NInRandVars), 2):
					# G(X+dh) = val1
					for eachVar in varId:
						varVal[eachVar] = matEvalPts[eachLine, varId[eachVar]]

					# Test if limstate is a string or a function
					if(type(self.limstate) == str):
						val1 = eval(self.limstate, globals(), varVal)
					else:
						val1 = self.limstate(**varVal)
						
					# G(X-dh) = val2
					for eachVar in varId:
						varVal[eachVar] = matEvalPts[eachLine+1, varId[eachVar]]
					
					# Test if limstate is a string or a function
					if(type(self.limstate) == str):
						val2 = eval(self.limstate, globals(), varVal)
					else:
						val2 = self.limstate(**varVal)

					gradGx[curId] = (val1-val2)/(2*dh*vecMean[curId])
					curId += 1

			elif diff == 'forward':
				for eachLine in range(1, (1+NInRandVars)):
					# G(X+dh) = val1
					for eachVar in varId:
						varVal[eachVar] = matEvalPts[eachLine, varId[eachVar]]
					
					# Test if limstate is a string or a function
					if(type(self.limstate) == str):
						val1 = eval(self.limstate, globals(), varVal)
					else:
						val1 = self.limstate(**varVal)

					gradGx[curId] = (val1-valG)/(dh*vecMean[curId])
					curId += 1

			elif diff == 'backward':
				for eachLine in range(1, (1+NInRandVars)):
					# G(X-dh) = val1
					for eachVar in varId:
						varVal[eachVar] = matEvalPts[eachLine, varId[eachVar]]
					
					# Test if limstate is a string or a function
					if(type(self.limstate) == str):
						val1 = eval(self.limstate, globals(), varVal)
					else:
						val1 = self.limstate(**varVal)

					gradGx[curId] = (valG-val1)/(dh*vecMean[curId])
					curId += 1

			#---------------------------------------------------------------

			#---------------------------------------------------------------
			# Get gradG (of y) by gradGx and Jxy
			#
			gradG = (Jxy.T).dot(gradGx)
			#---------------------------------------------------------------


			#---------------------------------------------------------------
			# Print gradient and alpha
			#
			self._PrintR('Gradient and direction cosines (alpha):')
			self._PrintR('         VarName | grad(g(X_i)) | alpha(i)')
			idx = 0
			absgrad = math.sqrt(gradG.dot(gradG))
			for eachVar in self.variableDesPt:
				self._PrintR(' %15s | %+12.5E | %+9.5E' % (eachVar, gradG[idx], gradG[idx]/absgrad))
				idx += 1
			self._PrintR(' ')
			#---------------------------------------------------------------


			#---------------------------------------------------------------
			# the min valu of grad must be greater than 0, if isn't print a warning 
			if min(abs(gradG)) == 0:
				self._PrintR('\n WARNING: One or more terms from gradG are equal to zero. \n Please pay attention to that!\n')
			#-------------------------------------------------------------------


			#-------------------------------------------------------------------
			# Get reduced uncorrelated points vector
			#

			# reduced uncorrelated points
			curVecRedPts = Jyx.dot(vecPts[:NInRandVars]-vecMean)

			# Current beta
			if cycle == 1:
				self.results['Beta'].append(math.sqrt(curVecRedPts.dot(curVecRedPts)))

			curBeta = self.results['Beta'][cycle]
			#-------------------------------------------------------------------


			#-------------------------------------------------------------------
			# FORM Method
			#
			if meth in ['HLRF', 'iHLRF', 'rHLRF']:
				#-------------------------------------------------------------------
				# HLRF - Rackwitz and Fiessler recursive method - Not Improved
				#
				# New reduced design point:
				newVecRedPts = gradG*(gradG.dot(curVecRedPts)-valG)*1/gradG.dot(gradG)


				#-------------------------------------------------------------------
				# Verify the convergence
				#
				if(abs(gradG.dot(curVecRedPts)) > 0.0):
					schwarzYgradG = abs(gradG.dot(curVecRedPts))/math.sqrt(gradG.dot(gradG)*curVecRedPts.dot(curVecRedPts))
				else:
					schwarzYgradG = 0.0

				self._PrintR('Schwarz inequality betwen y* and gradG = %f (it must be next to 1).' % schwarzYgradG)
				if lastcycle == True:
					if abs(valG) < self.controls['tolLS'] and (1-schwarzYgradG) < self.controls['tolRel']:
						#self._PrintR('\nFinal design point found on cycle %d.' % cycle)
						#self._PrintR('Performing a last cycle with final values.')
						#lastcycle = True
						self._PrintR('Convergence criterias checked.')
						self._PrintR('That\'s all folks.')
						self._stnumb = 0
						break
					else:
						self._PrintR('Convergence criterias not checked.')
						self._PrintR('Sorry, it wasn\'t the last cycle...')
						lastcycle = False

				#-------------------------------------------------------------------
				# If not converged yet and meth is rHLRF or iHLRF
				#
				if meth == 'rHLRF':
					dk = newVecRedPts - curVecRedPts
					newVecRedPts = curVecRedPts + self._options['rHLRF_relax']*dk

				if meth == 'iHLRF':
					#-------------------------------------------------------------------
					# iHLRF - improved Rackwitz and Fiessler recursive method
					#

					# parameters
					par_a = self._options['iHLRF_par_a']
					par_b = self._options['iHLRF_par_b']

					# direction
					dk = newVecRedPts - curVecRedPts

					# Line search
					self._PrintR('Starting line search for iHLRF step.')
					# find ck
					#	ck by Beck 2019
					val1 = math.sqrt(curVecRedPts.dot(curVecRedPts)/gradG.dot(gradG))
					val2 = 0.00

					if abs(valG) >= self.controls['tolLS']:
						# yk+dk = newVecRedPts
						val2 = 1/2*newVecRedPts.dot(newVecRedPts)/abs(valG)

					ck = max(val1, val2)*self._options['iHLRF_prod_ck'] + self._options['iHLRF_add_ck']

					# As ck must be greater than ||y*||/||gradG(y*)||
					while ck < val1:
						ck = 2*ck

					# Ck must always grow
					ck = max(ck, self._last_ck)
					self._last_ck = ck
					self._PrintR('ck value: %f.' % ck)


					#-----------------------------------------------------------
					# Find lambdk:
					#

					# Initial calcs
					#	m(y) = 1/2*y.y^T + c*|g(y)|

					# 	gradMy = grad(m(y))
					gradMy = np.zeros(NInRandVars)
					gradMy = curVecRedPts + ck*valG/abs(valG)*gradG

					#--
					# Now we search for lambdk
					#
					# 	myk = m(y)
					myk = 1/2*curVecRedPts.dot(curVecRedPts) + ck*abs(valG)

					# Find maxnk
					maxdk = max(abs(dk))
						# If b**n*max(dk) is less than tolRel, y ~= y+b**n*dk
					maxnk = math.ceil(math.log(self.controls['tolRel']/maxdk)/math.log(par_b))
						# I'm not a good guy and so we will do less!
					maxnk += 0
						# But if limit state value doesn't change anymore it will stop!
					stepnk = self._options['iHLRF_step_lambdk_test']

					# We need max(b**n), since 0<b<1 and n start as 0, where
					#	the first value that the condition is true is the value,
					#	so we can use stepnk to avoid evaluate like 100 ANSYS
					#	simulations and use the first...
					nk = 0
					valG_nk = 99999999.0
					done = False
					forcenk = False

					self._PrintR('iHLRF step range from %f to %f (%.0f steps), being test step size %d.' % (par_b**nk, par_b**maxnk, maxnk, stepnk))

					for step in range(math.ceil(maxnk/stepnk)):
						# current lenght
						curlen = min(stepnk, maxnk-(step)*stepnk)

						# points
						#matEvalPts = np.zeros([(1+2*NInRandVars+NInConstVars), (NInRandVars+NInConstVars+NOutVars)])
						matEvalPts = np.zeros([curlen, (NInRandVars+NInConstVars+NOutVars)])
						matEvalPts[:, 0:(NInRandVars+NInConstVars)] = vecPts
						lambdks = np.zeros(curlen)

						for eachnk in range(curlen):
							lambdks[eachnk] = par_b**nk
							matEvalPts[eachnk, 0:NInRandVars] = vecMean + Jxy.dot(curVecRedPts + lambdks[eachnk]*dk)
							nk += 1

						# pass ANSYSvars to ANSYS
						if self._ANSYS:
							self.ansys.ClearValues()
							self.ansys.SetLength(curlen)

							# Send from the list to ANSYS varInValues
							for eachVar in self.ansys.varInNames:
								eachVar = eachVar.lower()
								self.ansys.SetVarInValues(eachVar, matEvalPts[:, varId[eachVar]])

							# Run ANSYS
							self.ansys.Run()

							# Get results from ANSYS
							resANSYS = self.ansys.GetVarOutValues()

							# If control APDL debug is true!
							if self._options['APDLdebug'] == True:
								self._PrintR('---\nPrinting \'resANSYS\' dict:')
								self._PrintR(resANSYS)
								self._PrintR('---\n')

							# Add ANSYS results to matEvalPts
							for eachVar in self.ansys.varOutNames:
								matEvalPts[:, varId[eachVar.lower()]] = resANSYS[eachVar]


						#print(matEvalPts)
						#print(matEvalPts.shape)

						# evaluate the limit state function for eachnk
						for eachnk in range(curlen):
							# dic of values in each evaluation
							varVal = {}

							# Save last valG_nk to compare it
							valG_nk_old = valG_nk

							# Eval Limit State
							for eachVar in varId:
								varVal[eachVar] = matEvalPts[eachnk, varId[eachVar]]
							#valG_nk = eval(self.limstate, globals(), varVal)
							# Test if limstate is a string or a function
							if(type(self.limstate) == str):
								varVal['userf'] = self._userf
								valG_nk = eval(self.limstate, globals(), varVal)
							else:
								valG_nk = self.limstate(**varVal)

							# Verify the condition
							lambdk = lambdks[eachnk]
							curYk = curVecRedPts + lambdk*dk
							mynk = 1/2*curYk.dot(curYk)+ck*abs(valG_nk)
							# Correct way - REAL ARMIJO RULE
							target = -par_a*lambdk*gradMy.dot(dk)
							# Some people talk about this: gradMy.dot(gradMy), but ??
							#target = -par_a*lambdk*gradMy.dot(gradMy)

							if (mynk - myk) <= target:
								self._PrintR('iHLRF step size is %f.' % lambdk)
								done = True
								break

							if valG_nk_old/valG_nk == 1:
							#if valG_nk_old/valG_nk >= (1-self.controls['tolRel']):
								forcenk = True
								break

						if done == True or forcenk == True:
							break

					if done == False or forcenk == True:

						if self._options['iHLRF_forced_lambdk'] == 'auto':
							lambdk = 1 - schwarzYgradG
						else:
							lambdk = self._options['iHLRF_forced_lambdk']

						self._PrintR('iHLRF step not found, forcing to %f.' % lambdk)

					# Save new reduced point
					newVecRedPts = curVecRedPts + lambdk*dk
			else:
				exception = Exception('FORM method \"%s\" is not implemented.' % meth)
				raise exception


			#-------------------------------------------------------------------
			# Transform reduced points back to real space
			#
			newBeta = math.sqrt(newVecRedPts.dot(newVecRedPts))
			self.results['Beta'].append(newBeta)

			NewVecPts = vecMean + Jxy.dot(newVecRedPts)

			# Save it to self.variableDesPt
			for eachVar in self.variableDesPt:
				self.variableDesPt[eachVar] = NewVecPts[varId[eachVar]]

			#-------------------------------------------------------------------


			#-------------------------------------------------------------------
			# Tell user how is it going
			#
			self._PrintR(' ')
			self._PrintR('  Current Beta = %2.3f.' % curBeta)
			self._PrintR(' ')
			self._PrintR('         VarName | Next Design Point')
			idx = 0
			for eachVar in self.variableDesPt:
				self._PrintR(' %15s | %8.5E' % (eachVar, self.variableDesPt[eachVar]))
				idx += 1
			self._PrintR(' ')
			#-------------------------------------------------------------------


			#-------------------------------------------------------------------
			# Verify the convergence
			#
			# Here we have a problem: I don't know from where Beck get this schwarzYgradG verificarion!
			# I've used that, but now I changed to the old fashion way.
			# To use that again you need to change the next verification and remove the
			# conditional lastcycle below too. 
			#
			#### Sometimes it didn't converge because if this crap
			if np.linalg.norm(newVecRedPts) > 0:
				relErrorPoint = max(abs((curVecRedPts-newVecRedPts)/newVecRedPts))
			else:
				relErrorPoint = 1.0
			absErrorBeta = abs(curBeta-newBeta)
			self._PrintR('Maximum relative error on design point = %1.4f.' % relErrorPoint)
			self._PrintR('Absolute error betwen current and next beta = %1.4f.' % absErrorBeta)

			#### if abs(valG) < self.controls['tolLS'] and (1-schwarzYgradG) < self.controls['tolRel']:
			#### if abs(valG) < self.controls['tolLS'] and relErrorPoint < self.controls['tolRel']:

			if abs(valG) < self.controls['tolLS']:
				# limstate value is ok
				if absErrorBeta < self.controls['tolRel']:
					# Converged by absolute difference betwen betas
					self._PrintR('\nFinal design point found on cycle %d by absolute difference betwen two betas.' % cycle)
					self._stnumb = 0
					lastcycle = True
					break
				elif (1-schwarzYgradG) < self.controls['tolRel']:
					lastcycle = True
					self._PrintR('\nFinal design point was probably found on cycle %d by Schwarz inequality betwen y* and gradG.' % cycle)
					self._PrintR('A new cycle will be started to confirm the limit state value and it\'s gradient.')
				else:
					lastcycle = False
			else:
				lastcycle = False

			#-------------------------------------------------------------------


		#-----------------------------------------------------------------------

		#-----------------------------------------------------------------------
		# After CYCLEs
		#
		timef = time.time()
		self.results['ElapsedTime'] = ((timef-timei)/60)

		# Save last cycles
		self._cycles = cycle

		# Get results
		#self.results['Beta'] = math.sqrt(newVecRedPts.dot(newVecRedPts))
		self.results['Pf'] = scipy.stats.norm.cdf(-self.results['Beta'][-1])

		# Put grad and alpha in self.results
		self.results['grad'] = {}
		self.results['alpha'] = {}
		absgrad = math.sqrt(gradG.dot(gradG))
		idx = 0
		for eachVar in self.variableDistrib:
			self.results['grad'][eachVar] = gradG[idx]
			self.results['alpha'][eachVar] = gradG[idx]/absgrad
			idx += 1

		self._PrintR('\n\n============================================================================\n')
		# Verify if stopped without convergence
		if cycle == self.controls['maxIter']:
			self._stnumb = 1
			self._PrintR(' WARNING:')
			self._PrintR('  The process was finished after reach the limit of iterations without')
			self._PrintR('  reach the error tolerance.\n')

		self._PrintR(' Total of iterations: %3.3E' % (self._cycles))
		self._PrintR(' Probability of failure (Pf): %2.4E' % (self.results['Pf']))
		self._PrintR(' Reliability index (Beta): %2.3f' % (self.results['Beta'][-1]))
		self._PrintR(' Elapsed time: %4.2f minutes.' % (self.results['ElapsedTime']))

		self._PrintR(' Final values:')
		self._PrintR('         VarName | D. Point    | grad(g(X_i)) | alpha(i) ')
		for eachVar in self.variableDesPt:
			self._PrintR(' %15s | %8.5E | %+12.5E | %+9.5E' % (eachVar, \
				self.variableDesPt[eachVar], self.results['grad'][eachVar], self.results['alpha'][eachVar]))
			#self._PrintR(' %s = %3.5E' % (eachVar, self.variableDesPt[eachVar]))
		self._PrintR('\n============================================================================\n\n')

		#-------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Send the return
		#

		# put all in a dic
		finret = {}
		finret['status'] = self._stnumb
		finret['Pf'] = self.results['Pf']
		finret['Beta'] = self.results['Beta'][-1]
		finret['DesignPoint'] = self.variableDesPt
		finret['gradG'] = self.results['grad']
		finret['alpha'] = self.results['alpha']
		finret['cycles'] = self._cycles

		return finret

		#-----------------------------------------------------------------------


	def ExportDataCSV(self, filename, description=None):
		"""
		Exports process data to a CSV file.

		Parameters
		----------
		filename : str, obligatory
			Name of file that will receive the values, doesn't need the
			extension ".csv", it will be placed automatically.

		description : str, optional
			A string that will be write in the beggining of the file.
		"""

		# Open file
		filename = filename+'.csv'
		try:
			f = open(filename, 'wt')
		except:
			exception = Exception('Unable to open the \"%s\" file for write data.' % filename)
			raise exception
		else:
			# Starts with sep=, for Microsoft Excel
			f.write('sep=,\n')

			# Description
			if description != None:
				f.write('%s\n\n' % description)

			f.write('Input data:\n')

			# Simulation Controllers:
			f.write('Process Controllers:\n')
			f.write(',Limit of iterations:,%d\n' % self.controls['maxIter'])
			f.write(',Relative error tolerance:,%2.3E\n' % self.controls['tolRel'])
			f.write(',Absolute LS error tolerance:,%2.3E\n' % self.controls['tolLS'])
			f.write(',deltah (for derivatives):,%2.3E\n' % self.controls['dh'])
			f.write(',Finite difference method:,%s\n' % self.controls['diff'])
			f.write(',FORM Method:,%s\n' % self.controls['meth'])

			# ANSYS Properties
			if self._ANSYS:
				f.write('\n')
				f.write('ANSYS Properties:\n')
				f.write(',ANSYS Model:\n')
				f.write(',,Model:,%s\n' % self.ansys.Model['inputname'])
				f.write(',,Extra files:,%s\n' % self.ansys.Model['extrafiles'])
				f.write(',,Input dir.:,%s\n' % self.ansys.Model['directory'])

				f.write(',ANSYS Input variables:\n')
				for eachVar in self.ansys.varInNames:
					f.write(',,%s\n' % eachVar)

				f.write(',ANSYS Output variables:\n')
				for eachVar in self.ansys.varOutNames:
					f.write(',,%s\n' % eachVar)

			f.write('\n')

			# Random variables
			f.write('Random variables:\n')
			f.write(',Name,Distribution,Mean,Standard Deviation,CV,Par1,Par2\n')
			for eachVar in self.variableDistrib:
				values = self.variableDistrib[eachVar]
				cmd = ',%s,%s' % (eachVar, values[0])
				for eachVal in values[1:]:
					cmd = '%s,%8.5E' % (cmd, eachVal)
				f.write('%s\n' % cmd)
			f.write('\n')

			# Constant variables
			f.write('Constant variables:\n')
			f.write(',Name,Value\n')
			for eachVar in self.variableConst:
				cmd = ',%s,%8.5E' % (eachVar, self.variableConst[eachVar])
				f.write('%s\n' % cmd)
			f.write('\n')

			# Correlation Matrix
			f.write('Correlation matrix:\n')
			# First line with Varnames:
			cmd = ','
			idx = 0
			for eachVar in self.varId:
				# Just random variables!
				if idx >= len(self.variableDistrib):
					break
				cmd = '%s,%s' % (cmd, eachVar)
				idx += 1
			cmd = '%s\n' % cmd
			f.write(cmd)
			# Matrix lines with first column as Varname
			idx = 0
			for eachLine in self.correlMat:
				cmd = ',%s' % list(self.varId.keys())[idx] # WTF DID I DO? But it works very well... xD
					# GO HORSE, GO GO GO GO GO!
				for eachVal in eachLine:
					cmd = '%s,%f' % (cmd, eachVal)
				cmd = '%s\n' % cmd
				f.write(cmd)
				idx += 1
			f.write('\n')


			# Limit state
			f.write('Limit State:,"%s"\n' % (self.limstate))
			f.write('\n')

			# Initial point
			f.write('Initial design point:\n')
			for eachVar in self.variableStartPt:
				f.write(',%s,%8.5E\n' % (eachVar, self.variableStartPt[eachVar]))
			f.write('\n')



			# Results part
			f.write('\nResults:\n')

			# Final results
			f.write('FORM results:\n')
			f.write(',Exit status:,%d,\n' % self._stnumb)
			f.write(',Total of iterations:,%d\n' % (self._cycles))
			f.write(',Probability of failure (Pf):,%2.4E\n' % self.results['Pf'])
			f.write(',Reliability Index (Beta):,%2.3f\n' % self.results['Beta'][-1])
			f.write(',Elapsed time (minutes):,%4.3f\n' % self.results['ElapsedTime'])
			f.write('\n')

			# Final Sampling Point
			f.write('Final values:\n')
			f.write(',Variable,D. Point,grad(g(X_i)),alpha(i)\n')
			for eachVar in self.variableDesPt:
				f.write(',%s,%8.5E,%11.5E,%9.5E\n' % (eachVar, \
					self.variableDesPt[eachVar], self.results['grad'][eachVar], self.results['alpha'][eachVar]))
			f.write('\n')

			# Beta values
			f.write('Beta indexes of cycles:\n')
			f.write(',Cycle,Beta\n')
			idx = 1
			for eachBeta in self.results['Beta'][1:]:
				f.write(',%d,%2.3f\n' % (idx, eachBeta))
				idx += 1

			# End
			f.close()

			self._PrintR('Simulation data exported to "%s".' % filename)
