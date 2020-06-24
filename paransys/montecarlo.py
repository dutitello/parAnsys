# -*- coding: UTF-8 -*-
"""
Module to use ANSYS by Probabilistic Design System (PDS) as a FEM tool for
running Monte Carlo simulations in Python.

This module works based on ParAnsys.ANSYS.

Docs are also available in https://dutitello.github.io/paransys/ 

"""

import os
import time
from paransys.ansys import ANSYS
import numpy as np
import scipy.stats
import math
from math import *

class MonteCarlo(object):
	"""
	This class performns Monte Carlo simulations inside Python using ParAnsys as
	a connection with ANSYS for evaluate FEM models.

	It is possible to run Monte Carlo simulations without using ANSYS, just
	defining the limit state equation and all variables.

	This code was made following the ideia of ANSYS being a tool for getting the
	ultimate load of the structure. This works applying a displacement in the
	loaded node, and then getting the biggest reaction force on that node,
	following this way the limit state defined here is 'R-S', where R are the
	values get from ANSYS and S the values generated in Python. It's also
	possible to work applying the true load on ANSYS, it's just necessary to
	formulate a valid limit state equation.

	|

	ATTENTION: When using ANSYS the weight of results from ANSYS variables
	are determined using the weights of all ANSYS input variables.
	
	|

	**To do**

	1) When structure has more than one limit state the PDF of sampling
	   distribution is the sum of all limit states sampling distributions vs their
	   sampling weights (h(x) = w1.h1(x) + w2.h2(x) + hi.wi(x)...)
	   It's already done the division of simulations for each cycle with the
	   limit state weights.

	2) When sampling distribution is different of real distribution Pf is going
	   wrong, so it's not able to be used, for now.
		
	|
	|

	**Class methods:**

	
	"""

	#---------------------------------------------------------------------------



	def __init__(self):
		"""

		"""
		# Definition of self varibles that just exists on MonteCarlo
		# Dictionary with variables distrib and parameters
		self.variableDistrib = {}
		self.samplingDistrib = {}
		self.variableConst = {}
		self.corlist = {}
		self.controls = {}
		self.ansysvars = []
		self.limstates = {}

		#
		self.PrintR = False

		# Control if ANSYS is being used
		self._ANSYS = False
		pass


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
			extrafiles = ['temps.txt', 'model1.ans', 'file.db']

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


	# Internal for variable distrib definition
	def _VarDistrib(self, type, name, distrib, mean, std, cv, par1, par2, limst=0):
		"""
		Internal function that verify the values for CreateVar
		and SetRandomVarSampl.

		Parameters
		----------
		type : int, obligatory
			type = 0: for CreateVar   -> variableDistrib
			type = 1: for SetRandomVarSampl -> samplingDistrib

		limst : integer, obligatory for type=1 with more than 1 limit state
			Limit state ID that will use current sampling distribution.

		"""
		# Verify the variable name
		if name in ['', None, ' ']:
			exception = Exception('You must define the name of the variable that will receive the values.')
			raise exception
		else:
			name = name.lower()
			# Sampling distribution not created yet
			if type is 1 and name not in self.variableDistrib:
				exception = Exception('Variable %s is not declared yet. '+
						'Before set sampling distribution you must create it '+
						'with CreateVar().')
				raise exception

		# Set distribution name as lower case
		distrib = distrib.lower()

		# CV or STD?
		if distrib not in ['constant', 'const', 'cons', 'c']:
			if cv is not None:
				std = cv*mean
			elif mean == 0.0:
				cv = 1e99
			else:
				cv = std/mean

		# Verify the distribution and then determine the parameters
		# Gaussian variable
		if distrib in ['gauss', 'gaus', 'gaussian', 'normal', 'norm']:
			# Type 0 = Random Variable distribution
			if type is 0:
				self.variableDistrib[name] = ['gauss', mean, std, cv]
				return self._PrintR('Variable %s defined as Gaussian with mean=%f, std. dev.=%f, CV=%f.'
									% (name, mean, std, cv))

			# Type 1 = Sampling distribution
			elif type is 1:
				# Vrify of samplingDistrib[limst] exists
				try:
					self.samplingDistrib[limst]
				except:
					# Not created yet
					self.samplingDistrib[limst] = {}
					self.samplingDistrib[limst][name] = ['gauss', mean, std, cv]
				else:
					# Already created
					self.samplingDistrib[limst][name] = ['gauss', mean, std, cv]

				return self._PrintR('Variable %s sampling defined as Gaussian with mean=%f, std. dev.=%f, CV=%f for limit state %d.'
									% (name, mean, std, cv, limst))

		# Lognormal distribution
		elif distrib in ['lognormal', 'logn', 'ln', 'log', 'lognorm']:
			# Distribution parameters are qsi, related do std dev, and
			#	lmbd, related to mean.
			#qsi = math.sqrt(math.log(1 + (cv)**2))
			#lmbd = math.log(mean) - 0.5*qsix**2

			# Type 0 = Random Variable distribution
			if type is 0:
				self.variableDistrib[name] = ['logn', mean, std, cv]
				return self._PrintR('Variable %s defined as LogNormal with mean=%f, std. dev.=%f, CV=%f.'
									% (name, mean, std, cv))

			# Type 1 = Sampling distribution
			elif type is 1:
				# Vrify of samplingDistrib[limst] exists
				try:
					self.samplingDistrib[limst]
				except:
					# Not created yet
					self.samplingDistrib[limst] = {}
					self.samplingDistrib[limst][name] = ['logn', mean, std, cv]
				else:
					# Already created
					self.samplingDistrib[limst][name] = ['logn', mean, std, cv]

				return self._PrintR('Variable %s sampling defined as LogNormal with mean=%f, std. dev.=%f, CV=%f for limit state %d.'
									% (name, mean, std, cv, limst))


		# Gumbel distribution
		elif distrib in ['gumbel', 'gumb', 'type1']:

			# Type 0 = Random Variable distribution
			if type is 0:
				self.variableDistrib[name] = ['gumbel', mean, std, cv]
				return self._PrintR('Variable %s defined as Gumbel with mean=%f, std. dev.=%f, CV=%f.'
									% (name, mean, std, cv))

			# Type 1 = Sampling distribution
			elif type is 1:
				# Vrify of samplingDistrib[limst] exists
				try:
					self.samplingDistrib[limst]
				except:
					# Not created yet
					self.samplingDistrib[limst] = {}
					self.samplingDistrib[limst][name] = ['gumbel', mean, std, cv]
				else:
					# Already created
					self.samplingDistrib[limst][name] = ['gumbel', mean, std, cv]

				return self._PrintR('Variable %s sampling defined as Gumbel with mean=%f, std. dev.=%f, CV=%f for limit state %d.'
									% (name, mean, std, cv, limst))

		# Constant value
		elif distrib in ['constant', 'const', 'cons', 'c']:
			# Type 0 = Random Variable distribution
			if type is 0:
				self.variableConst[name] = mean
				return self._PrintR('Variable %s defined as Constant with value=%f'
									% (name, mean))

			# Type 1 = Sampling distribution
			elif type is 1:
				# Not possible, break
				exception = Exception('Error on %s: You can not use a constant value for sampling!'
										% (distrib.upper()))
				raise exception

		else:
			exception = Exception('Distribution %s set on variable %s is not recognized.'
									% (distrib.upper(), name))
			raise exception



	# Setting distribution of variables
	def CreateVar(self, name, distrib, mean, std=0, cv=None, par1=None, par2=None):
		"""
		Create a Random Variable

		If it's used on ANSYS it need to be told, so after this use:

		>>> mc.SetANSYSVar(name)

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

		self._VarDistrib(type=0, name=name, limst=None, distrib=distrib, mean=mean, std=std, cv=cv, par1=par1, par2=par2)



	# Sampling definition
	def SetRandomVarSampl(self, name, limst, distrib, mean, std=0, cv=None, par1=None, par2=None):
		"""
		Sets the sampling distribution of a variable to performn Importance
		Sampling the simulations.

		ATTENTION: When using ANSYS the weight of results from ANSYS variables
		are determined using the weights of all ANSYS input variables.

		Parameters
		----------
		name : str, obligatory
			Name of variable.

		limst : integer, obligatory for type=1 with more than 1 limit state
			Limit state ID that will use current sampling distribution.

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

		if limst > (len(self.limstates)-1):
			exception = Exception('Limit state %d is not created yet. Please create it before set sampling distribution.' % limst)
			raise exception


		self._VarDistrib(type=1, name=name, limst=limst, distrib=distrib, mean=mean, std=std, cv=cv, par1=par1, par2=par2)


	
	
	def SetCorrel(self, var1, var2, correl):
		"""
		Set the correlation betwen two variables. 
		
		The values will be transformed by the Nataf process before running.


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
                         'Before set the correlation you must create the random ' +
                         'variable with CreateVar().')
			raise exception

		if var2 not in self.variableDistrib:
			exception = Exception('Variable "%s" is not declared yet. ' % var2 +
                         'Before set the correlation you must create the random ' +
                         'variable with CreateVar().')
			raise exception

		if var1 == var2:
			exception = Exception(
			    'You cannot change the correlation from a variable with itself.')
			raise exception

		if correl < -1 or correl > 1:
			exception = Exception('Correlation must be a value betwen -1 and +1.')
			raise exception

		# Store correlations on correlation list self.corlist
		self.corlist[var1, var2] = correl
		self.corlist[var2, var1] = correl

		self._PrintR('Correlation betwen \"%s\" and \"%s\" set as %f.' %
		             (var1, var2, correl))



	def _SetControls(self, Ns, Nmaxcycles, CVPf=0.00, tolAdPt=False):
		"""
		Set the controls of simulation process.

		Parameters
		----------
		Ns : integer, obligatory
			Number of simulations performed on each cycle.
			After each cycle the convergence of simualtion is verified.
			When using Importance Sampling with Adaptive Sampling, after each
			cycle the new sampling point will be determined.

		Nmaxcycles : integer, obligatory
			Maximum number of cycles to be performed, if CVPf is not reached on
			Nmaxcycles the simulation will be interrupted.

		CVPf : float, optional
			Target value of Probability Failure Coefficient of Variation, when
			reached the simulation stops.

		tolAdPt : float or False, optional
			Maximum relative tolerance for adaptive sampling point search.
			If the value is "False" it disable adaptive sampling, simulations
			will use always the user set point.

		"""
		if Ns >= 0 and Nmaxcycles >= 0 and CVPf >= 0 and (tolAdPt is False or tolAdPt >= 0):
			# Save controls variable
			self.controls['Ns'] = Ns
			self.controls['Nmaxcycles'] = Nmaxcycles
			self.controls['CVPf'] = CVPf
			self.controls['tolAdPt'] = tolAdPt

			self._PrintR('Process controls set as:')
			self._PrintR('   Simulations per cycle: %d' % Ns)
			self._PrintR('   Maximum number of cycles: %d' % Nmaxcycles)
			self._PrintR('   Total maximum number of simulations: %2.3E' % (Ns*Nmaxcycles))
			self._PrintR('   Target CVPf: %2.3E' % CVPf)
			if tolAdPt is not False:
				self._PrintR('   Maximum relative tolerance for adaptive sampling point search: %2.3E' % tolAdPt)
		else:
			exception = Exception('Error while setting simulation controls. Please verify the set values.')
			raise exception


	def CreateLimState(self, equat, weight=1.00, userf=None):
		"""
		Create and Set a new limit state.*

		The number ID of LimitStates are generated automatically starting at 0
		for the first.

		*** Current version supports only one limit state!**

		ATTENTION: When using ANSYS the weight of results from ANSYS variables
		are determined using the weights of all ANSYS input variables.

		Parameters
		----------
		equat : str, obligatory
			String with the equation of the limit state. It must be write as a
			function of defined variables (In and Out).

		weight : float, obligatory only with more than 1 limit state
			The weight of current limit state, it determines how the simulations
			are distributed betwen all the limit states.
			The sum of all limit states must be 1.00, so, if there is just one
			limit state it's weight should be 1.00

		userf : function, optional
			An user defined function that could be used inside the limit state
			equation, called inside equat as ``userf()``. Each limit state has it's
			own userf, but you can use the same Python function for all limit states.
			For example, you can create a complex Python function with loops, ifs
			and whatever for evaluate the R part of your limit state function
			for a concrete beam. An example is showed after.

		First example: if ANSYS returns the maximum load on a truss as variable
		FxMAX, and applied loads to be tested are ``(g+q)*sin(theta)``, where
		``g``, ``q``, theta are defined random variables created with ``CreateVar()``.

		.. code-block:: python

			mc.CreateLimState(equat='FxMAX-(g+q)*sin(theta)', weight=1.00)

		Note that you can use math expressions as ``sin()``, ``cos()``, ``tan()``, ``sqrt()``
		from Python math module inside the equation.

		|

		Second example: you have a steel bar in tension that hasn't hardening.
		It's stress is a function of ``(def, fy, E)``, where ``def`` is current
		deformation, ``fy`` is yield stress and ``E`` the elastic moduli,
		you can create inside your code	an function like:

		.. code-block:: python

			def stress(def, fy, E):
				if def > fy/E:
					return fy
				else:
					return def*E

		And now defining ``userf=stress`` we can:

		.. code-block:: python

			mc.CreateLimState(equat='userf(def,fy,E)-q', weight=1.00, userf=stress)

		where ``def``, ``fy``, ``E`` and ``q`` are random variables.
		Note that the function inside the limit state equation should be
		called as ``userf()`` with the parameters from ``stress``.

		Or we can do the same using the functions instead of the string:

		.. code-block:: python

			mc.CreateLimState(equat=stress)

		"""
		if weight <= 0 or weight > 1.00:
			exception = Exception('The weigth of limit state must be greater than zero and less equal to 1.00.')
			raise exception

		# Add current limit state to the list
		act = len(self.limstates)

		#---------------
		#
		#	TAKE A LOOK HERE !
		#
		#
		#
		#
		#
		#
		#
		# For now just 1 limit state with weight=1
		# 	this part prevents more than 1 limit states.
		#
		#
		#
		#
		weight = 1
		if act >= 1:
			exception = Exception('Current version of ParAnsys does not support more than one limit states.')
			raise exception
		#else:
		#	self._PrintR('\n\nCurrent version of ParAnsys supports only one limit state, do not try to set another.\n\n')
		#
		#
		#
		#
		#
		#
		#---------------
		
		if(type(equat) is str):
			# String equation
			# Change equation to lowcase
			equat = equat.lower()
		else:
			self._userf = None
	
		# Fourth value limstates[act][3] will be the Ns for LS
		self.limstates[act] = [equat, weight, userf, 0]

		#
		try:
			self.samplingDistrib[act]
		except:
			self.samplingDistrib[act] = {}

		# Verify the sum off all limit states until now
		#	itsn't possible to verify the equantion since we dont have ANSYS out values
		tempsum = 0.00
		for each in self.limstates:
			tempsum += self.limstates[each][1]

		if tempsum > 1.00:
			exception = Exception('The sum of limit states weights exceeded 1.00, please verify the inserted values.')
			raise exception

		return self._PrintR('Limit state %d defined as "%s" with weigth %f, the sum of weights untill now is %f.' %(act, equat, weight, tempsum))



	def SetANSYSVar(self, name):
		"""
		Mark a Random variable as ANSYS variable.

		ATTENTION: When using ANSYS the weight of results from ANSYS variables
		are determined using the weights of all ANSYS input variables.

		Parameters
		----------
		name : str, obligatory
			Name of variable.

		"""

		# Verify ANSYS object
		if self._ANSYS is False:
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
				exception = Exception('This variable name is not declared. '+
						'Only Random variables can be set as ANSYS variables.\n'+
						'Please use CreateVar() to declare it.')
				raise exception

		# Declare it on ansys object
		self.ansys.CreateVarIn(name)

		# If exists and isn't an ANSYS variable yet append to ansysvars
		#if name not in self.ansysvars:
		#	self.ansysvars.append(name)
		#	return self._PrintR('Variable %s set as an ANSYS variable.' % name)
		#else:
		#	return self._PrintR('Variable %s is already set as an ANSYS variable.' % name)

	def _GenRandomVW(self, varDistrib, sampDist, NCValues_h):
		"""
		Internal function for Generate Random Values and it's Weights
			(when sampDist not False)


		Parameters
		----------
		varDistrib : float, obligatory
			A list with the real distribution parameters,
			[distrib, mean, std, par1, par2]

		sampDist : float, obligatory
			A list with the sampling distribution parameters, like varDistrib.
			It also could be "False", in this case there is no sampling distribution.

		NCValues_h : 1D np.array
			INPUT Normalized correlated (or not) values from sampling distirbuition (h_X(x))


		Returns
		-------
		[randval, weights, NCValues_f] where:
			* randval : float
			  Random values generated.

			* weights : float
			  Weight of each value, when sampDist being used.

			* NCValues_f : 1D np.array
			  OUTPUT Normalized correlated (or not) values for real distribution (f_X(x))
		"""

		# Initial values
		Nsi = len(NCValues_h)
		randval = np.zeros(Nsi)
		weights = np.zeros(Nsi)
		NCValues_f = np.zeros(Nsi)

		# Verify the existence of sampling distrib
		if sampDist is False:
			# Current distribution for generate values is varDistrib
			curDist = varDistrib
		else:
			curDist = sampDist

		# Generate it
		# Gauss distribution
		if curDist[0] is 'gauss':
			# Sampling distribution values (X_h)
			randval = curDist[1] + curDist[2]*NCValues_h

			# Weights
			# Not sampling
			if sampDist is False:
				NCValues_f = NCValues_h.copy()
				weights[:] = 1

			# Real distrib is gaussian:
			elif varDistrib[0] is 'gauss':
				# Real normal reduced correlated values (Z_f)
				NCValues_f = (randval - varDistrib[1])/varDistrib[2]

				# Weights 
				weights = (scipy.stats.norm.pdf(randval, varDistrib[1], varDistrib[2]) / scipy.stats.norm.pdf(NCValues_f, 0, 1)) \
						/ (scipy.stats.norm.pdf(randval, curDist[1], curDist[2]) / scipy.stats.norm.pdf(NCValues_h, 0, 1))

			# Real distrib is logn:
			#elif varDistrib[0] is 'logn':
				# eq for varDistrib
			#	qsiv = math.sqrt(math.log(1 + (varDistrib[2]/varDistrib[1])**2))
			#	lambdav = math.log(varDistrib[1]) - 0.5*qsiv**2

			#	weights = (1/(qsiv*randval)*scipy.stats.norm.pdf(np.log(randval), lambdav, qsiv)) \
			#			 / scipy.stats.norm.pdf(randval, curDist[1], curDist[2])

			else:
				exception = Exception('For now it is not possible to use different distributions in the same var.')
				raise exception

		# Lognormal distribution
		elif curDist[0] is 'logn':
			# equiv std dev
			qsic = math.sqrt(math.log(1 + (curDist[3])**2))
			# equiv mean
			lambdac = math.log(curDist[1]) - 0.5*qsic**2

			# Random values
			randval = np.exp(lambdac + qsic*NCValues_h)
			#np.random.lognormal(lambdac, qsic, Nsi)

			# Weights
			# Not sampling
			if sampDist is False:
				NCValues_f = NCValues_h.copy()
				weights[:] = 1

			# Real distrib is logn:
			elif varDistrib[0] is 'logn':
				# eq for varDistrib
				qsiv = math.sqrt(math.log(1 + (varDistrib[3])**2))
				lambdav = math.log(varDistrib[1]) - 0.5*qsiv**2

				# Real normal reduced correlated values (Z_f)
				NCValues_f = (lambdac-lambdav+qsic*NCValues_h)/qsiv

				weights = (scipy.stats.norm.pdf(np.log(randval), lambdav, qsiv)/qsiv/scipy.stats.norm.pdf(NCValues_f)) \
						/ (scipy.stats.norm.pdf(np.log(randval), lambdac, qsic)/qsic/scipy.stats.norm.pdf(NCValues_h))

				# old and wrong! who is s???
				# scipy.stats.lognorm.pdf(x, s, loc, scale)
				#weights = scipy.stats.lognorm.pdf(randval, s, lambdav, qsiv) \
				#		/  scipy.stats.lognorm.pdf(randval, s, lambdac, qsic)

			# Real distrib is gaussian:
			#elif varDistrib[0] is 'gauss':
			#	# weights = pdfgauss/pdflogn
			#	weights = scipy.stats.norm.pdf(randval, varDistrib[1], varDistrib[2]) \
			#			/ (scipy.stats.norm.pdf(np.log(randval), lambdac, qsic)/(qsic*randval))

			else:
				exception = Exception('For now it is not possible to use different distributions in the same var.')
				raise exception

		# Gumbel distribution
		elif curDist[0] is 'gumbel':
			# scale parm - eqv to std
			sclc = math.sqrt(6)*curDist[2]/math.pi

			# loc parm - equiv to mean
			locc = curDist[1] - 0.57721*sclc

			# Random values
			randval = scipy.stats.gumbel_r.ppf(scipy.stats.norm.cdf(NCValues_h), loc=locc, scale=sclc)
			#randval = np.random.gumbel(locc, sclc, Nsi)

			# Weights
			# Not sampling
			if sampDist is False:
				NCValues_f = NCValues_h.copy()
				weights[:] = 1

			# Real distrib is gumbel
			elif varDistrib[0] is 'gumbel':
				# eq for varDistrib
				sclv = math.sqrt(6)*varDistrib[2]/math.pi
				locv = varDistrib[1] - 0.57721*sclv

				# Real normal reduced correlated values (Z_f)
				NCValues_f = scipy.stats.norm.ppf(scipy.stats.gumbel_r.cdf(randval, loc=locv, scale=sclv))

				weights = (scipy.stats.gumbel_r.pdf(randval, loc=locv, scale=sclv)/scipy.stats.norm.pdf(NCValues_f)) \
						/ (scipy.stats.gumbel_r.pdf(randval, loc=locc, scale=sclc)/scipy.stats.norm.pdf(NCValues_h))

			else:
				exception = Exception('For now it is not possible to use different distributions in the same var.')
				raise exception

		# Return
		return [randval, weights, NCValues_f]


	def Run(self, Ns, Nmaxcycles, CVPf=0.00, tolAdPt=False):
		"""
		Run the Monte Carlo simulation.

		Parameters
		----------
		Ns : integer, obligatory
			Number of simulations performed on each cycle.
			After each cycle the convergence of simualtion is verified.
			When using Importance Sampling with Adaptive Sampling, after each
			cycle the new sampling point will be determined.

		Nmaxcycles : integer, obligatory
			Maximum number of cycles to be performed, if CVPf is not reached on
			Nmaxcycles the simulation will be interrupted.

		CVPf : float, optional
			Target value of Probability Failure Coefficient of Variation, when
			reached the simulation stops.

		tolAdPt : float or False, optional
			Maximum relative tolerance for adaptive sampling point search.
			If the value is "False" it disable adaptive sampling, simulations
			will use always the user set point.


		**Returns a dictionary with:**

			* stnumb : integer
			  Status of solution, values can be found after this list.

			* Pf : float
			  Probability of failure.

			* Beta : float
			  Reliability index.

			* CVPf : float
			  Coefficient of Variation of Probability of failure

			* {SamplingPoints} : dictionary of dictionaries
			  Dictionary with sampling points used, or founded in case of
			  adaptive sampling, for each Variable on each Limit State.
			  (SamplingPoints[eachLS][eachVar])

			* cycles : int
			  Number of cycles performed to obtain the solution.

			* distparms : dictionary of dictionaries
			  Return mean (gMean) and standart deviation (gStd) of each limit state function.  


		**Status values:**

		* 0: no problem;
		* 1: warning, maximum of cycles reached with no convergence of CVPf;
		* 99: undefined error!

		"""
		#---------------
		#
		#
		#
		#	TAKE A LOOK HERE !
		#
		#
		#
		#
		# 	For now just 1 limit state with weight=1
		#
		#
		#
		#
		#
		#
		#
		#---------------


		#-----------------------------------------------------------------------
		# Set the controls

		self._SetControls(Ns=Ns, Nmaxcycles=Nmaxcycles, CVPf=CVPf, tolAdPt=tolAdPt)

		#if self.controls == {}:
		#	exception = Exception('Before Run the Monte Carlo simulation you '+
		#						  'must define the controls of simulation '+
		#						  'process with SetControls().')
		#	raise exception

		# Adaptive condition
		if tolAdPt is not False:
			adapt = True
		else:
			adapt = False

		# Time in the beggining
		timei = time.time()

		# It starts with undefined error
		stnumb = 99

		self._PrintR('\nStarting simulation cycles.')
		#-----------------------------------------------------------------------


		# Before all we have to set current
		#
		#
		# For each cycle we need to generate all data, pass ANSYS data for
		# self.ansys object, run it, get results, process monte carlo results,
		# calibrate new point (if necessary)
		#

		#=======================================================================
		# Before first Run:
		#	Create variables;
		#	Set initial values;
		#	Verify controls of process;
		#
		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Verify if ANSYS is being used
		if self._ANSYS:
			# Set Length
			self.ansys.SetLength(Ns)

			# Verify if the model is set
			if self.ansys.Model == {}:
				exception = Exception('Before running ANSYS you must define the model that will be analysed with SetANSYSModel().')
				raise exception
		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Sampling Point For each Limit State (dictionary of dictionaries)
		#	and Sampling Point properties (dic. of dics.)
		#
		#
		#
		self.SPForLS = {}
		self.SPprops = {}


		#-----------------------------------------------------------------------

		#-----------------------------------------------------------------------
		# Monte Carlo simulation controls of each LIMIT STATE
		#	(dictionary of dictionaries)
		#
		# These has one for each limit state on each cycle
		self.MCControl_LS = {}

		# Failure count for each limit state on each cycle
		self.MCControl_LS['Nfi'] = {}
		# Square Failure count for each limit state on each cycle
		self.MCControl_LS['Nfi2'] = {}
		# Cumulative Pf of each limit state on each cycle
		self.MCControl_LS['Pf'] = {}
		# CV of cumulative Pf for each cycle (Pf)
		self.MCControl_LS['CVPf'] = {}
		# Mean values of LS g(x)
		self.MCControl_LS['gMean'] = {}
		# Std.Dev values of LS g(x)
		self.MCControl_LS['gStd'] = {}
		# S1 and S2 for means and std devs (internal use)
		gS1 = {}
		gS2 = {}
		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Monte Carlo simulation controls GENERAL
		#	(dictionary of arrays)
		#
		# These has one for each limit state on each cycle
		self.MCControl = {}

		# Failure count on each cycle
		self.MCControl['Nfi'] = []
		# Square Failure count on each cycle
		self.MCControl['Nfi2'] = []
		# Cumulative Pf on each cycle
		self.MCControl['Pf'] = []
		# Reliability index Beta
		self.MCControl['Beta'] = []
		# CV of cumulative Pf for each cycle (Pf)
		self.MCControl['CVPf'] = []
		#-----------------------------------------------------------------------



		#-----------------------------------------------------------------------
		# Run all limit states looking for settings and completing some
		#	dictionaries.
		#
		# Temp sum for verifying
		sumNsi = 0
		for eachLS in self.limstates:
			# Start an empty dic.
			self.SPForLS[eachLS] = {}
			# Weigth values on W of dic.
			self.SPForLS[eachLS]['W'] = 0
			# Point on Pt of dict.
			self.SPForLS[eachLS]['pt'] = {}
			# Weight zero for initial point

			# Set Ns for each limit state on self.limstates[eachLS]
			Nsi = round(self.limstates[eachLS][1]*Ns)
			sumNsi += Nsi
			self.limstates[eachLS][3] = Nsi

			# Empty arrays
			self.MCControl_LS['Nfi'][eachLS] = []
			self.MCControl_LS['Nfi2'][eachLS] = []
			self.MCControl_LS['Pf'][eachLS] = []
			self.MCControl_LS['CVPf'][eachLS] = []
			self.MCControl_LS['gMean'][eachLS] = []
			self.MCControl_LS['gStd'][eachLS] = []

			# Start empty!
			gS1[eachLS] = 0
			gS2[eachLS] = 0

			# Run all variables declared on each limit state to get Sampling Point
			for eachVar in self.samplingDistrib[eachLS]:
				# SPForLS it's a copy of sampling distribution to keep values in case of adaptive sampling.
				self.SPForLS[eachLS]['pt'][eachVar] = list.copy(self.samplingDistrib[eachLS][eachVar])
					## self.samplingDistrib[limst][name] = [distrib, mean, std]

		# If sumNsi < Ns add simulations to last LS
		if sumNsi < Ns:
			self.limstates[eachLS][3] = Ns - (sumNsi - Nsi)
			print('Is not possible to split Ns correctly betwen current limit states, '
				+ 'last limit state received some extra simulations.')

		# If sumNsi > Ns remove simulations from last LS
		if sumNsi > Ns:
			extra = sumNsi - Ns
			self.limstates[eachLS][3] = self.limstates[eachLS][3] - extra
			print('Is not possible to split Ns correctly betwen current limit states, '
				+ 'last limit lost %d simulations.' % extra)

		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Variables and Constants lengths
		NInRandVars = len(self.variableDistrib)
		#NInConstVars = len(self.variableConst)
		#-----------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Equivalent Correlation matrix 
		
		# Create var id list
		varId = {}
		idx = 0
		# For random variables
		for eachVar in self.variableDistrib:
			varId[eachVar] = idx
			idx += 1

		for eachVar in self.variableConst:
			varId[eachVar] = idx
			idx += 1

		self.varId = varId

		# Initial correlation Matrix
		self.correlMat = np.eye(NInRandVars)

		for each in self.corlist:
			i = varId[each[0]]
			j = varId[each[1]]
			cor = self.corlist[each]

			# Apply Nataf to transform the correlation
			var1props = self.variableDistrib[each[0]]
			var2props = self.variableDistrib[each[1]]

			# Both are gauss
			if (var1props[0] is 'gauss' and var2props[0] is 'gauss'):
				cor = cor

			# Both are LN
			elif var1props[0] is 'logn' and var2props[0] is 'logn':
				cv1 = var1props[3]
				cv2 = var2props[3]
				cor = cor*(math.log(1+cor*cv1*cv2) /
				           (cor*math.sqrt(math.log(1+cv1**2)*math.log(1+cv2**2))))

			# Both are Gumbel
			elif var1props[0] is 'gumbel' and var2props[0] is 'gumbel':
				cor = cor*(1.064 - 0.069*cor + 0.005*cor**2)

			# One is gauss and other is logn
			elif (var1props[0] is 'gauss' and var2props[0] is 'logn') \
				or (var2props[0] is 'gauss' and var1props[0] is 'logn'):

				# who is logn?
				if var1props[0] is 'logn':
					cv = var1props[3]
				else:
					cv = var2props[3]

				# cor is
				cor = cor*cv/math.sqrt(math.log(1+cv**2))

			# One is gauss and other is gumbel
			elif (var1props[0] is 'gauss' and var2props[0] is 'gumbel') \
				or (var2props[0] is 'gauss' and var1props[0] is 'gumbel'):
				cor = 1.031*cor

			# One is logn and other is gumbel
			elif (var1props[0] is 'logn' and var2props[0] is 'gumbel') \
				or (var2props[0] is 'logn' and var1props[0] is 'gumbel'):
				# who is logn?
				if var1props[0] is 'logn':
					cv = var1props[3]
				else:
					cv = var2props[3]

				# cor is
				cor = cor*(1.029 + 0.001*cor + 0.014*cv + 0.004*cor**2 + 0.233*cv**2 - 0.197*cor*cv)

			# Forbiden zone
			else:
				exception = Exception('When applying NATAF on variables \"%s\" and \"%s\" the variables ' +
                                    'conditions wasn\'t possible.')
				raise exception

			# Save it!
			self.correlMat[i, j] = cor
			self.correlMat[j, i] = cor

		# Lower Cholesky matrix
		matL = np.linalg.cholesky(self.correlMat)
		matCorInv = np.linalg.inv(self.correlMat)

		#-----------------------------------------------------------------------



		#-----------------------------------------------------------------------
		#
		# 	Cycle loop's
		#
		#
		for cycle in range(1, 1+self.controls['Nmaxcycles']):
			# 1+Ncycle because Python counts from 0

			#-------------------------------------------------------------------
			# Beggining of cycle
			#
			self._PrintR('\n---\nSimulation cycle %d of %d(cycle limit).'
						 % (cycle, self.controls['Nmaxcycles']))
			#-------------------------------------------------------------------

			#-------------------------------------------------------------------
			# Generate Random values
			#

			# Store variables Values and Weights
			#	(dictionary of dictionaries)
			varsValues = {}

			# Generate random normal matrix with all values
			normalValuesMatrix = np.random.normal(0.0, 1.0, (NInRandVars, Nsi))
			# Apply correlation with x_c=L.x
			normalCorrelatedMatrix_h = matL.dot(normalValuesMatrix)
			normalCorrelatedMatrix_f = normalCorrelatedMatrix_h.copy()

			# Run for each variable each limit state, of limit state get limit
			# 	state size (Ns), if there is SamplDistrb generate Ns from
			#	sampling, else generate from RandomDistrib.
			for eachVar in self.variableDistrib:
				eachVar = eachVar.lower()
				# For each variable has values and weights
				varsValues[eachVar] = {}
				# All values as 0 in the beggining
				varsValues[eachVar]['values'] = np.zeros(Ns)
				# All weights are 1 in the beggining
				varsValues[eachVar]['weights'] = 1+np.zeros(Ns)

				# Initial position on values array (of limit state values)
				posi = 0
				# Final position on values array (of limit state values)
				posf = 0

				for eachLS in self.limstates:
					# Get posf
					Nsi = self.limstates[eachLS][3]
					posf = posi + Nsi

					# Get True Random Distribution
					varDistrib = self.variableDistrib[eachVar]

					# Verify the existence of SamplingDistrib/SPForLS, if exists get it
					if eachVar in self.SPForLS[eachLS]['pt']:
						# Have sampling distribution
						sampDist = self.SPForLS[eachLS]['pt'][eachVar]
					else:
						# Have no sampling distribution
						sampDist = False


					# Call _GenRandomVW to get values and weights
					[varsValues[eachVar]['values'][posi:posf], varsValues[eachVar]['weights'][posi:posf], \
						normalCorrelatedMatrix_f[varId[eachVar], posi:posf]] = \
						self._GenRandomVW(varDistrib, sampDist, normalCorrelatedMatrix_h[varId[eachVar], posi:posf])

					# Set posf as next posi
					posi = posf

			# Joint distributions weights
			varsValues['__joint_all_w__'] = 1+np.zeros(Ns)
			# Ratio of joint normal distributions (fX/hX) for Nataf process - if not using importance sampling it will be 1 !
			for each in range(Ns):
				# matCorInv
				Z_f = normalCorrelatedMatrix_f[:, each]
				Z_h = normalCorrelatedMatrix_h[:, each]
				varsValues['__joint_all_w__'][each] = math.exp(-1/2*(Z_f.T).dot(matCorInv.dot(Z_f))+1/2*(Z_h.T).dot(matCorInv.dot(Z_h)))

			# Now for constats variables
			for eachVar in self.variableConst:
				eachVar = eachVar.lower()
				# For each variable has values and weights
				varsValues[eachVar] = {}
				# All values as 0 in the beggining
				varsValues[eachVar]['values'] = Ns*[self.variableConst[eachVar]]
				varsValues[eachVar]['weights'] = Ns*[1]
			#-------------------------------------------------------------------


			#-------------------------------------------------------------------
			# If ANSYS is being used:
			#	Set variables values;
			#	Run ANSYS;
			#	Import results;
			#	Put results on varsValues
			#
			if self._ANSYS:
				# Initialize results dic.
				ansysRes = {}

				# Set weigths as 1 for afeter do the product of all weights of
				#	ANSYS variables
				ansysRes['weights'] = np.zeros(Ns)
				ansysRes['weights'][:] = 1

				# Runs ansys object looking for ANSYS variables
				for eachVar in self.ansys.varInNames:
					eachVar = eachVar.lower()
					# Send value to ansys object
					self.ansys.SetVarInValues(eachVar, varsValues[eachVar]['values'])

					# Join the weight of variables to ansysRes['weights']
					ansysRes['weights'] = ansysRes['weights'] * varsValues[eachVar]['weights']

				# Run ANSYS
				self.ansys.Run()

				# Import results
				ansysRes['values'] = self.ansys.GetVarOutValues()

				# Put results and weights on varsValues
				for eachVar in self.ansys.varOutNames:
					#eachVar = eachVar.lower()
					varsValues[eachVar.lower()] = {}
					varsValues[eachVar.lower()]['values'] = ansysRes['values'][eachVar.upper()]
					varsValues[eachVar.lower()]['weights'] = ansysRes['weights']

			#-------------------------------------------------------------------
			#print('ansysRes:')
			#print(ansysRes)
			#print('*varsValues:')
			#print(varsValues)



			#-------------------------------------------------------------------
			# Eval each Limit State function
			#

			self._PrintR('Evaluating limit state functions.')

			# Run all the limit states list (self.limstates[act] = [equat, weight, [], 0])
			#	Store positions of eachLS
			posi = 0
			posf = 0

			# Vectors for store failures weights
			#	Igw: = 0,	   if not failed
			#		 = weight, if failed
			Igw = np.zeros(Ns)


			for eachLS in self.limstates:
				# Get range of LS
				Nsi = self.limstates[eachLS][3]
				posf = posi + Nsi

				# For eachLS run eachSim(ulation) and evaluate equat from LS
				for eachSim in range(posi, posf):
					# Dictionary of values in this simulation and LS (used on eval)
					varVal = {}

					# Weigth of simulation (starts as 1)
					#simWei = 1.00
					simWei = varsValues['__joint_all_w__'][eachSim]

					# Run all variables from variableDistrib, if it has SamplingDistirb
					#	it will be erased in next step...
					for eachVar in self.variableDistrib:
						eachVar = eachVar.lower()
						varVal[eachVar] = varsValues[eachVar]['values'][eachSim]

					# Run all variables from ANSYSOut and SamplingDistirbs
						# Run eachVar of eachLS in eachSim getting value and weight
						#for eachVar in self.limstates[eachLS][2]:
					for eachVar in self.samplingDistrib[eachLS]:
						eachVar = eachVar.lower()
						varVal[eachVar] = varsValues[eachVar]['values'][eachSim]
						simWei = simWei * varsValues[eachVar]['weights'][eachSim]


					for eachVar in self.variableConst:
						eachVar = eachVar.lower()
						varVal[eachVar] = varsValues[eachVar]['values'][eachSim]

					#----------------------------------------------------------
					# if running ansys:
					# Attention here: the ANSYS output variables has the weight
					# of all input variables, so it mustn't be used, if used
					# some weights will be applyed 2 times!
					#
					if self._ANSYS:
						for eachVar in self.ansys.varOutNames:
							eachVar = eachVar.lower()
							varVal[eachVar] = varsValues[eachVar]['values'][eachSim]
							# DO NOT REMOVE NEXT COMMENT!
							##simWei = simWei * varsValues[eachVar]['weights'][eachSim]

					# Evaluate each simulation with its LS equation
						#Igw[eachSim] = eval(self.limstates[eachLS][0], globals(), varVal)
					# Old solution:
					#curSLSvalue = eval(self.limstates[eachLS][0], globals(), varVal)

					# Test if limstate is a string or a function
					if(type(self.limstates[eachLS][0]) is str):
						varVal['userf'] = self.limstates[eachLS][2]
						curSLSvalue = eval(self.limstates[eachLS][0], globals(), varVal)
					else:
						curSLSvalue = self.limstates[eachLS][0](**varVal)

					# Count failures
					if curSLSvalue <= 0:
						Igw[eachSim] = simWei
					else:
						Igw[eachSim] = 0

					# add curSLSvalue to gS1 and gS2 for mean and std
					gS1[eachLS] += curSLSvalue**2
					gS2[eachLS] += curSLSvalue

				# Determine current mean and std.dev for LS
				self.MCControl_LS['gMean'][eachLS].append(gS2[eachLS]/(Nsi*cycle))
				self.MCControl_LS['gStd'][eachLS].append(math.sqrt((Nsi*cycle*gS1[eachLS]-gS2[eachLS]**2)/(Nsi*cycle*(Nsi*cycle-1))))

				# Convergence for each Limit State
				# Nf and Nf**2 (weighted)
				csIgw = sum(Igw[posi:posf])
				self.MCControl_LS['Nfi'][eachLS].append(csIgw)

				csIgw2 = sum(Igw[posi:posf]**2)
				self.MCControl_LS['Nfi2'][eachLS].append(csIgw2)

				# Current Pf of limit state
				cPfi = sum(self.MCControl_LS['Nfi'][eachLS])/(cycle*Nsi)
				self.MCControl_LS['Pf'][eachLS].append(cPfi)

				# Current CVPf of Limit state
				sumNfi  = sum(self.MCControl_LS['Nfi'][eachLS])
				sumNfi2 = sum(self.MCControl_LS['Nfi2'][eachLS])
				cCVPf = 1/cPfi * 1/math.sqrt((cycle*Nsi)*(cycle*Nsi-1)) \
							 *(sumNfi2 - 1/(cycle*Nsi)*(sumNfi)**2)**0.50
				self.MCControl_LS['CVPf'][eachLS].append(cCVPf)

				#self._PrintR('**Pf=%3.8E; CVPf=%f\n' % (cPfi, cCVPf))

				# Set next LS range
				posi = posf

			# Convergence for entire simulation
			# Nf and Nf**2 (weighted)
			csIgw = sum(Igw)
			self.MCControl['Nfi'].append(csIgw)

			csIgw2 = sum(Igw**2)
			self.MCControl['Nfi2'].append(csIgw2)

			# Current Pf of limit state
			cPfi = sum(self.MCControl['Nfi'])/(cycle*Ns)
			self.MCControl['Pf'].append(cPfi)

			# Current Beta
			self.MCControl['Beta'].append(-scipy.stats.norm.ppf(cPfi))

			# Current CVPf of Limit state
			sumNfi  = sum(self.MCControl['Nfi'])
			sumNfi2 = sum(self.MCControl['Nfi2'])
			cCVPf = 1/cPfi * 1/math.sqrt((cycle*Ns)*(cycle*Ns-1)) \
							*(sumNfi2 - 1/(cycle*Ns)*(sumNfi)**2)**0.50
			self.MCControl['CVPf'].append(cCVPf)

			########################################
			# It was used for debug/control
			#print('varsValues:')
			#print(varsValues)
			#print('Igw:')
			#print(Igw)
			#print('csIgw:')
			#print(csIgw)
			#print('csIgw2:')
			#print(csIgw2)
			#print('cPfi:')
			#print(cPfi)
			#print('cCVPf:')
			#print(cCVPf)
			#print('---------------')

			# Print main results
			self._PrintR('\nSolution on Cycle %d with %3.3E simulations:' % (cycle, cycle*Ns))
			self._PrintR('Pf=%2.4E \nBeta=%2.3f \nCVPf=%2.3f \n' % (cPfi, -scipy.stats.norm.ppf(cPfi), cCVPf))
			
			# Print limit states mean and std dev
			self._PrintR('Limit states means and standard deviations:')
			for eachLS in self.limstates:
				self._PrintR('  Limit state %d: mean=%.4E, std.dev.=%.4E.' % (eachLS, self.MCControl_LS['gMean'][eachLS][-1], self.MCControl_LS['gStd'][eachLS][-1]))

			# Verify the convergence criteria for CVPf after 3rd cycle
			# Avoid CVPf < 1E-5!
			if cCVPf <= self.controls['CVPf'] and cycle > 3 and cCVPf > 1E-5:
				self._PrintR('CVPf convergence criteria reached on cycle %d with %3.3E simulations.' % (cycle, cycle*Ns))
				self._PrintR('Finalizing process.')
				# Set status as 0 (no problem)
				stnumb = 0
				break


			#-------------------------------------------------------------------



			#-------------------------------------------------------------------
			#	Find new design point
			#

			if adapt is True:
				self._PrintR('Evaluating new design point.')
				# Initial Igw position
				posi = 0

				# Initial max relative error is 0, in the end will have the max
				#	value from adaption
				maxrelerror = 0

				# eachLS has its sampling distrib
				for eachLS in self.limstates:
					# Get range of LS
					Nsi = self.limstates[eachLS][3]
					# Final Igw position of LS
					posf = posi + Nsi

					# Copy of current sample variables (that will be the old sampling)
					#oldSamp = {}
					#oldSamp = self.SPForLS[eachLS]['pt']

					# Old and Current simualtions weights
					oldW = self.SPForLS[eachLS]['W']
					curW = sum(Igw[posi:posf])

					# Calibrate one-by-one each variable of each LS
					self._PrintR('New design point for Limit State %d:' % eachLS)
					for eachVar in self.SPForLS[eachLS]['pt']:
						# Old and Current simulations point coord * simulations weights
						oldPW = (self.SPForLS[eachLS]['pt'][eachVar][1] * self.SPForLS[eachLS]['W'])
						curPW = sum(varsValues[eachVar]['values'][posi:posf] * Igw[posi:posf])

						# New point
						NP = (oldPW+curPW)/(oldW+curW)

						# Print
						self._PrintR(' %s = %3.5E' % (eachVar, NP))

						# Relative error of current point abs(new-old)/new
						relerror = abs(NP-self.SPForLS[eachLS]['pt'][eachVar][1])/NP
						# The biggest will be maxrelerror
						maxrelerror = max(maxrelerror, relerror)

						# Save New Point
						self.SPForLS[eachLS]['pt'][eachVar][1] = NP

					# Future old weight
					self.SPForLS[eachLS]['W'] = oldW + curW

				# Print max point error of LS
				self._PrintR('Max. relative error on sampling point search is %f.' % maxrelerror)
				if maxrelerror <= self.controls['tolAdPt'] and cycle > 3:
					self._PrintR('Sampling point search converged on cycle %d.' % cycle)
					adapt = False

					# Current LS Igw is Igw[posi:posf]
					#Igw[posi:posf]
					# loop on eachvar of eachlimstate
					# X_i+1,j=sum(X_i,j*W_i)/sum(W_i,j)


		#-------------------------------------------------------------------


		#-------------------------------------------------------------------
		# After CYCLEs
		#
		timef = time.time()
		# Save last cycles
		self.cycles = cycle
		self.MCControl['ElapsedTime'] = ((timef-timei)/60)

		self._PrintR('\n\n=======================================================================\n')
		# Verify if stopped with no convergence
		if cycle == self.controls['Nmaxcycles']:
			stnumb = 1
			self._PrintR(' WARNING:')
			self._PrintR('  The process was finished after reach the limit of cycles without')
			self._PrintR('  reach the expected CVPf.\n')

		self._PrintR(' Total of simulations: %3.3E' % (Ns*cycle))
		self._PrintR(' Probability of failure (Pf): %2.4E' % (self.MCControl['Pf'][-1]))
		self._PrintR(' Reliability index (Beta): %2.3f' % (self.MCControl['Beta'][-1]))
		self._PrintR(' CV of Prob. of failure (CVPf): %2.3f' % (self.MCControl['CVPf'][-1]))
		self._PrintR(' Elapsed time: %f minutes.' % (self.MCControl['ElapsedTime']))
		# Print limit states mean and std dev
		self._PrintR(' Final Limit States means and standard deviations:')
		for eachLS in self.limstates:
			self._PrintR('   Limit state %d: mean=%.4E, std.dev.=%.4E.' % (eachLS, self.MCControl_LS['gMean'][eachLS][-1], self.MCControl_LS['gStd'][eachLS][-1]))
		self._PrintR('\n=======================================================================\n\n')
		#-------------------------------------------------------------------


		#-----------------------------------------------------------------------
		# Send the return
		#

		# rsampdict is the dictionary with sampling points in each limit state
		rsampdict = {}
		# distparms is the dictionary with g(X) mean and std dev for each ls
		distparms = {}
		for eachLS in self.limstates:
			rsampdict[eachLS] = {}
			for eachVar in self.SPForLS[eachLS]['pt']:
				rsampdict[eachLS][eachVar] = self.SPForLS[eachLS]['pt'][eachVar][1]
			distparms[eachLS] = {}
			distparms[eachLS]['gMean'] = self.MCControl_LS['gMean'][eachLS][-1]
			distparms[eachLS]['gStd'] = self.MCControl_LS['gStd'][eachLS][-1]

		# put all in a dict
		self._stnumb = stnumb

		finret = {}
		finret['status'] = self._stnumb
		finret['Pf'] = self.MCControl['Pf'][-1]
		finret['Beta'] = self.MCControl['Beta'][-1]
		finret['CVPf'] = self.MCControl['CVPf'][-1]
		finret['SamplingPoints'] = rsampdict
		finret['cycles'] = cycle
		finret['distparms'] = distparms

		return finret

		#-----------------------------------------------------------------------


	def GetSolutionControl(self, thing):
		"""
		Return values of Monte Carlo solution controllers.

		Parameters
		----------
		thing: str, obligatory
			Control that will be returned. Available things are listed below.

		Available things
		----------------
		In function of N:
			* 'N_Pf' = Probability of failure
			* 'N_Beta' = Reliability index
			* 'N_CVPf' = CV of Probability of failure

		Returns
		-------
		2D numpy array of floats:
			Each line has simulation number and requested value on this simulation.

		"""

		# Number of simulations on each cycle and cycles
		Ns = self.controls['Ns']
		cycles = self.cycles


		if thing is 'N_Pf':
			# N vs Pf
			result = np.zeros([cycles, 2])
			result[:, 0] = range(Ns, Ns*(cycles+1), Ns)
			result[:, 1] = self.MCControl['Pf']

		elif thing is 'N_Beta':
			# N vs Beta
			result = np.zeros([cycles, 2])
			result[:, 0] = range(Ns, Ns*(cycles+1), Ns)
			result[:, 1] = self.MCControl['Beta']

		elif thing is 'N_CVPf':
			# N vs CVPf
			result = np.zeros([cycles, 2])
			result[:, 0] = range(Ns, Ns*(cycles+1), Ns)
			result[:, 1] = self.MCControl['CVPf']

		else:
			exception = Exception('Error while getting values of Monte Carlo solution control. '+
				'"%s" is not recognized as a valid "thing".' % (thing))
			raise exception

		# return
		self._PrintR('Returning values of "%s".' % thing)
		return result

	def ExportDataCSV(self, filename, description=None):
		"""
		Exports Simulation data to a CSV file.

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
			exception = Exception('Unable to open the %s file for write data.' % filename)
			raise exception
		else:
			# Starts with sep=, for Microsoft Excel
			f.write('sep=,\n')

			# Description
			if description is not None:
				f.write('%s\n\n' % description)

			f.write('Input data:\n')

			# Simulation Controllers:
			f.write('Simulation Controllers:\n')
			f.write(',Ns/Cycle:,%d\n' % self.controls['Ns'])
			f.write(',MaxCycles:,%d\n' % self.controls['Nmaxcycles'])
			f.write(',CVPf target:,%2.4f\n' % self.controls['CVPf'])
			f.write(',tol. Adapt.:,%s\n' % str(self.controls['tolAdPt']))
			f.write('\n')

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
			f.write(',Name, Distribution, Mean, Standard Deviation, CV, Par1, Par2\n')
			for eachVar in self.variableDistrib:
				values = self.variableDistrib[eachVar]
				cmd = ',%s,%s' % (eachVar, values[0])
				for eachVal in values[1:]:
					cmd = '%s,%f' % (cmd, eachVal)
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
				# WTF DID I DO? But it works very well... xD
				# GO HORSE! GO!!
				cmd = ',%s' % list(self.varId.keys())[idx]
				for eachVal in eachLine:
					cmd = '%s,%f' % (cmd, eachVal)
				cmd = '%s\n' % cmd
				f.write(cmd)
				idx += 1
			f.write('\n')

			# Limit states and theirs Sampling Distributions
			f.write('Limit States:\n')
			for eachLS in self.limstates:
				f.write(',Limit State %d:,"%s"\n' % (eachLS, self.limstates[eachLS][0]))
				f.write(',,Weight:,%f\n' % self.limstates[eachLS][1])
				f.write(',,Ns/Cycle:,%f\n' % self.limstates[eachLS][3])
				f.write(',,Starting Sampling Distributions:\n')
				f.write(',,,Name, Distribution, Mean, Standard Deviation, CV, Par1, Par2\n')
				for eachVar in self.samplingDistrib[eachLS]:
					values = self.samplingDistrib[eachLS][eachVar]
					cmd = ',,,%s,%s' % (eachVar, values[0])
					for eachVal in values[1:]:
						cmd = '%s,%f' % (cmd, eachVal)
					f.write('%s\n' % cmd)
				f.write('\n')

			# Results part
			f.write('\nResults:\n')

			# Final results
			f.write('Monte Carlo results:\n')
			f.write(',Exit status:,%d,\n' % self._stnumb)
			f.write(',Total of simulations:,%3.3E\n' % (self.cycles*self.controls['Ns']))
			f.write(',Probability of failure (Pf):,%2.4E\n' % self.MCControl['Pf'][-1])
			f.write(',Reliability Index (Beta):,%2.3f\n' % self.MCControl['Beta'][-1])
			f.write(',CV of Prob. of failure (CVPf):,%2.3f\n' % self.MCControl['CVPf'][-1])
			f.write(',Elapsed time (minutes):,%4.3f\n' % self.MCControl['ElapsedTime'])
			f.write('\n')

			# Final Sampling Point
			if self.controls['tolAdPt'] is not False:
				f.write('Final sampling point:\n')
				for eachLS in self.SPForLS:
					f.write(',Limit State %d:\n' % (eachLS))
					for eachVar in self.SPForLS[eachLS]['pt']:
						f.write(',,%s,%f\n' % (eachVar, self.SPForLS[eachLS]['pt'][eachVar][1]))
					f.write('\n')

			# Convergence on each Cycle
			f.write('\nProcess Convergence:\n')
			result = np.zeros([self.cycles, 4])
			result[:, 0] = range(self.controls['Ns'], self.controls['Ns']*(self.cycles+1), self.controls['Ns'])
			result[:, 1] = self.MCControl['Pf']
			result[:, 2] = self.MCControl['Beta']
			result[:, 3] = self.MCControl['CVPf']

			f.write(',N,Pf,Beta,CVPf\n')
			for eachLine in result:
				f.write(',%d, %2.4E, %2.3f, %2.3f\n' % \
					(eachLine[0],eachLine[1],eachLine[2],eachLine[3]))
			# End
			f.close()

			self._PrintR('Simulation data exported to "%s".' % filename)


	def Graph(self, things, show=True, savefile=False):
		"""
		Generate graphics of Monte Carlo solution controllers.

		Things can be a list of data that will be ploted in the same figure, the
		figure doesn't need to be opened, it could be just saved, or just opened.

		Parameters
		----------
		things : list of strings, obligatory
			List of data that will be ploted in the same figure. Available things
			are listed below.

		show : bool, optional
			Sinalize if figure should be opened.
			Defaults to True.

		savefile : str/bool, optional
			If it's False doesn't save anything.
			If it's a string it will be used as directory+name that figure will
			have, it shouldn't have extension, since it will be SVG.
			Defaults to False.


		Available things
		----------------
		With N as horizontal axis:
			* 'N_Pf' = Probability of failure
			* 'N_Beta' = Reliability index
			* 'N_CVPf' = CV of Probability of failure

		"""
		import matplotlib.pyplot as plt

		# Verify if things is a list, if isn't put on a list
		if type(things) == list:
			nplots = len(things)
			things = things
		else:
			nplots = 1
			things = [things]

		idx = 0

		# Create figure and plots
		plots = plt.subplots(nplots, 1)
		plt.figure(num=1, figsize=[4,8])
		plt.subplots_adjust(hspace = 0.45, left = 0.15, right = 0.95, \
							top = 0.95, bottom = 0.12)

		# Now run all on the list
		for eachPlot in things:
			idx += 1
			if eachPlot == 'N_CVPf':
				xlab = 'Simulations'
				ylab = 'CVPf'
			elif eachPlot == 'N_Pf':
				xlab = 'Simulations'
				ylab = 'Probability of failure'
			elif eachPlot == 'N_Beta':
				xlab = 'Simulations'
				ylab = r'Reliability index ($\beta$)'
			else:
				exception = Exception('Error while getting values of Monte Carlo solution control. '+
					'"%s" is not recognized as a valid "thing".' % (eachPlot))
				raise exception

			if nplots > 1:
				plots[1][idx-1].plot(self.GetSolutionControl(eachPlot)[:,0], self.GetSolutionControl(eachPlot)[:,1])
				plots[1][idx-1].set_ylabel(ylab)
				plots[1][idx-1].grid(True)
				# Just last has xlabel
				plots[1][-1].set_xlabel(xlab)
			else:
				plots[1].plot(self.GetSolutionControl(eachPlot)[:,0], self.GetSolutionControl(eachPlot)[:,1])
				plots[1].set_ylabel(ylab)
				plots[1].grid(True)
				plots[1].set_xlabel(xlab)

		# Save, or not
		if savefile is not False:
			file = savefile+'.svg'
			self._PrintR('Saving figure as "%s".' % file)
			plt.savefig(file, transparent=True)

		# Show, or not
		if show is True:
			plt.show()
