# -*- coding: UTF-8 -*-
"""
This class uses the ANSYS class to performn a Monte Carlo analysis on ANSYS.


"""
import os
from paransys.ansys import ANSYS

# child class MonteCarloANSYS: all the commands from class ANSYS if not declared
# here will be imported from there if called outside. Here can be called by
# ANSYS.command(self,...)

class MonteCarloANSYS(ANSYS):
	"""
	This is a child class of ANSYS, they are really similar.
	It's __init__ is the same of ANSYS class, an most commands are the same,
	the difference here are some commands specified for Monte Carlo analysis.

	CreateVarIn, SetVarInDistrib, SetVarInSampl, CreateVarOut are related to the
	parameters for ANSYS analysis.


	This code was made following the idea of ANSYS be an tool for getting the
	ultimate load of a structure. To do that the APDL script apply a big
	displacement in the structure and then put the maximum reaction force on the
	displaced node in the PARAMETER variable that will be used as control here.
	The real loads for reliability test are generated inside Python and compared
	with the ultimate load of structure (G=R-S).

	When applying displacement and getting maximum force on the structure
	soliciting force is generated on


	Parameters: (same as ANSYS)
	---
	exec_loc : str, obligatory
		Location of ANSYS executable file.

	run_location : str, optional
		ANSYS working directory. Recomended to be a separated directory.
		Defaults to current directory.

	jobname : str, optional
		ANSYS jobname. Defaults to 'file'.

	nproc : int, optional
		Number of processors. Defaults to 2.

	override : bool, optional
		Attempts to delete the *.lock file at working directory.
		It's useful when ANSYS was interrupted.
		Defaults to False

	cleardir : bool, optional
		Delete all the files from ANSYS working directory when call the Run command.
		Defaults to False

	add_flags : str, optional
		Additional flags to be called with ANSYS.
		If it's an academic version use add_flags='-aa_r'
		Do not use '-b -i -o'
		Flags can be found at https://www.sharcnet.ca/Software/Ansys/16.2.3/en-us/help/ans_ope/Hlp_G_OPE3_1.ht

	"""

	#---------------------------------------------------------------------------
	# Definition of self varibles that just exists on MonteCarlo
	# Dictionary with variables distrib and parameters
	variabledistrib = {}
	samplingdistrib = {}
	controls = {}
	ansysvars = []
	limstates = {}



	# Here the values from Input variables aren't set, their distribuiton is set.
	# Avoid to set the values by parent command SetVarInValues()
	def SetVarInValues(self, name, values):
		"""
		This command couldn't be used with Monte Carlo.
		"""
		exception = Exception('Using class MonteCarloANSYS the values are not' +
							  'set by user, instead use SetVarInDistrib() to '+
							  'set variable distribution.')
		raise exception


	# Internal for variable distrib definition
	def _VarInDistrib(self, type, name, distrib, mean, cv, par1, par2, limst=0):
		"""
		Internal function that verify the values for SetVarInDistrib
		and SetVarInSampl.

		Parameters:
		---
		type : int, obligatory
			type = 0: for SetVarInDistrib -> variabledistrib
			type = 1: for SetVarInSampl   -> samplingdistrib

		limst : integer, obligatory for type=1 with more than 1 limit state
			Limit state ID that will use current sampling distribution.

		"""
		# Verify the variable name
		if name in ['', None, ' ']:
			exception = Exception('You must define the name of the variable that will receive the values.')
			raise exception
		else:
			name = name.upper()
			if name not in self.varInNames:
				exception = Exception('This variable name is not declared.\n'+
								'Please use CreateVarIn() to declare it.')
				raise exception


		# Set distribution name as lower case
		distrib = distrib.lower()


		# Verify the distribution and then determine the parameters

		# Gaussian variable
		if distrib in ['gauss', 'gaus', 'gaussian', 'normal']:
			# If mean=0 the standard deviation will be 0 to, to avoid this when
			#	mean=0 CV will be the standard deviation.
			if mean is 0:
				std = cv
			else:
				std = mean * cv


			# Type 0 = Random Variable distribution
			if type is 0:
				self.variabledistrib[name] = [distrib, mean, std]
				return self._PrintR('Variable %s defined as Gaussian with mean=%f and std. dev.=%f'
									% (name, mean, std))

			# Type 1 = Sampling distribution
			elif type is 1:
				self.samplingdistrib[name] = [limst, distrib, mean, std]
				return self._PrintR('Variable %s sampling defined as Gaussian with mean=%f and std. dev.=%f for limit state %d.'
									% (name, mean, std, limst))

		# Lognormal distribution
		elif distrib in ['lognormal', 'logn', 'ln', 'log', 'lognorm']:
			# Distribution parameters are qsi, related do std dev, and
			#	lmbd, related to mean.
			qsi = math.sqrt(math.log(1 + (cv)**2))
			lmbd = math.log(mean) - 0.5*qsix**2



			# Type 0 = Random Variable distribution
			if type is 0:
				self.variabledistrib[name] = [distrib, lmbd, qsi]
				return self._PrintR('Variable %s defined as LogNormal with lmbd=%f and qsi=%f'
									% (name, lmbd, qsi))

			# Type 1 = Sampling distribution
			elif type is 1:
				self.samplingdistrib[name] = [limst, distrib, lmbd, qsi]
				return self._PrintR('Variable %s sampling defined as LogNormal with lmbd=%f and qsi=%f for limit state %d.'
									% (name, lmbd, qsi, limst))


		# Gumbel distribution
		elif distrib in ['gumbel', 'gumb', 'type1']:
			pass

		# Constant value
		elif distrib in ['constant', 'const', 'cons', 'c']:
			# Type 0 = Random Variable distribution
			if type is 0:
				self.variabledistrib[name] = [distrib, mean]
				return self._PrintR('Variable %s defined as Constant with value=%f'
									% (name, mean))

			# Type 1 = Sampling distribution
			elif type is 1:
				self.samplingdistrib[name] = [limst, distrib, mean]
				return self._PrintR('Variable %s sampling defined as Constant with value=%f for limit state %d.'
									% (name, mean, limst))

		else:
			exception = Exception('Distribution %s set on variable %s is not recognized.'
									% (distrib.upper(), name))
			raise exception



	# Setting distribution of variables
	def SetVarInDistrib(self, name, distrib, mean, cv=0, par1=None, par2=None):
		"""
		Sets the distribution of a variable to performn the simulations.

		Parameters:
		---
		name : str, obligatory
			Name of variable declared in CreateVarIn()

		distrib : str, obligatory
			Probabilistic variable distribution type.
			For all distributions Mean and CV are related to Normal distribution
			(the code determines the parameters for the desired distribution)

			Available types are:
			- gaussian (or gauss, normal) (If mean=0 CV will be Stand. Deviation)
			- lognormal (or log, logn, ln, lognorm)
			- gumbel (or gumb, type1)
			- constant (or const) - Constant value (doesn't need CV)
			-

		mean : float, obligatory
			Standard mean of variable values.

		cv : float, obligatory
			Coefficient of variation* of variable values.
			(CV = StandardDeviation/Mean)

			* If mean is 0 the value set on CV will be the Standard Deviation!

		par1 and par2 : float, optional
			Parameters for future implementations.

		"""

		self._VarInDistrib(type=0, name=name, limst=None, distrib=distrib, mean=mean, cv=cv, par1=par1, par2=par2)



	# Sampling definition
	def SetVarInSampl(self, name, limst, distrib, mean, cv=0, par1=None, par2=None):
		"""
		Sets the sampling distribution of a variable to performn Importance
		Sampling the simulations.

		Parameters:
		---
		name : str, obligatory
			Name of variable declared in CreateVarIn()

		limst : integer, obligatory for type=1 with more than 1 limit state
			Limit state ID that will use current sampling distribution.

		distrib : str, obligatory
			Probabilistic variable distribution type.
			For all distributions Mean and CV are related to Normal distribution
			(the code determines the parameters for the desired distribution)

			Available types are:
			- gaussian (or gauss, normal) (If mean=0 CV will be Stand. Deviation)
			- lognormal (or log, logn, ln, lognorm)
			- gumbel (or gumb, type1)
			- constant (or const) - Constant value (doesn't need CV)

		mean : float, obligatory
			Standard mean of variable values.

		cv : float, obligatory
			Coefficient of variation* of variable values.
			(CV = StandardDeviation/Mean)

			* If mean is 0 the value set on CV will be the Standard Deviation!

		par1 and par2 : float, optional
			Parameters for future implementations.

		"""
		self._VarInDistrib(type=1, name=name, limst=limst, distrib=distrib, mean=mean, cv=cv, par1=par1, par2=par2)



	def SetLength(self, length):
		"""
		This command couldn't be used with Monte Carlo.
		"""
		exception = Exception('When using class MonteCarloANSYS the Length of '+
							  'data is defined in the controls of simulation.'+
							  'Instead of SetLength() use SetControls().')
		raise exception


	def SetControls(self, Ns, Nmaxcycles, cvpf=0.00):
		"""
		Set the control parameters of simulation.

		Parameters:
		---
		Ns : integer, obligatory
			Number of simulations performed on each cycle.
			After each cycle the convergence of simualtion is verified.
			When using Importance Sampling with Adaptive Sampling, after each
			cycle the new sampling point will be determined.

		Nmaxcycles : integer, obligatory
			Maximum number of cycles to be performed, if cvpf is not reached on
			Nmaxcycles the simulation will be interrupted.

		cvpf : float, optional
			Target value of Probability Failure Coefficient of Variation, when
			reached the simulation stops.

		"""
		if Ns >= 0 and Nmaxcycles >= 0 and cvpf >= 0:
			self.length = Ns
			self.controls['Ns'] = Ns
			self.controls['Nmaxcycles'] = Nmaxcycles
			self.controls['cvpf'] = cvpf

			return _PrintR('Simulation control set. Ns=%d, NMaxCycles=%d, NMax=%d, CVPf_targ=%f'
							% (Ns, Nmaxcycles, Ns*Nmaxcycles, cvpf))
		else:
			exception = Exception('Error while setting simulation controls. Please verify the set values.')
			raise exception


	def SetLimitState(self, equat, varlist, weight=1.00):
		"""
		Create and Set a new limit state
		The number ID of LimitStates are generated automatically starting at 0
		for the first.

		Parameters:
		---
		equat : str, obligatory
			String with the equation of the limit state. It must be write as a
			function of defined variables (In and Out).

		varlist : list of strings, obligatory
			A list with the name of variables used in equat

		weight : float, obligatory only with more than 1 limit state
			The weight of current limit state, it determines how the simulations
			are distributed betwen all the limit states.
			The sum of all limit states must be 1.00, so, if there is just one
			limit state it's weight should be 1.00

		For example, if ANSYS returns the maximum load on a truss as variable
		FxMAX, and applied loads to be tested are '(g+q)*sin(theta)', where
		g, q, theta are defined variables created with CreateVarIn(...), and
		this is the unique limite state for the analysis:
		SetLimitState(equat='FxMAX-(g+q)*sin(theta)', varlist=['FxMAX', 'g', 'q', 'theta'], weight=1.00)

		"""
		if weight <= 0 or weight > 1:
			exception = Exception('The weigth of limit state must be greater than zero and less equal to 1.00.')
			raise exception

		# Add current limit state to the list
		act = len(limstates)
		limstates[act] = [equat, varlist, weight]

		# Verify the sum off all limit states until now
		#	itsn't possible to verify the equantion since we dont have ANSYS out values
		tempsum = 0
		for each in limstates:
			tempsum += limstates[each][2]

		if tempsum > 1.00:
			exception = Exception('The sum of limit states weights exceeded 1.00, please verify the inserted values.')
			raise exception

		return _PrintR('Limit state %d defined, the sum of weights untill now is %f.' %(each, act))


	def SetANSYSVar(self, name):
		"""
		Set an Input variable as ANSYS variable

		Parameters:
		---
		name : str, obligatory
			Name of variable.

		"""

		# Verify the variable name
		if name in ['', None, ' ']:
			exception = Exception('You must define the name of the variable.')
			raise exception
		else:
			name = name.upper()
			if name not in self.varInNames:
				exception = Exception('This variable name is not declared.\n'+
								'Please use CreateVarIn() to declare it.')
				raise exception

		# If exists append to ansysvars
		ansysvars.append(name)


	def Run(self):

		# Performn the conditional clear
		self._ClearForRun()

		# Write the pdsrun.inp file
		self._writePDSfile()

		# Write the sample file
		self._writeSAMPfile()

		# Run
		self._Run()


	#def SetVarIn
