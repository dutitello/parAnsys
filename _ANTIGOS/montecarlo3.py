# -*- coding: UTF-8 -*-
"""
This class uses the ANSYS class to performn a Monte Carlo analysis on ANSYS.


"""
import os
from paransys.ansys import ANSYS

# child class MonteCarloANSYS: all the commands from class ANSYS if not declared
# here will be imported from there if called outside. Here can be called by
# ANSYS.command(self,...)

class MonteCarloANSYS(object):
	"""
	This is a child class of ANSYS, they are really similar.
	It's __init__ is the same of ANSYS class, an most commands are the same,
	the difference here are some commands specified for Monte Carlo analysis.

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

	def __init__(self, exec_loc=None, run_location=os.getcwd(), jobname='file',nproc=2, override=False, cleardir=False, add_flags=''):
		"""
		Define the properties of ANSYS execution.
		(As defined in the main class ANSYS)
		"""

		# Create ANSYS object
		self.ansys = ANSYS(exec_loc, run_location, jobname, nproc, override, cleardir, add_flags)




	# Here the values from Input variables aren't set, their distribuiton is set.

	# Avoid to try set the values by parent command SetVarInValues()
	def SetVarInValues(self, name, values):
		"""
		This command couldn't be used with Monte Carlo.
		"""
		exception = Exception('Using class MonteCarloANSYS the values are not' +
							  'set by user, instead use SetVarInDistrib() to '+
							  'set variable distribution.')
		raise exception

	def SetVarInDistrib(self, name, distrib, mean, cv, par1=None, par2=None):
		"""
		Sets the distribution of a variable to performn the simulations.

		Parameters:
		---
		name : str, obligatory
			Name of variable declared in CreateVarIn()

		distrib : str, obligatory
			Probabilistic variable distribution type.

			Available types are:
			- Gaussian (or gauss, normal)
			- Lognormal (or log, logn, lognorm)
			- Gumbel (or gumb, type1)
			- Constant (or const) - Constant value (not random!)
			-


		"""

		# Verify the name
		if name in ['', None, ' ']:
			exception = Exception('You must define the name of the variable that will receive the values.')
			raise exception
		else:
			name = name.upper()
			if name not in self.varInNames:
				exception = Exception('This variable name is not declared.\n'+
								'Please use CreateVarIn("name") to declare it.')
				raise exception

		# Verify the distribution and determine the parameters
