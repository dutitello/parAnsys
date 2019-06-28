# -*- coding: UTF-8 -*-
"""
This module calls ANSYS by Probabilistic Design System (PDS) as an FEM tool for
running a couple of parametric analysis and get some variables back to Python.

Built using ANSYS documentation from https://www.sharcnet.ca/Software/Ansys/

This module makes no claim to own any rights to ANSYS, it's just an interface
to call the program owned by ANSYS.

"""

import os
import numpy as np
import time


class ANSYS(object):
	"""
	This class creates a couple of script files with some parameters defined here
	by the user, copy an APDL script file created by the user with the model to be
	analysed, run everything on ANSYS and get the results from some defined
	parameters back to Python.

	Parameters
	----------
	exec_loc : str, obligatory
		Location of ANSYS executable file.

	run_location : str, optional
		ANSYS working directory. Recomended to be a separated directory.
		Defaults to ansys_anl on current directory.

	jobname : str, optional
		ANSYS jobname. Defaults to \'file\'.

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


	def __init__(self, exec_loc=None, run_location=os.getcwd()+'\\ansys_anl\\', jobname='file',
			     nproc=2, override=False, cleardir=False, add_flags=''):
		"""
		Define the properties of ANSYS exection.

		(As defined in the main class ANSYS)
		"""

		# Verify if the exec_loc is defined and real
		if exec_loc is None:
			exception = Exception('Undefined ANSYS executable location. \n You must define ANSYS executable location.')
			raise exception
		else:
			# Verify if the file exists
			if not os.path.isfile(exec_loc):
				exception = Exception('Invalid ANSYS executable file as \"%s\"' % exec_loc)
				raise exception

		# Save the options to the main object
		self.ANSYSprops = {}
		self.ANSYSprops['exec_loc']  = exec_loc
		self.ANSYSprops['run_location'] = run_location
		self.ANSYSprops['jobname']   = jobname
		self.ANSYSprops['nproc']     = nproc
		self.ANSYSprops['override']  = override
		self.ANSYSprops['cleardir']  = cleardir
		self.ANSYSprops['add_flags'] = add_flags

		# Create lists and set initial values

		# Variables names
		self.varInNames = []
		self.varOutNames = []

		# Variables values
		self.varInValues = {}
		self.varOutValues = {}

		# Analysis length
		self.length = 0

		# Model properties
		self.Model = {}

		# Defalts to print always
		self.PrintR = True


		self._PrintR('ANSYS properties defined as:')
		self._PrintR('   Executable file: \"%s\".' % self.ANSYSprops['exec_loc'])
		self._PrintR('   Working directory: \"%s\".' % self.ANSYSprops['run_location'])
		self._PrintR('   Jobname: \"%s\".' % self.ANSYSprops['jobname'])
		self._PrintR('   Number of processors used: \"%s\".' % self.ANSYSprops['nproc'])
		self._PrintR('   Override lock file: \"%s\".' % self.ANSYSprops['override'])
		self._PrintR('   Clear working directory: \"%s\".' % self.ANSYSprops['cleardir'])
		self._PrintR('   Additional flags: \"%s\".' % self.ANSYSprops['add_flags'])


	def Info(self, act=False):
		"""
		Turn on/off the return of the commands to Python.

		Parameters
		----------
		act : bool, obligatory
			True turn On and False turn Off the return of the commands to Python.
		"""

		self.PrintR = act
		self._PrintR('Now the commands will send a return to Python (like this).')

	def _PrintR(self, value):
		"""
		Internal function to print or not the return of commands based on command Info()
		"""
		if self.PrintR:
			print(value, flush=True)
		else:
			pass


	def _writePDSfile(self):
		"""
		Internal function to create/write the APDL file with the PDS analysis. (pdsrun.inp)
		"""
		# Verify if self.length >0
		if self.length <= 0:
			exception = Exception('Before run you must define the length of the analysis with SetLength()'+
									'and set the all the variables values.')
			raise exception

		# Run all the varInValues[var] looking for an empty
		for each in self.varInNames:
			if self.varInValues[each] is None:
				exception = Exception('Input variable \"%s\" has no defined values.' % each)
				raise exception

		# Try to open the file for writing
		try:
			runfile = '%s\\pdsrun.inp' % self.ANSYSprops['run_location']
			f = open(runfile, 'wt')
		except:
			exception = Exception('Unable to open the pdsrun.inp file for writting.')
			raise exception

		# Now try to write
		try:
			self._PrintR('Writing the pdsrun.inp file.')
			# some ANSYS initial commands
			f.write('FINISH\n')
			f.write('/GOPR\n')
			f.write('KEYW,PR_SET,1\n')
			f.write('KEYW,PR_STRUC,1\n')

			# Open PDS module
			f.write('/PDS\n')

			# Clear the PDS database
			f.write('PDCLR,ALL\n')

			# Here we put the from SetModel(..) that is named here as "current.inp"
			f.write('PDANL,current,inp\n')

			# Input variables
			for each in self.varInNames:
				# Variables are declared with uniform distribution betwen 1.5*maxval and 0.5*minval
				cmd = 'PDVAR,%s,UNIF,%f,%f\n' % (each, 0.5*min(self.varInValues[each]), 1.5*max(self.varInValues[each]))
				f.write(cmd)

			# Output/control variables
			for each in self.varOutNames:
				# Just the varname and distrib=RESP
				cmd = 'PDVAR,%s,RESP\n' % (each)
				f.write(cmd)

			# We will use Monte Carlo simulation with user sampling, that way we can
			# simulate the points declared in each variable easily and without problems.
			f.write('PDMETH,MCS,USER\n')
			# The sample points file is named here as "current.samp"
			f.write('PDUSER,current,samp\n')

			# Runs the Analysis - The analysis is named as current so the results file
			# will be $jobname%_current.pdrs
			f.write('PDEXE,current\n')

			# Create a delay of 2s on ANSYS just to make sure that all were done
			f.write('/WAIT,2\n')

			# Close the file
			self._PrintR('Saving the pdsrun.inp file.')
			f.close()

		except:
			exception = Exception('Error while writting the pdsrun.inp file.')
			raise exception
		else:
			pass



	def _writeSAMPfile(self):
		"""
		Internal function to create/write the sample points file to run the PDS. (current.samp)
		The sample file format has the PDEXE name in the first line, and the second line is:

		ITER CYCL LOOP $INPUT VAR NAMES$

		With a new column for each variable.
		The values for ITER and CYCL are always 1, LOOP is the number of solution range(0 to length)
		"""

		#
		self._PrintR('Writing the current.samp file.')

		# Create initial the header string
		header = 'current\nITER CYCL LOOP'

		# Create the initial format string
		format = '%d %d %d'

		# Create initial values array with ones
		samparray = np.zeros([self.Length, 3+len(self.varInNames)]) + 1

		# Fill the loop column
		samparray[:, 2] = range(1, self.Length+1)

		# Each variable will add his name at header and a %15.8E in the format
		# A counter from 2 (3rd column)
		i = 2
		for each in self.varInNames:
			# Add strings
			header = '%s %s' % (header, each)
			format = '%s %s' % (format, '%15.8E')

			# Add one in counter
			i += 1

			# Place values
			try:
				samparray[:, i] = self.varInValues[each]
			except:
				exception = Exception('Error while passing the values of \"%s\" to the current.samp file.' % each)
				raise exception

		# Save the file
		try:
			sampfile = '%s\\current.samp' % self.ANSYSprops['run_location']
			np.savetxt(sampfile, samparray, delimiter=' ', newline='\n', header=header, comments='', fmt=format)
		except:
			exception = Exception('Error while passing the values of \"%s\" to the current.samp file.' % each)
			raise exception
		pass


	def SetModel(self, inputname, extrafiles=[], directory=os.getcwd()):
		"""
		Set the input script file to be used and extra files that should be copied together.
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

		# Verify if the input script file exists
		if not os.path.isfile(str(directory+'\\'+inputname)):
			exception = Exception('Current input script file (\"%s\") does not exists in (\"%s\") to be copied.' % (inputname, directory))
			raise exception

		# Copy the script file to workdirectory with the name 'current.inp'
		errcopy = os.system('copy /Y %s %s' % (str(directory+'\\'+inputname), str(self.ANSYSprops['run_location']+'\\current.inp')))
		if errcopy is not 0:
			exception = Exception('It was not possible to copy the input script file. (\"%s\")' % str(directory+'\\'+inputname))
			raise exception

		# Copy the extra files
		for each in extrafiles:
			errcopy = os.system('copy /Y %s %s' % (str(directory+'\\'+each), str(self.ANSYSprops['run_location']+'\\'+each)))
			if errcopy is not 0:
				exception = Exception('It was not possible to copy an extra file (\"%s\").' % each)
				raise exception

		# Clear last Model properties and set new
		self.Model = {}
		self.Model['inputname'] = inputname
		self.Model['extrafiles'] = extrafiles
		self.Model['directory'] = directory

		self._PrintR('Input script file and extra files copied to working directory.')
		self._PrintR('   Main APDL script: \"%s\"' % self.Model['inputname'])
		self._PrintR('   Extra model files: \"%s\"' % self.Model['extrafiles'])
		self._PrintR('   Input directory: \"%s\"' % self.Model['directory'])


	def _ClearForRun(self):
		"""
		Internal function to clear the entire directory or just the lock file before running
		"""
		# Verify the clear condition
		if self.ANSYSprops['cleardir']:
			# Try to clear the working directory
			self._PrintR('Cleaning the files from ANSYS working directory (\"%s\")' % self.ANSYSprops['run_location'])
			delhand = os.system('del /q %s\\*' % self.ANSYSprops['run_location'])
			if delhand is not 0:
				exception = Exception('Unable to clear the ANSYS working directory.')
				raise exception

			# Put the model back in the working directory
			self.SetModel(self.Model['inputname'], self.Model['extrafiles'], self.Model['directory'])

		else:
			# Verify the override condition
			lockfile = '%s\\%s.lock' % (self.ANSYSprops['run_location'], self.ANSYSprops['jobname'])
			if self.ANSYSprops['override'] and os.path.isfile(lockfile):
				self._PrintR('Deleting lock file.')
				delhand = os.system('del /q %s' % lockfile)
				if delhand is not 0:
					exception = Exception('Unable to delete lock file (\"%s\").' % lockfile)
					raise exception

			# Before execution ALWAYS erase the $jobname$.err file
			self._PrintR('Deleting old error log file.')
			errfile = '%s\\%s.err' % (self.ANSYSprops['run_location'], self.ANSYSprops['jobname'])
			delhand = os.system('del /q %s' % errfile)
		pass


	def _Run(self):
		# ANSYS Run Parameters https://www.sharcnet.ca/Software/Ansys/16.2.3/en-us/help/ans_ope/Hlp_G_OPE3_1.html
		# Ansys PDS commands = pdsrun.inp and out=pdsout.out
		#
		# DO NOT REMOVE THE SPACES BEFORE -np !
		#
		cmd = 'start "ANSYS" /d "%s" /min /wait /b "%s" "   -np %d -j %s -b -i pdsrun.inp -o pdsout.out %s" ' % (self.ANSYSprops['run_location'],
				self.ANSYSprops['exec_loc'], self.ANSYSprops['nproc'], self.ANSYSprops['jobname'], self.ANSYSprops['add_flags'])

		self._PrintR('Running ANSYS.')

		# Get time before run ANSYS
		timei = time.time()

		ansyshand = os.system(cmd)
		if ansyshand is not 0:
			exception = Exception('ANSYS exited with error id=%d. Please verify the output file (\"%s\").\n\n' %(ansyshand, self.ANSYSprops['run_location']+'\\pdsout.out'))
			raise exception

		# verify if it has an error on $jobname$.err
		errfile = '%s\\%s.err' % (self.ANSYSprops['run_location'], self.ANSYSprops['jobname'])
		f = open(errfile).read()
		#if '*** ERROR ***' in f:
		#	exception = Exception('ANSYS exited with an error. Please verify the output file (%s).\n\n' % (self.ANSYSprops['run_location']+'\\pdsout.out'))
		#	raise exception

		# Time after ANSYS
		timef = time.time()

		self._PrintR('Solution is done. It took %f minutes.' % ((timef-timei)/60))


	def Run(self):
		"""
		Execute the analysis on ANSYS.

		"""
		# Verify if model was set
		if self.Model == {}:
			exception = Exception('Before running ANSYS you must define the model that will be analysed with SetModel().')
			raise exception

		# Verify if there are output variables
		if self.varOutNames == []:
			exception = Exception('Before running ANSYS you must define output variables.')
			raise exception

		# Performn the conditional clear
		self._ClearForRun()

		# Write the pdsrun.inp file
		self._writePDSfile()

		# Write the sample file
		self._writeSAMPfile()

		# Run
		self._Run()


	@property
	def Length(self):
		"""
		Return the number of analysis to be executed and the length of INPUT arrays.

		"""
		return self.length


	def SetLength(self, length):
		"""
		Define the number of analysis to be executed and the length of INPUT arrays.
		Must be set before the variables.

		Parameters
		----------
		length: int, obligatory
			Number of analysis.

		"""

		# control variable
		valuesset = False

		# if at least one InValue is set can't change de length
		for each in self.varInNames:
			if self.varInValues[each] is not None:
				valuesset = True
				exception = Exception('At least one value were already set to INPUT variables. \n'+
							'To change the length you have to clear the variables values with ClearValues().')
				raise exception
				break

		# if there is no value setted and length > 0 can change
		if valuesset is False:
			if length <= 0:
				exception = Exception('The analysis length must be greater than 0.')
				raise exception
			else:
				self.length = length
				self._PrintR('Analysis lenght set to %d.' % self.length)


	def CreateVarIn(self, name):
		"""
		Create an INPUT variable.

		Parameters
		----------
		name : str, obligatory
			Variable name.

			Do not use spaces in the name!

			The same name cannot be used with INPUT and OUTPUT variables.

			It's not case sensitivy, in fact it will be saved in uppercase
			because of ANSYS.

		"""
		if name in ['', None, ' ']:
			exception = Exception('You must define an unique name to each variable.')
			raise exception
		else:
			name = name.upper()
			if name in self.varInNames or name in self.varOutNames:
				exception = Exception('The same name cannot be used with INPUT and OUTPUT variables.')
				raise exception
			else:
				self.varInNames.append(name)
				self.varInValues[name] = None
				self._PrintR('ANSYS input variable \"%s\" created.' % name)


	def CreateVarOut(self, name):
		"""
		Create an OUTPUT variable.

		Parameters
		----------
		name : str, obligatory
			Variable name.

			Do not use spaces in the name!

			The same name cannot be used with INPUT and OUTPUT variables.

			It's not case sensitivy, in fact it will be saved in uppercase
			because of ANSYS.
		"""
		if name in ['', None, ' ']:
			exception = Exception('You must define an unique name to each variable.')
			raise exception
		else:
			name = name.upper()
			if name in self.varInNames or name in self.varOutNames:
				exception = Exception('The same name cannot be used with INPUT and OUTPUT variables.')
				raise exception
			else:
				self.varOutNames.append(name)
				self.varOutValues[name] = None
				self._PrintR('Variable \"%s\" declared as ANSYS output variable.' % name)


	def SetVarInValues(self, name, values):
		"""
		Set the values of an INPUT variable

		Parameters
		----------
		name : str, obligatory
			Input variable name that will receive the values.

		values : 1D np.array of floats, obligatory
			A 1D numpy.array() with the length of this class and the values to be analysed.
			If the array is not 1D just the first column (0) will be used.
		"""

		# Verify the set length
		if self.length <= 0:
			 exception = Exception('Before set the values you must define the length of the analysis with SetLength().')
			 raise exception

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

		# Verify if it's a list, in this case transf to an np.array
		if type(values) is list:
			values = np.array(values)

		# Verify if it has more than one column
		try:
			values.shape[1]
		except:
			pass
		else:
			self._PrintR('The values that are trying to be set to \"%s\" has more than one column, just the first (0) will be used.' % (name))
			values = np.copy(values[:,0])

		# Verify if length is GT or LT the expected value
		if values.shape[0] < self.length:
			exception = Exception('The length of values that are trying to be set to \"%s\" is less than the defined length (%d).' % (name, self.length))
			raise exception
		elif values.shape[0] > self.length:
			self._PrintR('The length of values that are trying to be set to \"%s\" is greater than the defined length (%d).\
				  \nJust the first %dth values will be used.' % (name, self.length, self.length))

			values = values[:self.length]

		# Apply
		try:
			self.varInValues[name] = values
		except:
			exception = Exception('Error setting the values of \"%s\"' % name)
			raise exception
		else:
			self._PrintR('Values of \"%s\" were set.' % name)


	def GetVarOutValues(self):
		"""
		Return the values of the Results file from ANSYS for all the OUTPUT variables

		"""

		# Verify the existence of results file "$jobname$_current.pdrs"
		resultsfile = '%s\\%s_current.pdrs' % (self.ANSYSprops['run_location'], self.ANSYSprops['jobname'])
		if os.path.isfile(resultsfile):

			# Import the results file
			self._PrintR('Importing results from PDS results file (\"%s\").' % resultsfile)
			resultsALL = np.genfromtxt(resultsfile, names=True, skip_header=1)

			# Verify the existence of ERRORS
			errorcount = resultsALL['ERR'].sum()
			if errorcount > 0:
				self._PrintR('\n\n\nATTENTION:\nThe 4th column in the results file called ERR warns about simulations that results should not be used.\n'+
					  'Current file has %d errors, you can acces this in the column ERR returned with the results.\n\n\n' % errorcount)
				_ = input('The execution is paused, press ENTER to continue.\n')

			# Get all the columns in varOutNames and varInNames the specified length
			results = {}
			#results['ERR'] = np.copy(resultsALL['ERR'][:self.length])
			results['ERR'] = np.copy(resultsALL['ERR'])

			#for each in self.varInNames:
			#	results[each] = np.copy(resultsALL[each])

			for each in self.varOutNames:
				results[each] = np.copy(resultsALL[each])

			for each in results:
				if results[each].size > self.length:
					results[each] = results[each][:self.length]

			return results
		else:
			# There is no results
			exception = Exception('There is no results file in current ANSYS working directory. \n'+
								  'Please verify if the analysis was run.')
			raise exception


	def ClearAll(self):
		"""
		Clear all the properties (not from ANSYS object)
		"""
		# Clear all
		self.Model = {}
		self.varInNames = []
		self.varOutNames = []
		self.varInValues = {}
		self.varOutValues = {}
		self.length = 0
		#self.PrintR = False

		self._PrintR('All the properties were cleared (not from ANSYS object).')

	def ClearValues(self):
		"""
		Clear the values of all variables.
		"""
		# Clear each input variable value
		for each in self.varInNames:
			self.varInValues[each] = None

		# Clear each output variable value
		for each in self.varOutNames:
			self.varOutValues[each] = None

		self._PrintR('The values of all variables were cleared. Now you can change the length parameter.')
