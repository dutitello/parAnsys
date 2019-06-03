# -*- coding: UTF-8 -*-
import pyansyspds
import os

class example():
	"""
	This example use ANSYS to get the ultimate force on a plane truss.

	The truss model will be on the running directory with the name "current.inp",
	basically when the material yields the structure displacement will grow without
	growing the reaction force on the node.
	Current script will test the truss for 2 Elastic Modulis (E) and 3 yield
	tensions (fy) and return the maximum force on the node that displacement is
	applied (ultimate force) on variable FxMAX.

	"""

	def __init__(self):
		print('For running current example we need the ANSYS executable location \
			  (with ansys170.exe, for example for version 17.0), and an directory \
			  for ANSYS work.')

		# Get necessary values
		ansysexe = input('Please insert ANSYS executable location: ')
		ansyswdir = input('Please insert the directory that ANSYS will work: ')


		# Create ANSYS object
		ansys = pyansyspds.AnsysPDS(exec_loc=ansysexe, run_location=ansyswdir)

		# Define the model
		f = open('current.inp', 'wr')
		f.write('')
		f.close()

		ansys.SetModel(inputname='current.inp', extrafiles=[], directory=os.getcwd())

		ansys.Info(True)

		# Declarando variaveis de entrada
		ansys.CreateVarIn('Fy')
		#ansys.CreateVarIn('E')

		# Declarando variaveis de saida
		ansys.CreateVarOut('FxMAX')
		#ansys.CreateVarOut('fymax')

		# Definindo tamanho das analises
		ansys.SetLength(5)

		# Tenta puxar algum resultado
		#res = ansys.GetVarOutValues()
		ansys.SetVarInValues('fy', np.array([10,20,25,30,35]))

		#ansys.SetVarInValues('E', [20])


		# Run ANSYS
		ansys.Run()

		# Apaga variaveis
		#ansys.ClearValues()

		# Tenta puxar algum resultado
		res = ansys.GetVarOutValues()
