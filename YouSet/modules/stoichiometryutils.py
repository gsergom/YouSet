import re, sys

def ReadStoichiometry(dict_equivalence):
	"""Allows the user to enter the desired stoichiometry of the complex. Reads the equivalences of chain IDs and asks the user to give the absolute frequency for each type of chain.

		Arguments:
			dict_equivalence -- dictionary of chain equivalences

		Returns:
			stoich -- dictionary with chain IDs as keys and absolute frequencies as values

		Exceptions:
			Only positive integers are accepted as values. Otherwise, it prints an error message.
	"""

	stoich = {}
	print("Provide the number of times that each chain appears in the complex and press ENTER.")
	for element in dict_equivalence:
		print_string = "Chain " + element + " (these chains are considered the same: " + ", ".join(dict_equivalence[element]) + "): "
		number = input(print_string)
		while True:
			try:
				number = int(number)
			except ValueError:
				print("This is not a valid value. Please enter an integer.")
				number = input(print_string)
			else:
				if number < 1:
					print("This is not a valid value. Please enter a positive integer.")
					number = input(print_string)
					continue
				else:
					break
		stoich[element] = stoich.setdefault(element,int(number))
	print("Stoichiometry has been readed correctly.")

	return(stoich)

def check_stoich(result_build_complex, seq_alignments):
	"""Checks the stoichiometry of the final complex model and prints it to standard error

		Arguments:
		result_build_complex -- tuple containing: (1) list of the chains that are in the final complex, (2) dictionary of chain equivalences
		seq_alignments -- tuple containing: (1) a dictionary with the chain ID as keys and the equivalent IDs as values, (2) a list of the chain IDs that correspond to nucleic acids.
"""
	stoichiometry = {}
	for element in result_build_complex[0]:
		for chain in result_build_complex[2]:
			if (element in seq_alignments[0][chain]):
				stoichiometry[chain] = stoichiometry.setdefault(chain, 0) + 1

	stoichiometry_string = ''
	stoichiometry_str_nucleic_acids = ''
	total_nucleic_acids = 0
	for element in stoichiometry:
		if element not in seq_alignments[1]:
			stoichiometry_string += element+str(stoichiometry[element])
		else:
			stoichiometry_str_nucleic_acids += element+str(stoichiometry[element])
			total_nucleic_acids += stoichiometry[element]
	print("The stoichiometry of the final model is %s + %s nucleic acid(s): %s." %(stoichiometry_string, total_nucleic_acids, stoichiometry_str_nucleic_acids), file = sys.stderr)
