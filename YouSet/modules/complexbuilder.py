from Bio.PDB import *

def CheckingNewChains(chains_movil, chains_fixed, equivalences, stoich):
	"""Takes the list of chains that are already in the complex and two chains that are evaluating for being added to the complex. Return a

	Arguments:
		chains_movil -- list of chains objects that are evaluating for being added to the complex
		chains_fixed -- list of chains objects that are already in the complex
		equivalences -- dictionary with chain IDs as keys and the equivalent IDs as values
		stoich -- dictionary with chain IDs as keys and absolute frequencies as values

	Returns:
		Dictionary with the IDs of chains_movil as keys, and as values strings indicating if they are 'Not yet in the complex', if they are 'Accepted' of if they are 'Not accepted'

"""
	possible_chains = {}
	for x in chains_movil:
		possible_chains[x.get_id()] = possible_chains.setdefault(x.get_id(), "Not yet in the complex")
		for y in chains_fixed:
			for element in equivalences:
				if y.get_id() in equivalences[element]:
					equivalents_chains_fixed = equivalences[element]
					break

			if (x.get_id() in equivalents_chains_fixed):
				if stoich[equivalents_chains_fixed[0]] > 0:
					possible_chains[x.get_id()] = "Accepted"
					break
				else:
					possible_chains[x.get_id()] = "Not accepted"
					break
	return possible_chains

def PossibleChains(equivalences, fixed_names, chain_name):
	"""Searches for all the identical chains to chain_name that are already included in the complex.

	Arguments:
		equivalences -- dictionary with chain IDs as keys and the equivalent IDs as values
		fixed_names -- list of chains IDs that are already included in the complex
		chain_name -- string representing the chain ID

	Returns:
		List containing the chain IDs

"""

	possible_chains = []
	for chain in equivalences:
		if chain_name in equivalences[chain]:
			for element in equivalences[chain]:
				if element in fixed_names:
					possible_chains.append(element)
	return possible_chains

def UsedIdentifiers(equivalences):
	"""Reads all the identifiers that have already been used by the user or have been newly added to the complex.

	Arguments:
		equivalences -- dictionary with chain IDs as keys and the equivalent IDs as values

	Returns:
		List with the identifiers
	"""

	str_chains = []
	int_chains = []
	for element in equivalences:
		for x in equivalences[element]:
			if x in ['1','2','3','4','5','6','7','8','9']:
				int_chains.append(x)
			else:
				str_chains.append(x)
	str_chains = sorted(str_chains)
	int_chains = sorted(int_chains)
	for element in int_chains:
		str_chains.append(element)
	return str_chains

def FindClashes(chains_fixed, re_chain_movil):
	"""Searches for clashes between pairs of chains.

	Arguments:
		chains_fixed -- list of residues
		re_chain_movil -- list of chains

	Returns:
		True if it finds one clash between the residues of chains_fixed and any residue belonging to the chains of re_chain_movil
		False if it does not find any clash
	"""

	aa=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
	for chain in chains_fixed:
		re_chain_fixed = [i for i in chain]
		for r1 in re_chain_fixed:
			for r2 in re_chain_movil:
				if r1.get_resname() in aa and r2.get_resname() in aa:
					distance = r1['CA'] - r2['CA']
					if distance < 2:
						return True
				else:
					for atom1 in r1.get_atoms():
						for atom2 in r2.get_atoms():
							distance = atom1 - atom2
							if distance < 1:
								return True
	return False


def BuildComplex(protein, out_path, equivalences, stoich, recursion = 0, base_structure = None):
	"""Build a complex by superimposition with a given set of interactions.

		Mandatory arguments:
			protein -- dictionary with chain ID as keys and a list of interactions in which that chain participates as value
			out_path -- path where to store the output pdb file
			equivalences -- dictionary with chain IDs as keys and the equivalent IDs as values
			stoich -- dictionary with chain IDs as keys and absolute frequencies as values

		Optional arguments:
			recursion -- integer indicating the recursivity level. By default, it is set to 0.
			base_structure -- structure object representing the final complex. By default, it is set to None.

		Returns:
			Tuple containing: (1) list of the chains that are in the final complex, (2) dictionary of chain equivalences

		Exceptions:
			- UnboundLocalError: when there is an error in the equivalences dictionary, it prints an error message
			- PDBExceptions.PDBConstructionException: when it tries to add a chain with an ID that is already in the complex, it changes the ID to a not used ID. If no more IDs are availabe, it passes.
	"""

	protein_sorted = sorted(protein.items())
	not_used_chains_dict = {}
	alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456789" #possible identifiers for the chains
	identifiers = UsedIdentifiers(equivalences)

	if (base_structure != None):
		first = 1
		fixed = base_structure
		chains_fixed = [i for i in fixed.get_chains()]
		chains_fixed_names = [i.get_id() for i in fixed.get_chains()]

	else:
		first = 0
		chains_fixed = []

	for chain_name in protein_sorted:
		for structure in chain_name[1]:
			if first == 0:
				fixed = structure.copy()
				first = 1
				chains_fixed = [i for i in structure.get_chains()]
				chains_fixed_names = [i.get_id() for i in fixed.get_chains()]
				for element in equivalences:
					for initial_chain in chains_fixed_names:
						if initial_chain in equivalences[element]:
							stoich[element] -= 1
			else:
				movil = structure.copy()
				chains_movil = [i for i in movil.get_chains()]
				chains_movil_names = [i.get_id() for i in movil.get_chains()]
				try:
					new_chains_check = CheckingNewChains(chains_movil, chains_fixed, equivalences, stoich)
				except UnboundLocalError:
					print("ERROR! Some sequence in the fasta file may not have a correspondence with any structure in the pdb files.")
					exit()
				count = 0
				repeated = 0

				for element in new_chains_check:
					if new_chains_check[element] == "Not accepted":
						count += 1
						common_chain = element
					elif new_chains_check[element] == "Not yet in the complex":
						new_chain = element
						repeated += 1
					elif new_chains_check[element] == "Accepted":
						new_chain = element

				if (count == 2): #we cannot add any chain due to stoichoimetry
					continue

				elif (repeated == 2): #there is no chain to use as "fixed" in superimposition
					not_used_chains_dict[chain_name[0]] = not_used_chains_dict.setdefault(chain_name[0], [])
					not_used_chains_dict[chain_name[0]].append(structure)
					continue

				elif count == 1 or count == 0: #we can add at least 1 chain
					if count == 0:
						for x in chains_movil:
							common_chain = x.get_id()
							break

					possible_common_chains = PossibleChains(equivalences, chains_fixed_names, common_chain)
					added = 0
					no_useful_structures = 0
					for possible_chain in possible_common_chains:
						if added == 0 and no_useful_structures < len(possible_common_chains):
							common_chain_fixed = fixed[0][possible_chain]
							common_chain_movil = movil[0][common_chain]
							new_chain_movil = movil[0][new_chain]

							#### START SUPERIMPOSITION ####
							sup = Superimposer()
							#Getting ist of atoms
							at_chain_fixed = [i for i in common_chain_fixed.get_atoms()]
							at_chain_movil = [i for i in common_chain_movil.get_atoms()]
							#Calculate minimum lenght and cutting if it is necessary
							minim_length = min(len(at_chain_fixed), len(at_chain_movil))
							at_chain_fixed = at_chain_fixed[:minim_length]
							at_chain_movil = at_chain_movil[:minim_length]
							#Transforming atomic coordenates
							sup.set_atoms(at_chain_fixed,at_chain_movil)
							for model in movil:
								for chain in model:
									for residue in chain:
										for atom in residue:
											atom.transform(sup.rotran[0],sup.rotran[1])

							#### CHECKING CLASHES ######
							re_chain_movil = [i for i in new_chain_movil]
							clashes = FindClashes(chains_fixed, re_chain_movil)
							if clashes == True:
								no_useful_structures += 1
							else:
								try:
									fixed[0].add(movil[0][new_chain])
								except PDBExceptions.PDBConstructionException:
									old_chain = new_chain
									for letter in alphabet:
										if letter not in identifiers:
											new_chain = letter
											break
									movil[0][old_chain].id = new_chain
									identifiers.append(new_chain)
									for element in equivalences:
										if old_chain in equivalences[element]:
											equivalences[element].append(new_chain)
									try:
										fixed[0].add(movil[0][new_chain])
									except PDBExceptions.PDBConstructionException:
										pass
								added = 1
								chains_fixed = [i for i in fixed.get_chains()]
								chains_fixed_names = [i.get_id() for i in fixed.get_chains()]
								for element in equivalences:
									if new_chain in equivalences[element]:
										stoich[element] -= 1
						elif (no_useful_structures == len(possible_common_chains)):
							not_used_chains_dict[chain_name[0]] = not_used_chains_dict.setdefault(chain_name[0], [])
							not_used_chains_dict[chain_name[0]].append(structure)

	chains_fixed = [i for i in fixed.get_chains()]
	chains_fixed_names = [i.get_id() for i in fixed.get_chains()]
	if (len(not_used_chains_dict) != 0) and recursion <10:
		recursion += 1
		return BuildComplex(not_used_chains_dict, out_path, equivalences, stoich, recursion, fixed)
	else:
		#Create and store PDB file
		io = PDBIO()
		io.set_structure(fixed)
		outfile = open(out_path, "wt")
		io.save(outfile)
		outfile.close()
		return(chains_fixed_names, not_used_chains_dict, equivalences)
