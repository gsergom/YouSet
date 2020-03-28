#!/usr/bin/env python3

import sys
from Bio.PDB import *

aa=['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
# create parser
parser = PDBParser(PERMISSIVE=1,QUIET=1)
inpu = sys.argv[1]
# read structure from file
structure = parser.get_structure(inpu,inpu+'.pdb')
model = structure[0]
chains = {} 
for chain in model:
	chains[chain.get_id()] = chain
#print(chains)
distances=[]
for element in chains:
	chain1 = chains[element]
	for element2 in chains:
		if element2 != element:
			chain2 = chains[element2]
			for residue1 in chain1:
				for residue2 in chain2:
					if residue1.get_resname() in aa and residue2.get_resname() in aa:
						while True:
							try:
								distance = residue1['CA'] - residue2['CA']
								if distance < 8:
									merged_chains = chain1.get_id()
									merged_chains += chain2.get_id()
									if merged_chains not in distances and ''.join(reversed(merged_chains)) not in distances:
										distances.append(merged_chains)
								break
							except:
								break
					elif residue1.get_resname() in aa or residue2.get_resname() in aa: 
						if residue1.get_resname() in aa: 
							for atom in residue2.get_atoms():
								while True:
									try:
										distance = residue1['CA'] - atom
										if distance < 8:
											merged_chains = chain1.get_id()
											merged_chains += chain2.get_id()
											if merged_chains not in distances and ''.join(reversed(merged_chains)) not in distances:
												distances.append(merged_chains)
										break
									except:
										break
						if residue2.get_resname() in aa: 
							for atom in residue1.get_atoms():
								while True:
									try:
										distance = atom - residue2['CA']
										if distance < 8:
											merged_chains = chain1.get_id()
											merged_chains += chain2.get_id()
											if merged_chains not in distances and ''.join(reversed(merged_chains)) not in distances:
												distances.append(merged_chains)
										break
									except:
										break
					else:
						for atom1 in residue1.get_atoms():
							for atom2 in residue2.get_atoms():
								while True:
									try:
										distance = atom1 - atom2
										if distance < 8:
											merged_chains = chain1.get_id()
											merged_chains += chain2.get_id()
											if merged_chains not in distances and ''.join(reversed(merged_chains)) not in distances:
												distances.append(merged_chains)
										break
									except:
										break

for element in distances:
	splitted_chains = [element[0:1], element[1:2]]

	class ChainSelect(Select):
	    def accept_chain(self, chain):
	        if chain.get_id()== splitted_chains[0] or chain.get_id()== splitted_chains[1]:
	            return 1
	        else:
	            return 0

	io = PDBIO()
	io.set_structure(structure)
	output = inpu +element+'.pdb'
	#print(output)
	io.save(output, ChainSelect())
