from Bio import Align, SeqIO, pairwise2
from Bio.Seq import Seq
import re, sys
from Bio.Alphabet import IUPAC, ProteinAlphabet

def SequenceAlignment(fasta_file):
	"""Perfoms pairwise sequence alignments to compare % of identity.

		Arguments:
			fasta_file -- path to the fasta file

		Returns:
			Tuple containing: (1) a dictionary with the chain ID as keys and the equivalent IDs as values, (2) a list of the chain IDs that correspond to nucleic acids.

		Exceptions:
			- When the format of the file does not correspond a fasta file, it prints an error message.
			- AttributeError: When some header does not match the regex, it prints an error message.
	"""

	first = 0
	chains_equivalence = {}
	chains_included = []
	dna_letters = 'GATC'
	rna_letters = 'GAUC'
	no_proteins = []
	with open(fasta_file, "rU") as handle:
		record = list(SeqIO.parse(handle, "fasta"))
		if not record:
			print("ERROR! The file %s is not a fasta file. Please provide a file in fasta format." %(fasta_file), file = sys.stderr)
			exit()
		record2 = record.copy()
	for sequence in record:
		seq1 = Seq(str(sequence.seq))
		protein_count = 0
		name1 = re.match('....:(.)', sequence.id)
		try:
			 name_chain1 = name1.group(1)
		except AttributeError:
			print("ERROR! The sequence header %s in the file %s does not have the correct format.\nCorrect header format: '>' + 4 characters corresponding to the name of the complex + ':' + 1 character corresponding to the ID of the chain.\nExample: >2C7A:A" %(sequence.id,fasta_file), file = sys.stderr)
			exit()

		for letter in str(sequence.seq):
			if (letter not in dna_letters) and (letter not in rna_letters):
				protein_count += 1

		if	protein_count == 0:
			no_proteins.append(name_chain1)

		if name_chain1 not in chains_included:
			chains_equivalence[name_chain1] =  [name_chain1]
			chains_included.append(name_chain1)
			for seq in record2:
				if sequence.id != seq.id:
					seq1 = Seq(str(sequence.seq))
					seq2 = Seq(str(seq.seq))
					alignments = pairwise2.align.globalxx(seq1,seq2)
					for alignment in alignments:
						score = alignment[2]/alignment[4]
						if (score > 0.95):
							name2 = re.match('....:(.)', seq.id)
							try:
								name_chain2 = name2.group(1)
							except AttributeError:
								print("ERROR! The sequence header %s in the file %s does not have the correct format.\nCorrect header format: '>' + 4 characters corresponding to the name of the complex + ':' + 1 character corresponding to the ID of the chain.\nExample: >2C7A:A" %(seq.id,fasta_file), file = sys.stderr)
								exit()

							if name_chain2 not in chains_equivalence[name_chain1]:
								chains_equivalence[name_chain1].append(name_chain2)
								chains_included.append(name_chain2)
		record2.remove(sequence)

	return(chains_equivalence, no_proteins)
