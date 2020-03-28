import re, sys, os
from Bio.PDB import *
def get_interactions_from_dir(input_dir):
    """Searches for files containing interactions between two PDB chains.
		Arguments:
			input_dir -- path to the directory containing the pdb interaction files.

		Returns:
			List with the input files that are identified as interaction files.

		Exceptions:
			- When any of the files names' does not correspond to the specified interaction filename, it prints an error message.
"""

    if os.path.isdir(input_dir):
        dir_files = os.listdir(input_dir)
        input_files = [element for element in dir_files if re.match('(....([A-Z 0-9 a-z]+))\.pdb', element)]
    else:
        print("%s is not a valid argument, must provide a directory containing pdb interaction files. The name of the files must be in the following format: 4 characters corresponding to the name of the complex + 2 characters corresponding the chains IDs + '.pdb'" %options.indir, file= sys.stderr)
        exit()

    return(input_files)


def check_interactions(input_files, input_dir):
    """Associates input files with their corresponding Biopython's Structure class objects and searches for interactions in which a chain participates.

		Arguments:
			input_files -- List with the input files that are identified as interaction files.

		Returns:
			Dictionary with chain ID as keys and a list of interactions in which that chain participates as value.
"""

    interaction = {}
    parser = PDBParser(PERMISSIVE=1,QUIET=1)
    for file in input_files:
        m = re.match('(....([A-Z 0-9 a-z]+))\.pdb', file)
        path = input_dir + file
        interaction[m.group(2)] = parser.get_structure(m.group(1), path)

    my_complex = {}
    for element in interaction:
        my_complex[element[0:1]] = my_complex.setdefault(element[0:1], [interaction[element]])
        if interaction[element] not in my_complex[element[0:1]]:
            my_complex[element[0:1]].append(interaction[element])

    return(my_complex)

def check_outdir(outdir, indir):
	"""Checks if output directory given by user is really a directory and sets an output name depending on input directory.

		Arguments:
			outdir -- path of the directory where to store results.
			indir -- path to the directory containing the input files.
		Returns:
			Path of the main script output_file.
		Exceptions:
			- When the given path does not correspond to a directory, it prints an error message.
"""
	output_file = indir.split('/')
	if os.path.isdir(outdir):
	    if '/' in outdir:
	        output_path = outdir + output_file[len(output_file)-2] + "_model.pdb"
	    else:
	        output_path = outdir + "/" +  output_file[len(output_file)-2] + "_model.pdb"
	else:
	    print("%s is not a valid argument, must provide a directory to store files" %options.outdir, file= sys.stderr)
	    exit()

	return(output_path)
