#!/usr/bin/env python3

import argparse, sys
from modules import *

#####SETTING ARGUMENTS FOR THE SCRIPT#####

arguparser = argparse.ArgumentParser(description = "Superimposition macrocomplex builder by Paula Lopez, Gerard Serrano and Laura Vila for the SBI subject")

required_args = arguparser.add_argument_group('mandatory arguments')

required_args.add_argument('-d', '--pdbdir',
                    dest = 'indir',
                    action = 'store',
                    required = True,
                    default = None,
                    help = 'Path to directory with all the base pdb files to use')

required_args.add_argument('-f', '--fasta',
                    dest = 'infasta',
                    action = 'store',
                    default = None,
                    required = True,
                    help = 'Path to fasta file with the sequences of the complex')

required_args.add_argument('-o', '--output',
                    dest = 'outdir',
                    action = 'store',
                    required = True,
                    default = None,
                    help = 'Path to output directory where to store the results')

arguparser.add_argument('-v', '--verbose',
                    dest = 'verbose',
                    action = 'store_true',
						  default = False,
                    help = 'Print log to stderr')

arguparser.add_argument('-s', '--stoichiometry',
                    dest = 'stoich',
                    action = 'store_true',
						  default = False,
                    help = 'Optional stoichiometry of the macrocomplex')

options = arguparser.parse_args()


##### GETTING INPUT FILES AND INTERACTIONS FROM DIR #####
if options.verbose:
    print('Checking files in input directory...', file = sys.stderr)

input_files = get_interactions_from_dir(options.indir)

if options.verbose:
    print("Done! %i files found." %len(input_files), file=sys.stderr)

if options.verbose:
    print('Associating chain interactions...', file = sys.stderr)

my_complex = check_interactions(input_files, options.indir)

if options.verbose:
    print('Done!', file = sys.stderr)

#### CHECKING OUTPUT DIRECTORY
output_path = check_outdir(options.outdir, options.indir)

##### CHECKING SEQUENCE EQUIVALENCES #####
if options.verbose:
    print('Checking for same chains with different names...', file = sys.stderr)

seq_alignments = SequenceAlignment(options.infasta)

if options.verbose:
    print('Done!', file = sys.stderr)

#### READING STOICHIOMETRY FROM THE USER ####
if options.stoich:
	print('%s different chains were found.' %(len(seq_alignments[0])), file= sys.stderr)
	output_stoich = ReadStoichiometry(seq_alignments[0])
else:
	output_stoich = {}
	for element in seq_alignments[0]:
		output_stoich[element] = output_stoich.setdefault(element, len(seq_alignments[0][element]))

##### RECONSTRUCT THE MACROCOMPLEX #####
if options.verbose:
	print("Reconstructing the macrocomplex...", file=sys.stderr)

result_build_complex = BuildComplex(my_complex, output_path, seq_alignments[0], output_stoich)

if options.verbose:
	print("Done!", file=sys.stderr)


#### Print if some interactions cannot be added to the complex
if len(result_build_complex[1]) != 0:
	for element in result_build_complex[1]:
		for interaction in result_build_complex[1][element]:
			print("The interaction %s could not be added to the final complex." %(interaction.get_id()), file=sys.stderr)

##### Stoichiometry_check #####
if options.verbose:
    check_stoich(result_build_complex, seq_alignments)
    print("Final model is %s" %output_path, file = sys.stderr)

exit(0)
