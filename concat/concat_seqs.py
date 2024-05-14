#!/usr/bin/env python3 
import os, sys
import re
import argparse
from glob import glob 
from collections import Counter

FASTA_REGEX = re.compile(r"\.fas?t?a?$", re.IGNORECASE)

# modify this if 
SEGMENTS = "PB2 PB1 PA HA NP NA M NS".split()
SEGMENTS_REGEX = [re.compile(r"[-_\|]" + s) for s in SEGMENTS]

def init_parser():
	parser = argparse.ArgumentParser(description="Concatenate (or extract specific segment from) influenza FASTA files.")
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("-s", "--segment", default='NONE', type=str, help="Named argument to extract a specific segment from all FASTA files. This works as a regular expression search.")
	group.add_argument("-c", "--concat", action='store_true', help="Flag to concatenate all segments together")

	parser.add_argument("-i", "--input",  required=True, help="Either single FASTA file, or folder of files, containing 8 segments per file.")
	parser.add_argument('-o', "--output", required=True, help="If using file input, output FASTA file name. If using folder input, output folder name.")
	parser.add_argument('-n','--header', default=None, help='Concatenated sequence header name. Only relevant in single file concat mode. Otherwise, default is the first element when splitting filename by both "_" and "."')

	return parser

def parse_fasta(filepath):
    """
    Parses a FASTA file and returns a list of tuples, where each tuple contains the header and sequence.

    Args:
        filepath (str): The path to the FASTA file.

    Returns:
        list: A list of tuples where each tuple contains the header and sequence.
    """
    # Initialize empty list to store the sequences.
    seqlist = []
    
    # Initialize empty variables to store the current header and sequence.
    header = ''
    seq = ''
    
    # Open the FASTA file and read each line.
    with open(filepath, 'r') as handle:
        for line in handle.readlines():
            # If the line starts with '>', it is a new header.
            if line[0] == '>':
                # If there is a previous header and sequence, add them to the list.
                if header != '':
                    seqlist.append((header, seq))
                
                # Update the header and clear the sequence.
                header = line.strip().lstrip('>')
                seq = ''
            else:
                # Add the sequence to the current sequence.
                seq += line.strip()
        
        # Add the last header and sequence to the list.
        seqlist.append((header, seq))
    
    # Return the list of sequences.
    return seqlist


def write_fasta(seqs, outpath):
	"""
    Writes sequences to a FASTA file.
    Args:
        seqs (list): List of tuples where each tuple contains a header and a sequence.
        outpath (str): Output file path.
    Returns:
        None
    """
	with open(outpath, 'w') as outfile:
		for header, seq in seqs:
			outfile.write(f">{header}\n{seq}\n")
	

def concatenate_sequences(filepath, header):
	"""
    Concatenates sequences from a FASTA file.

    Args:
        filepath (str): Path to the input FASTA file.
        header (str): Header for the concatenated sequence.

    Returns:
        tuple: A tuple containing the header and concatenated sequence.
              If the FASTA file contains less than 8 sequences, None is returned.
    """
	# Parse the FASTA file
	sequences = parse_fasta(filepath)

	# Check if the FASTA file contains 8 sequences
	if len(sequences) != 8:
		print(f"ERROR: Skipping {filepath}. There are {len(sequences)} sequences present instead of 8.")
		return None

	# Sort the sequences according to the SEGMENTS order
	sequences_sorted = []
	for segment, expr in zip(SEGMENTS, SEGMENTS_REGEX):
		# Search for sequences that match the current segment
		result = [seq for seq in sequences if expr.search(seq[0])]

		# Catch specific case where there are off target regex matches
		counts = [len(expr.search(seq[0]).groups())  for seq in sequences if expr.search(seq[0])]
		counter = Counter(counts)
		most_common = counter.most_common(1)[0]

		# Find best match amidst off target regex matches
		if len(result) == 8 and len(counter) > 1 and most_common[1] == 1:
			N_best_match = [n for n, x in enumerate(counts) if x == most_common[0]][0]
			sequences_sorted.append(result[N_best_match])
		
		# Normal case, no off-target matches, one result returned
		elif len(result) == 1:
			sequences_sorted.append(result[0])
		
		# Incorrect number of results, return None
		else:
			print(f"ERROR: Could not isolate sequence for segment {segment}. {len(result)} results found. Exiting.")
			return None

	
	# Concatenate the sequences
	concatenated_seq = ''
	for _ , seq in sequences_sorted:
		concatenated_seq += seq

	# Return the concatenated sequence with the provided header
	return (header, concatenated_seq)

def extract_segment(filepath, target_segment):
	"""
	Extracts a single sequence from a FASTA file that matches the target segment.

	Args:
		filepath (str): Path to the input FASTA file.
		target_segment (str): Regular expression to match the target segment.

	Returns:
		tuple: A tuple containing the header and sequence that matches the target segment.
		      If no sequence is found, None is returned.
	"""
	# Parse the FASTA file
	sequences = parse_fasta(filepath)

	# Search for sequences that match the target segment
	search_segments = [x for x in sequences if re.search(target_segment, x[0])] 

	# Check if exactly one sequence is found
	if len(search_segments) != 1: 
		# If no sequence is found, print an error message
		print("ERROR: Could not locate a singular sequence that matches the segment of interest. You may need to specify a more advanced / specific regex using the --segment argument.")
		print("Search results: ", search_segments)
		return None 

	# Return the single sequence that matches the target segment
	return search_segments[0]

def process_fasta(infile, outfile, segment, concat, header):
	# Set the output file name
	header = header if header else re.split(r"_|\.", os.path.basename(infile))[0]

	# Concatenate the sequences if the --concat flag is provided
	# Otherwise, extract the sequence that matches the target segment
	if concat:
		seq = concatenate_sequences(infile, header)
	else:
		seq = extract_segment(infile, segment)
	
	if seq:
		# Write the sequence to the output file
		write_fasta([seq], outfile)


def main(args):
	"""
	Main function to handle the execution of the script.

	Args:
		args (argparse.Namespace): The parsed command line arguments.
	"""

	# If the input is a file, process it
	if os.path.isfile(args.input):
		print("Processing: " + os.path.basename(args.input))
		# Call the process_fasta function with the input file, output file, segment, concat, and header arguments
		process_fasta(args.input, args.output, args.segment, args.concat, args.header)
	
	# If the input is a directory, process all FASTA files in the directory
	elif os.path.isdir(args.input):
		
		# Create the output directory if it does not exist
		os.makedirs(args.output, exist_ok=True)

		# List all files in the input directory that match the FASTA regex
		file_list = glob(os.path.join(args.input, '*'))
		file_list = [x for x in file_list if FASTA_REGEX.search(x)]

		# Iterate over each file
		for f in file_list:
			# Print the name of the file being processed
			print("Processing: " + os.path.basename(f))
			# Set the output file path
			outfile = os.path.join(args.output, os.path.basename(f))
			# Call the process_fasta function with the input file, output file, segment, concat, and header arguments
			process_fasta(f, outfile, args.segment, args.concat, args.header)
	# If the input file/folder does not exist, print an error message
	else:
		print("ERROR: Input FASTA file/folder does not exist", file=sys.stderr)


if __name__ == "__main__":
	args = init_parser().parse_args()
	main(args)