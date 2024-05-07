import os, sys, pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import argparse


def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--fasta',  required=True, help='FASTA file containing all query sequences in correct order. Assuming GenBank format starting with accession number followed by a space.')	
	parser.add_argument('-g', '--gffpath', required=True, help='Input path containing GFF3 files. Created using "fetch_gff3.py".')	
	parser.add_argument('-o', '--outgff', required=True, help='Output path for single combined GFF3 file')	
	parser.add_argument('-O', '--outfasta', required=True, help='Output path for single combined FASTA file')	
	parser.add_argument('-n', '--name', default='NC_000000.0', help='The name for the concatenated reference sequence.')	
	parser.add_argument('-F', '--force', action='store_true', help='Overwrite outputs.')	
	
	return parser

def parse_gff_header(filepath, concat_name):
	header = []

	with open(filepath, 'r') as infile:
		for line in infile.readlines():
			if line.startswith('#'):
				header.append(line)

	# modify the header line containing the accession number & length information 
	header = [f"##sequence-region {concat_name} 1 {{}}\n" if x.startswith('##sequence-region') else x for x in header]

	return "".join(header)

def parse_gff(filepath):
	cols = 'seqid source type start end score strand phase attributes'.split()
	return pd.read_csv(filepath, sep='\t', comment='#', names=cols)

def concat_fasta_seq(filepath, concat_name):

	concat_seq = SeqRecord('', id=concat_name, name=concat_name, description=concat_name)
	
	for s in SeqIO.parse(filepath, 'fasta'):
		concat_seq.seq += s.seq

	return concat_seq


def build_combined_gff3(seq_path, gff_path, concat_name):

	gff_list = []
	position = 0 
	header = ''

	for n , s in enumerate(SeqIO.parse(seq_path, 'fasta')):

		if n == 0:
			header = parse_gff_header(os.path.join(gff_path, s.id + '.gff3'), concat_name)
		
		df = parse_gff(os.path.join(gff_path, s.id + '.gff3'))

		df['start'] += position
		df['end'] += position

		gff_list.append(df)
		position += len(s)

	combined = pd.concat(gff_list).reset_index(drop=True)
	combined['seqid'] = concat_name
	return combined, header


def main(args):
	# concatenate FASTA sequences together
	concat_seq = concat_fasta_seq(args.fasta, args.name)

	# Build concatenated GFF3 file
	gff_combined, gff_header = build_combined_gff3(args.fasta, args.gffpath, args.name)
	
	# output concatenated GFF3 file 
	if os.path.isfile(args.outgff) or os.path.isfile(args.outfasta):
		if args.force:
			os.remove(args.outgff)
			os.remove(args.outfasta)
		else:
			print("ERROR: Outputs already exist. Avoiding overwrite. If overwrite desired, specify the -F force flag.")
			sys.exit(1)

	# write the new GFF3 file
	with open(args.outgff, 'a') as outfile:
		# write the GFF header and load in the new concatenated sequence length
		outfile.write(gff_header.format(str(gff_combined['end'].max())))

		# write the GFF dataframe contents after the header
		gff_combined.to_csv(outfile, sep='\t', index=False, header=False)

	# write the concatenated FASTA file 
	SeqIO.write(concat_seq, args.outfasta, 'fasta')

if __name__ == "__main__":
	args = init_parser().parse_args()
	main(args)