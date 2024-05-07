import os
from Bio import Entrez
from Bio import SeqIO
import argparse

def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument( 'fasta',  help='FASTA file containing all query sequences in correct order. Assuming GenBank header format starting with accession number followed by a space.')	
	parser.add_argument( 'outpath',  help='Output path for GFF3 files. Script will make directory if it doesnt exist.')	
	
	return parser

def fetch_gff3_files(accession_list, output_dir):

	for accession in accession_list:
		# Fetch the GFF3 file
		print(accession)
		handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gff3", retmode="text")
		gff3_data = handle.read()

		# Write the data to a file in the specified output directory
		with open(os.path.join(output_dir, f"{accession}.gff3"), "w") as file:
			file.write(gff3_data)

		handle.close()

def main(args):
	# optionally load the user's NCBI API key 
	if os.environ['NCBI_EMAIL'] and os.environ['NCBI_API_KEY']:
		Entrez.email = os.environ['NCBI_EMAIL']
		Entrez.api_key = os.environ['NCBI_API_KEY']

	# Read accession numbers from a file or list them here
	accession_numbers = [x.id for x in SeqIO.parse(args.fasta, 'fasta')]

	# Make the output directory if it does not exist 
	os.makedirs(args.outpath, exist_ok=True)

	# Fetch and write GFF3 files
	fetch_gff3_files(accession_numbers, args.outpath)

if __name__ == "__main__":
	args = init_parser().parse_args()
	main(args)

