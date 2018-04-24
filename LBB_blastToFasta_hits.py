
from Bio import SeqIO
import argparse

def get_args():
	parser = argparse.ArgumentParser(description="This will allow us to make two files, one with and one without hits", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-fi', '--fasta', type=str, help='peptide', required=True)
	parser.add_argument('-blast', '--blast', type=str, help='blastp', required=True)
	parser.add_argument('-out', '--out', help='Prefix', required=True)

	args = parser.parse_args()
	FASTA = args.fasta
	BLAST = args.blast
	PREFIX = args.out
	
	return FASTA, BLAST, PREFIX
FASTA, BLAST, PREFIX = get_args()

HITS_FILE = open(PREFIX + '.hits.fasta', 'w')
NOHITS_FILE = open(PREFIX + '.nohits.fasta', 'w')

FILES = open(BLAST).read() 
for RECORD in SeqIO.parse(FASTA, 'fasta'):
	if RECORD.id in FILES:
		SeqIO.write(RECORD, HITS_FILE, 'fasta')
	else: 
		SeqIO.write(RECORD, NOHITS_FILE, 'fasta')

HITS_FILE.close()
NOHITS_FILE.close()
