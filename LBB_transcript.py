#this script uses blastx and protein file
#	outputs statistics: 
#		Total number of basepairs
#		Number of transcripts
#		Number of "unigenes"
#		Mean and median transcript length
#		N50 explanation here

from Bio import SeqIO
import statistics as stat
import pandas as pd
import argparse
def get_args():
	parser = argparse.ArgumentParser(description="this parses transcripts", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-in', '--infile', type=str, help='Insert the fasta pep file', required=True)
	parser.add_argument('-blast', '--blastx', type=str, help='Insert the blastx file', required=True)
	
	args = parser.parse_args()
	INFILE = args.infile
	BLASTX = args.blastx
	
	return INFILE, BLASTX
INFILE, BLASTX = get_args()

ID_LIST = []
LENLIST = []
IDLIST = []

with open('stats.txt', 'w') as STATS: 

	for RECORD in SeqIO.parse(INFILE, 'fasta'):
			LENLIST.append(len(RECORD.seq))
	#this calculates the total basepairs in length
	print('length = ' + str(sum(LENLIST)) + '.')
	STATS.write('full lengths = ' + str(sum(LENLIST)) + '.\n')
	
	#this calculates the leght of the transcript 
	print ('transcripts = ' + str(len(LENLIST)) + '.')
	STATS.write('transcripts = ' + str(len(LENLIST)) + '.\n')
	
	#this calculates the mean of the length of the transcript
	print('mean transcripts length = ' + str(stat.mean(LENLIST)) + '.')
	STATS.write('mean transcritps length = ' + str(stat.mean(LENLIST)) + '.\n')
	
	#this calculates the N50, this is the mean or median of legths
	print('N50 of transcripts = ' + str(stat.median(LENLIST)) + '.\n')
	
	for RECORD in SeqIO.parse(INFILE, 'fasta'):
		IDLIST.append(RECORD.id.rsplit('_',1))
		
IDS = pd.DataFrame(IDLIST, columns = ['uni', 'iso'])
COUNTUNI = len(IDS['uni'].unique())
print('we have ' + str(COUNTUNI) + ' unigenes in the file.')
