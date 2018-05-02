#this script uses blastx and protein file
#	outputs statistics: 
#		Total number of basepairs
#		Number of transcripts
#		Number of "unigenes"
#		Mean and median transcript length
#		N50 explanation here
#		calculates the # of transcripts that have blast hits


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

ID_LIS = []
LEN_LIS = []
IDLIS = []

with open('stats.txt', 'w') as STATS: 

	for RECORD in SeqIO.parse(INFILE, 'fasta'):
			LEN_LIS.append(len(RECORD.seq))
	#THIS CALCULATES THE TOTAL BP LENGTH
	print('length = ' + str(sum(LEN_LIS)) + '.')
	STATS.write('full lengths = ' + str(sum(LEN_LIS)) + '.\n')
	
	#THIS CALCULATES THE LENGTH OF THE TRANSCRIPT 
	print ('transcripts = ' + str(len(LEN_LIS)) + '.')
	STATS.write('transcripts = ' + str(len(LEN_LIS)) + '.\n')
	
	#THIS CALCULATES THE MEAN OF THE LENGTH OF THE TRANSCRIPTS
	print('mean transcripts length = ' + str(stat.mean(LEN_LIS)) + '.')
	STATS.write('mean transcritps length = ' + str(stat.mean(LEN_LIS)) + '.\n')
	
	#THIS CALCULATES THE N50, THIS WAS COPIED FROM DARAY'S GITHUB
	LENLIS = sorted(LEN_LIS)
	HALF = sum(LEN_LIS)/2
	TOTALENTRIES = len(LEN_LIS)
	for ENTRY in range(1, TOTALENTRIES):
		RUNNINGTOTAL = sum(LEN_LIS[0:ENTRY])
		if RUNNINGTOTAL >= HALF:
			print('N50 = ' + str(LEN_LIS[ENTRY-1]) + '.')
			STATS.write('N50 = ' + str(LEN_LIS[ENTRY-1]) + '.\n')
			break
	for RECORD in SeqIO.parse(INFILE, 'fasta'):
		ID_LIS.append(RECORD.id.rsplit('_',1))
	###
	##THIS CALCULATES THE UNIGENES FROM THE FILES
	###
	IDS = pd.DataFrame(ID_LIS, columns = ['uni', 'iso'])
	COUNTUNI = len(IDS['uni'].unique())
	print('unigenes = ' + str(COUNTUNI))
	STATS.write('unigenes ' + str(COUNTUNI) + '.\n')

####
##PART TWO OF THE HW. THE IDEA HERE IS TO COUNT THE SEQS THAT HAVE HITS IN THE BLAST
###
	#reading the file as a table and use the columns to count the hits
	TRANSCRIPT_HITS = pd.read_table(BLASTX, sep='\t', usecols=[1,2], names=['ID', 'db'])
	HITS = len(TRANSCRIPT_HITS.ID.unique())
	print('UNIQUE HITS = ' + str(HITS) + '.')
	
#####
##PART THREE, NUMBER OF UNIQUE PROTEINS IN THE REFERENCE PROTEOME WITH BLAST HIT
#####
	##count the number of times there is a hit in the blast file