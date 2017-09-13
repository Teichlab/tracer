import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', help='GENCODE FASTA file to strip long names from and make a gene map out of')
args = parser.parse_args()

transcripts = []

with open(args.file,'r') as fid1:
	with open('transcripts.fasta','w') as fid2:
		with open('genemap.tsv','w') as fid3:
			for line in fid1:
				if line[0]=='>':
					#we have a transcript name on our hands, time to get parsing
					#only interested in ENST and ENSG, which are the first two fields
					line = line[1:].split('|')[:2]
					#not interested in the dot version number
					line[0] = line[0].split('.')[0]
					line[1] = line[1].split('.')[0]
					#so is this a duplicate of an already registered ENST as that happens sometimes
					if line[0] in transcripts:
						toggle = 0
					else:
						toggle = 1
						transcripts.append(line[0])
						#now we can make our cleaned up transcript name and our gene map entry
						fid2.write('>'+line[0]+'\n')
						fid3.write(line[0]+'\t'+line[1]+'\n')
				else:
					#copy over actual transcript sequence, unless we're skipping the transcript
					if toggle == 1:
						fid2.write(line)
