import glob
import os

def save_extensions(samples):
	# read in the lines from each sample file
	for sample in samples:
		fasta = os.path.join(sample, '4-MITObim', 'mt-candidate.fasta')
		with open(fasta, 'U') as fh:
			lines = fh.readlines()
			for l1, l2 in zip(lines[0::2], lines[1::2]):
				# if line starts with >, read gene name before '_'
				if l1.find('>')!=-1:					
					gene = l1.split('_')[0] 
					#pop off > for gene name
					gene = gene[1:]
					# save sample name, gene name, newline, sequence in file under that gene
					filename = str(gene) + '.txt'
					with open(filename, 'a') as f:
						f.write('>'+str(sample.replace('_mt_candidate.fasta','_')) + '_' + str(gene)+'\n'+ str(l2))

        



if __name__ == '__main__':
	samples = sorted(glob.glob("P00*"))
	save_extensions(samples)
    
