from __future__ import division
import glob
import os

# Place this file in the same directory as the results and seed folders
# (i.e., the folders starting with P00 that were created by automate_mitobim.py
# and the folder called seeds shoould be in the same directory as this program)
# To run: $ python make_results_table.py
# The script makes two files, seedLengths.txt and extensionLengths.txt,
# with seed and extended seed lengths respectively. 



def sorted_contig_lens(filename):
    with open(filename, 'U') as fh:
        lines = fh.readlines()
    nameslens = [(l1.strip().strip('>').split('_')[0],
                  len(l2.strip().replace('X', '').replace('N', '')))
                 for l1, l2 in zip(lines[0::2], lines[1::2])]
    nameslens = sorted(nameslens)
    return zip(*nameslens)


def count_seqs(filename):
    with open(filename, 'U') as fh:
        return len(fh.readlines()) // 4


def save_table(samples, is_seeds, filename):
    with open(filename, 'w') as fh:

        for sample in samples:
            if is_seeds:
                fasta = sample
                fastq = None
                sample = os.path.splitext(os.path.basename(sample))[0]
            else:
                fasta = os.path.join(sample, '4-MITObim', 'mt-candidate.fasta')
                fastq = os.path.join(sample, '4-MITObim', 'mt-candidate-readpool.fastq')

            names, lens = sorted_contig_lens(fasta)

            pool_size = count_seqs(fastq) if fastq is not None else 'seed'

            if sample in samples[0]:
                fh.write('\t'.join(['SampleName', 'Pool'] + list(names))+'\n')
            fh.write('\t'.join([sample, str(pool_size)] + [str(seqlen) for seqlen in lens])+'\n')


if __name__ == '__main__':
    import sys
    # if len(sys.argv) > 1:
    #     samples = sorted(sys.argv[1:])
    #     seeds = True
    # else:
    #     samples = sorted(glob.glob('P00*'))
    #     seeds = False
    samples = sorted(glob.glob("seeds/*.seeds"))
    save_table(samples, is_seeds=True, filename="seedLengths.txt")
    samples = sorted(glob.glob("P00*"))
    save_table(samples, is_seeds=False, filename="extensionLengths.txt")
