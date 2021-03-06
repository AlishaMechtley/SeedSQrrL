SeedSQrrL: A new tool for mitochondrial Genome Reconstruction in Non-Model Organisms
====================================================================================
The tool consists of three programs: 1) MitoDBMaker, and 2) MitoDBRelativeMaker, and 3) MitoDBExtractor. They are to be run in that order. The first program creates an SQLite database, populates it with the gene sequences and metadata (from the NCBI nucleotide and taxonomy databases) for a list of organisms that you provide then generates a list of genes for those organisms in which no match was found. If a list of genes is not provided via command line, a default list is used (12S, 16S, COX1, CYTB, ND2, ND5).
The next program, mitoDBRelativeMaker, will find relatives for the resulting list of unmatched genes using the NCBI Taxonomy DB as a guide and then populate the database with the genes for those relatives. The final program, mitoDBExtractor, extracts those genes into individual files to be run with mitobim. From there, the genes are extended using mitobim and contigs can be placed together (e.g. using Geneious). I have modified the jupyter notebook provided on the mitobim github page to be run with the files in the Anchor Lab in a script called automate_mitobim.py that is to be run after the extractor program.


### How to Install

Use Anaconda Python 2.7 distribution.

Install the necessary bs4 and biopython (for mitobim) modules using the following command in a terminal:
conda install beautifulsoup4 biopython


Copy the three SeedSQrrL programs and the automate_mitobim.py program into the folder where you are going to run them.


Install mitobim and its necessary packages (mira, flash, perl if not already installed but probably is). For mitobim, just copy MITObim_1.8.pl to /usr/local/bin or elsewhere in your $PATH and chmod a+x. For mira, you can download the platform-appropriate compiled binaries (mira_4.0.2_darwin…tar.bz2 for OSX, mira_4.0.2_linux-gnu...tar.bz2 for linux), and copy everything from the tar.bz2 file’s bin/ folder into the same place as mitobim above (e.g. /usr/local/bin). Copy everything from the mira library (lib) folder into usr/local/lib


Notes: use ls -l to make sure that the permissions are correct (-rwxr-xr-x). If there is an @ at the end of the permissions, use xattr -d com.apple.quarantine /usr/local/bin/filename on the file. This will likely need to be done on the mira executable and the library files. If you are using a recent version of mira (e.g. greater than 4.0.2) then mitobim needs to specify the bait and output files with a command. So you must add these to line 290 of the mitobim.pl file.
Example, replace


	@output = qx($mirabait -k $k_bait -n 1 temp_baitfile.fasta $readpool $strainname-readpool-it$currentiteration);



With:
	
	@output = qx($mirabait -k $k_bait -n 1 -b temp_baitfile.fasta $readpool -o $strainname-readpool-it$currentiteration);









### How to run

To Run SeedSQrrL:


1. Run mitoDBMaker
```$python mitoDBmaker.py [samplelist].csv [optional genelist]```


Example, let’s say you wanted ND4 in addition to the six default genes:
```python mitoDBmaker.py SampleList.csv [\'CO1\',\'ND2\',\'12S\',\'16S\',\'COX1\',\'ND5\',\'ND4\']```


2. Run mitoDBRelativeMaker
```$python mitoDBRelativemaker.py [samplelist]NeedReference.csv```


3. Create folder called “seeds”. Run mitoDBExtractor to Generate seed files.


 ![Figure1](https://raw.githubusercontent.com/AlishaMechtley/SeedSQrrL/master/images/InputOutput.png)

*Figure 1: SeedSQrrL Input/Output Workflow*




#### To run Mitobim
Organize your samples in a folder called RawReads, with individual samples underneath, in folders named named Sample_PXXXX_FG_IXXXX. When FG are any two letters. 
The (still gzipped) sample files go inside. For each sample, if there is a folder called redo inside, samples therein will be used instead.

Make sure all seeds are in a folder called seeds, and are named by the sample id (without the Sample_ part, should start with 'P00').

Run ```$python automate_mitobim.py```







### Troubleshooting

Be sure there are no spaces between items in the sample list.

If an organism is not found in the taxonomy database or a gene is not found, check that the spelling is correct.

Currently gene synonyms for six default genes are written into the program using information found in a human gene database. If there are a different set of genes you wish to use, you should add potential synonyms into the mitoDBmaker and mitoDBRelativeMaker programs. I would like to eventually automate this but I have to find a suitable database (e.g. not restricted to human mitochondrial gene synonyms).
Example:
if gene == "ND5":
gene_synonym = ["Mitochondrially Encoded NADH:Ubiquinone Oxidoreductase Core Subunit 5", "NADH Dehydrogenase Subunit 5", "EC 1.6.5.3", "MTND5", "ND5", "Mitochondrially Encoded NADH Dehydrogenase 5","NADH Dehydrogenase", Subunit 5 (Complex I)", "NADH-Ubiquinone Oxidoreductase Chain 5", "Complex I ND5 Subunit", "NADH Dehydrogenase 5", "NADH5"]


Notes on improvements:
Want options to extend one seed at a time (different files), select the highest rank, and restrict the use of bad references for seeds by id.

