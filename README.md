# EPSplorer
Program which takes .fasta files of genomes (nucleotides) and detects exopolysaccharide (EPS) gene clusters, visualising them using gene arrows from gggenes. 

Clone the program, then add the genomes you want analysed to the genomes/ directory. Go intro scripts/EPSplorer.sh and edit the working directory (WD, line 3) to the EPSplorer directory and the amount of threads you want to use (threads, line 4). Now you can run the pipeline by writing "bash *path to EPSplorer.sh*".

The only part of any script that should be edited by the user are the two lines mentioned above.

You also need to add a magstats.tsv file, with the columns "bin" and "midas4tax". Under "bin" should be the names of your genome fasta files (without .fasta), and under "midas4tax" you should be the names you want them to have in the figures.
