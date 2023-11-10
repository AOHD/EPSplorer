# EPSplorer

Program which takes .fasta files of genomes (nucleotides) and detects exopolysaccharide (EPS) gene clusters, visualising them using gene arrows from gggenes. 

## How to install
The program requires Conda and InterProScan/5.38-76.0-foss-2020b
Clone the program, then add the genomes you want analysed to the genomes/ directory. 
Load Conda and InterProScan. Run the following commands:
conda env create --name R_env --file *path to R_env.yml*
conda create -n prokka_env prokka=1.14.6
conda create -n blast_env blast=2.12

Modify the magstats.tsv file. Under "bin" should be the names of your genome fasta files (without .fasta), and under "midas4tax" you should be the names you want them to have in the figures.

## How to run
Write source *path to EPSplorer.sh* -w *path to EPSplorer director* -t *amount of threads to be allocated (default 20)*
