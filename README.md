# EPSplorer

Program which takes .fasta files of genomes (nucleotides) and detects exopolysaccharide (EPS) gene clusters, visualising them using gene arrows from gggenes. 

## Installation

This program requires Conda and InterProScan/5.38-76.0-foss-2020b.

Clone the program, then add .fasta files of the genomes you want analysed to ```genomes/```.

Load Conda.

Run the following commands:

```
conda env create --name R_env --file *path to R_env.yml*
conda create -n prodigal_env -c bioconda prodigal=2.6.3
conda create -n blast_env -c bioconda blast=2.12
```

Modify the magstats.tsv file. Under "bin" should be the names of your genome fasta files (without .fasta), and under "midas4tax" should be the names you want them to have in the figures.

## Running the Program

```
source *path to EPSplorer.sh* -w <working_director (EPSplorer/)y> [-t <threads>] [-i <InterProScan>] [-h]"
```

## Output

This program detects putative EPS gene clusters using PSI-BLAST, displaying them in gene arrow plots. If the ```-i``` flag is supplied, InterProScan's Pfam library is used to highlight domains that are related to polysaccharide production. An overview of the different kinds of EPS gene clusters detected in the analysed genomes is supplied in an excel file. An Upset plot is also generated, depicting PSI-BLAST hits shared between different EPS queries. 
