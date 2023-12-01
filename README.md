# EPSplorer

Program which detects exopolysaccharide (EPS) gene clusters in bacterial genomes. 

## Installation

This program requires Conda and InterProScan/5.38-76.0-foss-2020b.

Clone the program, then add .fasta files of the genomes you want analysed to ```genomes/```. Please only submit .fasta files with a single sequence in them.

Make sure InterProScan and Conda are installed and loaded.

Run the following commands:

```
conda env create --name R_env --file *path to R_env.yml*
conda create -n prodigal_env -c bioconda prodigal=2.6.3
conda create -n blast_env -c bioconda blast=2.12
```

Create the magstats.tsv file, which should have two columns, "bin" and midas4_tax. Under "bin" should be the names of your genome fasta files (without .fasta), and under "midas4tax" should be the names you want them to have in the figures. 

To help create magststs.tsv, you can use the ```scripts/magstats.sh ``` script, which will create a magstats file based on the filenames in the ```genomes/``` folder. Be sure to run it from the ```EPSplorer/``` directory.

## Running the Program

```
source *path to EPSplorer.sh* -w <working_director (EPSplorer/)y> [-t <threads>] [-i <InterProScan>] [-h]"
```

## Output

This program detects putative EPS gene clusters using PSI-BLAST, displaying them in gene arrow plots. If the ```-i``` flag is supplied, InterProScan's Pfam library is used to highlight domains that are related to polysaccharide production. 

An overview of the different kinds of EPS gene clusters detected in the analysed genomes is supplied as an excel file. An Upset plot is also generated, depicting PSI-BLAST hits shared between different EPS queries. 

If one wishes to do downstream analysis of the detected gene clusters, useful information can be found in ```data/output_proximity_filtration```. Here, .faa files for all detected EPS genes can be found (```fasta_output/```), as well as .tsv files containing information about individual gene clusters and their genes  (```psi_operon_full/```)
