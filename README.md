# EPSplorer

Program which detects exopolysaccharide (EPS) gene clusters in bacterial genomes. 

## Installation

This program requires Conda and InterProScan/5.38-76.0-foss-2020b.

Clone the program, then add .fasta files of the genomes you want analysed to ```genomes/```. Please only submit .fasta files with a single sequence in them.

<<<<<<< HEAD
Make sure InterProScan and Conda are installed and loaded.
=======
Make sure InterProScan/5.38-76.0-foss-2020b and Conda  is installed and loaded.
>>>>>>> b1a084a4e46fab78f8b7379540aed56ee6ef9829

Run the following commands:

```
conda env create --name R_env --file *path to R_env.yml*
conda create -n prodigal_env -c bioconda prodigal=2.6.3
conda create -n blast_env -c bioconda blast=2.12
```

Once the conda environments are created, you do not have to create them on subsequent runs.

## Running the Program

```
source *path to EPSplorer.sh* -w <working_director (EPSplorer/)y> [-t <threads>] [-i <InterProScan>] [-h]"
```

## Output

This program detects putative EPS gene clusters using PSI-BLAST, displaying them in gene arrow plots. If the ```-i``` flag is supplied, InterProScan's Pfam library is used to highlight domains that are related to polysaccharide production. 

An overview of the different kinds of EPS gene clusters detected in the analysed genomes is supplied as an excel file. An Upset plot is also generated, depicting PSI-BLAST hits shared between different EPS queries. 

<<<<<<< HEAD
If one wishes to do downstream analysis of the detected gene clusters, useful information can be found in ```data/output_proximity_filtration```. Here, .faa files for all detected EPS genes can be found (```fasta_output/```), as well as .tsv files containing information about individual gene clusters and their genes  (```psi_operon_full/```)
=======
If one wishes to do downstream analysis of the detected gene clusters, useful information can be found½ in ```data/output_proximity_filtration```. Here, .faa files for all detected EPS genes can be found (```fasta_output/```), as well as .tsv files containing information about individual gene clusters and their genes  (```psi_operon_full/```)
>>>>>>> b1a084a4e46fab78f8b7379540aed56ee6ef9829
