# EPSplorer

Program which detects exopolysaccharide (EPS) gene clusters in bacterial genomes. 

## Installation

This program requires Mamba.

Clone the program, then add .fasta files of the genomes you want analysed to ```genomes/```. Please only submit .fasta files with a single sequence in them.

Make sure Mamba is installed and loaded.

Run the following commands:

```
mamba env create --name R_env --file *path to R_env.yml*
mamba create -n prodigal_env -c bioconda prodigal=2.6.3
mamba create -n blast_env -c bioconda blast=2.12
mamba create -n ips_v5.59_91 -c bioconda interproscan 
```

Once the conda environments are created, you do not have to create them on subsequent runs. The InterProScan conda environment does not work out the gate, and requires you to run the script described by andradejon in this thread to work:
https://github.com/ebi-pf-team/interproscan/issues/305

In the scripts/ folder there is an .sh file configured for slurm to run the andradejon script with sbatch on the CMC servers.

## Running the Program

```
source *path to EPSplorer.sh* -w <working_directory> [-t <threads>] [-i <InterProScan>] [-h]"
```

## Output

This program detects putative EPS gene clusters using PSI-BLAST, displaying them in gene arrow plots. If the ```-i``` flag is supplied, InterProScan's Pfam library is used to highlight domains that are related to polysaccharide production. 

An overview of the different kinds of EPS gene clusters detected in the analysed genomes is supplied as an excel file. An Upset plot is also generated, depicting PSI-BLAST hits shared between different EPS queries. 

If one wishes to do downstream analysis of the detected gene clusters, useful information can be found in ```data/output_proximity_filtration```. Here, .faa files for all detected EPS genes can be found (```fasta_output/```), as well as .tsv files containing information about individual gene clusters and their genes  (```psi_operon_full/```)
