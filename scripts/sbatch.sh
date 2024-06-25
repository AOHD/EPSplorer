#!/usr/bin/bash -l
#SBATCH --job-name=EPSplorer_Eawag
#SBATCH --output=/home/bio.aau.dk/zs85xk/sbatch/job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=default-op
#SBATCH --cpus-per-task=192
#SBATCH --mem=288G
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aohd@bio.aau.dk

# Exit on first error and if any variables are unset
set -eu
module purge

ips="TRUE"
# load software modules or environments
#ips="FALSE"
# run one or more commands as part a full pipeline script,
# or call scripts from elsewhere. Make sure you use the same
# number of max threads as requested 
cd /home/bio.aau.dk/zs85xk/projects/EPSplorer
source EPSplorer.sh -w . -t 192 -i