#!/usr/bin/env bash 
start_time=$(date +%s)

while getopts ":w:t:ih" opt; do
  case ${opt} in
    w)
      WD="$OPTARG"
      ;;
    t)
      threads="$OPTARG"
      ;;
    i)
      ips="TRUE"
      ;;
    h)
      echo "Usage: $(basename $0) -w <working_directory> [-t <threads>] [-i <InterProScan>] [-h]"
      return
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

# Check if WD is provided
if [ -z "$WD" ]; then
  echo "Error: Missing required option -w <working directory>"
  exit 1
fi

# Default threads if not provided
if [ -z "$threads" ]; then
  threads=20
fi

# Default threads if not provided
if [ -z "$ips" ]; then
  ips="FALSE"
fi





exec 2>"$ERROR_LOG"
cd $WD

if [ -e $WD/"temp" ]; then
  rm -r $WD/"temp"
fi

if [ -e $WD/"figures" ]; then
  rm -r $WD/"figures"
fi

if [ -e $WD/"data" ]; then
  rm -r $WD/"data"
fi

if [ -e $WD/"temp" ]; then
  rm -r $WD/"temp"
fi

ERROR_LOG=$WD/"error.log"
mkdir $WD/"temp"
TMPDIR=$WD/"temp"
mkdir $WD/"data"
mkdir $WD/"data/prodigal"
mkdir $WD/"data/databases"
mkdir $WD/"data/interproscan_results"
mkdir $WD/"data/output_proximity_filtration"
mkdir $WD/"data/psiblast_results"
mkdir $WD/"figures"
mkdir $WD/"figures/operons"

##Command which takes all the .faa files in the genomes/ folder and removes all underscores from the file
for file in $WD/genomes/*.fasta; do sed -i "s/_//g" $file; done

source $WD/scripts/magstats.sh

##Command which runs prodigal on all .fasta files of $WD/prodigal/
echo Annotating genomes with prodigal
conda activate prodigal_env
for file in $WD/genomes/*.fasta; do mkdir $WD/data/prodigal/$(basename $file .fasta)/; prodigal -a $WD/data/prodigal/$(basename $file .fasta)/$(basename $file .fasta).faa -c -f gff -o $WD/data/prodigal/$(basename $file .fasta)/$(basename $file .fasta).gff -i $file; done
conda deactivate

##Command which takes all the .faa files in the prodigal output folder and changes the > + 8 letter string before each prodigal ID to the name of the .fasta file
for file in $WD/genomes/*.fasta; do sed -i "s/>[^_]*_/>$(basename $file .fasta)_/g" $WD/data/prodigal/$(basename $file .fasta)/*.faa; done

##Command which takes all the .faa files in the prodigal output folder and removes all asterisks from the file
for file in $WD/genomes/*.fasta; do sed -i "s/\*//g" $WD/data/prodigal/$(basename $file .fasta)/*.faa; done

##Command which takes all the .gff files in the prodigal output folder and changes the ID= + 8 letter string before each prodigal ID to the name of the .fasta file
for file in $WD/genomes/*.fasta; do sed -i "s/ID=[^_]*_/ID=$(basename $file .fasta)_/g" $WD/data/prodigal/$(basename $file .fasta)/*.gff; done

##Makeblastdb command which makes a blast database for each .fasta file in $WD/data/prodigal/, 
##based on the .faa file in the prodigal output folder. Each database should be in its own folder in $WD/databases/
echo Making blast databases
conda activate blast_env
for file in $WD/genomes/*.fasta; do makeblastdb -in $WD/data/prodigal/$(basename $file .fasta)/*.faa -dbtype prot -out $WD/data/databases/$(basename $file .fasta)/$(basename $file .fasta); done


#Execute psiblast.sh
echo Running PSI-BLAST
data="$WD/data/databases/*"

##Query path definitions, non-MSA queries
operon_fasta=$WD/queries_psiblast/*.faa

## PSI-BLAST of all operons against all databases in $WD/queries_psiblast/databases/
for database in ${data[@]}; do
mkdir $WD/data/psiblast_results/$(basename $database)
for operon in ${operon_fasta[@]}; do
psiblast -query $operon -db $database/$(basename $database) -out $WD/data/psiblast_results/$(basename $database)/$(basename $operon) -evalue 0.0001 -qcov_hsp_perc 20 -max_hsps 10 -max_target_seqs 100000 -outfmt 6 -num_iterations 20 -comp_based_stats 1 -num_threads $threads
done
echo $database completed
done
conda deactivate

##Move all .faa files in $WD/data/prodigal/*/* to $WD/data/prodigal/
for file in $WD/data/prodigal/*/*.faa; do
cp $file $WD/data/prodigal
echo $file
done


##Concatenate all .faa files with the same name in $WD/data/psiblast_results/*/*.faa
echo Creating concatenated GFF file
mkdir $WD/data/psiblast_results/concatenated
for file in $WD/data/psiblast_results/*/*.faa; do cat $file >> $WD/data/psiblast_results/concatenated/$(basename $file .faa).faa; done

##Run $WD/scripts/generate_gff.R
conda activate R_env
Rscript $WD/scripts/generate_gff.R "$WD"


##Run $WD/scripts/proximity_main.R
echo Running proximity filtration
Rscript $WD/scripts/proximity_main.R "$WD"

conda deactivate

#If data/output_proximity_filtration/fasta_output is empty, exit with error
if [ -z "$(ls -A $WD/data/output_proximity_filtration/fasta_output)" ]; then
  echo "Error: No polysaccharides found in proximity filtration"
  exit 1
fi

##InterProScan analysis
if [ "$ips" = "TRUE" ]; then
echo Running InterProScan
# Folder with all subsetted fasta files
FASTA=$WD/data/output_proximity_filtration/fasta_output
# Results folder
RESULTS=$WD/data/interproscan_results

for j in $FASTA/*; do
# For each polysaccharide
j=`basename $j`
echo $j
mkdir -p $RESULTS/$j
mkdir -p $RESULTS/$j\_fasta_removed
for i in $FASTA/$j/*; do
	# For each MAG id
	# Removing .faa from the name
	i=`basename $i | rev | cut -c 5- | rev`
	echo $i
	if [[ ! -s $RESULTS/$j\_fasta_removed/$i.gff3 ]]; then
	# Running interproscan
	## Only Pfam database
	## Only GFF3 file is written out
	source interproscan.sh -dp -T $TMPDIR -cpu $threads -appl Pfam -f GFF3 -o $RESULTS/$j/$i.gff3 -i $FASTA/$j/$i.faa
	
	# Remove trailing fasta sequence in gff3 file and leading three lines
	F=$RESULTS/$j/$i.gff3
	sed '/^##FASTA$/,$d' ${F} | sed '1,/^##interproscan*/d' > $RESULTS/$j\_fasta_removed/$i.gff3
	fi
done
done
fi

conda activate R_env

##Run $WD/scripts/ips_main.R
echo Generating gene arrow plots
Rscript $WD/scripts/ips_main.R "$WD" "$(printf "%q" "$ips")"

echo Generating overview excel file
Rscript $WD/scripts/overview.R "$WD"

cd $WD

exec 2>&1

conda deactivate
unset ips
unset threads


end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo Finished in $elapsed_time seconds! Have a nice day!



